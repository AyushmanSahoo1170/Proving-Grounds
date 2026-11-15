#!/usr/bin/env python3
"""Proving Grounds: rank Earth terrain windows against Moon/Mars targets, backtest, export for the map.
  python proving_grounds.py selftest                          # synthetic terrain; tests the logic; writes a bannered demo
  python proving_grounds.py run                               # real run: downloads and clips everything it needs (internet + rasterio)
Targets live in data/moon.tif, mars.tif, mars2.tif (elevation in metres, ~25 km windows). Existing files are reused; --moon/--mars/--mars2 override."""
import argparse, json, math, os, sys, time, urllib.request, warnings
import numpy as np
from scipy.ndimage import uniform_filter
from scipy.stats import wasserstein_distance
warnings.filterwarnings("ignore", category=RuntimeWarning)

WIN_KM, MIN_RES, THRESH, MIN_BARE = 25.0, 90.0, 0.7, 0.5
MOON_URL = "https://pgda.gsfc.nasa.gov/data/LOLA_5mpp/Site01/Site01_final_adj_5mpp_surf.tif"   # NASA PGDA, Connecting Ridge, 5 m LOLA
MOLA_URL = "https://planetarymaps.usgs.gov/mosaic/Mars_MGS_MOLA_DEM_mosaic_global_463m.tif"      # USGS, MOLA 463 m global
BINS = [0, .5, 1, 2, 3, 5, 8, 12, 20, 90]; CEN = np.array([.25, .75, 1.5, 2.5, 4, 6.5, 10, 16, 40.0])
CRIT = ["Slope distribution", "Roughness", "Relief"]
TARGETS = {  # weights are assumptions set from the task, not fitted to known sites
    "moon": {"name": "Shackleton Connecting Ridge (Moon)", "activity": "Rover mobility", "w": [.5, .3, .2], "R": 1737.4e3},
    "mars": {"name": "Arcadia Planitia (Mars)", "activity": "Ice access", "w": [.3, .3, .4], "R": 3396.19e3}}
TARGETS["moon"]["at"] = None; TARGETS["mars"]["at"] = (47.2, 184.3)   # Arcadia Planitia centre (USGS gazetteer)
TARGETS["mars2"] = {"name": "Deuteronilus Mensae (Mars)", "activity": "Ice access", "w": [.3, .3, .4], "R": 3396.19e3, "at": (43.8, 29.4)}  # HiRISE-imaged lobate debris aprons
KNOWN = {"Atacama (Yungay)": (-24.1, -69.9), "Death Valley": (36.5, -116.9), "Haughton Crater, Devon Island": (75.4, -89.7),
         "Mars Desert Research Station, Utah": (38.4, -110.8), "Mauna Kea, Hawaii": (19.8, -155.5),
         "McMurdo Dry Valleys": (-77.5, 162.0), "Askja, Iceland": (65.0, -16.8), "Meteor Crater, Arizona": (35.0, -111.0)}
CAND = {"Dasht-e Lut, Iran": (30.5, 59.0), "Danakil Depression, Ethiopia": (14.2, 40.3), "Namib Desert, Namibia": (-24.5, 15.5),
        "Tibesti, Chad": (21.0, 17.5), "Altiplano-Puna, Argentina": (-25.5, -67.0), "Craters of the Moon, Idaho": (43.2, -113.5),
        "Lonar crater, India": (19.98, 76.51), "Changthang plateau, Ladakh": (33.5, 78.5), "Qaidam Basin, China": (37.0, 93.5),
        "Pinacate volcanic field, Mexico": (31.8, -113.5), "Svalbard": (78.2, 16.0), "Axel Heiberg Island, Canada": (79.4, -90.8),
        "Rub' al Khali, Saudi Arabia": (20.0, 50.0), "Amboy lava field, California": (34.55, -115.78)}
LIMITS = ["Terrain only (slope, roughness, relief). Thermal, mineral and illumination layers are not included.",
          "Windows with under 50% bare or sparse ground (ESA WorldCover 2021) are dropped. Where WorldCover has no tile, windows are kept unmasked.",
          "Mars targets use the 463 m global MOLA model, so Earth is compared with them at that coarser scale.",
          "Earth cannot match vacuum, low gravity, radiation or extreme cold.",
          "Activity weights are assumptions set from the task, not fitted to known sites.",
          "The backtest ranks about 8 known analogs among desert-like regions only, so it shows consistency, not proof.",
          "Scores rank sites for field checking. They are candidates, not verified analogs."]

def agg(z, fx, fy):
    h, w = (z.shape[0] // fy) * fy, (z.shape[1] // fx) * fx
    return np.nanmean(z[:h, :w].reshape(h // fy, fy, w // fx, fx), axis=(1, 3))

def agg_to(z, dx, dy, res):
    fx, fy = max(1, round(res / dx)), max(1, round(res / dy))
    return agg(z, fx, fy), dx * fx, dy * fy

def feats(z, dx, dy):
    z = np.where(np.isnan(z), np.nanmean(z), z)
    gy, gx = np.gradient(z, dy, dx)
    sl = np.degrees(np.arctan(np.hypot(gx, gy)))
    h = np.histogram(sl, bins=BINS)[0].astype(float); h /= h.sum()
    p5, p95 = np.percentile(z, [5, 95])
    return {"sh": h, "sm": float(sl.mean()), "rough": float((z - uniform_filter(z, 3)).std()), "relief": float(p95 - p5)}

def parts(a, b, tau):  # three 0-1 similarities: slope distribution (Earth Mover's distance), roughness ratio, relief ratio
    ratio = lambda x, y: (min(x, y) + .1) / (max(x, y) + .1)
    return [math.exp(-wasserstein_distance(CEN, CEN, a["sh"], b["sh"]) / tau), ratio(a["rough"], b["rough"]), ratio(a["relief"], b["relief"])]

def backtest(W, n=2000, seed=0):
    """Region-level test. Where do known analogs rank among regions, and how does that compare with chance?
    (Comparing a region's BEST window with all windows is biased: the best of ~10 windows lands near the 90th percentile by luck.)"""
    regs = {}
    for w in W: regs.setdefault(w["roi"], []).append(w["score"])
    names = list(regs); best = np.array([max(regs[k]) for k in names])
    known = {w["roi"] for w in W if w["kind"] == "known"}; kidx = [i for i, k in enumerate(names) if k in known]
    if not kidx or len(names) < 3: return None
    rank = lambda b: (-b).argsort(kind="stable").argsort(kind="stable") + 1
    obs, rng = float(np.median(rank(best)[kidx])), np.random.default_rng(seed)
    allsc, cuts = np.array([w["score"] for w in W]), np.cumsum([len(regs[k]) for k in names])[:-1]
    null = np.array([np.median(rank(np.array([q.max() for q in np.split(rng.permutation(allsc), cuts)]))[kidx]) for _ in range(n)])
    return {"median_rank": obs, "n_regions": len(names), "n_known": len(kidx), "expected_rank": round((len(names) + 1) / 2, 1), "p": round(float((null <= obs).mean()), 3)}

def analyse(rois, tz):
    out, rng = {}, np.random.default_rng(0)
    for b in tz:
        T = TARGETS[b]; z, dx, dy = tz[b]; res = max(MIN_RES, (dx + dy) / 2)
        Zt, ex, ey = agg_to(z, dx, dy, res); tf = feats(Zt, ex, ey); W = []
        tau = max(0.5, 0.35 * tf["sm"])   # slope tolerance in degrees scales with how steep the target is
        for r in rois:
            Z, ex, ey = agg_to(r["z"], r["dx"], r["dy"], res); fx, fy = ex / r["dx"], ey / r["dy"]
            B = agg_to(r["bare"].astype("float32"), r["dx"], r["dy"], res)[0] if r.get("bare") is not None else None
            wp = max(16, int(WIN_KM * 1000 / ((ex + ey) / 2)))
            for i in range(Z.shape[0] // wp):
                for j in range(Z.shape[1] // wp):
                    s = Z[i * wp:(i + 1) * wp, j * wp:(j + 1) * wp]
                    if np.isnan(s).mean() > .1: continue
                    bare = None if B is None else float(np.nanmean(B[i * wp:(i + 1) * wp, j * wp:(j + 1) * wp]))
                    if bare is not None and bare < MIN_BARE: continue
                    f = feats(s, ex, ey)
                    W.append({"roi": r["name"], "kind": r["kind"], "truth": r.get("truth"), "bare": None if bare is None else round(bare, 2),
                              "lat": round(r["top"] - (i * wp + wp / 2) * fy * r["dlat"], 3),
                              "lon": round(r["left"] + (j * wp + wp / 2) * fx * r["dlon"], 3), "f": f, "p": parts(tf, f, tau)})
        S, w0 = np.array([w["p"] for w in W]), np.array(T["w"])
        base = S @ w0 / w0.sum()
        Wm = w0 * rng.uniform(.8, 1.2, (200, 3))                       # weights varied by +/-20 percent
        ranks = (-(S @ Wm.T / Wm.sum(1))).argsort(0).argsort(0) + 1
        lo, hi = np.percentile(ranks, [5, 95], axis=1)
        for k, w in enumerate(W):
            w.update(score=round(float(base[k]), 3), agree=int((S[k] >= THRESH).sum()), lo=int(lo[k]), hi=int(hi[k]))
        best = {}
        for w in W: best[w["roi"]] = max(best.get(w["roi"], 0), w["score"])
        order = sorted(best, key=best.get, reverse=True)
        known = [{"name": n, "best": best[n], "pct": round(float((base < best[n]).mean() * 100), 1), "rank": order.index(n) + 1}
                 for n in order if any(w["roi"] == n and w["kind"] == "known" for w in W)]
        out[b] = {"windows": W, "validation": {"fair": backtest(W), "n_roi": len(order), "n_win": len(W), "known": known,
                  "median_pct": round(float(np.median([k["pct"] for k in known])), 1) if known else None,
                  "hit25": round(float(np.mean([k["pct"] >= 75 for k in known])), 2) if known else None},
                  "target": {"name": T["name"], "activity": T["activity"], "w": T["w"], "res_m": round(res), "tau": round(tau, 2),
                             "f": tf}}
    return out

def clean(o):
    if isinstance(o, dict): return {k: clean(v) for k, v in o.items() if k != "truth"}
    if isinstance(o, list): return [clean(v) for v in o]
    if isinstance(o, np.ndarray): return [round(float(x), 3) for x in o]
    if isinstance(o, (float, np.floating)): return round(float(o), 3)
    return o

def save(res, synthetic, out):
    doc = {"synthetic": synthetic, "generated": time.strftime("%Y-%m-%d %H:%M"), "win_km": WIN_KM, "crit": CRIT, "limits": LIMITS, "results": clean(res)}
    s = json.dumps(doc, separators=(",", ":"))
    open(out, "w").write(s); open(out.replace(".json", ".js"), "w").write("window.RESULTS=" + s + ";")
    print("wrote", out, "and", out.replace(".json", ".js"))

# ---------- real data ----------
def load_target(path, body, cyl=False):
    import rasterio
    with rasterio.open(path) as s:
        z, t, nd = s.read(1).astype("float32"), s.transform, s.nodata
        wkt = s.crs.to_wkt().lower() if s.crs else ""
        geo = s.crs.is_geographic if s.crs else abs(t.a) < 1   # no CRS: pixel size under 1 means degrees, otherwise metres
    if nd is not None: z[z == nd] = np.nan
    R = TARGETS[body]["R"]
    if geo:
        lat = t.f + t.e * z.shape[0] / 2
        return z, abs(t.a) * math.pi / 180 * R * math.cos(math.radians(lat)), abs(t.e) * math.pi / 180 * R
    if cyl or any(k in wkt for k in ("cylindrical", "equirect", "eqc", "plate_carree")):  # east-west pixels shrink away from the equator
        lat = math.degrees((t.f + t.e * z.shape[0] / 2) / R)
        return z, abs(t.a) * math.cos(math.radians(lat)), abs(t.e)
    return z, abs(t.a), abs(t.e)

def bare_mask(la, lo, shape):  # ESA WorldCover 10 m 2021 (3-degree COGs, public bucket): True where bare, sparse or snow/ice
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.windows import from_bounds
    la3, lo3, fl, fo = 3 * math.floor(la / 3), 3 * math.floor(lo / 3), math.floor(la), math.floor(lo)
    n = f"ESA_WorldCover_10m_2021_v200_{'N' if la3 >= 0 else 'S'}{abs(la3):02d}{'E' if lo3 >= 0 else 'W'}{abs(lo3):03d}_Map.tif"
    with rasterio.open(f"/vsicurl/https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/{n}") as s:
        lc = s.read(1, window=from_bounds(fo, fl, fo + 1, fl + 1, s.transform), out_shape=shape, resampling=Resampling.nearest)
    return np.isin(lc, (60, 70))

def load_roi(name, kind, lat, lon):  # Copernicus GLO-30 COG tile (public AWS bucket, no login)
    import rasterio
    la, lo = math.floor(lat), math.floor(lon)
    n = f"Copernicus_DSM_COG_10_{'N' if la >= 0 else 'S'}{abs(la):02d}_00_{'E' if lo >= 0 else 'W'}{abs(lo):03d}_00_DEM"
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR"):
        with rasterio.open(f"/vsicurl/https://copernicus-dem-30m.s3.amazonaws.com/{n}/{n}.tif") as s:
            z, t = s.read(1).astype("float32"), s.transform
    z[z < -500] = np.nan
    z[z == 0] = np.nan   # GLO-30 stores ocean as exactly 0 m; masking it keeps coastlines from looking flat
    bare = None
    try: bare = bare_mask(lat, lon, z.shape)
    except Exception as e: print("   (no land-cover mask here:", type(e).__name__ + ")")
    return {"name": name, "kind": kind, "z": z, "bare": bare, "dx": abs(t.a) * 111320 * math.cos(math.radians(lat)), "dy": abs(t.e) * 110574,
            "top": t.f, "left": t.c, "dlat": abs(t.e), "dlon": abs(t.a)}

def clip_mars(lat, lon, out, half_km=12.5):  # reads only a small window of the 2 GB USGS file over HTTP
    import rasterio
    from rasterio.windows import Window
    R = 3396000.0
    with rasterio.open("/vsicurl/" + MOLA_URL) as s:
        t = s.transform
        if abs(s.bounds.left + math.pi * R) > 5000: print("  warning: unexpected file extent", s.bounds, "- check the Mars clip")
        cx, cy, hy = R * math.radians(lon - 180), R * math.radians(lat), half_km * 1000
        hx = hy / math.cos(math.radians(lat))
        c0, r0 = int((cx - hx - t.c) / t.a), int((cy + hy - t.f) / t.e)
        c1, r1 = int((cx + hx - t.c) / t.a) + 1, int((cy - hy - t.f) / t.e) + 1
        w = Window(c0, r0, c1 - c0, r1 - r0)
        z, prof = s.read(1, window=w), s.profile.copy()
        prof.update(height=z.shape[0], width=z.shape[1], transform=s.window_transform(w), tiled=False, compress="deflate")
        for k in ("blockxsize", "blockysize"): prof.pop(k, None)
    with rasterio.open(out, "w", **prof) as d: d.write(z, 1)

def fetch(a=None):
    os.makedirs("data", exist_ok=True)
    if not os.path.exists("data/moon.tif"):
        print("Downloading the Moon target (NASA PGDA Site 01, Connecting Ridge)...")
        try:
            with urllib.request.urlopen(urllib.request.Request(MOON_URL, headers={"User-Agent": "proving-grounds/1.0"})) as r, open("data/moon.tmp", "wb") as f:
                n = 0
                while chunk := r.read(1 << 20): f.write(chunk); n += len(chunk); print(f"\r  {n >> 20} MB", end="", flush=True)
            os.replace("data/moon.tmp", "data/moon.tif"); print()
        except Exception as e: print("\n  Moon download failed:", e, "\n  Save it by hand as data/moon.tif from", MOON_URL)
    for b in ("mars", "mars2"):
        if not os.path.exists(f"data/{b}.tif"):
            print("Clipping", TARGETS[b]["name"], "from the USGS MOLA model...")
            try: clip_mars(*TARGETS[b]["at"], f"data/{b}.tif")
            except Exception as e: print("  automatic clip failed:", e, "\n  Use the manual QGIS or gdal_translate steps in the README.")

def run(a):
    fetch(); tz = {}
    for b in TARGETS:
        p = getattr(a, b, None) or f"data/{b}.tif"
        if os.path.exists(p): tz[b] = load_target(p, b, a.mars_cyl and b != "moon"); print("target ready:", b)
        else: print("skipping target", b, "(no file)")
    if not tz: sys.exit("No target files found. See README, Real run.")
    rois = []
    for kind, d in (("known", KNOWN), ("candidate", CAND)):
        for n, (la, lo) in d.items():
            try: rois.append(load_roi(n, kind, la, lo)); print("loaded", n)
            except Exception as e: print("skipped", n, "->", e)
    if not rois: sys.exit("No Earth tiles could be downloaded. Check your internet connection.")
    save(analyse(rois, tz), False, a.out)

def revalidate(a):  # recompute the backtest from an existing results file (no downloads, takes seconds)
    d = json.load(open(a.out))
    for b, r in d["results"].items():
        r["validation"]["fair"] = backtest(r["windows"]); print(b, r["validation"]["fair"])
    d["limits"] = LIMITS; s = json.dumps(d, separators=(",", ":"))
    open(a.out, "w").write(s); open(a.out.replace(".json", ".js"), "w").write("window.RESULTS=" + s + ";"); print("updated", a.out)

# ---------- self-test on synthetic terrain ----------
def fractal(n, beta, amp, rng):
    f = np.fft.fftfreq(n); k = np.hypot(*np.meshgrid(f, f)); k[0, 0] = 1
    z = np.fft.ifft2(np.fft.fft2(rng.standard_normal((n, n))) * k ** (-beta / 2)).real
    return amp * (z - z.mean()) / z.std()

CLASSES = {"plain": (3.2, 6), "hills": (2.8, 120), "rugged": (2.2, 500)}

def selftest(a):
    rng, names = np.random.default_rng(1), list(CLASSES)
    tz = {"moon": (fractal(208, *CLASSES["rugged"], rng), 118., 118.), "mars": (fractal(208, *CLASSES["plain"], rng), 120., 120.), "mars2": (fractal(208, *CLASSES["hills"], rng), 120., 120.)}
    rois = []
    for i, (n, (la, lo), kind) in enumerate([(n, c, "known") for n, c in KNOWN.items()] + [(n, c, "candidate") for n, c in CAND.items()]):
        cls = names[i % 3]
        rois.append({"name": "SYNTHETIC " + n, "kind": kind, "truth": cls, "z": fractal(832, *CLASSES[cls], rng), "dx": 120., "dy": 120.,
                     "top": la + .5, "left": lo - .5, "dlat": 120 / 110574, "dlon": 120 / (111320 * math.cos(math.radians(la)))})
    res = analyse(rois, tz)
    for b, want in (("moon", "rugged"), ("mars", "plain"), ("mars2", "hills")):
        top = sorted(res[b]["windows"], key=lambda w: -w["score"])[:8]
        hits = sum(w["truth"] == want for w in top); print(b, "top-8 windows that are the matching class:", hits, "/ 8")
        assert hits >= 6, "scoring failed to recover the matching terrain class"
    print("self-test passed"); save(res, True, a.out)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest="cmd", required=True)
    for c in ("selftest", "run", "fetch", "revalidate"):
        p = sp.add_parser(c); p.add_argument("--out", default="web/results.json")
    r = sp.choices["run"]
    for b in TARGETS: r.add_argument("--" + b, help=f"override data/{b}.tif")
    r.add_argument("--mars-cyl", action="store_true", help="force the simple-cylindrical correction for Mars DEMs")
    a = ap.parse_args(); {"run": run, "fetch": fetch, "selftest": selftest, "revalidate": revalidate}[a.cmd](a)
