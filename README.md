# Proving Grounds

**Live viewer:** https://proving-grounds-gamma.vercel.app/

Finds places on Earth whose terrain matches NASA's Moon Base target (Shackleton Connecting Ridge) and two Mars ice-access targets (Arcadia Planitia and Deuteronilus Mensae), then backtests the method against known analog sites. Built for NASA Space Apps 2026: *Identify Earth Locations that Analog the Permanent Moon Base Locations and Mars*.

The results in `web/` are real. `python proving_grounds.py run` generated them on 15 November 2026 (22:39). The viewer shows no warning banner because `web/results.json` has `"synthetic": false`.

## What is in this folder
- `proving_grounds.py`  the pipeline: downloads data, scores, backtests, exports
- `web/index.html`  the map viewer (Leaflet; Esri World Imagery and OpenStreetMap basemaps, switchable in the layer control)
- `web/results.js`, `web/results.json`  the real results from the 15 November 2026 run
- `requirements.txt`  numpy, scipy, rasterio
- `vercel.json`, `.vercelignore`  deploy only the static viewer in `web/`
- `data/`  created on the first real run (target elevation files), not committed

## Results
Every 25 km Earth window gets a 0-1 score against each target. 22 Earth regions are defined (8 known analog sites and 14 candidates). 17 regions (Moon) and 16 regions (Mars) produced at least one scored window.

| Target | Windows scored | Best new candidates (score, rank range when weights vary by 20%) |
|---|---|---|
| Shackleton Connecting Ridge (Moon), rover mobility | 187 | Altiplano-Puna, Argentina (0.502, rank 1-4); Dasht-e Lut, Iran (0.493, rank 1-12) |
| Arcadia Planitia (Mars), ice access | 169 | Pinacate volcanic field, Mexico (0.754, 2-3); Qaidam Basin, China (0.746, 2-4); Amboy lava field, California (0.735, 3-4) |
| Deuteronilus Mensae (Mars), ice access | 169 | Pinacate volcanic field, Mexico (0.793, 2-3); Qaidam Basin, China (0.780, 3-4); Danakil Depression, Ethiopia (0.771, 4-5) |

- Meteor Crater, Arizona (a known analog) is the top window for both Mars targets (0.815 and 0.837, rank range 1-1).
- No Earth window is a close match for the Moon ridge: the best score is 0.50 and it agrees on only 1 of the 3 measures. The Mars targets are matched much better.

### Backtest against known analogs
Known analogs are ranked among regions and compared with chance by a permutation test (2,000 shuffles).

| Target | Known analogs tested | Median rank of known analogs (expected by chance) | p | Best known analog |
|---|---|---|---|---|
| Moon ridge | 7 | 7 of 17 (9.0) | 0.091 | Mauna Kea, Hawaii (rank 3) |
| Arcadia Planitia | 6 | 9.5 of 16 (8.5) | 0.369 | Meteor Crater, Arizona (rank 1) |
| Deuteronilus Mensae | 6 | 9.5 of 16 (8.5) | 0.356 | Meteor Crater, Arizona (rank 1) |

How to read this: for the Moon, known analogs rank somewhat better than chance, but not significantly. For Mars, only Meteor Crater stands out; the other known analogs rank no better than chance. Treat the backtest as a consistency check, not proof. The percentiles shown in the viewer compare each region's best window with all windows, which is biased upward, so the region-level test above is the fair one.

## Quick start
1. `pip install -r requirements.txt`
2. Open `web/index.html` (or the live link) to explore the committed results.
3. Optional logic check: `python proving_grounds.py selftest --out selftest.json`

**Warning:** `selftest` writes a synthetic demo to its `--out` path. Without `--out` that path is `web/results.json` (and `web/results.js`), which replaces the real results and brings back the synthetic banner. Always pass `--out`, as above. If you overwrite them by accident, restore them with `git checkout -- web/results.json web/results.js`.

## Real run
`python proving_grounds.py run`

It does everything in one go and reuses files already in `data/`:
1. Downloads the Moon target: NASA PGDA's 5 m LOLA model of Site 01, the Connecting Ridge (saved as `data/moon.tif`).
2. Clips two 25 km Mars targets out of the USGS 463 m MOLA model by reading only those windows over HTTP: Arcadia Planitia (47.2 N, 184.3 E) and Deuteronilus Mensae (43.8 N, 29.4 E).
3. Downloads one Copernicus GLO-30 tile per Earth region and the matching ESA WorldCover land-cover window (used to drop windows that are less than half bare or sparse ground).
4. Scores, backtests and writes `web/results.json` and `web/results.js`.

Then open `web/index.html`. To publish, deploy to Vercel (see below).

Other commands (all accept `--out`, default `web/results.json`):
- `python proving_grounds.py fetch`  steps 1-2 only: download or clip the target files into `data/`
- `python proving_grounds.py revalidate`  recompute the backtest from an existing results file in seconds, with no downloads
- `python proving_grounds.py selftest`  synthetic terrain; tests the scoring logic (see the warning above)

`run` also accepts `--moon`, `--mars` and `--mars2` to use your own file instead of `data/<name>.tif`, and `--mars-cyl` to force the simple-cylindrical correction for Mars DEMs.

If the automatic Mars clip fails, clip by hand and save as `data/mars.tif` (and `data/mars2.tif`):
`gdal_translate -projwin 236469 2810109 273264 2785109 /vsicurl/https://planetarymaps.usgs.gov/mosaic/Mars_MGS_MOLA_DEM_mosaic_global_463m.tif data/mars.tif`
(check first that `gdalinfo` shows Upper Left near -10668849, 5334424). In QGIS use Raster > Extraction > Clip Raster by Extent with the same numbers. If Mars slopes look too gentle, add `--mars-cyl`.

If `pip install rasterio` is hard on your laptop, run it in Google Colab (`!pip install rasterio`), then download `web/results.js` and `web/results.json` back into the folder.

## Deploy to Vercel
The viewer in `web/` is a static site, so there is no build step. `vercel.json` points Vercel at `web/`, and `.vercelignore` keeps the Python pipeline and `data/` out of the upload. The live viewer is deployed this way.

- **From GitHub:** push this folder to a repo, then in Vercel choose Add New > Project, import the repo and press Deploy. Leave the settings as they are.
- **From the CLI:** `npm i -g vercel`, then run `vercel` in this folder (`vercel --prod` for the production URL).

The pipeline does not run on Vercel. Run `python proving_grounds.py run` on your machine, then redeploy with the new `web/results.js` and `web/results.json`.

## Viewer features
- Target switcher (Moon ridge, Arcadia Planitia, Deuteronilus Mensae)
- Weight sliders that re-rank everything live
- Rank-stability chart, and a site profile card with "where Earth fails"
- CSV download of the top 10, and a field validation plan per site
- Backtest panel against known analogs
- A "next layers, not used yet" list that links to candidate datasets for future work

## Method
Each target and each 25 km Earth window is described by three terrain measures:
- **Slope distribution:** a histogram of slope in 9 bins from 0 to 90 degrees.
- **Roughness:** the standard deviation of elevation after removing a 3 by 3 local mean.
- **Relief:** the 5th to 95th percentile elevation range.

Earth elevation is block-averaged to the target's pixel size (at least 90 m), so all comparisons share one scale. Each measure gets a 0-1 similarity. Slope uses Earth Mover's distance, `exp(-distance / tau)`, where `tau` is 35% of the target's mean slope in degrees (minimum 0.5). Roughness and relief use min/max ratios. The score is the weighted mean of the three. Task weights are editable in `TARGETS` (slope / roughness / relief: Moon rover mobility 0.5 / 0.3 / 0.2, Mars ice access 0.3 / 0.3 / 0.4). Windows with more than 10% missing elevation, or less than 50% bare or sparse ground, are dropped.

Each window also gets an agreement count (measures at 70% or better) and a rank range: the weights are varied by 20% in 200 random draws, and the 5th to 95th percentile of the window's rank is reported.

The backtest works at region level. It takes the best window of each region, ranks the regions, and compares the median rank of the known analogs with a permutation null (2,000 shuffles of window scores that keep region sizes). The reported p is the share of shuffles that do at least as well as the real ranking.

## Reproducibility
- The weight draws and the permutation test use fixed seeds, so `revalidate` reproduces the backtest from `web/results.json`.
- The self-test uses its own fixed seed on generated terrain and fails if the scoring cannot recover the matching terrain class.
- Results depend on the data versions used: ESA WorldCover 2021 (v200), Copernicus GLO-30, the PGDA LOLA Site 01 model and the USGS MOLA 463 m mosaic.

## Honest limits
- Terrain only (slope, roughness, relief). No thermal, mineral, illumination, cave or accessibility layers yet.
- Land cover comes from ESA WorldCover 2021. Windows with under 50% bare or sparse ground are dropped. Where WorldCover has no tile (for example Antarctica), windows are kept unmasked.
- Mars targets use the 463 m global MOLA model, so Earth is compared with them at that coarser scale.
- Earth cannot match vacuum, low gravity, radiation or extreme cold.
- Weights are assumptions set from the task, not fitted to known sites. Region coordinates in `KNOWN` and `CAND` are approximate starting points; each region is the 1 degree tile that contains the point. Edit them freely.
- The backtest uses 6 or 7 known analogs among desert-like regions only, and its p-values (0.09 to 0.37) are not significant. It shows consistency, not proof.
- Rankings are candidates for field checking, not verified analogs.

## Data sources and credits
- Moon target: NASA PGDA, 5 m LOLA model of Site 01 (Connecting Ridge): `https://pgda.gsfc.nasa.gov/data/LOLA_5mpp/Site01/Site01_final_adj_5mpp_surf.tif`
- Mars targets: USGS Astrogeology, MOLA 463 m global mosaic: `https://planetarymaps.usgs.gov/mosaic/Mars_MGS_MOLA_DEM_mosaic_global_463m.tif`
- Earth elevation: Copernicus GLO-30 DEM (30 m), read from the public `copernicus-dem-30m` bucket
- Land cover: ESA WorldCover 2021 (v200), read from the public `esa-worldcover` bucket
- Viewer: Leaflet 1.9.4, with Esri World Imagery and OpenStreetMap basemaps

Follow each provider's attribution and terms when you reuse their data.
