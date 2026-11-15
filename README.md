# Proving Grounds

Finds places on Earth whose terrain matches NASA's Moon Base target (Shackleton Connecting Ridge) and a Mars ice-access target (Arcadia Planitia), then backtests the method against known analog sites. Built for NASA Space Apps 2026: *Identify Earth Locations that Analog the Permanent Moon Base Locations and Mars*.

## What is in this folder
- `proving_grounds.py`  the pipeline: downloads data, scores, backtests, exports
- `web/index.html`  the map viewer (Leaflet; basemap is NASA Blue Marble, switchable in the layer control)
- `web/results.js`, `web/results.json`  **SYNTHETIC test data** so the viewer opens straight away. A warning banner shows until you replace it with real results.
- `data/`  created on the first real run (target elevation files)

## Quick test (2 minutes)
1. `pip install -r requirements.txt`
2. `python proving_grounds.py selftest` checks the scoring on generated terrain for all three targets.
3. Open `web/index.html`.

## Real run
`python proving_grounds.py run`

It does everything in one go and reuses files already in `data/`:
1. Downloads the Moon target: NASA PGDA's 5 m LOLA model of Site 01, the Connecting Ridge (saved as `data/moon.tif`).
2. Clips two 25 km Mars targets out of the USGS 463 m MOLA model by reading only those windows over HTTP: Arcadia Planitia (47.2 N, 184.3 E) and Deuteronilus Mensae (43.8 N, 29.4 E).
3. Downloads one Copernicus GLO-30 tile per Earth region and the matching ESA WorldCover land-cover window (used to drop windows that are less than half bare or sparse ground).
4. Scores, backtests and writes `web/results.json` and `web/results.js`.

Then open `web/index.html`. The banner disappears when the data is real. To publish, deploy to Vercel (see below).

If the automatic Mars clip fails, clip by hand and save as `data/mars.tif` (and `data/mars2.tif`):
`gdal_translate -projwin 236469 2810109 273264 2785109 /vsicurl/https://planetarymaps.usgs.gov/mosaic/Mars_MGS_MOLA_DEM_mosaic_global_463m.tif data/mars.tif`
(check first that `gdalinfo` shows Upper Left near -10668849, 5334424). In QGIS use Raster > Extraction > Clip Raster by Extent with the same numbers. If Mars slopes look too gentle, add `--mars-cyl`.

If `pip install rasterio` is hard on your laptop, run it in Google Colab (`!pip install rasterio`), then download `web/results.js` and `web/results.json` back into the folder.

## Deploy to Vercel
The viewer in `web/` is a static site, so there is no build step. `vercel.json` points Vercel at `web/`, and `.vercelignore` keeps the Python pipeline and `data/` out of the upload.

- **From GitHub:** push this folder to a repo, then in Vercel choose Add New > Project, import the repo and press Deploy. Leave the settings as they are.
- **From the CLI:** `npm i -g vercel`, then run `vercel` in this folder (`vercel --prod` for the production URL).

The pipeline does not run on Vercel. Run `python proving_grounds.py run` on your machine, then redeploy with the new `web/results.js` and `web/results.json`.

## Viewer features
- Target switcher (Moon ridge, Arcadia Planitia, Deuteronilus Mensae)
- Weight sliders that re-rank everything live
- Rank-stability chart, site profile card with "where Earth fails"
- CSV download of the top 10, and a field validation plan per site
- Backtest panel against known analogs

## Method
Each target and each 25 km Earth window is described by three terrain measures: slope distribution, roughness and relief. Earth elevation is block-averaged to the target's pixel size so all comparisons share one scale. Similarity uses Earth Mover's distance for slope and ratios for roughness and relief. The score combines them with task weights (editable in `TARGETS`). Each window also gets an agreement count (criteria at 70% or better) and a rank range when weights vary by 20%. The backtest asks where the best window of each known analog ranks among all windows.

## Honest limits
- Terrain only. No thermal, mineral, illumination, cave or accessibility layers yet.
- Land cover comes from ESA WorldCover 2021. Where it has no tile (for example Antarctica), windows are kept unmasked.
- Weights are assumptions, not fitted to known sites. Region coordinates in `KNOWN` and `CAND` are approximate starting points; each region is the 1 degree tile that contains the point. Edit them freely.
- Rankings are candidates for field checking, not verified analogs.
- The targets use NASA LOLA and MOLA data. List them, and any AI tools you used, on the submission form.
