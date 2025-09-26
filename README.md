# globargo

Repository supporting the manuscript on biogeochemical Argo (BGC-Argo) detection of
the eddy subduction pump. It contains the detection algorithms, figure generation
code, notebooks, and processed data products used in the analysis.

## Repository layout

| Path | Purpose |
| --- | --- |
| `scripts/` | Production detection algorithms for physical and carbon subduction (plus supporting functions). |
| `results/` | Reproducible scripts that regenerate the manuscript-quality figures once data products are in place. |
| `data/` | Processed inputs and intermediate outputs (see `data/README.md` for sourcing instructions). |
| `pubfig/` | Exported figures ready for publication. |
| `notebooks/` | Jupyter notebooks used for exploratory plots (e.g., histogram & ESP calculations). |
| `OLD_script/` | Historical prototypes and helper utilities kept for provenance. |
| `renv/`, `renv.lock` | Reproducible R environment captured with `renv`. |
| `globargo_repo.Rproj` | RStudio project file for interactive development. |

## Prerequisites

- **R**: ≥ 4.3 (tested with 4.3.3 on Ubuntu 24.04).
- **System libraries**: installing `sf`, `terra`, and `rnaturalearth` typically requires GDAL/GEOS/PROJ development headers; consult your OS package manager if `renv::restore()` reports missing system dependencies.
- **R packages** (automatically handled by `renv`): `tidyverse`, `gsw`, `oce`, `robustbase`, `zoo`, `ggpubr`, `segmented`, `pracma`, `mgcv`, `gratia`, `sf`, `ggspatial`, `rnaturalearth`, `patchwork`, `viridis`, `conflicted`, and others captured in `renv.lock`.
- **OneArgo toolbox**: the detection scripts rely on `oneArgo::load_float_data()` to fetch BGC-Argo profiles. Install it from <https://github.com/ArgoCanada/oneArgo> and configure the data cache (e.g., set `options(oneArgo.data = "path/to/cache")`).
- **Python (optional)**: the notebook in `notebooks/` expects a recent Python 3 environment with `pandas`, `xarray`, `netCDF4`, and plotting libraries if you plan to re-run it.

## Environment management with `renv`

1. Install the `renv` bootstrapper in your user library if necessary:
   ```r
   install.packages("renv")
   ```
2. From the repository root, restore the project library to match the versions recorded in `renv.lock`:
   ```r
   renv::restore()
   ```
   This step installs all R dependencies into a project-local library under `renv/library/`.
3. When adding new packages or updating existing ones, run `renv::snapshot()` to refresh `renv.lock` before committing.

## Running the detection algorithms

### 1. Prepare supporting inputs

- **Manual classifications**: export the latest manual review spreadsheets (e.g., `manual_classification.csv`, `carbon_cat1_class.csv`, `carbon_cat2_class.csv`) into `data/`. These identify which WMO/cycle combinations should be processed.
- **Float list (`wmolist`)**: create an RDS file containing the numeric vector of float WMO identifiers you want to analyse and place it under `data/wmolist.rds` (or update the script to point at your preferred location).
- **OneArgo configuration**: confirm that `oneArgo::load_float_data()` can download or access the required profiles before running the algorithms.

### 2. Physical subduction detection

1. Edit `scripts/detection_algorithm_physical_subduction.R` if needed to point `wmolist <- readRDS("fill_in")` at your `data/wmolist.rds`.
2. Optional: adjust the hyper-parameters (`cutoff`, `resolution`, `window`) near the top of the script.
3. Execute the algorithm from the repository root:
   ```bash
   Rscript scripts/detection_algorithm_physical_subduction.R
   ```
4. The script downloads the required float profiles via `oneArgo`, computes anomalies, and saves the consolidated detections to `data/detected_physical_subd_events.csv`.

### 3. Carbon subduction detection

1. Place the manual classification CSV mentioned above in `data/manual_classification.csv` (the script expects at least `WMO`, `CYCLE_NUMBER`, and `Category` columns).
2. Edit `scripts/detection_algorithm_carbon_subduction.R` so that the `classified_df <- read.csv(...)` and `wmolist` inputs point at your files and adjust output paths if you want them to land inside this repository (the historical absolute paths under `/data/GLOBARGO/` can be replaced with `data/` paths).
3. Run the carbon detector (it loops through category 1 and category 2 events):
   ```bash
   Rscript scripts/detection_algorithm_carbon_subduction.R
   ```
4. The script augments the manual classifications with recalculated anomalies/POC and writes per-category outputs (update the `write_csv()` paths to land in `data/` if required).

Both detection scripts are long-running: they fetch large NetCDF files and perform robust statistics per cycle. Run them on a machine with adequate RAM/disk and expect several hours for the full float list.

## Regenerating figures

Once the detection outputs and derived data sets listed in `data/README.md` are available:

1. Ensure the renv environment is active (`renv::restore()` completed successfully).
2. From the repository root, run:
   ```bash
   Rscript results/script_figures_generation.R
   ```
   The script sets the working directory to `scripts/` to reuse helper functions, loads the processed data from `data/`, and saves figure panels into `pubfig/`.
3. Optional: open `notebooks/04HistogramFigureAndESP.ipynb` to regenerate the ESP histograms after updating `data/histogram_data_full_for_python.csv`.

Refer to `data/README.md` for the provenance of each intermediate required by the figure script.

## Additional notes

- The `OLD_script/` folder documents exploratory analyses and is useful when reproducing historical data products (e.g., GAM surfaces, sensitivity experiments). They often contain hard-coded absolute paths—update them to point to your local clone before re-running.
- Large upstream data products (e.g., OCIM-2 sequestration grids, manual classification spreadsheets) are not tracked in Git. Coordinate with the manuscript authors to obtain the raw inputs.
