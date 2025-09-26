## Quick purpose
This repository contains R scripts used to detect and classify eddy subduction events from Argo float data (physical and carbon detection). The goal of these instructions is to give an AI code assistant the contextual knowledge needed to make safe, useful edits and run the code locally.

## High-level architecture (what to read first)
- scripts/ contains the runnable detection workflows: `detection_algorithm_physical_subduction.R` and `detection_algorithm_carbon_subduction.R`.
- Each detection script sources a corresponding annex file in the same folder (`*_annex_fun.R`) which exposes the small function library (e.g. `downscale_data_fun_wo_out`, `downscale_data_fun`, `perform_checks`, `compute_derivatives_vectorized`). Read the annex files to understand core data transforms.
- data/ holds input RDS/CSV used by the scripts (e.g. `WMO_list.rds`, `manual_classification*.csv`). Scripts expect specific filenames — check the scripts before renaming or moving files.

## Key data flow and patterns
- Load: scripts call an external `load_float_data()` to fetch a single-float dataframe. That function is NOT in this repo — it comes from the OneArgo toolbox or the runtime R environment. Ensure it is available before running the scripts.
- Preprocess: each script computes ABS_SAL, CONS_TEMP, SAT_DOXY, AOU, SIGMA0 using `gsw` and other packages.
- Resampling: the annex functions downscale profiles into regular pressure bins (20 m and 40 m) via `downscale_data_fun_wo_out()` and `downscale_data_fun()` and then compute robust residuals (rolling median, trimmed mean TM_9).
- Detection: residuals are scaled (IQR-normalized), thresholded (default cutoff 1.96) and candidate events are validated with `perform_checks()` which uses first and second derivative tests.

## Where to run / working-directory gotchas (important)
- The scripts use relative paths and assume they are executed with working directory set to `scripts/`. Examples that demonstrate this pattern:
  - `source("detection_algorithm_physical_subduction_annex_fun.R")` (no path — annex file expected in current dir)
  - `wmolist <- readRDS("../data/WMO_list.rds")` (reads from repo-root/data when executed from `scripts/`)
  - Some write operations use `write_csv("data/detected_physical_subd_events.csv")` (note the missing `..`) — verify the intended output path before running.

## How to run (examples)
Run from the repository root by changing into the scripts directory first (recommended):

```bash
cd scripts
Rscript detection_algorithm_physical_subduction.R
# or
Rscript detection_algorithm_carbon_subduction.R
```

If you run the scripts from a different working directory, update the relative paths (the scripts rely on `../data/` for inputs and local annex files being in the current folder).

## Important functions & small examples (copyable for edits)
- downscale_data_fun_wo_out(df, bin_width = 40)
  - Groups by pressure bins, averages variables, and returns a downsampled profile with PRES_ADJUSTED set to bin centers.
- perform_checks(profile_data, target_level, variable_name, second_deriv, window)
  - Uses `compute_derivatives_vectorized()` to compute d/dP and d2/dP2 and returns boolean checks used to validate candidate anomalies.

## Local mock loader & runner (for broken OneArgo)
- If OneArgo's `load_float_data()` is unavailable, a lightweight mock and runner were added for local testing:
  - `scripts/mock_load_float_data.R` - defines `load_float_data()` that returns synthetic profile data.
  - `scripts/run_detection_with_mock.R` - runs a minimal detection pass using the mock loader and the existing annex functions.
  - Run locally from repository root:

```bash
Rscript scripts/run_detection_with_mock.R
```

Note: this environment does not include R so we could not execute the runner here — run the command on your machine with R and the packages listed in the Dependencies section.

## Dependencies and environment
- The scripts assume an R environment with at least: tidyverse, robustbase, gsw, zoo, oce, ggpubr, segmented, conflicted, pracma. The `load_float_data()` helper is external (OneArgo or equivalent).
- If a package is missing, use `install.packages()` or your environment manager (renv/packrat) — this repo does not include an environment lockfile.

## Common pitfalls and things to verify before editing or running
- Filenames in `data/` may differ from those referenced in scripts (e.g. `manual_classification_raw.csv` exists in the repo but scripts refer to `../data/manual_classification.csv`). Confirm exact filenames and update scripts or data accordingly.
- Relative paths are fragile: prefer running `cd scripts` first or modify scripts to build absolute paths programmatically (e.g. via here::here()).
- `load_float_data()` must be available in the R session (not present in this repo). If you cannot access OneArgo, mock this function during tests with a small dataframe matching expected columns (PRES_ADJUSTED, PSAL_ADJUSTED, TEMP_ADJUSTED, DOXY, etc.).

## Editing guidelines for AI agents
- Small, local changes only: follow the project pattern of putting helpers in `*_annex_fun.R` and keeping main workflows in `detection_algorithm_*.R`.
- When changing IO paths, update both read and write locations consistently and add a short comment explaining why you changed working directory assumptions.
- Preserve numeric hyperparameters unless a user asks to tune them; they appear at the top of the scripts (e.g. `cutoff <- 1.96`, `window <- 60`, `resolution <- 40`).

## Files to open first when debugging
- `scripts/detection_algorithm_physical_subduction.R`
- `scripts/detection_algorithm_carbon_subduction.R`
- `scripts/detection_algorithm_physical_subduction_annex_fun.R`
- `scripts/detection_algorithm_carbon_subduction_annex_fun.R`
- `data/WMO_list.rds`, `data/README.md`, and any `data/manual_classification*.csv` for input expectations

If anything above is unclear, tell me which part you want expanded (paths, function signatures, example data shapes) and I will iterate.
