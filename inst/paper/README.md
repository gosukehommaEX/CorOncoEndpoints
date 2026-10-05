# Figures and tables of the article

This folder contains the code that creates every figure and table of the
article on CorOncoEndpoints (main text and Supporting Information). It is
excluded from the package build (`.Rbuildignore`).

## Folders

| Folder | Content |
|---|---|
| `functions/` | Functions used by the scripts (one function per file), loaded by `load_functions.R` |
| `data_generation/` | Scripts that compute or simulate the results and save them in `data/` |
| `data/` | Results (`*.rds`); pilot runs are saved in `data/pilot/` |
| `figures_and_tables/` | Scripts that read `data/` only and write `output/` |
| `output/figures/` | Figures (`.eps` for the article, `.pdf` for viewing) |
| `output/tables/` | Complete LaTeX table environments (`\input{}` in the manuscript; need booktabs and threeparttable) |
| `output/numbers/`, `output/numbers.csv` | Every number shown in the figures and tables, with its source file |
| `validation/` | Validation of the functions against survival, stats, BuyseTest and independent Python code |

## Required packages

survival, BuyseTest (validation); ggplot2 (>= 3.5.0) and patchwork (figures).
These are not listed in DESCRIPTION.

## Order of execution

Open `CorOncoEndpoints.Rproj` (working directory = package root), then run:

```r
devtools::load_all()
source("inst/paper/validation/validate_example_functions.R")
# optional pilot runs (1% of the trials; prints the projected run time)
pilot <- TRUE
source("inst/paper/data_generation/02_design_simulation.R")
source("inst/paper/data_generation/03_two_in_one_simulation.R")
source("inst/paper/data_generation/04_win_statistics_simulation.R")
rm(pilot)
# full runs
source("inst/paper/data_generation/01_model_quantities.R")
source("inst/paper/data_generation/02_design_simulation.R")
source("inst/paper/data_generation/03_two_in_one_simulation.R")
source("inst/paper/data_generation/04_win_statistics_simulation.R")
source("inst/paper/data_generation/05_computing_time.R")
source("inst/paper/figures_and_tables/make_all.R")
```

The tables S2 and S3 read the outputs of `inst/reproduce/reproduce_published.R`
and `inst/validation/validate_generator.R`, which are not run again.

## Figures and tables

| Item | Script | Data |
|---|---|---|
| Figure 1 | `figure1_expexp_hazards.R` | `model_quantities.rds` |
| Figure 2 | `figure2_os_distribution.R` | `model_quantities.rds` |
| Figure 3 | `figure3_correlations.R` | `model_quantities.rds` |
| Figure 4 | `figure4_two_in_one_type1.R` | `two_in_one_simulation.rds` |
| Table 1 | `table1_calibration.R` | `model_quantities.rds` |
| Table 2 | `table2_required_events.R` | `design_simulation.rds` |
| Table 3 | `table3_win_statistics.R` | `win_statistics_simulation.rds` |
| Table S1 | `tableS1_inputs.R` | none |
| Table S2 | `tableS2_reproduction.R` | `inst/reproduce/output/reproduce_results.csv` |
| Table S3a-c | `tableS3_generator_validation.R` | `inst/validation/output/validate_generator.csv` |
| Table S4 | `tableS4_function_validation.R` | `validation/output/validate_example_functions.csv` |
| Table S5 | `tableS5_computing_time.R` | `computing_time.rds` |
| Table S6 | `tableS6_two_in_one.R` | `two_in_one_simulation.rds` |
| Table S7 | `tableS7_win_statistics.R` | `win_statistics_simulation.rds` |
| Figure S1 | `figureS1_survival_by_response.R` | `model_quantities.rds` |
| Figure S2 | `figureS2_death_proportion.R` | `model_quantities.rds` |
| Figure S3 | `figureS3_two_in_one_alternative.R` | `two_in_one_simulation.rds` |
| Figure S4 | `figureS4_two_in_one_g1.R` | `two_in_one_simulation.rds` |

## Seeds

`rOncoEndpoints()` uses the dqrng generator; `set.seed()` has no effect. Each
batch of 1,000 trials uses the seed 10^6 x (script number) + 10^3 x (scenario
number) + (batch number) (`functions/seed_of.R`); the expansion part of the
2-in-1 trials adds 500,000.
