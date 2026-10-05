# Create all figures and tables of the article and combine the recorded numbers
#
# Run from the package root after the data generation scripts and the
# validation script:
#   source("inst/paper/figures_and_tables/make_all.R")
# Writes inst/paper/output/figures, inst/paper/output/tables,
# inst/paper/output/numbers/*.csv and inst/paper/output/numbers.csv.

scripts <- c("figure1_expexp_hazards.R", "figure2_os_distribution.R",
             "figure3_correlations.R", "figure4_two_in_one_type1.R",
             "table1_calibration.R", "table2_required_events.R", "table3_win_statistics.R",
             "tableS1_inputs.R", "tableS2_reproduction.R", "tableS3_generator_validation.R",
             "tableS4_function_validation.R", "tableS5_computing_time.R",
             "tableS6_two_in_one.R", "tableS7_win_statistics.R",
             "figureS1_survival_by_response.R", "figureS2_death_proportion.R",
             "figureS3_two_in_one_alternative.R", "figureS4_two_in_one_g1.R")
for (s in scripts) {
  message("Running ", s)
  source(file.path("inst", "paper", "figures_and_tables", s), local = new.env())
}
num_dir <- file.path("inst", "paper", "output", "numbers")
num <- do.call(rbind, lapply(list.files(num_dir, pattern = "\\.csv$", full.names = TRUE),
                             utils::read.csv, stringsAsFactors = FALSE))
utils::write.csv(num, file.path("inst", "paper", "output", "numbers.csv"), row.names = FALSE)
message("numbers.csv: ", nrow(num), " numbers from ", length(unique(num$item)), " figures and tables")
