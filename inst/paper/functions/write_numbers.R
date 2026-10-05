# Record the numbers shown in a figure or table
#
# Writes output/numbers/<item>.csv with one row per number. The future script
# that checks the numbers in the manuscript against the R results reads these
# files (combined into output/numbers.csv by make_all.R).
#
# Arguments
#   item        figure or table name, e.g. "table1_calibration"
#   key         unique key of each number within the item
#   value       numeric value (unrounded)
#   description short description of each number
#   source      data file the number comes from
#   dir         output directory
# Value
#   the data frame written (invisibly)
write_numbers <- function(item, key, value, description, source,
                          dir = file.path("inst", "paper", "output", "numbers")) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  df <- data.frame(item = item, key = key, value = as.numeric(value),
                   description = description, source = source,
                   stringsAsFactors = FALSE)
  if (anyDuplicated(df$key)) stop("Duplicated keys in ", item, call. = FALSE)
  utils::write.csv(df, file.path(dir, paste0(item, ".csv")), row.names = FALSE)
  invisible(df)
}
