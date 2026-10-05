# Information stored with every data file of the article
#
# Arguments
#   t_start value of proc.time() at the start of the script
#   pilot   TRUE for a pilot run
# Value
#   list with the creation time, elapsed seconds, pilot flag, R and package
#   versions and platform
run_info <- function(t_start, pilot) {
  list(created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
       elapsed_sec = unname((proc.time() - t_start)[["elapsed"]]),
       pilot = isTRUE(pilot), r_version = R.version.string,
       package_version = as.character(utils::packageVersion("CorOncoEndpoints")),
       platform = R.version$platform)
}
