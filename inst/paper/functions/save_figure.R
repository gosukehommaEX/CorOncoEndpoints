# Save a ggplot2 figure as .eps and .pdf
#
# The .eps file (used in the article) is written with cairo_ps() and the .pdf
# file with cairo_pdf(); both at 800 dpi, the resolution used for any part
# that has to be rasterized.
#
# Arguments
#   plot   ggplot object
#   name   file name without extension
#   width, height size in inches
#   dir    output directory
# Value
#   the two file paths (invisibly)
save_figure <- function(plot, name, width, height,
                        dir = file.path("inst", "paper", "output", "figures")) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  f_eps <- file.path(dir, paste0(name, ".eps"))
  f_pdf <- file.path(dir, paste0(name, ".pdf"))
  ggplot2::ggsave(f_eps, plot, device = grDevices::cairo_ps, width = width,
                  height = height, units = "in", dpi = 800,
                  fallback_resolution = 800)
  ggplot2::ggsave(f_pdf, plot, device = grDevices::cairo_pdf, width = width,
                  height = height, units = "in", dpi = 800,
                  fallback_resolution = 800)
  invisible(c(f_eps, f_pdf))
}
