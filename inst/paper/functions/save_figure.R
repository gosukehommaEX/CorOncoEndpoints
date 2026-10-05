# Save a ggplot2 figure as .eps and .pdf
#
# The .eps file (used in the article) is written with grDevices::postscript()
# and the .pdf file with grDevices::pdf(), both with the Helvetica family.
# These devices use the standard PostScript fonts (Greek letters from plotmath
# expressions come from the Symbol font), so that no Type 3 font appears when
# the .eps file is converted to PDF. cairo_ps() was not used because its .eps
# files gave a Type 3 font after conversion with Ghostscript. All labels must
# therefore be ASCII or plotmath expressions. dpi = 800 is passed for any part
# that would be rasterized (there is none).
#
# Arguments
#   plot   ggplot or patchwork object
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
  ggplot2::ggsave(f_eps, plot, device = "eps", width = width, height = height,
                  units = "in", dpi = 800, family = "Helvetica")
  ggplot2::ggsave(f_pdf, plot, device = grDevices::pdf, width = width, height = height,
                  units = "in", dpi = 800, family = "Helvetica", useDingbats = FALSE)
  invisible(c(f_eps, f_pdf))
}
