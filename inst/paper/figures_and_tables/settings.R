# Common settings of the figure and table scripts for the article
#
# Sourced at the top of each script in inst/paper/figures_and_tables. The
# scripts only read the files in inst/paper/data (and the outputs of
# inst/reproduce, inst/validation and inst/paper/validation); they do not
# simulate. Run from the package root.
#
# Statistics in Medicine: figures as .eps (800 dpi for graphs), no tints
# (shaded areas), keys inside the artwork, compound figures as one file.

for (p in c("ggplot2", "patchwork")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop("The ", p, " package is required for the figures.", call. = FALSE)
  }
}
if (utils::packageVersion("ggplot2") < "3.5.0") {
  stop("ggplot2 >= 3.5.0 is required (legend.position.inside).", call. = FALSE)
}
if (!exists("OncoArm")) library(CorOncoEndpoints)
source(file.path("inst", "paper", "load_functions.R"))
library(ggplot2)
library(patchwork)
data_dir <- file.path("inst", "paper", "data")
fig_width <- 6.5      # inches; the manuscript scales figures to \textwidth
z_975 <- stats::qnorm(0.975)

# Okabe-Ito colours (distinguishable with colour vision deficiency), used
# together with line types so that the figures remain readable in grey scale
pal <- c("#000000", "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")
theme_paper <- theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", colour = "black"),
        legend.key = element_rect(fill = "white", colour = NA),
        legend.background = element_rect(fill = "white", colour = NA),
        plot.tag = element_text(face = "bold"))
