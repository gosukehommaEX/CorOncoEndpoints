# Write a complete LaTeX table environment
#
# The file contains \begin{table} ... \end{table} with the caption, label,
# a booktabs tabular and the notes (threeparttable and tablenotes), so that
# the manuscript only needs \input{}. Requires the LaTeX packages booktabs and
# threeparttable.
#
# Arguments
#   name    file name without extension
#   caption caption text (LaTeX)
#   label   LaTeX label
#   align   column specification of the tabular, e.g. "llrr"
#   header  character vector of header rows (each ending with \\)
#   body    character vector of body rows (each ending with \\, or \midrule)
#   notes   character vector of notes (one \item each), or NULL
#   size    font size command such as "\small", or NULL
#   dir     output directory
# Value
#   the file path (invisibly)
write_table_tex <- function(name, caption, label, align, header, body,
                            notes = NULL, size = NULL,
                            dir = file.path("inst", "paper", "output", "tables")) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  lines <- c("\\begin{table}[htbp]",
             "\\centering",
             paste0("\\caption{", caption, "}"),
             paste0("\\label{", label, "}"),
             "\\begin{threeparttable}",
             if (!is.null(size)) size,
             paste0("\\begin{tabular}{", align, "}"),
             "\\toprule",
             header,
             "\\midrule",
             body,
             "\\bottomrule",
             "\\end{tabular}",
             if (!is.null(notes)) c("\\begin{tablenotes}", "\\footnotesize",
                                    paste0("\\item ", notes), "\\end{tablenotes}"),
             "\\end{threeparttable}",
             "\\end{table}")
  f <- file.path(dir, paste0(name, ".tex"))
  writeLines(lines, f, useBytes = TRUE)
  invisible(f)
}
