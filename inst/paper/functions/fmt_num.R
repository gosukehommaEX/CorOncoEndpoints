# Format numbers with a fixed number of decimals for tables
#
# NA is shown as "--". A minus sign is written in math mode.
#
# Arguments
#   x      numeric vector
#   digits number of decimals
fmt_num <- function(x, digits) {
  out <- formatC(x, format = "f", digits = digits)
  out[is.na(x)] <- "--"
  neg <- !is.na(x) & substr(out, 1, 1) == "-"
  out[neg] <- paste0("$-$", substring(out[neg], 2))
  out
}
