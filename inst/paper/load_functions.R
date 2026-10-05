# Load the functions used by the scripts in inst/paper (run from the package
# root, after devtools::load_all() or library(CorOncoEndpoints))

for (f in list.files(file.path("inst", "paper", "functions"), pattern = "\\.R$",
                     full.names = TRUE)) {
  source(f)
}
rm(f)
