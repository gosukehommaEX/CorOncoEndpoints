# Seed of a batch of simulated trials
#
# seed = 10^6 * script + 10^3 * scenario + batch, so that every batch of every
# scenario of every data generation script has its own seed. A second call of
# rOncoEndpoints() for the same batch (the expansion part of the 2-in-1 design)
# adds 500000.
#
# Arguments
#   script   number of the data generation script (1 to 9)
#   scenario scenario number (1 to 499)
#   batch    batch number (1 to 999)
seed_of <- function(script, scenario, batch) {
  stopifnot(script >= 1, script <= 9, scenario >= 1, scenario <= 499,
            batch >= 1, batch <= 999)
  as.integer(1e6 * script + 1e3 * scenario + batch)
}
