## .................................................................................
## Purpose: Run Newman Approximations on Large Simulation Models
##
## Author: Nick Brazeau
##
## Date: 22 June, 2023
##
## Notes:
## .................................................................................
library(fomes)
library(igraph)

#............................................................
# read in my Nets
#...........................................................
snin <- list.files("simulations/00_snakeinput_data/", full.names = T)
snin <- lapply(snin, readRDS) %>%
  dplyr::bind_rows()

#............................................................
# fomes wrapper
#...........................................................
integralget_newman_mean_final_epidemic_size <- function(network, durationI, betaI,
                                                        network_manip, param, val, path) {
  out <- fomes:::get_newman_mean_final_epidemic_size(graph = network,
                                                     taui = durationI,
                                                     rij = betaI,
                                             transmApprox = FALSE,
                                             initu = 0.5, iters = 1e6, tol = 1e-5)
  return(out)
}

#............................................................
# Run Newman Approx
#...........................................................
sninout <- snin %>%
  dplyr::mutate(NewmanApprox = purrr::pmap_dbl(., integralget_newman_mean_final_epidemic_size))

#............................................................
# out
#...........................................................
sninout %>%
  dplyr::mutate(dplyr::mutate(basenet = stringr::str_extract(basename(unique(ret$netpath)), "(?<=_)[0-9]+(?=\\.RDS$)"))) %>%
  dplyr::select(-c("param", "path", "network")) %>%
  saveRDS(., file = "analyses/zbackend-large_simulation_sirpent_runs/newman_approximation_from_inputs.RDS")
