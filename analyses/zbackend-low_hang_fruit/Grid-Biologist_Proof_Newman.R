## .................................................................................
## Purpose: Biologist Proof Large
##
## Author: Nick Brazeau
##
## Date: 12 June, 2023
##
## Notes:
## .................................................................................
source("R/Newman_Eqs.R")
library(tidyverse)
library(igraph)
library(fomes)

#............................................................
# wrapper functions
#...........................................................
fomes_ct_wrapper <- function(N, rij, taui, conmat) {
  N <- 1e3
  iters <- 100
  final_size <- rep(NA, iters)
  for (i in 1:iters) {
    sim <- fomes::sim_Gillespie_nSIR(Iseed = 1,
                                     N = N,
                                     beta = rep(rij, N),
                                     dur_I = taui,
                                     rho = 1e-27,
                                     init_contact_mat = conmat,
                                     term_time = Inf)
    final_size[i] = sum(sim$Event_traj == "transmission") + 1 # for init
  }
  return(final_size)
}


fomes_discrete_wrapper <- function(N, rij, taui, conmat) {
  N <- 1e3
  iters <- 100
  dicfinal_size <- rep(NA, iters)
  for (i in 1:iters) {
    dicsim <- fomes:::sim_DTDC_nSIR(Iseed = 1,
                                    N = N,
                                    beta = rep(rij, N),
                                    dur_I = taui,
                                    init_contact_mat = conmat,
                                    time_steps = Inf)
    dicfinal_size[i] = N - dicsim$Susc[length(dicsim$Susc)]
  }
  return(dicfinal_size)
}

approxget_newman_mean_final_epidemic_size <- function(graph, taui, rij) {
  out <- get_newman_mean_final_epidemic_size(graph, taui, rij,
                                             transmApprox = TRUE,
                                             initu = 0.5, iters = 1e6, tol = 1e-5)
  return(out)
}

integralget_newman_mean_final_epidemic_size <- function(graph, taui, rij) {
  out <- get_newman_mean_final_epidemic_size(graph, taui, rij,
                                             transmApprox = FALSE,
                                             initu = 0.5, iters = 1e6, tol = 1e-5)
  return(out)
}

#............................................................
# magic numbers
#...........................................................
set.seed(42)
# constants
N <- 1e3
pconn <- 5 * N

g <- igraph::erdos.renyi.game(n = N, p.or.m = pconn,
                              type = "gnm", directed = F) # using same network
hist(degree(g))
conmat <- igraph::as_adjacency_matrix(g, sparse = T)
# varying params
rij <- seq(0.01, 0.1, by = 0.0025)
taui <- c(2:10)

#............................................................
# grid run
#...........................................................
pgrid <- tidyr::expand_grid(rij, taui) %>%
  dplyr::mutate(graph = list(g),
                conmat = list(conmat))


#......................
# run newman
#......................
pgrid$newmanApprox <- purrr::pmap_dbl(pgrid[,c("graph", "rij", "taui")],
                                      approxget_newman_mean_final_epidemic_size)
pgrid$newmanIntegral <- purrr::pmap_dbl(pgrid[,c("graph", "rij", "taui")],
                                        integralget_newman_mean_final_epidemic_size)
hist(pgrid$newmanIntegral)
# let's subset to more interesting ones
pgrid <- pgrid %>%
  dplyr::filter(newmanIntegral < 0.95 & newmanIntegral > 0.05)


#......................
# run fomes
#......................
pgrid$fomesCT <- purrr::pmap(pgrid[,c("conmat", "rij", "taui")],
                             fomes_ct_wrapper, .progress = T)
pgrid$fomesDT <- purrr::pmap(pgrid[,c("conmat", "rij", "taui")],
                             fomes_discrete_wrapper, .progress = T)

#............................................................
# out
#...........................................................
saveRDS(pgrid, "analyses/zbackend-low_hang_fruit/NewmanBioProof.RDS")

