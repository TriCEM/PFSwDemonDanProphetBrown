## .................................................................................
## Purpose: Backend script for looking at Mass Action Divergence
##
## Author: Nick Brazeau
##
## Date: 22 June, 2023
##
## Notes:
## .................................................................................
library(tidyverse)
library(fomes)
#............................................................
# Functions
#...........................................................
# mass action wrapper
MAsim_Gillespie_SIR_wrapper <- function(reps, N, beta, dur_I) {
  ret <- fomes:::sim_Gillespie_SIR(
    Iseed = 1, N = N,
    beta = beta * N,
    dur_I = dur_I,
    term_time = Inf)
  fs <- N - ret[nrow(ret), "Susc"]
  return(fs)
}

# vary degree distribution network wrapper
netddsim_Gillespie_nSIR <- function(reps, N, beta, dur_I, conmat, dd) {
  ret <- fomes::sim_Gillespie_nSIR(Iseed = 1,
                            N = N,
                            beta = rep(beta, N),
                            dur_I = dur_I,
                            rho = 1e-27,
                            init_contact_mat = conmat,
                            term_time = Inf)
  out <- summary(ret)$FinalEpidemicSize
  return(out)
}

# vary rho network wrapper
netrhosim_Gillespie_nSIR <- function(reps, N, beta, dur_I, conmat, rho) {
  ret <- fomes::sim_Gillespie_nSIR(Iseed = 1,
                                   N = N,
                                   beta = rep(beta, N),
                                   dur_I = dur_I,
                                   rho = rho,
                                   init_contact_mat = conmat,
                                   term_time = Inf)
  out <- summary(ret)$FinalEpidemicSize
  return(out)
}




#............................................................
# Magic Numbers
#...........................................................
N <- 100
dur_I <- 5
Iseed <- 1
betaind <- 0.0035 # see Keeling & Grenfell 2000 PMID: 10677276 for justification of bN in MA
reps <- 100

#............................................................
# Mass Action
#...........................................................
# will reuse
matab <- tidyr::expand_grid(reps = 1:reps, N = N, beta = betaind, dur_I = dur_I)
matab <- matab %>%
  dplyr::mutate(fs = purrr::pmap_dbl(., MAsim_Gillespie_SIR_wrapper, .progress = T))

#............................................................
# Network for Degree Distribution
#...........................................................
#......................
# make contact networks
#......................
probofcon <- seq(0.6, 1, 0.1)
edges <- floor(N*probofcon)
edges[length(edges)] <- 99 # one less
conmats <- tibble::tibble(dd = edges) %>%
  dplyr::mutate(conmat = purrr::map(dd, function(x, N){igraph::as_adjacency_matrix( igraph::degree.sequence.game(out.deg = rep(x, N), method = "vl") )}, N = N)
                )
nettabs <- tidyr::expand_grid(reps = 1:reps, N = N, beta = betaind, dur_I = dur_I, dd = edges)
nettabs <- dplyr::left_join(nettabs, conmats, by = "dd")
nettabs <- nettabs %>%
  dplyr::mutate(fs = purrr::pmap_dbl(., netddsim_Gillespie_nSIR, .progress = T))


#............................................................
# Network for Varying Rho
#...........................................................
#......................
# make contact networks
#......................
rhovary <- sapply(c(5,3,0,-1,-2,-3), function(x){10^(-x)})
conmat <- igraph::as_adjacency_matrix( igraph::degree.sequence.game(out.deg = rep(20, N), method = "vl") )
rhotabs <- tidyr::expand_grid(reps = 1:reps, N = N, beta = betaind, dur_I = dur_I,
                              conmat = list(conmat), rho = rhovary)
rhotabs <- rhotabs %>%
  dplyr::mutate(fs = purrr::pmap_dbl(., netrhosim_Gillespie_nSIR, .progress = T))


#............................................................
# out
#...........................................................
out <- list(nettabs = nettabs, matab = matab, rhotabs = rhotabs)
saveRDS(out, paste0(here::here(),
                    "/analyses/zbackend-low_hang_fruit/mass_action_divergence.RDS"))
