## .................................................................................
## Purpose: Generate `fomes` stochastic runs for analysis
##
## Author: Nick Brazeau
##
## Date: 11 March, 2023
##
## Notes: Using the outline of PFS Manuscript
## .................................................................................
library(tibble)
library(igraph)

#++++++++++++++++++++++++++++++++++++++++++
### Source Functions and Immutables        ####
#++++++++++++++++++++++++++++++++++++++++++
source("R/network_manipulations.R")
set.seed(42) # meaning of life
# population size
N <- 1e3

#++++++++++++++++++++++++++++++++++++++++++
### Maestro Setup                      ####
#++++++++++++++++++++++++++++++++++++++++++
# Based on discussion want to have the following structure:
#  * 5 beta values
#  * 5 duration of illness
#  * 10 "base" ER networks that will be adjusted by network mechanisms (default is p = 0.3)
#    * Network mechanisms:
#      * 10 * 5 degree distributions
#      * 10 modularity models
#      * 10 clustering
#      * 10 neighbor exchange rate

#......................
# STEP 0: generate beta and duration of illness
# save these out for later call and expansion for snakemake
# separating this out for optimal batching
#......................
betaI <-seq(0.1, 1, by = 0.2)
durationI <- seq(3, 15, by = 3)
sirparams <- tidyr::expand_grid(betaI = betaI,
                                durationI = durationI)
SIRpath <- "data/raw_data/SIRparams/sirparam_vals.RDS"
dir.create(sub("sirparam_vals.RDS", "", SIRpath), recursive = T)
saveRDS(sirparams, SIRpath)
#..........
# make Mass Action Model
#..........
outdir <- "data/raw_data/mass_action_network/"
dir.create(outdir, recursive = T)
maoupath <- paste0(outdir, "massactionnetwork.RDS")
manet <- matrix(data = 1, nrow = N, ncol = N)
diag(manet) <- 0
manet <- igraph::graph_from_adjacency_matrix(manet)
saveRDS(manet, file = maoupath)
massactiondf <- tibble::tibble(
  network_manip = "massaction",
  param = "ma",
  val = as.character(0),
  path = maoupath)

#......................
# STEP 1: generate 10 base networks
#......................
outdir <- "data/raw_data/base_networks/"
dir.create(outdir, recursive = T)
basenetpaths <- sapply(1:10, function(x) paste0(outdir, "base_b0_", x, ".RDS"))

# create networks with igraph degree sequence game
base_networks <- replicate(10,
                           igraph::degree.sequence.game(out.deg = rep(50, N),
                             method = "vl"))

# save out
mapply(function(x,y){saveRDS(x, file = y)},
       x = base_networks, y = basenetpaths)

#++++++++++++++++++++++++++++++++++++++++++
### Network Mechanisms        ####
#++++++++++++++++++++++++++++++++++++++++++
# Performing network manipulations on our 10 base models
maestro_dfbase <- tibble::tibble(
  network_manip = "base",
  param = "b",
  val = as.character(0),
  path = basenetpaths) %>%
  dplyr::mutate(num = stringr::str_extract(path, "[0-9]*.RDS"),
                num = stringr::str_replace(num, ".RDS", ""),
                num = as.numeric(num)) %>%
  dplyr::arrange(num) %>%
  dplyr::select(-c("num"))

#++++++++++++++++++++++++++++++++++++++++++
#### Degree Distribution        ####
#++++++++++++++++++++++++++++++++++++++++++
# expand out grid for degree distributions for homogenous and heterogenous variance
degprobdist <- c(0.05, 0.1, 0.15, 0.2, 0.25)
degvardist <- c(0, 1, 5, 10, 25)
dfdegdist <- tidyr::expand_grid(basenetpaths, degprobdist, degvardist) %>%
  magrittr::set_colnames(c("basenetpath", "degprob", 'degvar')) %>%
  tibble::as_tibble()

#......................
# perform manipulations
#......................
# run manipulations via search using wrapper for mem optim
dfdegdist <- dfdegdist %>%
  dplyr::mutate(new_network = purrr::pmap(.,
                                          wrapper_manip_degdist,
                                          .progress = TRUE))

#......................
# write these out
#......................
outdir <- "data/raw_data/degreedist_networks/"
dir.create(outdir, recursive = T)
maestro_dfdegdist <- dfdegdist %>%
  dplyr::mutate(basenetcnt = stringr::str_extract(basenetpath, "[0-9]*.RDS"),
                basenetcnt = stringr::str_replace(basenetcnt, ".RDS", ""),
                basenetcnt = as.numeric(basenetcnt),
                network_manip = "degreedist",
                param = "mean-var",
                val = purrr::map2_chr(degprob, degvar, function(x, y){paste(x, y, sep = "-")}),
                path = paste0(outdir,
                              network_manip, "_",
                              param, "_", val, "_",
                              basenetcnt, ".RDS")) %>%
  dplyr::arrange(param, val, basenetcnt) %>%
  dplyr::select(c("network_manip", "param", "val", "path", "new_network"))

# write out networks out of scope
mapply(function(x, y) { saveRDS(x, file = y) },
       x = maestro_dfdegdist$new_network,
       y = maestro_dfdegdist$path)
# drop for mem
maestro_dfdegdist <- maestro_dfdegdist %>%
  dplyr::select(-c("new_network"))


#++++++++++++++++++++++++++++++++++++++++++
#### Modularity       ####
#++++++++++++++++++++++++++++++++++++++++++
# expand out grid for edges to remove to make modular communities
edgedrop <- c(5, 10, 15, 25, 50, 100, 250, 500, 750, 1000)
dfmodularity <- tidyr::expand_grid(basenetpaths, edgedrop) %>%
  magrittr::set_colnames(c("basenetpath", "edge_rm_num"))


#......................
# perform manipulations
#......................
dfmodularity <- dfmodularity %>%
  dplyr::mutate(new_network = purrr::pmap(.,
                                          wrapper_manip_modular_rmedges,
                                          .progress = TRUE))

#......................
# write these out
#......................
outdir <- "data/raw_data/modularity_networks/"
dir.create(outdir, recursive = T)
maestro_dfmodularity <- dfmodularity %>%
  dplyr::mutate(val = as.character(edge_rm_num)) %>%
  dplyr::mutate(basenetcnt = stringr::str_extract(basenetpath, "[0-9]*.RDS"),
                basenetcnt = stringr::str_replace(basenetcnt, ".RDS", ""),
                basenetcnt = as.numeric(basenetcnt),
                network_manip = "modularity",
                param = "edgerm",
                path = paste0(outdir,
                              network_manip, "_",
                              param, "_", val, "_",
                              basenetcnt, ".RDS")) %>%
  dplyr::arrange(param, val, basenetcnt) %>%
  dplyr::select(c("network_manip", "param", "val", "path", "new_network"))

# write out networks out of scope
mapply(function(x, y) { saveRDS(x, file = y) },
       x = maestro_dfmodularity$new_network,
       y = maestro_dfmodularity$path)
# drop for mem
maestro_dfmodularity <- maestro_dfmodularity %>%
  dplyr::select(-c("new_network"))


#++++++++++++++++++++++++++++++++++++++++++
#### Unity       ####
#++++++++++++++++++++++++++++++++++++++++++
# expand out grid for edges to remove to make modular communities
edgeadd <- c(5, 10, 15, 25, 50, 100, 250, 500, 750, 1000)
dfunity <- tidyr::expand_grid(basenetpaths, edgeadd) %>%
  magrittr::set_colnames(c("basenetpath", "edge_add_num"))

#......................
# perform manipulations
#......................
dfunity <- dfunity %>%
  dplyr::mutate(new_network = purrr::pmap(.,
                                          wrapper_manip_unity_addedges,
                                          .progress = TRUE))

#......................
# write these out
#......................
outdir <- "data/raw_data/unity_networks/"
dir.create(outdir, recursive = T)
maestro_dfunity <- dfunity %>%
  dplyr::mutate(val = as.character(edge_add_num)) %>%
  dplyr::mutate(basenetcnt = stringr::str_extract(basenetpath, "[0-9]*.RDS"),
                basenetcnt = stringr::str_replace(basenetcnt, ".RDS", ""),
                basenetcnt = as.numeric(basenetcnt),
                network_manip = "unity",
                param = "edgeadd",
                path = paste0(outdir,
                              network_manip, "_",
                              param, "_", val, "_",
                              basenetcnt, ".RDS")) %>%
  dplyr::arrange(param, val, basenetcnt) %>%
  dplyr::select(c("network_manip", "param", "val", "path", "new_network"))

# write out networks out of scope
mapply(function(x, y) { saveRDS(x, file = y) },
       x = maestro_dfunity$new_network,
       y = maestro_dfunity$path)
# drop for mem
maestro_dfunity <- maestro_dfunity %>%
  dplyr::select(-c("new_network"))


#++++++++++++++++++++++++++++++++++++++++++
#### Clustering       ####
#++++++++++++++++++++++++++++++++++++++++++
# expand out grid to add edges between dyad pairs and make new clusters (triangles)
clusted <- c(0.01,seq(0.05, 0.25, by = 0.025))
dfclusted <- tidyr::expand_grid(basenetpaths, clusted) %>%
  magrittr::set_colnames(c("basenetpath", "new_transitivity_prob"))


#......................
# perform manipulations
#......................
dfclusted <- dfclusted %>%
  dplyr::mutate(new_network = purrr::pmap(.,
                                          wrapper_manip_clust_edges,
                                          .progress = TRUE))

#......................
# write these out
#......................
outdir <- "data/raw_data/cluster_networks/"
dir.create(outdir, recursive = T)
maestro_dfclusted <- dfclusted %>%
  dplyr::mutate(val = as.character(new_transitivity_prob)) %>%
  dplyr::mutate(basenetcnt = stringr::str_extract(basenetpath, "[0-9]*.RDS"),
                basenetcnt = stringr::str_replace(basenetcnt, ".RDS", ""),
                basenetcnt = as.numeric(basenetcnt),
                network_manip = "cluster",
                param = "transitivity",
                path = paste0(outdir,
                              network_manip, "_",
                              param, "_", val, "_",
                              basenetcnt, ".RDS")) %>%
  dplyr::arrange(param, val, basenetcnt) %>%
  dplyr::select(c("network_manip", "param", "val", "path", "new_network"))

# write out networks out of scope
mapply(function(x, y) { saveRDS(x, file = y) },
       x = maestro_dfclusted$new_network,
       y = maestro_dfclusted$path)
# drop for mem
maestro_dfclusted <- maestro_dfclusted %>%
  dplyr::select(-c("new_network"))


#++++++++++++++++++++++++++++++++++++++++++
#### Neighbor Exchange       ####
#++++++++++++++++++++++++++++++++++++++++++
nexchange_rate <- sapply(seq(9, 0), function(x){10^(-x)})
# this function will capitalize on the internal functionality of `fomes` and
# does not need a static network manipulation as above
maestro_dfNEdyn <- tidyr::expand_grid(basenetpaths, nexchange_rate) %>%
  magrittr::set_colnames(c("path", "val")) %>%
  dplyr::mutate(network_manip = "NEdynamicity",
                param = "exchrate",
                val = as.character(val)) %>%
  dplyr::mutate(basenetcnt = stringr::str_extract(path, "[0-9]*.RDS"),
                basenetcnt = stringr::str_replace(basenetcnt, ".RDS", ""),
                basenetcnt = as.numeric(basenetcnt)) %>%
  dplyr::arrange(param, val, basenetcnt) %>%
  dplyr::select(c("network_manip", "param", "val", "path"))


#++++++++++++++++++++++++++++++++++++++++++
### Bring Together       ####
#++++++++++++++++++++++++++++++++++++++++++
# combine network rows
maestro <- dplyr::bind_rows(massactiondf,
                            maestro_dfbase,
                            maestro_dfdegdist,
                            maestro_dfmodularity,
                            maestro_dfunity,
                            maestro_dfclusted,
                            maestro_dfNEdyn)

# add in outpath
maestro_out <- maestro %>%
  dplyr::mutate(SIRParampath = SIRpath) %>%
  dplyr::mutate(outpath = purrr::pmap_chr(., function(SIRParampath, network_manip, param, val, path){
    paste0("SimRet_", network_manip, param, val, "-", sub(".RDS", "", basename(path)),".RDS")}
  ))


#......................
# save out map file for snakemake
#......................
readr::write_tsv(maestro_out, file = "data/raw_data/maestro_map.tsv")
