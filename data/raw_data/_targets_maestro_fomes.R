## .................................................................................
## Purpose: Generate `fomes` stochastic runs for analysis
##
## Author: Nick Brazeau
##
## Date: 11 March, 2023
##
## Notes: Using the outline of PFS Manuscript
## .................................................................................
library(targets)
library(tarchetypes)
library(tibble)
library(igraph)

#++++++++++++++++++++++++++++++++++++++++++
### Magic Number Source        ####
#++++++++++++++++++++++++++++++++++++++++++
source("data/raw_data/magic_numbers.R")
set.seed(seed)

#++++++++++++++++++++++++++++++++++++++++++
### Setup Target Specifications        ####
#++++++++++++++++++++++++++++++++++++++++++
# Set target options:
tar_option_set(
  packages = c("tidyverse", "fomes"), # package dependencies
  format = "rds",
  memory = "transient",
  garbage_collection = TRUE, # optimize mem
  storage = "worker") # offload main script from bottlenecking
options(clustermq.scheduler = "SLURM",
        clustermq.template = "/nas/longleaf/home/nfb/clustermq/clustermq_template.slurm") # my location for clustermq

# Run the R scripts in the R/ folder with your custom functions:
tar_source("R/")

#++++++++++++++++++++++++++++++++++++++++++
### Maestro Setup                      ####
#++++++++++++++++++++++++++++++++++++++++++
# Based on discussion want to have the following structure:
#  * 25 beta values
#  * 25 duration of illness
#  * 10 "base" ER networks that will be adjusted by network mechanisms (default is p = 0.3)
#    * Network mechanisms:
#      * 10 degree distributions
#      * 10 modularity models
#      * 10 clustering
#      * 10 neighbor exchange rate

#......................
# STEP 0: generate beta and duration of illness
#......................
betaI <- seq(0.1, 1, length.out = 25)
durationI <- seq(2,50, by = 2)

#......................
# STEP 1: generate 10 base networks
#......................
netbasenames <- sapply(1:10, function(x) paste0("BN", x))
outdir <- "data/raw_data/base_networks/"
dir.create(outdir, recursive = T)
# simple function to save out nets
mk_base_nets_sout <- function(nm, n, p.or.m, type, outdir) {
  # ER game
  out <- igraph::erdos.renyi.game(
    n = n,
    p.or.m = p.or.m,
    type = type )
  # save out
  saveRDS(out, paste0(outdir, nm, ".RDS"))
  }

# save out to dir not w/in scope
lapply(netbasenames, mk_base_nets_sout,
       n = N, p.or.m = baseERprob, type = "gnp", outdir = outdir)

#......................
# make base maestro
#......................
basenetpaths <- list.files(outdir, pattern = ".RDS", full.names = T)
base_maestro <- tibble::as_tibble(expand.grid(betaI, durationI, basenetpaths)) %>%
  magrittr::set_colnames(c("betaI", "durI", "basenetpaths")) %>%
  dplyr::mutate(basenetpaths = as.character(basenetpaths),
                basenetnames = sub(".RDS", "", basename(basenetpaths)))

#++++++++++++++++++++++++++++++++++++++++++
### Network Mechanisms        ####
#++++++++++++++++++++++++++++++++++++++++++
# Performing network manipulations on our 10 base models

#++++++++++++++++++++++++++++++++++++++++++
#### Degree Distribution        ####
#++++++++++++++++++++++++++++++++++++++++++
degdist <- seq(0.1, 1, length.out = 10)
dfdegdist <- expand.grid(basenetpaths, degdist) %>%
  magrittr::set_colnames(c("basenetpath", "degdistit"))

