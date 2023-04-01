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
### Magic Number Source        ####
#++++++++++++++++++++++++++++++++++++++++++
source("data/raw_data/magic_numbers.R")
set.seed(seed)

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
netbasenames <- sapply(1:10, function(x) paste0("base_b0_", x))
outdir <- "data/raw_data/base_networks/"
dir.create(outdir, recursive = T)

# save out to dir not w/in scope
lapply(netbasenames, mk_base_nets_sout,
       n = N, p.or.m = baseERprob, type = "gnp", outdir = outdir)

#++++++++++++++++++++++++++++++++++++++++++
### Network Mechanisms        ####
#++++++++++++++++++++++++++++++++++++++++++
# Performing network manipulations on our 10 base models
basenetpaths <- list.files(path = outdir, pattern = ".RDS", full.names = T)
maestro_dfbase <- tibble::tibble(
  network_manip = "base",
  param = "b",
  val = 0,
  path = basenetpaths) %>%
  dplyr::mutate(num = stringr::str_extract(path, "[0-9]*.RDS"),
                num = stringr::str_replace(num, ".RDS", ""),
                num = as.numeric(num)) %>%
  dplyr::arrange(num) %>%
  dplyr::select(-c("num"))

#++++++++++++++++++++++++++++++++++++++++++
#### Degree Distribution        ####
#++++++++++++++++++++++++++++++++++++++++++
# expand out grid for degree distributions for homogenous
degprobdist <- seq(0.1, 1, length.out = 10)
dfdegdist_homogen <- expand.grid(basenetpaths, degprobdist,
                                 stringsAsFactors = F) %>%
  magrittr::set_colnames(c("basenetpath", "degprob")) %>%
  dplyr::mutate(degvar = 0)

# expand out grid for degree distributions for heterogenous
degvardist <- c(0.5, 1, 5, 10, 25)
dfdegdist_hetero <- expand.grid(basenetpaths, degprobdist, degvardist,
                                stringsAsFactors = F) %>%
  magrittr::set_colnames(c("basenetpath", "degprob", 'degvar'))

# combine these and setup wrapper for search function
dfdegdist <- dplyr::bind_rows(dfdegdist_homogen, dfdegdist_hetero) %>%
  dplyr::mutate(searches = 1e3) %>%
  tibble::as_tibble()

#......................
# perform manipulations
#......................
# run manipulations via search using wrapper for mem optim
dfdegdist <- dfdegdist %>%
  dplyr::mutate(new_network = purrr::pmap(.,
                                          manip_degdist_wrapper))



#......................
# write these out
#......................
outdir <- "data/raw_data/degreedist_networks/"
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


#++++++++++++++++++++++++++++++++++++++++++
#### Clustering       ####
#++++++++++++++++++++++++++++++++++++++++++




#++++++++++++++++++++++++++++++++++++++++++
#### Neighbor Exchange       ####
#++++++++++++++++++++++++++++++++++++++++++


#++++++++++++++++++++++++++++++++++++++++++
### Bring Together       ####
#++++++++++++++++++++++++++++++++++++++++++
# sanity check
ls()[ grepl(pattern = "maestro_*", ls()) ]
