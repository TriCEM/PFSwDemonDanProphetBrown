## .................................................................................
## Purpose: Wrapper script for running fomes from command line
##
## Author: Nick Brazeau
##
## Notes:
## .................................................................................

#++++++++++++++++++++++++++++++++++++++++++
### dependencies     ####
#++++++++++++++++++++++++++++++++++++++++++
deps <- c("fomes", "goodegg", "optparse", "dplyr", "tibble", "magrittr", "igraph")
deps <- !sapply(deps, function(x){x %in% installed.packages()[,1]} ) # note this a named vector

# catch fomes remote
if(deps["fomes"]) {
  if (!"remotes" %in% installed.packages()[,1]){
    install.packages("remotes")
  }
  remotes::install_github("TriCEM/fomes")
  deps <- deps[names(deps) != "fomes"]
}

# catch goodegg remote
if(deps["goodegg"]) {
  if (!"remotes" %in% installed.packages()[,1]){
    install.packages("remotes")
  }
  remotes::install_github("nickbrazeau/goodegg")
  deps <- deps[names(deps) != "goodegg"]
}

# rest of deps
if (any(deps)) {
  install.packages(names(deps)[deps])
}

#......................
# call dependencies
#......................
library(fomes)
library(goodegg)
library(igraph)
library(optparse)
library(dplyr)
library(tibble)
library(magrittr)

#++++++++++++++++++++++++++++++++++++++++++
### parse CL inputs     ####
#++++++++++++++++++++++++++++++++++++++++++
option_list=list(

  make_option(c("-m", "--mod"),
              type = "character", default = NULL,
              help = paste("Model type"),
              metavar = "character"),


  make_option(c("-b", "--betaI"),
              type = "character", default = NULL,
              help = paste("File path for beta values to iterate over"),
              metavar = "character"),

  make_option(c("-d", "--dur"),
              type = "character", default = NULL,
              help = paste("File path for duration of infection values to iterate over"),
              metavar = "character"),

  make_option(c("-n", "--netpath"),
              type = "character", default = NULL,
              help = paste("File path for network to consider"),
              metavar = "character"),

  make_option(c("-v", "--val"),
              type = "numeric", default = NULL,
              help = paste("Value of NE for specific Gillespie SIR-NE Model Simulation"),
              metavar = "character"),

  make_option(c("-r", "--reps"),
              type = "integer", default = NULL,
              help = paste("Number of reps to consider"),
              metavar = "character"),

  make_option(c("-o", "--output"),
              type = "character", default = NULL,
              help = paste("Output filename to write result of Gillespie SIR-NE Model Simulation"),
              metavar = "character"),

  make_option(c("-x", "--outdir"),
              type = "character", default = NULL,
              help = paste("Output directory to write result of Gillespie SIR-NE Model Simulation"),
              metavar = "character"),

  make_option(c("-s", "--Rdir"),
              type = "character", default = NULL,
              help = paste("Root directory of R project"),
              metavar = "character")

)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

#++++++++++++++++++++++++++++++++++++++++++
### Unpack from CL        ####
#++++++++++++++++++++++++++++++++++++++++++
beta <- readRDS(paste0(opt$Rdir, opt$beta))
durI <- readRDS(paste0(opt$Rdir, opt$dur))
netgraph <- readRDS(paste0(opt$Rdir, opt$netpath))
conmat <- igraph::as_adjacency_matrix(netgraph, sparse = F)
reps <- 1:opt$reps

#++++++++++++++++++++++++++++++++++++++++++
### Make Map Dataframe to iterate over        ####
#++++++++++++++++++++++++++++++++++++++++++
simmap <- tibble::as_tibble(expand.grid(opt$mod, beta, durI, opt$val, reps,
                                        stringsAsFactors = F)) %>%
  dplyr::mutate(conmat = list(conmat)) %>%
  magrittr::set_colnames(c("mod", "beta", "durI", "val", "rep", "conmat"))


#++++++++++++++++++++++++++++++++++++++++++
### Wrapper Function        ####
#++++++++++++++++++++++++++++++++++++++++++
wrap_sim_fomes <- function(mod, beta, durI, val, reps, conmat) {

  #......................
  # setup (const, storage, etc)
  #......................
  N <- nrow(conmat)
  beta <- rep(beta, N)

  #......................
  # core
  #......................
  if (mod == "NE") {
    ret <- fomes::sim_Gillespie_SIR(Iseed = 1,
                             N = N,
                             beta = beta,
                             dur_I = durI,
                             init_contact_mat = conmat,
                             rho = as.numeric(val),
                             term_time = Inf,
                             return_contact_matrices = FALSE)

  } else {
    ret <- fomes::sim_Gillespie_SIR(Iseed = 1,
                             N = N,
                             beta = beta,
                             dur_I = durI,
                             init_contact_mat = conmat,
                             rho = .Machine$double.xmin,
                             term_time = Inf,
                             return_contact_matrices = FALSE)
  }

  #......................
  # out
  #......................
  return(ret)
}


#++++++++++++++++++++++++++++++++++++++++++
### Run        ####
#++++++++++++++++++++++++++++++++++++++++++
simout <- simmap %>%
  dplyr::mutate(simout = purrr::pmap(., wrap_sim_fomes,
                                     .progress = TRUE))

#++++++++++++++++++++++++++++++++++++++++++
### Save        ####
#++++++++++++++++++++++++++++++++++++++++++
saveRDS(simout, file = paste0(opt$outdir, opt$output))
