## .................................................................................
## Purpose: Wrapper script for running fomes from command line
##
## Author: Nick Brazeau
## .................................................................................
# Temporarily suppress warning from RNG and "expect" val to be numeric for strict snakemake bash mode
defaultwarnings <- getOption("warn")
options(warn = -1)


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
### Functions        ####
#++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++
#### Fomes Wrapper Function        #####
#++++++++++++++++++++++++++++++++++++++++++
wrap_sim_fomes <- function(seed, mod, beta, durI, val, reps, conmat) {
  # set seed
  set.seed(seed)
  #......................
  # setup (const, storage, etc)
  #......................
  N <- nrow(conmat)
  beta <- rep(beta, N)

  #......................
  # core
  #......................
  if (mod == "NE") {
    ret <- fomes::sim_Gillespie_nSIR(Iseed = 1,
                                     N = N,
                                     beta = beta,
                                     dur_I = durI,
                                     init_contact_mat = conmat,
                                     rho = as.numeric(val),
                                     term_time = Inf,
                                     return_contact_matrices = FALSE)

  } else if (mod == "massaction") {
    maconmat <- matrix(1, N, N)
    diag(maconmat) <- 0
    ret <- fomes::sim_Gillespie_nSIR(Iseed = 1,
                                     N = N,
                                     beta = beta*N,
                                     dur_I = durI,
                                     init_contact_mat = maconmat,
                                     rho = as.numeric(val),
                                     term_time = Inf,
                                     return_contact_matrices = FALSE)



  } else {
    ret <- fomes::sim_Gillespie_nSIR(Iseed = 1,
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
#### Observation Bias Function Wrapper #####
#++++++++++++++++++++++++++++++++++++++++++
#' @title Introduce Bias
#' @inheritParams wrap_sim_fomes
#' @description NO adaptive SA - just single value
#' @details
#' @returns
#' @export

sim_observation_bias <- function(mod, beta, durI, val, reps, conmat,
                                 bias = c(0.01, 0.05, 0.1)) {
  #............................................................
  # checks
  #............................................................
  goodegg::assert_single_string(mod)
  goodegg::assert_single_numeric(beta)
  goodegg::assert_single_numeric(durI)
  goodegg::assert_single_string(val)
  goodegg::assert_single_int(reps)
  goodegg::assert_square_matrix(as.matrix(conmat))
  goodegg::assert_numeric(bias)

  #............................................................
  # setup (const, storage, etc)
  #............................................................
  ### Call Seed for Reproducibility
  RNGkind(sample.kind = "Rounding")
  # each beta, dur, and network - and now bias - will have it's own set of #rep independent seed
  seeds <- sample(1:1e3, size = length(bias)*reps, replace = F)


  #............................................................
  # core
  #............................................................
  #......................
  # introduce bias
  #......................
  biasedconmat <- tibble::tibble(bias = bias,
                                 conmat = list(conmat))

  biasconmatfunx <- function(conmat, bias) {
    newconmat <- conmat <- as.matrix(conmat)
    # lift to just lower tri for symmetry
    conmat <- conmat[lower.tri(conmat)]
    #
    p <- sum(conmat)/length(conmat) # flip w/ orig deg dist density to not alter underlying network properties
    v <- sample(1:length(conmat), floor(length(conmat) * bias), replace = F)
    conmat[v] <- rbinom(n = length(v), size = 1, prob = c(p, 1-p)) # potentially change values and introduce bias
    # now make changes and make it symmetric
    newconmat <- matrix(NA, nrow = nrow(newconmat), ncol = ncol(newconmat))
    newconmat[lower.tri(newconmat)] <- conmat
    newconmat[upper.tri(newconmat)] <- t(newconmat)[upper.tri(newconmat)]
    diag(newconmat) <- 0
    return(newconmat)
  }
  biasedconmat <- biasedconmat %>%
    dplyr::mutate(biasconmat = purrr::pmap(., biasconmatfunx)) %>%
    dplyr::select(-c("conmat")) %>%
    dplyr::rename(conmat = biasconmat)


  #......................
  # expand out simmap
  #......................
  reps <- 1:reps
  simmap <- tidyr::expand_grid(mod, val, reps, beta, durI, bias) %>%
    dplyr::mutate(seed = seeds) %>%
    dplyr::left_join(., biasedconmat, by = "bias")
  # not wrap_fomes doesn't have bias
  simmap$biasfinalsize <- purrr::pmap(simmap[, c("seed", "mod", "beta",
                                                 "durI", "val", "reps", "conmat")],
                                      wrap_sim_fomes)
  #............................................................
  # out
  #............................................................
  return(simmap)
}




#++++++++++++++++++++++++++++++++++++++++++
### Command Line     ####
#++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++
#### parse CL inputs     #####
#++++++++++++++++++++++++++++++++++++++++++
option_list=list(
  make_option(c("-n", "--input"),
              type = "character", default = NULL,
              help = paste("Input filename to read simulation guides from"),
              metavar = "character"),

  make_option(c("-o", "--output"),
              type = "character", default = NULL,
              help = paste("Output filename to write result of Gillespie SIR-NE Model Simulation"),
              metavar = "character"),

  make_option(c("-r", "--reps"),
              type = "integer", default = NULL,
              help = paste("Number of reps to consider"),
              metavar = "character"),

  make_option(c("-b", "--bias"),
              type = "logical", default = NULL,
              help = paste("Whether we are doing bias simulations"),
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


#++++++++++++++++++++++++++++++++++++++++++
#### Unpack from CL        #####
#++++++++++++++++++++++++++++++++++++++++++
simguide <- readRDS(opt$input)
betaI <- simguide$betaI
durationI <- simguide$durationI
mod <- simguide$network_manip
val <- simguide$val
netgraph <- simguide$network[[1]]
conmat <- igraph::as_adjacency_matrix(netgraph, sparse = T)
reps <- opt$reps
output <- opt$output


### Call Seed for Reproducibility
RNGkind(sample.kind = "Rounding")
# each beta, dur, and network - and now bias - will have it's own set of #rep independent seed
fomesseeds <- sample(1:reps*1e2, size = reps, replace = F)



#++++++++++++++++++++++++++++++++++++++++++
### Observation Bias ####
#++++++++++++++++++++++++++++++++++++++++++
# run
biasout <- sim_observation_bias(
  mod = mod,
  beta = betaI,
  durI = durationI,
  val = val,
  reps = reps,
  conmat = list(conmat),
  bias = c(0.01, 0.05, 0.1))
# adjust ouput
biasoutput <- output
biasoutput <- stringr::str_replace(biasoutput, ".RDS", "-BIAStesting.RDS")
# save
saveRDS(biasout, file = biasoutput)


#++++++++++++++++++++++++++++++++++++++++++
### Run What you Brung Main ####
#++++++++++++++++++++++++++++++++++++++++++
Start <- Sys.time()
# mk tbl
runmaestro <- tibble::tibble(
  beta = betaI,
  durI = durationI,
  val = val,
  reps = reps,
  conmat = list(conmat),
  output = output)
# get seends
theseseeds <- fomesseeds[1:nrow(runmaestro)]
# run
runmaestro <- runmaestro %>%
  dplyr::mutate(seed = theseseeds) %>%
  dplyr::relocate(seed) %>%
  dplyr::mutate(fomesout = purrr::pmap(., wrap_sim_fomes))
Sys.time() - Start
# tidy up
runmaestro <- runmaestro %>%
  dplyr::mutate(net = opt$netpath) %>%
  dplyr::select(c("seed", "beta", "durI", "mod", "val", "net", "fomesout"))


saveRDS(runmaestro,
        file = output)

# turn warnings back to default
options(warn = defaultwarnings)
