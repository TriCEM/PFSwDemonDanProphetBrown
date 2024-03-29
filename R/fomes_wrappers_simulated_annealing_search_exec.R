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
  retsumm <- summary(ret)
  return(retsumm)
}


#++++++++++++++++++++++++++++++++++++++++++
#### Adaptive Simulated Annealer Functions #####
#++++++++++++++++++++++++++++++++++++++++++
#' @title Cost
#' @param
#' @description
#' @details
#' @returns
#' @export
cost <- function(finalsizes) {
  cst <- 1/var(finalsizes)
  if (is.infinite(cst)) {
    cst <- 1e10 # this is 1/1e-10 which is near the vector limit for R and its ability to calculate variance: `var(c(rep(0,1e9),1))`
  }
  return(cst)
}



#' @title Propose
#' @param
#' @description
#' @details
#' @returns
#' @export
proposal <- function(currB, currD, iters) {
  # since we need this to be a neighborhood of solutions, will only change beta or duration for each prop
  if (rbinom(1,1,0.5)) {
    newB <- round(rnorm(n = 1, mean = currB, sd = (currB * (1-0.23)/sqrt(iters))), digits = 2)
    newB <- ifelse(newB <= 0, 1e-2, newB)
    newD <- currD
  } else {
    newB <- currB
    newD <- round(rnorm(n = 1, mean = currD, sd = (currD * (1-0.23)/sqrt(iters))), digits = 0)
    newD <- ifelse(newD <= 0, 1, newD)
  }
  out <- c(newB, newD)
  return(out)
}

#' @title Run SA
#' @param
#' @description
#' @details Note, in scope behavior of writing out functions
#' @returns Tibble: NB, the chain is the simulated annealing accepted chain versus the "all"
#' which contain the proposed values (whether or not they were accepted)
#' @export


adaptive_sim_anneal <- function(maxIter, Iters, AddOnIters, coolingB = 1e-3, Temp = 1,
                                mod, beta, durI, val, reps, conmat,
                                output) {
  #......................
  # checks
  #......................
  goodegg::assert_single_int(maxIter)
  goodegg::assert_single_int(Iters)
  goodegg::assert_single_int(AddOnIters)
  goodegg::assert_single_numeric(coolingB)
  goodegg::assert_single_numeric(Temp)
  goodegg::assert_single_string(mod)
  goodegg::assert_single_numeric(beta)
  goodegg::assert_single_numeric(durI)
  goodegg::assert_single_string(val)
  goodegg::assert_single_int(reps)
  goodegg::assert_square_matrix(as.matrix(conmat))
  goodegg::assert_single_string(output)

  #......................
  ### Call Seed for Reproducibility
  #......................
  RNGkind(sample.kind = "Rounding")
  # each beta, dur, and network will have it's own set of #rep independent seed
  seeds <- sample(1:1e6, size = (maxIter * reps), replace = F)
  # split up by reps
  seedrun <- split(seeds, factor(sort(rep(1:maxIter, times = reps))))

  #......................
  # setup (const, storage, etc)
  #......................
  costall <- rep(NA, maxIter)
  costchain <- rep(NA, maxIter)
  temprun <- rep(NA, maxIter)
  betaall <- rep(NA, maxIter)
  betachain <- rep(NA, maxIter)
  durIall <- rep(NA, maxIter)
  durIchain <- rep(NA, maxIter)
  SimOutRun <- lapply(1:maxIter, function(x){return(NULL)})
  reps <- 1:reps

  #......................
  # init
  #......................
  initTemp <- Temp # store initial temp for reignition
  i <- 1
  simmap <- tidyr::expand_grid(mod, val, reps) %>%
    dplyr::mutate(beta = beta, # for c/w below code
                  durI = durI,
                  seed = seedrun[[i]],
                  conmat = list(conmat))
  simout <- simmap %>%
    dplyr::mutate(simout = purrr::pmap(., wrap_sim_fomes))
  finalsize <- simout %>%
    dplyr::mutate(simoutlite = purrr::map(simout, fomes::tidyout)) %>%
    dplyr::mutate(finalsize = purrr::map_dbl(simoutlite, "FinalEpidemicSize")) %>%
    dplyr::pull(finalsize)
  currcost <- cost(finalsize)
  currpos <- c(beta, durI)
  # store initial
  costchain[1] <- costall[1] <- currcost # init is same for chain and accept
  temprun[1] <- Temp
  betachain[1] <- betaall[1] <- currpos[1]
  durIchain[1] <- durIall[1] <- currpos[2]
  SimOutRun[[1]] <- finalsize

  #......................
  # core
  # NB, must make an interactive for loop with while loop because vector allocation prespecified
  #......................
  for (i in 2:maxIter){ # R 1 based and we init above
    # PROPOSE
    newpos <- proposal(currB = currpos[1],
                       currD = currpos[2],
                       iters = Iters)
    # CALC COST
    simmap <- tidyr::expand_grid(mod, val, reps) %>%
      dplyr::mutate(beta = newpos[1],
                    durI = newpos[2],
                    seed = seedrun[[i]],
                    conmat = list(conmat))
    simout <- simmap %>%
      dplyr::mutate(simout = purrr::pmap(., wrap_sim_fomes))
    finalsize <- simout %>%
      dplyr::mutate(simoutlite = purrr::map(simout, fomes::tidyout)) %>%
      dplyr::mutate(finalsize = purrr::map_dbl(simoutlite, "FinalEpidemicSize")) %>%
      dplyr::pull(finalsize)
    newcost <- cost(finalsize)

    # update ALL
    costall[i] <- newcost
    betaall[i] <- currpos[1]
    durIall[i] <- currpos[2]
    SimOutRun[[i]] <- finalsize

    # ACCEPT MOVE?
    u <- runif(1)
    p <- min( exp(-(newcost-currcost))/Temp, 1 )
    accept <- u <= p

    ## UPDATES
    # update current position and cost
    currpos <- if(accept){newpos}else{currpos}
    currcost <- if(accept){newcost}else{currcost}
    # REIGNITE: Increase temperature if we move from local minima to a new potentially better area
    # DAMPENING: update temp with exponential decay function
    if (accept) { # dynamically increase search and reignite temp
      Iters <- Iters + AddOnIters
      Temp <- max(Temp, initTemp * (1 - i/maxIter))
    } else {
      Temp <- Temp * exp(-coolingB * i)
    }
    # STORAGE items
    costchain[i] <- currcost
    temprun[i] <- Temp
    betachain[i] <- currpos[1]
    durIchain[i] <- currpos[2]

    # DETERMINE SEARCH, stop if we hit iters mark and are done reigniting
    if (i == Iters) {
      break
    }
  }

  #......................
  # out
  #......................
  seedrun <- seedrun[1:min(Iters, maxIter)]   # drop since we autofilled seeds
  SimOutFinSz <- SimOutRun[!unlist(lapply(SimOutRun, is.null))]
  # bring together in tibble
  out <- tibble::tibble(
    costall = costall,
    costchain = costchain,
    temp = temprun,
    betaall = betaall,
    betachain = betachain,
    durIall = durIall,
    durIchain = durIchain) %>%
    dplyr::filter(!is.na(costall)) # drop if we don't get to max iters

  out$seed <- seedrun
  out$finalsize <- SimOutFinSz
  return(out)
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

  make_option(c("-z", "--maxIter"),
              type = "integer", default = NULL,
              help = paste("Max iters for adaptive simulated annealing search"),
              metavar = "character"),

  make_option(c("-i", "--Iters"),
              type = "integer", default = NULL,
              help = paste("Start iters for adaptive simulated annealing search"),
              metavar = "character"),

  make_option(c("-a", "--AddOnIters"),
              type = "integer", default = NULL,
              help = paste("Additional iters to add for adaptive simulated annealing search"),
              metavar = "character"),

  make_option(c("-c", "--coolingB"),
              type = "numeric", default = NULL,
              help = paste("Cooling factor for slow decrease in adaptive simulated annealing"),
              metavar = "character"),

  make_option(c("-t", "--Temperature"),
              type = "numeric", default = NULL,
              help = paste("Initial temperature for adaptive simulated annealing"),
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
maxIter <- opt$maxIter
Iters <- opt$Iters
AddOnIters <-opt$AddOnIters
coolingB <- opt$coolingB
Temp <- opt$Temperature



#++++++++++++++++++++++++++++++++++++++++++
### Run What you Brung Main ####
#++++++++++++++++++++++++++++++++++++++++++
Start <- Sys.time()
SAmaestro <- tibble::tibble(
  beta = betaI,
  durI = durationI,
  maxIter = maxIter,
  Iters = Iters,
  AddOnIters = AddOnIters,
  coolingB = coolingB,
  Temp = Temp,
  mod = mod,
  val = val,
  reps = reps,
  conmat = list(conmat),
  output = output) %>%
  dplyr::mutate(SAout = purrr::pmap(., adaptive_sim_anneal))
Sys.time() - Start
# tidy up
SAmaestro <- SAmaestro %>%
  dplyr::select(c("beta", "durI", "SAout", "mod", "val")) %>%
  dplyr::mutate(net = opt$netpath)


saveRDS(SAmaestro,
        file = output)

# turn warnings back to default
options(warn = defaultwarnings)
