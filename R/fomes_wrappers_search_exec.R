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
mod <- opt$mod
val <- opt$val
beta <- as.numeric(opt$beta)
durI <- as.numeric(opt$dur)
netgraph <- readRDS(paste0(opt$Rdir, opt$netpath))
conmat <- igraph::as_adjacency_matrix(netgraph, sparse = F)
reps <- 1:opt$reps
output <- opt$output
outdir <- opt$outdir

#++++++++++++++++++++++++++++++++++++++++++
### Fomes Wrapper Function        ####
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
### Adaptive Simulated Annealer Functions ####
#++++++++++++++++++++++++++++++++++++++++++
#' @title Cost
#' @param
#' @description
#' @details
#' @returns
#' @export
cost <- function(mod, beta, durI, val, reps, conmat,
                 outdir, output, itername){

  #......................
  # Running on SimMap
  #......................
  simmap <- tidyr::expand_grid(mod, beta, durI, val, reps) %>%
    dplyr::mutate(conmat = list(conmat)) %>%
    magrittr::set_colnames(c("mod", "beta", "durI", "val", "rep", "conmat"))
  simout <- simmap %>%
    dplyr::mutate(simout = purrr::pmap(., wrap_sim_fomes))
  simoutlite <- simout %>%
    dplyr::mutate(simoutlite = purrr::map(simout, fomes::tidyout)) %>%
    dplyr::select(-c("conmat", "simout"))

  #......................
  # Save Every Run
  # this is not being returned in scope of the function
  #......................
  saveRDS(simout, file = paste0(outdir, "SAiter", itername, "-", output))
  saveRDS(simoutlite, file = paste0(outdir, "lite-SAiter", itername, "-", output))

  #......................
  # CALC COST
  #......................
  costvec <- simoutlite %>%
    dplyr::mutate(finalsize = purrr::map_dbl(simoutlite, "FinalEpidemicSize")) %>%
    dplyr::pull(finalsize)
  costvar <- 1/var(costvec)
  return(costvar)
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
#' @returns
#' @export


adaptive_sim_anneal <- function(maxIter, Iters, AddOnIters, coolingB = 1e-3, initTemp = 1,
                                mod, beta, durI, val, reps, conmat,
                                outdir, output) {
  #......................
  # checks
  #......................
  #goodegg::assert

  #......................
  # setup (const, storage, etc)
  #......................
  costrun <- rep(NA, maxIter)
  temprun <- rep(NA, maxIter)
  betarun <- rep(NA, maxIter)
  durIrun <- rep(NA, maxIter)
  Temp <- initTemp
  currcost <- cost(mod = mod, beta = beta, durI = durI,
                   val = val, reps = reps, conmat = conmat,
                   outdir = outdir, output = output, itername = 0)
  currpos <- c(beta, durI)
  #......................
  # core
  # NB, must make an interactive for loop with while loop because vector allocation prespecified
  #......................
  srch <- TRUE
  i <- 1
  while(srch) {
    # PROPOSE
    newpos <- proposal(currB = currpos[1],
                       currD = currpos[2],
                       iters = Iters)
    # CALC COST
    newcost <- cost(mod = mod, beta = newpos[1], durI = newpos[2],
                    val = val, reps = reps, conmat = conmat,
                    outdir = outdir, output =  output, itername = i)
    # ACCEPT MOVE?
    u <- runif(1)
    p <- min( exp(-(newcost-currcost))/Temp, 1 )
    accept <- u <= p

    ## UPDATES
    # update current position and cost
    currpos <- if(accept){newpos}else{currpos}
    currcost <- if(accept){newcost}else{currcost}
    # Increase temperature if we move from local minima to a new potentially better area
    # DAMPENING: update temp with slow decrease: https://towardsdatascience.com/optimization-techniques-simulated-annealing-d6a4785a1de7
    if (accept) { # dynamically increase search and reignite temp
      Iters <- Iters + AddOnIters
      Temp <- max(Temp, initTemp * (1 - i/maxIter))
    } else {
      Temp <- Temp/(1 + Temp*coolingB)
    }
    # STORAGE items
    costrun[i] <- currcost
    temprun[i] <- Temp
    betarun[i] <- currpos[1]
    durIrun[i] <- currpos[2]
    # DETERMINE SEARCH
    if (i == maxIter) {
      srch <- FALSE
    } else if (i == Iters) {
      srch <- FALSE
    }
    # update i
    i <- i+1
  }

  #......................
  # out
  #......................
  out <- list(
    costrun = costrun[!is.na(costrun)],
    temprun = temprun[!is.na(temprun)],
    betarun = betarun[!is.na(betarun)],
    durIrun = durIrun[!is.na(durIrun)]
  )
  return(out)
}


#++++++++++++++++++++++++++++++++++++++++++
### Run What you Brung ####
#++++++++++++++++++++++++++++++++++++++++++
SAout <- adaptive_sim_anneal(maxIter = 100,
                             Iters = 25,
                             AddOnIters = 5,
                             coolingB = 1e-3,
                             initTemp = 1,
                             mod = mod,
                             beta = beta,
                             durI = durI,
                             val = val,
                             reps = reps,
                             conmat = conmat,
                             outdir = outdir,
                             output = output
                             )
saveRDS(SAout,
        file = paste0(outdir, "SAchain-", output))
