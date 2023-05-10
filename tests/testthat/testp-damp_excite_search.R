## .................................................................................
## Purpose: test script for a dampening-excitation search algorithm
##
## Author: Nick Brazeau
##
## Date: 14 April, 2023
##
## Notes:
## .................................................................................

#............................................................
# create data that is difficult to find a solution
#...........................................................
#set.seed(42)
# paramdim <- 1e2
# poor_landscape <- matrix(runif(n = (paramdim+1)^2, min = 100, max = 500), paramdim+1, paramdim+1)
# # going to make region in interest at center
# poor_landscape[40:60,40:60] <- 100
# poor_landscape[45:55,45:55] <- 75
# poor_landscape[48:52,48:52] <- 50
# poor_landscape[49:50,49:50] <- 10
#image(poor_landscape)
paramdim <- 1e3
poor_landscape <- matrix(runif(n = (paramdim+1)^2, min = 100, max = 500), paramdim+1, paramdim+1)
# going to make region in interest at center
poor_landscape[480:520,480:520] <- 50
poor_landscape[499:500,499:500] <- 10
image(poor_landscape)




#............................................................
# simulated annealing function
#...........................................................
cost <- function(i,j,landscape){
  return(landscape[i,j])
}

proposal <- function(curri, currj, Boundi, Boundj) {
  # since we need this to be a neighborhood of solutions, will only change i or j for each prop
  if (rbinom(1,1,0.5)) {
    newi <- round(runif(1, min = 1, max = Boundi), digits = 0)
    newj <- currj
  } else {
    newi <- curri
    newj <- round(runif(1, min = 1, max = Boundj), digits = 0)
  }
  out <- c(newi, newj)
  return(out)
}

#' @title
#' @param
#' @description
#' @details
#' @returns
#' @export


adaptive_sim_anneal <- function(Initi, Initj, Boundi, Boundj, landscape,
                                maxIter, Iters, AddOnIters, coolingB = 1e-3, initTemp = 1) {
  #......................
  # checks
  #......................
  #goodegg::assert

  #......................
  # setup (const, storage, etc)
  #......................
  costrun <- rep(NA, maxIter)
  temprun <- rep(NA, maxIter)
  irun <- rep(NA, maxIter)
  jrun <- rep(NA, maxIter)
  Temp <- initTemp
  currcost <- cost(i = Initi, j = Initj, landscape)
  currpos <- c(Initi, Initj)
  #......................
  # core
  # NB, must make an interactive for loop with while loop because vector allocation prespecified
  #......................
  srch <- TRUE
  i <- 1
  while(srch) {
    # PROPOSE
    newpos <- proposal(curri = currpos[1],
                       currj = currpos[2],
                       Boundi = Boundi, Boundj = Boundj)
    # CALC COST
    newcost <- cost(i = newpos[1], j = newpos[2], landscape)
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
    irun[i] <- currpos[1]
    jrun[i] <- currpos[2]
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
    irun = irun[!is.na(irun)],
    jrun = jrun[!is.na(jrun)]
  )
  return(out)
}



#............................................................
# run
#...........................................................
#debug(adaptive_sim_anneal)
out <- adaptive_sim_anneal(Initi = 1, Initj = 1,
                           Boundi = paramdim, Boundj = paramdim,
                           landscape = poor_landscape,
                           maxIter = 1e3,
                           Iters = 1e2,
                           AddOnIters = 1e1,
                           coolingB = 1e-3, initTemp = 1)
length(out$costrun)


#............................................................
# plot
#...........................................................
library(tidyverse)
library(gganimate)
i <- rep(1:1001, 1001)
j <- sort(rep(1:1001, 1001))
val <- as.vector(poor_landscape)
# landscape
landscapeplotdf <- tibble::tibble(i = i, j = j, val = val)
# add init
initdf <- tibble::tibble(rep = 0, i  = 1, j = 1, temp = 1, cost = NA)

# search
searchdf <- tibble::tibble(rep = 1:length(out$temprun),
                           i = out$irun,
                           j = out$jrun,
                           temp = out$temprun,
                           cost = out$costrun)
searchdf <- dplyr::bind_rows(initdf, searchdf)

searchdf <- searchdf %>%
  dplyr::mutate(iend = lead(i),
                jend = lead(j))




# plot
plotObj <- ggplot() +
  geom_tile(data = landscapeplotdf,
            aes(x = i, y = j, fill = val)) +
  geom_segment(data = searchdf,
            aes(x = i, xend = iend, y = j, yend = jend,
                color = temp),
            alpha = 0.5, linewidth = 1,
            arrow = arrow(type = "closed",
                          ends = "first",
                          length = unit(0.07, "cm"))) +
  scale_fill_viridis_c() +
  scale_color_gradient("Search Temp",
                       low = "#fee0d2", high = "#cb181d") +
  theme(plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right",
        legend.title = element_text(family = "Helvetica", face = "bold", vjust = 0.85, size = 12),
        legend.text = element_text(family = "Helvetica", hjust = 0.5, vjust = 0.5, size = 10),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank())

#plotObj
plotAnim <- plotObj + gganimate::transition_reveal(rep)

#gganimate::animate(anim1, fps = 1)
gganimate::anim_save(animation = plotAnim,
                     filename = "~/Desktop/adaptive_SA_search.gif")

# prop vs acceptned




