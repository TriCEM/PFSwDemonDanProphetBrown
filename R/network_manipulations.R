## .................................................................................
## Purpose: Functions for creating the network manipulations used within maestro
##
## Notes:
## .................................................................................

#' @title Generate Base Networks with Erdos Renyi Algorithm
#' @inheritParams igraph erdos.renyi.game
#' @param nm numeric; base network name
#' @param outdir char; path base for base network
#' @description simple function to save out nets
#' @importFrom igraph erdos.renyi.game
#' @returns igraph network and writes it out to outdir
#
mk_base_nets_sout <- function(nm, n, p.or.m, type, outdir) {
  # ER game
  out <- igraph::erdos.renyi.game(
    n = n,
    p.or.m = p.or.m,
    type = type )
  # save out
  saveRDS(out, paste0(outdir, nm, ".RDS"))
}


#' @title Use Random Search to Manipulate the Degree Distribution of a Network
#' @param graph_network input graph (class igraph)
#' @param degprob numeric; edge density probability per node
#' @description Function to find most overlaps when manipulating the degree
#' distribution of an original network
#' @details
#' @importFrom truncnorm rtruncnorm
#' @importFrom igraph degree.sequence.game, intersection, ecount
#' @returns list containing dataframe of searches for networks and the best network

manip_degdist <- function(graph_network, degprob = 0.5, degvar = 0, searches = 1e3) {
  #......................
  # checks
  #......................
  goodegg::assert_eq("igraph", class(graph_network),
                     message = "The graph_network object must be have the igraph class (i.e. generate network with igraph)")


  #......................
  # setup (const, storage, etc)
  #......................
  # new degree sequence based on mean and variance
  nodecount <- vcount(graph_network)
  new_edge_density <- sapply(1:nodecount, function(x, m, v){
    out <- truncnorm::rtruncnorm(n = 1, a = 0, b = 1, mean = m, sd = sqrt(v))
    return(out)
  }, m = degprob, v = degvar)
  new_edge_density <- round( new_edge_density * nodecount )
  # generate random numbers for search
  seednum <- sample(1:searches*1e2, size = searches, replace = F)

  #......................
  # core
  #......................
  ## FUNCTION for doing search
  # identify new potential graphs based on overlaps to curr graph
  identify_potent_deggraphs <- function(seednum, new_edge_density,
                                        graph_network) {
    # seed
    set.seed(seednum)
    # make new graph
    new_graph_network <- igraph::degree.sequence.game(out.deg = new_edge_density,
                                                      method = "vl")
    # calculate overlaps --> cost is based on edges that are shared
    overlaps <- igraph::graph.intersection(graph_network, new_graph_network)
    overlaps <- igraph::ecount(overlaps)

    # out
    return(overlaps)
  }

  # put together in dfmap
  pot_graphs_df <- tibble::tibble(seednum = seednum,
                                  new_edge_density = list(new_edge_density),
                                  graph_network = list(graph_network))
  # do search
  pot_graphs_df <- pot_graphs_df %>%
    dplyr::mutate(overlaps = purrr::pmap_dbl(., identify_potent_deggraphs)) %>%
    dplyr::select(c("seednum", "overlaps"))


  # identify best graph
  bst_sd <- pot_graphs_df %>%
    dplyr::filter(overlaps == max(overlaps)) %>%
    dplyr::pull(seednum)

  # catch multiple maxes
  if (length(bst_sd) > 1) {
    rw <- sample(1:length(bst_sd), size = 1)
    bst_sd <- bst_sd[rw]
  }
  # set seed
  set.seed(bst_sd)
  bst_ovlp_grph <- igraph::degree.sequence.game(out.deg = new_edge_density,
                                                method = "vl")
  #......................
  # out
  #......................
  # tidy up
  out <- list(best_overlapping_graph = bst_ovlp_grph,
              potential_graphs_df = pot_graphs_df)
  return(out)
}

#' @title Wrapper for `manip_degdist`
#' @noMd
#' @return igraph network
wrapper_manip_degdist <- function(basenetpath, degprob, degvar, searches) {
  grphnet <- readRDS(basenetpath)
  out <- manip_degdist(graph_network = grphnet,
                       degprob = degprob,
                       degvar = degvar,
                       searches = searches)
  out <- out$best_overlapping_graph
  return(out)
}



#' @title Increase Network Modularity by Removing Connected Edges
#' @inheritParams manip_degdist
#' @param edge_rm_num int; Number of edges to remove from current graph
#' @description To increase modularity of a given network, a number of well-connected
#' edges are removed. Edge connectedness is determined by the \code{igraph::edge_betweenness}
#' function. The process is deterministic with the most connected edges being removed
#' sequentially
#' @returns network graph (class igraph) with new modularity

manip_modular_rmedges <- function(graph_network, edge_rm_num) {
  #......................
  # checks
  #......................
  goodegg::assert_eq("igraph", class(graph_network),
                     message = "The graph_network object must be have the igraph class (i.e. generate network with igraph)")
  goodegg::assert_single_int(edge_rm_num)

  #......................
  # core
  #......................
  # identify most connected edges by betweeness
  btwnconn <- igraph::edge_betweenness(graph_network, directed = F)
  # drop num edges specified based on most connected
  btwnord <- rev( order(btwnconn) )[1:edge_rm_num]
  # drop edges
  new_mod <- igraph::delete_edges(graph = graph_network,
                                  edges = E(graph_network)[btwnord])
  #......................
  # out
  #......................
  return(new_mod)
}




#' @title Wrapper for `manip_modular_rmedges`
#' @noMd
#' @return igraph network
wrapper_manip_modular_rmedges <- function(basenetpath, edge_rm_num) {
  grphnet <- readRDS(basenetpath)
  out <- manip_modular_rmedges(graph_network = grphnet,
                               edge_rm_num = edge_rm_num)
  return(out)
}



#' @title Identify Potential Edges to Create Triangles (Clustering) in Network
#' @inheritParams manip_degdist
#' @description Identifies potential edges that would create triangles from
#' the input network. Returns those edges as an A(i,j) listing, where A is the adjacency
#' matrix
#' @details Based on premise of squaring adjacency matrices to identify nodes with
#' path lengths of two
#' @returns list of potential edges that would introduce new triangles in the graph

get_potential_triangle_edges <- function(graph_network) {
  #......................
  # checks
  #......................
  goodegg::assert_eq("igraph", class(graph_network),
                     message = "The graph_network object must be have the igraph class (i.e. generate network with igraph)")

  #......................
  # setup
  #......................
  # storage
  newtri_conns <- list()
  # squared adjacency matrix gives all paths of 2 (Gilbert Strang, pg 78)
  gnet_adjmat <- igraph::as_adjacency_matrix(graph_network, sparse = F)
  gnet_adjmatsq <- gnet_adjmat %*% gnet_adjmat
  # drop out elements we don't need to avoid redundancies (matrix is symmetric)
  diag(gnet_adjmatsq) <- 0
  gnet_adjmatsq[upper.tri(gnet_adjmatsq)] <- 0

  #......................
  # core
  #......................
  # corner case and normal case
  if ( all(gnet_adjmatsq[lower.tri(gnet_adjmatsq)] %in% c(0,1)) ) { # catch corner case of line or perimeter
    for(i in 1:nrow(gnet_adjmatsq)) {
      for (j in 1:ncol(gnet_adjmatsq)) {
        if(gnet_adjmatsq[i,j] == 1) {
          newtri_conns <- append(newtri_conns, list(c(i,j)))
        }
      }
    } # end nested for loop corner case

  } else { # "normal" case

    for(i in 1:nrow(gnet_adjmatsq)) {
      for (j in 1:ncol(gnet_adjmatsq)) {
        if(gnet_adjmatsq[i,j] == 2) {
          newtri_conns <- append(newtri_conns, list(c(i,j)))
        }
      }
    } # end nested for loop
  } # end if/else corner case catch

  #......................
  # out
  #......................
  return(newtri_conns)
}


#' @title Add New Edges to Induce Clustering (Triangles)
#' @inheritParams manip_degdist
#' @param edge_add_num int; Number of edges to add
#' @description Identify a subset of dyad pairs that could have a new edge
#' introduced between them to induce a cluster, or new triangle.
#' @details Adds edges sequentially based on node numbering
#' @importFrom PFSwDemonDanProphetBrown get_potential_triangle_edges
#' @returns network graph (igraph class) with new triangle clusters

manip_clust_addedges <- function(graph_network, edge_add_num) {
  #......................
  # checks
  #......................
  goodegg::assert_eq("igraph", class(graph_network),
                     message = "The graph_network object must be have the igraph class (i.e. generate network with igraph)")
  goodegg::assert_single_int(edge_add_num)

  #......................
  # setup (const, storage, etc)
  #......................
  new_graph_network <- graph_network
  #......................
  # core
  #......................
  potedges <- get_potential_triangle_edges(graph_network)
  # catch
  if (length(potedges) > edge_add_num) {
    stop("You have requested to add more edges than there are potential dyad pairs.
         Either submit a new graph network or decrease the number of requested
         edges")
  }
  # add edges in place
  for (i in 1:edge_add_num) {
    new_graph_network <- igraph::add_edges(graph = new_graph_network,
                                           edges = potedges[[i]])
  }

  #......................
  # out
  #......................
  return(new_graph_network)
}

#' @title Wrapper for `manip_clust_addedges`
#' @noMd
#' @return igraph network
wrapper_manip_clust_addedges <- function(basenetpath, edge_add_num) {
  grphnet <- readRDS(basenetpath)
  out <- manip_clust_addedges(graph_network = grphnet,
                              edge_add_num = edge_add_num)
  return(out)
}
