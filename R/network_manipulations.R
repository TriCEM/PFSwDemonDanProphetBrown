## .................................................................................
## Purpose: Functions for creating the network manipulations used within maestro
##
## Notes:
## .................................................................................

#' @title Generate Base Networks with Degree Sequence Algorithm
#' @inheritParams igraph degree.sequence.game
#' @param n integer; Node count
#' @param outdir char; path base for base network
#' @description simple function to save out nets
#' @importFrom igraph erdos.renyi.game
#' @returns igraph network and writes it out to outdir
#
mk_base_nets_sout <- function(n, out.deg, outdir) {
  # ER game
  out <- igraph::degree.sequence.game(out.deg = rep(degden, n),
                                      method = "vl", )
  # save out
  saveRDS(out, paste0(outdir, nm, ".RDS"))
}



#' @title Finite CRP for Network Code
#' @param graph_network input graph (class igraph)
#' @param edge_delta integer; number of edges to add or remove
#' @param theta numeric; concentration parameter
#' @description
#' @details Assumes that there is already a base categorical distribution
#' and that all of the finite tables have been identified, and we are just
#' adding individuals to those existing tables
#' @returns New Degree Distribution

finite_crp_network <- function(graph_network, edge_delta, theta) {
  #......................
  # checks
  #......................
  goodegg::assert_eq("igraph", class(graph_network),
                     message = "The graph_network object must be have the igraph class (i.e. generate network with igraph)")
  goodegg::assert_single_numeric(theta)
  goodegg::assert_single_pos(theta)
  goodegg::assert_single_int(edge_delta)

  #......................
  # setup (const, storage, etc)
  #......................
  x <- degree(graph_network)
  # catch if removing edges versus adding
  if (edge_delta < 0) {
    x <- x * -1
  }
  #......................
  # core
  #......................
  for (i in 1:edge_delta) {
    prd <- x/(i - 1 + theta)
    td <- sample(1:length(x), 1, prob = prd)
    x[td] <- x[td] + 1
  }
  #......................
  # out
  #......................
  return(abs(x))
}


#' @title Randomly Adds or Removes Edges to Manipulate the Degree Distribution of a Network
#' @param graph_network input graph (class igraph)
#' @param new_degprob numeric; edge density probability per node
#' @param new_degvar numeric; edge density probability per node
#' @details
#' @importFrom truncnorm rtruncnorm
#' @importFrom igraph degree.sequence.game, intersection, ecount
#' @returns list containing dataframe of searches for networks and the best network

manip_degdist <- function(graph_network, new_degprob = 0.5,
                          new_degvar = 5) {
  #......................
  # checks
  #......................
  goodegg::assert_eq("igraph", class(graph_network),
                     message = "The graph_network object must be have the igraph class (i.e. generate network with igraph)")


  #......................
  # setup (const, storage, etc)
  #......................
  # identify edges that are needed
  nodecount <- igraph::vcount(graph_network)
  curr_edges <- igraph::ecount(graph_network)
  new_edges <- round( new_degprob * (nodecount * (nodecount - 1) / 2 ))
  need_ed <- new_edges - curr_edges
  # plan to manipulate in place
  new_graph_network <- graph_network
  # identify all potential edges
  pot_edges <- igraph::simplify(igraph::complementer(graph_network))

  # catch zero
  if (new_degvar == 0) {
    new_degvar <- .Machine$double.xmin
  }
  #......................
  # core
  #......................
  # account for edge dispersion based on variance
  node_degchanges <- finite_crp_network(graph_network = graph_network,
                                        edge_delta = need_ed,
                                        theta = 1/new_degvar)

  # add or delete edges accordingly
  node_degchanges <- node_degchanges - degree(graph_network)
  for (i in 1:length(node_degchanges)) {
    #......................
    # ADDING edges
    #......................
    if (node_degchanges[i] > 0) {
      # identify edges with node of interest
      edges_with_node <- igraph::incident(graph = pot_edges, v = node_degchanges[i])
      # catch if we need to add more edges than is possible
      if (length(edges_with_node) < node_degchanges[i]) {
        edges_to_add <- edges_with_node
      } else {
        edges_to_add <- sample(x = edges_with_node, size = node_degchanges[i])
      }
      # class liftover
      edges_to_add <- igraph::ends(pot_edges, edges_to_add)
      # make changes
      graph_network <- igraph::add_edges(graph = graph_network,
                                         edges = edges_to_add)

      #......................
      # Removing edges
      #......................
    } else {
      # identify edges with node of interest
      edges_with_node <- igraph::incident(graph = graph_network, v = node_degchanges[i])
      # catch if we need to add more edges than is possible
      if (length(edges_with_node) < node_degchanges[i]) {
        edges_to_rm <- edges_with_node
      } else {
        edges_to_rm <- sample(x = edges_with_node, size = node_degchanges[i])
      }
      # class liftover
      edges_to_rm <- igraph::ends(pot_edges, edges_to_rm)
      # make changes
      graph_network <- igraph::delete_edges(graph = graph_network,
                                            edges = edges_to_rm)
    }
  }

  #......................
  # out
  #......................
  return(graph_network)
}




#' @title Wrapper for `manip_degdist`
#' @noMd
#' @return igraph network
wrapper_manip_degdist <- function(basenetpath, new_degprob, new_degvar) {
  grphnet <- readRDS(basenetpath)
  out <- manip_degdist(graph_network = grphnet,
                       new_degprob = degprob,
                       new_degvar = degvar)
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
