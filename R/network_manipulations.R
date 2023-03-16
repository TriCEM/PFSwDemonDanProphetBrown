#' @title Use Random Search to Manipulate the Degree Distribution of a Network
#' @param graph_network input graph (class igraph)
#' @param degprob numeric; edge density probability per node
#' @description
#' @details
#' @importFrom truncnorm rtruncnorm
#' @importFrom igraph
#' @returns

manip_degdist <- function(graph_network, degprob = 0.5, degvar = 0,
                          searches = 1e3) {
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
  rands <- sample(1:.Machine$integer.max, size = searches, replace = F)


  #......................
  # core
  #......................
  # identify new potential graphs
  # cost is based on edges that are shared
  sapply(rands, function(seednum, new_edge_density, graph_network){
    seed(seednum)
    new_graph_network <- igraph::degree.sequence.game(out.deg = new_edge_density,
                                                      method = "vl")
    overlaps <- igraph::graph.intersection(graph_network, new_graph_network)

    # out
    out <- tibble(randnum = seednum, overlaps = overlaps)
    return(out)
  }, new_edge_density = new_edge_density, graph_network = graph_network)


  #......................
  # out
  #......................
  return()
}





#' @title Identify Potential Edges to Create Triangles in Network
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
  gnet_adjmat <- igraph::as_adjacency_matrix(gnet, sparse = F)
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
