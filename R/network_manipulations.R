#' @title Identify Potential Edges to Create Triangles in Network
#' @param graph_network input graph (class igraph)
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
