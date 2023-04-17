## .................................................................................
## Purpose: Functions for creating the network manipulations used within maestro
##
## Notes:
## .................................................................................

#' @title Stick Breaking Process to Spread Edges
#' @inheritParams manip_degdist
#' @param edge_delta integer; Number of edges to spread out
#' @description Use a stick-breaking process from a beta distribution
#' with alpha = 1 and beta = new_degvar
#' @importFrom truncnorm rtruncnorm
#' @returns Vector of integers that correspond to edge changes
#' @export

stickbreak_network <- function(graph_network,
                               edge_delta,
                               new_degvar) {
  #......................
  # checks
  #......................
  goodegg::assert_single_int(edge_delta)
  goodegg::assert_single_pos(new_degvar)

  #......................
  # setup (const, storage, etc)
  #......................
  nodecount <- igraph::vcount(graph_network)
  remain_stick_length <- 1
  wi <- rep(0, nodecount)

  #......................
  # core
  #......................
  if (new_degvar == 0) {
    edge_changes <- round(edge_delta/nodecount)
    edge_changes <- rep(edge_changes, nodecount)

  } else {
    #......................
    # draw weights through stick break
    #......................
    # run stick break
    for (i in 1:nodecount) {
      tsamp <- rbeta(1, 1, new_degvar)
      wi[i] <- remain_stick_length * tsamp
      remain_stick_length <- remain_stick_length * (1 - tsamp)
    } # end for loop
    # spread out edges
    # NB, by rounding, we may end up with more or fewer edges than
    # edge_delta, this should be minimal
    edge_changes <- round( edge_delta * wi )
  } # end else

  #......................
  # out
  #......................
  return(edge_changes)
}





#' @title Randomly Adds or Removes Edges to Manipulate the Degree Distribution of a Network
#' @param graph_network input graph (class igraph)
#' @param new_degprob numeric; edge density probability per node
#' @param new_degvar numeric; edge density probability per node
#' @details
#' @return graph network updated for degree distribution
#' @export

manip_degdist <- function(graph_network, new_degprob = 0.5,
                          new_degvar = 5) {
  #......................
  # checks
  #......................
  goodegg::assert_eq("igraph", class(graph_network),
                     message = "The graph_network object must be have the igraph class (i.e. generate network with igraph)")
  goodegg::assert_single_pos(new_degprob)
  goodegg::assert_single_pos(new_degvar)


  #......................
  # setup (const, storage, etc)
  #......................
  # identify edges that are needed
  nodecount <- igraph::vcount(graph_network)
  curr_edges <- igraph::ecount(graph_network)
  new_edges <- round( new_degprob * (nodecount * (nodecount - 1) / 2 ))
  need_ed <- new_edges - curr_edges

  #......................
  # core
  #......................
  # account for edge dispersion based on variance
  node_degchanges <- stickbreak_network(graph_network = graph_network,
                                        edge_delta = need_ed,
                                        new_degvar = new_degvar)

  # identify all potential edges that we could add
  pot_graph <- igraph::simplify(igraph::complementer(graph_network))
  # add or delete edges accordingly
  for (i in 1:length(node_degchanges)) {
    # identify edges with node of interest
    edges_with_node <- igraph::incident(graph = pot_graph, v = i)
    #......................
    # ADDING edges
    #......................
    if (node_degchanges[i] > 0) {
      # catch if we need to add more edges than is possible
      if (length(edges_with_node) < abs(node_degchanges[i])) {
        edges_to_add <- edges_with_node
      } else {
        edges_to_add <- sample(x = edges_with_node, size = abs(node_degchanges[i]))
      }
      # class liftover
      edges_to_add <- igraph::ends(pot_graph, edges_to_add)
      edges_to_add <- unlist(apply(edges_to_add, 1, as.vector, simplify = F))
      # make changes
      graph_network <- igraph::add_edges(graph = graph_network,
                                         edges = edges_to_add)
      #......................
      # Removing edges
      #......................
    } else if (node_degchanges[i] < 0) {
      # don't need complement here, as we are removing existing edges
      # identify edges with node of interest
      edges_with_node <- igraph::incident(graph = graph_network, v = i)
      # catch if we need to rm more edges than is possible
      if (length(edges_with_node) < abs(node_degchanges[i])) {
        edges_to_rm <- edges_with_node
      } else {
        edges_to_rm <- sample(x = edges_with_node, size = abs(node_degchanges[i]))
      }
      # class liftover
      edges_to_rm <- igraph::ends(graph_network, edges_to_rm)
      edges_to_rm <- apply(edges_to_rm, 1, function(x){paste(x, collapse = "|")})
      # make changes
      graph_network <- igraph::delete_edges(graph = graph_network,
                                            edges = edges_to_rm)
    } else {
      next
    } # end ifelse
  }

  #......................
  # out
  #......................
  return(graph_network)
}




#' @title Wrapper for `manip_degdist`
#' @noMd
#' @return igraph network
#' @export

wrapper_manip_degdist <- function(basenetpath, degprob, degvar) {
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
#' @export

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
#' @export

wrapper_manip_modular_rmedges <- function(basenetpath, edge_rm_num) {
  grphnet <- readRDS(basenetpath)
  out <- manip_modular_rmedges(graph_network = grphnet,
                               edge_rm_num = edge_rm_num)
  return(out)
}

#' @title Increase Network Unity by Adding Connected Edges
#' @inheritParams manip_degdist
#' @param edge_add_num int; Number of edges to add from current graph
#' @description To increase unity of a given network, a number of well-connected
#' edges that could exist in the saturated graph (as identified with \code{igraph::complementer`}
#' are added. Edge connectedness is determined by the \code{igraph::edge_betweenness}
#' function. The process is deterministic with the most connected edges being added
#' sequentially
#' @returns network graph (class igraph) with new modularity
#' @export

manip_unity_addedges <- function(graph_network, edge_add_num) {
  #......................
  # checks
  #......................
  goodegg::assert_eq("igraph", class(graph_network),
                     message = "The graph_network object must be have the igraph class (i.e. generate network with igraph)")
  goodegg::assert_single_int(edge_add_num)
  #......................
  # setup (const, storage, etc)
  #......................
  new_mod <- graph_network
  comp_graph_network <- igraph::complementer(graph_network)

  #......................
  # core
  #......................
  # identify most connected edges by betweeness from the complemented graph
  btwnconn <- igraph::edge_betweenness(comp_graph_network, directed = F)
  # add num edges specified based on most connected
  btwnord <- rev( order(btwnconn) )[1:edge_add_num]
  # because this is a different graph than original, we need to do perform additional unpacking/descriptive work
  new_edges_to_add <- igraph::ends(comp_graph_network, igraph::E(comp_graph_network)[btwnord])
  # add edges
  # modify graph in place (clearer/slightly safer than doing this w/ vector form in igraph)
  for (i in 1:nrow(new_edges_to_add)) {
    new_mod <- igraph::add_edges(graph = new_mod,
                                 edges = c(new_edges_to_add[i,1], # from
                                           new_edges_to_add[i,2])) # to
  }

  #......................
  # out
  #......................
  return(new_mod)
}




#' @title Wrapper for `manip_unity_addedges`
#' @noMd
#' @return igraph network
#' @export

wrapper_manip_unity_addedges <- function(basenetpath, edge_add_num) {
  grphnet <- readRDS(basenetpath)
  out <- manip_unity_addedges(graph_network = grphnet,
                              edge_add_num = edge_add_num)
  return(out)
}





#' @title Identify Potential Edges to Create Triangles (Clustering) in Network
#' @inheritParams manip_degdist
#' @description Identifies potential edges that would create triangles from
#' the input network. Returns those edges as an A(i,j) listing, where A is the adjacency
#' matrix
#' @details Based on premise of the cross product of an adjacency matrix with itself
#' will identify common neighbors, which can then be used to identify nodes
#' that could have potential triangular relationships. NB, this relationship is only
#' true for undirected, symmetric network (otherwise the transpose is needed)
#' @returns list of potential edges that would introduce new triangles in the graph
#' @export

get_potential_triangle_edges <- function(graph_network) {
  #......................
  # checks
  #......................
  goodegg::assert_eq("igraph", class(graph_network),
                     message = "The graph_network object must be have the igraph class (i.e. generate network with igraph)")

  #......................
  # setup
  #......................
  # squared adjacency matrix gives common neighbors and path lengths (Gilbert Strang, pg 78)
  gnet_adjmat <- igraph::as_adjacency_matrix(graph_network, sparse = F)
  gnet_adjmatsq <- gnet_adjmat %*% gnet_adjmat
  # drop out elements we don't need to avoid redundancies (matrix is symmetric)
  diag(gnet_adjmatsq) <- 0
  gnet_adjmatsq[upper.tri(gnet_adjmatsq)] <- 0

  #......................
  # core
  #......................
  # do nodes i,j share a common neighbor
  #   common neighbors mean potential for transitivity
  newijs <- which(gnet_adjmatsq > 0)
  i <- newijs %% nrow(gnet_adjmatsq) # identify ij based on R order vector from matrix
  i <- ifelse(i == 0, nrow(gnet_adjmatsq), i)
  j <- ceiling( newijs /ncol(gnet_adjmatsq) )
  newtri_conns <- mapply(function(x,y){return(c(x,y))},
                         x = i, y = j, SIMPLIFY = FALSE)

  #......................
  # out
  #......................
  return(newtri_conns)
}


#' @title Add or Remove Edges to Change Clustering (Triangles)
#' @inheritParams manip_degdist
#' @param new_transitivity_prob numeric; The probability that adjacent vertices of a
#' node are connected (i.e. the clustering coefficient)
#' @description Method for adding or removing edges to change the probability of transitivity,
#' or clustering coefficient, of the current graph
#' @details In cases were the desired probability of transitivity is greater than the
#' current probability of transitivity, a subset of dyad pairs that could have a new edge
#' introduced between them to induce a new triangle (cluster) are identified and created.
#' In instances where the desired probability of transitivity is lower than the current,
#' edges that make up triangles are removed.
#' @importFrom PFSwDemonDanProphetBrown get_potential_triangle_edges
#' @returns network graph (igraph class) with new clustering coefficient
#' @export

manip_clust_edges <- function(graph_network, new_transitivity_prob) {
  #......................
  # checks
  #......................
  goodegg::assert_eq("igraph", class(graph_network),
                     message = "The graph_network object must be have the igraph class (i.e. generate network with igraph)")
  goodegg::assert_single_numeric(new_transitivity_prob)

  #......................
  # setup (const, storage, etc)
  #......................
  new_graph_network <- graph_network
  start_trans <- igraph::transitivity(graph_network)
  #......................
  # core
  #......................
  # NB the number of triangles that we induce or remove will be not be a 1:1 correlation
  # with edges removed due to interdependencies
  # therefore, will loop through and calculate transitivity iteratively
  # there is a potential to run out of edges to sample from (either in adding or removing) but have not added catch for that

  if ( new_transitivity_prob > start_trans ) { # adding clusters
    # identify edges that we could potentially add
    pot_add_edges <- get_potential_triangle_edges(graph_network)

    while( new_transitivity_prob > igraph::transitivity(new_graph_network) ) { # transitivity being updated iteratively, so needs to be called iteratively
      #......................
      # sampling efficiency yields approx transitivity, not exact
      #......................
      diffprob <- new_transitivity_prob - igraph::transitivity(new_graph_network)
      nedge <- ceiling( diffprob * igraph::ecount(new_graph_network) )
      new_edge_index <-  sample(1:length(pot_add_edges), size = nedge)
      # pick edge to add and add it
      new_edges <- unlist(pot_add_edges[ new_edge_index ])
      new_graph_network <- igraph::add_edges(new_graph_network,
                                             edges = new_edges)
      # drop used edge
      pot_add_edges <- pot_add_edges[!(1:length(pot_add_edges) %in% new_edge_index)]
    } # end while
  } else if (new_transitivity_prob < start_trans ) { # removing clusters
    # identify edges that we could potentially remove
    pot_rm_edges <- igraph::triangles(graph_network) # per igraph docs: For efficiency, all triangles are returned in a single vector. The first three vertices belong to the first triangle, etc.
    # per igraph docs need to split this
    pot_rm_edges <- igraph::as_ids(pot_rm_edges)
    grps <- factor( sort(rep(1:(length(pot_rm_edges)/3), 3)) )
    pot_rm_edges <- split(pot_rm_edges, f = grps)

    while (new_transitivity_prob < igraph::transitivity(new_graph_network) ) {

      #......................
      # sampling efficiency yields approx transitivity, not exact
      # but given dampening effect, should be close...
      #......................
      diffprob <- abs( new_transitivity_prob - igraph::transitivity(new_graph_network) )
      nedge <- ceiling( diffprob * igraph::ecount(new_graph_network) )
      del_edge_index <- sample(1:length(pot_add_edges), size = nedge)

      # pick left to right vs right to left
      if(rbinom(1)) {
        pot_rm_edges[[del_edge_index]] <- paste(pot_rm_edges[[del_edge_index]][1:2], collapse = "|")
      } else {
        pot_rm_edges[[del_edge_index]] <- paste(pot_rm_edges[[del_edge_index]][2:3], collapse = "|")
      }
      # update graph
      del_edge <- unlist(pot_rm_edges[del_edge_index])
      new_graph_network <- igraph::delete_edges(new_graph_network,
                                                edges = del_edge) # r a character vector containing the IDs or names of the source and target vertices, separated by |
      # drop used edge
      pot_rm_edges <- pot_rm_edges[!(1:length(pot_rm_edges) %in% del_edge_index)]
    } # end while
  } else {
    next
  }

  #......................
  # out
  #......................
  return(new_graph_network)
}

#' @title Wrapper for `manip_clust_addedges`
#' @noMd
#' @return igraph network
#' @export

wrapper_manip_clust_edges <- function(basenetpath, new_transitivity_prob) {
  grphnet <- readRDS(basenetpath)
  out <- manip_clust_edges(graph_network = grphnet,
                           new_transitivity_prob = new_transitivity_prob)
  return(out)
}
