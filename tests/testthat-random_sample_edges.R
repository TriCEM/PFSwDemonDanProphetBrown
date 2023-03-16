## .................................................................................
## Purpose: Demonstrate randomly sampling of edges from base network behaves as expected
##
## Author: Nick Brazeau
##
## Date: 13 March, 2023
##
## Notes:
## .................................................................................
library(igraph)

#............................................................
# play net
#...........................................................
# make edge mat and network
edge_matrix <- matrix(c(1,2, 2,3, 3,4, 4,5, 5,1),
                      ncol = 2, byrow = TRUE)
gnet <- graph_from_edgelist(edge_matrix,
                            directed = F)
# add node names
V(gnet)$name <- LETTERS[1:5]

#............................................................
# accessing edges
#...........................................................
# edge IDs are 1-based and go from 1 to num edges
igraph::E(gnet) # all edges
igraph::ecount(gnet) # count of edges

# table of edges
igraph::E(gnet)[[]]

# table of vertices
igraph::V(gnet)[[]]

# get the vertix IDs of the second edge in the graph
E(gnet)[2]



# vertices are similar 1-based and go to num vertices but also have node names that we can store (from what I did above)
igraph::get.edges(gnet, 1) # return vertices ids between that edge
igraph::ends(gnet, 1, names = T) # return vertices names between specified edge
igraph::ends(gnet, 1, names = F) # return vertices ids between specified edge
# ^^ functions are equivalent
igraph::get.edgelist(gnet, names = F) # all edges between nodes
igraph::get.edgelist(gnet, names = T) # all edges between nodes


#............................................................
# adding and deleting edges
#...........................................................
# add AD edge
gnet_add <- add_edges(gnet,
                      edges = c(1,4))
plot(gnet_add,
     vertex.label = V(gnet)$name,
     vertex.label.dist = 3,
     vertex.label.font = 2)


# delete AB edge
gnet_del <- igraph::delete.edges(gnet, edges = 1)
plot(gnet_del,
     vertex.label = V(gnet)$name,
     vertex.label.dist = 3,
     vertex.label.font = 2)

# subset edges
# TODO




#............................................................
# identifying paths with the most betweeness or connectedness
#...........................................................
# make new obvious between net
edge_matrix <- matrix(c(1,2, 2,3, 2,4, 2,5, 5,6, 5,7, 5,8),
                      ncol = 2, byrow = TRUE)
gnet_betw <- graph_from_edgelist(edge_matrix,
                                 directed = F)
V(gnet_betw)$name <- LETTERS[1:8] # add node names

# viz
plot(gnet_betw,
     vertex.label = V(gnet_betw)$name,
     vertex.label.dist = 3,
     vertex.label.font = 2)

# identify betweenness (see if expected)
igraph::edge_betweenness(gnet_betw, directed = F)
E(gnet_betw)[[]]

#......................
# works as expected but ?unwanted? result for a graph like below
#......................
set.seed(24)
gran <- igraph::sample_gnp(10, 3 / 10)
plot(gran)
igraph::edge_betweenness(gran)
E(gran)[[]]
# 1--5 has the highest value, which is correct but not the
# wanted behavior because we don't want to end up with unconnected nodes
# TODO qualification for isolation?



#............................................................
# clustering
#...........................................................
# squares of adjacency matrices for off diagonal elements give the
# number of "walks" of length 2 (and A^K gives K walks - Gilbert Strang, pg 78)
# for triangles, this is a value of 1
# for potential dyads A-B B-C that could form new triangle A-C
# the value is 2
# note, these are distances, so triangle becomes symmetric
#
# Although the walk distance is correct, when the network is a perimeter
# or line, it correctly says that all nodes that are adjacent have a single
# path of 2, when for our purposes any edge could be added to those nodes
# in order to induce a triangle

#......................
# corner case
#......................
plot(gnet,
     vertex.label = V(gnet)$name,
     vertex.label.dist = 3,
     vertex.label.font = 2)
adj_mat <- igraph::as_adjacency_matrix(gnet, sparse = F)
adj_matsq <- adj_mat %*% adj_mat

#......................
# normal case
#......................
plot(gnet_add,
     vertex.label = V(gnet_add)$name,
     vertex.label.dist = 3,
     vertex.label.font = 2)
adj_mat <- igraph::as_adjacency_matrix(gnet_add, sparse = F)
adj_matsq <- adj_mat %*% adj_mat

