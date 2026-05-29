is_complete_k_partite <- function(graph) {

  while (vcount(graph) > 1) {
    # the vertex set
    vertices <- V(graph)
    # the first vertex
    v <- vertices[1]
    # the first vertex belongs to this parition
    partition <- c(v)
    # the neighbours of the first vertex
    neighbours_of_v <- neighbors(graph, v)
    for (i in 2:vcount(graph)) {
      # for each other vertex, either it is a neighbour of v
      # or it has the same neighbours as v
      u <- vertices[i]
      if (!are_adjacent(graph, v, u)) {
        if (setequal(neighbors(graph, u), neighbours_of_v)) {
          partition <- c(partition, u)
        } else {
          return (FALSE)
        }
      }
    }
    # remove the vertices from partition
    graph <- delete_vertices(graph, partition)
  }

  # after the while loop, return true
  return (TRUE)
}

# not complete k-partite example
is_complete_k_partite(make_graph(edges = c(1,2, 2,4, 2,5, 1,3, 3,5, 5,6, 6,7, 5,7, 4,7), n = 7, directed = FALSE))
# complete 4-partite example, (1,5), (2,4), (3,7), (6)
is_complete_k_partite(make_graph(edges = c(1,2, 1,3, 1,4, 1,6, 1,7, 2,3, 2,5, 2,6, 2,7, 3,4, 3,5, 3,6, 4,5, 4,6, 4,7, 5,6, 5,7, 6,7), n = 7, directed = FALSE))
