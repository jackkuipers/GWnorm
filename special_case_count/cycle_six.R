is_6_cycle_complement <- function(graph) {
  # return false if the graph does not have 6 vertices
  if (vcount(graph) != 6) {
    return (FALSE)
  }
  
  # the graph is C6 if there are 6 edges, connected, all vertices have degree 2
  if (ecount(graph) == 6 && is_connected(graph) && sum(degree(graph) == 2) == 6) {
    return (TRUE)
  }
  
  # complement
  graph <- complementer(graph)
  return (ecount(graph) == 6 && is_connected(graph) && sum(degree(graph) == 2) == 6)
}

# C6 example 5 1 3 2 4 6
is_6_cycle_complement(make_graph(edges = c(1,3, 1,5, 2,3, 2,4, 5,6, 6,4), n = 6, directed = FALSE))
# C6 complement example
is_6_cycle_complement(make_graph(edges = c(1,2, 1,4, 1,6, 2,5, 2,6, 3,4, 3,5, 3,6, 4,5), n = 6, directed = FALSE))
# not C6 examples
is_6_cycle_complement(make_graph(edges = c(1,3, 1,5, 2,3, 2,4, 5,6, 6,4), n = 7, directed = FALSE))
is_6_cycle_complement(make_graph(edges = c(1,4, 1,5, 2,3, 2,4, 5,6, 6,4), n = 7, directed = FALSE))
