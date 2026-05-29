#' Find triangle contains two missing edges # do not export
#'
#' This function checks whether there is a triangle in the graph(graph2)
#' contains at least two edges from the missing edges (graph1).
#'
#' @param graph1 An igraph object, corresponding to the missing edges
#' @param graph2 An igraph object, corresponding to the super graph
#' @return TRUE if there is such a triangle
#' @export
form_triangle <- function(graph1, graph2) {
  degree <- igraph::degree(graph1)

  p <- igraph::vcount(graph1)
  for (v in 1:p) {
    if (degree[v] >= 2) {
      d <- degree[v]
      neighbours <- igraph::neighbors(graph1, v) # at least 2 neighbours
      for (i in 1:(d-1)) {
        for (j in (i+1):d) {
          if (igraph::are_adjacent(graph2, i, j)) {
            return (TRUE)
          }
        }
      }
    }
  }

  return (FALSE)
}
