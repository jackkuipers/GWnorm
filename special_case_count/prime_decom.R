#' Prime Decomposition # do not export
#'
#' This function returns the prime components and the minimal separators of a
#' connected graph.
#'
#' @param graph An igraph object, a connected graph
#' @examples
#' graph <- igraph::make_graph(edges = c(1,2, 2,4, 2,5, 1,3, 3,5, 5,6, 6,7, 5,7, 4,7), n = 7, directed = FALSE)
#' prime_decom(graph)
#' @export
prime_decom <- function(graph) {

  prime_components <- list()
  min_sep <- list()
  subgraph_stack <- igraph::groups(igraph::components(graph))

  while (length(subgraph_stack) >= 1) {
    found = FALSE # a flag

    # pop the last graph from stack
    curr_graph_v <- subgraph_stack[[length(subgraph_stack)]] # get
    subgraph_stack <- subgraph_stack[-length(subgraph_stack)] # remove
    curr_graph_v <- sort(curr_graph_v)
    curr_graph <- igraph::induced_subgraph(graph, curr_graph_v, impl = c("auto"))
    all_cliques <- igraph::cliques(curr_graph)
    # sort the cliques by length to make sure we look at minimal clique
    all_cliques <- all_cliques[order(sapply(all_cliques, length))]

    for (cli in all_cliques) {
      # obtain the induced subgraph by removing cli
      remaining_v <- setdiff(igraph::V(curr_graph), cli)
      removed <- igraph::induced_subgraph(curr_graph, remaining_v, impl = c("auto"))

      # get connected components of the induced subgraph
      clu <- igraph::components(removed)
      if (clu$no > 1) {
        found = TRUE
        connected_components <- igraph::groups(clu)

        # append subgraph_stack
        length(subgraph_stack) <- length(subgraph_stack) + clu$no
        for (i in 1:clu$no) {
          vertices <- union(remaining_v[connected_components[[i]]], cli)
          subgraph_stack[[length(subgraph_stack) - i + 1]] <- curr_graph_v[vertices]
        }

        # append min_sep
        length(min_sep) <- length(min_sep) + (clu$no - 1)
        for (i in 1:(clu$no - 1)) {
          min_sep[[length(min_sep) - i + 1]] <- curr_graph_v[cli]
        }

        # break the for loop
        break
      }
    }

    if (!found) {
      # it is a prime component
      prime_components <- c(prime_components, list(curr_graph_v))
    }
  }

  return (list(components = prime_components, separators = min_sep))
}
