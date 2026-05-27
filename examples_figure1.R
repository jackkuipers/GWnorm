
library(GWnorm)

set.seed(42) # for BDgraph comparisons
beta <- 9 # corresponds to delta = 20


# Figure 1f: 5-cycle
p <- 5
adj <- matrix(0, nrow = p, ncol = p)
adj[1,5] <- adj[1,2] <- adj[2,3] <- adj[3,4] <- adj[4,5] <- 1 # 5-cycle
graph <- igraph::graph_from_adjacency_matrix(adj, mode = c("max"))
(IG_val <- I_G_special(graph, beta))
# convert to G_C
I_GtoC_G(graph, IG_val, beta)
# compare with BDgraph gnorm()
I_G_BD(graph, beta, diag(p), n_samp = 1e5, C_G_flag = TRUE)


# Figure 1d: 6-cycle
p <- 6
adj <- matrix(0, nrow = p, ncol = p)
adj[1,6] <- adj[1,2] <- adj[2,3] <- adj[3,4] <- adj[4,5] <- adj[5,6] <- 1 # 6-cycle
graph <- igraph::graph_from_adjacency_matrix(adj, mode = c("max"))
(IG_val <- I_G_special(graph, beta))
# convert to G_C
I_GtoC_G(graph, IG_val, beta)
# compare with BDgraph gnorm()
I_G_BD(graph, beta, diag(p), n_samp = 1e5, C_G_flag = TRUE)


# Figure 1a: 6-cycle complement
graph <- igraph::complementer(graph) # complement of 6-cycle
(IG_val <- I_G_special(graph, beta))
# convert to G_C
I_GtoC_G(graph, IG_val, beta)
# compare with BDgraph gnorm()
I_G_BD(graph, beta, diag(p), n_samp = 1e5, C_G_flag = TRUE)


# Figure 1c: nearly complete graph with three separate missing edges
adj[upper.tri(adj)] <- 1
adj[2,3] <- adj[1,4] <- adj[5,6] <- 0 # complete graph with some edges missing
graph <- igraph::graph_from_adjacency_matrix(adj, mode = c("max"))
(IG_val <- I_G_special(graph, beta))
# convert to G_C
I_GtoC_G(graph, IG_val, beta)
# compare with BDgraph gnorm()
I_G_BD(graph, beta, diag(p), n_samp = 1e5, C_G_flag = TRUE)


# Figure 1b: graph with three pairs of missing edges
p <- 8
adj <- matrix(0, nrow = p, ncol = p)
adj_add1 <- adj_add2 <- adj_add3 <- adj
adj[1, c(2:5,7)] <- adj[2, 4:8] <- adj[3, 5:8] <- adj[4, 5:6] <- 1
graph <- igraph::graph_from_adjacency_matrix(adj, mode = c("max"))
adj_add1[c(1,7), 8] <- 1 # blue missing edges
adj_add2[c(1,5), 6] <- 1 # red missing edges
adj_add3[2, 3] <- adj_add3[3, 4] <- 1 # dark blue missing edges
graph_star <- igraph::graph_from_adjacency_matrix(adj + adj_add1 + adj_add2 + adj_add3, 
                                                    mode = c("max"))
graph1 <- igraph::graph_from_adjacency_matrix(adj + adj_add2 + adj_add3, 
                                                  mode = c("max"))
graph2 <- igraph::graph_from_adjacency_matrix(adj + adj_add1 + adj_add3, 
                                              mode = c("max"))
graph3 <- igraph::graph_from_adjacency_matrix(adj + adj_add1 + adj_add2, 
                                              mode = c("max"))
IG1 <- I_G_special(graph1, beta)
IG2 <- I_G_special(graph2, beta)
IG3 <- I_G_special(graph3, beta)
IG_star <- I_G_special(graph_star, beta)
(IG_val <- IG1 + IG2 + IG3 - 2*IG_star)
# convert to G_C
I_GtoC_G(graph, IG_val, beta)
# compare with BDgraph gnorm()
I_G_BD(graph, beta, diag(p), n_samp = 1e5, C_G_flag = TRUE)


# Figure 1e:
p <- 7
adj <- matrix(0, nrow = p, ncol = p)
adj[1,7] <- adj[1,4] <- adj[2,4] <- adj[2,5] <- adj[5,6] <- adj[6,7] <- 1 # cycle around edge
adj[3,4] <- adj[3,6] <- adj[5,7] <- 1 # other edges
graph <- igraph::graph_from_adjacency_matrix(adj, mode = c("max"))
int_f <- function(t) {
  gsl::hyperg_U(1/2, 1, qchisq(t, 2*beta + 5)/2)^3
}
int_result <- stats::integrate(int_f, lower = 0, upper = 1)
(IG_val <- 3*log(pi)/2 + CholWishart::lmvgamma(beta + 5/2, 4) + 3*lgamma(beta + 2) + log(int_result$value))
# convert to G_C
I_GtoC_G(graph, IG_val, beta)
# compare with BDgraph gnorm()
I_G_BD(graph, beta, diag(p), n_samp = 1e5, C_G_flag = TRUE)


