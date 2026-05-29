library(igraph)

source('prime_decom.R')
source('form_triangle.R')
source('cycle_six.R')
source('complete_partite.R')

# graph is an igraph object
# output is a vector of dimensions
# 0 means exact formula involving Gamma functions
# 1 means one dimensional integral
integrate_dim <- function(graph) {
  # if the input is an empty graph
  if (vcount(graph) == 0) {
    return (c())
  }
  
  # prime decomposition
  prime_components <- prime_decom(graph)$components
  dim <- rep(0, length(prime_components))
  
  # find the dimension for each prime component
  for (i in 1:length(prime_components)) {
    # get the induced subgraph of the prime component
    subgraph <- induced_subgraph(graph, prime_components[[i]], impl = c("auto"))
    
    # (C3) if the subgraph is C6 or its complement, then 3 dimensional integral
    if (is_6_cycle_complement(subgraph)) {
      dim[i] <- 1
      next
    }
    
    # (B3) if it is complete k-partite, then exact formula
    if (is_complete_k_partite(subgraph)) {
      next
    }
    
    # chordal information of the prime component
    chordal_subgraph <- is_chordal(subgraph, fillin = TRUE, newgraph = TRUE)
    
    # (B2) if it is chordal or has fill in 1 or 2, then exact formula
    if (chordal_subgraph$chordal || length(chordal_subgraph$fillin) <= 4) {
      next
    }
      
    # (C2) if three missing edges form a triangle, then 3 dimensional integral
    if (length(chordal_subgraph$fillin) == 6) {
      if (length(unique(chordal_subgraph$fillin)) == 3) {
        dim[i] <- 1
        next
      }
    }
        
    # (B1) if no two missing edges form a triangle in chordal, then exact formula
    chordal_graph <- chordal_subgraph$newgraph
    missing_edges_graph <- difference(chordal_graph, graph)
    if (!form_triangle(missing_edges_graph, chordal_graph)) {
      next
    }
    
    # for other subgraphs, we use Theorem 1
    dim[i] <- length(chordal_subgraph$fillin) / 2
  }
  
  return (dim)
}

integrate_dim(make_graph(edges = c(1,2, 2,4, 2,5, 1,3, 3,5, 5,6, 6,7, 5,7, 4,7), n = 7, directed = FALSE))

js <- 1:5
ns <- 10*js
ks <- 1:4
ds <- 0.25 + 0.25*ks

n_reps <- 4e4

frac_mat <- matrix(NA, length(ks), length(js))

set.seed(123)
for (j in js) {
  n <- ns[j]
  print(n)
  for (k in ks) {
    d <- ds[k]
    print(d)
    p <- d*2/(n-1) # for reasonable n?!

    fracs <- rep(NA, n_reps)
    for (i in 1:n_reps) {
      graph <- sample_gnp(n,p)
      #print(mean(degree(graph)))
      comps <- integrate_dim(graph)
      fracs[i] <- (sum(comps > 1) < 1)
    }
    frac_mat[k, j] <- mean(fracs)
  }
}

frac_mat

save(paste0("frac_mat_",n_reps,".Rdata"), frac_mat)
