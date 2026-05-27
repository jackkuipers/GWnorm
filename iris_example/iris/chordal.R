library(igraph)
library(matrixcalc)
library(CholWishart)
library(complexplus)

## input a chordal graph, relabel it according to MCS,
## find cliques and separators
chordalMCS <- function(graph) {
  if (is_chordal(graph)$chordal) {
    ## relabel the vertices PEOMCS
    PEO_alpha <- as.vector(max_cardinality(graph)$alpha) ## $alpham1 : PEO from MCS
    graph <- permute(graph, PEO_alpha)
    graph$PEO <- order(PEO_alpha)
    
    n <- vcount(graph)
    cliques <- list()
    kdeg <- vector()
    for (i in 1:(n)) {
      neighbours <- as.vector(graph[[i]][[1]])
      cliques[[i]] <- c(i,intersect(neighbours,c(i:n)))
      kdeg[i] <- length(cliques[[i]]) - 1
    }
    graph$cliques <- cliques
    graph$kdeg <- kdeg
    
    ## find perfect clique sequence
    perfectCliques <- list(cliques[[1]])
    representitives <- c(1)
    for (i in 2:n) {
      if (kdeg[i] >= kdeg[i-1]) {
        perfectCliques <- append(perfectCliques, list(cliques[[i]]))
        representitives <- c(representitives, i)
      }
    }
    perfectCliques <- rev(perfectCliques)
    representitives <- rev(representitives)
    graph$perfectCliques <- perfectCliques
    graph$representitives <- representitives
    
    ## find the separators
    separators <- list()
    unions <- vector()
    if (length(representitives) > 1) {
      for (i in 1:(length(representitives) - 1)) {
        unions <- union(unions, perfectCliques[[i]])
        separators[[i]] <- sort(intersect(unions, perfectCliques[[i+1]]))
      }
    }
    graph$separators <- separators
    return(graph)
  } else {
    print("The input graph is not chordal.")
  }
}


## find the ratio of the normalising constant
ratioChordal <- function(graph, scaleMatrix) {
  100 * log(2) + lgamma(27) + lgamma(27) + lgamma(27) - lgamma(2)- lgamma(2)- lgamma(2)+ lgamma(26.5)  - lgamma(1.5) -27 * log(det(B[c(1,3),c(1,3)])) - 27 *log(det(B[c(4,2),c(4,2)])) -27 * log(det(B[c(1,4),c(1,4)]))  +26.5 * log(B[1,1]) +26.5 * log(B[4,4]) 
  
  scaleMatrix <- scaleMatrix[graph$PEO, graph$PEO]
  output <- 100 * log(2)
  
  ### gamma
  for (i in 1:4) {
    deg <- graph$kdeg[i]
    output <- output + lgamma(25.5 + (deg + 2) / 2) - lgamma(0.5 + (deg + 2) / 2)
  }
  
  ### determinant (cliques)
  for (i in 1:length(graph$perfectCliques)) {
    clique <- graph$perfectCliques[[i]]
    size <- length(clique)
    output <- output - (25.5 + (size + 1) / 2) * log(det(as.matrix(scaleMatrix[clique,clique])))
  }
  
  ### determinant (separators)
  if (length(graph$separators) > 0) {
    for (i in 1:length(graph$separators)) {
      if (length(graph$separators[[i]]) > 0) {
        separator <- graph$separators[[i]]
        size <- length(separator)
        output <- output + (25.5 + (size + 1) / 2) * log(det(as.matrix(scaleMatrix[separator,separator])))
      }
    }
  }
  
  return (exp(output))
}


