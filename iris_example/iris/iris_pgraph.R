
source("H2m_fn.R")
source("chordal.R")

data(iris)
p <- 4
data <- iris[101:150, 1:p]
U <- cov(data)*(nrow(data) - 1)
b <- nrow(data) + p-1
beta <- (b + p-1) / 2
B <- U + diag(p)
A <- cov2cor(B)
options(digits = 10)

#############
# all adjacency matrices on 4 vertices
allAdjMat <- list()

edgeCombinations <- as.matrix(expand.grid(replicate(6, c(0, 1), simplify = FALSE)))
vectors <- matrix(0, nrow = 64, ncol = 16)
vectors[,2:4] = edgeCombinations[,1:3]
vectors[,7:8] = edgeCombinations[,4:5]
vectors[,12] = edgeCombinations[,6]

for (i in 1:64) {
  allAdjMat[[i]] <- t(matrix(vectors[i,], nrow = 4))
}
############

############
## get all the graphs
graphs <- list()
chordal <- vector()
for (i in 1:64) {
  graph <- graph_from_adjacency_matrix(allAdjMat[[i]], mode = c("undirected"))
  V(graph)$name <- c("1SL", "2SW", "3PL", "4PW")
  chordal[i] <- is.chordal(graph)$chordal
  graphs[[i]] <- graph
}

cycle4Indices <- which(chordal == FALSE) # 31, 46, 52

# find cliques and separators
graphsMCS <- list()
for (i in 1:64) {
  graphsMCS[[i]] <- chordalMCS(graphs[[i]])
}

#####

#####
## get the ratios for the 61 chordal graphs
ratios <- vector()
for (i in 1:64) {
  if (i == 31 | i == 46 | i == 52) {
    ratios[i] <- 0
  }
  else {
    ratios[i] <- ratioChordal(graphsMCS[[i]], B)
  }
}

### G31
ratios[31] <- exp(H2m_fn(beta, ds = matrix(c(A[2,3:4], A[1,4:3]), ncol = 2, byrow=TRUE)) + lgamma(28) + lgamma(27.5) + lgamma(27) + lgamma(26.5) - log(pi) + 100 * log(2) - 27.5 * log(prod(diag(B))) - 3 * lgamma(2.5) - lgamma(2) - lgamma(1.5) + lgamma(3))

### G46
ratios[46] <- exp(H2m_fn(beta, ds = matrix(c(A[3,4], A[2,3], A[1,2], A[1,4]), ncol = 2, byrow=TRUE)) + lgamma(28) + lgamma(27.5) + lgamma(27) + lgamma(26.5) - log(pi) + 100 * log(2) - 27.5 * log(prod(diag(B))) - 3 * lgamma(2.5) - lgamma(2) - lgamma(1.5) + lgamma(3))

### G52
ratios[52] <- exp(H2m_fn(beta, ds = matrix(c(A[1,2:3], A[3:2,4]), ncol = 2, byrow=TRUE)) + lgamma(28) + lgamma(27.5) + lgamma(27) + lgamma(26.5) - log(pi) + 100 * log(2) - 27.5 * log(prod(diag(B))) - 3 * lgamma(2.5) - lgamma(2) - lgamma(1.5) + lgamma(3))

prob <- ratios / sum(ratios)

weighted_graph <- matrix(0, p, p)
for (i in 1:length(prob)) {
  G <- allAdjMat[[i]]
  weighted_graph <- weighted_graph + prob[i]*G
}

weighted_graph_true <- weighted_graph

save(weighted_graph_true, file = "Iris_graph.Rdata")

