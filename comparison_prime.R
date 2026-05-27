
library(igraph)
library(GWnorm)

#for (converge_flag in c(TRUE, FALSE)) {
for (converge_flag in c(TRUE)) {

if (converge_flag) {
  p_vec <- c(10)
} else {
  p_vec <- c(10, 15, 20)
}

for (p in p_vec) {

if (converge_flag) {
# sample size for GWnorm, BDgraph gets twice that
  n_loop <- 100
  n_seeds <- 40
  n_samp_vec_GW <- 2^c(0:5)*1e4
  n_samp_vec_BD <- 2*n_samp_vec_GW
} else { # for different dimensions
  n_loop <- 40
  n_seeds <- 50
  n_samp_vec_GW <- 5e5
  if (p == 10) { # try to match times
    n_samp_vec_BD <- 1e6
  }
  if (p == 15) { # try to match times
    n_samp_vec_BD <- 6e5
  }
  if (p == 20) { # try to match times
    n_samp_vec_BD <- 3e5
  }
}

out_data <- NULL

for (seed_number in 1:n_seeds){

set.seed(100+seed_number)
d <- 2 # density of 2 gives good chance of connected prime graph
edge_prob <- d*2/(p-1) # for reasonable p at least?!

counter <- 0
good_graph <- FALSE
while (!good_graph) {
  graph <- sample_gnp(p, edge_prob)
  good_graph <- check_prime_connected(graph)
  counter <- counter + 1
}

beta <- 1
adj <- as.matrix(as_adjacency_matrix(graph, type = "upper"))
prec_mat <- BDgraph::rgwish(adj = adj)
sigma <- cov2cor(solve(prec_mat))
data <- BDgraph::rmvnorm(n = p^2/5, sigma = sigma) # need at least the dimension plus extras
D <- cor(data) # keep it in correlation form

tau <- I_Gnorm(graph, beta, D, err_flag = TRUE)$tau

for (ii in 1:length(n_samp_vec_GW)) {

start.time <- Sys.time()
BD_outs <- rep(NA, n_loop)
for(k in 1:n_loop){
 BD_outs[k] <- I_G_BD(graph, beta, D, n_samp = n_samp_vec_BD[ii])
}
end.time <- Sys.time()
time.taken <- (end.time - start.time)/n_loop
BD_outs[which(BD_outs == -Inf)] <- NA # remove any -Inf errors
out_data_local <- data.frame(seed = seed_number, p = p, method = "BD", tau = tau, n_samp = n_samp_vec_BD[ii], time = time.taken, 
                             mean = mean(BD_outs, na.rm = TRUE), sd = sd(BD_outs, na.rm = TRUE), na_count = sum(is.na(BD_outs)))
print(BD_outs)
print(out_data_local)
out_data <- rbind(out_data, out_data_local)

start.time <- Sys.time()
GW_outs <- rep(NA, n_loop)
for(k in 1:n_loop){
  GW_outs[k] <- I_G_MC(graph, beta, D, n_samp = n_samp_vec_GW[ii])
}
end.time <- Sys.time()
time.taken <- (end.time - start.time)/n_loop
out_data_local <- data.frame(seed = seed_number, p = p, method = "GW", tau = tau, n_samp = n_samp_vec_GW[ii], time = time.taken, 
                             mean = mean(GW_outs, na.rm = TRUE), sd = sd(GW_outs, na.rm = TRUE), na_count = sum(is.na(GW_outs)))
print(GW_outs)
print(out_data_local)
out_data <- rbind(out_data, out_data_local)

}

}

save(out_data, file = paste0("./sim_results/g", p, "prime", length(n_samp_vec_GW), ".rData"))

}
}
