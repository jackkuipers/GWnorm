
library(igraph)
library(GWnorm)

for (p in c(10, 15, 20)) {

for (density_level in 1:3) {

for (beta in 10^(0:2)){

  n_samp_GW <- 5e5

  if (p == 10) { # try to match times
    n_samp_BD <- 1e6
  }
  if (p == 15) { # try to match times
    n_samp_BD <- 6e5
  }
  if (p == 20) { # try to match times
    n_samp_BD <- 3e5
    if (density_level > 1) {
      n_samp_BD <- 6e5
    } 
  }

n_loop <- 40

out_data <- NULL

if (!file.exists(paste0("./sim_results/g", p, "d", density_level, "beta", beta, ".rData"))) {

for (seed_number in 1:50){

set.seed(100+seed_number)
d <- 2
if (density_level == 1) {
  edge_prob <- d*2/(p-1) # for reasonable p at least?!
}
if (density_level == 2) {
  edge_prob <- 1/2 # middle density
}
if (density_level == 3) {
  edge_prob <- 1 - d*2/(p-1) # complement graph to low density
}

graph <- sample_gnp(p, edge_prob)

adj <- as.matrix(as_adjacency_matrix(graph, type = "upper"))
prec_mat <- BDgraph::rgwish(adj = adj)
sigma <- cov2cor(solve(prec_mat))
data <- BDgraph::rmvnorm(n = p^2/5, sigma = sigma) # need at least the dimension plus extras
D <- cor(data) # keep it in correlation form

start.time <- Sys.time()
BD_outs <- rep(NA, n_loop)
for(k in 1:n_loop){
 BD_outs[k] <- I_Gnorm(graph, beta, D, n_samp = n_samp_BD, useBDgraph = TRUE)
}
end.time <- Sys.time()
time.taken <- (end.time - start.time)/n_loop
BD_outs[which(BD_outs == -Inf)] <- NA # remove any -Inf errors
out_data_local <- data.frame(seed = seed_number, p = p, method = "BD", beta = beta, dl = density_level, n_samp = n_samp_BD, time = time.taken, 
                             mean = mean(BD_outs, na.rm = TRUE), sd = sd(BD_outs, na.rm = TRUE), na_count = sum(is.na(BD_outs)))
print(BD_outs)
print(out_data_local)
out_data <- rbind(out_data, out_data_local)

start.time <- Sys.time()
GW_outs <- rep(NA, n_loop)
for(k in 1:n_loop){
  GW_outs[k] <- I_Gnorm(graph, beta, D, n_samp = n_samp_GW)
}
end.time <- Sys.time()
time.taken <- (end.time - start.time)/n_loop
out_data_local <- data.frame(seed = seed_number, p = p, method = "GW", beta = beta, dl = density_level, n_samp = n_samp_GW, time = time.taken, 
                             mean = mean(GW_outs, na.rm = TRUE), sd = sd(GW_outs, na.rm = TRUE), na_count = sum(is.na(GW_outs)))
print(GW_outs)
print(out_data_local)
out_data <- rbind(out_data, out_data_local)

}

save(out_data, file = paste0("./sim_results/g", p, "d", density_level, "beta", beta, ".rData"))

}

}
  
}
  
}
