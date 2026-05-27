
#seednumber <- 101
#top_i <- 5

N <- nrow(data)
b <- 3

nits_base <- 5.5e5

rmses <- rep(NA, top_i)
times <- rep(NA, top_i)

for(ii in 1:top_i) {
  
  set.seed(seednumber)
  
  n_its <- 2^(ii-1)*nits_base

start.time <- Sys.time()

if (method_val == "BDmpl_bd") {
  BDgraph_out <- BDgraph::bdgraph.mpl(data, iter = n_its, burnin = n_its/4, algorithm = "bdmcmc", g.prior = 0.5)
} else {
  n_its <- n_its*3
  BDgraph_out <- BDgraph::bdgraph.mpl(data, iter = n_its, burnin = n_its/4, algorithm = "rjmcmc", g.prior = 0.5)
}

# total time taken 
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units = "secs")
time.taken

times[ii] <- time.taken

load("./iris/Iris_graph.Rdata")

G_av <- BDgraph_out$p_links

rmses[ii] <- sqrt(mean((weighted_graph_true[upper.tri(weighted_graph_true)] - G_av[upper.tri(G_av)])^2))

}

res <- data.frame(method = method_val, seed = seednumber, is = 1:top_i, time = times, rmse = rmses)


