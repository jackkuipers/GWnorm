
nits_base <- 7e4

rmses <- rep(NA, top_i)
times <- rep(NA, top_i)

for(ii in 1:top_i) {
  
  set.seed(seednumber)
  
  n_its <- 2^(ii-1)*nits_base

start.time <- Sys.time()  
  
  G_av <- WWA(data, n_its)
# time take for sampler
end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units = "secs")

times[ii] <- time.taken

load("./iris/Iris_graph.Rdata")

rmses[ii] <- sqrt(mean((weighted_graph_true[upper.tri(weighted_graph_true)] - G_av[upper.tri(G_av)])^2))

}

res <- data.frame(method = "WWA", seed = seednumber, is = 1:top_i, time = times, rmse = rmses)
