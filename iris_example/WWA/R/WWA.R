
#' @export
WWA <- function(data, n_its, burn_in = 0.25, b = 3) {
  
  df_0 <- b  # Degrees of freedom of the G-Wishart prior for the supergraph
  # data properties
  n <- nrow(data)
  p <- ncol(data)
  
  U <- stats::cov(data)*(nrow(data) - 1) # scatter matrix

  burnin_it <- burn_in*n_its
  
  G <- matrix(0L, nrow = p, ncol = p)
  G_av <- matrix(0, nrow = p, ncol = p)
  
  for (it in 1:n_its) {
    G <- update_G_WWA(G, 0.5, 3, U, n, 1)
    if (it > burnin_it) {
      G_av <- G_av + G
    }
  }
  G_av <- G_av/(n_its - burnin_it)
  G_av
}


#' MCMC update for a Gaussian graphical model
#' 
#' @description The MCMC update is a marginal update on the graph in the
#'  Gaussian graphical model with a G-Wishart prior. The rate matrix of the
#'  G-Wishart prior is chosen as the identity matrix. The edges are assumed to
#'  be independent a-priori with equal prior inclusion probability.
#' 
#' @param adj The symmetric adjacency matrix of the graph with zeroes on the
#'   diagonal.
#' @param edge_prob The prior edge inclusion probability.
#' @param df_0 The degrees of freedom of the G-Wishart prior.
#' @param U The scatter matrix of the observations.
#' @param n The number of observations constituting `U`.
#' @return The symmetric adjacency matrix of the graph after the MCMC update.
#' @export
update_G_WWA <- function (adj, edge_prob, df_0, U, n, n_edge) {
  p <- nrow(adj)
  
  return (tryCatch(
    expr = update_G_Rcpp(
      adj = adj,
      
      edge_prob_mat = matrix(edge_prob, nrow = p, ncol = p),
      df = df_0 + n,
      df_0 = df_0,
      rate = diag(p) + U,
      n_edge = n_edge,  # Number of single edge updates attempted in this MCMC update
      seed = sample.int(n = .Machine$integer.max, size = 1L),
      loc_bal = FALSE
    ), error = function (e) {
      warning(paste("`update_G_Rcpp` failed:", e))
      return (adj)
    }
  ))
}
