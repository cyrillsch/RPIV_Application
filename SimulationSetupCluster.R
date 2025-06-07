## R script that is sourced in SimulationsH0Cluster.R and SimulationsHACluster.R
################################################################################

# to install the RPIV package:
# devtools::install_github("cyrillsch/RPIV")
library(RPIV)
library(ivreg)
library(parallel)


# n_clusters is the number of clusters
# cluster_size is the size of clusters
# s is the amount of dependence within clusters, s = 0 is i.i.d. data, s = 1 means
# all observations in the same cluster are identical.

generate_data <- function(n_clusters, cluster_size, s, overidentified, heteroskedastic, gamma, violation = "Z_squared"){
  clustering <- rep(1:n_clusters, each = cluster_size)
  clnorm <- function(){
    return(s * rep(rnorm(n_clusters), each = cluster_size) + (1-s) * rnorm(n_clusters * cluster_size))
  }
  C <- cbind(clnorm(), clnorm())
  Z <- 0.5 * C[, 1] - 0.5 * C[, 2] +  clnorm()
  if(overidentified){
    Z <- cbind(Z, clnorm())
  }
  Z <- as.matrix(Z)
  H <- clnorm()
  del <- sign(H) + 0.5 * clnorm()
  X <- Z[, 1] -  0.5 * C[, 1] + del
  if(overidentified){
    X <- X - 0.5 * Z[, 2]
  }
  eps <- H + 0.5 * clnorm()
  if(heteroskedastic){
    eps <- eps * Z[, 1]^2
  }
  Y <- 2 - X - 0.5 * C[, 1] + C[, 2] + eps
  if(violation == "Z_squared"){
    Y <-  Y + gamma * Z[, 1]^2
  }
  if(violation == "sign_Z"){
    Y <-  Y + gamma * sign(Z[, 1])
  }
  if(violation == "misspec_squared"){
    Y <-  Y + gamma * (- X - 0.5 * C[, 1] + C[, 2])^2
  }
  if(violation == "misspec_sign"){
    Y <-  Y + gamma * sign(- X - 0.5 * C[, 1] + C[, 2])
  }
  return(list(X = X, Y = Y, Z = Z, C = C, clustering = clustering))
}


# returns 4 p_values for a given dataset: cluster, heteroskedastic, homoskedastic and J.
# If we are in just-identified setting, we add Z^2 as instrument for the J-test

get_p_values <- function(dat){
  # J-test p-value
  if(NCOL(dat$Z) == NCOL(dat$X)){
    mod_ivreg <- ivreg(Y ~ X + C | Z + I(Z^2) + C, data = dat)
  } else {
    mod_ivreg <- ivreg(Y ~ X + C | Z + C, data = dat)
  }
  p_J <- summary(mod_ivreg)$diagnostics[2 + NCOL(dat$X), 4]
  # use less trees and less tuning to make the simulations faster
  mod_RP <- RPIV_test(Y = dat$Y, X = dat$X, C = dat$C, Z=dat$Z, 
                      variance_estimator = c("homoskedastic", "heteroskedastic", "cluster"),
                      clustering = dat$clustering,
                      regr_par = list(num.trees = 200, num_min.node.size = 15))
  p_RP_hom <- mod_RP$homoskedastic$p_value
  p_RP_het <- mod_RP$heteroskedastic$p_value
  p_RP_clu <- mod_RP$cluster$p_value
  return(c(p_RP_clu, p_RP_het, p_RP_hom, p_J))
}

calc_rej_rates <- function(nrep, n_cores, n_clusters, cluster_size, s, overidentified, heteroskedastic, gamma, violation, sig_level){
  one_rep <- function(i){
    dat <- generate_data(n_clusters, cluster_size, s, overidentified, heteroskedastic, gamma, violation)
    return(get_p_values(dat))
  }
  list_rep <- mclapply(1:nrep, one_rep, mc.cores = n_cores)
  mat_rep <- do.call(rbind, list_rep)
  rej_rates <- apply(mat_rep, 2, function(v){mean(v<=sig_level)})
  return(rej_rates)
}


