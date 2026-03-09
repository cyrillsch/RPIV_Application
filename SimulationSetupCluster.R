## Simulation setup under clustered dependence structure
################################################################################

# to install the RPIV package:
# devtools::install_github("cyrillsch/RPIV")
library(RPIV)
library(parallel)


# n_clusters is the number of clusters
# cluster_size is the size of clusters
# s is the amount of dependence within clusters, s = 0 is i.i.d. data, s = 1 means
# all observations in the same cluster are identical.


generate_data <- function(n_clusters, cluster_size, s, n_control, n_IV, iv_strength, heteroskedastic, gamma, violation = "Z_squared"){
  clustering <- rep(1:n_clusters, each = cluster_size)
  clnorm <- function(){
    return(sqrt(s) * rep(rnorm(n_clusters), each = cluster_size) + sqrt(1-s) * rnorm(n_clusters * cluster_size))
  }
  Z <- replicate(n_IV, clnorm())
  if(n_control > 0){
    C <- replicate(n_control, clnorm())
    Z[, 1:min(c(n_control, n_IV))] <- (Z[, 1:min(c(n_control, n_IV))] + C[, 1:min(c(n_control, n_IV))])/sqrt(2)
  } else {
    C <- NULL
  }
  H <- clnorm()
  del <- -H + 0.3 * clnorm()
  if(n_control > 0){
    X <- iv_strength * tanh(rowSums(Z)/sqrt(n_IV)) + 0.3 * C[, 1] + del
  } else {
    X <- iv_strength * tanh(rowSums(Z)/sqrt(n_IV)) + del
  }
  
  eps <- H + 0.3 * clnorm()
  if(heteroskedastic){
    eps <- eps * abs(Z[, 1])
  }
  if(n_control > 0){
    linear_part <- -X + 0.5 * rowSums(C)
  } else {
    linear_part <- -X
  }
  Y <- 2 + linear_part + eps
  
  if(violation == "Z_squared"){
    Y <-  Y + gamma * Z[, 1]^2
  }
  if(violation == "sign_Z"){
    Y <-  Y + gamma * sign(Z[, 1])
  }
  
  if(violation == "misspec_squared"){
    Y <-  Y + gamma * (linear_part)^2
  }
  if(violation == "misspec_sign"){
    Y <-  Y + gamma * sign(linear_part)
  }
  return(list(X = X, Y = Y, Z = Z, C = C, clustering = clustering))
}




# returns 6 p_values for a given dataset: RP clust./hom./het. and weak RP clust./hom./het.


get_p_values <- function(dat){
  mod_RP <- RPIV_test(Y = dat$Y, X = dat$X, C = dat$C, Z=dat$Z, 
                      variance_estimator = c("homoskedastic", "heteroskedastic", "cluster"),
                      clustering = dat$clustering,
                      regr_par = list(num.trees = 200, num_min.node.size = 15))
  p_RP_hom <- mod_RP$homoskedastic$p_value
  p_RP_het <- mod_RP$heteroskedastic$p_value
  p_RP_clu <- mod_RP$cluster$p_value
  
  beta_grid <- seq(-10, 10, length.out = 200)
  calc_p_weak_RPIV <- weak_RPIV_test(Y = dat$Y, X = dat$X, C = dat$C, Z = dat$Z, variance_estimator = c("homoskedastic", "heteroskedastic", "cluster"),
                                     clustering = dat$clustering,
                                     regr_par = list(num.trees = 200, num_min.node.size = 20))
  Ts_weak_RPIV <- sapply(beta_grid, calc_p_weak_RPIV, type = "fit")
  p_weakRP_hom <- max(1-pnorm(Ts_weak_RPIV[1,]))
  p_weakRP_het <- max(1-pnorm(Ts_weak_RPIV[2,]))
  p_weakRP_clu <- max(1-pnorm(Ts_weak_RPIV[3,]))
  return(c(p_RP_clu, p_RP_hom, p_RP_het, p_weakRP_clu, p_weakRP_hom, p_weakRP_het))
}

calc_rej_rates <- function(nrep, n_cores, n_clusters, cluster_size, s, n_control, n_IV, iv_strength, heteroskedastic, gamma, violation = "Z_squared", sig_level = 0.05){
  one_rep <- function(i){
    dat <- generate_data(n_clusters, cluster_size, s, n_control, n_IV, iv_strength, heteroskedastic, gamma, violation)
    return(get_p_values(dat))
  }
  list_rep <- mclapply(1:nrep, one_rep, mc.cores = n_cores)
  mat_rep <- do.call(rbind, list_rep)
  rej_rates <- apply(mat_rep, 2, function(v){mean(v<=sig_level)})
  return(rej_rates)
}


