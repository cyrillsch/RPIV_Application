## R script that is sourced in SimulationsH0.R and SimulationsHA.R
##################################################################

# to install the RPIV package:
# devtools::install_github("cyrillsch/RPIV")
library(RPIV)

library(ivreg)
library(parallel)

source("AdditionalImplementations.R")


# n is the sample size
# n_control is the number of control variables
# n_IV is the number of instruments
# iv_strength is strength of IV
# heteroskedastic is boolean
# gamma denotes the amount of violation that is added to Y
# violation either "Z_squared", "sign_Z", "misspec_squared" or "misspec_sign"
# Instruments Z ~ N(0,1) and (optionally) controls C ~ N(0,1)
# are partially correlated by averaging overlapping columns and scaling by 1/sqrt(2).
# The endogenous regressor X is driven by the tanh of ascaled sum of instruments,
# iv_strength * tanh(rowSums(Z) / sqrt(n_IV)), plus an error correlated with the outcome error.
# Allows for heteroskedasticity and several forms of violations of the IV assumptions
generate_data <- function(n, n_control, n_IV, iv_strength, heteroskedastic, gamma, violation = "Z_squared"){
  Z <- matrix(rnorm(n * n_IV), ncol = n_IV)
  if(n_control > 0){
    C <- matrix(rnorm(n * n_control), ncol = n_control)
    Z[, 1:min(c(n_control, n_IV))] <- (Z[, 1:min(c(n_control, n_IV))] + C[, 1:min(c(n_control, n_IV))])/sqrt(2)
  } else {
    C <- NULL
  }
  H <- rnorm(n)
  del <- -H + 0.3 * rnorm(n)
  if(n_control > 0){
    X <- iv_strength * tanh(rowSums(Z)/sqrt(n_IV)) + 0.3 * C[, 1] + del
  } else {
    X <- iv_strength * tanh(rowSums(Z)/sqrt(n_IV)) + del
  }
  
  eps <- H + 0.3 * rnorm(n)
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
  return(list(X = X, Y = Y, Z = Z, C = C))
}




# Returns 8 p-values for the given dataset:
# - RPIV homoskedastic and heteroskedastic
# - J-test (if we are in the just-identified setting, we add Z^2)
# - Delgado et al.'s (2006) smooth test both with asymptotic and bootstrap p-value
# - weak_RPIV test, maximal p-value over a grid of betas (homoskedastic and heteroskedastic). We use the type = "fit" option.
# - Antoine and Lavergne's (2023) ICM test, maximal p-value over a grid of betas

get_p_values <- function(dat){
  # J-test p-value
  if(NCOL(dat$Z) == 1){
    if(is.null(dat$C)){
      mod_ivreg <- ivreg(dat$Y ~ dat$X | dat$Z + I(dat$Z^2))
    } else {
      mod_ivreg <- ivreg(dat$Y ~ dat$X + dat$C | dat$Z + I(dat$Z^2) + dat$C)
    }
  } else {
    if(is.null(dat$C)){
      mod_ivreg <- ivreg(dat$Y ~ dat$X | dat$Z)
    } else {
      mod_ivreg <- ivreg(dat$Y ~ dat$X + dat$C | dat$Z + dat$C)
    }
  }
  p_J <- summary(mod_ivreg)$diagnostics[2 + NCOL(dat$X), 4]
  # use less trees and less tuning to make the simulations faster
  mod_RP <- RPIV_test(Y = dat$Y, X = dat$X, C = dat$C, Z=dat$Z, 
                      variance_estimator = c("homoskedastic", "heteroskedastic"),
                      regr_par = list(num.trees = 200, num_min.node.size = 20))
  p_RP_hom <- mod_RP$homoskedastic$p_value
  p_RP_het <- mod_RP$heteroskedastic$p_value
  p_smooth <- unname(smooth_Test(Y = dat$Y, X = dat$X, C = dat$C, Z = dat$Z))
  p_smooth_asymp <- p_smooth[1]
  p_smooth_boot <- p_smooth[2]
  
  beta_grid <- seq(-10, 10, length.out = 200)
  calc_p_ICM <- ICM_test(Y = dat$Y, X = dat$X, C = dat$C, Z = dat$Z)
  ps_ICM <- sapply(beta_grid, calc_p_ICM)
  p_ICM <- max(ps_ICM)
  
  calc_p_weak_RPIV <- weak_RPIV_test(Y = dat$Y, X = dat$X, C = dat$C, Z = dat$Z, variance_estimator = c("homoskedastic", "heteroskedastic"),
                                     regr_par = list(num.trees = 200, num_min.node.size = 20))
  Ts_weak_RPIV <- sapply(beta_grid, calc_p_weak_RPIV, type = "fit")
  p_weakRP_hom <- max(1-pnorm(Ts_weak_RPIV[1,]))
  p_weakRP_het <- max(1-pnorm(Ts_weak_RPIV[2,]))
  
  return(c(p_RP_hom = p_RP_hom, p_RP_het = p_RP_het, p_J = p_J, p_smooth_asymp = p_smooth_asymp, p_smooth_boot = p_smooth_boot,
           p_weakRP_hom = p_weakRP_hom, p_weakRP_het = p_weakRP_het, p_ICM = p_ICM))
}


calc_rej_rates <- function(nrep, n_cores, n, n_control, n_IV, iv_strength, heteroskedastic, gamma, violation = "Z_squared", sig_level = 0.05){
  one_rep <- function(i){
    dat <- generate_data(n, n_control, n_IV, iv_strength, heteroskedastic, gamma, violation)
    return(get_p_values(dat))
  }
  list_rep <- mclapply(1:nrep, one_rep, mc.cores = n_cores)
  mat_rep <- do.call(rbind, list_rep)
  rej_rates <- apply(mat_rep, 2, function(v){mean(v<=sig_level)})
  return(rej_rates)
}


# Calculates rejection rates for weak_RPIV tests and Antoine and Lavergne's ICM test for a vector of candidate beta values
calc_rej_rates_weak <- function(nrep, n_cores, n, n_control, n_IV, iv_strength, heteroskedastic, gamma, violation = "Z_squared", sig_level = 0.05, beta_vec){
  one_rep <- function(i){
    dat <- generate_data(n, n_control, n_IV, iv_strength, heteroskedastic, gamma, violation)
    calc_p_weak_RPIV <- weak_RPIV_test(Y = dat$Y, X = dat$X, C = dat$C, Z = dat$Z, variance_estimator = c("homoskedastic", "heteroskedastic"),
                                       regr_par = list(num.trees = 200, num_min.node.size = 20))
    Ts_weak_RPIV <- sapply(beta_vec, calc_p_weak_RPIV, type = "fit")
    ps_weakRP_hom <- 1-pnorm(Ts_weak_RPIV[1,])
    ps_weakRP_het <- 1-pnorm(Ts_weak_RPIV[2,])
    calc_p_ICM <- ICM_test(Y = dat$Y, X = dat$X, C = dat$C, Z = dat$Z)
    ps_ICM <- sapply(beta_vec, calc_p_ICM)
    return(list(ps_weakRP_hom, ps_weakRP_het,ps_ICM))
  }
  list_rep <- mclapply(1:nrep, one_rep, mc.cores = n_cores)
  p_weakRP_hom <- do.call(rbind, lapply(list_rep, `[[`, 1))
  p_weakRP_het <- do.call(rbind, lapply(list_rep, `[[`, 2))
  p_ICM <- do.call(rbind, lapply(list_rep, `[[`, 3))
  c_rr <- function(v){mean(v <= sig_level)}
  rr_weakRP_hom <- apply(p_weakRP_hom, 2, c_rr)
  rr_weakRP_het <- apply(p_weakRP_het, 2, c_rr)
  rr_ICM <- apply(p_ICM, 2, c_rr)
  return(list(rr_weakRP_hom = rr_weakRP_hom, rr_weakRP_het = rr_weakRP_het, rr_ICM = rr_ICM, beta_vec = beta_vec))
}

