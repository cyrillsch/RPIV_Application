## Implementation of "Smooth test" by Delgado, Dominguez and Lavegne, 2006 in the setting of linear IV models
smooth_Test <- function(Y, X, C = NULL, Z,  h, B = 999){
  n <- length(Y)
  Xbar <- cbind(rep(1,n), X, C)
  Zbar <- cbind(Z, C)
  Xhat <- lm(Xbar ~ Zbar)$fitted.values
  beta <- lm(Y ~ -1 + Xhat)$coefficients
  psi <- Y - Xhat %*% beta
  D2 <- as.matrix(dist(Zbar))^2
  q <- ncol(Zbar)
  K <- (2 * pi)^(-q / 2) * exp(-D2 / (2 * h^2))
  diag(K) <- 0
  calc_stat <- function(psi0){
    Kpsi0 <- K %*% psi0
    Tn0 <- sum(psi0 * Kpsi0)/(h^q * n * (n-1))
    return(Tn0)
  }
  Tn <- calc_stat(psi)
  w_probs <- c((sqrt(5) + 1)/2/sqrt(5), (sqrt(5) - 1)/2/sqrt(5))
  w_points <- c(-(sqrt(5) -1)/2, (sqrt(5) + 1)/2)
  weights <- matrix(sample(w_points, n * B, replace = TRUE, prob = w_probs), nrow = n)
  TBoot <- apply(weights, 2, function(wrow){calc_stat(psi * wrow)})
  K2psi2 <- K^2 %*% psi^2
  Vn <- 2 * sum(psi^2 * K2psi2)/(h^q * n * (n-1))
  tn <- n * h^(q/2) * Tn / Vn^(0.5)
  tBoot <- n * h^(q/2) * TBoot / Vn^(0.5)
  p_asymp <- 1-pnorm(tn)
  p_boot <- (1 + sum(tBoot > tn))/(B + 1)
  return(c(p_asymp = p_asymp, p_boot = p_boot))
}



## implementation of approach from Antoine and Lavergne 2023

ICM_test <- function(Y, X, C = NULL, Z, B = 499){
  ## Influence of exogenous control variables is projected out according to paper
  if(is.null(C)){
    Y <- lm(Y~1)$residuals
    X <- lm(X~1)$residuals
    Z <- lm(Z~1)$residuals
  } else {
    Y <- lm(Y~C)$residuals
    X <- lm(X~C)$residuals
    Z <- lm(Z~C)$residuals
  }
  n <- length(Y)
  ## Authors "recommend to scale the exogenous instruments by a measure of dispersion, such
  ## as their ampirical standard deviation", although I am not sure if they are applying this in
  ## their simulations.
  Z <- apply(as.matrix(Z), 2, scale, center = FALSE)
  d <- NCOL(Z)
  ## first step: estimate conditional variance Omega(z).
  ## Paper uses Gaussian kernel with rule-of-thumb bandwidth in simulations.
  ## from scaling, all components of Z have standard deviation 1, hence
  ## Silverman's Rule of thumb is
  h <- (4/(d+2))^(1/(d+4))/n^(1/(d+4))
  D2 <- as.matrix(dist(Z)^2)
  K <- exp(-D2 / (2 * h^2))
  fhat <- as.numeric(K %*% rep(1, n) / (n * h))
  bigY <- cbind(Y, X)
  bigYhat  <-  K %*% bigY / (n * h) / fhat
  bigY_center <- bigY - bigYhat
  p <- NCOL(bigY)
  Omega_cond <- array(NA, dim = c(p, p, n))
  
  ## Omega_cond[, , j] should contain Omega(Z_j)
  for(j in 1:n){
    Omega_cond[, , j] <- t(bigY_center) %*% diag(K[j, ]) %*% bigY_center/ (n * h) / fhat[j]
  }
  Omega <- apply(Omega_cond, c(1,2), sum) / n
  
  ## Second step: calculate test statistic.
  ## We need the function w(). The authors only specify that w() is a triangle density
  ## At a different place, they specify that it is imposed that the squared integral of w() equals one. 
  ## But the normalization is just scaling and not important. More important is the choice of the support of w,
  ## but there is no indication about this. Also, we assume that a triangle density
  ## in multiple dimensions is just the product of triangle densities.
  ## Hence, we use the following choice (with support [-1,1])
  w <- function(z){
    prod(pmax(1-abs(z), 0))
  }
  W <- outer(1:n, 1:n, Vectorize(function(i,j){
    w(Z[i,] - Z[j,])/n
  }))
  
  ICM <- function(beta){
    b0 <- c(1, -beta)
    bigYb0 <- bigY %*% b0
    return(t(bigYb0) %*%  W %*% bigYb0 / t(b0) %*% Omega %*% b0)
  }
  VarS <- function(beta){
    b0 <- c(1, -beta)
    num <- sapply(1:n, function(i){t(b0) %*% Omega_cond[, , i] %*% b0})
    denom <- t(b0) %*% Omega %*% b0
    return(num / as.numeric(denom))
  }
  G0 <- matrix(rnorm(n * B), nrow = n)
  calc_p <- function(beta){
    G <- G0 * sqrt(VarS(beta))
    a <- rowSums((t(G) %*% W) * t(G))
    return((sum(a > as.numeric(ICM(beta))) + 1)/(B + 1))
  }
  return(calc_p)
}
