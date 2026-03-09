source("SimulationSetup.R")

nrep <- 1000
n_cores <- 70


############################################################################
# Simulation to assess how often H_0(beta) is rejected as a function of beta
############################################################################


n_vec <- c(100, 300)
n_control_vec <- c(0, 2)
n_IV_vec <- c(1, 2, 5)
iv_strength_vec <- c(0, 0.5, 1)
heteroskedastic_vec <- c(FALSE, TRUE)

n_beta <- 50

beta_vec <- seq(-4, 2, length.out = n_beta)

data_H0 <- data.frame(matrix(NA, ncol = 5 + 3 * length(beta_vec), nrow = length(n_vec) * length(n_control_vec) * length(n_IV_vec) * length(iv_strength_vec) * length(heteroskedastic_vec)))
colnames(data_H0) <- c("n", "n_control", "n_iv", "iv_strength", "heteroskedastic", paste0("weakRP_hom_beta", 1:n_beta), paste0("weakRP_het_beta", 1:n_beta), paste0("ICM_beta", 1:n_beta))



RNGkind("L'Ecuyer-CMRG")
set.seed(915)
row_count <- 1

for(n in n_vec){
  for(n_control in n_control_vec){
    for(n_IV in n_IV_vec){
      for(iv_strength in iv_strength_vec){
        for(heteroskedastic in heteroskedastic_vec){
          data_H0$n[row_count] <- n
          data_H0$n_control[row_count] <- n_control
          data_H0$n_iv[row_count] <- n_IV
          data_H0$iv_strength[row_count] <- iv_strength
          data_H0$heteroskedastic[row_count] <- heteroskedastic
          rej_rates <- calc_rej_rates_weak(nrep, n_cores, n, n_control, n_IV, iv_strength, heteroskedastic, gamma = 0, violation = "Z_squared", sig_level = 0.05, beta_vec = beta_vec)
          data_H0[row_count, 5 + 1:n_beta] <- rej_rates$rr_weakRP_hom
          data_H0[row_count, 5 + n_beta + 1:n_beta] <- rej_rates$rr_weakRP_het
          data_H0[row_count, 5 + 2*n_beta + 1:n_beta] <- rej_rates$rr_ICM
          row_count <- row_count + 1
        }
      }
      # save intermediate results
      save(data_H0, file = paste0("SimulationResults/WeakPower/row_", row_count, ".RData"))
    }
  }
}



