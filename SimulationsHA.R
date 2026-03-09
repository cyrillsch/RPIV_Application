# Simulations under the alternative hypothesis
##############################################


source("SimulationSetup.R")

nrep <- 1000
n_cores <- 70


############################################
# Simulation with varying violation_strength
############################################


n_vec <- c(300)
n_control_vec <- c(0, 2)
n_IV_vec <- c(1, 2)
iv_strength_vec <- c(1)
heteroskedastic_vec <- c(FALSE, TRUE)
gamma_mat <- cbind(seq(-1, 1, length.out = 25), seq(-3, 3, length.out = 25),
                   seq(-4, 4, length.out = 25), seq(-8, 8, length.out = 25))
violation_vec <- c("Z_squared", "sign_Z", "misspec_squared", "misspec_sign")

data_HA <- data.frame(matrix(NA, ncol = 15, 
                             nrow = length(n_vec) * length(n_control_vec) * length(n_IV_vec) * length(iv_strength_vec) * length(heteroskedastic_vec) * nrow(gamma_mat) * length(violation_vec)))
colnames(data_HA) <- c("n", "gamma", "violation", "n_control", "n_iv", "iv_strength", "heteroskedastic", "RP_hom", "RP_het", "J", "smooth_asymp", "smooth_boot", "weakRP_hom", "weakRP_het", "ICM")



RNGkind("L'Ecuyer-CMRG")
set.seed(915)
row_count <- 1

for(n in n_vec){
  for(n_control in n_control_vec){
    for(n_IV in n_IV_vec){
      for(iv_strength in iv_strength_vec){
        for(heteroskedastic in heteroskedastic_vec){
          for(l in 1:length(violation_vec)){
            violation <- violation_vec[l]
            for(k in 1:nrow(gamma_mat)){
              gamma <- gamma_mat[k, l]
              data_HA$gamma[row_count] <- gamma
              data_HA$violation[row_count] <- violation
              data_HA$n[row_count] <- n
              data_HA$n_control[row_count] <- n_control
              data_HA$n_iv[row_count] <- n_IV
              data_HA$iv_strength[row_count] <- iv_strength
              data_HA$heteroskedastic[row_count] <- heteroskedastic
              rej_rates <- calc_rej_rates(nrep, n_cores, n, n_control, n_IV, iv_strength, heteroskedastic, gamma = gamma, violation = violation, sig_level = 0.05)
              data_HA[row_count, 8:15] <- rej_rates
              row_count <- row_count + 1
            }
            # save intermediate results
            save(data_HA, file = paste0("SimulationResults/HA_vary_gamma/row_", row_count, ".RData"))
          }
        }
      }
    }
  }
}


#####################################
# Simulation with varying IV strength
#####################################


n_vec <- c(300)
n_control_vec <- c(0, 2)
n_IV_vec <- c(1, 2)
iv_strength_vec <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8, 1)
heteroskedastic_vec <- c(FALSE, TRUE)
gamma_mat <- cbind(c(1), c(3), c(4), c(8))
violation_vec <- c("Z_squared", "sign_Z", "misspec_squared", "misspec_sign")

data_HA <- data.frame(matrix(NA, ncol = 15, 
                             nrow = length(n_vec) * length(n_control_vec) * length(n_IV_vec) * length(iv_strength_vec) * length(heteroskedastic_vec) * nrow(gamma_mat) * length(violation_vec)))
colnames(data_HA) <- c("n", "gamma", "violation", "n_control", "n_iv", "iv_strength", "heteroskedastic", "RP_hom", "RP_het", "J", "smooth_asymp", "smooth_boot", "weakRP_hom", "weakRP_het", "ICM")



RNGkind("L'Ecuyer-CMRG")
set.seed(915)
row_count <- 1

for(n in n_vec){
  for(n_control in n_control_vec){
    for(n_IV in n_IV_vec){
      for(iv_strength in iv_strength_vec){
        for(heteroskedastic in heteroskedastic_vec){
          for(l in 1:length(violation_vec)){
            violation <- violation_vec[l]
            for(k in 1:nrow(gamma_mat)){
              gamma <- gamma_mat[k, l]
              data_HA$gamma[row_count] <- gamma
              data_HA$violation[row_count] <- violation
              data_HA$n[row_count] <- n
              data_HA$n_control[row_count] <- n_control
              data_HA$n_iv[row_count] <- n_IV
              data_HA$iv_strength[row_count] <- iv_strength
              data_HA$heteroskedastic[row_count] <- heteroskedastic
              rej_rates <- calc_rej_rates(nrep, n_cores, n, n_control, n_IV, iv_strength, heteroskedastic, gamma = gamma, violation = violation, sig_level = 0.05)
              data_HA[row_count, 8:15] <- rej_rates
              row_count <- row_count + 1
            }
          }
          # save intermediate results
          save(data_HA, file = paste0("SimulationResults/HA_vary_iv_strength/row_", row_count, ".RData"))
        }
      }
    }
  }
}


#######################################
# Simulation with varying number of IVs
#######################################


n_vec <- c(300)
n_control_vec <- c(0, 2)
n_IV_vec <- c(5, 10, 15, 20, 25)
iv_strength_vec <- c(1)
heteroskedastic_vec <- c(FALSE, TRUE)
gamma_mat <- cbind(c(1), c(3), c(4), c(8))
violation_vec <- c("Z_squared", "sign_Z", "misspec_squared", "misspec_sign")

data_HA <- data.frame(matrix(NA, ncol = 15, 
                             nrow = length(n_vec) * length(n_control_vec) * length(n_IV_vec) * length(iv_strength_vec) * length(heteroskedastic_vec) * nrow(gamma_mat) * length(violation_vec)))
colnames(data_HA) <- c("n", "gamma", "violation", "n_control", "n_iv", "iv_strength", "heteroskedastic", "RP_hom", "RP_het", "J", "smooth_asymp", "smooth_boot", "weakRP_hom", "weakRP_het", "ICM")



RNGkind("L'Ecuyer-CMRG")
set.seed(915)
row_count <- 1

for(n in n_vec){
  for(n_control in n_control_vec){
    for(n_IV in n_IV_vec){
      for(iv_strength in iv_strength_vec){
        for(heteroskedastic in heteroskedastic_vec){
          for(l in 1:length(violation_vec)){
            violation <- violation_vec[l]
            for(k in 1:nrow(gamma_mat)){
              gamma <- gamma_mat[k, l]
              data_HA$gamma[row_count] <- gamma
              data_HA$violation[row_count] <- violation
              data_HA$n[row_count] <- n
              data_HA$n_control[row_count] <- n_control
              data_HA$n_iv[row_count] <- n_IV
              data_HA$iv_strength[row_count] <- iv_strength
              data_HA$heteroskedastic[row_count] <- heteroskedastic
              rej_rates <- calc_rej_rates(nrep, n_cores, n, n_control, n_IV, iv_strength, heteroskedastic, gamma = gamma, violation = violation, sig_level = 0.05)
              data_HA[row_count, 8:15] <- rej_rates
              row_count <- row_count + 1
            }
          }
          # save intermediate results
          save(data_HA, file = paste0("SimulationResults/HA_vary_n_iv/row_", row_count, ".RData"))
        }
      }
    }
  }
}



