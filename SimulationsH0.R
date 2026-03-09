# Simulations under the null hypothesis
#######################################

source("SimulationSetup.R")

nrep <- 1000
n_cores <- 70


###########################
# Simulation with varying n
###########################


n_vec <- c(50, 100, 150, 200, 250, 300, 400, 500)
n_control_vec <- c(0, 2)
n_IV_vec <- c(1, 2)
iv_strength_vec <- c(1)
heteroskedastic_vec <- c(FALSE, TRUE)

data_H0 <- data.frame(matrix(NA, ncol = 13, nrow = length(n_vec) * length(n_control_vec) * length(n_IV_vec) * length(iv_strength_vec) * length(heteroskedastic_vec)))
colnames(data_H0) <- c("n", "n_control", "n_iv", "iv_strength", "heteroskedastic", "RP_hom", "RP_het", "J", "smooth_asymp", "smooth_boot", "weakRP_hom", "weakRP_het", "ICM")



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
          rej_rates <- calc_rej_rates(nrep, n_cores, n, n_control, n_IV, iv_strength, heteroskedastic, gamma = 0, violation = "Z_squared", sig_level = 0.05)
          data_H0[row_count, 6:13] <- rej_rates
          row_count <- row_count + 1
        }
      }
    }
  }
  # save intermediate results
  save(data_H0, file = paste0("SimulationResults/H0_vary_n/row_", row_count, ".RData"))
}



#####################################
# Simulation with varying IV strength
#####################################


n_vec <- c(300)
n_control_vec <- c(0, 2)
n_IV_vec <- c(1, 2)
iv_strength_vec <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8, 1)
heteroskedastic_vec <- c(FALSE, TRUE)

data_H0 <- data.frame(matrix(NA, ncol = 13, nrow = length(n_vec) * length(n_control_vec) * length(n_IV_vec) * length(iv_strength_vec) * length(heteroskedastic_vec)))
colnames(data_H0) <- c("n", "n_control", "n_iv", "iv_strength", "heteroskedastic", "RP_hom", "RP_het", "J", "smooth_asymp", "smooth_boot", "weakRP_hom", "weakRP_het", "ICM")



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
          rej_rates <- calc_rej_rates(nrep, n_cores, n, n_control, n_IV, iv_strength, heteroskedastic, gamma = 0, violation = "Z_squared", sig_level = 0.05)
          data_H0[row_count, 6:13] <- rej_rates
          row_count <- row_count + 1
        }
        # save intermediate results
        save(data_H0, file = paste0("SimulationResults/H0_vary_iv_strength/row_", row_count, ".RData"))
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

data_H0 <- data.frame(matrix(NA, ncol = 13, nrow = length(n_vec) * length(n_control_vec) * length(n_IV_vec) * length(iv_strength_vec) * length(heteroskedastic_vec)))
colnames(data_H0) <- c("n", "n_control", "n_iv", "iv_strength", "heteroskedastic", "RP_hom", "RP_het", "J", "smooth_asymp", "smooth_boot", "weakRP_hom", "weakRP_het", "ICM")



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
          rej_rates <- calc_rej_rates(nrep, n_cores, n, n_control, n_IV, iv_strength, heteroskedastic, gamma = 0, violation = "Z_squared", sig_level = 0.05)
          data_H0[row_count, 6:13] <- rej_rates
          row_count <- row_count + 1
        }
      }
      # save intermediate results
      save(data_H0, file = paste0("SimulationResults/H0_vary_n_iv/row_", row_count, ".RData"))
    }
  }
}

