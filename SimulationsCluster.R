source("SimulationSetupCluster.R")

nrep <- 1000
n_cores <- 70


## H0

cluster_size <- 4
s_vec <- seq(0, 1, by = 0.1)

n_vec <- c(100, 200, 400, 800)
n_IV <- 1
n_control <- 2
iv_strength <- 1
heteroskedastic <- FALSE


gamma <- 0

data_H0_clust <- data.frame(matrix(
  NA,
  ncol = 14,
  nrow = length(n_vec) * length(s_vec)
))
colnames(data_H0_clust) <- c(
  "n", "n_clusters", "cluster_size", "s",
  "n_control", "n_iv", "iv_strength", "heteroskedastic",
  "RP_clu", "RP_hom", "RP_het",
  "weakRP_clu", "weakRP_hom", "weakRP_het"
)

RNGkind("L'Ecuyer-CMRG")
set.seed(915)
row_count <- 1

for (n in n_vec) {
  n_clusters <- n/cluster_size
  for (s in s_vec) {
    data_H0_clust$n[row_count] <- n
    data_H0_clust$n_clusters[row_count] <- n_clusters
    data_H0_clust$cluster_size[row_count] <- cluster_size
    data_H0_clust$s[row_count] <- s
    
    data_H0_clust$n_control[row_count] <- n_control
    data_H0_clust$n_iv[row_count] <- n_IV
    data_H0_clust$iv_strength[row_count] <- iv_strength
    data_H0_clust$heteroskedastic[row_count] <- heteroskedastic
    
    rej_rates <- calc_rej_rates(
      nrep = nrep, n_cores = n_cores,
      n_clusters = n_clusters, cluster_size = cluster_size, s = s,
      n_control = n_control, n_IV = n_IV, iv_strength = iv_strength,
      heteroskedastic = heteroskedastic, gamma = gamma
    )
    data_H0_clust[row_count, c("RP_clu","RP_hom","RP_het","weakRP_clu","weakRP_hom","weakRP_het")] <- rej_rates
    row_count <- row_count + 1
    save(data_H0_clust, file = paste0("SimulationResults/Cluster_H0_vary_s/row_", row_count, ".RData"))
  }
}



## HA

# Design
cluster_size <- 4
n <- 300
n_clusters <- n / cluster_size


s_vec <- seq(0, 1, by = 0.1)

n_IV <- 1
n_control <- 2
iv_strength <- 1
heteroskedastic <- FALSE


violation_vec <- c("Z_squared", "sign_Z", "misspec_squared", "misspec_sign")
gamma_map <- c(Z_squared = 1, sign_Z = 3, misspec_squared = 4, misspec_sign = 8)


data_HA_clust <- data.frame(matrix(
  NA,
  ncol = 16,
  nrow = length(s_vec) * length(violation_vec)
))
colnames(data_HA_clust) <- c(
  "n", "n_clusters", "cluster_size", "s",
  "n_control", "n_iv", "iv_strength", "heteroskedastic",
  "violation", "gamma",
  "RP_clu", "RP_hom", "RP_het",
  "weakRP_clu", "weakRP_hom", "weakRP_het"
)

RNGkind("L'Ecuyer-CMRG")
set.seed(915)
row_count <- 1

for (violation in violation_vec) {
  gamma <- unname(gamma_map[violation])
  
  for (s in s_vec) {
    data_HA_clust$n[row_count] <- n
    data_HA_clust$n_clusters[row_count] <- n_clusters
    data_HA_clust$cluster_size[row_count] <- cluster_size
    data_HA_clust$s[row_count] <- s
    
    data_HA_clust$n_control[row_count] <- n_control
    data_HA_clust$n_iv[row_count] <- n_IV
    data_HA_clust$iv_strength[row_count] <- iv_strength
    data_HA_clust$heteroskedastic[row_count] <- heteroskedastic
    
    data_HA_clust$violation[row_count] <- violation
    data_HA_clust$gamma[row_count] <- gamma
    
    rej_rates <- calc_rej_rates(
      nrep = nrep, n_cores = n_cores,
      n_clusters = n_clusters, cluster_size = cluster_size, s = s,
      n_control = n_control, n_IV = n_IV, iv_strength = iv_strength,
      heteroskedastic = heteroskedastic,
      gamma = gamma, violation = violation
    )
    
    data_HA_clust[row_count, c("RP_clu","RP_hom","RP_het","weakRP_clu","weakRP_hom","weakRP_het")] <- rej_rates
    row_count <- row_count + 1
    
    save(data_HA_clust, file = paste0("SimulationResults/Cluster_HA_vary_s/row_", row_count, ".RData"))
  }
}