## R script to reproduce Figure 4 in Section A.2.1. in the Supplementary Material
#################################################################################

source("SimulationSetupCluster.R")

nrep <- 1000
n_cores <- 18

n_vec <- c(100, 200, 400, 800)
cluster_size <- 4
n_clusters_vec <- n_vec/cluster_size
s_vec <- seq(0, 1, length.out = 11)

data_H0 <- data.frame(matrix(NA, ncol = 10, nrow = length(n_vec) * length(s_vec)))
colnames(data_H0) <- c("n", "overidentified", "heteroskedastic", "rr_RP_clu", "rr_RP_het", "rr_RP_hom", "rr_J", "s", "cluster_size", "n_clusters")

heteroskedastic <- FALSE
overidentified <- FALSE


RNGkind("L'Ecuyer-CMRG")
set.seed(915)
row_count <- 1
for(n_clusters in n_clusters_vec){
  for(s in s_vec){
    n <- n_clusters * cluster_size
    data_H0[row_count, 4:7] <- calc_rej_rates(nrep, n_cores, n_clusters, cluster_size, s, overidentified, heteroskedastic, gamma = 0, violation = "Z_squared", sig_level = 0.05)
    data_H0[row_count, 1] <- n
    data_H0[row_count, 2] <- overidentified
    data_H0[row_count, 3] <- heteroskedastic
    data_H0[row_count, 8] <- s
    data_H0[row_count, 9] <- cluster_size
    data_H0[row_count, 10] <- n_clusters
    row_count <- row_count + 1
  }
  # save intermediate results
  save(data_H0, file = paste("SimulationResultsCluster/H0_n_", n, ".RData", sep = ""))
}


# plot
library(ggplot2)
library(dplyr)
library(tidyr)

load("SimulationResultsCluster/H0_n_800.RData")

r025 <- qbinom(0.025, size = nrep, prob = 0.05)/nrep
r975 <- qbinom(0.975, size = nrep, prob = 0.05)/nrep

data_H0 <- data_H0 %>%
  mutate(
    panel = factor(n, labels = paste0("n = ", n_vec))
  )

data_long <- data_H0 %>%
  select(s, rr_RP_clu, rr_RP_het, rr_RP_hom, rr_J, panel) %>%
  pivot_longer(cols = c(rr_RP_clu, rr_RP_het, rr_RP_hom, rr_J),
               names_to = "method",
               values_to = "rejection_rate") %>%
  mutate(method = factor(method,
                         levels = c("rr_RP_clu", "rr_RP_het", "rr_RP_hom", "rr_J"),
                         labels = c("RP Cluster", "RP Het.", "RP Hom.", "Overid. J")))

line_types <- c("RP Cluster" = "dotdash", "RP Het." = "solid", "RP Hom." = "dashed", "Overid. J" = "dotted")

ggplot(data_long, aes(x = s, y = rejection_rate, linetype = method)) +
  geom_hline(yintercept = 0.05, color = "grey") +
  geom_hline(yintercept = r025, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = r975, linetype = "dashed", color = "grey") +
  geom_line(linewidth = 0.5) +
  scale_linetype_manual(values = line_types) +
  facet_wrap(~ panel, scales = "free_y") +
  labs(x = "s", y = "Rejection Rate", linetype = "Method") +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)
  )

ggsave("Plots/H0_Cluster.pdf", width = 5.5, height = 4)
