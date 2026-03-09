# R script to reproduce the figures in Section F of the supplement
################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)


## H0

load("SimulationResults/Cluster_H0_vary_s/row_45.RData")

method_levels <- c("RP_clu","RP_hom","RP_het","weakRP_clu","weakRP_hom","weakRP_het")
method_labels <- c("RP clust.","RP hom.","RP het.","weak RP clust.","weak RP hom.","weak RP het.")

line_types <- c(
  "RP clust." = "dotdash",
  "RP hom."   = "solid",
  "RP het."   = "dashed",
  "weak RP clust."   = "dotdash",
  "weak RP hom."   = "solid",
  "weak RP het." = "dashed"
)

line_colors <- c(
  "RP clust." = "blue",
  "RP hom."   = "blue",
  "RP het."   = "blue",
  "weak RP clust."   = "darkgreen",
  "weak RP hom."   = "darkgreen",
  "weak RP het." = "darkgreen"
)

n_vec <- c(100, 200, 400, 800)
s_vec <- seq(0, 1, by = 0.1)

nrep <- 1000

r025 <- qbinom(0.025, size = nrep, prob = 0.05)/nrep
r975 <- qbinom(0.975, size = nrep, prob = 0.05)/nrep

data_long <- data_H0_clust %>%
  pivot_longer(cols = all_of(method_levels),
               names_to = "method", values_to = "rejection_rate") %>%
  mutate(method = factor(method, levels = method_levels, labels = method_labels),
         panel = factor(paste0("n == ", n), levels = paste0("n == ", n_vec)))

p <- ggplot(data_long, aes(x = s, y = rejection_rate, linetype = method, color = method)) +
  geom_hline(yintercept = 0.05, color = "grey") +
  geom_hline(yintercept = r025, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = r975, linetype = "dashed", color = "grey") +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ panel, ncol = 4, labeller = label_parsed) +
  scale_linetype_manual(values = line_types) +
  scale_color_manual(values = line_colors) +
  scale_x_continuous(breaks = seq(0,1,by=0.2)) +
  labs(x = expression(s[clust]), y = "Rejection rate") +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave("Plots/Cluster_H0_vary_s.pdf", p, width = 8, height = 3.4)


## HA

load("SimulationResults/Cluster_HA_vary_s/row_45.RData")

violation_labs <- c(
  Z_squared = "z squared",
  sign_Z = "sign(z)",
  misspec_squared = "misspec. squared",
  misspec_sign = "misspec. sign"
)

data_long <- data_HA_clust %>%
  mutate(violation = factor(violation, levels = names(violation_labs), labels = unname(violation_labs))) %>%
  pivot_longer(cols = all_of(method_levels),
               names_to = "method", values_to = "rejection_rate") %>%
  mutate(method = factor(method, levels = method_levels, labels = method_labels))

p <- ggplot(data_long, aes(x = s, y = rejection_rate, linetype = method, color = method)) +
  geom_hline(yintercept = 0.05, color = "grey") +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ violation, ncol = 4, scales = "free") +
  scale_linetype_manual(values = line_types) +
  scale_color_manual(values = line_colors) +
  scale_x_continuous(breaks = seq(0,1,by=0.2)) +
  labs(x = expression(s[clust]), y = "Rejection rate") +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave("Plots/Cluster_HA_vary_s.pdf", p, width = 8, height = 3.4)
