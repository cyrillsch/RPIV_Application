## R script to reproduce Figure 1 in Section 4.1.1.
###################################################

source("SimulationSetup.R")

nrep <- 1000
n_cores <- 18

n_vec <- c(25, 50, 75, 100, 150, 200, 250, 300, 400, 500)

data_H0 <- data.frame(matrix(NA, ncol = 6, nrow = length(n_vec) * 4))
colnames(data_H0) <- c("n", "overidentified", "heteroskedastic", "rr_RP_het", "rr_RP_hom", "rr_J")


RNGkind("L'Ecuyer-CMRG")
set.seed(915)
row_count <- 1
for(n in n_vec){
  for(overidentified in c(FALSE, TRUE)){
    for(heteroskedastic in c(FALSE, TRUE)){
      data_H0[row_count, 4:6] <- calc_rej_rates(nrep, n_cores, n, overidentified, heteroskedastic, gamma = 0, violation = "Z_squared", sig_level = 0.05)
      data_H0[row_count, 1] <- n
      data_H0[row_count, 2] <- overidentified
      data_H0[row_count ,3] <- heteroskedastic
      row_count <- row_count + 1
    }
  }
  # save intermediate results
  save(data_H0, file = paste("SimulationResults/H0_n_", n, ".RData", sep = ""))
}


# plot
library(ggplot2)
library(dplyr)
library(tidyr)

load("SimulationResults/H0_n_500.RData")

r025 <- qbinom(0.025, size = nrep, prob = 0.05)/nrep
r975 <- qbinom(0.975, size = nrep, prob = 0.05)/nrep

data_H0 <- data_H0 %>%
  mutate(
    id_type = ifelse(overidentified, "overidentified", "just-identified"),
    het_type = ifelse(heteroskedastic, "heteroskedastic", "homoskedastic"),
    panel = factor(paste(id_type, het_type, sep = ",\n"),
                   levels = c("just-identified,\nhomoskedastic",
                              "overidentified,\nhomoskedastic",
                              "just-identified,\nheteroskedastic",
                              "overidentified,\nheteroskedastic"))
  )

data_long <- data_H0 %>%
  select(n, rr_RP_het, rr_RP_hom, rr_J, panel) %>%
  pivot_longer(cols = c(rr_RP_het, rr_RP_hom, rr_J),
               names_to = "method",
               values_to = "rejection_rate") %>%
  mutate(method = factor(method,
                         levels = c("rr_RP_het", "rr_RP_hom", "rr_J"),
                         labels = c("RP Het.", "RP Hom.", "Overid. J")))

line_types <- c("RP Het." = "solid", "RP Hom." = "dashed", "Overid. J" = "dotted")

ggplot(data_long, aes(x = n, y = rejection_rate, linetype = method)) +
  geom_hline(yintercept = 0.05, color = "grey") +
  geom_hline(yintercept = r025, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = r975, linetype = "dashed", color = "grey") +
  geom_line(linewidth = 0.5) +
  scale_linetype_manual(values = line_types) +
  facet_wrap(~ panel, scales = "free_y") +
  labs(x = "n", y = "Rejection Rate", linetype = "Method") +
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

ggsave("Plots/H0.pdf", width = 5.5, height = 4)
