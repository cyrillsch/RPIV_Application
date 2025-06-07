## R script to reproduce Figures 2 and 3 in Section 4.1.2.
## and Figures 7 and 8 in Section C in the Supplementary Material
#################################################################

source("SimulationSetup.R")

nrep <- 1000
n_cores <- 30

n <- 400


gamma_mat <- cbind(seq(-0.75, 0.75, length.out = 45), seq(-1.5, 1.5, length.out = 45),
                   seq(-0.5, 0.5, length.out = 45), seq(-3, 3, length.out = 45))
violation_vec <- c("Z_squared", "sign_Z", "misspec_squared", "misspec_sign")



data_HA <- data.frame(matrix(NA, ncol = 7, nrow = NROW(gamma_mat) * length(violation_vec) * 4))
colnames(data_HA) <- c("gamma", "violation", "overidentified", "heteroskedastic", "rr_RP_het", "rr_RP_hom", "rr_J")


RNGkind("L'Ecuyer-CMRG")
set.seed(915)
row_count <- 1

for(overidentified in c(FALSE, TRUE)){
  for(heteroskedastic in c(FALSE, TRUE)){
    for(l in 1:length(violation_vec)){
      violation <- violation_vec[l]
      for(k in 1:NROW(gamma_mat)){
        gamma <- gamma_mat[k, l]
        data_HA[row_count, 5:7] <- calc_rej_rates(nrep, n_cores, n, overidentified, heteroskedastic, gamma = gamma, violation = violation, sig_level = 0.05)
        data_HA[row_count, 1] <- gamma
        data_HA[row_count, 2] <- violation
        data_HA[row_count, 3] <- overidentified
        data_HA[row_count ,4] <- heteroskedastic
        row_count <- row_count + 1
      }
      # save intermediate results
      save(data_HA, file = paste("SimulationResults/HA_", row_count, ".RData", sep = ""))
    }
      
  }
}



# plot
library(ggplot2)
library(dplyr)
library(tidyr)

load("SimulationResults/HA_721.RData")

data_HA$violation_label <- recode(data_HA$violation,
                                  "Z_squared" = "Z squared",
                                  "sign_Z" = "sign(Z)",
                                  "misspec_squared" = "misspec. squared",
                                  "misspec_sign" = "misspec. sign")

data_HA$violation_label <- factor(data_HA$violation_label,
                                  levels = c("Z squared", "sign(Z)", "misspec. squared", "misspec. sign"))
data_long <- data_HA %>%
    pivot_longer(
    cols = c(rr_RP_het, rr_RP_hom, rr_J),
    names_to = "method",
    values_to = "rejection_rate"
  ) %>%
  mutate(method = factor(method,
                         levels = c("rr_RP_het", "rr_RP_hom", "rr_J"),
                         labels = c("RP Het.", "RP Hom.", "Overid. J")))

line_types <- c("RP Het." = "solid", "RP Hom." = "dashed", "Overid. J" = "dotted")

for (overidentified in c(FALSE, TRUE)) {
  for (heteroskedastic in c(FALSE, TRUE)) {
    
    plot_title <- paste0(ifelse(overidentified, "overidentified, ", "just-identified, "),
                         ifelse(heteroskedastic, "heteroskedastic", "homoskedastic"))
    
    subset_data <- data_long %>%
      filter(overidentified == !!overidentified,
             heteroskedastic == !!heteroskedastic)
    
    p <- ggplot(subset_data, aes(x = gamma, y = rejection_rate, linetype = method)) +
      geom_hline(yintercept = 0.05, color = "grey50", linetype = "solid") +
      geom_line(size = 0.5) +
      facet_wrap(~ violation_label, nrow = 1, scales = "free_x") +  
      scale_linetype_manual(values = line_types) +
      labs(x = "Violation Strength", y = "Rejection Rate", linetype = "Method") +
      theme_bw() +
      theme(
        strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.spacing = unit(1, "lines"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 10)
      ) +
      ggtitle(plot_title)
    
    outname <- paste0("Plots/HA_", ifelse(overidentified, "overidentified_", "justidentified_"),
                      ifelse(heteroskedastic, "heteroskedastic", "homoskedastic"), ".pdf")
    ggsave(outname, plot = p, width = 6, height = 3)
  }
}

