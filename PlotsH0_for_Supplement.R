# R script to reproduce the figures in Section D.3 of the supplement
################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)


nrep <- 1000
r025 <- qbinom(0.025, size = nrep, prob = 0.05)/nrep
r975 <- qbinom(0.975, size = nrep, prob = 0.05)/nrep


# sim is either "vary_n", "vary_iv_strength" or "vary_n_iv"
plot_results <- function(sim){
  filename <- switch(sim,
                     vary_n = "SimulationResults/H0_vary_n/row_65.RData",
                     vary_iv_strength = "SimulationResults/H0_vary_iv_strength/row_57.RData",
                     vary_n_iv = "SimulationResults/H0_vary_n_iv/row_21.RData")
  load(filename)
  if(sim %in% c("vary_n", "vary_iv_strength")){
    data_H0 <- data_H0 %>%
      mutate(
        het_type = ifelse(heteroskedastic, "heteroskedastic", "homoskedastic"),
        panel = factor(paste0("n_iv_", n_iv, "_n_control_", n_control, "_het_type_", het_type))
      )
  }
  if(sim == "vary_n_iv"){
    data_H0 <- data_H0 %>%
      mutate(
        het_type = ifelse(heteroskedastic, "heteroskedastic", "homoskedastic"),
        panel = factor(paste0("n_control_", n_control, "_het_type_", het_type))
      )
  }
  data_H0 <- data_H0 %>%
    mutate(
      panel = factor(
        panel,
        levels = c(
          sort(grep("homoskedastic", unique(panel), value = TRUE)),
          sort(grep("heteroskedastic", unique(panel), value = TRUE))
        )
      )
    )
  
  data_long <- data_H0 %>%
    select(n, n_iv, iv_strength, RP_hom, RP_het, J, smooth_asymp, smooth_boot, weakRP_hom, weakRP_het, ICM, panel) %>%
    pivot_longer(cols = c(RP_hom, RP_het, J, smooth_asymp, smooth_boot, weakRP_hom, weakRP_het, ICM),
                 names_to = "method",
                 values_to = "rejection_rate") %>%
    mutate(method = factor(method,
                           levels = c("RP_hom", "RP_het", "weakRP_hom", "weakRP_het", "J", "ICM", "smooth_asymp", "smooth_boot"),
                           labels = c("RP hom.", "RP het.", "weak RP hom.", "weak RP het.", "overid. J", "ICM", "smooth asymp.", "smooth boot.")))
  
  line_types <- c("RP hom." = "dashed", "RP het." = "solid", "overid. J" = "dotted", "smooth asymp." = "dotted", "smooth boot." = "dotted", "weak RP hom." = "dashed", "weak RP het." = "solid", "ICM" = "dotted")
  line_colors <- c("RP hom." = "blue", "RP het." = "blue", "overid. J" = "orange", "smooth asymp." = "purple", "smooth boot." = "magenta", "weak RP hom." = "darkgreen", "weak RP het." = "darkgreen", "ICM" = "red")
  x_var <- switch(sim,
                  vary_n = "n",
                  vary_iv_strength = "iv_strength",
                  vary_n_iv = "n_iv")
  gsub1 <- switch(sim,
                  vary_n = "n_iv_(\\d+)_n_control_(\\d+)_het_type_(.*)",
                  vary_iv_strength = "n_iv_(\\d+)_n_control_(\\d+)_het_type_(.*)",
                  vary_n_iv = "n_control_(\\d+)_het_type_(.*)")
  gsub2 <- switch(sim,
                  vary_n =
                    "atop(n[IV] == \\1*','~n[c] == \\2, \\3)",
                  vary_iv_strength =
                    "atop(n[IV] == \\1*','~n[c] == \\2, \\3)",
                  vary_n_iv =
                    "atop(n[c] == \\1, \\2)"
  )
  
  x_breaks <- switch(sim,
                     vary_n = c(100, 200, 300, 400, 500),
                     vary_iv_strength = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                     vary_n_iv = c(5, 10, 15, 20, 25))
  x_lab <- switch(sim,
                  vary_n = "n",
                  vary_iv_strength = expression("IV strength " * pi),
                  vary_n_iv = expression(n[IV]))
  title <- switch(sim, 
                  vary_n = expression("Rejection rates under " * H[0] * ", varying n"),
                  vary_iv_strength = expression("Rejection rates under " * H[0] * ", varying IV strength " * pi),
                  vary_n_iv = expression("Rejection rates under " * H[0] * ", varying " * n[IV]))
  ggplot(data_long, aes(x = .data[[x_var]], y = rejection_rate, linetype = method, color = method)) +
    geom_hline(yintercept = 0.05, color = "grey") +
    geom_hline(yintercept = r025, linetype = "dashed", color = "grey") +
    geom_hline(yintercept = r975, linetype = "dashed", color = "grey") +
    geom_line(linewidth = 0.6) +
    scale_linetype_manual(values = line_types, name = "Method") +
    scale_color_manual(values = line_colors, name = "Method") +
    scale_x_continuous(breaks = x_breaks) +
    facet_wrap(
      ~ panel,
      scales = if(sim == "vary_n_iv") "fixed" else "free_y",
      ncol = 4,
      labeller = labeller(
        panel = function(x) {
          gsub(gsub1, gsub2, x)
        },
        .default = label_parsed
      )
    ) +
    labs(x = x_lab, y = "Rejection Rate", linetype = "Method", title = title) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 10),
      legend.position = "bottom",
      legend.key.size = unit(0.8, "lines"),       
      legend.key.width = unit(1.2, "lines"),        
      legend.spacing.x = unit(0.2, "lines"),      
      legend.spacing.y = unit(0.2, "lines"),      
      legend.text = element_text(size = 9),
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.spacing = unit(0.8, "lines")
    ) + 
    {if(sim == "vary_n_iv") coord_cartesian(ylim = c(0, 0.3)) else NULL}
  if(sim == "vary_n"){
    ggsave("Plots/H0_vary_n.pdf", width = 8.5, height = 5)
  }
  if(sim == "vary_iv_strength"){
    ggsave("Plots/H0_vary_iv_strength.pdf", width = 8.5, height = 5)
  }
  if(sim == "vary_n_iv"){
    ggsave("Plots/H0_vary_n_iv.pdf", width = 8.5, height = 3.4)
  }
}

plot_results("vary_n")
plot_results("vary_iv_strength")
plot_results("vary_n_iv")

