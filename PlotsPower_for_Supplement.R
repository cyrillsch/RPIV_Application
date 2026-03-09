# R script to reproduce the figures in Section D.5 of the supplement
################################################################################

load("SimulationResults/WeakPower/row_73.RData")

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)


n_beta <- 50
beta_vec <- seq(-4, 2, length.out = n_beta)


# helper to reshape one row into beta-long format
reshape_beta_data <- function(data_H0, beta_vec){
  beta_df <- data_H0 %>%
    mutate(row_id = row_number()) %>%
    pivot_longer(
      cols = starts_with("weakRP_hom_beta"),
      names_to = "beta_id",
      values_to = "weakRP_hom"
    ) %>%
    mutate(beta_index = as.integer(str_extract(beta_id, "\\d+"))) %>%
    select(n, n_control, n_iv, iv_strength, heteroskedastic, row_id, weakRP_hom, beta_index) %>%
    left_join(
      data_H0 %>%
        mutate(row_id = row_number()) %>%
        pivot_longer(
          cols = starts_with("weakRP_het_beta"),
          names_to = "beta_id",
          values_to = "weakRP_het"
        ) %>%
        mutate(beta_index = as.integer(str_extract(beta_id, "\\d+"))) %>%
        select(row_id, beta_index, weakRP_het),
      by = c("row_id", "beta_index")
    ) %>%
    left_join(
      data_H0 %>%
        mutate(row_id = row_number()) %>%
        pivot_longer(
          cols = starts_with("ICM_beta"),
          names_to = "beta_id",
          values_to = "ICM"
        ) %>%
        mutate(beta_index = as.integer(str_extract(beta_id, "\\d+"))) %>%
        select(row_id, beta_index, ICM),
      by = c("row_id", "beta_index")
    ) %>%
    mutate(beta = beta_vec[beta_index])
  return(beta_df)
}

plot_beta_power <- function(iv_strength_value,
                            data_H0,
                            beta_vec){
  
  data_long <- reshape_beta_data(data_H0, beta_vec) %>%
    filter(iv_strength == iv_strength_value) %>%
    mutate(
      het_type = ifelse(heteroskedastic, "heteroskedastic", "homoskedastic"),
      panel_raw = paste0(
        "n_", n,
        "_iv_", n_iv,
        "_ctrl_", n_control,
        "_", het_type
      )
    )
  panel_levels <- data_long %>%
    distinct(n, n_iv, n_control, heteroskedastic, panel_raw) %>%
    arrange(
      n,
      heteroskedastic,
      n_iv,
      n_control
    ) %>%
    pull(panel_raw)
  
  data_long <- data_long %>%
    mutate(
      panel = factor(panel_raw, levels = panel_levels)
    ) %>%
    pivot_longer(
      cols = c(weakRP_hom, weakRP_het, ICM),
      names_to = "method",
      values_to = "rejection_rate"
    ) %>%
    mutate(
      method = factor(
        method,
        levels = c("weakRP_hom", "weakRP_het", "ICM"),
        labels = c("weak RP hom.", "weak RP het.", "ICM")
      )
    )
  
  line_types <- c(
    "weak RP hom." = "dashed",
    "weak RP het." = "solid",
    "ICM" = "dotted"
  )
  
  line_colors <- c(
    "weak RP hom." = "darkgreen",
    "weak RP het." = "darkgreen",
    "ICM" = "red"
  )
  
  title <- switch(as.character(iv_strength_value),
                  "0" = expression("Rejection rate of " * H[0](beta) * ", IV strength " * pi == 0),
                  "0.5" = expression("Rejection rate of " * H[0](beta) * ", IV strength " * pi == 0.5),
                  "1" = expression("Rejection rate of " * H[0](beta) * ", IV strength " * pi == 1))
  p <- ggplot(
    data_long,
    aes(x = beta, y = rejection_rate,
        linetype = method, color = method)
  ) +
    geom_hline(yintercept = 0.05, color = "grey") +
    geom_vline(xintercept = -1, color = "black", linetype = "dashed", linewidth = 0.1) +
    geom_line(linewidth = 0.6) +
    facet_wrap(
      ~ panel,
      ncol = 4,
      labeller = labeller(
        panel = function(x) {
          gsub(
            "n_(\\d+)_iv_(\\d+)_ctrl_(\\d+)_(.*)",
            "atop(n == \\1*','~n[IV] == \\2*','~n[c] == \\3, \\4)",
            x
          )
        },
        .default = label_parsed
      )
    ) +
    scale_linetype_manual(values = line_types) +
    scale_color_manual(values = line_colors) +
    labs(
      x = expression(beta[0]),
      y = "Rejection Rate",
      title = title
    ) +
    scale_x_continuous(breaks = seq(floor(min(beta_vec)),
                                    ceiling(max(beta_vec)),
                                    by = 2)) +
    theme_bw() +
    theme(
      strip.text = element_text(face = "bold", size = 9),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.spacing = unit(1, "lines"),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 10)
    )
  
  ggsave(
    filename = paste0("Plots/WeakPower_beta_ivstrength_", iv_strength_value, ".pdf"),
    plot = p,
    width = 8,
    height = 12
  )
}



plot_beta_power(0,   data_H0, beta_vec)
plot_beta_power(0.5, data_H0, beta_vec)
plot_beta_power(1,   data_H0, beta_vec)



