# R script to reproduce the figures in Section D.4 of the supplement
################################################################################


library(ggplot2)
library(dplyr)
library(tidyr)

# Methods settings
method_levels <- c(
  "RP_hom", "RP_het", "weakRP_hom", "weakRP_het",
  "J", "ICM", "smooth_asymp", "smooth_boot"
)

method_labels <- c(
  "RP hom.", "RP het.", "weak RP hom.", "weak RP het.",
  "overid. J", "ICM", "smooth asymp.", "smooth boot."
)

violation_labels <- c("Z_squared", "sign_Z", "misspec_squared", "misspec_sign")

line_types <- c(
  "RP hom." = "dashed",
  "RP het." = "solid",
  "weak RP hom." = "dashed",
  "weak RP het." = "solid",
  "overid. J" = "dotted",
  "ICM" = "dotted",
  "smooth asymp." = "dotted",
  "smooth boot." = "dotted"
)

line_colors <- c(
  "RP hom." = "blue",
  "RP het." = "blue",
  "weak RP hom." = "darkgreen",
  "weak RP het." = "darkgreen",
  "overid. J" = "orange",
  "ICM" = "red",
  "smooth asymp." = "purple",
  "smooth boot." = "magenta"
)



plot_HA_results <- function(sim){
  filename <- switch(
    sim,
    vary_gamma       = "SimulationResults/HA_vary_gamma/row_801.RData",
    vary_iv_strength = "SimulationResults/HA_vary_iv_strength/row_225.RData",
    vary_n_iv        = "SimulationResults/HA_vary_n_iv/row_81.RData"
  )
  load(filename)
  data_HA <- data_HA %>%
    mutate(
      het_type = ifelse(heteroskedastic, "heteroskedastic", "homoskedastic")
    )
  if(sim == "vary_n_iv"){
    data_HA <- data_HA %>%
      mutate(panel_raw = paste0(
        "ctrl_", n_control,
        "_", het_type
      ))
  } else {
    data_HA <- data_HA %>%
      mutate(panel_raw = paste0(
        "iv_", n_iv,
        "_ctrl_", n_control,
        "_", het_type
      ))
  }
  if(sim == "vary_n_iv"){
    panel_levels <- data_HA %>%
      distinct(n_control, heteroskedastic, panel_raw) %>%
      arrange(
        heteroskedastic,
        n_control
      ) %>%
      pull(panel_raw)
  } else {
    panel_levels <- data_HA %>%
      distinct(n_iv, n_control, heteroskedastic, panel_raw) %>%
      arrange(
        heteroskedastic,
        n_iv,
        n_control
      ) %>%
      pull(panel_raw)
  }
  
  
  # Convert to factor
  data_HA <- data_HA %>%
    mutate(panel = factor(panel_raw, levels = panel_levels))
  
  data_HA <- data_HA %>%
    mutate(panel = factor(panel_raw, levels = panel_levels))
  
  # Reshape to long format
  data_long <- data_HA %>%
    pivot_longer(
      cols = all_of(method_levels),
      names_to = "method",
      values_to = "rejection_rate"
    ) %>%
    mutate(
      violation = as.character(violation),
      method = factor(method, levels = method_levels, labels = method_labels)
    )
  
  # X-axis settings
  x_var <- switch(
    sim,
    vary_gamma       = "gamma",
    vary_iv_strength = "iv_strength",
    vary_n_iv        = "n_iv"
  )
  
  x_lab <- switch(sim,
                  vary_gamma = expression("Violation strength " * s[viol]),
                  vary_iv_strength = expression("IV strength " * pi),
                  vary_n_iv = expression(n[IV]))
  
  x_breaks <- switch(
    sim,
    vary_gamma       = NULL,
    vary_iv_strength = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    vary_n_iv        = c(5, 10, 15, 20, 25)
  )
  
  
  for (v in violation_labels){
    v_type = switch(v,
                    Z_squared = "z squared",
                    sign_Z = "sign(z)",
                    misspec_squared = "misspec. squared",
                    misspec_sign = "misspec. sign")
    title <- switch(sim,
                    vary_gamma = bquote("Rejection rates under " * H[A] * ", varying " * s[viol] * ", violation type" ~ .(v_type)),
                    vary_iv_strength = bquote("Rejection rates under " * H[A] * ", varying IV strength " * pi * ", violation type" ~ .(v_type)),
                    vary_n_iv = bquote("Rejection rates under " * H[A] * ", varying " * n[IV] * ", violation type" ~ .(v_type)))
    
    df_v <- data_long %>% filter(violation == v)
    if (sim == "vary_gamma"){
      x_breaks <- seq(range(df_v$gamma)[1], range(df_v$gamma)[2], length.out = 5)
    }
    panel_labeller <- if (sim == "vary_n_iv") {
      labeller(
        panel = function(x) {
          gsub(
            "ctrl_(\\d+)_(.*)",
            "atop(n[c] == \\1, \\2)",
            x
          )
        },
        .default = label_parsed
      )
    } else {
      labeller(
        panel = function(x) {
          gsub(
            "iv_(\\d+)_ctrl_(\\d+)_(.*)",
            "atop(n[IV] == \\1*','~n[c] == \\2, \\3)",
            x
          )
        },
        .default = label_parsed
      )
    }
    p <- ggplot(
      df_v,
      aes(
        x = .data[[x_var]],
        y = rejection_rate,
        color = method,
        linetype = method
      )
    ) +
      geom_line(linewidth = 0.6) +
      facet_wrap(~ panel, ncol = 4, scales = "free_y", labeller = panel_labeller) +
      scale_color_manual(values = line_colors) +
      scale_linetype_manual(values = line_types) +
      scale_x_continuous(breaks = x_breaks) +
      labs(x = x_lab, y = "Rejection rate", title = title) +
      geom_hline(yintercept = 0.05, color = "grey") +
      theme_bw() +
      theme(
        strip.text = element_text(face = "bold", size = 10),
        legend.position = "bottom",
        legend.key.size = unit(0.8, "lines"),       
        legend.key.width = unit(1.2, "lines"),        
        legend.spacing.x = unit(0.2, "lines"),      
        legend.spacing.y = unit(0.2, "lines"),      
        legend.text = element_text(size = 9),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing = unit(0.8, "lines")
      ) 
    outfile <- paste0("Plots/HA_", sim, "_", v, ".pdf")
    if(sim == "vary_n_iv"){
      ggsave(outfile, p, width = 8.5, height = 3.4)
    } else {
      ggsave(outfile, p, width = 8.5, height = 5)
    }
  }
}

# Run
plot_HA_results("vary_gamma")
plot_HA_results("vary_iv_strength")
plot_HA_results("vary_n_iv")
