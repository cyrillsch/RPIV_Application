# R script to reproduce the figures in Section 5.2 and 5.3 in the paper
################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)

method_levels <- c(
  "RP_hom", "RP_het", "weakRP_hom", "weakRP_het",
  "J", "ICM", "smooth_asymp", "smooth_boot"
)

method_labels <- c(
  "RP hom.", "RP het.", "weak RP hom.", "weak RP het.",
  "overid. J", "ICM", "smooth asymp.", "smooth boot."
)

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

violation_facet_labels <- c(
  "Z_squared" = "z squared",
  "sign_Z" = "sign(z)",
  "misspec_squared" = "misspec. squared",
  "misspec_sign" = "misspec. sign"
)

nrep <- 1000

r025 <- qbinom(0.025, size = nrep, prob = 0.05)/nrep
r975 <- qbinom(0.975, size = nrep, prob = 0.05)/nrep


################################################################################
## H0, vary n
################################################################################

load("SimulationResults/H0_vary_n/row_65.RData")  

data_H0 <- data_H0 %>%
  mutate(
    het_type = ifelse(heteroskedastic, "heteroskedastic", "homoskedastic"),
    panel_raw = paste0("iv_", n_iv, "_ctrl_", n_control, "_", het_type)
  )


data_H0 <- data_H0 %>%
  filter(
    n_control == 2,
    n_iv %in% c(1, 2)
  )

panel_levels <- data_H0 %>%
  distinct(n_iv, n_control, heteroskedastic, panel_raw) %>%
  arrange(heteroskedastic, n_iv, n_control) %>%
  pull(panel_raw)

data_H0 <- data_H0 %>%
  mutate(panel = factor(panel_raw, levels = panel_levels))

data_long <- data_H0 %>%
  pivot_longer(
    cols = all_of(method_levels),
    names_to = "method",
    values_to = "rejection_rate"
  ) %>%
  mutate(
    method = factor(method, levels = method_levels, labels = method_labels)
  )


panel_labeller <- labeller(
  panel = function(x) {
    gsub(
      "iv_(\\d+)_ctrl_(\\d+)_(.*)",
      "atop(n[IV] == \\1*','~n[c] == \\2, \\3)",
      x
    )
  },
  .default = label_parsed
)


p <- ggplot(
  data_long,
  aes(
    x = n,
    y = rejection_rate,
    color = method,
    linetype = method
  )
) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ panel, ncol = 4, scales = "free_y", labeller = panel_labeller) +
  scale_color_manual(values = line_colors) +
  scale_linetype_manual(values = line_types) +
  scale_x_continuous(breaks = c(100, 200, 300, 400, 500)) +
  labs(
    x = expression(n),
    y = "Rejection rate",
    title = bquote("Rejection rates under " * H[0] * ", varying n")
  ) +
  geom_hline(yintercept = 0.05, color = "grey") +
  geom_hline(yintercept = r025, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = r975, linetype = "dashed", color = "grey") +
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


ggsave("Plots/Main_H0.pdf", p, width = 8, height = 3.4)


################################################################################
## H0: vary n_IV
################################################################################

load("SimulationResults/H0_vary_n_iv/row_21.RData")


data_H0_niv <- data_H0 %>%
  mutate(
    het_type = ifelse(heteroskedastic, "heteroskedastic", "homoskedastic"),
    panel_raw = paste0("ctrl_", n_control, "_", het_type)
  ) %>%
  filter(n_control == 2)


panel_levels <- data_H0_niv %>%
  distinct(heteroskedastic, panel_raw) %>%
  arrange(heteroskedastic) %>%  
  pull(panel_raw)

data_H0_niv <- data_H0_niv %>%
  mutate(panel = factor(panel_raw, levels = panel_levels))


data_long_niv <- data_H0_niv %>%
  pivot_longer(
    cols = all_of(method_levels),
    names_to = "method",
    values_to = "rejection_rate"
  ) %>%
  mutate(
    method = factor(method, levels = method_levels, labels = method_labels)
  )

panel_labeller_niv <- labeller(
  panel = function(x) {
    gsub(
      "ctrl_(\\d+)_(.*)",
      "atop(n[c] == \\1, \\2)",
      x
    )
  },
  .default = label_parsed
)

p_niv <- ggplot(
  data_long_niv,
  aes(
    x = n_iv,
    y = rejection_rate,
    color = method,
    linetype = method
  )
) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ panel, ncol = 2, scales = "fixed", labeller = panel_labeller_niv) +
  scale_color_manual(values = line_colors) +
  scale_linetype_manual(values = line_types) +
  scale_x_continuous(breaks = c(5, 10, 15, 20, 25)) +
  labs(
    x = expression(n[IV]),
    y = "Rejection rate",
    title = bquote("Rejection rates under " * H[0] * ", varying " * n[IV])
  ) +
  geom_hline(yintercept = 0.05, color = "grey") +
  geom_hline(yintercept = r025, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = r975, linetype = "dashed", color = "grey") +
  coord_cartesian(ylim = c(0, 0.3)) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "right",
    legend.key.size = unit(0.8, "lines"),
    legend.key.width = unit(1.2, "lines"),
    legend.spacing.x = unit(0.2, "lines"),
    legend.spacing.y = unit(0.2, "lines"),
    legend.text = element_text(size = 9),
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(0.8, "lines")
  )

ggsave("Plots/Main_H0_vary_n_iv.pdf", p_niv, width = 5, height = 2.6)




################################################################################
## HA, vary violation strength s
################################################################################


load("SimulationResults/HA_vary_gamma/row_801.RData")

data_HA <- data_HA %>%
  mutate(
    het_type = ifelse(heteroskedastic, "heteroskedastic", "homoskedastic")
  )


data_HA_sub <- data_HA %>%
  filter(
    heteroskedastic == FALSE,
    n_iv == 1,
    n_control == 2
  )


data_long <- data_HA_sub %>%
  pivot_longer(
    cols = all_of(method_levels),
    names_to = "method",
    values_to = "rejection_rate"
  ) %>%
  mutate(
    violation = as.character(violation),
    method = factor(method, levels = method_levels, labels = method_labels)
  ) %>%
  mutate(
    violation = factor(
      violation,
      levels = c("Z_squared", "sign_Z", "misspec_squared", "misspec_sign")
    )
  )



p <- ggplot(
  data_long,
  aes(
    x = gamma,
    y = rejection_rate,
    color = method,
    linetype = method
  )
) +
  geom_line(linewidth = 0.6) +
  facet_wrap(
    ~ violation,
    ncol = 4,
    scales = "free",
    labeller = labeller(violation = violation_facet_labels)
  ) +
  scale_color_manual(values = line_colors) +
  scale_linetype_manual(values = line_types) +
  labs(
    x = expression("Violation strength " * s[viol]),
    y = "Rejection rate",
    title = bquote("Rejection rates under " * H[A] * ", varying violation strength, homoskedastic, " * n[IV] == 1 * "," ~ n[c] == 2)
  ) +
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


ggsave("Plots/Main_HA.pdf", p, width = 8, height = 3.4)



################################################################################
## HA, vary IV strength pi
################################################################################

load("SimulationResults/HA_vary_iv_strength/row_225.RData")

data_HA_pi <- data_HA %>%
  mutate(
    het_type = ifelse(heteroskedastic, "heteroskedastic", "homoskedastic")
  ) %>%
  filter(
    heteroskedastic == FALSE,
    n_iv == 1,
    n_control == 2
  )

data_long_pi <- data_HA_pi %>%
  pivot_longer(
    cols = all_of(method_levels),
    names_to = "method",
    values_to = "rejection_rate"
  ) %>%
  mutate(
    violation = as.character(violation),
    method = factor(method, levels = method_levels, labels = method_labels),
    violation = factor(
      violation,
      levels = c("Z_squared", "sign_Z", "misspec_squared", "misspec_sign")
    )
  )

p_pi <- ggplot(
  data_long_pi,
  aes(
    x = iv_strength,      
    y = rejection_rate,
    color = method,
    linetype = method
  )
) +
  geom_line(linewidth = 0.6) +
  facet_wrap(
    ~ violation,
    ncol = 4,
    scales = "free",
    labeller = labeller(violation = violation_facet_labels)
  ) +
  scale_color_manual(values = line_colors) +
  scale_linetype_manual(values = line_types) +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  labs(
    x = expression("IV strength " * pi),
    y = "Rejection rate",
    title = bquote("Rejection rates under " * H[A] * ", homoskedastic, " *
                     n[IV] == 1 * "," ~ n[c] == 2 * ", varying IV strength " * pi)
  ) +
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

ggsave("Plots/Main_HA_vary_iv_strength.pdf", p_pi, width = 8, height = 3.4)





