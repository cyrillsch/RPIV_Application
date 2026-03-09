# R script to reproduce the real data analyses in Section 5.4 of the paper
# and Section E of the supplement
#####################################################

library(ivdoctr)
library(wooldridge)
library(ivreg)
library(parallel)

library(dplyr)
library(tidyr)
library(ggplot2)

source("AdditionalImplementations.R")

# to install the RPIV package:
# devtools::install_github("cyrillsch/RPIV")
library(RPIV)



## Card dataset
# Analysis according to Wooldridge, J. M. (2013),
# Introductory Econometrics: A Modern Approach, 5 edn, South-Western Cengage Learning, Mason

data(card)

# with expersq
fit1 <- ivreg(lwage ~ educ + exper + expersq  + black + smsa + south +  smsa66 +  reg662 + reg663 + reg664 + reg665 + reg666 + reg667 + reg668 + reg669|nearc4 + exper + expersq + black + smsa + south +  smsa66 + reg662 + reg663 + reg664 + reg665 + reg666 + reg667 + reg668 + reg669, data = card)
summary(fit1)
# this replicates table 15.1 in Wooldridge (2013)

# without expersq
fit0 <- ivreg(lwage ~ educ + exper  + black + smsa + south +  smsa66 +  reg662 + reg663 + reg664 + reg665 + reg666 + reg667 + reg668 + reg669|nearc4 + exper + black + smsa + south +  smsa66 + reg662 + reg663 + reg664 + reg665 + reg666 + reg667 + reg668 + reg669, data = card)
summary(fit0)



Y <-  card$lwage
X <- card$educ
Z <- card$nearc4
# without expersq
C0 <- with(card, cbind(exper, black, smsa, south, smsa66, reg662, reg663, reg664, reg665, reg666, reg667, reg668, reg669))
# with expersq
C1 <- with(card, cbind(exper, expersq, black, smsa, south, smsa66, reg662, reg663, reg664, reg665, reg666, reg667, reg668, reg669))

# for replicability of parallelization
RNGkind("L'Ecuyer-CMRG")
set.seed(915)

ncores <- 50
B <- 100

beta_vec <- seq(-0.5, 1, length.out = 100)

smooth0 <- smooth_Test(Y, X, C0, Z, B = 1999)
smooth1 <- smooth_Test(Y, X, C1, Z, B = 1999)

RPIV0_hom <- mclapply(1:B, function(i){RPIV_test(Y = Y, X = X, C = C0, Z = Z, variance_estimator = "homoskedastic")$p_value}, mc.cores = ncores)
RPIV0_het <- mclapply(1:B, function(i){RPIV_test(Y = Y, X = X, C = C0, Z = Z, variance_estimator = "heteroskedastic")$p_value}, mc.cores = ncores)

RPIV1_hom <- mclapply(1:B, function(i){RPIV_test(Y = Y, X = X, C = C1, Z = Z, variance_estimator = "homoskedastic")$p_value}, mc.cores = ncores)
RPIV1_het <- mclapply(1:B, function(i){RPIV_test(Y = Y, X = X, C = C1, Z = Z, variance_estimator = "heteroskedastic")$p_value}, mc.cores = ncores)


ICM0 <- ICM_test(Y, X, C0, Z, B = 1999)
ICM1 <- ICM_test(Y, X, C1, Z, B = 1999)

ICM0 <- sapply(beta_vec, ICM0)
ICM1 <- sapply(beta_vec, ICM1)

get_weakRPIV <- function(C){
  temp_weakRP <- weak_RPIV_test(Y, X, C, Z, variance_estimator = c("homoskedastic", "heteroskedastic"))
  return(sapply(beta_vec, temp_weakRP, type = "fit"))
}

weakRPIV0 <- mclapply(1:B, function(i){get_weakRPIV(C0)}, mc.cores = ncores)
weakRPIV1 <- mclapply(1:B, function(i){get_weakRPIV(C1)}, mc.cores = ncores)




save(smooth0, smooth1, RPIV0_hom, RPIV0_het, RPIV1_hom, RPIV1_het, ICM0, ICM1, weakRPIV0, weakRPIV1, file = "RealDataResults/Card/p_values.RData")

load("RealDataResults/Card/p_values.RData")


RPIV0_het <- do.call(c, RPIV0_het)
RPIV0_hom <- do.call(c, RPIV0_hom)
RPIV1_het <- do.call(c, RPIV1_het)
RPIV1_hom <- do.call(c, RPIV1_hom)

weakRPIV0_hom <- 1-pnorm(do.call(rbind, lapply(weakRPIV0, function(m){m[1,]})))
weakRPIV0_het <- 1-pnorm(do.call(rbind, lapply(weakRPIV0, function(m){m[2,]})))

weakRPIV1_hom <- 1-pnorm(do.call(rbind, lapply(weakRPIV1, function(m){m[1,]})))
weakRPIV1_het <- 1-pnorm(do.call(rbind, lapply(weakRPIV1, function(m){m[2,]})))

weakRPIV0_hom <- 2 * apply(weakRPIV0_hom, 2, median)
weakRPIV0_het <- 2 * apply(weakRPIV0_het, 2, median)
weakRPIV1_hom <- 2 * apply(weakRPIV1_hom, 2, median)
weakRPIV1_het <- 2 * apply(weakRPIV1_het, 2, median)

RPIV0_het <- 2 * median(RPIV0_het)
RPIV0_hom <- 2 * median(RPIV0_hom)
RPIV1_het <- 2 * median(RPIV1_het)
RPIV1_hom <- 2 * median(RPIV1_hom)




n_beta <- length(beta_vec)

df_plot <- data.frame(beta = beta_vec, weakRP_hom = c(weakRPIV0_hom, weakRPIV1_hom), weakRP_het = c(weakRPIV0_het, weakRPIV1_het),
                      ICM = c(ICM0, ICM1), smooth_asymp = rep(c(smooth0[1], smooth1[1]), each = n_beta),
                      smooth_boot = rep(c(smooth0[2], smooth1[2]), each = n_beta),
                      RP_hom = rep(c(RPIV0_hom, RPIV1_hom), each = n_beta), RP_het = rep(c(RPIV0_het, RPIV1_het), each = n_beta),
                      setting = rep(c("without expersq", "with expersq"), each = length(beta_vec)))





# Methods settings
method_levels <- c(
  "RP_hom", "RP_het", "weakRP_hom", "weakRP_het",
   "ICM", "smooth_asymp", "smooth_boot"
)

method_labels <- c(
  "RP hom.", "RP het.", "weak RP hom.", "weak RP het.",
   "ICM", "smooth asymp.", "smooth boot."
)


line_types <- c(
  "RP hom." = "dashed",
  "RP het." = "solid",
  "weak RP hom." = "dashed",
  "weak RP het." = "solid",
  "ICM" = "dotted",
  "smooth asymp." = "dotted",
  "smooth boot." = "dotted"
)

line_colors <- c(
  "RP hom." = "blue",
  "RP het." = "blue",
  "weak RP hom." = "darkgreen",
  "weak RP het." = "darkgreen",
  "ICM" = "red",
  "smooth asymp." = "purple",
  "smooth boot." = "magenta"
)


df_long <- df_plot %>%
  pivot_longer(cols = all_of(method_levels),
               names_to = "method", values_to = "p_value") %>%
  mutate(
    method = factor(method, levels = method_levels, labels = method_labels)
  )

# replace 0 by very small p-value
df_long$p_value <- ifelse(df_long$p_value == 0, 1e-12, df_long$p_value)


p <- ggplot(
  df_long,
  aes(x = beta, y = p_value, color = method, linetype = method)
) +
  geom_hline(yintercept = 0.05, color = "grey") +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ setting, ncol = 2) +
  scale_linetype_manual(values = line_types) +
  scale_color_manual(values = line_colors) +
  labs(
    x = expression(beta[0]),
    y = "p-value"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.key.size = unit(0.8, "lines"),
    legend.key.width = unit(1.2, "lines"),
    legend.spacing.x = unit(0.2, "lines"),
    legend.spacing.y = unit(0.2, "lines"),
    legend.text = element_text(size = 9),
    panel.spacing = unit(0.8, "lines")
  ) +
  scale_y_log10(
    limits = c(5e-9, NA),
    breaks = c(5e-2, 5e-4, 5e-6, 5e-8),
    labels = scales::label_scientific(),
    oob = scales::squish
  )

ggsave(
  filename = paste0("Plots/Card.pdf"),
  plot = p,
  width = 5,
  height = 3
)






################################################################################

## Protestantism and Literacy Rate from 
# Becker, S. O. & Woessmann, L. (2009), ‘Was Weber wrong? A human capital theory of
# protestant economic history’, The Quarterly Journal of Economics 124(2), 531–596.

# Data is avaliable from library(ivdoctr)
data(weber)

Y <- weber$f_rw
X <- weber$f_prot
Z <- weber$kmwittenberg
C <- with(weber, cbind(f_young, f_jew, f_fem, f_ortsgeb, f_pruss, hhsize, lnpop, gpop, f_miss, f_blind, f_deaf, f_dumb))

fit0 <- ivreg(Y ~ X + C| Z + C)
summary(fit0)
# this replicates the middle column Table III in Becker and Woessmann (2009)
# (with very small differences in the standard error estimates)

# include Z^2 as instrument and look at the Sargan test
fit1 <- ivreg(Y ~ X + C | Z + I(Z^2) + C)
summary(fit1)

# p-value from J-test with Z^2 is
J <- summary(fit1)$diagnostics[3,4]


# for replicability of parallelization
RNGkind("L'Ecuyer-CMRG")
set.seed(915)

ncores <- 50
B <- 100

beta_vec <- seq(-0.5, 0.5, length.out = 100)

smooth <- smooth_Test(Y, X, C, Z)


RPIV_hom <- mclapply(1:B, function(i){RPIV_test(Y = Y, X = X, C = C, Z = Z, variance_estimator = "homoskedastic")$p_value}, mc.cores = ncores)
RPIV_het <- mclapply(1:B, function(i){RPIV_test(Y = Y, X = X, C = C, Z = Z, variance_estimator = "heteroskedastic")$p_value}, mc.cores = ncores)



ICM <- ICM_test(Y, X, C, Z)


ICM <- sapply(beta_vec, ICM)


get_weakRPIV <- function(C){
  temp_weakRP <- weak_RPIV_test(Y, X, C, Z, variance_estimator = c("homoskedastic", "heteroskedastic"))
  return(sapply(beta_vec, temp_weakRP, type = "fit"))
}

weakRPIV <- mclapply(1:B, function(i){get_weakRPIV(C)}, mc.cores = ncores)

save(J, smooth, RPIV_hom, RPIV_het, ICM, weakRPIV, file = "RealDataResults/Weber/p_values.RData")


load("RealDataResults/Weber/p_values.RData")

RPIV_het <- do.call(c, RPIV_het)
RPIV_hom <- do.call(c, RPIV_hom)

RPIV_het <- 2 * median(RPIV_het)
RPIV_hom <- 2 * median(RPIV_hom)


weakRPIV_hom <- 1 - pnorm(do.call(rbind, lapply(weakRPIV, function(m) m[1, ])))
weakRPIV_het <- 1 - pnorm(do.call(rbind, lapply(weakRPIV, function(m) m[2, ])))

weakRPIV_hom <- 2 * apply(weakRPIV_hom, 2, median)
weakRPIV_het <- 2 * apply(weakRPIV_het, 2, median)

n_beta <- length(beta_vec)

df_plot <- data.frame(
  beta = beta_vec,
  weakRP_hom   = weakRPIV_hom,
  weakRP_het   = weakRPIV_het,
  ICM          = ICM,
  smooth_asymp = rep(smooth[1], n_beta),
  smooth_boot  = rep(smooth[2], n_beta),
  RP_hom       = rep(RPIV_hom, n_beta),
  RP_het       = rep(RPIV_het, n_beta),
  J            = J
)

method_levels <- c(
  "RP_hom", "RP_het", "weakRP_hom", "weakRP_het",
  "J", "ICM", "smooth_asymp", "smooth_boot"
)

method_labels <- c(
  "RP hom.", "RP het.", "weak RP hom.", "weak RP het.",
  "overid. J", "ICM", "smooth asymp.", "smooth boot."
)

line_types <- c(
  "RP hom."        = "dashed",
  "RP het."        = "solid",
  "weak RP hom."   = "dashed",
  "weak RP het."   = "solid",
  "overid. J"      = "dotted",
  "ICM"            = "dotted",
  "smooth asymp."  = "dotted",
  "smooth boot."   = "dotted"
)

line_colors <- c(
  "RP hom."        = "blue",
  "RP het."        = "blue",
  "weak RP hom."   = "darkgreen",
  "weak RP het."   = "darkgreen",
  "overid. J"      = "orange",
  "ICM"            = "red",
  "smooth asymp."  = "purple",
  "smooth boot."   = "magenta"
)


df_long <- df_plot %>%
  pivot_longer(
    cols = all_of(method_levels),
    names_to = "method",
    values_to = "p_value"
  ) %>%
  mutate(
    method = factor(method, levels = method_levels, labels = method_labels)
  )

# replace 0 p-values by minimum nonzero p_value
min_p_value <- min(df_long$p_value[df_long$p_value > 0])

# replace 0 by very small p-value
df_long$p_value <- ifelse(df_long$p_value == 0, min_p_value, df_long$p_value)


p <- ggplot(
  df_long,
  aes(x = beta, y = p_value, color = method, linetype = method)
) +
  geom_hline(yintercept = 0.05, color = "grey") +
  geom_line(linewidth = 0.6) +
  scale_linetype_manual(values = line_types) +
  scale_color_manual(values = line_colors) +
  labs(
    x = expression(beta[0]),
    y = "p-value"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.size = unit(0.8, "lines"),
    legend.key.width = unit(1.2, "lines"),
    legend.spacing.x = unit(0.2, "lines"),
    legend.spacing.y = unit(0.2, "lines"),
    legend.text = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)
  ) +
  scale_y_log10(
    breaks = c(5e-2, 5e-6, 5e-10, 5e-14),
    labels = scales::label_scientific()
  )


ggsave(
  filename = "Plots/Weber.pdf",
  plot = p,
  width = 5,
  height = 3.5
)



