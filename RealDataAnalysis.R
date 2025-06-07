## R script to reproduce the tables from Section 4.2.
#####################################################

library(ivdoctr)
library(wooldridge)
library(ivreg)
library(parallel)

# to install the RPIV package:
# devtools::install_github("cyrillsch/RPIV")
library(RPIV)

## Card dataset
# Analysis according to Wooldridge (2013), Introductory Econometrics, A Modern Approach

data(card)

# with expersq
fit1 <- ivreg(lwage ~ educ + exper + expersq  + black + smsa + south +  smsa66 +  reg662 + reg663 + reg664 + reg665 + reg666 + reg667 + reg668 + reg669|nearc4 + exper + expersq + black + smsa + south +  smsa66 + reg662 + reg663 + reg664 + reg665 + reg666 + reg667 + reg668 + reg669, data = card)
summary(fit1)
# this replicates table 15.1 in Wooldridge (2013)

# without expersq
fit0 <- ivreg(lwage ~ educ + exper  + black + smsa + south +  smsa66 +  reg662 + reg663 + reg664 + reg665 + reg666 + reg667 + reg668 + reg669|nearc4 + exper + black + smsa + south +  smsa66 + reg662 + reg663 + reg664 + reg665 + reg666 + reg667 + reg668 + reg669, data = card)
summary(fit0)

ncores <- 15
B <- 50

# for replicability of parallelization
RNGkind("L'Ecuyer-CMRG")
set.seed(915)

Y = card$lwage
X <- card$educ
Z <- card$nearc4

# without expersq
C0 <- with(card, cbind(exper, black, smsa, south, smsa66, reg662, reg663, reg664, reg665, reg666, reg667, reg668, reg669))
p0_Card_hom <- mclapply(1:B, function(i){RPIV_test(Y = Y, X = X, C = C0, Z = Z, variance_estimator = "homoskedastic")$p_value}, mc.cores = ncores)
p0_Card_het <- mclapply(1:B, function(i){RPIV_test(Y = Y, X = X, C = C0, Z = Z, variance_estimator = "heteroskedastic")$p_value}, mc.cores = ncores)

# with expersq
C1 <- with(card, cbind(exper, expersq, black, smsa, south, smsa66, reg662, reg663, reg664, reg665, reg666, reg667, reg668, reg669))
p1_Card_hom <- mclapply(1:B, function(i){RPIV_test(Y = Y, X = X, C = C1, Z = Z, variance_estimator = "homoskedastic")$p_value}, mc.cores = ncores)
p1_Card_het <- mclapply(1:B, function(i){RPIV_test(Y = Y, X = X, C = C1, Z = Z, variance_estimator = "heteroskedastic")$p_value}, mc.cores = ncores)

save(p0_Card_hom, p0_Card_het, p1_Card_hom, p1_Card_het, file = "RealDataResults/Card.RData")

load("RealDataResults/Card.RData")

# aggregate p-values by using twice the median

p0_Card_hom <- do.call(c, p0_Card_hom)
p0_Card_het <- do.call(c, p0_Card_het)

p1_Card_hom <- do.call(c, p1_Card_hom)
p1_Card_het <- do.call(c, p1_Card_het)

2*median(p0_Card_hom)
2*median(p0_Card_het)

2*median(p1_Card_hom)
2*median(p1_Card_het)


################################################################################

## Protestantism and Literacy Rate from Becker and Woessmann (2009), 
# Was weber wrong? a human capital theory of protestant economic history. 

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
summary(fit1)$diagnostics[3,4]

ncores <- 15
B <- 50

# for replicability of parallelization
RNGkind("L'Ecuyer-CMRG")


set.seed(915)
p_Becker_hom <- mclapply(1:B, function(i){RPIV_test(Y = Y, X = X, C = C, Z = Z, variance_estimator = "homoskedastic")$p_value}, mc.cores = ncores)
p_Becker_het <- mclapply(1:B, function(i){RPIV_test(Y = Y, X = X, C = C, Z = Z, variance_estimator = "heteroskedastic")$p_value}, mc.cores = ncores)

save(p_Becker_het, p_Becker_hom, file = "RealDataResults/Becker.RData")

load("RealDataResults/Becker.RData")

# aggregate p-values by using twice the median

p_Becker_hom <- do.call(c, p_Becker_hom)
p_Becker_het <- do.call(c, p_Becker_het)

2*median(p_Becker_hom)
2*median(p_Becker_het)

