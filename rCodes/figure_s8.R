#
# figure S8
#

library(mediation)
library(dplyr)

# ------------------------------------ LOAD DATA ---------------------------------------%

load('data/data.RData')
source('rCodes/statFuncs.R')
source('rCodes/figureFuncs.R')

# ---------------------------------- PREPARE DATA ---------------------------------------%

# generate data used for logistic regressions
out <- getSeroconversionData(conversion, hi_paired, st)

medEff <- getMedEff(out, interact = FALSE)
medEffInteract <- getMedEff(out, interact = TRUE)


# ------------------------------------ SAVE PLOT ---------------------------------------%

jpeg('figures/figure_s8.jpeg',
     width = 15*300, height = 4*300, res = 300)

layout(matrix(1:5, nrow = 1, ncol = 5, byrow = TRUE),
       widths = c(1.5, rep(1.5, 4)), heights = 2)
sapply(0:4, function(p) medEffPlot(medEff, medEffInteract, panel = p, st))

dev.off()
