#
# figure S10
#

library(dplyr)
library(mgcv)

# ------------------------------------ LOAD DATA ---------------------------------------%

load('data/data.RData')
source('rCodes/statFuncs.R')
source('rCodes/figureFuncs.R')

# ---------------------------------- PREPARE DATA ---------------------------------------%

out <- getSeroconversionData(conversion, hi_paired, st)

# ------------------------------------ SAVE PLOT ---------------------------------------%

jpeg('figures/figure_s.jpeg',
     width = 12*300, height = 12*300, res = 300)

layout(matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE),
       widths = rep(1, 2), heights = rep(1, 2))
sapply(1:4, getSplinePlotForMetric, out = out)

dev.off()
