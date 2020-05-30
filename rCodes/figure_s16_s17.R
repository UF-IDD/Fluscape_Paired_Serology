#
# figure S16 and S17
#

library(dplyr)
library(mgcv)

# ------------------------------------ LOAD DATA ---------------------------------------%

load('data/data.RData')
source('rCodes/statFuncs.R')
source('rCodes/figureFuncs.R')

# ---------------------------------- PREPARE DATA ---------------------------------------%

out <- getSeroconversionData(conversion, hi_paired, st)

pred <- getPred(out, spline = FALSE)
pred.spline <- getPred(out, spline = TRUE)

# ------------------------------------ SAVE PLOT ---------------------------------------%

jpeg('figures/figure_s16.jpeg',
     width = 10*300, height = 10*300, res = 300)

layout(matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE),
       widths = rep(1, 2), heights = rep(1, 2))
sapply(1:4, predPlot, pred = pred)

dev.off()


jpeg('figures/figure_s17.jpeg',
     width = 10*300, height = 10*300, res = 300)

layout(matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE),
       widths = rep(1, 2), heights = rep(1, 2))
sapply(1:4, predPlot, pred = pred.spline)

dev.off()
