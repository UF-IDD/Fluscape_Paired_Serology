#
# figure S15
#

library(dplyr)

# ------------------------------------ LOAD DATA ---------------------------------------%

load('data/data.RData')
source('rCodes/statFuncs.R')
source('rCodes/figureFuncs.R')

# ---------------------------------- PREPARE DATA ---------------------------------------%

coeff.uni <- getUniCoeff(conversion, st)

coeff.mul <- getMulCoeff(hi_paired, conversion, st)

# ------------------------------------ SAVE PLOT ---------------------------------------%

jpeg('figures/figure_s15.jpeg',
     width = 15*300, height = 7*300, res = 300)

layout(matrix(1:2, nrow = 1, ncol = 2),
       widths = rep(1.5, 2), heights = 1)
coefHeatmap(coeff.uni, uni = TRUE)
coefHeatmap(coeff.mul, uni = FALSE)

dev.off()
