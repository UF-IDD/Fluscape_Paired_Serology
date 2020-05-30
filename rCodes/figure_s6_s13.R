#
# figure S6 and S13
#

library(mgcv)
library(dplyr)

# ------------------------------------ LOAD DATA ---------------------------------------%

load('data/data.RData')
source('rCodes/statFuncs.R')
source('rCodes/figureFuncs.R')

# ---------------------------------- PREPARE DATA ---------------------------------------%

# data used for regression on seroconversion 
out <- getSeroconversionData(conversion, hi_paired, st)


# ------------------------------------ SAVE PLOT ---------------------------------------%

jpeg('figures/figure_s6.jpeg',
     width = 15*300, height = 7*300, res = 300)

AICComparePlot(dat = out, spline = FALSE)

dev.off()


jpeg('figures/figure_s13.jpeg',
     width = 15*300, height = 7*300, res = 300)

AICComparePlot(dat = out, spline = TRUE)

dev.off()
