#
# figure 3
#

# ------------------------------------ LOAD DATA ---------------------------------------%

library(mgcv)
load('data/data.RData')
source('rCodes/figureFuncs.R')


# ------------------------------------ SAVE PLOT ---------------------------------------%


jpeg('figures/figure_3.jpeg',
     width = 15*300, height = 6*300, res = 300)

layout(matrix(c(1:3, 7, 4:7), nrow = 2, ncol = 4, byrow = TRUE),
       w = c(1,1,1,0.5), h = rep(1,2))

# age at sampling
sapply(1:3, splineTiterOnAgePlotCombined, 
       st = st,
       data = conversion, 
       age.type = "sampling")

# age at circulation
sapply(1:3, splineTiterOnAgePlotCombined, 
       st = st,
       data = conversion, 
       age.type = "isolation")

# legend
splineTiterOnAgePlotLegend(st)

dev.off()
