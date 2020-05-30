#
# figure 2
#

library(reshape2)
library(mgcv)

# ------------------------------------ LOAD DATA ---------------------------------------%

load('data/data.RData')
source('rCodes/figureFuncs.R')
source('rCodes/statFuncs.R')

# ------------------------------------ SAVE PLOT ---------------------------------------%

# jpeg('figures/figure_2.jpeg',
#      width = 15*300, height = 11*300, res = 300)
# 
# 
# layout(matrix(c(1:9,
#                 10,10,11), 
#               nrow = 3, ncol = 4, 
#               byrow = FALSE),
#        widths = c(3.3, 2.2, 2, 0.7), heights = rep(2,3))

jpeg('figures/figure_2.jpeg',
     width = 18*300, height = 13*300, res = 300)


layout(matrix(c(1:9,
                10,10,11), 
              nrow = 3, ncol = 4, 
              byrow = FALSE),
       widths = c(4.2, 2.2, 2, 0.7), 
       heights = rep(2, 3))

# heatmap
sapply(1:2, function(x) titerAgeHeatmap(hi_paired, 
                                        st = st, 
                                        panel = x))
converAgeHeatmap(conversion, st)

# titer distribution
sapply(1:3, titerDistribution, conversion = conversion)

# contour plot
sapply(1:3, contourPlot, hi_paired = hi_paired, conversion = conversion)

# scale plot
scalePlot2()
converScalePlot()

dev.off()