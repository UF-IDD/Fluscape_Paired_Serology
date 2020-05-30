#
# figure S2
#

# ------------------------------------ LOAD DATA ---------------------------------------%

library(mgcv)
load('data/data.RData')
source('rCodes/statFuncs.R')
source('rCodes/figureFuncs.R')

# ------------------------------------ PREPARE DATA ---------------------------------------%

# auc
stat_auc <- aucData(data = hi_paired)
models_auc <- aucGam(stat_auc)

# w40
stat_w40 <- widthData(data = hi_paired, z = 3, normalize = TRUE)
models_w40 <- widthGam(stat_w40)

# w10
stat_w10 <- widthData(data = hi_paired, z = 1, normalize = TRUE)
models_w10 <- widthGam(stat_w10)

# aty
stat_aty <- atyData(hi_paired = hi_paired ,data = conversion)
models_aty <- atyGam(stat_aty)


# ------------------------------------ SAVE PLOT ---------------------------------------%


jpeg('figures/figure_s2.jpeg',
     width = 9*300, height = 14*300, res = 300)

layout(matrix(1:12, 
              nrow = 4, ncol = 3, 
              byrow = TRUE),
       widths = rep(1.5, 3), heights = rep(1,4))

# auc
sapply(1:3, function(panel) aucGamDataPlot(stat = stat_auc, 
                                           m = models_auc, 
                                           panel = panel, 
                                           all = TRUE, 
                                           combinedPlot = TRUE))

# width
sapply(1:3, function(panel) widthGamDataPlot(stat = stat_w40, 
                                             m = models_w40, 
                                             panel = panel, 
                                             cutoff = 3))

sapply(1:3, function(panel) widthGamDataPlot(stat = stat_w10, 
                                             m = models_w10, 
                                             panel = panel, 
                                             cutoff = 1))

# aty
sapply(1:3, function(panel) atyGamDataPlot(stat = stat_aty, 
                                           m = models_aty, 
                                           panel = panel))


dev.off()
