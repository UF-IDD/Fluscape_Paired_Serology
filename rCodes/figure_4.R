#
# figure 4
#

library(dplyr)
library(mgcv)

# ------------------------------------ LOAD DATA ---------------------------------------%

load('data/data.RData')
source('rCodes/figureFuncs.R')
source('rCodes/statFuncs.R')

# ------------------------------------ PREPARE DATA ---------------------------------------%

stat <- birthAdjustedMetrics(hi_paired, conversion)

models_auc <- aucGam(stat)
models_aty <- atyGam(stat)

tem = stat[ ,c('participant_id', 'age_at_sampling', 
               'width_V1_40', 'width_V4_40', 'width_delta_40')]
colnames(tem) = c('participant_id', 'age_at_sampling', 
                  'width_V1', 'width_V4', 'width_delta')
models_w40 <- widthGam(tem)

# ------------------------------------ SAVE PLOT ---------------------------------------%

jpeg('figures/figure_4.jpeg',
     width = 14*300, height = 9*300, res = 300)

layout(matrix(1:12, 
              nrow = 3, ncol = 4, 
              byrow = TRUE),
       widths = c(2, rep(1.5, 3)), heights = rep(1,3))

# auc 
conceptualPlotStatistics(hi_paired = hi_paired, 
                         conversion = conversion,
                         panel = 1)
sapply(1:3, function(panel) aucGamDataBirthAdjustedPlot(stat = stat, 
                                                        m = models_auc, 
                                                        panel = panel,
                                                        legend = LETTERS[2:4]))

# width 1:40
conceptualPlotStatistics(hi_paired = hi_paired, 
                         conversion = conversion,
                         panel = 2)
sapply(1:3, function(panel) widthGamDataBirthAdjustedPlot(stat = tem, 
                                                          m = models_w40, 
                                                          panel = panel, 
                                                          cutoff = 3, 
                                                          spp = FALSE,
                                                          legend = LETTERS[6:8]))

# aty
conceptualPlotStatistics(hi_paired = hi_paired, 
                         conversion = conversion,
                         panel = 3)

sapply(1:3, function(panel) atyGamDataBirthAdjustedPlot(stat = stat, 
                                                        m = models_aty, 
                                                        panel = panel,
                                                        legend = LETTERS[10:12]))

dev.off()