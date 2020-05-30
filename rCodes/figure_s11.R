#
# figure S11
#

library(dplyr)
library(mgcv)

# ------------------------------------ LOAD DATA ---------------------------------------%

load('data/data_never_vaccinated.RData')
source('rCodes/statFuncs.R')
source('rCodes/figureFuncs.R')

# ------------------------------------ PREPARE DATA ---------------------------------------%

stat <- birthAdjustedMetrics(hi_paired, conversion)

models_auc <- aucGam(stat)
models_aty <- atyGam(stat)

w40 = stat[ ,c('participant_id', 'age_at_sampling', 
               'width_V1_40', 'width_V4_40', 'width_delta_40')]
colnames(w40) = c('participant_id', 'age_at_sampling', 
                  'width_V1', 'width_V4', 'width_delta')
models_w40 <- widthGam(w40)


w10 = stat[ ,c('participant_id', 'age_at_sampling', 
               'width_V1_10', 'width_V4_10', 'width_delta_10')]
colnames(w10) = c('participant_id', 'age_at_sampling', 
                  'width_V1', 'width_V4', 'width_delta')
models_w10 <- widthGam(w10)


# ------------------------------------ SAVE PLOT ---------------------------------------%


jpeg('figures/figure_s11.jpeg',
     width = 9*300, height = 14*300, res = 300)

layout(matrix(1:12, 
              nrow = 4, ncol = 3, 
              byrow = TRUE),
       widths = rep(1.5, 3), heights = rep(1,4))

# auc
sapply(1:3, function(panel) aucGamDataBirthAdjustedPlot(stat = stat, 
                                                        m = models_auc, 
                                                        panel = panel,
                                                        legend = LETTERS[1:3]))

# width
sapply(1:3, function(panel) widthGamDataBirthAdjustedPlot(stat = w40, 
                                                          m = models_w40, 
                                                          panel = panel, 
                                                          cutoff = 3, 
                                                          spp = FALSE,
                                                          legend = LETTERS[4:6]))

sapply(1:3, function(panel) widthGamDataBirthAdjustedPlot(stat = w10, 
                                                          m = models_w10, 
                                                          panel = panel, 
                                                          cutoff = 1, 
                                                          spp = FALSE,
                                                          legend = LETTERS[7:9]))

# aty
sapply(1:3, function(panel) atyGamDataBirthAdjustedPlot(stat = stat, 
                                                        m = models_aty, 
                                                        panel = panel,
                                                        legend = LETTERS[10:12]))


dev.off()
