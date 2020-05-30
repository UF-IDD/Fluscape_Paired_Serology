#
# figure s3
#

library(dplyr)
library(mgcv)

# ------------------------------------ LOAD DATA ---------------------------------------%

load('data/data.RData')
source('rCodes/figureFuncs.R')
source('rCodes/statFuncs.R')

# ------------------------------------ PREPARE DATA ---------------------------------------%

stat <- birthAdjustedMetrics(hi_paired, conversion)

# w10
w10 = stat[ ,c('participant_id', 'age_at_sampling', 
               'width_V1_10', 'width_V4_10', 'width_delta_10')]
colnames(w10) = c('participant_id', 'age_at_sampling', 
                  'width_V1', 'width_V4', 'width_delta')
models_w10 <- widthGam(w10)

# w40
w40 = stat[ ,c('participant_id', 'age_at_sampling', 
               'width_V1_40', 'width_V4_40', 'width_delta_40')]
colnames(w40) = c('participant_id', 'age_at_sampling', 
                  'width_V1', 'width_V4', 'width_delta')
models_w40 <- widthGam(w40)

# ------------------------------------ SAVE PLOT ---------------------------------------%

jpeg('figures/figure_s3.jpeg',
     width = 10*300, height = 6*300, res = 300)

layout(matrix(1:6, 
              nrow = 2, ncol = 3, 
              byrow = TRUE),
       widths = rep(1.5, 3), heights = rep(1,3))

# width 1:10
sapply(1:3, function(panel) widthGamDataBirthAdjustedPlot(stat = w10, 
                                                          m = models_w10, 
                                                          panel = panel, 
                                                          cutoff = 1, 
                                                          spp = TRUE,
                                                          legend = LETTERS[1:3]))

# width 1:40
sapply(1:3, function(panel) widthGamDataBirthAdjustedPlot(stat = w40, 
                                                          m = models_w40, 
                                                          panel = panel, 
                                                          cutoff = 3, 
                                                          spp = TRUE,
                                                          legend = LETTERS[4:6]))


dev.off()