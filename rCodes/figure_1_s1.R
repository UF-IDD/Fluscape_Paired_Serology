#
# figure 1 & s1
#

library(dplyr)

# ------------------------------------ LOAD DATA ---------------------------------------%

load('data/data.RData')
source('rCodes/figureFuncs.R')
source('rCodes/statFuncs.R')

# ------------------------------------ PREPARE DATA ---------------------------------------%

paired <- fig1Data(hi_paired)
stat <- birthAdjustedMetrics(hi_paired, conversion)


# ------------------------------------ SAVE PLOT ---------------------------------------%

#ids <- selID(paired)

#These are randomly selected ids using the above function
ids <- c("L42H20P02", "L50H15P03", "L03H19P03", 
         "L45H01P01", "L27H16P04", "L24H10P01")

jpeg('figures/figure_1.jpeg',
      width = 20*300, height = 11.5*300, res = 300)

layout(matrix(1:6, 
              nrow = 2, ncol = 3, 
              byrow = TRUE),
       widths = rep(1.5, 3), 
       heights = rep(1, 2))

mat <- rbind(matrix(1:6, 
                    nrow = 2, ncol = 3, 
                    byrow = FALSE),
             matrix(7:12, 
                    nrow = 2, ncol = 3, 
                    byrow = FALSE))
layout(mat,
       widths = rep(1.5, 3), 
       heights = c(1, 0.4, 1, 0.4))

sapply(1:length(ids), function(i) {
  
  individualPlot(i = i, 
                 ids = ids, 
                 paired = paired)
  individualMetricValuePlot(i = i,
                            ids = ids,
                            stat = stat)
  
})

dev.off()


# # random generation
# ids.midage <- c(sample(stat$participant_id[stat$age_at_sampling %in% c(41:50)],
#                        3, replace = FALSE),
#                 sample(stat$participant_id[stat$age_at_sampling %in% c(51:60)],
#                        3, replace = FALSE))

ids.midage <- c("L22H16P03", "L20H22P01", "L21H10P01",
                "L32H09P01", "L15H08P01", "L08H10P02")


jpeg('figures/figure_s1.jpeg',
     width = 20*300, height = 11.5*300, res = 300)

layout(matrix(1:6, 
              nrow = 2, ncol = 3, 
              byrow = TRUE),
       widths = rep(1.5, 3), 
       heights = rep(1, 2))

mat <- rbind(matrix(1:6, 
                    nrow = 2, ncol = 3, 
                    byrow = FALSE),
             matrix(7:12, 
                    nrow = 2, ncol = 3, 
                    byrow = FALSE))
layout(mat,
       widths = rep(1.5, 3), 
       heights = c(1, 0.4, 1, 0.4))

sapply(1:length(ids.midage), function(i) {
  
  individualPlot(i = i, 
                 ids = ids.midage, 
                 paired = paired)
  individualMetricValuePlot(i = i,
                            ids = ids.midage,
                            stat = stat)
  
})

dev.off()