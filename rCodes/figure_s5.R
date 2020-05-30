#
# figure s5
#

library(tidyr)

# ------------------------------------ LOAD DATA ---------------------------------------%

load('data/data.RData')
source('rCodes/statFuncs.R')
source('rCodes/figureFuncs.R')

# ------------------------------------ PREPARE DATA ---------------------------------------%

changes <- getDataForTiterChanges(conversion)

# ------------------------------------ SAVE PLOT ---------------------------------------%

jpeg('figures/figure_s5.jpeg',
     width = 12*450, height = 6*450, res = 450)

layout(matrix(c(1, 2), 
              nrow = 1, ncol = 2, 
              byrow = TRUE), 
       widths = rep(1, 2), 
       heights = rep(1, 2))

# changes in four strains
plotChangesInFourStrains(changes, byStrain = FALSE)
plotChangesInFourStrains(changes, byStrain = TRUE)

dev.off()
