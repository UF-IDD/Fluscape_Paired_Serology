#
# figure 5
#

library(tidyr)

# ------------------------------------ LOAD DATA ---------------------------------------%

load('data/data.RData')
source('rCodes/statFuncs.R')
source('rCodes/figureFuncs.R')


# ------------------------------------ SAVE PLOT ---------------------------------------%

jpeg('figures/figure_5.jpeg',
     width = 14*450, height = 10*450, res = 450)

layout(matrix(c(1, 2), 
              nrow = 2, ncol = 1), 
       widths = c(4, 4), 
       heights = c(1.5, 3))

# gmt
propChangePlot(conversion, st, all = TRUE, panel = 1 , supp = FALSE)

# distribution of strains by change
strainDistrByTiterChanges(data = conversion, 
                          st = st,
                          panel = 2)

dev.off() 

