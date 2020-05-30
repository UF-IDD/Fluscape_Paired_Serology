#
# figure S12
#

library(dplyr)

# ------------------------------------ LOAD DATA ---------------------------------------%

load('data/data.RData')
source('rCodes/statFuncs.R')
source('rCodes/figureFuncs.R')


# ---------------------------------- PREPARE PLOT ---------------------------------------%

uniData <- getDataForUniTiter(conversion[conversion$age_at_isolation > 0, ], 
                              hi_paired, st)

uniEstAdj <- sapply(1:4, function(s) getStrainCoeffForSeroconv(data = uniData, 
                                                               strain = s,
                                                               includeStraini = TRUE))

uniEst <- sapply(1:4, function(s) getStrainCoeffForSeroconv(data = uniData, 
                                                            strain = s,
                                                            includeStraini = FALSE))

# ------------------------------------ SAVE PLOT ---------------------------------------%

jpeg('figures/figure_s12.jpeg',
     width = 12*300, height = 12*300, res = 300)

layout(matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE),
       widths = rep(1, 4), heights = rep(1.5, 2))

sapply(1:4, function(s) ORSingleStrainPlot(uniEst, 
                                           uniEstAdj, 
                                           strain = s, 
                                           st))

dev.off()