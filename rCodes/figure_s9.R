#
# figure S9
#

# ------------------------------------ LOAD DATA ---------------------------------------%

load('data/data.RData')
source('rCodes/statFuncs.R')
source('rCodes/figureFuncs.R')


# ---------------------------------- PREPARE DATA --------------------------------------%

ageDistr <- getAgeTiterChanges(conversion)
ci <- boostrap(ageDistr, n = 2000)

# ------------------------------------ SAVE PLOT ---------------------------------------%

jpeg('figures/figure_s9.jpeg',
     width = 8*300, height = 8*300, res = 300)

nr <- 2; nc <- 2
layout(matrix(1:(nr*nc), nrow = nr, ncol = nc, byrow = TRUE),
       widths = rep(1, nc), heights = rep(1, nr))

sapply(1:(nc*nr), function(panel) histPlot(stat = ageDistr, 
                                      panel = panel))

dev.off()