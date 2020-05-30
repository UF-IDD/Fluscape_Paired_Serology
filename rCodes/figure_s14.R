#
# figure S14
#

# ------------------------------------ LOAD DATA ---------------------------------------%

load('data/data.RData')
source('rCodes/figureFuncs.R')

# ------------------------------------ SAVE PLOT ---------------------------------------%

jpeg('figures/figure_s14.jpeg',
     width = 18*300, height = 16*300, res = 300)

layout(matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE), widths = rep(1.3,2), heights = rep(1,2))
strainDistrByTiterChangesCeiling(data = conversion, st = st,
                                 inset = TRUE, 
                                 panel = 1,
                                 legend = FALSE)
strainDistrByTiterChangesCeiling(data = conversion[conversion$log_hi_baseline <= 3, ], 
                                 st = st,
                                 inset = TRUE, 
                                 panel = 2,
                                 legend = FALSE)
strainDistrByTiterChangesCeiling(data = conversion[conversion$log_hi_baseline > 3, ], 
                                 st = st,
                                 inset = TRUE, 
                                 panel = 3,
                                 legend = FALSE)

dev.off()