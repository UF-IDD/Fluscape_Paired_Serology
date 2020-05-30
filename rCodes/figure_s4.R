#
# figure S4
#

# ------------------------------------ LOAD DATA ---------------------------------------%

load('data/data.RData')
source('rCodes/statFuncs.R')
source('rCodes/figureFuncs.R')

# ---------------------------------- PREPARE DATA ---------------------------------------%

data <-  conversion
data$strain <- factor(data$strain, 
                      levels = st$strain)
# fit logistic regression
m <- glm(seroconversion ~ age_baseline + log_hi_baseline + strain,
         family = binomial, data = data)

# ------------------------------------ SAVE PLOT ---------------------------------------%

jpeg('figures/figure_s4.jpeg',
     width = 8*300, height = 5*300, res = 300)

seroconvORPlot(m)

dev.off()