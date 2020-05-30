#
# table s9
#

library(Rarity)
library(dplyr)
library(mgcv)

# ------------------------------ LOAD DATA AND FUNCTIONS ---------------------------------------%

load('data/data_never_vaccinated.RData')
source('rCodes/statFuncs.R')

# ------------------------------------ ANALYSIS --------------------------------------------%

# generate data used for logistic regressions
out <- getSeroconversionData(conversion, hi_paired, st)


# ----------------------------------- SAVE TABLE ---------------------------------------%

# results shown in table 14
table_s9 <- summaryRegressionFull(out, range = 2)
write.csv(table_s9, 'tables/table_s9.csv', row.names = TRUE)