#
# table 1, s3 and s4
#

library(Rarity)
library(dplyr)
library(mgcv)

# ------------------------------ LOAD DATA AND FUNCTIONS ---------------------------------------%

load('data/data.RData')
source('rCodes/statFuncs.R')

# ------------------------------------ ANALYSIS --------------------------------------------%

# generate data used for logistic regressions
out <- getSeroconversionData(conversion, hi_paired, st)


# ----------------------------------- SAVE TABLE ---------------------------------------%

# results shown in table 3 and s4
table_1 <- summaryRegressionFull(out, range = 2)
write.csv(table_1, 'tables/table_1_s4.csv', row.names = TRUE)

# results shown in table s3 and s4
table_s3 <- summaryRegressionFull(out, range = 1)
write.csv(table_s3, 'tables/table_s3_s4.csv', row.names = TRUE)
