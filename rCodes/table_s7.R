#
# table s7
#

library(Rarity)
library(dplyr)
library(mgcv)

# ------------------------------ LOAD DATA AND FUNCTIONS ---------------------------------------%

load('data/data.RData')
source('rCodes/statFuncs.R')

# ------------------------------------ FUNCTIONS --------------------------------------------%

##' 
##' Function to fit univariable logistic regression of seroconversion 
##' to four recent strains
##' 
##' @param out list of dataset for seroconversion, returned from getSeroconversionData()
##' 
##' @return table s7
##' 

summaryUnivariableEachother <- function(out){
  
  dgt <- 2
  
  tab <- matrix(nrow = 5*4 + 4, ncol = 5)
  colnames(tab) <-c("Age at sampling",
                    "Titer for strain i",
                    "Titer for strain i - 1",
                    "AUC",
                    "Width, cut off 1:40")
  
  strains <- c(sapply(1:4, function(i) out[[i]]$outcome_strain[1]))
  
  rownames(tab) <- c(strains[1],
                     colnames(tab),
                     strains[2],
                     colnames(tab),
                     strains[3],
                     colnames(tab),
                     strains[4],
                     colnames(tab))
  
  terms <- c('age_at_sampling',
             'titer_outcome_strain_v1',
             'titer_previous_strain_v1',
             'AUC_from_birth_adj',
             'W40_from_birth_adj')
  
  for(i in 1:4){ # strain
    
    for(outcome in 1:5){
      
      r <- (i-1)*6 + outcome + 1
      tab[r, outcome] <- "--"
      
      for(predictor in c(1:5)[-outcome]){
        
        fml <- as.formula(paste(terms[outcome], 
                                terms[predictor], 
                                sep = ' ~ '))
        
        m <- lm(fml, data = out[[i]])
        
        coeff <- round(coef(m)[2], dgt)
        ci <- round(confint(m)[2, ], dgt)
        
        p <- summary(m)$coefficients[2, 4]
        star <- ""
        star[p < 0.05] <- "*"
        
        r <- (i-1)*6 + predictor + 1
        tab[r, outcome] <- paste(sprintf("%.02f", coeff),
                                 " (", 
                                 sprintf("%.02f", ci[1]),
                                 ", ", 
                                 sprintf("%.02f", ci[2]),
                                 ")", star, sep = "")
        
      }
      
    }
    
  }
  
  tab[is.na(tab)] <- ""
  
  return(tab)
  
}


# ----------------------------------- SAVE TABLE ---------------------------------------%

# generate data used for logistic regressions
out <- getSeroconversionData(conversion, hi_paired, st)

# results shown in table s10
table_s7 <- summaryUnivariableEachother(out)
write.csv(table_s7, 'tables/table_s7.csv', row.names = TRUE)

