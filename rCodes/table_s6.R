#
# table s6
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
##' @return table s6
##' 

summaryUnivariableRegression <- function(out){
  
  dgt <- 2
  
  tab <- matrix(nrow = 3 + 4*2, ncol = 4)
  rownames(tab) <- c("Age at sampling",
                     "Titer for strain i",
                     "Titer for strain i - 1",
                     "AUC, exposed strains",
                     "ATY, exposed strains",
                     "W10, exposed strains",
                     "W40, exposed strains",
                     "AUC, all strains",
                     "ATY, all strains",
                     "W10, all strains",
                     "W40, all strains")
  colnames(tab) <- c(sapply(1:4, function(i) out[[i]]$outcome_strain[1]))
  
  for(i in 1:4){ # strain
    
    # age
    base <- as.formula(seroconversion ~ age_at_sampling)
    m <- glm(base, data = out[[i]], family = binomial())
    
    coeff <- round(exp(coef(m)[2]), dgt)
    ci <- round(exp(confint(m)[2, ]), dgt)
    
    p <- summary(m)$coefficients[2, 4]
    star <- ""
    star[p < 0.05] <- "*"
    
    tab[1, i] <- paste(sprintf("%.02f", coeff),
                       " (", 
                       sprintf("%.02f", ci[1]),
                       ", ", 
                       sprintf("%.02f", ci[2]),
                       ")", star, sep = "")
    
    # titer i
    base <- as.formula(seroconversion ~ titer_outcome_strain_v1)
    m <- glm(base, data = out[[i]], family = binomial())
    
    coeff <- round(exp(coef(m)[2]), dgt)
    ci <- round(exp(confint(m)[2, ]), dgt)
    
    p <- summary(m)$coefficients[2, 4]
    star <- ""
    star[p < 0.05] <- "*"
    
    tab[2, i] <- paste(sprintf("%.02f", coeff),
                       " (", 
                       sprintf("%.02f", ci[1]),
                       ", ", 
                       sprintf("%.02f", ci[2]),
                       ")", star, sep = "")
    
    # titer i-1
    base <- as.formula(seroconversion ~ titer_previous_strain_v1)
    m <- glm(base, data = out[[i]], family = binomial())
    
    coeff <- round(exp(coef(m)[2]), dgt)
    ci <- round(exp(confint(m)[2, ]), dgt)
    
    p <- summary(m)$coefficients[2, 4]
    star <- ""
    star[p < 0.05] <- "*"
    
    tab[3, i] <- paste(sprintf("%.02f", coeff),
                       " (", 
                       sprintf("%.02f", ci[1]),
                       ", ", 
                       sprintf("%.02f", ci[2]),
                       ")", star, sep = "")
    
    base <- "seroconversion ~ "
    
    for(range in 1:2){
      
      for(k in 1:4){
        
        metric <- c("AUC", "ATY", "W10", "W40")[k]
        rg <- c('1968', 'birth', 'adj', 'recent')[c(3,1)[range]]
        
        used <- colnames(out[[i]])[grepl(metric, colnames(out[[i]]))]
        used <- used[!(grepl("residual", used))]
        used <- used[grepl(rg, used)]
        
        
        fml <- as.formula(paste(base, used, sep = " "))
        m <- glm(fml, data = out[[i]], family = binomial())
        
        coeff <- round(exp(coef(m)[2]), dgt)
        ci <- round(exp(confint(m)[2, ]), dgt)
        
        p <- summary(m)$coefficients[2, 4]
        star <- ""
        star[p < 0.05] <- "*"
        
        r <- (range-1)*4 + 3 + k
        tab[r, i] <- paste(sprintf("%.02f", coeff),
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

# results shown in table s11
table_s6 <- summaryUnivariableRegression(out)
write.csv(table_s6, 'tables/table_s6.csv', row.names = TRUE)

