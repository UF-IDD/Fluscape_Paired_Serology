#
# table s5
#

library(Rarity)
library(dplyr)
library(mgcv)

# ------------------------------ LOAD DATA AND FUNCTIONS ---------------------------------------%

load('data/data.RData')
source('rCodes/statFuncs.R')

# ------------------------------------ FUNCTIONS --------------------------------------------%

##' 
##' Function to fit logistic regressions of seroconversion to four recent strains
##' spline term on age
##' 
##' @param out list of dataset for seroconversion, returned from getSeroconversionData()
##' @param range 1 = all strains
##'              2 = post-birth strains
##' 
##' 
##' @return table s5
##' 

summaryRegressionSplineDeviance <- function(out, range){
  
  dgt <- 2
  
  tab <- matrix(nrow = 5 + 3 + 4*4, ncol = 4)
  rownames(tab) <- c("Model 1",
                     "Titer for strain i",
                     "Titer for strain i - 1",
                     "% of deviance explained",
                     "Model 2",
                     "Titer for strain i",
                     "Titer for strain i - 1",
                     "AUC",
                     "% of deviance explained",
                     "Model 3",
                     "Titer for strain i",
                     "Titer for strain i - 1",
                     "ATY", 
                     "% of deviance explained",
                     "Model 4",
                     "Titer for strain i",
                     "Titer for strain i - 1",
                     "Width, cut off at 1:10",
                     "% of deviance explained",
                     "Model 5",
                     "Titer for strain i",
                     "Titer for strain i - 1",
                     "Width, cut off at 1:40",
                     "% of deviance explained")
  
  colnames(tab) <- c(sapply(1:4, function(i) out[[i]]$outcome_strain[1]))
  
  for(i in 1:4){ # strain
    
    # null model
    m0 <- gam(seroconversion ~ 1, data = out[[i]], family = binomial())
    d0 <- deviance(m0)
    
    # model 1
    base <- as.formula(seroconversion ~ s(age_at_sampling) + titer_outcome_strain_v1 + 
                         titer_previous_strain_v1)
    m <- gam(base, data = out[[i]], family = binomial())
    
    coeff <- coef(m)[2:3]
    se <- summary(m)$se[2:3]
    ci <- cbind(coeff - 1.96*se,
                coeff + 1.96*se)
    
    coeff <- round(exp(coeff), dgt)
    ci <- round(exp(ci), dgt)
    
    p <- summary(m)$p.pv[2:3]
    star <- rep("", 2)
    star[p < 0.05] <- "*"
    
    tab[2:3, i] <- paste(sprintf("%.02f", coeff),
                         " (", 
                         sprintf("%.02f", ci[ ,1]),
                         ", ", 
                         sprintf("%.02f", ci[ ,2]),
                         ")", star, sep = "")
    
    pDev <- (d0 - deviance(m))/d0*100
    tab[4, i] <- sprintf("%.01f", pDev)
    
    base <- "seroconversion ~ s(age_at_sampling) + titer_outcome_strain_v1 + 
    titer_previous_strain_v1"
    
    for(k in 1:4){
      
      metric <- c("AUC", "ATY", "W10", "W40")[k]
      rg <- c('1968', 'adj')[range]
      
      used <- colnames(out[[i]])[grepl(metric, colnames(out[[i]]))]
      used <- used[grepl(rg, used)]
      
      fml <- as.formula(paste(base, used, sep = " + "))
      m <- gam(fml, data = out[[i]], family = binomial())
      
      coeff <- coef(m)[2:4]
      se <- summary(m)$se[2:4]
      ci <- cbind(coeff - 1.96*se,
                  coeff + 1.96*se)
      
      coeff <- round(exp(coeff), dgt)
      ci <- round(exp(ci), dgt)
      
      p <- summary(m)$p.pv[2:4]
      star <- rep("", 3)
      star[p < 0.05] <- "*"
      
      r1 <- (k-1)*5+6
      r2 <- (k-1)*5+8
      tab[r1:r2, i] <- paste(sprintf("%.02f", coeff),
                             " (", 
                             sprintf("%.02f", ci[ ,1]),
                             ", ", 
                             sprintf("%.02f", ci[ ,2]),
                             ")", star, sep = "")
      
      pDev <- (d0 - deviance(m))/d0*100
      tab[r2+1, i] <- sprintf("%.01f", pDev)
      
    }
    
    
  }
  
  tab[is.na(tab)] <- ""
  
  return(tab)
  
  }


# ----------------------------------- SAVE TABLE ---------------------------------------%

# generate data used for logistic regressions
out <- getSeroconversionData(conversion, hi_paired, st)

# results shown in table s4
table_s5 <- summaryRegressionSplineDeviance(out, range = 2)
write.csv(table_s5, 'tables/table_s5.csv', row.names = TRUE)

