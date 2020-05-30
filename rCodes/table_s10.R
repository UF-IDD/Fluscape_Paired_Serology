#
# table s9
#

library(Rarity)
library(dplyr)
library(mgcv)

# ------------------------------ LOAD DATA AND FUNCTIONS ---------------------------------------%

load('/Volumes/GoogleDrive/My Drive/R projects/Fluscape/r codes/spatial_pattern_serology/Fluscape_serology/data/data.RData')
source('/Volumes/GoogleDrive/My Drive/R projects/Fluscape/r codes/spatial_pattern_serology/Fluscape_serology//rCodes/statFuncs.R')
source("../../source/R/GeneralUtility.r")
source("/Volumes/GoogleDrive/My Drive/R projects/Fluscape/r codes/spatial_pattern_serology/Fluscape_serology/rCodes/GeneralUtility_BY.r")


# ------------------ FUNCTIONS TO CLEAN SAMPLE TIME AND REGRESSION ----------------------------%

getSampleTime <- function(ids){
  
  setwd('/Users/bingyiyang/Documents/Fluscape/fluscape/trunk/manuscripts/paired_titer')
  
  topdir <- "../../"
  
  hh <- load.and.merge.household.data.V1.V2.V3.V4(missing.to.na = TRUE,
                                                  topdir = topdir,
                                                  date.format = 3)
  
  hh.short <- hhVisitDate(topdir)
  
  sample.time <- t(sapply(ids, function(i){
    
    l <- substr(i, 1, 3)
    h <- substr(i, 4, 6)
    
    hh.short[hh.short$HH_ID == h &
               hh.short$LOC_ID == l &
               hh.short$Visit %in% c('V1', 'V4'), 
             'HH_VisitA1DateTime']
    
  }))
  
  sample.time <- data.frame(id = ids,
                            v1 = sample.time[ ,1],
                            v4 = sample.time[ ,2],
                            stringsAsFactors = FALSE)
  
  sample.time[ ,2] <- as.Date(sample.time[ ,2],
                              '%Y-%m-%d')
  sample.time[ ,3] <- as.Date(sample.time[ ,3],
                              '%Y-%m-%d')
  
  sample.time$v1[sample.time$v1 == as.Date('2002-11-02',
                                           '%Y-%m-%d')] <- as.Date('2010-11-02',
                                                                   '%Y-%m-%d')
  
  sample.time$v1.day <- as.numeric(sample.time$v1 - as.Date('2009-01-01', '%Y-%m-%d'))
  sample.time$v4.day <- as.numeric(sample.time$v4 - as.Date('2014-01-01', '%Y-%m-%d'))
  
  sample.time
  
}

summaryRegressionFull <- function(out, range){
  
  dgt <- 2
  
  tab <- matrix(nrow = 5 + 4 + 5*4, ncol = 4)
  rownames(tab) <- c("Model 1",
                     "Age at sampling",
                     "Titer for strain i",
                     "Titer for strain i - 1",
                     "% of deviance explained",
                     "Model 2",
                     "Age at sampling",
                     "Titer for strain i",
                     "Titer for strain i - 1",
                     "AUC",
                     "% of deviance explained",
                     "Model 3",
                     "Age at sampling",
                     "Titer for strain i",
                     "Titer for strain i - 1",
                     "ATY",
                     "% of deviance explained",
                     "Model 4",
                     "Age at sampling",
                     "Titer for strain i",
                     "Titer for strain i - 1",
                     "Width, cut off at 1:10",
                     "% of deviance explained",
                     "Model 5",
                     "Age at sampling",
                     "Titer for strain i",
                     "Titer for strain i - 1",
                     "Width, cut off at 1:40",
                     "% of deviance explained")
  
  colnames(tab) <- c(sapply(1:4, function(i) out[[i]]$outcome_strain[1]))
  
  for(i in 1:4){ # strain
    
    # null model
    m0 <- glm(seroconversion ~ 1, data = out[[i]], family = binomial())
    d0 <- deviance(m0)
    
    # model 1
    base <- as.formula(seroconversion ~ age_at_sampling + titer_outcome_strain_v1 + 
                         titer_previous_strain_v1 + s(v1.day, v4.day))
    m <- gam(base, data = out[[i]], family = binomial())
    
    
    coeff <- coef(m)[2:4]
    se <- summary(m)$se[2:4]
    ci <- cbind(coeff - 1.96*se,
                coeff + 1.96*se)
    
    coeff <- round(exp(coeff), dgt)
    ci <- round(exp(ci), dgt)
    
    p <- summary(m)$p.pv[2:4]
    star <- rep("", length(coeff))
    star[p < 0.05] <- "*"
    
    tab[2:4, i] <- paste(sprintf("%.02f", coeff),
                         " (", 
                         sprintf("%.02f", ci[ ,1]),
                         ", ", 
                         sprintf("%.02f", ci[ ,2]),
                         ")", star, sep = "")
    
    pDev <- (d0 - deviance(m))/d0*100
    tab[5, i] <- sprintf("%.01f", pDev)
    
    base <- "seroconversion ~ age_at_sampling + titer_outcome_strain_v1 + 
    titer_previous_strain_v1 + s(v1.day, v4.day)"
    
    for(k in 1:4){ # summary metric
      
      metric <- c("AUC", "ATY", "W10", "W40")[k]
      rg <- c('1968', 'adj')[range]
      
      used <- colnames(out[[i]])[grepl(metric, colnames(out[[i]]))]
      used <- used[grepl(rg, used)]
      
      fml <- as.formula(paste(base, used, sep = " + "))
      m <- gam(fml, data = out[[i]], family = binomial())

      
      coeff <- coef(m)[2:5]
      se <- summary(m)$se[2:5]
      ci <- cbind(coeff - 1.96*se,
                  coeff + 1.96*se)
      
      coeff <- round(exp(coeff), dgt)
      ci <- round(exp(ci), dgt)
      
      p <- summary(m)$p.pv[2:5]
      star <- rep("", length(coeff))
      star[p < 0.05] <- "*"
      
      r1 <- (k-1)*6 + 7
      r2 <- (k-1)*6 + 10
      tab[r1:r2, i] <- paste(sprintf("%.02f", coeff),
                             " (", 
                             sprintf("%.02f", ci[ ,1]),
                             ", ", 
                             sprintf("%.02f", ci[ ,2]),
                             ")", star, sep = "")
      
      pDev <- (d0 - deviance(m))/d0*100
      tab[r2 + 1, i] <- sprintf("%.01f", pDev)
      
      
    }
    
    
  }
  
  tab[is.na(tab)] <- ""
  
  return(tab)
  
}


# ------------------------------------ ANALYSIS --------------------------------------------%

# generate data used for logistic regressions
out <- getSeroconversionData(conversion, hi_paired, st)

# get sample time
ids <- unique(conversion$participant_id)
sample.time <- getSampleTime(ids)

# merge the two data sets
for(i in 1:length(out)){
  
  out[[i]] <- left_join(out[[i]], sample.time)
  
}

# ----------------------------------- SAVE TABLE ---------------------------------------%

# results shown in table 3 and s4
table_sx <- summaryRegressionFull(out, range = 2)
write.csv(table_sx, '/Volumes/GoogleDrive/My Drive/R projects/Fluscape/r codes/spatial_pattern_serology/Fluscape_serology/tables/table_sx.csv', row.names = TRUE)


d <- seq(as.Date('2009-01-01',
                 '%Y-%m-%d'),
         as.Date('2015-12-31',
                 '%Y-%m-%d'),
         1)
x <- sapply(d, function(x) sum(sample.time$v1 == x,
                          na.rm = TRUE) + 
         sum(sample.time$v4 == x,
             na.rm = TRUE))
