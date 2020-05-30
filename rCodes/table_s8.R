#
# table s8
#

library(Rarity)
library(tidyr)
library(mgcv)

# ------------------------------ LOAD DATA AND FUNCTIONS ---------------------------------------%

load('data/data.RData')
conversion_all <- conversion

load('data/data_never_vaccinated.RData')
conversion_no_vacc <- conversion

participant <- readRDS('data/participant.rds')

# ------------------------------------ FUNCTIONS --------------------------------------------%

##' 
##' Function to generate table to compare demographic characteristics
##' of participants who self-reported to have not been vaccinated against influenza
##' 
##' @param conversion for all participants
##' @param conversion_no_vacc for participants reported never vaccinated against influenza
##' @param part, list of participants 
##' 
##' @return table s13
##' 

compTable <- function(conversion, conversion_no_vacc, part){
  
  ids <- list(unique(conversion$participant_id),
              unique(conversion_no_vacc$participant_id))
  
  tab <- data.frame(matrix(nrow = 20,
                           ncol = 3))
  
  colnames(tab) <- c("All",
                     "Never vaccinated",
                     "P")
  
  rownames(tab) <- c("N",
                     "Sex", "Male", "Female",
                     "Age group, yrs", "< 10", "10-19", "20-29", 
                     "30-39", "40-49", "50-59", "60 +",
                     "Employment status", 
                     "Full Time", "Self Employed", "Retired", 
                     "Student", "Homemaker", "unemployed", "other")
  
  for(i in 1:2){# all or never vaccinated
    
    df <- part$Blood_available
    df <- df[df$participant_id %in% ids[[i]], ]
    
    # total # of participants
    tab[1, i] <- length(ids[[i]])
    
    # sex
    tab[3:4, i] <- table(df$sex)
    
    # age group
    df$ageGroup <- floor(df$age_at_sampling/10)
    df$ageGroup[df$ageGroup>=6 & !is.na(df$ageGroup)] <- 6
    
    tab[6:12, i] <- table(df$ageGroup)
    
    # employment status
    tab[14:19, i] <- table(df$occupation)[c(1,2,4,5,6,8)]
    tab[20, i] <- sum(table(df$occupation)[c(3,7,9,10)])
    
  }
  
  # chi-sq test
  p <- rep(NA, 3)
  
  tst <- chisq.test(tab[3:4, 1:2])
  p[1] <- sprintf("%.02f", tst$p.value)
  
  tst <- chisq.test(tab[6:12, 1:2])
  p[2] <- sprintf("%.02f", tst$p.value)
  
  tst <- chisq.test(tab[14:20, 1:2])
  p[3] <- sprintf("%.02f", tst$p.value)
  
  p[p == '0.00'] <- '<0.01'
  
  # add percentage
  ps <- sapply(1:2, function(i) round(tab[2:20,i]/tab[1, i]*100, 1))
  
  for(i in 1:2){
    
    tab[2:20 ,i] <- paste(tab[2:20 ,i],
                          " (",
                          sprintf("%.01f", ps[,i]),
                          ")",
                          sep = "")
    
  }
  
  tab[c(3, 6, 14), 3] <- p
  
  tab[tab=="NA (NA)"] <- ""
  tab[is.na(tab)] <- ""
  
  return(tab)
  
  
}


# ----------------------------------- SAVE TABLE ---------------------------------------%

table_s8 <- compTable(conversion = conversion_all, 
                      conversion_no_vacc = conversion_no_vacc,
                      part = participant)
write.csv(table_s8, 'tables/table_s8.csv', row.names = TRUE)


