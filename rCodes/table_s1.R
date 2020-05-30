#
# table s1
#

# ------------------------------------ LOAD DATA ---------------------------------------%

load('data/data.RData')

participant <- readRDS('data/participant.rds')

# ------------------------------------ FUNCTION ---------------------------------------%

generateTableS1 <- function(hi_paired, participant){
  
  part <- list()
  part[[1]] <- participant$Blood_available
  part[[2]] <- participant$Baseline[!(participant$Baseline$participant_id %in% participant$Blood_available$participant_id), ]
  part[[3]] <- participant$`Follow-up`[!(participant$`Follow-up`$participant_id %in% participant$Blood_available$participant_id), ]
  part[[3]] <- part[[3]][!is.na(part[[3]]$sex), ]
  
  tab <- data.frame(matrix(nrow = 27,
                           ncol = 6))
  
  colnames(tab) <- c("Serum Available",
                     "Unavailable for baseline",
                     "P-value",
                     "Serum Available",
                     "Unavailable for follow-up",
                     "P-value")
  
  rownames(tab) <- c("N",
                     "Sex", "Male", "Female",
                     "Age group, yrs", "< 10", "10-19", "20-29", 
                     "30-39", "40-49", "50-59", "60 +",
                     "Employment status", 
                     "Full Time", "Self Employed", "Retired", 
                     "Student", "Homemaker", "unemployed", "other",
                     "Years since last influenza vaccination",
                     "Less than 1", "1", "2-5", "More than 5",
                     "Never", "Unknown/unsure")
  
  for(i in 1:2){
    
    df <- list(part[[1]],
               part[[i + 1]])
    
    col.used <- ((i-1)*3+1):(i*3)
    
    for(j in 1:2){
      
      # total # of participants
      tab[1, col.used[j]] <- nrow(df[[j]])
      
      # sex
      tab[3:4, col.used[j]] <- table(df[[j]]$sex)
      
      # age group
      df[[j]]$ageGroup <- floor(df[[j]]$age_at_sampling/10)
      df[[j]]$ageGroup[df[[j]]$ageGroup>=6 & !is.na(df[[j]]$ageGroup)] <- 6
      
      tab[6:12, col.used[j]] <- table(df[[j]]$ageGroup)
      
      # employment status
      for(k in 1:6){
        tab[k+13, col.used[j]] <- sum(df[[j]]$occupation == c(1,2,4,5,6,8)[k], na.rm = TRUE)
      }
      tab[20, col.used[j]] <- nrow(df[[j]]) - sum(tab[14:19, col.used[j]])
      
      # vaccine group
      tab[22:26, col.used[j]] <- table(df[[j]]$last_flu_vaccine)[c(2:5, 1)]
      tab[27, col.used[j]] <- sum(table(df[[j]]$last_flu_vaccine)[c(6:7)], na.rm = TRUE) + sum(is.na(df[[j]]$last_flu_vaccine))
      
    }
    
    # chisq.test
    # sex
    tab[3, col.used[3]] <- sprintf('%.02f', chisq.test(tab[3:4, col.used[1:2]])$p.value)
    
    # age
    tab[6, col.used[3]] <- sprintf('%.02f', chisq.test(tab[6:12, col.used[1:2]])$p.value)
    
    # employment status
    tab[14, col.used[3]] <- sprintf('%.02f', chisq.test(tab[14:20, col.used[1:2]])$p.value)
    
    # vaccine group
    if(i == 1){
      tab[22, col.used[3]] <- sprintf('%.02f', chisq.test(tab[22:27, col.used[1:2]])$p.value)
    } else {
      tab[22, col.used[3]] <- sprintf('%.02f', chisq.test(tab[22:27, col.used[1:2]])$p.value)
    }
    
    
  }
  
  
  # add percentage
  ps <- sapply(c(1,2,4,5), function(i) round(tab[2:27,i]/tab[1, i]*100, 1))
  
  for(i in 1:4){
    
    tab[2:27 ,c(1,2,4,5)[i]] <- paste(tab[2:27 ,c(1,2,4,5)[i]],
                                " (",
                                sprintf("%.01f", ps[ ,i]),
                                ")",
                                sep = "")
    
                          
  }
  
  tab[tab=="NA (NA)"] <- ""
  tab[is.na(tab)] <- ""
  tab[tab == '0.00'] <- '<0.01'
  
  tab <- tab[ ,-4]
  
  return(tab)
  
}

table_s1 <- generateTableS1(hi_paired, participant)

write.csv(table_s1, 'tables/table_s1.csv', row.names = TRUE)