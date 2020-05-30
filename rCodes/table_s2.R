#
# table s2
#

# ------------------------------------ LOAD DATA ---------------------------------------%

load('data/data.RData')

# ------------------------------------ FUNCTION ----------------------------------------%

getGMTTab <- function(conversion, st){
  
  tab <- data.frame(strain = st$strain,
                    allBaseline = NA,
                    allFollowup = NA,
                    postBirthBaseline = NA,
                    postBirthFollowup = NA,
                    preBirthBaseline = NA,
                    preBirthFollowup = NA,
                    stringsAsFactors = FALSE)
  
  tab2 <- data.frame(strain = 'Overall',
                     allBaseline = NA,
                     allFollowup = NA,
                     postBirthBaseline = NA,
                     postBirthFollowup = NA,
                     preBirthBaseline = NA,
                     preBirthFollowup = NA,
                     stringsAsFactors = FALSE)
  
  
  for(i in 1:3){# all strains, post-birth strains or pre-birth strains
    
    dat <- conversion
    
    if(i == 2){
      dat <- conversion[conversion$age_at_isolation >= 0, ]
    }
    if(i == 3){
      dat <- conversion[conversion$age_at_isolation < 0, ]
    }
    
    for(j in 1:2){# baseline/follow-up
      
      col.used <- c('log_hi_baseline',
                    'log_hi_followup')[j]
      
      # strain-specific
      gmt <- sapply(st$strain, function(s){
        
        
        tem <- dat[dat$strain == s & !is.na(dat[ ,col.used]), 
                   col.used]
        
        if(length(tem)>1){
          est <- t.test(tem)
          
          a <- c(est$estimate, est$conf.int)
          round(2^a*5, 1)
          
        } else {rep(NA, 3)}
        
      })
      
      tab[ ,(i-1)*2+j+1] <- paste(sprintf('%.01f', gmt[1, ]),
                                  ' (', sprintf('%.01f', gmt[2, ]), 
                                  ', ', sprintf('%.01f', gmt[3, ]), ')',
                                  sep = '')
      
      # overall gmt
      tem <- dat[!is.na(dat[ ,col.used]), 
                 col.used]
      
      est <- t.test(tem)
      
      a <- c(est$estimate, est$conf.int)
      gmt <- round(2^a*5, 1)
      
      tab2[ ,(i-1)*2+j+1] <- paste(sprintf('%.01f', gmt[1]),
                                   ' (', sprintf('%.01f', gmt[2]), 
                                   ', ', sprintf('%.01f', gmt[3]), ')',
                                   sep = '')
      
    }
    
  }
  
  
  tab <- rbind(tab, tab2)
  
  return(tab)
  
}

# ----------------------------------- SAVE TABLE ---------------------------------------%

table_s2 <- getGMTTab(conversion, st)

write.csv(table_s2, 'tables/table_s2.csv', row.names = TRUE)
