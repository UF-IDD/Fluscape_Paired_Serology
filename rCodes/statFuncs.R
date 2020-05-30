#
# 2019-10-16
# Functions for summary statistics
#

##' Function to fit GAM with joint spline term of age at sampling and year of isolation
##' 
##' @param hi_paired
##' @param conversion
##' @param panel, 1 = baseline
##'               2 = follow-up
##'               3 = changes between visits
##' 
##' @return fitted gam of titer on joint spline term of year of isolation and age at sampling
##

fitGam <- function(hi_paired, conversion, panel){
  
  # data
  if(panel %in% c(1,2)){
    
    if(panel==1){ # visit 1
      dat <- hi_paired[hi_paired$visit=="V1", ]
    }
    if(panel==2){ # visit 4
      dat <- hi_paired[hi_paired$visit=="V4", ]
    }
    
    model <- gam(log_hi_titer ~ s(age_at_sampling, year_of_isolation),
                 data = dat)
  }
  
  
  if(panel==3){ # delta
    
    dat <- conversion
    
    model <- gam(fold_of_change ~ s(age_baseline, year_of_isolation),
                 data = dat)
  }
  
  return(model)
  
}



##' Function to format data to calculate AUC
##' @param data hi_paired or its subset that contains post-birth strains
##' 
##' @return dataframe contains individual unnormalized AUC
##' 

aucData <- function(data){
  
  stat <- unique(data[data$visit=="V1" ,c("participant_id",
                                          "age_at_sampling")])
  stat <- data.frame(stat, stringsAsFactors = FALSE)
  
  tem <- t(sapply(stat$participant_id, function(i) aucCal(i, data)))
  colnames(tem) <- paste("AUC_", c("V1", "V4", "delta"), sep="")
  
  stat <- cbind(stat, tem)
  
  return(stat)
}


##' Function to calculate AUC
##' 
##' @param i, individual id
##' @param data, hi_paired or its subset that contains post-birth strains
##' 
##' @return dataframe with birth adjusted metrics
##' 

aucCal <- function(i, data){
  
  #id <- ids[i]
  
  df <- data[data$participant_id==i, ]
  
  df <- df[order(df$visit, df$year_of_isolation), ]
  
  df <- list(df[df$visit=="V1", ],
             df[df$visit=="V4", ])
  
  auc <- unlist(lapply(df, function(f){
    
    year <- f$year_of_isolation
    titer <- f$log_hi_titer
    
    len <- nrow(f)
    
    y1 <- titer[1:(len-1)]
    y2 <- titer[2:len]
    t1 <- year[1:(len-1)]
    t2 <- year[2:len]
    
    
    sum((y1 + y2)/2*(t2 - t1))
    
  }))
  
  return(c(auc, diff(auc)))
  
}


##' Function to fit univariable GAM with AUC
##' 
##' @param stat, dataframe contains individual AUC, derived from aucData
##' 
##' @return fitted GAM for baseline, follow-up and changes
##' 

aucGam <- function(stat){
  
  m <- list(gam(AUC_V1 ~ s(age_at_sampling), data=stat),
            gam(AUC_V4 ~ s(age_at_sampling), data=stat),
            gam(AUC_delta ~ s(age_at_sampling), data=stat))
  
  return(m)
}


##' Function to format data to calculate ATY
##' @param data conversion or its subset that contains post-birth strains
##' 
##' @return dataframe contains individual unnormalized AUC
##'

atyData <- function(hi_paired, data){
  
  stat <- unique(hi_paired[hi_paired$visit=="V1" ,c("participant_id",
                                                    "age_at_sampling")])
  
  stat <- data.frame(stat, stringsAsFactors = FALSE)
  
  
  tem <- t(sapply(stat$participant_id, function(i) atyCal(i, data)))
  colnames(tem) <- paste("ATY_", c("V1", "V4", "delta"), sep="")
  
  stat <- cbind(stat, tem)
  
  stat <- stat[!is.na(stat$ATY_delta), ]
  
  return(stat)
}


##' Function to calculate ATY
##' 
##' @param i, individual id
##' @param data, hi_paired or its subset that contains post-birth strains
##' 
##' @return dataframe with birth adjusted metrics
##' 

atyCal <- function(i, data){
  
  df <- data[data$participant_id==i, ]
  
  if(nrow(df) == 0){
    aty <- rep(NA, 3)
    return(aty)
  } else {
   
     df <- df[order(df$year_of_isolation), ]
    
    aty <- unlist(lapply(c("log_hi_baseline",
                           "log_hi_followup",
                           "fold_of_change"), function(x){
                             
                             year <- df$year_of_isolation
                             titer <- df[, x]
                             
                             len <- nrow(df)
                             
                             y1 <- titer[1:(len-1)]
                             y2 <- titer[2:len]
                             t1 <- year[1:(len-1)]
                             t2 <- year[2:len]
                             
                             auc <- sum((y1 + y2)/2*(t1 - t2), na.rm=TRUE)
                             aty <- sum((y1 + y2)*(t1 + t2)/4*(t1 - t2), na.rm=TRUE)/auc
                             if(is.infinite(aty)){aty <- NA}
                             
                             return(aty)
                             
                           }))
  }
  
  return(aty)
  
}


##' Function to fit univariable GAM with ATY
##' 
##' @param stat, dataframe contains individual ATY, derived from aucData
##' 
##' @return fitted GAM for baseline, follow-up and changes
##'

atyGam <- function(stat){
  
  m <- list(gam(ATY_V1 ~ s(age_at_sampling), data=stat),
            gam(ATY_V4 ~ s(age_at_sampling), data=stat),
            gam(ATY_delta ~ s(age_at_sampling), data=stat))
  
  return(m)
}


##' Function to prepare width data
##' 
##' @param data hi_paired or its subset that contains post-birth strains
##' @param z, cut-off 1 = 1:10,
##'                   3 = 1:40
##' @param normalize, TRUE or FALSE
##' 
##' @return dataframe of width for baseline, follow-up and difference between the two visits

widthData <- function(data, z, normalize){
  
  stat <- unique(data[data$visit=="V1" ,c("participant_id",
                                          "age_at_sampling")])
  stat <- data.frame(stat, stringsAsFactors = FALSE)
  
  tem <- t(sapply(stat$participant_id, function(i) widthCal(i, z, data)))
  colnames(tem) <- paste("width_", c("V1", "V4", "delta"), sep="")
  
  if(normalize){
    # normalized widths
    width.max <- 2014 - 1968
    width.min <- 0
    for(i in 1:2){
      tem[,i] <- (tem[,i]-width.min)/(width.max - width.min)
    }
    
    # normalized diff width
    width.min <- width.max*(-1)
    tem[,3] <- (tem[,3]-width.min)/(width.max - width.min)
  }
  
  stat <- cbind(stat, tem)
  
  return(stat)
  
}


##' Function to calculate the width of titer above threshold
##' 
##' @param i participant ID
##' @param z, cut-off 1 = 1:10,
##'                   3 = 1:40
##' @param data hi_paired or its subset that contains post-birth strains
##' 
##' @return width for visit 1, visit 4 and difference between the two visits
##' 

widthCal <- function(i, z, data){
  
  df <- data[data$participant_id==i, ]
  
  if(nrow(df) <= 2){
    return(rep(NA, 3))
  } else {
    df <- df[order(df$year_of_isolation), ]
    
    len <- length(unique(df$strain))
    
    tem <- unlist(lapply(c("V1", "V4"), function(x) {
      
      tem <- df[df$visit==x, ]
      
      y.1 <- tem$log_hi_titer[1:(len-1)]
      y.2 <- tem$log_hi_titer[2:len]
      t.1 <- tem$year_of_isolation[1:(len-1)]
      t.2 <- tem$year_of_isolation[2:len]
      
      t.1GreaterThanZ <- y.1>=z
      t.2GreaterThanZ <- y.2>=z
      
      width <- rep(0, len-1)
      
      # condition 1
      con <- t.1GreaterThanZ & t.2GreaterThanZ # both greater than z
      width[con] <- t.2[con] - t.1[con]
      
      # condition 2
      con <- !(t.1GreaterThanZ) & t.2GreaterThanZ # yi+1 greater than z & yi less than z
      width[con] <- t.2[con] - (((z-y.2)*t.1 + (y.1-z)*t.2)/(y.1-y.2))[con]
      
      # condition 3
      con <- t.1GreaterThanZ & !(t.2GreaterThanZ) # yi+1 less than z & yi greater than z
      width[con] <- (((z-y.2)*t.1 + (y.1-z)*t.2)/(y.1-y.2))[con] - t.1[con]
      
      return(sum(width))
      
    }))
    
    return(c(tem, diff(tem)))
    
  } 
  
}


##' Function to fit univariable GAM with width
##' 
##' @param stat, dataframe contains individual ATY, derived from aucData
##' 
##' @return fitted GAM for baseline, follow-up and changes
##'

widthGam <- function(stat){
  
  m <- list(gam(width_V1 ~ s(age_at_sampling), data=stat),
            gam(width_V4 ~ s(age_at_sampling), data=stat),
            gam(width_delta ~ s(age_at_sampling), data=stat))
  
  return(m)
}


##' Function to calculate the birth adjusted metrics
##' 
##' @param data subset of hi_paired which contains titers to post-birth strains
##' @param conversion dataframe
##' 
##' @return dataframe with birth adjusted metrics
##' 

birthAdjustedMetrics <- function(hi_paired, conversion){
  
  hi_paired$year_of_birth <- hi_paired$year_of_sampling - hi_paired$age_at_sampling
  data <- hi_paired[hi_paired$year_of_isolation > hi_paired$year_of_birth, ] # select post-birth strains
  conversion_used <- conversion[conversion$age_at_isolation > 0 & !is.na(conversion$age_at_isolation), ]
  
  # get summary metrics
  auc <- aucData(data)
  aty <- atyData(hi_paired, data = conversion_used)
  w10 <-  widthData(data, z = 1, normalize = FALSE)
  w40 <-  widthData(data, z = 3, normalize = FALSE)
  colnames(w10)[3:5] <- paste(colnames(w10)[3:5], "_10", sep = "")
  colnames(w40)[3:5] <- paste(colnames(w40)[3:5], "_40", sep = "")
  
  out <- left_join(auc, aty)
  out <- left_join(out, w10)
  out <- left_join(out, w40)
  
  tem <- t(sapply(out$participant_id, function(i){
    
    used <- data[data$participant_id == i, ]
    
    numSt <- length(unique(used$strain))
    
    uwATY <- mean(used$year_of_isolation)
    
    y1 <- max(used$year_of_birth, 1968)
    widthMax <- 2014 - y1
    
    return(c(numSt, uwATY, widthMax))
  }))
  colnames(tem) <- c("num_St_Exposed", "unweighted_ATY", "width_Max")
  
  out <- cbind(out, tem)
  
  for(i in 3:5){
    out[ ,i] <- out[ ,i]/out$num_St_Exposed
  }
  
  # for(i in 6:8){
  #   out[ ,i] <- out[ ,i] - out$unweighted_ATY
  # }
  
  for(i in 9:14){
    out[ ,i] <- out[ ,i] / out$width_Max
  }
  
  return(out)
  
}


##' Function to prepare data for seroconversion analysis
##' 
##' @param conversion dataframe
##' @param hi_paired subset of hi_paired which contains titers to post-birth strains
##' @param st dataframe contains strain informations
##' 
##' 
##' @return a list contains dataset for each of the four recent strains
##' 

getSeroconversionData <- function(conversion, hi_paired, st){
  
  id <- unique(conversion$participant_id)
  
  # four recent strains
  out <- lapply(1:4, function(s){
    
    out <- data.frame(matrix(nrow = length(id),
                             ncol = 2),
                      stringsAsFactors = FALSE)
    colnames(out) <- c("participant_id", "outcome_strain")
    
    out$participant_id <- id
    out$outcome_strain <- rep(st$strain[17 + s], length(id))
    
    
    # pre-titer for outcome strain
    st.used <- st$strain[17 + s]
    used.outcome <- conversion[conversion$strain == st.used, ]
    
    out <- left_join(out, 
                     used.outcome[ ,c("participant_id", 
                                      "age_baseline", 
                                      "seroconversion", 
                                      "log_hi_baseline")])
    
    # pre-titer for previous strain
    st.pre <- st$strain[16 + s]
    if(s == 2){
      st.pre <- st$strain[17]
    }
    used.pre <- conversion[conversion$strain == st.pre, ]
    colnames(used.pre)[colnames(used.pre) == 'log_hi_baseline'] <- "log_hi_baseline_pre" 
    out <- left_join(out, used.pre[ ,c("participant_id", "log_hi_baseline_pre")])
    
    
    # auc
    auc <- aty <- w10 <- w40 <- list()
    if(s == 2){
      st.auc <- st$strain[1:16]
    } else {
      st.auc <- st$strain[1:(15+s)] # this has been revised
    }
    
    yr.used <- st$year_of_isolation[17 + s]
    
    st.recent <- st.auc[(length(st.auc)-3):length(st.auc)]
    
    # auc from 1968 
    data <- hi_paired[hi_paired$strain %in% st.auc, ]
    auc[[1]] <- aucData(data)
    colnames(auc[[1]])[1] <- "participant_id"
    out <- left_join(out, auc[[1]][ ,c("participant_id", "AUC_V1")])
    
    # auc from birth
    data$year_of_birth <- data$year_of_sampling - data$age_at_sampling
    data <- data[data$year_of_birth <= data$year_of_isolation, ]
    auc[[2]] <- aucData(data)
    n <- sapply(auc[[2]]$participant_id, function(i){
      length(data$strain[data$participant_id == i & data$visit == "V1"])
    })
    auc[[2]]$AUC_V1_adj <-  auc[[2]]$AUC_V1/n
    colnames(auc[[2]])[c(1, 3, 6)] <- c("participant_id", "AUC_V1_birth", "AUC_V1_birth_adj")
    out <- left_join(out, auc[[2]][ ,c("participant_id", "AUC_V1_birth_adj")])
    
    
    # aty
    # aty from 1968 
    data <- conversion[conversion$strain %in% st.auc, ]
    aty[[1]] <- atyData(hi_paired, data)
    colnames(aty[[1]])[1] <- "participant_id"
    out <- left_join(out, aty[[1]][ ,c("participant_id", "ATY_V1")])
    
    # aty from birth
    data$year_of_birth <- 2010 - data$age_baseline
    data <- data[data$year_of_birth <= data$year_of_isolation, ]
    aty[[2]] <- atyData(hi_paired, data)
    colnames(aty[[2]])[c(1, 3)] <- c("participant_id", "ATY_V1_birth_adj")
    out <- left_join(out, aty[[2]][ ,c("participant_id", "ATY_V1_birth_adj")])
    
    
    # width > 10
    # width > 10 from 1968 
    data <- hi_paired[hi_paired$strain %in% st.auc, ]
    w10[[1]] <- widthData(data, 1, normalize = TRUE)
    colnames(w10[[1]])[1] <- "participant_id"
    out <- left_join(out, w10[[1]][ ,c("participant_id", "width_V1")])
    
    data$year_of_birth <- data$year_of_sampling - data$age_at_sampling
    data <- data[data$year_of_birth <= data$year_of_isolation, ]
    w10[[2]] <- widthData(data, 1, normalize = FALSE)
    width.max <- sapply(w10[[2]]$participant_id, function(i){
      used <- data[data$participant_id == i, ]
      y <- max(used$year_of_birth, 1968)
      yr.used - y
    })
    w10[[2]]$width_V1_birth_adj <- w10[[2]]$width_V1 / width.max
    colnames(w10[[2]])[1] <- "participant_id"
    out <- left_join(out, w10[[2]][ ,c("participant_id", "width_V1_birth_adj")])
    
    
    # width > 40
    # width > 40 from 1968 
    data <- hi_paired[hi_paired$strain %in% st.auc, ]
    w40[[1]] <- widthData(data, 3, normalize = TRUE)
    colnames(w40[[1]])[c(1, 3)] <- c("participant_id", "width40_V1")
    out <- left_join(out, w40[[1]][ ,c("participant_id", "width40_V1")])
    
    data$year_of_birth <- data$year_of_sampling - data$age_at_sampling
    data <- data[data$year_of_birth <= data$year_of_isolation, ]
    w40[[2]] <- widthData(data, 3, normalize = FALSE)
    width.max <- sapply(w40[[2]]$participant_id, function(i){
      used <- data[data$participant_id == i, ]
      y <- max(used$year_of_birth, 1968)
      # 2014 - y
      yr.used - y
    })
    w40[[2]]$width40_V1_birth_adj <- w40[[2]]$width_V1 / width.max
    colnames(w40[[2]])[1] <- "participant_id"
    out <- left_join(out, w40[[2]][ ,c("participant_id", "width40_V1_birth_adj")])
    

    colnames(out) <- c("id",                      "outcome_strain",           "age_at_sampling",
                       "seroconversion",          "titer_outcome_strain_v1",  "titer_previous_strain_v1",
                       "AUC_from_1968",           "AUC_from_birth_adj",       
                       "ATY_from_1968",           "ATY_from_birth_adj",
                       "W10_from_1968",           "W10_from_birth_adj", 
                       "W40_from_1968",           "W40_from_birth_adj")
    
    out <- out[rowSums(is.na(out)) == 0, ]
    
    return(out)
    
  })
  
  return(out)
}


##' Function to prepare data for distribution of changes for 4 recent strains
##' 
##' @param conversion dataframe
##' 
##' 
##' @return a list contains dataset of
##'         1) Distribution of number of strains with increase and 
##'            the changes among each category
##'         2) Distribution of changes for each of the four strains
##'         3) Chi-squared test
##' 

getDataForTiterChanges <- function(conversion){
  
  dat <- conversion[conversion$year_of_isolation >= 2009, ]
  ids <- unique(dat$participant_id)
  
  out <- spread(dat[, c("participant_id", 
                        "strain",
                        "fold_of_change")],
                strain, fold_of_change)
  for(i in 2:5){
    out[out[,i] < 0, i] <- -1 # decrease
    out[out[,i] > 1, i] <- 2 # seroconversion
  }
  
  out <- out[ ,c("participant_id",
                 "A/Perth/2009",
                 "A/Victoria/2009",
                 "A/Texas/2012",
                 "A/HongKong/2014")]
  
  out$NumStrainsIncrease <- rowSums(out[,2:5] > 0)
  
  
  byNumber <- matrix(nrow = 5, ncol = 5)
  byStrain <- matrix(nrow = 4, ncol = 5)
  
  pValues <- rep(NA, 2)
  
  # by number
  xs <- matrix(nrow = 5, ncol = 4)
  for(i in 0:4){
    
    used <- out[out$NumStrainsIncrease==i, ]
    
    N <- nrow(used)
    
    used <- c(used$`A/Perth/2009`, used$`A/Victoria/2009`, 
              used$`A/Texas/2012`, used$`A/HongKong/2014`)
    x <- sapply(-1:2, function(t) sum(used==t))
    n <- sum(x)
    p <- x/n
    
    byNumber[i+1, ] <- c(N, p)
    xs[i+1, ] <- x
    
  }
  
  est <- chisq.test(xs)
  if(est$p.value < 0.01){
    pValues[1] <- '<0.01'
  } else {pValues[1] <- sprintf('%.02f', est$p.value)}
  
  # by strain
  xs <- matrix(nrow = 4, ncol = 4)
  for(i in 1:4){
    
    used <- out[ ,i+1]
    x <- sapply(-1:2, function(t) sum(used==t))
    n <- sum(x)
    p <- x/n
    ci <- sapply(x, function(t) binom.test(t, n)$conf.int)
    
    N <- nrow(out)
    
    byStrain[i, ] <- c(N, p)
    
    xs[i, ] <- x
    
  }
  
  est <- chisq.test(xs)
  if(est$p.value < 0.01){
    pValues[2] <- '<0.01'
  } else {pValues[2] <- sprintf('%.02f', est$p.value)}
  
  
  return(list(byNumber = byNumber,
              byStrain = byStrain,
              est = est))
  
}


##' 
##' Function to fit and summary logistic regressions of seroconversion on fuller exposure history
##' results used in Table 1, Table S3-4
##' 
##' @param out list of dataset for each of the four recent strains, output from function 'getSeroconversionData'
##' @param range strains included 
##'              1 = all tested strains
##'              2 = post-birth strains
##' 
##' @return a table contains results from logistic regressions
##'

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
                         titer_previous_strain_v1)
    m <- glm(base, data = out[[i]], family = binomial())
    
    coeff <- round(exp(coef(m)[2:4]), dgt)
    ci <- round(exp(confint(m)[2:4, ]), dgt)
    
    p <- summary(m)$coefficients[2:4, 4]
    star <- rep("", 3)
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
    titer_previous_strain_v1"
    
    for(k in 1:4){ # summary metric
      
      metric <- c("AUC", "ATY", "W10", "W40")[k]
      rg <- c('1968', 'adj')[range]
      
      used <- colnames(out[[i]])[grepl(metric, colnames(out[[i]]))]
      used <- used[grepl(rg, used)]
      
      fml <- as.formula(paste(base, used, sep = " + "))
      m <- glm(fml, data = out[[i]], family = binomial())
      
      coeff <- round(exp(coef(m)[2:5]), dgt)
      ci <- round(exp(confint(m)[2:5, ]), dgt)
      
      p <- summary(m)$coefficients[2:5, 4]
      star <- rep("", 4)
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


##' 
##' Function to summary decreased/unchanged/increased titer by individual
##' 
##' @param conversion dataset contains changes of titers
##' 
##' @return matrix contains age, number of decrease/increase/seroconversion
##'

getAgeTiterChanges <- function(conversion){
  
  ids <- unique(conversion$participant_id)
  
  out <- t(sapply(ids, function(i){
    
    used <- conversion$fold_of_change[conversion$participant_id == i]
    
    age <- unique(conversion$age_baseline[conversion$participant_id == i])[1]
    
    d <- sum(used < 0, na.rm = TRUE)
    
    u <- sum(used == 0, na.rm = TRUE)
    
    a <- sum(used > 0, na.rm = TRUE)
    
    b <- sum(used >= 2, na.rm = TRUE)
    
    return(c(age, d, u, a, b))
    
  }))
  
  out <- data.frame(out, stringsAsFactors = FALSE)
  colnames(out) <- c("age", "numDecrease", "numNoChange", "numIncrease", "numFourFold")
  
  return(out)
  
}


##' 
##' Function to boostrap the age distribution of the studied population
##' 
##' @param stat dataset contains individual age; derived from getAgeTiterChanges()
##' @param n number of boostraps
##' 
##' @return matrix contains ci for each age groups
##'

boostrap <- function(stat, n){
  
  age <- stat$age
  p <- table(floor(age/10))/length(age)
  p <- p/sum(p)
  
  tem <- lapply(1:n, function(i){
    a.new <- sample(0:8, 777, replace = TRUE, prob = p)
    n.used <- sapply(0:8, function(j) sum(a.new == j))
    return(n.used/777)
  })
  tem <- do.call("rbind", tem)
  
  t(sapply(1:length(p), function(x) quantile(tem[,x], c(.025, .975))))  
}



##' 
##' Function to calculate proportion of changes in GMT
##' 
##' @param conversion
##' @param st
##' @param population
##' 
##' @return matrix contains ci for each age groups
##'

getPropChange <- function(conversion, st, population){
  
  res <- lapply(1:2, function(i){
    
    dat <- conversion
    
    if(i == 1){ # post-birth strain only
      dat <- conversion[conversion$age_at_isolation >= 0, ]
    }
    
    sc <- 2.5
    
    titers <- lapply(st$strain, function(s){
      
      t1 <- dat[dat$strain == s & !is.na(dat[ ,'log_hi_baseline']), 
                'log_hi_baseline']
      t2 <- dat[dat$strain == s & !is.na(dat[ ,'log_hi_followup']), 
                'log_hi_followup']
      t3 <- dat[dat$strain == s & !is.na(dat[ ,'fold_of_change']), 
                'fold_of_change']
      
      return(cbind(t1, t2, t3))
      
    })
    
    out <- lapply(1:2, function(j){
      
      sapply(1:21, function(s){
        
        t <- titers[[s]][ ,j]
        t <- 2^t*5
        t <- log2(t/sc)
        
        est <- t.test(t)
        a <- c(est$estimate, est$conf.int)
        return(a)
        
      })
      
    })
    
    out[[3]] <- sapply(1:21, function(s){
      
      t1 <- titers[[s]][ ,1]
      t1 <- 2^t1*5
      t1 <- log2(t1/sc)
      
      t2 <- titers[[s]][ ,2]
      t2 <- 2^t2*5
      t2 <- log2(t2/sc)
      
      t <- (t2 - t1)/t1
      
      est <- t.test(t)
      a <- c(est$estimate, est$conf.int)
      return(a)
      
    })
    
    return(out)
    
  })
  
  return(res)
  
}


##' 
##' Function to derive 95% CI of strain-specific of hi titer on delta titer
##' 
##' @param dat convserion
##' 
##' @return dataframe of point estimates and 95% CI
##'

lmInteractCI <- function(dat, st){
  
  n <- nrow(st)
  
  strain.factor <- data.frame(matrix(dat$strain,
                                     nrow = nrow(dat), ncol = n,
                                     byrow = FALSE),
                              stringsAsFactors = FALSE)
  for(i in 1:n){
    
    st.used <- c(st$strain[i],
                 st$strain[setdiff(1:n, i)])
    
    strain.factor[,i] <- factor(strain.factor[,i], 
                                levels = st.used)
    
  }
  
  colnames(strain.factor) <- paste("strain", 1:n, sep="")
  dat <- cbind(dat, strain.factor)
  
  est <- matrix(nrow = n, ncol = 3)
  rownames(est) <- st$strain
  
  for(i in 1:n){
    
    formular <- paste("fold_of_change ~ log_hi_baseline + age_baseline + strain",
                      i, 
                      "+ log_hi_baseline*strain", i,
                      sep = "")
    
    m <- do.call(lm, list(as.formula(formular), 
                          data=dat))
    
    # extract point and confidence interval
    p <- m$coefficients[2]
    ci <- confint(m)[2, ]
    
    est[i, ] <- c(p, ci)
    
  }
  
  
  return(est)
  
}


##' 
##' Function to prepare data used for univariable logistic regression 
##' on seroconversion to four recent strains
##' 
##' @param conversion
##' @param hi_paired
##' @param st
##' 
##' @return list of dataset used for univariable logistic regression 
##'

getDataForUniTiter <- function(conversion, hi_paired, st){
  
  id <- unique(conversion$participant_id)
  
  st$strain.new <- gsub('/', '_', st$strain)
  
  # four recent strains
  out <- lapply(1:4, function(s){
    
    out <- data.frame(matrix(nrow = length(id),
                             ncol = 2),
                      stringsAsFactors = FALSE)
    colnames(out) <- c("participant_id", "outcome_strain")
    
    out$participant_id <- id
    out$outcome_strain <- rep(st$strain[17 + s], length(id))
    
    
    # pre-titer for outcome strain
    st.used <- st$strain[17 + s]
    used.outcome <- conversion[conversion$strain == st.used, ]
    
    out <- left_join(out, used.outcome[ ,c("participant_id", 
                                           "age_baseline", 
                                           "seroconversion", 
                                           "log_hi_baseline")])
    
    # pre-titer for previous strain
    for(k in 1:(16 + s)){
      
      st.pre <- st$strain[k]
      used.pre <- conversion[conversion$strain == st.pre, ]
      col.name <- paste('titer', st$strain.new[k], sep = '_') 
      colnames(used.pre)[colnames(used.pre) == 'log_hi_baseline'] <- col.name
      out <- left_join(out, used.pre[ ,c("participant_id", col.name)])
      
      
    }
    colnames(out)[1:5] <- c("id",                      "outcome_strain",           "age_at_sampling",
                            "seroconversion",          "titer_outcome_strain_v1")
    
    
    return(out)
    
  })
  
  
  return(out)
  
}


##' 
##' Function to get coefficients from univariable logistic regression 
##' on seroconversion to four recent strains
##' results used for comparing results with pre-existing titer to strain i
##' 
##' @param data, uniData returned from getDataForUniTiter()
##' @param strain ith outcome strain
##'               1 = A/Perth/2009
##'               2 = A/Victoria/2009
##'               3 = A/Texas/2012
##'               4 = A/HongKong/2014
##' @param includeStraini whether or not incluede pre-existing titer to strain i
##' 
##' @return list of coefficients from univariable logistic regression 
##'

getStrainCoeffForSeroconv <- function(data, strain, includeStraini){
  
  dat.used <- data[[strain]]
  col.used <- colnames(dat.used)[6:ncol(dat.used)]
  nTotal <- length(col.used)
  dgt <- 2
  
  if(includeStraini){
    base <- 'seroconversion ~ age_at_sampling + titer_outcome_strain_v1 + '
    n <- 4
  } else {
    base <- 'seroconversion ~ '
    n <- 2
  }
  
  
  est <- sapply(1:nTotal, function(k){
    fml <- as.formula(paste(base,
                            col.used[k],
                            sep = ''))
    m <- glm(fml, data = dat.used, family = binomial())
    coeff <- round(exp(coef(m))[n], dgt)
    ci <- round(exp(confint(m)[n, ]), dgt)
    return(c(coeff, ci))
  })
  
  colnames(est) <- col.used
  
  return(est)
  
}



##' 
##' Function to get coefficients from univariable logistic regression 
##' on seroconversion to four recent strains
##' results used for heatmap plot
##' 
##' @param conversion
##' @param st
##' 
##' @return list of coefficients from univariable logistic regression 
##'

getUniCoeff <- function(conversion, st){
  
  strains <- st$strain
  
  coeff <- lapply(strains, function(s.out){
    
    t <- conversion$fold_of_change[conversion$strain == s.out]
    
    c.out <- rep(0, length(t))
    c.out[t >= 2] <- 1
    
    sapply(strains, function(s.pred){
      
      pred <- conversion$log_hi_baseline[conversion$strain == s.pred]
      
      m <- summary(glm(c.out ~ pred, family = binomial()))
      
      round(m$coefficients[2, c(1, 4)], 3)
      
    })
    
  })
  
  coeff <- do.call('rbind', coeff)
  
}



##' 
##' Function to get coefficients from multivariable logistic regression 
##' on seroconversion to four recent strains
##' results used for heatmap plot
##' 
##' @param hi_paired
##' @param conversion
##' @param st
##' 
##' @return list of coefficients from multivariable logistic regression 
##'

getMulCoeff <- function(hi_paired, conversion, st){
  
  strains <- st$strain
  
  hi_paired_wide <- spread(hi_paired[,c("participant_id",
                                        "visit",
                                        "age_at_sampling",
                                        "strain",
                                        "log_hi_titer")],
                           strain, log_hi_titer
  )
  hi_paired_wide <- hi_paired_wide[hi_paired_wide$visit == 'V1', ]
  hi_paired_wide <- hi_paired_wide[ ,-2]
  colnames(hi_paired_wide)[1] <- 'participant_id'
  hi_paired_wide <- cbind(hi_paired_wide[ ,1:2], hi_paired_wide[ ,st$strain])
  
  colnames(hi_paired_wide)[3:23] <- gsub('/', '_', colnames(hi_paired_wide)[3:23])
  
  coeff.m <- lapply(strains, function(s.out){
    
    t <- conversion[conversion$strain == s.out, 
                    c('participant_id', "seroconversion")]
    
    tem <- left_join(hi_paired_wide, t)
    
    predictors <- paste(colnames(tem)[2:23], collapse = ' + ')
    fml <- as.formula(paste('seroconversion', predictors, sep = ' ~ '))
    
    m <- summary(glm(fml, family = binomial(), data = tem))
    
    round(t(m$coefficients[3:23, c(1, 4)]), 3)
    
  })
  
  coeff.m <- do.call('rbind', coeff.m)
  
}


##' 
##' Function to generate tables comparing AICs/BICs
##' 
##' @param out data used to assess seroconversion, returned from getSeroconversionData()
##' 
##' @return table contains AICs/BICs
##'

AICCompareTable <- function(out){
  
  tab <- matrix(nrow = ncol(out[[1]]) - 2, ncol = 8)
  rownames(tab) <- c('NULL',
                     'Age',
                     "Age + Titer i",
                     "Age + Titer i + Titer i-1",
                     paste("Age + Titer i + Titer i-1",
                           " + ",
                           colnames(out[[1]])[7:ncol(out[[1]])],
                           sep = ""))
  colnames(tab) <- c(sapply(1:4, function(i) paste(c('AIC', 'BIC'), 
                                                   out[[i]]$outcome_strain[1], 
                                                   sep = "_")))
  
  for(i in 1:4){ # strain
    
    # null
    fml <- as.formula(seroconversion ~ 1)
    m <- glm(fml, data = out[[i]], family = binomial())
    tab[1, (i-1)*2+1] <- AIC(m)
    tab[1, i*2] <- BIC(m)
    
    # age only
    fml <- as.formula(seroconversion ~ age_at_sampling)
    m <- glm(fml, data = out[[i]], family = binomial())
    tab[2, (i-1)*2+1] <- AIC(m)
    tab[2, i*2] <- BIC(m)
    
    # age + titer
    fml <- update(fml, . ~ . + titer_outcome_strain_v1)
    m <- glm(fml, data = out[[i]], family = binomial())
    tab[3, (i-1)*2+1] <- AIC(m)
    tab[3, i*2] <- BIC(m)
    
    # age + titer + pre titer
    fml <- update(fml, . ~ . + titer_previous_strain_v1)
    m <- glm(fml, data = out[[i]], family = binomial())
    tab[4, (i-1)*2+1] <- AIC(m)
    tab[4, i*2] <- BIC(m)
    
    # age + titer + pre titer + others
    base <- "seroconversion ~ age_at_sampling + titer_outcome_strain_v1 + titer_previous_strain_v1"
    for(k in 1:(nrow(tab)-4)){
      
      used <- colnames(out[[1]])[k + 6]
      fml <- as.formula(paste(base, used, sep = " + "))
      m <- glm(fml, data = out[[i]], family = binomial())
      tab[k+4, (i-1)*2+1] <- AIC(m)
      tab[k+4, i*2] <- BIC(m)
      
    }
    
  }
  
  tab <- data.frame(tab, stringsAsFactors = FALSE)
  
  rnames <- rownames(tab)
  metric <- c('AUC', "ATY", 'W10', 'W40')
  rg <- c('1968', 'adj')
  tab$metric <- tab$range <-  'NA'
  for(k in 1:4){
    tab$metric[grepl(metric[k], rnames)] <- metric[k]
  }
  for(k in 1:2){
    tab$range[grepl(rg[k], rnames)] <- rg[k]
  }
  
  tab <- tab[c(1:6,
               11:12,
               9:10,
               7:8), ]

  
  return(tab)
  
}


##' 
##' Function to generate tables comparing AICs/BICs with spline on age
##' 
##' @param out data used to assess seroconversion, returned from getSeroconversionData()
##' 
##' @return table contains AICs/BICs from models with spline on age
##'

AICCompareTableSpline <- function(out){
  
  tab <- matrix(nrow = ncol(out[[1]]) - 2, ncol = 8)
  rownames(tab) <- c('NULL',
                     'Age',
                     "Age + Titer i",
                     "Age + Titer i + Titer i-1",
                     paste("Age + Titer i + Titer i-1",
                           " + ",
                           colnames(out[[1]])[7:ncol(out[[1]])],
                           sep = ""))
  colnames(tab) <- c(sapply(1:4, function(i) paste(c('AIC', 'BIC'), 
                                                   out[[i]]$outcome_strain[1], 
                                                   sep = "_")))
  
  for(i in 1:4){
    
    # null
    fml <- as.formula(seroconversion ~ 1)
    m <- glm(fml, data = out[[i]], family = binomial())
    tab[1, (i-1)*2+1] <- AIC(m)
    tab[1, i*2] <- BIC(m)
    
    # age only
    fml <- as.formula(seroconversion ~ s(age_at_sampling))
    m <- gam(fml, data = out[[i]], family = binomial())
    tab[2, (i-1)*2+1] <- AIC(m)
    tab[2, i*2] <- BIC(m)
    
    # age + titer
    fml <- update(fml, . ~ . + titer_outcome_strain_v1)
    m <- gam(fml, data = out[[i]], family = binomial())
    tab[3, (i-1)*2+1] <- AIC(m)
    tab[3, i*2] <- BIC(m)
    
    # age + titer + pre titer
    fml <- update(fml, . ~ . + titer_previous_strain_v1)
    m <- gam(fml, data = out[[i]], family = binomial())
    tab[4, (i-1)*2+1] <- AIC(m)
    tab[4, i*2] <- BIC(m)
    
    # age + titer + pre titer + others
    base <- "seroconversion ~ s(age_at_sampling) + titer_outcome_strain_v1 + titer_previous_strain_v1"
    for(k in 1:(nrow(tab)-4)){
      
      used <- colnames(out[[1]])[k + 6]
      fml <- as.formula(paste(base, used, sep = " + "))
      m <- gam(fml, data = out[[i]], family = binomial())
      tab[k+4, (i-1)*2+1] <- AIC(m)
      tab[k+4, i*2] <- BIC(m)
      
    }
    
  }
  
  tab <- data.frame(tab, stringsAsFactors = FALSE)
  
  
  rnames <- rownames(tab)
  metric <- c('AUC', "ATY", 'W10', 'W40')
  rg <- c('1968', 'adj')
  tab$metric <- tab$range <-  'NA'
  for(k in 1:4){
    tab$metric[grepl(metric[k], rnames)] <- metric[k]
  }
  for(k in 1:2){
    tab$range[grepl(rg[k], rnames)] <- rg[k]
  }
  
  tab <- tab[c(1:6,
               11:12,
               9:10,
               7:8), ]
  
  return(tab)
  
}


##' 
##' Function to estimate mediation effect of previous exposure on seroconversion
##' 
##' @param out data used to assess seroconversion, returned from getSeroconversionData()
##' @interact whether or not to include interaction between previous exposure 
##' and pre-existing titers to the outcome strain
##' 
##' @return 
##'

getMedEff <- function(out, interact){
  
  eff.auc <- eff.w40 <- matrix(nrow = 9, ncol = 4)
  
  if(interact){
    
    for(i in 1:4){
      
      # auc
      model.0 <- glm(seroconversion ~ age_at_sampling + AUC_from_birth_adj,
                     data = out[[i]], family = binomial("probit"))
      
      model.M <- lm(titer_outcome_strain_v1 ~ age_at_sampling + AUC_from_birth_adj,
                    data = out[[i]])
      model.Y <- glm(seroconversion ~ age_at_sampling + titer_outcome_strain_v1 + AUC_from_birth_adj + 
                       titer_outcome_strain_v1 * AUC_from_birth_adj,
                     data = out[[i]], family = binomial("probit"))
      
      results <- mediate(model.M, model.Y, 
                         treat='AUC_from_birth_adj', 
                         mediator='titer_outcome_strain_v1',
                         boot=TRUE, sims = 1000,
                         control.value = 0, treat.value = 1)
      
      eff.auc[ ,i] <- c(results$d.avg, results$d.avg.ci,
                        results$z.avg, results$z.avg.ci,
                        results$tau.coef, results$tau.ci)
      
      rm(results)
      
      # width 40
      model.0 <- glm(seroconversion ~ age_at_sampling + W40_from_birth_adj,
                     data = out[[i]], family = binomial("probit"))
      
      model.M <- lm(titer_outcome_strain_v1 ~ age_at_sampling + W40_from_birth_adj,
                    data = out[[i]])
      model.Y <- glm(seroconversion ~ age_at_sampling + titer_outcome_strain_v1 + W40_from_birth_adj +
                       titer_outcome_strain_v1 * W40_from_birth_adj,
                     data = out[[i]], family = binomial("probit"))
      
      results <- mediate(model.M, model.Y, 
                         treat='W40_from_birth_adj', 
                         mediator='titer_outcome_strain_v1',
                         boot=TRUE, sims = 1000,
                         control.value = 0, treat.value = 1)
      
      eff.w40[ ,i] <- c(results$d.avg, results$d.avg.ci,
                        results$z.avg, results$z.avg.ci,
                        results$tau.coef, results$tau.ci)
      
      rm(results)
      
      
    }
    
  }
  
  if(!interact){
    
    for(i in 1:4){
      
      # auc
      model.0 <- glm(seroconversion ~ age_at_sampling + AUC_from_birth_adj,
                     data = out[[i]], family = binomial("probit"))
      
      model.M <- lm(titer_outcome_strain_v1 ~ age_at_sampling + AUC_from_birth_adj,
                    data = out[[i]])
      model.Y <- glm(seroconversion ~ age_at_sampling + titer_outcome_strain_v1 + AUC_from_birth_adj,
                     data = out[[i]], family = binomial("probit"))
      
      results <- mediate(model.M, model.Y, 
                         treat='AUC_from_birth_adj', 
                         mediator='titer_outcome_strain_v1',
                         boot=TRUE, sims = 1000,
                         control.value = 0, treat.value = 1)
      
      eff.auc[ ,i] <- c(results$d.avg, results$d.avg.ci,
                        results$z.avg, results$z.avg.ci,
                        results$tau.coef, results$tau.ci)
      
      rm(results)
      
      
      # width 40
      model.0 <- glm(seroconversion ~ age_at_sampling + W40_from_birth_adj,
                     data = out[[i]], family = binomial("probit"))
      
      model.M <- lm(titer_outcome_strain_v1 ~ age_at_sampling + W40_from_birth_adj,
                    data = out[[i]])
      model.Y <- glm(seroconversion ~ age_at_sampling + titer_outcome_strain_v1 + W40_from_birth_adj,
                     data = out[[i]], family = binomial("probit"))
      
      results <- mediate(model.M, model.Y, 
                         treat='W40_from_birth_adj', 
                         mediator='titer_outcome_strain_v1',
                         boot=TRUE, sims = 1000,
                         control.value = 0, treat.value = 1)
      
      eff.w40[ ,i] <- c(results$d.avg, results$d.avg.ci,
                        results$z.avg, results$z.avg.ci,
                        results$tau.coef, results$tau.ci)
      
      rm(results)
      
      
    }
    
  }
  
  return(list(eff.auc, eff.w40))
  
}


##' 
##' Function to get predictions from seroconversion models
##' 
##' @param out data used to assess seroconversion, returned from getSeroconversionData()
##' @spline whether or not to include spline term on age
##' 
##' @return list of predictions
##

getPred <- function(out, spline){
  
  rtList <- list()
  
  for(i in 1:4){ # strain
    
    used <- out[[i]]
    
    rt <- data.frame(age.group = NA,
                     seroconversion = used$seroconversion,
                     pred.auc = NA,
                     pred.width40 = NA,
                     stringsAsFactors = FALSE)
    rt$age.group <- floor(used$age_at_sampling/10)
    rt$age.group[rt$age.group == 8] <- 7
    rt$age.group[rt$age.group == 0] <- 1
    
    if(!spline){
      
      # prediction from auc model
      fml <- as.formula("seroconversion ~ age_at_sampling + titer_outcome_strain_v1 + titer_previous_strain_v1 + AUC_from_birth_adj")
      m <- glm(fml, family = binomial(), data = used)
      rt$pred.auc <- plogis(predict(m))
      
      # prediction from w40 model
      fml <- as.formula("seroconversion ~ age_at_sampling + titer_outcome_strain_v1 + titer_previous_strain_v1 + W40_from_birth_adj")
      m <- glm(fml, family = binomial(), data = used)
      rt$pred.width40 <- plogis(predict(m))
      
    }
    
    if(spline){
      
      # prediction from auc model
      fml <- as.formula("seroconversion ~ s(age_at_sampling) + titer_outcome_strain_v1 + titer_previous_strain_v1 + AUC_from_birth_adj")
      m <- gam(fml, family = binomial(), data = used)
      rt$pred.auc <- plogis(predict(m))
      
      # prediction from w40 model
      fml <- as.formula("seroconversion ~ s(age_at_sampling) + titer_outcome_strain_v1 + titer_previous_strain_v1 + W40_from_birth_adj")
      m <- gam(fml, family = binomial(), data = used)
      rt$pred.width40 <- plogis(predict(m))
      
    }
    
    
    # predictions for age groups
    rtList[[i]] <- sapply(1:7, function(i){
      
      # observations
      x <- sum(rt$seroconversion[rt$age.group == i])
      n <- sum(rt$age.group == i)
      est <- binom.test(x, n)
      obs <- c(est$estimate, est$conf.int)
      
      # prediction by auc
      pred.auc <- quantile(rt$pred.auc[rt$age.group == i],
                           c(.5, .25, .75))
      pred.width40 <- quantile(rt$pred.width40[rt$age.group == i],
                               c(.5, .25, .75))
      
      return(c(obs, pred.auc, pred.width40))
      
      
    })
    
    
  }
  
  return(rtList)
  
}


