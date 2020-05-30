#
# Descriptions in main text
#

library(dplyr)

# ------------------------------------ LOAD DATA ---------------------------------------%

load('/Volumes/GoogleDrive/My Drive/R projects/Fluscape/r codes/spatial_pattern_serology/Fluscape_serology/data/data.RData')
source('/Volumes/GoogleDrive/My Drive/R projects/Fluscape/r codes/spatial_pattern_serology/Fluscape_serology/rCodes/figureFuncs.R')
source('/Volumes/GoogleDrive/My Drive/R projects/Fluscape/r codes/spatial_pattern_serology/Fluscape_serology/rCodes/statFuncs.R')

# ------------------------------------ PREPARE DATA ---------------------------------------%

stat <- birthAdjustedMetrics(hi_paired, conversion)

models_auc <- aucGam(stat)
models_aty <- atyGam(stat)

tem = stat[ ,c('participant_id', 'age_at_sampling', 
               'width_V1_40', 'width_V4_40', 'width_delta_40')]
colnames(tem) = c('participant_id', 'age_at_sampling', 
                  'width_V1', 'width_V4', 'width_delta')
models_w40 <- widthGam(tem)

# -------------------------------------- PREDICT ------------------------------------------%


getBootstrapOfAUCRatio <- function(stat, n){
  
  new = data.frame(age_at_sampling = c(20, 50),
                   stringsAsFactors = FALSE)
  nrows <- nrow(stat)
  
  out <- lapply(1:n, function(x){
    
    rNo <- sample(1:nrows, nrows, replace = TRUE)
    stat.new <- stat[rNo, ]
    
    m_auc <- aucGam(stat.new)
    
    pred.auc <- c(predict(m_auc[[1]], 
                          newdata = new),
                  predict(m_auc[[2]], 
                          newdata = new))
    
    return(pred.auc)
    
    
  })
  
  out <- do.call('rbind', out)
  
  colnames(out) <- c('auc_v1_20', 'auc_v1_50',
                     'auc_v4_20', 'auc_v4_50')
  
  return(out)
  
}

out <- getBootstrapOfAUCRatio(stat, n = 1000)

quantile(out[,1]/out[,2],
         c(.025, .975))

quantile(out[,3]/out[,4],
         c(.025, .975))


getBootstrapOfWidthRatio <- function(stat, n){
  
  new = data.frame(age_at_sampling = seq(0, 80, 0.2),
                   stringsAsFactors = FALSE)
  nrows <- nrow(stat)
  
  out <- lapply(1:n, function(x){
    
    rNo <- sample(1:nrows, nrows, replace = TRUE)
    stat.new <- stat[rNo, ]
    
    # w40
    tem = stat.new[ ,c('participant_id', 'age_at_sampling', 
                       'width_V1_40', 'width_V4_40', 'width_delta_40')]
    colnames(tem) = c('participant_id', 'age_at_sampling', 
                      'width_V1', 'width_V4', 'width_delta')
    m_w40 <- widthGam(tem)
    
    ratio_w40 <- sapply(1:2, function(i){
      
      pred <- predict(m_w40[[i]], 
                      newdata = new)
      peak <- max(pred)
      w40_50yr <- pred[new$age_at_sampling == 50]
      
      w40_50yr/peak
      
    })
    
    # w10
    tem = stat.new[ ,c('participant_id', 'age_at_sampling', 
                       'width_V1_10', 'width_V4_10', 'width_delta_10')]
    colnames(tem) = c('participant_id', 'age_at_sampling', 
                      'width_V1', 'width_V4', 'width_delta')
    m_w10 <- widthGam(tem)
    
    ratio_w10 <- sapply(1:2, function(i){
      
      pred <- predict(m_w10[[i]], 
                      newdata = new)
      peak <- max(pred)
      w10_50yr <- pred[new$age_at_sampling == 50]
      
      w10_50yr/peak
      
    })
    
    return(c(ratio_w40, ratio_w10))
    
    
  })
  
  out <- do.call('rbind', out)
  
  colnames(out) <- c('w40_v1', 'w40_v4',
                     'w10_v1', 'w10_v4')
  
  return(out)
  
}

ratio <- getBootstrapOfWidthRatio(stat, n = 1000)

r <- sapply(1:4, function(x) quantile(out[,x],
                                      c(.025, .975)))
round(r, 2)


# aty by seroconversion to recent strains

stat$seroconversion = sapply(stat$participant_id, function(i){
  
  x = conversion$seroconversion[conversion$participant_id == i&
                                  conversion$strain %in% c('A/Perth/2009',
                                                           'A/Victoria/2009',
                                                           'A/Texas/2012',
                                                           'A/HongKong/2014')]
  sum(x) > 0
  
})

ggplot(stat, aes(x = seroconversion,
                 y = ATY_V1)) +
  geom_boxplot()




