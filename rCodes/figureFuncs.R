#
# 2019-10-16
# Functions for plots
#

#------------------------------- Figure 1 & S1 ----------------------------------------%

##' Function to prepare data used for figure 1 
##' 
##' @param hi_paired
##' 
##' @return data that removed missing values in sex and year of sampling
##'

fig1Data <- function(hi_paired){
  
  del <- rowSums(is.na(hi_paired[ ,c("sex", "year_of_sampling")]))
  hi_paired <- hi_paired[del==0, ]
  hi_paired$year_of_birth <- hi_paired$year_of_sampling - hi_paired$age_at_sampling
  
  hi_paired <- hi_paired[order(hi_paired$visit,
                               hi_paired$year_of_birth, 
                               hi_paired$year_of_isolation), ]
  
  return(hi_paired)
  
}


##' Function to randomly choose id from age groups
##' 
##' @param paired, dataframe from "fig1Data" function
##' 
##' @return a vector contains randomly selected participant id
##'

selID <- function(paired){
  
  # age group 
  x <- list(c(0, 10), 
            c(11, 20), 
            c(21, 40),
            c(41, 60),
            c(61, 80),
            c(81, 100))
  
  ids <- rep(NA, 6)
  
  for(i in 1:6){
    
    x.used <- x[[i]][1]:x[[i]][2]
    
    con <- paired$age_at_sampling %in% x.used & paired$visit == 'V1'
    
    id.used <- unique(paired$participant_id[con])
    
    ids[i] <- sample(id.used, 1)
    
  }
  
  return(ids)
  
}


##' Function to plot individual antibody profile
##' 
##' @param i, panel i
##' @param ids, a vector of randomly selected participant id
##' @param paired, dataframe from "fig1Data" function
##' 
##' @return panel of Fig. 1
##'

individualPlot <- function(i, ids, paired){
  
  id <- ids[i]
  age <- unique(paired$age_at_sampling[paired$participant_id==id & paired$visit == 'V1'])
  
  # identify titer
  con.1 <- paired$participant_id==id & 
    paired$visit==c("V1")
  con.2 <- paired$participant_id==id & 
    paired$visit==c("V4")
  
  # individual data for plot
  data <- NULL
  
  data$strain <- unique(paired$strain)
  data$titer.1 <- sapply(data$strain, function(x) {
    a <- paired$log_hi_titer[con.1 & paired$strain==x]
    if(length(a)==0){a <- NA}
    return(a)
  })
  data$titer.2 <- sapply(data$strain, function(x) {
    a <- paired$log_hi_titer[con.2 & paired$strain==x]
    if(length(a)==0){a <- NA}
    return(a)
  })
  data$year <- sapply(data$strain, function(x){
    
    unique(paired$year_of_isolation[paired$strain==x])
    
  })
  
  data$delta.titer <- data$titer.2 - data$titer.1
  
  data <- data.frame(data)
  data <- data[!is.na(data$titer.1), ]
  
  # splined data
  
  sm.1 <- spline(data$year, data$titer.1, method = "natural")
  sm.2 <- spline(data$year, data$titer.2, method = "natural")
  
  
  data.spline <- data.frame(year = sm.1$x,
                            titer.1 = sm.1$y,
                            titer.2 = sm.2$y)
  
  # replace negative with 0
  data.spline$titer.1[data.spline$titer.1<0] <- 0
  data.spline$titer.2[data.spline$titer.2<0] <- 0
  
  data.spline$increase <- data.spline$titer.2 - data.spline$titer.1
  data.spline$increase[data.spline$increase<=0] <- 0
  data.spline$decrease <- data.spline$titer.2 - data.spline$titer.1
  data.spline$decrease[data.spline$decrease>=0] <- 0
  
  # year of birth
  hi_paired$year_of_birth <- hi_paired$year_of_sampling - hi_paired$age_at_sampling
  birth <- min(unique(hi_paired$year_of_birth[hi_paired$participant_id == id]))
  
  x.min <- 1965
  x.max <- 2015
  x.diff.1 <- 5
  x.diff.2 <- 1
  
  y.min <- 0
  y.max <- 8
  
  cols <- c("#3f93aa",
            "#fa8072") # blue and red
  
  cols.change <- c("#590059",
                   "#8BAD23") # purple and green
  
  par(mar = c(4,6,4,1))
  plot(NULL, 
       xlim = c(x.min, x.max), 
       ylim = c(y.min, y.max),
       las = 1,
       xlab = c("Year of isolation"),
       ylab = c("HAI titer"),
       main = paste(age, " years",
                    sep = ""),
       axes = FALSE,
       cex.lab = 1.5,
       cex.main = 1.5)
  
  # time period of visit
  rect(2009.923,
       y.min,
       2011.058,
       y.max,
       col = adjustcolor(cols[1], 
                         alpha.f = 0.2),
       border = FALSE) # visit 1
  
  rect(2014.458,
       y.min,
       2015.416,
       y.max,
       col = adjustcolor(cols[2], 
                         alpha.f = 0.2),
       border = FALSE) # visit 4
  
  
  # titer reference line
  segments(x.min, log2(10/5),
           x.max, log2(10/5),
           lty = 2)
  segments(x.min, log2(40/5),
           x.max, log2(40/5),
           lty = 3)
  
  # brown area
  rowMins <- sapply(1:nrow(data.spline), function(i){
    max(min(data.spline[i, 2:3]), 0)
  })
  
  polygon(x = c(data.spline$year, rev(data.spline$year)),
          y = c(rep(0, length(data.spline$year)), 
                rev(rowMins)),
          border = NA,
          col = adjustcolor("#867b7b",
                            alpha.f = 0.5))
  #col = "#c0b1b1")
  
  for(j in 1:2){
    
    x <- data.spline$year
    y <- data.spline$titer.1 + data.spline[ ,c("increase", "decrease")[j]]
    
    polygon(x = c(x, rev(x)),
            y = c(data.spline$titer.1, rev(y)),
            border = NA,
            col = adjustcolor(cols.change[j],
                              alpha.f = 0.7))
    
    points(data$year, data[ ,c("titer.1", "titer.2")[j]],
           pch = c(1,2)[j],
           col = cols[j],
           cex = 1.5,
           lwd = 1.5)
    
    lines(x, data.spline[ ,c("titer.1", "titer.2")[j]],
          lty = c(1,3)[j],
          col = cols[j],
          lwd = 2)
    
  }
  
  # time period before birth
  segments(birth, y.min, birth, y.max, lty = 4)
  
  
  # add axis
  axis(1, seq(x.min, x.max, x.diff.1), 
       pos = y.min - 0.5, cex.axis = 1.5)
  axis(1, seq(x.min, x.max, x.diff.2), 
       pos = y.min - 0.5, labels = NA,
       tcl = -0.2, col.ticks = 1, 
       col=NA, cex.axis = 1.5)
  
  axis(2, seq(y.min, y.max, 1), 
       pos = x.min - 0.5, 
       las=1,
       labels = c(0, 2^seq(y.min+1, y.max, 1)*5), 
       cex.axis = 1.5)
  
  mtext(# letters[i], 
    LETTERS[i],
    side = 2,
    line = 3.5, 
    at = y.max + (y.max - y.min)*0.15,
    outer = FALSE, 
    las = 1,
    cex = 1.5,
    font = 2)
  
  
}

individualMetricValuePlot <- function(i, ids, stat){
  
  # text metric values
  used <- stat[stat$participant_id == ids[i], ]
  
  x.lim <- c(0, 1)
  y.lim <- c(0, 1)
  
  
  par(mar = c(4,10,2,7))
  plot(NULL, 
       xlim = x.lim, 
       ylim = y.lim,
       las = 1,
       xlab = '',
       ylab = '',
       axes = FALSE,
       cex.lab = 1.5,
       cex.main = 1.5)
  
  y.pos <- seq(y.lim[2], y.lim[1], length.out = 4)
  x.pos <- seq(x.lim[1], x.lim[2], length.out = 5)
  
  text(x.pos[1],
       y.pos[2:4],
       c('Baseline',
         'Follow-up',
         'Delta'),
       cex = 1.5,
       adj = c(0.5, 0.5),
       xpd = TRUE)
  
  text(x.pos[2:5],
       y.pos[1],
       c('nAUC',
         'nW40',
         'nW10',
         'nATY'),
       cex = 1.5,
       adj = c(0.5, 0.5),
       xpd = TRUE)
  

  text(x.pos[2],
       y.pos[2:4],
       sprintf('%.01f',
               c(used$AUC_V1,
                 used$AUC_V4,
                 used$AUC_delta)),
       cex = 1.5,
       adj = c(0.5, 0.5),
       xpd = TRUE)
  
  text(x.pos[3],
       y.pos[2:4],
       sprintf('%.01f',
               c(used$width_V1_40,
                 used$width_V4_40,
                 used$width_delta_40)),
       cex = 1.5,
       adj = c(0.5, 0.5),
       xpd = TRUE)
  
  text(x.pos[4],
       y.pos[2:4],
       sprintf('%.01f',
               c(used$width_V1_10,
                 used$width_V4_10,
                 used$width_delta_10)),
       cex = 1.5,
       adj = c(0.5, 0.5),
       xpd = TRUE)
  
  text(x.pos[5],
       y.pos[2:4],
       sprintf('%.01f',
               c(used$ATY_V1,
                 used$ATY_V4,
                 used$ATY_delta)),
       cex = 1.5,
       adj = c(0.5, 0.5),
       xpd = TRUE)
  
  space <- 0.12
  rect(x.lim[1] - diff(x.lim)*space,
       y.lim[1] - diff(y.lim)*space,
       x.lim[2] + diff(x.lim)*space,
       y.lim[2] + diff(y.lim)*space,
       xpd = TRUE)
  
  segments(x.lim[1] - diff(x.lim)*space,
           mean(y.pos[1:2]),
           x.lim[2] + diff(x.lim)*space,
           mean(y.pos[1:2]),
           xpd = TRUE)
  
  segments(mean(x.pos[1:2]),
           y.lim[1] - diff(y.lim)*space,
           mean(x.pos[1:2]),
           y.lim[2] + diff(y.lim)*space,
           xpd = TRUE)
  
}



#------------------------------- Figure 2 ----------------------------------------%

##' Function to plot heatmap of titers at baseline and follow-up
##' 
##' @param data, hi_paired
##' @param st, st
##' @param panel, 1 = baseline
##'               2 = follow-up
##' 
##' @return panel A and B of Fig. 2
##'

titerAgeHeatmap <- function(data, st, panel){
  
  strain <- st[order(st$year_of_isolation), ]
  strain <- strain$strain
  
  df <- data[data$visit == c('V1', 'V4')[panel], ]
  
  # convert to a 21*777 matrix
  used <- acast(df,
                df$strain ~ df$participant_id,
                value.var = "log_hi_titer")
  
  # order participants by age
  ids <- data.frame(Participant_ID = colnames(used),
                    stringsAsFactors = FALSE)
  ids$Age <- sapply(1:nrow(ids), function(x){
    unique(data$age_at_sampling[data$participant_id == ids$Participant_ID[x] &
                                  data$visit == 'V1'])
  }) # use baseline age for all panels
  
  ids <- ids[order(ids$Age), ]
  
  ids$Age.group <- floor(ids$Age/10) + 1
  ids$Age.group[ids$Age.group>=8] <- 8 # >=70 yrs
  
  n.age <- table(ids$Age.group)
  
  # order matrix by age
  used <- used[ ,ids$Participant_ID]
  
  x.min <- 0
  x.max <- ncol(used)
  y.min <- 0
  y.max <- length(strain)
  
  cols <- colorRampPalette(c("#9FBDEF", "#efd19f", "#e5635a"))(9)
  
  margin <- c(4, 12.5, 1, 0)
  
  par(mar = margin,
      new = FALSE)
  plot(NULL, xlim=c(x.min, x.max),
       ylim=c(y.min, y.max),
       xlab = "", ylab = "",
       axes = FALSE)
  
  for(j in 1:y.max){ # string
    
    tt <- used[rownames(used)==strain[j], ]
    
    sapply(1:x.max, function(i){
      rect(xleft = i-1,
           # ybottom = y.min+j-0.9,
           ybottom = y.min+j-1,
           xright = i,
           # ytop = y.min+j-0.1,
           ytop = y.min+j,
           col = cols[tt[i]+2],
           border = NA)
    })
    
  }
  
  rect(x.min, y.min, x.max, y.max)
  
  strain[strain == 'X31'] <- 'X31 (1970)'
  axis(2,
       #at = seq(y.max-0.5, y.min, -1),
       at = seq(y.min+0.5, y.max, 1),
       labels = strain,
       pos= x.min,
       las = 1,
       cex.axis = 1.5)
  
  axis(1, cumsum(n.age),
       pos = y.min,
       labels = NA, tick = TRUE,
       tcl = 0.4, col.ticks = 1, col = NA,
       cex.axis = 1.5)
  
  axis(1, cumsum(n.age),
       pos = y.max,
       labels = NA, tick = TRUE,
       tcl = -0.4, col.ticks = 1, col = NA,
       cex.axis = 1.5)
  
  txt <- c("0~9",
           "10~19",
           "20~29",
           "30~39",
           "40~49",
           "50~59",
           "60~69",
           "70+")
  text(c(0,cumsum(n.age))[1:8]+n.age/2,
       rep(y.min - 0.5, 8),
       txt,
       xpd=TRUE,
       srt = 90,
       adj =c(1, 0.5),
       cex = 1.5)
  
  txt <- "A"
  if(panel == 2){txt <- "B"}
  
  # txt <- "a"
  # if(panel == 4){txt <- "b"}
  
  mtext(txt,
        side = 2,
        las = 1,
        at = y.max,
        line = 10.5,
        outer = FALSE,
        cex = 2,
        font = 2)
  
  
}


##' Function to plot heatmap of change in titers
##' 
##' @param data, conversion
##' @param st, st
##' @param panel, 1 = baseline
##'               2 = follow-up
##' 
##' @return panel C of Fig. 2
##'

converAgeHeatmap <- function(data, st){
  
  strain <- st[order(st$year_of_isolation), ]
  strain <- strain$strain
  
  df <- data[data$strain != "", ]
  
  used <- acast(df,
                df$strain ~ df$participant_id,
                value.var = "fold_of_change")
  used[used < -4] <- -4
  used[used > 4] <- 4
  
  # order participants by month
  ids <- data.frame(Participant_ID = colnames(used),
                    stringsAsFactors = FALSE)
  ids$Age <- sapply(1:nrow(ids), function(x) unique(df$age_baseline[df$participant_id==ids$Participant_ID[x]]))
  
  ids <- ids[order(ids$Age), ]
  
  ids$Age.group <- floor(ids$Age/10) + 1
  ids$Age.group[ids$Age.group>=8] <- 8 # >=70 yrs
  
  n.age <- table(ids$Age.group)
  
  n.month <- sapply(1:8, function(x) sum(ids$Age.group==x))
  
  # order matrix by age
  used <- used[ ,ids$Participant_ID]
  
  x.min <- 0
  x.max <- ncol(used)
  y.min <- 0
  y.max <- length(strain)
  
  cols <- c(colorRampPalette(c("#8BAD23",
                               "#e6f8da"))(4),
            # "#f2f2f2",
            "#E5E5E5",
            colorRampPalette(c("#cdb2cd",
                               "#590059"))(4))
  
  
  margin <- c(4, 12.5, 1, 0)
  
  par(mar = margin,
      new = FALSE)
  plot(NULL, xlim=c(x.min, x.max),
       ylim=c(y.min, y.max),
       xlab = "", ylab = "",
       axes = FALSE)
  
  for(j in 1:y.max){ # strainring
    
    tt <- used[rownames(used)==strain[j], ]
    
    sapply(1:x.max, function(i){
      rect(xleft = i-1,
           ybottom = y.min+j - 1,
           # ybottom = y.min+j-0.9,
           xright = i,
           ytop = y.min+j,
           # ytop = y.min+j-0.1,
           col = cols[tt[i]+5],
           border = NA)
    })
    
  }
  
  rect(x.min, y.min, x.max, y.max)
  strain[strain == 'X31'] <- 'X31 (1970)'
  
  axis(2,
       #at = seq(y.max-0.5, y.min, -1),
       at = seq(y.min+0.5, y.max, 1),
       labels = strain,
       pos= x.min,
       las = 1,
       cex.axis = 1.5)
  
  
  axis(1, cumsum(n.month),
       pos = y.min,
       labels = NA, tick = TRUE,
       tcl = 0.4, col.ticks = 1, col = NA,
       cex.axis = 1.5)
  
  axis(1, cumsum(n.month),
       pos = y.max,
       labels = NA, tick = TRUE,
       tcl = -0.4, col.ticks = 1, col = NA,
       cex.axis = 1.5)
  
  txt <- c("0~9",
           "10~19",
           "20~29",
           "30~39",
           "40~49",
           "50~59",
           "60~69",
           "70+")
  text(c(0,cumsum(n.month))[1:8]+n.month/2,
       rep(y.min - 0.5, 8),
       txt,
       xpd=TRUE,
       srt = 90,
       adj =c(1, 0.5),
       cex = 1.5)
  
  mtext(# 'c',
    "C",
    side = 2,
    las = 1,
    at = y.max,
    line = 10.5,
    outer = FALSE,
    cex = 2,
    font = 2)
  
}


##' Function to plot distribution of titers (across strains and participants)
##' 
##' @param conversion
##' @param panel, 1 = baseline
##'               2 = follow-up
##'               3 = changes between visits
##' 
##' @return panel E, E and F of Fig. 2
##'

titerDistribution <- function(conversion, panel){
  
  cols <- colorRampPalette(c("#9FBDEF", "#efd19f", "#e5635a"))(9)
  
  margin <- c(4, 6, 1, 0)
  x.min <- c(0, 0, -4)[panel]
  x.max <- c(8, 8,  4)[panel]
  y.min <- 0
  y.max <- 0.4
  
  x <- seq(x.min, x.max, 1)
  y <- sapply(x, function(i) sum(conversion[ ,c("log_hi_baseline", 
                                                "log_hi_followup", 
                                                "fold_of_change")[panel]]==i, 
                                 na.rm = TRUE))
  if(panel == 3){
    y[1] <- sum(conversion$fold_of_change <= x.min, na.rm = TRUE)
    y[length(y)] <- sum(conversion$fold_of_change >= x.max, na.rm = TRUE)
    
    cols <- c(colorRampPalette(c("#8BAD23",
                                 "#e6f8da"))(4),
              #"#f2f2f2",
              "#E5E5E5",
              colorRampPalette(c("#cdb2cd",
                                 "#590059"))(4))
  }
  p <- y/sum(y)
  
  par(mar = margin)
  plot(NULL,
       xlim = c(x.min, x.max+1),
       ylim = c(y.min, y.max),
       xlab = "",
       ylab = "",
       axes = FALSE,
       cex.lab = 1.5)
  
  rect(x+0.2,
       rep(0, length(y)),
       x + 0.8,
       p,
       col = cols,
       border = NA)
  
  axis(1, c(x, x.max+1), pos = y.min, labels = NA, cex.axis = 1.5)
  
  if(panel == 3){
    txt <- c(expression(""<="1/16"),
             expression('1/8'),
             expression('1/4'),
             expression('1/2'),
             expression("1"),
             expression("2"),
             expression("4"),
             expression("8"),
             expression("">="16"))
    text(x + 0.5, rep(y.min - (y.max - y.min)*0.05),
         txt, xpd = TRUE, cex = 1.5)
  } else {
    txt <- 2^x*5
    txt[1] <- 0
    text(x + 0.5, rep(y.min - (y.max - y.min)*0.05),
         txt, xpd = TRUE, cex = 1.5)
  }
  mtext(c("HAI titer",
          "HAI titer",
          "Fold of HAI titer changes")[panel], 1, line = 2)
  
  axis(2, seq(y.min, y.max, 0.1), 
       pos = x.min, las = 1,
       labels = seq(y.min, y.max, 0.1)*100, cex.axis = 1.5)
  mtext("Proportion, %", 2, line = 3)
  
  mtext(# c("d", "e", "f")[panel],
    c("D", "E", "F")[panel],
    side = 2,
    las = 1,
    at = y.max,
    line = 4.5,
    outer = FALSE,
    cex = 2,
    font = 2)
  
}


##' Function to generate countour plot
##' 
##' @param conversion
##' @param panel, 1 = baseline
##'               2 = follow-up
##'               3 = changes between visits
##' 
##' @return panel G, H and I of Fig. 2
##

contourPlot <- function(hi_paired, conversion, panel){
  
  model <- fitGam(hi_paired, conversion, panel)
  n.grid <- 200
  
  x.min <- 0
  x.max <- 90
  y.min <- 1965
  y.max <- 2015
  z.min <- c(0,0,-1)[panel]
  z.max <- c(6,6,3)[panel]
  len <- length(z.min:z.max)
  
  margin <- c(4, 6, 1, 0)
  
  if(panel %in% c(1,2)){
    len.cols <- length(0:8)
    cols <- colorRampPalette(c("#9FBDEF", "#efd19f", "#e5635a"))(len.cols)
  } else {
    len.cols <- length(-4:4)
    cols <- c(colorRampPalette(c("#8BAD23",
                                 "#e6f8da"))(4),
              "#E5E5E5",
              colorRampPalette(c("#cdb2cd",
                                 "#590059"))(4))
  }
  if(panel %in% c(1,2)){
    cols <- cols[1:6]
  } else cols <- cols[4:8]
  
  age.start <- seq(x.min, x.max, length = n.grid)
  CirculatingYear <- seq(y.min, y.max, length = n.grid)
  
  newd <- data.frame(matrix(0, n.grid * n.grid, 0))
  newd[,1] <- rep(age.start, n.grid)
  newd[,2] <- c(sapply(CirculatingYear, function(x) rep(x, n.grid)))
  
  if(panel%in%c(1,2)){
    names(newd) <- c("age_at_sampling", "year_of_isolation")
  } else {
    names(newd) <- c("age_baseline", "year_of_isolation")
  }
  
  
  fv <- predict.gam(model, newdata = newd, 
                    se.fit = TRUE, type = "response")
  z <- fv$fit
  if(panel %in% c(1,2)){
    z[z<0] <- 0 # make negative log titer to 0 -> undetectable
  }
  
  z <- matrix(z, n.grid, n.grid)
  
  par(mar = margin)
  plot(NULL, 
       xlim=c(x.min, x.max), 
       ylim=c(y.min, y.max), xlab="", ylab="", axes=F)
  .filled.contour(age.start, CirculatingYear, z,
                  # xlim = c(x.min, x.max),
                  # ylim = c(y.min, y.max),
                  # zlim = c(z.min, z.max),
                  #nlevels = len,
                  levels = z.min:z.max,
                  col = cols)
  lvls <- len*2
  if(panel==3){lvls <- len*5}
  contour(age.start, CirculatingYear, z,
          xlim = c(x.min, x.max),
          ylim = c(y.min, y.max),
          zlim = c(z.min, z.max),
          nlevels = lvls,
          lty = 1,
          add = TRUE,
          col = grey(0.6))
  
  axis(1, at = seq(x.min, x.max, 10), pos = y.min, cex.axis = 1.5)
  mtext("Age at sampling, years", 1, line = 2)
  axis(2, at = seq(y.min, y.max, 10), pos = x.min, las=1, cex.axis = 1.5)
  mtext("Year of isolation", 2, line = 3)
  
  segments(x.min, y.max, x.max, y.max)
  segments(x.max, y.min, x.max, y.max)
  
  segments(0, 2010, 45, 1965, lty=3, lwd=2)
  
  mtext(# c("g", "h", "i")[panel],
    c("G", "H", "I")[panel],
    side = 2,
    las = 1,
    at = y.max,
    line = 4.5,
    outer = FALSE,
    cex = 2,
    font = 2)
  
}


##' Function to plot color scales for titers
##' 
##' @return color scale for baseline and follow-up visit, Fig. 2
##

scalePlot2 <- function(){
  
  cols <- colorRampPalette(c("#9FBDEF", "#efd19f", "#e5635a"))(9)
  
  x.min <- 0
  x.max <- 1
  y.min <- 0
  y.max <- 2
  
  y.diff <- (y.max-y.min)/length(cols)
  y.pos <- seq(y.min, y.max, y.diff)
  
  x.pos <- c(x.min+0.2,
             x.max-0.7)
  
  margin <- c(3, 0, 1, 0)
  
  par(mar = margin,
      new = FALSE)
  plot(NULL, xlim=c(x.min, x.max),
       #ylim=c(y.min-y.max*1.5, y.max),
       ylim=c(y.min-y.max, y.max*2),
       xlab = "", ylab = "",
       axes = FALSE)
  
  rect(rep(x.pos[1], 9),
       y.pos[1:9],
       rep(x.pos[2], 9),
       y.pos[2:10],
       col = cols,
       border = NA)
  
  
  rect(x.pos[1],
       y.min,
       x.pos[2],
       y.max)
  
  axis(2, y.pos[2:9],
       pos =x.pos[1],
       labels = NA, tick = TRUE,
       tcl = 0.4, col.ticks = 1, col = NA,
       cex.axis = 1.5)
  
  axis(4, y.pos[2:9],
       pos = x.pos[2],
       labels = NA, tick = TRUE,
       tcl = 0.4, col.ticks = 1, col = NA,
       cex.axis = 1.5)
  
  txt <- 2^c(0:8)*5
  txt[1] <- 0
  text(rep(x.pos[2]+0.13, 9),
       c(y.pos+y.diff/2)[1:9],
       txt,
       xpd = TRUE,
       cex = 1.5,
       adj = c(0, 0.5))
  
  text(x.pos[2]+0.5,
       (y.min+y.max)/2,
       "HAI titer",
       srt = 90,
       font = 2,
       cex = 1.5,
       xpd = TRUE)
  
  #mtext("HI titer", 2, line = -2, font = 2)
  
  
}


##' Function to plot color scales for changes in titers between visits
##' 
##' @return color scale for changes, Fig. 2
##

converScalePlot <- function(){
  
  cols <- c(colorRampPalette(c("#8BAD23",
                               "#e6f8da"))(4),
            "#E5E5E5",
            colorRampPalette(c("#cdb2cd",
                               "#590059"))(4))
  
  
  x.min <- 0
  x.max <- 1
  y.min <- 0
  y.max <- 2
  
  y.diff <- (y.max-y.min)/length(cols)
  y.pos <- seq(y.min, y.max, y.diff)
  
  x.pos <- c(x.min+0.2,
             x.max-0.7)
  
  margin <- c(3, 0, 1, 0)
  
  par(mar = margin,
      new = FALSE)
  plot(NULL, xlim=c(x.min, x.max),
       ylim=c(y.min, y.max),
       xlab = "", ylab = "",
       axes = FALSE)
  
  rect(rep(x.pos[1], 9),
       y.pos[1:9],
       rep(x.pos[2], 9),
       y.pos[2:10],
       col = cols,
       border = NA)
  
  
  rect(x.pos[1],
       y.min,
       x.pos[2],
       y.max)
  
  axis(2, y.pos[2:9],
       pos =x.pos[1],
       labels = NA, tick = TRUE,
       tcl = 0.4, col.ticks = 1, col = NA,
       cex.axis = 1.5)
  
  axis(4, y.pos[2:9],
       pos = x.pos[2],
       labels = NA, tick = TRUE,
       tcl = 0.4, col.ticks = 1, col = NA,
       cex.axis = 1.5)
  
  txt <- c(expression(""<="1/16"),
           expression('1/8'),
           expression('1/4'),
           expression('1/2'),
           expression("1"),
           expression("2"),
           expression("4"),
           expression("8"),
           expression("">="16"))
  
  text(rep(x.pos[2]+0.12, 9),
       c(y.pos+y.diff/2)[1:9],
       txt,
       xpd = TRUE,
       cex = 1.5,
       adj = c(0, 0.5))
  
  text(x.pos[2] + 0.5,
       (y.min+y.max)/2,
       "Fold of HAI titer changes",
       srt = 90,
       font = 2,
       cex = 1.5,
       xpd = TRUE)
  
  #mtext("HI titer", 2, line = -2, font = 2)
  
  
}


#------------------------------- Figure 3 ----------------------------------------%

##'
##' Function to plot spline of titer on age stratifiy by strains
##'
##' @param data conversion
##' @param st st, dataset constains information on strains
##' @param panel 1 = baseline
##'              2 = follow-up
##'              3 = changes
##' @param age.type "sampling" or "isolation"             
##' 
##' @return Fig. 3
##'

splineTiterOnAgePlotCombined <- function(data, st, panel, age.type){
  
  # plot
  x.min.plot <- x.min <- 0
  x.max <- x.max.plot <- 80
  if(age.type=="isolation"){
    x.min.plot <- -20
    x.max.plot <- 80
  }
  y.min <- 0
  y.max <- 8
  
  # splines
  spl <- NULL
  for(i in 1:nrow(st)){
    
    used <- data[data$strain==st$strain[i], ]
    del <- rowSums(is.na(used))
    used <- used[del==0, ] 
    
    if(age.type == "sampling"){
      
      fml <- as.formula(paste(c("log_hi_baseline",
                                "log_hi_followup",
                                "fold_of_change")[panel],
                              " ~ s(age_baseline)"))
      m <- gam(fml, data = used)
      
      new.data <- data.frame(age_baseline = seq(x.min, x.max, 0.1))
    }
    
    if(age.type == "isolation"){
      
      x.min <- st$year_of_isolation[i] - 2010
      x.max <- min(max(used$age_at_isolation), x.max.plot)
      
      fml <- as.formula(paste(c("log_hi_baseline",
                                "log_hi_followup",
                                "fold_of_change")[panel],
                              " ~ s(age_at_isolation)"))
      m <- gam(fml, data = used)
      new.data <- data.frame(age_at_isolation = seq(x.min, x.max, 0.1))
    }
    
    spl[[i]] <- cbind(new.data[,1],
                      predict(m, newdata = new.data))
    
  }
  
  margin <- c(6,8,2,0)
  cols <- c('#800000', '#e6194b', '#fabebe', '#9a6324', '#f58231',
            '#ffd8b1', '#ffe119', '#808080', '#bfef45', '#3cb44b',
            '#aaffc3', '#469990', '#42d4f4', '#000000', '#000075',
            '#4363d8', '#911eb4', '#e6beff', '#f032e6', '#a9a9a9', '#996300')
  cols <- adjustcolor(cols, alpha.f = 0.8)
  
  if(age.type == "sampling"){
    mainTitle <- c("Baseline visit",
                   "Follow-up visit",
                   "Differences between two visits")[panel]
  } else {
    mainTitle <- ''
  }
  
  par(new = FALSE,
      mar = margin)
  
  plot(NULL, xlim=c(x.min.plot, x.max),
       ylim=c(y.min, y.max),
       xlab = paste("Age at ", age.type, ", years", sep = ""), 
       ylab = "HAI titer",
       axes = FALSE,
       main = mainTitle,
       cex.main = 2,
       cex.lab = 1.5)
  
  a <- sapply(1:nrow(st), function(i) 
    lines(spl[[i]][,1], spl[[i]][,2], lwd = 1.5, col = cols[i]))
  
  axis(2, seq(y.min, y.max, 2), las = 1,
       labels = 2^(seq(y.min, y.max, 2))*5, 
       cex.axis = 1.5, 
       pos = x.min.plot - (x.max.plot-x.min.plot)*0.05)
  axis(2, seq(y.min, y.max, 1), las = 1, 
       cex.axis = 1.5, 
       pos = x.min.plot - (x.max.plot-x.min.plot)*0.05,
       tcl = -0.2, labels = NA, xpd = TRUE,
       col = NA, col.ticks = "black")
  
  
  axis(1, at = seq(x.min.plot, x.max, 20), 
       cex.axis = 1.5, pos = y.min - (y.max-y.min)*0.05)
  axis(1, at = seq(x.min.plot, x.max.plot, 10), 
       cex.axis = 1.5, pos = y.min - (y.max-y.min)*0.05,
       tcl = -0.3, labels = NA, xpd = TRUE,
       col = NA, col.ticks = "black")
  axis(1, at = seq(x.min.plot, x.max.plot, 1), 
       cex.axis = 1.5, pos = y.min - (y.max-y.min)*0.05,
       tcl = -0.1, labels = NA, xpd = TRUE,
       col = NA, col.ticks = "black")
  
  if(age.type == "sampling"){
    txt <- LETTERS[panel]
    # txt <- letters[panel]
  } else {
    txt <- LETTERS[panel+3]
    # txt <- letters[panel+3]
  }
  
  mtext(txt, 
        side = 2,
        line = 4.5, 
        at = y.max + 0.5,
        outer = FALSE, 
        las = 1,
        cex = 2,
        font = 2)
  
}


##'
##' Function to plot legend of strains

##' @param st st, dataset constains information on strains
##' 
##' @return legend of Fig. 3
##'

splineTiterOnAgePlotLegend <- function(st){
  
  # plot
  x.min <- 0
  x.max <- 1
  y.min <- 0
  y.max <- 1
  
  margin <- c(4,0,4,0)
  
  cols <- c('#800000', '#e6194b', '#fabebe', '#9a6324', '#f58231',
            '#ffd8b1', '#ffe119', '#808080', '#bfef45', '#3cb44b',
            '#aaffc3', '#469990', '#42d4f4', '#000000', '#000075',
            '#4363d8', '#911eb4', '#e6beff', '#f032e6', '#a9a9a9', '#996300')
  cols <- adjustcolor(cols, alpha.f = 0.8)
  
  par(new = FALSE,
      mar = margin)
  
  plot(NULL, xlim=c(x.min, x.max),
       ylim=c(y.min, y.max),
       xlab = "", ylab = "",
       axes = FALSE)
  
  legend(x = "center",
         legend = st$strain,
         col = cols,
         lty = 1,
         bty = "n",
         xpd = TRUE,
         cex = 1.5,
         lwd = 1.5)
  
}


#------------------------------- Figure 4 ----------------------------------------%

##'
##' Example plot of summary statistics of antibody profile
##' 
##' @param hi_paired
##' @param conversion
##' @param panel, 1 = auc, 2 = width, 3 = aty
##' 
##' @return panel A, E and I of Fig. 4
##' 

conceptualPlotStatistics <- function(hi_paired, conversion, panel){
  
  id <- "L50H16P01"
  
  used <- hi_paired[hi_paired$participant_id==id & hi_paired$visit == "V4", 
                    c("year_of_isolation", "log_hi_titer")]
  used <- used[order(used$year_of_isolation), ]
  
  x.min <- 1965
  x.max <- 2015
  x.diff.1 <- 5
  x.diff.2 <- 1
  
  y.min <- 0
  y.max <- 8
  
  cols <- c("#7b896a", # green
            "#e6a953", # dark green
            "#0099cc") # blue
  
  cols <- c("#7b896a", # green
            "#c3bc1a", # dark green
            #"#545e49")
            "#5bd0c8") # blue
  
  
  # par(mar = c(4,6,2,0))
  par(mar = c(6,6,5,2),
      mgp = c(4, 1, 4))
  plot(NULL, 
       xlim = c(x.min, x.max), 
       ylim = c(y.min, y.max),
       las = 1,
       xlab = c("Year of isolation"),
       # ylab = c("HAI titer"),
       ylab = c("log HAI titer"),
       axes = FALSE,
       cex.lab = 1.5)
  
  # obserbed data
  x <- used$year_of_isolation
  y <- used$log_hi_titer
  
  lines(x, y, lwd = 2, col = "#999999", type = "b")
  
  if(panel == 1){ # auc
    
    polygon(x = c(x, rev(x)),
            y = c(rep(0, nrow(used)), 
                  rev(y)),
            border = NA,
            col = adjustcolor(cols[panel],
                              alpha.f = 0.5))
    len <- length(x)
    
    sapply(1:(len-1), function(i){
      
      polygon(x = c(x[i], x[i], x[i+1], x[i+1]),
              y = c(y[i],0,  0, y[i + 1]),
              col = NULL,
              border = "#999999",
              lty = 3)
      
    })
    
    # legend for auc
    points(1968, 7.5, pch = 15, col = cols[panel])
    
    segments(1967, 6.7, 1969, 6.7, col = "#999999", lty = 1)
    points(1968, 6.7, pch = 1, col = "#999999")
    
    # text(1972, c(7.5, 6.7),
    #      c('nAUC', 'Observation'),
    #      adj = c(0, 0.5))
    
    text(1972, c(7.5, 6.7),
         c('nAUC = 5.6', 'Observation'),
         adj = c(0, 0.5))

    
  }
  
  if(panel == 2){
    
    cutoff <- 3
    y.1 <- y
    y.1[y<cutoff] <- cutoff
    
    solid <- unlist(sapply(1:length(y), function(i){
      
      if(i == 1){
        rt <- all(y[i:(i+1)] >= 3)
      }
      
      if(i > 1 & i < length(y)){
        rt <- all(y[(i-1):i] >= 3) | all(y[i:(i+1)] >= 3)
      }
      
      if(i == length(y)){
        rt <- all(y[(i-1):i] >= 3)
      }
      
      rt
      
    }))
    
    points(x[solid],
           y[solid],
           col = '#b87624',
           pch = 19)
    
    # polygon(x = c(x, rev(x)),
    #         y = c(rep(cutoff, nrow(used)), 
    #               rev(y.1)),
    #         border = NA,
    #         col = adjustcolor(cols[panel],
    #                           alpha.f = 0.5))
    
    lines(x, rep(cutoff, nrow(used)), col = cols[panel], lty = 2, lwd = 2)
    
    rect(x[1],
         y.min - 0.15,
         x[21],
         y.min - 0.4,
         col = '#e9d5bd',
         border = NA)
    
    rect(x[1],
         y.min - 0.15,
         x[2],
         y.min - 0.4,
         col = '#b87624',
         border = NA)
    
    rect(x[14],
         y.min - 0.15,
         x[21],
         y.min - 0.4,
         col = '#b87624',
         border = NA)
    
    # legend for width
    points(1966, 7.5, pch = 15, col = '#b87624')
    text(1966.5, 7.5, '/')
    text(1967, 7.5, '(')
    points(1967.5, 7.5, pch = 15, col = '#b87624')
    text(1968.5, 7.5, '+')
    points(1969.5, 7.5, pch = 15, col = '#e9d5bd')
    text(1970, 7.5, ')')
    
    segments(1967, 6.7, 1969, 6.7, col = cols[panel], lty = 2)
    
    segments(1967, 5.9, 1969, 5.9, col = "#999999", lty = 1)
    points(1968, 5.9, pch = 1, col = "#999999")
    
    # text(1972, c(7.5, 6.7, 5.9),
    #      c('nWidth', 'Cutoff', 'Observation'),
    #      adj = c(0, 0.5))
    
    text(1972, c(7.5, 6.7, 5.9),
         c('nWidth = 0.3', 'Cutoff = 1:40', 'Observation'),
         adj = c(0, 0.5))
    
    
    # legend(x.min, y.max,
    #        legend = c("Width", "Cutoff", "Observation"),
    #        col = c(cols[panel], cols[panel], "#999999"),
    #        bty = "n",
    #        pch = c(15, NA, 1),
    #        lty = c(NA, 2, 1),
    #        xpd = TRUE)
    
    
  }  
  
  if(panel == 3){
    
    stat <- atyData(hi_paired, conversion)
    
    aty <- stat$ATY_V4[stat$participant_id == id]
    
    aty.mean <- mean(used$year_of_isolation)
    
    segments(aty, y.min, aty, y.max-2, col = cols[panel], lwd = 3)
    segments(aty.mean, y.min, aty.mean, y.max-2, col = cols[panel], lwd = 3, lty = 2)
    
    # legend for auc
    segments(1967, 7.5, 1969, 7.5, col = cols[panel], lty = 2)
    
    segments(1967, 6.7, 1969, 6.7, col = cols[panel], lty = 1)
    
    segments(1967, 5.9, 1969, 5.9, col = "#999999", lty = 1)
    points(1968, 5.9, pch = 1, col = "#999999")
    
    # text(1972, c(7.5, 6.7, 5.9),
    #      c('Unweighted average isolation year', 'nATY', 'Observation'),
    #      adj = c(0, 0.5))
    
    text(1972, c(7.5, 6.7, 5.9),
         c('Unweighted average isolation year = 1991.7', 
           'nATY = 1996.4', 
           'Observation'),
         adj = c(0, 0.5))
    
    # legend(x.min, y.max,
    #        legend = c("Unweighted average isolation year",
    #                   "ATY",
    #                   "Observation"),
    #        col = c(cols[panel], cols[panel], "#999999"),
    #        bty = "n",
    #        pch = c(NA, NA, 1),
    #        lty = c(2, 1, 1),
    #        lwd = 1,
    #        xpd = TRUE)
    
  }

  
  axis(1, seq(x.min, x.max, x.diff.1), 
       pos = y.min - 0.5,
       cex.axis = 1.5)
  axis(1, seq(x.min, x.max, x.diff.2), 
       pos = y.min - 0.5, labels = NA,
       tcl = -0.2, col.ticks = 1, col=NA,
       cex.axis = 1.5)
  
  axis(2, seq(y.min, y.max, 1), 
       pos = x.min - 0.5, 
       las=1,
       # labels = 2^seq(y.min, y.max, 1)*5,
       labels = seq(y.min, y.max, 1),
       cex.axis = 1.5)
  
  mtext(# c("a", "e", "i")[panel],
    c("A", "E", "I")[panel],
    side = 2,
    line = 4.5, 
    at = y.max+(y.max-y.min)*0.1,
    outer = FALSE, 
    las = 1,
    cex = 1.5,
    font = 2)
  
}


##'
##' Plot of auc of pre-birth strains
##' 
##' @param data, dataframe of conversion 
##' @param inset indicate whether to plot inset or not
##' @param panel, 1 = all titers, 2 = pre-titers <= 80, 3 = pre-titers > 80
##' 
##' @return panel B, C and D of Fig. 4
##'

aucGamDataBirthAdjustedPlot <- function(stat, m, panel, legend){
  
  stat <- stat[order(stat$age_at_sampling), ]
  
  x.min <- 0
  x.max <- 80
  y.min <- -5
  y.max <- 15
  
  cols <- c("#3f93aa",
            "#fa8072",
            "#590059")
  
  cols.lines <- c("#255866",
                  "#af594f",
                  "#350035")
  
  fit <- predict(m[[panel]], newdata = stat, se = TRUE)
  
  par(mar = c(6,6,5,2),
      mgp = c(4, 1, 4))
  plot(NULL, 
       xlim = c(x.min, x.max),
       ylim = c(y.min, y.max),
       xlab = "Age at baseline, years", 
       ylab = "nAUC",
       las = 1,
       type = "p",
       pch = 19,
       col = adjustcolor(#"#ff6666",
         cols[panel],
         alpha.f = 0.2),
       bty = "n",
       cex = 0.6,
       axes = FALSE,
       cex.lab = 1.5)
  
  # AUC data
  points(stat$age_at_sampling, 
         stat[ ,c("AUC_V1", "AUC_V4", "AUC_delta")[panel]],
         pch = 19,
         col = adjustcolor(#"#ff6666",
           cols[panel],
           alpha.f = 0.2),
         bty = "n",
         cex = 0.6)
  
  # AUC
  lines(stat$age_at_sampling, fit$fit, 
        col = cols.lines[panel],
        #col = grey(0.5), 
        pch = 19,
        lwd = 2.5)
  lines(stat$age_at_sampling, fit$fit-1.96*fit$se.fit, 
        col = cols.lines[panel],
        lty = 3,
        lwd = 2.5)
  lines(stat$age_at_sampling, fit$fit+1.96*fit$se.fit, 
        col = cols.lines[panel],
        lty = 3,
        lwd = 2.5)
  
  axis(1, at = seq(x.min, x.max, 20),
       pos = y.min - (y.max-y.min)*0.02,
       cex.axis = 1.5)
  axis(1, at = seq(x.min, x.max, 10),
       pos = y.min - (y.max-y.min)*0.02,
       tcl = -0.3, labels = NA, xpd = TRUE,
       col = NA, col.ticks = "black",
       cex.axis = 1.5)
  axis(1, at = seq(x.min, x.max, 1),
       pos = y.min - (y.max-y.min)*0.02,
       tcl = -0.1, labels = NA, xpd = TRUE,
       col = NA, col.ticks = "black",
       cex.axis = 1.5)
  
  
  mtext(# letters[panel + 1],
        # LETTERS[panel + 1],
        legend[panel], 
        side = 2,
        line = 4.5,
        at = y.max+(y.max-y.min)*0.1,
        cex = 1.5,
        font = 2,
        outer = FALSE,
        las = 1)
  
  
  axis(2, at = seq(y.min, y.max, 5),
       pos = x.min - (x.max-x.min)*0.02, las = 1,
       cex.axis = 1.5)
  
}


##' Function to plot width of pre-birth strains
##' 
##' @param stat dataframe of widths returned from widthData()
##' @param m list of fitted gam models returned from widthGam()
##' @param panel 
##' @param cutoff
##' @param spp, panel used in figure 3 or supp figure 
##'          
##' @return panel F, G and H of Fig. 4 
##'

widthGamDataBirthAdjustedPlot <- function(stat, m, panel, cutoff, spp, legend){
  
  stat <- stat[order(stat$age_at_sampling), ]
  
  x.min <- 0
  x.max <- 80
  y.min <- 0
  y.max <- 1
  
  cols <- c("#3f93aa",
            "#fa8072",
            "#590059")
  
  cols.lines <- c("#255866",
                  "#af594f",
                  "#350035")
  
  fit <- predict(m[[panel]], newdata = stat, se = TRUE)
  
  par(mar = c(6,6,5,2),
      mgp = c(4, 1, 4))
  plot(stat$age_at_sampling, 
       stat[ ,c("width_V1", "width_V4", "width_delta")[panel]], 
       xlim=c(x.min, x.max),
       ylim=c(y.min, y.max),
       xlab = "Age at baseline, years", 
       # ylab = paste("nWidth, cutoff at 1:", 2^cutoff*5, sep=""),
       ylab = paste("nWidth, 1:", 2^cutoff*5, sep=""),
       las = 1,
       type = "p",
       pch = 17,
       col = adjustcolor(#"#ff6666",
         cols[panel],
         alpha.f = 0.2),
       bty = "l",
       axes = FALSE,
       cex.lab = 1.5)
  
  lines(stat$age_at_sampling, fit$fit, 
        col = cols.lines[panel],
        lwd = 2.5)
  
  lines(stat$age_at_sampling, fit$fit-1.96*fit$se.fit, 
        col = cols.lines[panel],
        lty = 3,
        lwd = 2.5)
  
  lines(stat$age_at_sampling, fit$fit+1.96*fit$se.fit, 
        col = cols.lines[panel],
        lty = 3,
        lwd = 2.5)
  
  axis(1, at = seq(x.min, x.max, 20),
       pos = y.min - (y.max - y.min)*0.02,
       cex.axis = 1.5)
  axis(1, at = seq(x.min, x.max, 10),
       pos = y.min - (y.max - y.min)*0.02,
       tcl = -0.3, labels = NA, xpd = TRUE,
       col = NA, col.ticks = "black",
       cex.axis = 1.5)
  axis(1, at = seq(x.min, x.max, 1),
       pos = y.min - (y.max - y.min)*0.02,
       tcl = -0.1, labels = NA, xpd = TRUE,
       col = NA, col.ticks = "black",
       cex.axis = 1.5)
  
  axis(2, at = seq(y.min, y.max, 0.2),
       pos = x.min - (x.max - x.min)*0.02, las = 1,
       cex.axis = 1.5)
  # axis(2, at = seq(y.min, y.max, 0.1),
  #      pos = x.min - (x.max - x.min)*0.02, las = 1,
  #      tcl = -0.15, labels = NA, xpd = TRUE,
  #      col = NA, col.ticks = "black")
  
  if(spp){
    
    if(cutoff==1){
      mtext(# letters[panel],
            # LETTERS[panel],
            legend[panel],
            side = 2,
            line = 4.5,
            at = y.max+(y.max-y.min)*0.1,
            cex = 1.5,
            font = 2,
            outer = FALSE,
            las = 1)
    }
    if(cutoff==3){
      mtext(# letters[panel + 3],
            # LETTERS[panel + 3],
            legend[panel],
            side = 2,
            line = 4.5,
            at = y.max+(y.max-y.min)*0.1,
            cex = 1.5,
            font = 2,
            outer = FALSE,
            las = 1)
    }
    
  } else {
    mtext(# letters[panel + 9],
          # LETTERS[panel + 5],
          legend[panel],
          side = 2,
          line = 4.5,
          at = y.max+(y.max-y.min)*0.1,
          cex = 1.5,
          font = 2,
          outer = FALSE,
          las = 1)
  }
  
  
}


##' Function to plot aty of pre-birth strains
##' 
##' @param stat dataframe of widths returned from widthData()
##' @param m list of fitted gam models returned from widthGam()
##' @param panel panel of figure, visit 1, 4 or difference
##'          
##' @return panel J, K and L of Fig. 4
##'

atyGamDataBirthAdjustedPlot <- function(stat, m, panel, legend){
  
  stat <- stat[order(stat$age_at_sampling), ]
  
  x.min <- 0
  x.max <- 80
  y.min <- 1965
  y.max <- 2015
  
  cols <- c("#3f93aa",
            "#fa8072",
            "#590059")
  
  cols.lines <- c("#255866",
                  "#af594f",
                  "#350035")
  
  fit <- predict(m[[panel]], newdata = stat, se = TRUE)
  
  par(mar = c(6,6,5,2),
      mgp = c(4, 1, 4))
  plot(stat$age_at_sampling, 
       stat[ ,c("ATY_V1", "ATY_V4", "ATY_delta")[panel]],
       xlim=c(x.min, x.max),
       ylim=c(y.min, y.max),
       xlab = "Age at baseline, years", 
       ylab = "nATY",
       las = 1,
       bty = "l",
       type = "p",
       axes = FALSE,
       cex.lab = 1.5,
       pch = 18,
       col = adjustcolor(#"#ff6666",
         cols[panel],
         alpha.f = 0.2))
  
  
  lines(stat$age_at_sampling, fit$fit, 
        col = cols.lines[panel],
        lwd = 2.5)
  lines(stat$age_at_sampling, fit$fit-1.96*fit$se.fit, 
        col = cols.lines[panel],
        lty = 3,
        lwd = 2.5)
  lines(stat$age_at_sampling, fit$fit+1.96*fit$se.fit, 
        col = cols.lines[panel],
        lty = 3,
        lwd = 2.5)
  
  # unweighted average isolation year
  age <- 0:80
  uaty <- sapply(age, function(a){
    
    yob <- 2010 - a
    yr <- st$year_of_isolation[st$year_of_isolation >= yob]
    mean(yr)
    
  })
  
  lines(age, uaty, lty = 2, lwd = 2)
  
  # yob
  y.end <- 2010
  yrs <- 1965:y.end
  lines(y.end - yrs,
        yrs,
        lty = 3,
        lwd = 2)
  
  axis(1, at = seq(x.min, x.max, 20),
       pos = y.min - (y.max - y.min)*0.02,
       cex.axis = 1.5)
  axis(1, at = seq(x.min, x.max, 10),
       pos = y.min - (y.max - y.min)*0.02,
       tcl = -0.3, labels = NA, xpd = TRUE,
       col = NA, col.ticks = "black",
       cex.axis = 1.5)
  axis(1, at = seq(x.min, x.max, 1),
       pos = y.min - (y.max - y.min)*0.02,
       tcl = -0.1, labels = NA, xpd = TRUE,
       col = NA, col.ticks = "black",
       cex.axis = 1.5)
  
  axis(2, at = seq(y.min, y.max, 10),
       pos = x.min - (x.max - x.min)*0.02, las = 1,
       cex.axis = 1.5)
  axis(2, at = seq(y.min, y.max, 5),
       pos = x.min - (x.max - x.min)*0.02, las = 1,
       tcl = -0.3, labels = NA, xpd = TRUE,
       col = NA, col.ticks = "black",
       cex.axis = 1.5)
  axis(2, at = seq(y.min, y.max, 1),
       pos = x.min - (x.max - x.min)*0.02, las = 1,
       tcl = -0.1, labels = NA, xpd = TRUE,
       col = NA, col.ticks = "black",
       cex.axis = 1.5)
  
  
  mtext(# letters[panel + 5],
        # LETTERS[panel + 9],
        legend[panel],
        side = 2,
        line = 4.5,
        at = y.max+(y.max-y.min)*0.1,
        cex = 1.5,
        font = 2,
        outer = FALSE,
        las = 1)
  
}


#------------------------------- Figure 5 ----------------------------------------%

##' 
##' Function to plot changes in strain-specific GMT
##' 
##' @param conversion
##' @param st
##' @param all if all strains included
##' @param panel
##' @param supp if supp figure or not
##'          
##' @return panel A of Fig.5
##'


propChangePlot <- function(conversion, st, all, panel, supp){
  
  dat <- getPropChange(conversion, st)
  
  if(supp){
    
    if(all){
      dat.used <- dat[[2]]
      txt <- 'All strains'
    } else {dat.used <- dat[[1]]; txt <- 'Post-birth strains'}
    
  } else {
    
    dat.used <- dat[[2]]
    txt <- ''
    
  }
  
  x.lim <- c(0, 21)
  y.lim <- c(1, 6)
  space <- 0.05
  
  par(mar = c(0,6,3,3),
      mgp = c(3, 1, 0))
  plot(NULL,
       las = 1,
       xlim = x.lim,
       ylim = y.lim,
       xlab = '',
       ylab = 'GMT',
       axes = FALSE,
       cex.lab = 1.5)
  
  for(i in 1:2){ #visit
    
    used <- dat.used[[i]]
    
    cols <- c("#3f93aa",
              "#fa8072")[i]
    
    points(seq(0.5, x.lim[2], 1),
           used[1, ],
           pch = c(16,17)[i], 
           cex = 1,
           col = cols)
    
    lines(seq(0.5, x.lim[2], 1),
          used[1, ], 
          col = cols,
          lty = c(2,3)[i])
    
  }
  
  
  axis(1, 
       at = seq(0, x.lim[2], 1), 
       pos = y.lim[1],
       labels = NA,
       xpd = TRUE)
  
  y <- 2^seq(y.lim[1]+1, y.lim[2], 1)*2.5
  axis(2, 
       at = seq(1, y.lim[2], 1), 
       # pos = y.lim[1] - diff(y.lim)*space,
       pos = x.lim[1],
       labels = c(0, y),
       cex.axis = 1.5, las = 1)
  
  mtext(# letters[panel], 
    LETTERS[panel],
    side = 2,
    line = 4.5, 
    at = y.lim[2] + diff(y.lim)*0.08,
    outer = FALSE, 
    las = 1,
    cex = 2,
    font = 2)
  
}


##' 
##' Function to plot strain distribution by titer changes
##' 
##' @param data conversion
##' @param st dataframe contains strain information
##' @param panel, 1 = all titers
##'               2 = pre-existing titer <= 1:80
##'               3 = pre-existing titer > 1:80
##'          
##' @return panel B of Fig.5
##'

strainDistrByTiterChanges <- function(data, st, panel){
  
  data <- data[!is.na(data$fold_of_change), ]
  data$change <- 0
  data$change[data$fold_of_change < 0] <- -1
  data$change[data$fold_of_change > 0] <- 1
  data$change[data$seroconversion == 1] <- 2
  
  p <- t(sapply(st$strain, function(s){
    
    used <- data[data$strain == s, ]
    
    table(used$change)/nrow(used)
    
  }))
  
  p <- data.frame(p)
  colnames(p) <- c('decrease', 'unchanged', 'twoFoldIncrease', 'seroconversion')
  p$increase <- p$`twoFoldIncrease` + p$seroconversion
  
  
  y.lim <- c(0, 1)
  x.lim <- c(0, nrow(st))
  
  cols <- c('#8BAD23', grey(0.9), '#af84af', '#590059')
  cols.seg <- c('#8BAD23', grey(0.9), '#af84af', '#590059')
  line.lwd = 1.5
  
  par(mar = c(12,6,3,3),
      mgp = c(3, 1, 0))
  plot(NULL,
       xlim = x.lim, ylim = y.lim,
       axes = FALSE,
       xlab = "", 
       ylab = "Proportion, %",
       cex.lab = 1.5)
  
  
  # decrease
  rect(xleft = seq(x.lim[1] + 0.2, x.lim[2], 1),
       ybottom = rep(0, x.lim[2]),
       xright = seq(x.lim[1] + 0.8, x.lim[2], 1),
       ytop = p$decrease,
       col = cols[1],
       border = NA)
  
  # no change
  rect(xleft = seq(x.lim[1] + 0.2, x.lim[2], 1),
       ybottom = p$decrease,
       xright = seq(x.lim[1] + 0.8, x.lim[2], 1),
       ytop = p$decrease + p$unchanged,
       col = cols[2],
       border = NA)
  
  # increase
  rect(xleft = seq(x.lim[1] + 0.2, x.lim[2], 1),
       ybottom = p$decrease + p$unchanged,
       xright = seq(x.lim[1] + 0.8, x.lim[2], 1),
       ytop = p$decrease + p$unchanged + p$twoFoldIncrease,
       col = cols[3],
       border = NA)
  
  # four fold increase
  rect(xleft = seq(x.lim[1] + 0.2, x.lim[2], 1),
       ybottom = p$decrease + p$unchanged + p$twoFoldIncrease,
       xright = seq(x.lim[1] + 0.8, x.lim[2], 1),
       ytop = rep(1, x.lim[2]),
       col = cols[4],
       border = NA)
  
  strain <- st$strain
  strain[strain == 'X31'] <- 'X31 (1970)'
  text(seq(x.lim[1] + 0.5, x.lim[2], 1),
       rep(-0.1, y.lim[2]),
       strain,
       cex = 1.5,
       srt = 90,
       adj = c(1, 0.5),
       xpd = TRUE)
  
  mtext(# letters[panel], 
    LETTERS[panel],
    side = 2,
    line = 4.5, 
    at = y.lim[2] + diff(y.lim)*0.08,
    outer = FALSE, 
    las = 1,
    cex = 2,
    font = 2)
  
  
  axis(2, 
       at = seq(y.lim[1], y.lim[2], 0.1),
       pos = x.lim[1],
       labels = seq(y.lim[1], y.lim[2], 0.1)*100,
       cex.axis = 1.5, las = 1)
  
  axis(1, 
       at = seq(x.lim[1], x.lim[2], 1),
       pos = y.lim[1], labels = NA)
  
  
}


#------------------------------- Figure S2 ----------------------------------------%


##'
##' Plot of auc using all strains
##' 
##' @param stat, dataframe contains auc, returned from aucData()
##' @param m list of fitted GAM to auc, returned from aucGam()
##' @param panel, 1 = baseline 
##'               2 = follow-up
##'               3 = changes
##' @param all logical, if all strains included
##' @param group number of subgroups, 1 to 5, valid only when all = FALSE
##' @param combinedPlot logical, if belongs to combined plot
##' 
##' @return panel A to C of Fig. s2
##'

aucGamDataPlot <- function(stat, m, panel, all, group, combinedPlot){
  
  stat <- stat[order(stat$age_at_sampling), ]
  
  x.min <- 0
  x.max <- 80
  y.min <- c(-100, -100, -100)[panel]
  y.max <- c(300, 300, 300)[panel]
  if(all == FALSE){
    y.min <- -60
    y.max <- 60
  }
  
  cols <- c("#3f93aa",
            "#fa8072",
            "#590059")
  
  cols.lines <- c("#255866",
                  "#af594f",
                  "#350035")
  
  fit <- predict(m[[panel]], newdata = stat, se = TRUE)
  
  par(mar = c(6,6,5,2),
      mgp = c(4, 1, 4))
  plot(NULL, 
       xlim = c(x.min, x.max),
       ylim = c(y.min, y.max),
       xlab = "Age at sampling, years", 
       ylab = "AUC",
       las = 1,
       type = "p",
       pch = 19,
       col = adjustcolor(#"#ff6666",
         cols[panel],
         alpha.f = 0.2),
       bty = "n",
       cex = 0.6,
       axes = FALSE,
       cex.lab = 1.5)
  
  if(all == FALSE){
    
    # birth cohort
    y1 <- c(1968, 1982, 1995, 2007)[group]
    y2 <- c(1979, 1992, 2004, 2014)[group]
    
    if(panel == 2){
      rect(max(2014-y1, x.min),
           y.min,
           max(2014-y2, x.min),
           y.max,
           col = adjustcolor(grey(0.7),
                             alpha.f = 0.4),
           border = NA)
    } else {
      rect(max(2010-y1, x.min),
           y.min,
           max(2010-y2, x.min),
           y.max,
           col = adjustcolor(grey(0.7),
                             alpha.f = 0.4),
           border = NA)
    }
    
  }
  
  
  # AUC data
  points(stat$age_at_sampling, 
         stat[ ,c("AUC_V1", "AUC_V4", "AUC_delta")[panel]],
         pch = 19,
         col = adjustcolor(cols[panel],
                           alpha.f = 0.2),
         bty = "n",
         cex = 0.6)
  
  # AUC
  lines(stat$age_at_sampling, fit$fit, 
        col = cols.lines[panel],
        #col = grey(0.5), 
        pch = 19,
        lwd = 2.5)
  lines(stat$age_at_sampling, fit$fit-1.96*fit$se.fit, 
        col = cols.lines[panel],
        lty = 3,
        lwd = 2.5)
  lines(stat$age_at_sampling, fit$fit+1.96*fit$se.fit, 
        col = cols.lines[panel],
        lty = 3,
        lwd = 2.5)
  
  axis(1, at = seq(x.min, x.max, 20),
       pos = y.min - (y.max-y.min)*0.02,
       cex.axis = 1.5)
  axis(1, at = seq(x.min, x.max, 10),
       pos = y.min - (y.max-y.min)*0.02,
       tcl = -0.3, labels = NA, xpd = TRUE,
       col = NA, col.ticks = "black",
       cex.axis = 1.5)
  axis(1, at = seq(x.min, x.max, 1),
       pos = y.min - (y.max-y.min)*0.02,
       tcl = -0.1, labels = NA, xpd = TRUE,
       col = NA, col.ticks = "black",
       cex.axis = 1.5)
  
  
  if(all==TRUE){
    if(combinedPlot){
      mtext(# letters[panel],
            LETTERS[panel],
            side = 2,
            line = 3.5,
            at = y.max+(y.max-y.min)*0.1,
            cex = 1.5,
            font = 2,
            outer = FALSE,
            las = 1)
    } else {
      mtext(# letters[panel],
            LETTERS[panel],
            side = 2,
            line = 4.5,
            at = y.max+(y.max-y.min)*0.1,
            cex = 1.5,
            font = 2,
            outer = FALSE,
            las = 1)
    }
    
    
    axis(2, at = seq(y.min, y.max, 100),
         pos = x.min - (x.max-x.min)*0.02, las = 1,
         cex.axis = 1.5)
    
  } else {
    
    if(panel == 1){
      l.used <- group
      mtext(# letters[l.used],
            LETTERS[l.used],
            side = 2,
            line = 4.5,
            at = y.max+(y.max-y.min)*0.1,
            cex = 1.5,
            font = 2,
            outer = FALSE,
            las = 1)
      
    }
    if(panel == 2){
      txt <- c("Strains from 1968 to 1977",
               "Strains from 1979 to 1989",
               "Strains from 1992 to 2002",
               "Strains from 2004 to 2014")[group]
      title(txt, cex.main = 1.5)
    }
    
    axis(2, at = seq(y.min, y.max, 30),
         pos = x.min - (x.max-x.min)*0.02, las = 1,
         cex.axis = 1.5)
    
  }
  
}



##' 
##' Function to plot fitted gam of width
##' 
##' @param stat dataframe of widths returned from widthData()
##' @param m list of fitted GAM to width, returned from widthGam()
##' @param panel, 1 = baseline 
##'               2 = follow-up
##'               3 = changes
##' @param cutoff 1 = 1:10
##'               3 = 1:40
##'          
##' @return panel D to I of Fig. s2
##'

widthGamDataPlot <- function(stat, m, panel, cutoff){
  
  stat <- stat[order(stat$age_at_sampling), ]
  
  x.min <- 0
  x.max <- 80
  y.min <- 0
  y.max <- 1
  
  cols <- c("#3f93aa",
            "#fa8072",
            "#590059")
  
  cols.lines <- c("#255866",
                  "#af594f",
                  "#350035")
  
  fit <- predict(m[[panel]], newdata = stat, se = TRUE)
  
  par(mar = c(6,6,5,2),
      mgp = c(4, 1, 4))
  plot(stat$age_at_sampling, 
       stat[ ,c("width_V1", "width_V4", "width_delta")[panel]], 
       xlim=c(x.min, x.max),
       ylim=c(y.min, y.max),
       xlab = "Age at sampling, years", 
       ylab = paste("Width, cutoff at 1:", 2^cutoff*5, sep=""),
       las = 1,
       type = "p",
       pch = 17,
       col = adjustcolor(#"#ff6666",
         cols[panel],
         alpha.f = 0.2),
       bty = "l",
       axes = FALSE,
       cex.lab = 1.5)
  
  lines(stat$age_at_sampling, fit$fit, 
        col = cols.lines[panel],
        lwd = 2.5)
  
  lines(stat$age_at_sampling, fit$fit-1.96*fit$se.fit, 
        col = cols.lines[panel],
        lty = 3,
        lwd = 2.5)
  
  lines(stat$age_at_sampling, fit$fit+1.96*fit$se.fit, 
        col = cols.lines[panel],
        lty = 3,
        lwd = 2.5)
  
  axis(1, at = seq(x.min, x.max, 20),
       pos = y.min - (y.max - y.min)*0.02,
       cex.axis = 1.5)
  axis(1, at = seq(x.min, x.max, 10),
       pos = y.min - (y.max - y.min)*0.02,
       tcl = -0.3, labels = NA, xpd = TRUE,
       col = NA, col.ticks = "black",
       cex.axis = 1.5)
  axis(1, at = seq(x.min, x.max, 1),
       pos = y.min - (y.max - y.min)*0.02,
       tcl = -0.1, labels = NA, xpd = TRUE,
       col = NA, col.ticks = "black",
       cex.axis = 1.5)
  
  axis(2, at = seq(y.min, y.max, 0.2),
       pos = x.min - (x.max - x.min)*0.02, las = 1,
       cex.axis = 1.5)

  if(cutoff==1){
    mtext(letters[panel + 6],
          # LETTERS[panel + 6],
          side = 2,
          line = 4.5,
          at = y.max+(y.max-y.min)*0.1,
          cex = 1.5,
          font = 2,
          outer = FALSE,
          las = 1)
  }
  
  if(cutoff==3){
    mtext(letters[panel + 9],
          # LETTERS[panel + 3],
          side = 2,
          line = 4.5,
          at = y.max+(y.max-y.min)*0.1,
          cex = 1.5,
          font = 2,
          outer = FALSE,
          las = 1)
  }
  
  
  
}


##' 
##' Function to plot fitted gam of aty
##' 
##' @param stat dataframe of aty returned from atyData()
##' @param m list of fitted GAM to aty, returned from atyGam()
##' @param panel, 1 = baseline 
##'               2 = follow-up
##'               3 = changes
##'          
##' @return panel J to L of Fig. s2
##'

atyGamDataPlot <- function(stat, m, panel){
  
  stat <- stat[order(stat$age_at_sampling), ]
  
  x.min <- 0
  x.max <- 80
  y.min <- 1965
  y.max <- 2015
  
  cols <- c("#3f93aa",
            "#fa8072",
            "#590059")
  
  cols.lines <- c("#255866",
                  "#af594f",
                  "#350035")
  
  fit <- predict(m[[panel]], newdata = stat, se = TRUE)
  
  par(mar = c(6,6,5,2),
      mgp = c(4, 1, 4))
  plot(stat$age_at_sampling, 
       stat[ ,c("ATY_V1", "ATY_V4", "ATY_delta")[panel]], 
       xlim=c(x.min, x.max),
       ylim=c(y.min, y.max),
       xlab = "Age at sampling, years", 
       ylab = "Average titer years",
       las = 1,
       type = "p",
       pch = 18,
       col = adjustcolor(#"#ff6666",
         cols[panel],
         alpha.f = 0.2),
       bty = "l",
       axes = FALSE,
       cex.lab = 1.5)
  
  lines(stat$age_at_sampling, fit$fit, 
        col = cols.lines[panel],
        lwd = 2.5)
  lines(stat$age_at_sampling, fit$fit-1.96*fit$se.fit, 
        col = cols.lines[panel],
        lty = 3,
        lwd = 2.5)
  lines(stat$age_at_sampling, fit$fit+1.96*fit$se.fit, 
        col = cols.lines[panel],
        lty = 3,
        lwd = 2.5)
  
  axis(1, at = seq(x.min, x.max, 20),
       pos = y.min - (y.max - y.min)*0.02,
       cex.axis = 1.5)
  axis(1, at = seq(x.min, x.max, 10),
       pos = y.min - (y.max - y.min)*0.02,
       tcl = -0.3, labels = NA, xpd = TRUE,
       col = NA, col.ticks = "black",
       cex.axis = 1.5)
  axis(1, at = seq(x.min, x.max, 1),
       pos = y.min - (y.max - y.min)*0.02,
       tcl = -0.1, labels = NA, xpd = TRUE,
       col = NA, col.ticks = "black",
       cex.axis = 1.5)
  
  axis(2, at = seq(y.min, y.max, 10),
       pos = x.min - (x.max - x.min)*0.02, las = 1,
       cex.axis = 1.5)
  axis(2, at = seq(y.min, y.max, 5),
       pos = x.min - (x.max - x.min)*0.02, las = 1,
       tcl = -0.3, labels = NA, xpd = TRUE,
       col = NA, col.ticks = "black",
       cex.axis = 1.5)
  axis(2, at = seq(y.min, y.max, 1),
       pos = x.min - (x.max - x.min)*0.02, las = 1,
       tcl = -0.1, labels = NA, xpd = TRUE,
       col = NA, col.ticks = "black",
       cex.axis = 1.5)
  
  # age reference lines
  y.end <- 2010
  yrs <- 1965:y.end
  lines(y.end - yrs,
        yrs,
        lty = 3)
  
  # average reference lines
  
  y.used <- c(1990.556,
              1992.8,
              1992.8)[panel]
  segments(x.min, y.used, x.max, y.used, lty = 2)
  
  mtext(letters[panel + 9],
        # LETTERS[panel + 9],
        side = 2,
        line = 4.5,
        at = y.max+(y.max-y.min)*0.1,
        cex = 1.5,
        font = 2,
        outer = FALSE,
        las = 1)
  
}


#------------------------------- Figure S3 ----------------------------------------%


##' 
##' Function to plot OR of seroconversion for strains 
##' 
##' @param m fitted logistic regression model
##'          
##' @return figure s4
##'

seroconvORPlot <- function(m){
  
  coeff <- c(1, exp(coef(m)[4:23]))
  ci <- exp(confint(m)[4:23, ])
  ci <- rbind(c(NA, NA), ci)
  
  x.lim <- c(0, 21)
  y.lim <- c(0, 25)
  x.pos <- 0.5:21
  
  par(mar=c(8,5,2,2))
  plot(NULL,
       xlim = x.lim, ylim = y.lim,
       axes = FALSE,
       xlab = "", 
       ylab = "Odds ratio")
  
  segments(x.lim[1], 1, x.lim[2], 1, lty = 3)
  
  points(x.pos, coeff, pch = 19)
  segments(x.pos, ci[,1], x.pos, ci[,2])
  segments(x.pos-0.15, ci[,1], x.pos+0.15, ci[,1])
  segments(x.pos-0.15, ci[,2], x.pos+0.15, ci[,2])
  
  axis(1, at = seq(x.lim[1], x.lim[2]), labels = NA, pos = y.lim[1])
  axis(2, at = seq(y.lim[1], y.lim[2], 5), pos = x.lim[1], las = 1)
  
  text(x.pos, y.lim[1]-0.05,
       st$strain,
       srt = 90,
       xpd = TRUE,
       adj = c(1,0.5))
  
}


#------------------------------- Figure S5 ----------------------------------------%


##' 
##' Function to distribution of changes for four recent strains 
##' 
##' @param changes a list contains information of changes distribution
##'                returned from getDataForTiterChanges()
##' @param byStrain TRUE = plot distribution of changes by each strain (panel B)
##'                 FALSE = plot distribution of changes by # of strains with increase (panel A)            
##'          
##' @return figure s5
##'

plotChangesInFourStrains <- function(changes, byStrain){
  
  if(byStrain){
    dat <- changes$byStrain
    xlab.txt <- ""
  } else {
    dat <- changes$byNumber
    xlab.txt <- "Number of recent strains with increased titers"
  }
  
  cols <- c('#8BAD23', grey(0.8), '#CDB2CD', '#590059')
  
  
  x.lim <- c(0, nrow(dat))
  
  if(byStrain){
    
    y.lim <- c(0, 1)
    
    par(mar = c(12,6,3,3))
    
    plot(NULL,
         xlim = x.lim,
         ylim = y.lim,
         las = 1,
         xlab = xlab.txt,
         ylab = 'Proportion, %',
         axes = FALSE,
         cex.lab = 1.5,
         cex.main = 1.5)
    
    for(i in 1:x.lim[2]){
      
      used <- c(0, cumsum(dat[i, 2:5]))
      
      sapply(1:4, function(x){
        
        rect(i + 0.3 - 1,
             used[x],
             i + 0.7 -1,
             used[x+1],
             col = cols[x],
             border = NA)
        
      })
      
    }
    
    axis(1, seq(x.lim[1], x.lim[2], 1), 
         pos = y.lim[1] - 0.02, labels = NA)
    
    axis(2, seq(y.lim[1], y.lim[2], 0.2), 
         pos = x.lim[1], las = 1, cex.axis = 1.5)
    
    text(seq(x.lim[1] + 0.5, x.lim[2], 1), 
         rep(-0.1, x.lim[2]), 
         c("A/Perth/2009",
           "A/Victoria/2009",
           "A/Texas/2012",
           "A/HongKong/2014"),
         xpd = TRUE, cex = 1.5,
         srt = 90, c(1,0.5))
    
    mtext(# letters[2], 
      LETTERS[2],
      side = 2,
      # line = 6, 
      line = 4, 
      at = y.lim[2],
      outer = FALSE, 
      las = 1,
      cex = 2,
      font = 2)
    
  }
  
  if(!byStrain){
    
    y.lim <- c(0, 0.6)
    
    par(mar = c(12,6,3,3))
    
    plot(NULL,
         xlim = x.lim,
         ylim = y.lim,
         las = 1,
         xlab = xlab.txt,
         ylab = 'Proportion, %',
         axes = FALSE,
         cex.lab = 1.5,
         cex.main = 1.5)
    
    for(i in 1:x.lim[2]){
      
      used <- c(0, cumsum(dat[i, 2:5])) * dat[i, 1]/sum(dat[ ,1])
      
      sapply(1:4, function(x){
        
        rect(i + 0.3 - 1,
             used[x],
             i + 0.7 -1,
             used[x+1],
             col = cols[x],
             border = NA)
        
      })
      
    }
    
    axis(1, seq(x.lim[1], x.lim[2], 1), 
         pos = y.lim[1] - 0.02, labels = NA)
    
    axis(2, seq(y.lim[1], y.lim[2], 0.2), 
         pos = x.lim[1], las = 1, cex.axis = 1.5)
    
    text(seq(x.lim[1] + 0.5, x.lim[2], 1), 
         rep(-0.1, x.lim[2]), 
         x.lim[1]:(x.lim[2]-1),
         xpd = TRUE, cex = 1.5)
    
    mtext(# letters[1], 
      LETTERS[1],
      side = 2,
      line = 4, 
      at = y.lim[2],
      outer = FALSE, 
      las = 1,
      cex = 2,
      font = 2)
    
    
  }
  
  
  
}


##' 
##' Function to plot AIC/BIC
##' 
##' @param dat data used to assess seroconversion, returned from getSeroconversionData()
##'          
##' @return figure s6 and s7
##'
##'

AICComparePlot <- function(dat, spline){
  
  if(!spline){
    tab <- AICCompareTable(dat)
  }
  if(spline){
    tab <- AICCompareTableSpline(dat)
  }
  
  tab.sub <- tab[tab$range != "NA", ]
  
  tab.sub$metric <- factor(tab.sub$metric,
                           levels = c("AUC",  "W40", "W10", "ATY"))
  tab.sub$range <- factor(tab.sub$range,
                          levels = c("1968", "adj"))
  
  tab.sub <- tab.sub[order(tab.sub$metric), ]
  
  layout(matrix(1:10, 
                nrow = 2, ncol = 5,
                byrow = TRUE), 
         widths = rep(c(1.2, rep(1, 4)), 2),
         heights = c(3, 1))
  
  sapply(0:9, AICComparePlotShortPanel, used = tab.sub, tab = tab)
  
  
}


#------------------------------ Figure S6 & S13 -------------------------------------%

##' 
##' Function to plot panel in AIC/BIC figure
##' 
##' @param panel  0 = labels for metrics
##'               1 = A/Perth/2009 models including metrics
##'               2 = A/Victoria/2009 models including metrics
##'               3 = A/Texas/2012 models including metrics
##'               4 = A/HongKong/2014 models including metrics
##'               
##'               5 = table legend
##'               6 = A/Perth/2009 models without metrics
##'               7 = A/Victoria/2009 models without metrics
##'               8 = A/Texas/2012 models without metrics
##'               9 = A/HongKong/2014 models without metrics
##'               
##' @param used internal subset of table contains AIC/BIC
##' @param tab table contains AIC/BIC
##' 
##'          
##' @return panel in figure s6 and s13
##'
##'

AICComparePlotShortPanel <- function(panel, used, tab){
  
  
  x.lim <- c(c(960,   880,  900, 1000)[panel],
             c(1020,  940,  960, 1060)[panel])
  x.diff <- c( 30,    30,   30,  30)[panel]
  if(panel == 0){
    x.lim <- c(0, 1)
  }
  y.lim <- c(0, 8)
  
  y.pos <- rev(seq(y.lim[1] + 0.4, y.lim[2], 1))
  y.pos.2 <- rev(seq(y.lim[1] + 0.6, y.lim[2], 1))
  
  cols <- c("#f1bc50", "#005668")
  metric <- c("AUC",  "W40", "W10", "ATY")
  
  
  if(panel == 0){
    
    par(mar = c(5,1,8,1))
    plot(NULL,
         xlim = x.lim,
         ylim = y.lim,
         axes = FALSE,
         xlab = "",
         ylab = "")
    
    text(rep(1, nrow(used)/2),
         y.pos + 0.1,
         c("All strains",
           "Post-birth strains"),
         xpd = TRUE,
         cex = 1.5,
         adj = c(1, 0.5))
    
    segments(rep(0.15, 2),
             seq(y.lim[1] + 0.5, y.lim[2], 2),
             rep(0.15, 2),
             seq(y.lim[1] + 1.5, y.lim[2], 2))
    
    segments(rep(0.15, 2),
             seq(y.lim[1] + 0.5, y.lim[2], 2),
             rep(0.18, 2),
             seq(y.lim[1] + 0.5, y.lim[2], 2))
    
    segments(rep(0.15, 2),
             seq(y.lim[1] + 1.5, y.lim[2], 2),
             rep(0.18, 2),
             seq(y.lim[1] + 1.5, y.lim[2], 2))
    
    text(rep(0.1, 2),
         rev(seq(y.lim[1] + 1, y.lim[2], 2)),
         metric,
         xpd = TRUE,
         adj = c(1, 0.5),
         cex = 1.5)
    
  }
  
  if(panel > 0 & panel < 5){
    
    par(mar = c(5,1,8,1))
    plot(NULL,
         xlim = x.lim,
         ylim = y.lim,
         axes = FALSE,
         xlab = "",
         ylab = "")
    
    rect(x.lim[1],
         y.lim[1],
         x.lim[2],
         y.lim[1] + 2,
         col = grey(0.9),
         border = NA)
    
    rect(x.lim[1],
         y.lim[1] + 4,
         x.lim[2],
         y.lim[1] + 6,
         col = grey(0.9),
         border = NA)
    
    used <- used[used$range %in% c("1968", "adj"), ]
    
    # aic
    points(used[ ,panel*2-1], y.pos, pch = 19, col = cols[1], cex = 1.5, xpd = TRUE)
    # bic
    points(used[ ,panel*2], y.pos.2, pch = 17, col = cols[2], cex = 1.5, xpd = TRUE)
    
    segments(tab[3, panel*2-1],
             y.lim[1],
             tab[3, panel*2-1],
             y.lim[2],
             lty = 2,
             col = cols[1],
             lwd = 1.5)
    
    segments(tab[3, panel*2],
             y.lim[1],
             tab[3, panel*2],
             y.lim[2],
             lty = 2,
             col = cols[2],
             lwd = 1.5)
    
    segments(tab[4, panel*2-1],
             y.lim[1],
             tab[4, panel*2-1],
             y.lim[2],
             lty = 3,
             col = cols[1],
             lwd = 1.5)
    
    segments(tab[4, panel*2],
             y.lim[1],
             tab[4, panel*2],
             y.lim[2],
             lty = 3,
             col = cols[2],
             lwd = 1.5)
    
    
    axis(1, seq(x.lim[1], x.lim[2], x.diff), pos = y.lim[1], cex.axis = 1.5)
    axis(3, seq(x.lim[1], x.lim[2], x.diff), pos = y.lim[2], cex.axis = 1.5)
    
    axis(1, seq(x.lim[1], x.lim[2], 5), 
         pos = y.lim[1], cex.axis = 1.5,
         labels = NA, tick = TRUE,
         tcl = -0.2, col.ticks = 1, col = NA)
    axis(3, seq(x.lim[1], x.lim[2], 5), 
         pos = y.lim[2], cex.axis = 1.5,
         labels = NA, tick = TRUE,
         tcl = -0.2, col.ticks = 1, col = NA)
    
    mtext("AIC", 1, line = 2, cex = 1)
    mtext("BIC", 3, line = 2, cex = 1)
    
    mtext(st$strain[17+panel], 
          3, line = 4, cex = 1.5, font = 2)
  }
  
  if(panel >= 5){
    
    x.lim <- c(-0, 1)
    y.lim <- c(0, 1)
    
    y.pos <- seq(y.lim[2], y.lim[1]+0.2, -0.2)
    x.pos.1 <- 0.2
    x.pos.2 <- 0.6
    
    par(mar = c(1,1,1,1))
    plot(NULL,
         xlim = x.lim,
         ylim = y.lim,
         xlab = "",
         ylab = "",
         axes = FALSE)
    segments(0.2,
             mean(y.pos[1:2]),
             0.9,
             mean(y.pos[1:2]))
    
    if(panel == 5){
      text(rep(x.pos.1, 5),
           y.pos,
           c("Predictors",
             'None',
             "Age",
             "Age + Titer to strain i",
             "Age + Titer to strain i & i - 1"),
           xpd = TRUE,
           adj = c(0, 0.5),
           cex = 1.5,
           font = 2)
      
    } else {
      
      k <- panel - 5
      
      text(c(x.pos.1, x.pos.2),
           rep(y.pos[1], 2),
           c("AIC", "BIC"),
           xpd = TRUE,
           adj = c(0, 0.5),
           cex = 1.5,
           font = 2)
      text(rep(x.pos.1, 4),
           y.pos[2:5],
           sprintf("%.01f", tab[1:4 ,k*2-1]),
           xpd = TRUE,
           adj = c(0, 0.5),
           cex = 1.5)
      text(rep(x.pos.2, 4),
           y.pos[2:5],
           sprintf("%.01f", tab[1:4 ,k*2]),
           xpd = TRUE,
           adj = c(0, 0.5),
           cex = 1.5)
      
    }
    
    
  }
  
}


#------------------------------- Figure S8 ----------------------------------------%

##' 
##' Function to plot mediation analysis
##' 
##' @param medEff, results of mediation analysis that did not assume interaction
##'                between pre-existing titer and AUC/W40
##' @param medEffInteract, results of mediation analysis that assumes interaction
##'                        between pre-existing titer and AUC/W40
##' @param panel 0 = Legen texts
##'               1 = A/Perth/2009 models including metrics
##'               2 = A/Victoria/2009 models including metrics
##'               3 = A/Texas/2012 models including metrics
##'               4 = A/HongKong/2014 models including metrics
##' @param st, strain information              
##'          
##' @return figure s8
##'

medEffPlot <- function(medEff, medEffInteract, panel, st){
  
  y.pos <- list(7:5,
                3:1)
  
  if(panel == 0){
    
    x.lim <- c(0, 1)
    y.lim <- c(0, 7)
    space <- 0.05
    
    par(3.5, 0, 3, 0)
    plot(NULL,
         xlim = x.lim,
         ylim = y.lim,
         las = 1, 
         xlab = '',
         ylab = '',
         cex.lab = 1.5,
         axes = FALSE)
    
    text(1,
         c(y.pos[[1]], y.pos[[2]]),
         c('Indirect effect',
           'Direct effect',
           'Total effect'),
         xpd = TRUE,
         cex = 1.5,
         adj = c(1, 0.5))
    
    segments(0.25,
             c(y.pos[[1]][1], y.pos[[2]][1]),
             0.25,
             c(y.pos[[1]][3], y.pos[[2]][3]))
    
    segments(0.25,
             c(y.pos[[1]][1], y.pos[[2]][1]),
             0.28,
             c(y.pos[[1]][1], y.pos[[2]][1]))
    
    segments(0.25,
             c(y.pos[[1]][3], y.pos[[2]][3]),
             0.28,
             c(y.pos[[1]][3], y.pos[[2]][3]))
    
    text(rep(0.15, 2),
         c(y.pos[[1]][2], y.pos[[2]][2]),
         c('AUC',
           'Width, 1:40'),
         xpd = TRUE,
         adj = c(1, 0.5),
         cex = 1.5)
    
    
  } else {
    
    x.lim <- c(0.5, 1.5)
    y.lim <- c(0, 7)
    space <- 0.05
    
    rows <- list(seq(1, 9, 3),
                 seq(2, 9, 3),
                 seq(3, 9, 3))
    
    par(3.5, 2, 3, 2)
    plot(NULL,
         xlim = x.lim,
         ylim = y.lim,
         las = 1, 
         xlab = 'Odds ratio',
         ylab = '',
         cex.lab = 1.5,
         axes = FALSE)
    
    segments(1, y.lim[1], 1, y.lim[2], lty = 3)
    
    for(i in 1:2){
      
      points(exp(medEff[[i]][rows[[1]], panel]),
             y.pos[[i]] + 0.1,
             pch = c(15, 16)[i])
      
      segments(exp(medEff[[i]][rows[[2]], panel]),
               y.pos[[i]] + 0.1,
               exp(medEff[[i]][rows[[3]], panel]),
               y.pos[[i]] + 0.1)
      
      points(exp(medEffInteract[[i]][rows[[1]], panel]),
             y.pos[[i]] - 0.1,
             pch = c(0, 1)[i])
      
      segments(exp(medEffInteract[[i]][rows[[2]], panel]),
               y.pos[[i]] - 0.1,
               exp(medEffInteract[[i]][rows[[3]], panel]),
               y.pos[[i]] - 0.1,
               lty = 2)
      
    }
    
    axis(1,
         at = seq(x.lim[1], x.lim[2], 0.25),
         pos = y.lim[1] - diff(y.lim)*space,
         labels = c(0.5, '', 1, '', 1.5),
         cex.axis = 1.5)
    
    title(main = st$strain[17 + panel],
          line = 2,
          font = 2, cex.main = 1.5)
    
  }
  
}


#------------------------------- Figure S9 ----------------------------------------%


##' 
##' Function to plot of number of strains by titer changes
##' 
##' @param stat dataframe of aty returned from atyData()
##' @param panel, 1 = decreased 
##'               2 = unchanged
##'               3 = increased by 2 fold
##'               4 = seroconversion
##'          
##' @return panel A to D of Fig. s9
##'

histPlot <- function(stat, panel){
  
  y <- stat[ ,panel+1]
  y <- y[y > 0]
  
  n <- length(y)
  x <- seq(1, 21, 1)
  p <- sapply(x, function(i) sum(y==i))/n
  
  cols <- c('#8BAD23', grey(0.8), '#CDB2CD', '#590059')
  x.lim <- c(0, 21)
  # y.lim <- c(0, 
  #            c(0.4, 0.1, 0.1, 0.2)[panel])
  y.lim <- c(0, 0.4)
  margin <- c(6,6,1,1)
  
  par(new = FALSE, mar = margin)
  plot(NULL, xlim = x.lim, ylim = y.lim, 
       xlab = c("Number of decreased strains",
                "Number of unchanged strains",
                "Number of increased strains",
                "Number of 4-fold increased strains")[panel], 
       ylab = "Proportion",
       las = 1, axes = FALSE,
       cex.lab = 1.5)
  
  rect(x-0.9, rep(0, length(x)),
       x-0.1, p, col = cols[panel], border = NA)
  
  axis(1, seq(x.lim[1], x.lim[2], diff(x.lim)), 
       pos = y.lim[1] - diff(y.lim)*0.02,
       tcl = -0.3, labels = NA, xpd = TRUE, 
       cex.axis = 1.5)
  axis(1, seq(x.lim[1], x.lim[2], 5),
       pos = y.lim[1] - diff(y.lim)*0.02, 
       col = NA, col.ticks = "black", 
       tcl = -0.3, cex.axis = 1.5)
  axis(1, seq(x.lim[1], x.lim[2], 1), 
       pos = y.lim[1] - diff(y.lim)*0.02,
       tcl = -0.15, labels = NA, xpd = TRUE,
       col = NA, col.ticks = "black", 
       cex.axis = 1.5)
  
  # axis(2, 
  #      seq(y.lim[1], y.lim[2], c(0.1, 0.05, 0.05, 0.1)[panel]), 
  #      las = 1, pos = x.lim[1], cex.axis = 1.5)
  axis(2, 
       seq(y.lim[1], y.lim[2], 0.1), 
       las = 1, pos = x.lim[1], cex.axis = 1.5)
  
  text(x.lim[1] - diff(x.lim)*0.3,
       y.lim[2],
       LETTERS[panel],
       # letters[panel],
       xpd = TRUE, font = 2, cex = 2)
  
}


##' 
##' Function to plot age distribution by titer changes
##' 
##' @param stat dataframe of aty returned from atyData()
##' @param panel, 1 = decreased 
##'               2 = unchanged
##'               3 = increased by 2 fold
##'               4 = seroconversion
##' @param ci confidence interval from boostraps, returned from boostrap()
##'          
##' @return panel E to F of Fig. s9
##'

ageDistributionPlot <- function(stat, panel, ci){
  
  age <- stat$age
  y <- stat[ ,panel+1] # by change
  y <- stat$age[y>0]
  
  n <- length(y)
  x <- seq(0, 80, 10)
  
  x.used <- sapply(x, function(i) sum(y>=i & y<i+10, na.rm = TRUE))
  x.base <- sapply(x, function(i) sum(age>=i & age<i+10, na.rm = TRUE))
  p <- x.used/n
  p.base <- x.base/length(age)
  
  fisher.p <- fisher.test(cbind(x.used, x.base), simulate.p.value = TRUE)$p.value
  fisher.p <- round(fisher.p, 3)
  
  cols <- c('#8BAD23', grey(0.8), '#CDB2CD', '#590059')
  x.lim <- c(0, 90)
  y.lim <- c(0, 0.4)
  margin <- c(6,6,1,1)
  
  par(new = FALSE, mar = margin)
  plot(NULL, xlim = x.lim, ylim = y.lim, 
       xlab = "Age at sampling, years", 
       ylab = "Proportion",
       las = 1, axes = FALSE,
       cex.lab = 1.5)
  
  rect(x+1, rep(0, length(x)),
       x+9, p, col = cols[panel], border = NA)
  
  lines(c(x, 90), 
        c(p.base, p.base[length(p.base)]), 
        # col = grey(0.5), 
        col = 'black', 
        "s", lwd = 1.5)
  lines(c(x, 90), 
        c(ci[,1], ci[length(p.base), 1]), 
        # col = grey(0.5), 
        col = 'black', 
        "s", lwd = 1.5, lty = 3)
  lines(c(x, 90), 
        c(ci[,2], ci[length(p.base), 2]), 
        # col = grey(0.5), 
        col = 'black', 
        "s", lwd = 1.5, lty = 3)
  
  axis(1, seq(x.lim[1], x.lim[2], 10),
       pos = y.lim[1] - diff(y.lim)*0.02, 
       tcl = -0.3, cex.axis = 1.5)
  
  axis(2, 
       seq(y.lim[1], y.lim[2], 0.1), 
       las = 1, pos = x.lim[1], cex.axis = 1.5)
  
  text(mean(x.lim), y.lim[2],
       paste('P = ', sprintf('%.03f', fisher.p), sep = ''),
       cex = 1.5, xpd = TRUE)
  
  text(x.lim[1] - diff(x.lim)*0.3,
       y.lim[2],
       # LETTERS[panel+4],
       letters[panel+4],
       xpd = TRUE, font = 2, cex = 2)
  
} 



#------------------------------- Figure S10 ----------------------------------------%


##' 
##' Function to plot spline shape of age in logistic regression model for seroconversions
##' 
##' @param out data used to assess seroconversion, returned from getSeroconversionData()
##'          
##' @return figure s10
##'
##'

getSplinePlotForMetric <- function(out, index){
  
  
  base <- 'seroconversion ~ s(age_at_sampling) + titer_outcome_strain_v1 + 
  titer_previous_strain_v1'
  
  predList <- lapply(1:4, function(i){ # strain
    
    metric <- c("AUC", "W40", "W10", "ATY")[index]
    rg <- 'adj'
    
    used <- colnames(out[[i]])[grepl(metric, colnames(out[[i]]))]
    used <- used[!(grepl("residual", used))]
    used <- used[grepl(rg, used)]
    
    fml <- as.formula(paste(base, '+ s(', used, ')', sep = ""))
    m <- gam(fml, data = out[[i]], family = binomial())
    
    new <- data.frame(age_at_sampling = seq(5, 80, 0.2))
    new$titer_outcome_strain_v1 <- mean(out[[i]]$titer_outcome_strain_v1)
    new$titer_previous_strain_v1 <- mean(out[[i]]$titer_previous_strain_v1)
    new$used <- mean(out[[i]][ ,used])
    colnames(new)[4] <- used
    
    pred <- predict(m, newdata = new, se.fit = TRUE)
    
    ci <- cbind(pred$fit - 1.96*pred$se.fit,
                pred$fit + 1.96*pred$se.fit)
    
    rt <- cbind(pred$fit, ci)
    
    return(rt)
    
  })
  
  cols <- c('#0392cf', '#7bc043', '#f37736', '#ee4035')
  
  age <- seq(5, 80, 0.2)
  
  x.lim <- c(0, 80)
  y.lim <- c(0, 6)
  
  par(mar = c(5,5,4,4))
  plot(NULL,
       xlim = x.lim,
       ylim = y.lim,
       las = 1,
       xlab = 'Age at sampling',
       ylab = 'Odds ratio',
       main = c("AUC", "Width 40", "Width 10", "ATY")[index],
       cex.lab = 1.5,
       cex.axis = 1.5,
       cex.main = 1.5)
  
  sapply(1:4, function(i){
    
    # lines(age, exp(predList[[i]][,1]), lty = i, lwd = 2)
    
    est <- exp(predList[[i]][,1])
    
    polygon(c(age, rev(age)),
            c(exp(predList[[i]][,2]),
              rev(exp(predList[[i]][,3]))),
            col = adjustcolor(cols[i], alpha.f = 0.2),
            border = NA)
    lines(age, exp(predList[[i]][,1]), 
          # lty = i, 
          lty = 1,
          lwd = 2, 
          col = adjustcolor(cols[i], alpha.f = 0.8))
    
  })
  
  if(index == 4){
    
    legend(x = 'topright',
           legend = c("A/Perth/2009", 
                      "A/Victoria/2009", 
                      "A/Texas/2012", 
                      "A/HongKong/2014"),
           # lty = 1:4,
           lty = 1,
           bty = 'n',
           cex = 1.5,
           col = cols)
    
  }
  
}



#------------------------------- Figure S12 ----------------------------------------%


##' 
##' Function to plot coefficients from univariable logistic regression 
##' on seroconversion to four recent strains
##' 
##' @param uniEst esitmation returned from lgetDataForUniTiter()
##' @param uniEstAdj esitmation returned from lgetDataForUniTiter() adjusted
##'                  for pre-existing titers to strain i
##' @param strain ith outcome strain
##'               1 = A/Perth/2009
##'               2 = A/Victoria/2009
##'               3 = A/Texas/2012
##'               4 = A/HongKong/2014
##' @param st dataframe contains information on strains
##'          
##' @return figure s12
##'
##'

ORSingleStrainPlot <- function(uniEst, uniEstAdj, strain, st){
  
  usedEst <- uniEst[[strain]]
  usedEstAdj <- uniEstAdj[[strain]]
  
  nTotal <- ncol(usedEst)
  x.1 <- seq(0.3, nTotal, 1)
  x.2 <- seq(0.6, nTotal, 1)
  
  x.lim <- c(0, nTotal)
  y.lim <- c(0, 2)
  y.diff <- 1
  space <- 0.05
  
  par(mar = c(12,5,5,5))
  plot(NULL,
       xlim = x.lim,
       ylim = y.lim,
       las = 1,
       xlab = '',
       ylab = 'Odds ratio',
       axes = FALSE,
       cex.lab = 1.5)
  
  segments(x.lim[1], 1, x.lim[2], 1, lty = 3)
  
  points(x.1, usedEst[1, ], pch = c(1, 2)[1], cex = 1.2)
  segments(x.1, usedEst[2, ], x.1, usedEst[3, ], lty = 1.5)
  segments(x.1-0.1, usedEst[2, ], x.1+0.1, usedEst[2, ], lty = 1.5)
  segments(x.1-0.1, usedEst[3, ], x.1+0.1, usedEst[3, ], lty = 1.5)
  
  points(x.2, usedEstAdj[1, ], pch = c(1, 2)[2], cex = 1.2, col = 'red')
  segments(x.2, usedEstAdj[2, ], x.2, usedEstAdj[3, ], lty = 1.5, col = 'red')
  segments(x.2-0.1, usedEstAdj[2, ], x.2+0.1, usedEstAdj[2, ], lty = 1.5, col = 'red')
  segments(x.2-0.1, usedEstAdj[3, ], x.2+0.1, usedEstAdj[3, ], lty = 1.5, col = 'red')
  
  rect(xleft = x.lim[1] - diff(x.lim)*space,
       ybottom = y.lim[1] - diff(y.lim)*space,
       xright = x.lim[2] + diff(x.lim)*space,
       ytop = y.lim[2] + diff(y.lim)*space,
       xpd = TRUE)
  
  axis(1, 
       seq(x.lim[1], x.lim[2], 1), 
       pos = y.lim[1] - diff(y.lim)*space, 
       cex.axis = 1.5,
       labels = NA, 
       tick = TRUE,
       tcl = -0.4, col.ticks = 1, col = NA)
  text(seq(x.lim[1]+0.5, x.lim[2], 1),
       y.lim[1] - diff(y.lim)*space*1.2,
       st$strain[1:nTotal],
       srt = 90, cex = 1.5, xpd = TRUE, adj = c(1, 0.5))
  
  axis(2, 
       seq(y.lim[1], y.lim[2], y.diff), 
       pos = x.lim[1] - diff(x.lim)*space, 
       cex.axis = 1.5,
       labels = seq(y.lim[1], y.lim[2], y.diff), 
       tick = TRUE, las = 1,
       tcl = -0.4, col.ticks = 1, col = NA)
  
  title(main = st$strain[strain + 17],
        line = 3,
        font = 2,
        xpd = TRUE,
        cex.main = 1.5)
  
  if(strain == 1){
    
    legend(x = 'topleft',
           legend = c('Univariable',
                      'Multivariable'),
           bty = 'n',
           col = c('black', 'red'),
           pch = 1:2,
           adj = c(0, 0),
           cex = 1.5)
    
  }
  
  
}


#------------------------------- Figure S14 ----------------------------------------%


##' 
##' Function to plot strain distribution by titer changes
##' 
##' @param data conversion
##' @param st dataframe contains strain information
##' @param inset TRUE or FALSE
##' @param panel, 1 = all titers
##'               2 = pre-existing titer <= 1:80
##'               3 = pre-existing titer > 1:80
##'          
##' @return figure s14
##'

strainDistrByTiterChangesCeiling <- function(data, st, inset, panel, legend){
  
  used <- list(data[data$fold_of_change == 0, ],
               data[data$fold_of_change < 0, ],
               data[data$fold_of_change > 0 & data$seroconversion ==0, ],
               data[data$seroconversion == 1, ])
  
  out <- data.frame(strain = st$strain,
                    year = st$year_of_isolation,
                    nochange = NA,
                    decrease = NA,
                    increase = NA,
                    fourfold = NA,
                    stringsAsFactors = FALSE)
  
  for(i in 1:nrow(out)){
    tem <- out$strain[i]
    out[i, 3:6] <- unlist(lapply(used, function(x) sum(x$strain == tem, na.rm = TRUE)))
  }
  
  n <- colSums(out[ ,3:6])
  p <- sapply(3:6, function(i) out[,i]/sum(out[,i]))
  ci <- lapply(3:6, function(i) {
    x <- out[,i]
    n <- sum(out[,i])
    t(sapply(x, function(xx) binom.test(xx, n)$conf.int))
  })
  
  # underlying distribution
  base <- data.frame(strain = st$strain,
                     p = NA,
                     stringsAsFactors = FALSE)
  base$p <- sapply(out$strain, function(i) sum(data$strain == i, na.rm = TRUE))
  base$p <- base$p/sum(base$p)
  
  x.max <- 0.2
  if(panel == 3){x.max = 0.5}
  x.space <- x.max/2.5
  x.lim <- c(0-x.space-x.max, x.space + x.max)
  y.lim <- c(0, nrow(st))
  
  cols <- c(grey(0.5), '#8BAD23', '#af84af', '#590059')
  cols.seg <- c(grey(0.5), '#8BAD23', '#af84af', '#590059')
  line.lwd = 1.5
  
  par(mar = c(6,6,6,6))
  plot(NULL,
       xlim = x.lim, ylim = y.lim,
       axes = FALSE,
       xlab = "", ylab = "")
  
  # baseline
  rect(xleft = rep(x.space, y.lim[2]),
       ybottom = seq(y.lim[1] + 0.2, y.lim[2], 1),
       xright = c(base$p[1], base$p)+ x.space,
       ytop = seq(y.lim[1] + 0.8, y.lim[2], 1),
       col = '#f0eae7',
       border = NA)
  
  rect(xleft = rep(x.space*(-1), y.lim[2]),
       ybottom = seq(y.lim[1] + 0.2, y.lim[2], 1),
       xright = c(base$p[1], base$p)*(-1) - x.space,
       ytop = seq(y.lim[1] + 0.8, y.lim[2], 1),
       col = '#f0eae7',
       border = NA)
  
  # unchanged
  segments(ci[[1]][,1]*(-1)-x.space,
           seq(y.lim[1] + 0.65, y.lim[2], 1),
           ci[[1]][,2]*(-1)-x.space,
           seq(y.lim[1] + 0.65, y.lim[2], 1),
           col = cols.seg[1], lwd = line.lwd)
  segments(ci[[1]][,1]*(-1)-x.space,
           seq(y.lim[1] + 0.70, y.lim[2], 1),
           ci[[1]][,1]*(-1)-x.space,
           seq(y.lim[1] + 0.60, y.lim[2], 1),
           col = cols.seg[1], lwd = line.lwd)
  segments(ci[[1]][,2]*(-1)-x.space,
           seq(y.lim[1] + 0.70, y.lim[2], 1),
           ci[[1]][,2]*(-1)-x.space,
           seq(y.lim[1] + 0.60, y.lim[2], 1),
           col = cols.seg[1], lwd = line.lwd)
  points((p[,1]+x.space)*(-1),
         seq(y.lim[1] + 0.65, y.lim[2], 1),
         col = adjustcolor(cols[1], alpha.f = 0.5), 
         xpd = TRUE, pch = 15)
  
  
  # decrease
  segments(ci[[2]][,1]*(-1)-x.space,
           seq(y.lim[1] + 0.35, y.lim[2], 1),
           ci[[2]][,2]*(-1)-x.space,
           seq(y.lim[1] + 0.35, y.lim[2], 1),
           col = cols.seg[2], lwd = line.lwd)
  segments(ci[[2]][,1]*(-1)-x.space,
           seq(y.lim[1] + 0.40, y.lim[2], 1),
           ci[[2]][,1]*(-1)-x.space,
           seq(y.lim[1] + 0.30, y.lim[2], 1),
           col = cols.seg[2], lwd = line.lwd)
  segments(ci[[2]][,2]*(-1)-x.space,
           seq(y.lim[1] + 0.40, y.lim[2], 1),
           ci[[2]][,2]*(-1)-x.space,
           seq(y.lim[1] + 0.30, y.lim[2], 1),
           col = cols.seg[2], lwd = line.lwd)
  points((p[,2]+x.space)*(-1),
         seq(y.lim[1] + 0.35, y.lim[2], 1),
         col = adjustcolor(cols[2], alpha.f = 0.5), 
         xpd = TRUE, pch = 15)
  
  # increase
  segments(ci[[3]][,1] + x.space,
           seq(y.lim[1] + 0.35, y.lim[2], 1),
           ci[[3]][,2] + x.space,
           seq(y.lim[1] + 0.35, y.lim[2], 1),
           col = cols.seg[3], lwd = line.lwd)
  segments(ci[[3]][,1] + x.space,
           seq(y.lim[1] + 0.40, y.lim[2], 1),
           ci[[3]][,1] + x.space,
           seq(y.lim[1] + 0.30, y.lim[2], 1),
           col = cols.seg[3], lwd = line.lwd)
  segments(ci[[3]][,2] + x.space,
           seq(y.lim[1] + 0.40, y.lim[2], 1),
           ci[[3]][,2] + x.space,
           seq(y.lim[1] + 0.30, y.lim[2], 1),
           col = cols.seg[3], lwd = line.lwd)
  points(p[,3]+x.space,
         seq(y.lim[1] + 0.35, y.lim[2], 1),
         col = adjustcolor(cols[3], alpha.f = 0.5), 
         xpd = TRUE, pch = 15)
  
  # four fold increase
  segments(ci[[4]][,1] + x.space,
           seq(y.lim[1] + 0.65, y.lim[2], 1),
           ci[[4]][,2] + x.space,
           seq(y.lim[1] + 0.65, y.lim[2], 1),
           col = cols.seg[4], lwd = line.lwd)
  segments(ci[[4]][,1] + x.space,
           seq(y.lim[1] + 0.70, y.lim[2], 1),
           ci[[4]][,1] + x.space,
           seq(y.lim[1] + 0.60, y.lim[2], 1),
           col = cols.seg[4], lwd = line.lwd)
  segments(ci[[4]][,2] + x.space,
           seq(y.lim[1] + 0.70, y.lim[2], 1),
           ci[[4]][,2] + x.space,
           seq(y.lim[1] + 0.60, y.lim[2], 1),
           col = cols.seg[4], lwd = line.lwd)
  points(p[,4]+x.space,
         seq(y.lim[1] + 0.65, y.lim[2], 1),
         col = adjustcolor(cols[4], alpha.f = 0.5), 
         xpd = TRUE, pch = 15)
  
  
  text(rep(0, y.lim[2]),
       seq(y.lim[1]+0.5, y.lim[2], 1),
       out$strain,
       cex = 1.5)
  
  mtext("Decrease/No change", 3, line = 1, 
        at = (x.space+x.lim[2])/2*(-1),
        font = 2)
  mtext("Increase", 3, line = 1, 
        at = (x.space+x.lim[2])/2,
        font = 2)
  
  mtext("Proportion, %", 1, line = 2, 
        at = (x.space+x.lim[2])/2)
  mtext("Proportion, %", 1, line = 2, 
        at = (x.space+x.lim[2])/2*(-1))
  
  mtext(letters[panel], 
        # LETTERS[panel], 
        side = 2,
        line = 3, 
        at = y.lim[2] + diff(y.lim)*0.08,
        outer = FALSE, 
        las = 1,
        cex = 2,
        font = 2)
  
  
  if(panel < 3){
    axis(1, at = seq(0, x.max, 0.05)+x.space, 
         pos = y.lim[1],
         labels = seq(0, x.max, 0.05)*100,
         cex.axis = 1.5)
    
    axis(1, at = (seq(0, x.max, 0.05)+x.space)*(-1), 
         pos = y.lim[1],
         labels = seq(0, x.max, 0.05)*100,
         cex.axis = 1.5)
  } else {
    axis(1, at = seq(0, x.max, 0.1)+x.space, 
         pos = y.lim[1],
         labels = seq(0, x.max, 0.1)*100,
         cex.axis = 1.5)
    
    axis(1, at = (seq(0, x.max, 0.1)+x.space)*(-1), 
         pos = y.lim[1],
         labels = seq(0, x.max, 0.1)*100,
         cex.axis = 1.5)
  }
  
  axis(2, at = seq(y.lim[1], y.lim[2], 1),
       pos = x.space*0.95, labels = NA)
  
  axis(4, at = seq(y.lim[1], y.lim[2], 1),
       pos = x.space*(-0.95), labels = NA)
  
  
  if(inset == TRUE){
    
    par(new = TRUE)
    x.lim <- c(0, 3)
    y.lim <- c(0, 1)
    x.lim.plot <- c(x.lim[1]-x.lim[2]*6, x.lim[2])
    y.lim.plot <- c(y.lim[1]-y.lim[2]*0.8, y.lim[2]*4)
    
    pp <- n/sum(n)
    
    plot(NULL,
         xlim = x.lim.plot, ylim = y.lim.plot,
         axes = FALSE,
         xlab = "", ylab = "")
    
    rect(seq(x.lim[1]+0.2, x.lim[2], 1),
         rep(y.lim[1], x.lim[2]),
         seq(x.lim[1]+0.8, x.lim[2], 1),
         pp[1:3],
         col = cols[1:3],
         border = NA)
    
    rect(x.lim[1]+2.2,
         y.lim[1],
         x.lim[1]+2.8,
         pp[4],
         col = cols[4],
         border = NA)
    
    axis(1, seq(x.lim[1], x.lim[2], 1), labels = NA, pos = y.lim[1])
    axis(2, seq(y.lim[1], y.lim[2], 0.25), 
         labels = seq(y.lim[1], y.lim[2], 0.25)*100, 
         pos = x.lim[1], las = 1)
    
    text(x.lim[1] - x.lim[2]/2,
         mean(y.lim),
         "Proportion, %", srt = 90)
    
    rect(xleft = x.lim[1] - x.lim[2]/3*2,
         ybottom = -0.2, 
         xright = x.lim[2] + x.lim[2]/5, 
         ytop = 1.5, 
         xpd = TRUE)
    mtext(letters[panel + 3], 
          # LETTERS[panel + 3], 
          side = 2,
          line = -34, 
          at = y.lim[2]+diff(y.lim)*0.2,
          outer = FALSE, 
          las = 1,
          cex = 2,
          font = 2)
    
  }
  
  
  if(legend == TRUE){
    
    legend(x = 0.18,
           y = 15,
           legend = c('Baseline',
                      'Follow-up',
                      'Changes in GMT',
                      '',
                      'Decrease',
                      'No change',
                      '2-fold increase',
                      '4-fold increase',
                      '(Seroconversion)'),
           pch = c(16, 17, 0, rep(15, 6)),
           lty = c(2, 3, rep(1, 7)),
           col = c("#3f93aa",
                   "#fa8072",
                   '#590059',
                   'white',
                   '#8BAD23', 
                   grey(0.5), 
                   '#af84af',
                   '#590059',
                   'white'),
           xpd = TRUE,
           cex = 1.5,
           lwd = 1.5)
    
  }
  
}


#------------------------------- Figure S15 ----------------------------------------%


##' 
##' Function to plot heatmap for strain-specific coefficients from logistic regression 
##' on seroconversion to four recent strains
##' 
##' @param coeff coeff.uni = univariable logistic regression
##'              coeff.mul = multivariable logistic regression
##' @param uni whether or not it is univariable analysis
##'          
##' @return figure s15
##'

coefHeatmap <- function(coeff, uni){
  
  x.lim <- c(0, 21)
  y.lim <- c(0, 21)
  if(uni){
    z.lim <- c(-0.82, 0.223)
    leg <- c(0.5, 0.75, 1, 1.25)
    # txt <- 'A'
    txt <- 'a'
  } else {
    z.lim <- c(-2.1, 0.823)
    leg <- c(0.125, 0.25, 0.5, 1, 2)
    # txt <- 'B'
    txt <- 'b'
  }
  z.diff <- 0.001
  n1 <- length(seq(z.lim[1], 0 - z.diff, z.diff))
  n2 <- length(seq(z.diff, z.lim[2], z.diff))
  
  rows <- seq(1, nrow(coeff), 2)
  
  cols <- c(colorRampPalette(c('#293a4b', '#34495e', '#85919e', '#adb6be', '#eaecee'))(n1),
            '#ffffff',
            colorRampPalette(c('#fdf9e7', '#fcf3cf', '#f9e79f', '#f6db6f', '#f1c40f'))(n2))
  
  p <- coeff[seq(2, 42, 2), ]
  stars <- matrix('**', nrow = nrow(p), ncol = ncol(p))
  stars[p >= 0.01] <- '*'
  stars[p >= 0.05] <- ''
  
  
  par(mar = c(10, 10, 2, 5.5))
  
  plot(NULL,
       xlim = x.lim,
       ylim = y.lim,
       axes = FALSE,
       xlab = '',
       ylab = '')
  
  for(i in 1:x.lim[2]){ # seroconversion strain
    
    for(j in 1:y.lim[2]){ # titer strain
      
      used <- (coeff[rows[i], j] - z.lim[1])*1000 + 1
      
      rect(i - 1,
           j - 1,
           i, j,
           col = cols[used],
           border = NA)
      
      text(i - 0.5, j - 0.5, stars[i, j], cex = 0.5)
      
    }
    
  }
  
  rect(x.lim[1], y.lim[1], x.lim[2], y.lim[2])
  
  axis(1, 
       at = seq(x.lim[1], x.lim[2], 1), 
       pos = y.lim[1],
       labels = NA,
       tcl = -0.4, col.ticks = 1, 
       col=NA, cex.axis = 1.5)
  
  axis(2, 
       at = seq(y.lim[1], y.lim[2], 1), 
       pos = x.lim[1],
       labels = NA,
       tcl = -0.4, col.ticks = 1, 
       col=NA, cex.axis = 1.5)
  
  text(seq(0.5, x.lim[2], 1),
       y.lim[1] - diff(y.lim)*0.05,
       st$strain,
       srt = 90, cex = 1, xpd = TRUE,
       adj= c(1,0.5))
  
  text(x.lim[1] - diff(x.lim)*0.05,
       seq(0.5, y.lim[2], 1),
       st$strain, las = 1,
       cex = 1, xpd = TRUE,
       adj= c(1,0.5))
  
  mtext('Seroconversion', side = 1,
        line = 8, outer = FALSE, cex = 1.5, font = 2)
  mtext('Pre-existing titer', side = 2,
        line = 8, outer = FALSE, cex = 1.5, font = 2)
  
  # legend
  
  y <- seq(y.lim[1], y.lim[2], length.out = length(cols) + 1)
  rect(x.lim[2] + diff(x.lim)*0.06,
       y[1:length(cols)],
       x.lim[2] + diff(x.lim)*0.1,
       y[2:length(y)],
       col = cols,
       border = NA,
       xpd = TRUE)
  
  rect(x.lim[2] + diff(x.lim)*0.06,
       y.lim[1],
       x.lim[2] + diff(x.lim)*0.1,
       y.lim[2],
       xpd = TRUE)
  
  
  y.pos <- round(log(leg), 3)
  y.pos <- (y.pos - z.lim[1])*1000 + 1
  
  segments(x.lim[2] + diff(x.lim)*0.1,
           y[y.pos],
           x.lim[2] + diff(x.lim)*0.115,
           y[y.pos],
           xpd = TRUE)
  
  text(x.lim[2] + diff(x.lim)*0.13,
       y[y.pos],
       leg,
       xpd = TRUE,
       cex = 1,
       adj = c(0, 0.5))
  
  mtext('Odds ratio', side = 4,
        line = 3.5, outer = FALSE, cex = 1, font = 1)
  
  mtext(txt, side = 2,
        at = y.lim[2] + diff(y.lim)*0.1,
        line = 9, outer = FALSE, 
        cex = 2, font = 2, las = 1)
  
}

#------------------------------- Figure S16 & 17 --------------------------------------%


##' 
##' Function to plot predicted probability of seroconversion 
##' vs. observed proportion of seroconversion
##' 
##' @param pred predictions returned from getPred()
##' @param panel ith outcome strain
##'               1 = A/Perth/2009
##'               2 = A/Victoria/2009
##'               3 = A/Texas/2012
##'               4 = A/HongKong/2014
##'          
##' @return figure s16 and s17
##'
##'

predPlot <- function(pred, panel){
  
  y.lim <- x.lim <- c(0.3, 0.8)
  
  par(mar = c(5,5,5,5))
  plot(NULL,
       xlim = x.lim,
       ylim = y.lim,
       # axes = FALSE,
       xlab = 'Predicted probability of seroconversion',
       ylab = 'Observed proportion of seroconversion',
       las = 1,
       cex.lab = 1.5,
       cex.axis = 1.5)
  
  segments(x.lim[1], y.lim[1], x.lim[2], y.lim[2], lty = 2)
  
  points(pred[[panel]][1, ],
         pred[[panel]][4, ],
         pch = 19, col = 'grey')
  
  # range of observation
  segments(pred[[panel]][2, ],
           pred[[panel]][4, ],
           pred[[panel]][3, ],
           pred[[panel]][4, ],
           col = 'grey')
  
  # range of prediction
  segments(pred[[panel]][1, ],
           pred[[panel]][5, ],
           pred[[panel]][1, ],
           pred[[panel]][6, ],
           col = 'grey')
  
  title(main = c('A/Perth/2009',
                 'A/Victoria/2009',
                 'A/Texas/2012',
                 'A/HongKong/2014')[panel])
  
  
}

