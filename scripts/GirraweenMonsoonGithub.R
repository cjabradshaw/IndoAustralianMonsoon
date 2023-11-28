## To what degree are the core records correlated to each other in terms of trends?
## To what extent (within chronological uncertainty) are the Girraween records correlated to
    # (i) the ‘china isotope record
    # (ii) the ‘dole effect’ record, and
    # (iii) the relative rainfall records from LOVECLIM and my original attempt based on
            # insolation and distance to coast (insolation model rainfall sheet)

## Corey Bradshaw
## November 2023

rm(list = ls())

# libraries
library(spatstat)
library(gstat)
library(maps)
library(sp)
library(ape)
library(permute)
library(ggplot2)
library(dplyr)
library(boot)
library(tmvnsim)
library(wCorr)
library(truncnorm)
library(orcutt)
library(lmtest)

# functions
# Hildreth-Lu (named for Clifford Hildreth and John Y. Lu, is a method for adjusting
# a linear model in response to the presence of serial correlation in the error term)
hildreth.lu.func <- function(r, model){
  x <- model.matrix(model)[,-1]
  y <- model.response(model.frame(model))
  n <- length(y)
  t <- 2:n
  y <- y[t]-r*y[t-1]
  x <- x[t]-r*x[t-1]
  
  return(lm(y~x))
}

hildreth.lu.order.func <- function(r, model, order){
  x <- model.matrix(model)[,-1]
  y <- model.response(model.frame(model))
  n <- length(y)
  t <- (order+1):n
  y <- y[t]-r*y[-c((n-(order-1)):n)]
  x <- x[t]-r*x[-c((n-(order-1)):n)]
  
  return(lm(y~x))
}

# geometric mean
gmMean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

## sample avoiding consecutive values of order o
NonSeqSample <- function(x, size, ord=1, replace) {
  # x = vector to sample
  # size = size of resample
  # ord = order to consider (1 = no consecutive, 2 = no 2nd-order consecutive ...)
  # replace = with (TRUE) or without (FALSE) replacement
  vsmp <- sort(sample(x=x, size=size, replace=F))
  diff.vsmp <- diff(vsmp)
  vsmp <- vsmp[-(which(diff.vsmp <= ord)+1)]
  new.size <- size - length(vsmp)
  
  while(new.size >= 1) {
    vsmp.add <- sort(sample(x=x, size=new.size+1, replace=F))
    vsmp <- sort(c(vsmp, vsmp.add))
    diff.vsmp <- diff(vsmp)
    vsmp <- vsmp[-(which(diff.vsmp <= ord)+1)]
    new.size <- size - length(vsmp)
  }
  return(vsmp)
}

##########
## DATA ##
##########

# China speleothem
# combined O isotope record, interpreted as lower values, more rain in China,
# so usually y axis is shown reversed we simply assume that the chronology & values are both perfect, no errors
chspeleo <- read.csv("ChSpeleo2.csv", header=T)
head(chspeleo)
tail(chspeleo)

# dole effect
# modified dole effect based on the isotope composition of atmospheric oxygen,
# high values mean strong tropical hydroclimate, we simply assume that the chronology
# and values are both perfect, no errors
dole <- read.csv("dole2.csv", header=T)
head(dole)

# relative insolation
# we assume that errors on age and insolation are zero
insol <- read.csv("insol.csv", header=T)
head(insol)

# relative rainfall
# based on insolation effect, distance to coast controlling rainfall
# still trying to figure out how to get an error, but age error is assumed zero as it is based on insolation
rain <- read.csv("rainrel.csv", header=T)
head(rain)

# HadCM3 precip anomaly (northern Australia)
prcp.H <- read.csv("NTRegionClimate(0-150ka)_anomalies(Precipitation).csv", header=T)
head(prcp.H)

# HadCM3 precip anomaly (Giraween cell)
prcpG.H <- read.csv("HadCM3_Gironly(0-150ka)_anomalies(Precip).csv", header=T)
colnames(prcpG.H) <- c("age", "prcpG.H")
head(prcpG.H)

# HadCM3 relative precip (Giraween cell); corrected for distance to coast
prcpGWcor.H <- read.csv("HADCMS rel rainfall.csv", header=T)
colnames(prcpGWcor.H)[2] <- "prcpGWcor.H"
head(prcpGWcor.H)

# LOVECLIM precip anomaly (northern Australia)
prcp.L <- read.csv("LOVECLIM_NTRegionClimate(1-150ka)_anomalies(Precipitation).csv", header=T)
head(prcp.L)

# LOVECLIM precip anomaly (South Asia)
prcp.SA <- read.csv("LOVECLIM_SARegionClimate(1-150ka)_anomalies(Precipitation).csv", header=T)
head(prcp.SA)

# LOVECLIM precip anomaly (Giraween cell)
prcpG.L <- read.csv("LOVECLIM_Gironly(1-150ka)_anomalies(Precip).csv", header=T)
colnames(prcpG.L) <- c("age", "prcpG.L")
head(prcpG.L)

# LOVECLIM relative precip (Giraween cell); corrected for distance to coast
prcpGWcor.L <- read.csv("rainLoveClGWdistcor2.csv", header=T)
head(prcpGWcor.L)
tail(prcpGWcor.L)

# % total organic carbon
# with error being 2% of the measured value (interpreted as higher TOC = wetter conditions)
# age uncertainty is 2 standard deviations from rbacon model
toc <- read.csv("toc.csv", header=T)
head(toc)

# hydrogen isotope
# calculated hydrogen isotope composition of precipitation based on n-alkane isotope data,
# errors are deviations from the mean of replicate/triplicate analyses
# interpreted as lower values, more intense monsoon rain
# age uncertainty is 2 standard deviations from rbacon model
Hiso <- read.csv("Hiso.csv", header=T)
head(Hiso)

# % tree
# % of tree pollen in total pollen (so trees:grass)
# with a flat assumption of a 2% counting error (interpretation is more trees = wetter climate),
# while age uncertainty is 2 standard deviations provided by rbacon model
tree <- read.csv("tree.csv", header=T)
tree$treeSD <- sqrt(((tree$tree/100)*(1 - (tree$tree/100)))/tree$pollenN)
head(tree)

## DATA STANDARISATION (MAKE AGE INTERVALS, TIMESPAN & SCALING CONSISTENT)
# standardise to set age interval (min, max, interval)
agest <- seq(0,150,0.5)

# insolation
insol.approx <- approx(insol$age, insol$insol, xout = agest)
insol.sta <- data.frame(insol.approx$x, insol.approx$y)
colnames(insol.sta) <- c("age", "insol")
head(insol.sta)
tail(insol.sta)
plot(insol.sta$age, insol.sta$insol, type="l", xlab="", ylab="insolation")

# rain
rain.approx <- approx(rain$age, rain$rainrel, xout = agest)
rain.sta <- data.frame(rain.approx$x, rain.approx$y)
colnames(rain.sta) <- c("age", "rainrel")
head(rain.sta)
tail(rain.sta)
plot(rain.sta$age, rain.sta$rainrel, type="l", xlab="", ylab="relative rain")

# Giraween LOVECLIM prcp coast-distance-corrected
prcpGcor.L.approx <- approx(prcpGWcor.L$age, prcpGWcor.L$prcpGWcor.L, xout = agest)
prcpGcor.L.sta <- data.frame(prcpGcor.L.approx$x, prcpGcor.L.approx$y)
colnames(prcpGcor.L.sta) <- c("age", "prcpGcor.L")
head(prcpGcor.L.sta)
tail(prcpGcor.L.sta)
plot(prcpGcor.L.sta$age, prcpGcor.L.sta$prcpGcor.L, type="l", xlab="", ylab="dist-corr rel prcp Giraween LOVECLIM")

# LOVECLIM precip anomaly (Giraween cell)
prcpG.L.approx <- approx(prcpG.L$age/1000, prcpG.L$prcpG.L, xout = agest)
prcpG.L.sta <- data.frame(prcpG.L.approx$x, prcpG.L.approx$y)
colnames(prcpG.L.sta) <- c("age", "prcpG.L")
head(prcpG.L.sta)
tail(prcpG.L.sta)
plot(prcpG.L.sta$age, prcpG.L.sta$prcpG.L, type="l", xlab="", ylab="prcp Giraween LOVECLIM")

# Giraween HadCM3 prcp coast-distance-corrected
prcpGcor.H.approx <- approx(prcpGWcor.L$age, prcpGWcor.H$prcpGWcor.H, xout = agest)
prcpGcor.H.sta <- data.frame(prcpGcor.H.approx$x, prcpGcor.H.approx$y)
colnames(prcpGcor.H.sta) <- c("age", "prcpGcor.H")
head(prcpGcor.H.sta)
tail(prcpGcor.H.sta)
plot(prcpGcor.H.sta$age, prcpGcor.H.sta$prcpGcor.H, type="l", xlab="", ylab="dist-corr rel prcp Giraween HadCM3")

# HadCM3 precip anomaly (Giraween cell)
prcpG.H.approx <- approx(prcpG.H$age/1000, prcpG.H$prcpG.H, xout = agest)
prcpG.H.sta <- data.frame(prcpG.H.approx$x, prcpG.H.approx$y)
colnames(prcpG.H.sta) <- c("age", "prcpG.H")
head(prcpG.H.sta)
tail(prcpG.H.sta)
plot(prcpG.H.sta$age, prcpG.H.sta$prcpG.H, type="l", xlab="", ylab="prcp Giraween HadCM3")


noerr.dat <- data.frame(insol.sta$insol, rain.sta$rainrel, prcpGcor.L.sta$prcpGcor.L, prcpG.L.sta$prcpG.L,
                        prcpGcor.H.sta$prcpGcor.H, prcpG.H.sta$prcpG.H)
colnames(noerr.dat) <- c("insol", "rain", "prcpGcor.L", "prcpG.L", "prcpGcor.H", "prcpG.H")
head(noerr.dat)

noerr.cor <- cor(na.omit(noerr.dat), method="spearman")
diag(noerr.cor) <- NA
noerr.cor[upper.tri(noerr.cor)] <- NA
noerr.cor <- noerr.cor[-1,-dim(noerr.cor)[2]]
noerr.cor


# number of iterations
iter <- 1000
itdiv <- iter/10

# how many resamples to avoid temporal autocorrelation?
resample.n <- 50

# remove Heinrich events?
Heinrich.rem <- 0 # 1 = yes, 0 = no

Heinrich0 <- c(12,12.5); Heinrich1 <- c(15.5,18); Heinrich2 <- c(24,26); Heinrich3 <- c(30,32)
Heinrich4 <- c(38,40); Heinrich5 <- c(47,49); Heinrich5a <- c(54,56); Heinrich6 <- c(60,63)
Heinrich7 <- c(66,68); Heinrich8 <- c(85.5,88); Heinrich9 <- c(103,105); Heinrich10 <- c(109,110)
Heinrich11 <- c(129,136); Heinrich12 <- c(139.5,140.5)

Heinrich.dates <- rbind(Heinrich0, Heinrich1, Heinrich2, Heinrich3, Heinrich4, Heinrich5, Heinrich5a,
                        Heinrich6, Heinrich7, Heinrich8, Heinrich9, Heinrich10, Heinrich11, Heinrich12)
colnames(Heinrich.dates) <- c("en","st")
rownames(Heinrich.dates) <- c(paste("H",seq(0,5,1),sep=""), "H5a", paste("H",seq(6,12,1),sep=""))
Heinrich.dates

MaO18insol.slope <- MaO18insol.p <- MaO18insol.R2 <-
  MaO18prcpGcorL.slope <- MaO18prcpGcorL.p <- MaO18prcpGcorL.R2 <-
  MaO18prcpGL.slope <- MaO18prcpGL.p <- MaO18prcpGL.R2 <- 
  MaO18prcpL.slope <- MaO18prcpL.p <- MaO18prcpL.R2 <- 
  MaO18prcpGcorH.slope <- MaO18prcpGcorH.p <- MaO18prcpGcorH.R2 <-
  MaO18prcpGH.slope <- MaO18prcpGH.p <- MaO18prcpGH.R2 <- 
  MaO18prcpH.slope <- MaO18prcpH.p <- MaO18prcpH.R2 <- 
  MaO18dole.slope <- MaO18dole.p <- MaO18dole.R2 <- 
  MaO18tree.slope <- MaO18tree.p <- MaO18tree.R2 <- 
  MaO18Hiso.slope <- MaO18Hiso.p <- MaO18Hiso.R2 <- 
  doleinsol.slope <- doleinsol.p <- doleinsol.R2 <-
  doleprcpGcorL.slope <- doleprcpGcorL.p <- doleprcpGcorL.R2 <-
  doleprcpGL.slope <- doleprcpGL.p <- doleprcpGL.R2 <- 
  doleprcpL.slope <- doleprcpL.p <- doleprcpL.R2 <- 
  doleprcpGcorH.slope <- doleprcpGcorH.p <- doleprcpGcorH.R2 <-
  doleprcpGH.slope <- doleprcpGH.p <- doleprcpGH.R2 <- 
  doleprcpH.slope <- doleprcpH.p <- doleprcpH.R2 <- 
  doletree.slope <- doletree.p <- doletree.R2 <- 
  Hisoinsol.slope <- Hisoinsol.p <- Hisoinsol.R2 <-
  Hisodole.slope <- Hisodole.p <- Hisodole.R2 <-
  HisoprcpGcorL.slope <- HisoprcpGcorL.p <- HisoprcpGcorL.R2 <-
  HisoprcpGL.slope <- HisoprcpGL.p <- HisoprcpGL.R2 <- 
  HisoprcpL.slope <- HisoprcpL.p <- HisoprcpL.R2 <- 
  HisoprcpGcorH.slope <- HisoprcpGcorH.p <- HisoprcpGcorH.R2 <-
  HisoprcpGH.slope <- HisoprcpGH.p <- HisoprcpGH.R2 <- 
  HisoprcpH.slope <- HisoprcpH.p <- HisoprcpH.R2 <- 
  Hisotree.slope <- Hisotree.p <- Hisotree.R2 <- 
  insolprcpGcorL.slope <- insolprcpGcorL.p <- insolprcpGcorL.R2 <-
  insolprcpGL.slope <- insolprcpGL.p <- insolprcpGL.R2 <- 
  insolprcpL.slope <- insolprcpL.p <- insolprcpL.R2 <- 
  insolprcpGcorH.slope <- insolprcpGcorH.p <- insolprcpGcorH.R2 <-
  insolprcpGH.slope <- insolprcpGH.p <- insolprcpGH.R2 <- 
  insolprcpH.slope <- insolprcpH.p <- insolprcpH.R2 <- 
  insoltree.slope <- insoltree.p <- insoltree.R2 <- 
  treeprcpGcorL.slope <- treeprcpGcorL.p <- treeprcpGcorL.R2 <-
  treeprcpGL.slope <- treeprcpGL.p <- treeprcpGL.R2 <- 
  treeprcpL.slope <- treeprcpL.p <- treeprcpL.R2 <- 
  treeprcpGcorH.slope <- treeprcpGcorH.p <- treeprcpGcorH.R2 <-
  treeprcpGH.slope <- treeprcpGH.p <- treeprcpGH.R2 <- 
  treeprcpH.slope <- treeprcpH.p <- treeprcpH.R2 <- 
  prcpLprcpSA.slope <- prcpLprcpSA.p <- prcpLprcpSA.R2 <- rep(NA,iter)

MaO18insol.slope.rs <- MaO18insol.p.rs <- MaO18insol.R2.rs <- MaO18insol.DWp.rs <-
  MaO18prcpGcorL.slope.rs <- MaO18prcpGcorL.p.rs <- MaO18prcpGcorL.R2.rs <- MaO18prcpGcorL.DWp.rs <-
  MaO18prcpGL.slope.rs <- MaO18prcpGL.p.rs <- MaO18prcpGL.R2.rs <- MaO18prcpGL.DWp.rs <- 
  MaO18prcpL.slope.rs <- MaO18prcpL.p.rs <- MaO18prcpL.R2.rs <-  MaO18prcpL.DWp.rs <-
  MaO18prcpGcorH.slope.rs <- MaO18prcpGcorH.p.rs <- MaO18prcpGcorH.R2.rs <- MaO18prcpGcorH.DWp.rs <-
  MaO18prcpGH.slope.rs <- MaO18prcpGH.p.rs <- MaO18prcpGH.R2.rs <- MaO18prcpGH.DWp.rs <- 
  MaO18prcpH.slope.rs <- MaO18prcpH.p.rs <- MaO18prcpH.R2.rs <- MaO18prcpH.DWp.rs <- 
  MaO18dole.slope.rs <- MaO18dole.p.rs <- MaO18dole.R2.rs <- MaO18dole.DWp.rs <- 
  MaO18tree.slope.rs <- MaO18tree.p.rs <- MaO18tree.R2.rs <- MaO18tree.DWp.rs <- 
  MaO18Hiso.slope.rs <- MaO18Hiso.p.rs <- MaO18Hiso.R2.rs <- MaO18Hiso.DWp.rs <- 
  doleinsol.slope.rs <- doleinsol.p.rs <- doleinsol.R2.rs <- doleinsol.DWp.rs <-
  doleprcpGcorL.slope.rs <- doleprcpGcorL.p.rs <- doleprcpGcorL.R2.rs <- doleprcpGcorL.DWp.rs <-
  doleprcpGL.slope.rs <- doleprcpGL.p.rs <- doleprcpGL.R2.rs <- doleprcpGL.DWp.rs <- 
  doleprcpL.slope.rs <- doleprcpL.p.rs <- doleprcpL.R2.rs <- doleprcpL.DWp.rs <- 
  doleprcpGcorH.slope.rs <- doleprcpGcorH.p.rs <- doleprcpGcorH.R2.rs <- doleprcpGcorH.DWp.rs <-
  doleprcpGH.slope.rs <- doleprcpGH.p.rs <- doleprcpGH.R2.rs <- doleprcpGH.DWp.rs <- 
  doleprcpH.slope.rs <- doleprcpH.p.rs <- doleprcpH.R2.rs <- doleprcpH.DWp.rs <- 
  doletree.slope.rs <- doletree.p.rs <- doletree.R2.rs <- doletree.DWp.rs <- 
  Hisoinsol.slope.rs <- Hisoinsol.p.rs <- Hisoinsol.R2.rs <- Hisoinsol.DWp.rs <-
  Hisodole.slope.rs <- Hisodole.p.rs <- Hisodole.R2.rs <- Hisodole.DWp.rs <-
  HisoprcpGcorL.slope.rs <- HisoprcpGcorL.p.rs <- HisoprcpGcorL.R2.rs <- HisoprcpGcorL.DWp.rs <-
  HisoprcpGL.slope.rs <- HisoprcpGL.p.rs <- HisoprcpGL.R2.rs <- HisoprcpGL.DWp.rs <- 
  HisoprcpL.slope.rs <- HisoprcpL.p.rs <- HisoprcpL.R2.rs <- HisoprcpL.DWp.rs <- 
  HisoprcpGcorH.slope.rs <- HisoprcpGcorH.p.rs <- HisoprcpGcorH.R2.rs <- HisoprcpGcorH.DWp.rs <-
  HisoprcpGH.slope.rs <- HisoprcpGH.p.rs <- HisoprcpGH.R2.rs <- HisoprcpGH.DWp.rs <- 
  HisoprcpH.slope.rs <- HisoprcpH.p.rs <- HisoprcpH.R2.rs <- HisoprcpH.DWp.rs <- 
  Hisotree.slope.rs <- Hisotree.p.rs <- Hisotree.R2.rs <- Hisotree.DWp.rs <- 
  insolprcpGcorL.slope.rs <- insolprcpGcorL.p.rs <- insolprcpGcorL.R2.rs <- insolprcpGcorL.DWp.rs <-
  insolprcpGL.slope.rs <- insolprcpGL.p.rs <- insolprcpGL.R2.rs <- insolprcpGL.DWp.rs <- 
  insolprcpL.slope.rs <- insolprcpL.p.rs <- insolprcpL.R2.rs <- insolprcpL.DWp.rs <- 
  insolprcpGcorH.slope.rs <- insolprcpGcorH.p.rs <- insolprcpGcorH.R2.rs <- insolprcpGcorH.DWp.rs <-
  insolprcpGH.slope.rs <- insolprcpGH.p.rs <- insolprcpGH.R2.rs <- insolprcpGH.DWp.rs <- 
  insolprcpH.slope.rs <- insolprcpH.p.rs <- insolprcpH.R2.rs <- insolprcpH.DWp.rs <- 
  insoltree.slope.rs <- insoltree.p.rs <- insoltree.R2.rs <- insoltree.DWp.rs <- 
  treeprcpGcorL.slope.rs <- treeprcpGcorL.p.rs <- treeprcpGcorL.R2.rs <- treeprcpGcorL.DWp.rs <-
  treeprcpGL.slope.rs <- treeprcpGL.p.rs <- treeprcpGL.R2.rs <- treeprcpGL.DWp.rs <- 
  treeprcpL.slope.rs <- treeprcpL.p.rs <- treeprcpL.R2.rs <- treeprcpL.DWp.rs <- 
  treeprcpGcorH.slope.rs <- treeprcpGcorH.p.rs <- treeprcpGcorH.R2.rs <- treeprcpGcorH.DWp.rs <-
  treeprcpGH.slope.rs <- treeprcpGH.p.rs <- treeprcpGH.R2.rs <- treeprcpGH.DWp.rs <- 
  treeprcpH.slope.rs <- treeprcpH.p.rs <- treeprcpH.R2.rs <- treeprcpH.DWp.rs <-
  prcpLprcpSA.slope.rs <- prcpLprcpSA.p.rs <- prcpLprcpSA.R2.rs <- prcpLprcpSA.DWp.rs <- rep(NA,iter)

allcor.arr <- array(data=NA, dim=c(13,13,iter))
alldat.arr <- array(data=NA, dim=c((length(agest)), 13, iter))

LOVECLIM.NA.SAcor.arr <- array(data=NA, dim=c(2,2,iter))
LOVECLIM.NA.SAdat.arr <- array(data=NA, dim=c((length(agest)), 2, iter))

for (i in 1:iter) {
  
  # China speleothem
  age.it <- val.it <-  rep(NA,dim(chspeleo)[1])
  for (t in 1:dim(chspeleo)[1]) {
    age.it[t] <- rtruncnorm(1, a=0, b=Inf, chspeleo$age[t], chspeleo$ageSD[t])
    val.it[t] <- rnorm(1, chspeleo$MaO18[t], chspeleo$MaO18SD[t])
  } # end t
  chspeleo.it <- data.frame(age.it, val.it)
  
  chspeleo.it.approx <- approx(chspeleo.it$age.it, chspeleo.it$val.it, xout = agest)
  chspeleo.it.sta <- data.frame(chspeleo.it.approx$x, chspeleo.it.approx$y)
  colnames(chspeleo.it.sta) <- c("age", "MaO18")
  
  # dole
  age.it <- val.it <-  rep(NA,dim(dole)[1])
  for (t in 1:dim(dole)[1]) {
    age.it[t] <- rtruncnorm(1, a=0, b=Inf, dole$age[t], dole$ageSD[t])
    val.it[t] <- rnorm(1, dole$dDE[t], dole$dDESD[t])
  } # end t
  dole.it <- data.frame(age.it, val.it)
  
  dole.it.approx <- approx(dole.it$age.it, dole.it$val.it, xout = agest)
  dole.it.sta <- data.frame(dole.it.approx$x, dole.it.approx$y)
  colnames(dole.it.sta) <- c("age", "dDE")
  
  # HadCM3 precipitation anomaly northern Australia
  val.it <-  rep(NA,dim(prcp.H)[1])
  for (t in 1:dim(prcp.H)[1]) {
    val.it[t] <- runif(1, min=prcp.H$LowCI[t], max=prcp.H$UpCI[t])
  } # end t
  prcpH.it <- data.frame(prcp.H$Time, val.it)
  colnames(prcpH.it)[1] <- "age"
  
  prcpH.it.approx <- approx(prcpH.it$age/1000, prcpH.it$val.it, xout = agest)
  prcpH.it.sta <- data.frame(prcpH.it.approx$x, prcpH.it.approx$y)
  colnames(prcpH.it.sta) <- c("age", "prcpH")
  
  # LOVECLIM precipitation anomaly (northern Australia)
  val.it <-  rep(NA,dim(prcp.L)[1])
  for (t in 1:dim(prcp.L)[1]) {
    val.it[t] <- runif(1, min=prcp.L$LowCI[t], max=prcp.L$UpCI[t])
  } # end t
  prcpL.it <- data.frame(prcp.L$Time, val.it)
  colnames(prcpL.it)[1] <- "age"
  
  prcpL.it.approx <- approx(prcpL.it$age/1000, prcpL.it$val.it, xout = agest)
  prcpL.it.sta <- data.frame(prcpL.it.approx$x, prcpL.it.approx$y)
  colnames(prcpL.it.sta) <- c("age", "prcpL")
  
  # LOVECLIM precipitation anomaly (South Asia)
  val.it <-  rep(NA,dim(prcp.SA)[1])
  for (t in 1:dim(prcp.SA)[1]) {
    val.it[t] <- runif(1, min=prcp.SA$LowCI[t], max=prcp.SA$UpCI[t])
  } # end t
  prcpSA.it <- data.frame(prcp.SA$Time, val.it)
  colnames(prcpSA.it)[1] <- "age"
  
  prcpSA.it.approx <- approx(prcpSA.it$age/1000, prcpSA.it$val.it, xout = agest)
  prcpSA.it.sta <- data.frame(prcpSA.it.approx$x, prcpSA.it.approx$y)
  colnames(prcpSA.it.sta) <- c("age", "prcpSA")
  
  # toc
  age.it <- val.it <-  rep(NA,dim(toc)[1])
  for (t in 1:dim(toc)[1]) {
    age.it[t] <- rtruncnorm(1, a=0, b=Inf, toc$age[t], toc$ageSD[t]/2)
    val.it[t] <- rnorm(1, toc$toc[t], toc$tocSD[t])
  } # end t
  toc.it <- data.frame(age.it, val.it)
  
  toc.it.approx <- approx(toc.it$age.it, toc.it$val.it, xout = agest)
  toc.it.sta <- data.frame(toc.it.approx$x, toc.it.approx$y)
  colnames(toc.it.sta) <- c("age", "toc")
  
  # H isotope
  age.it <- val.it <-  rep(NA,dim(Hiso)[1])
  for (t in 1:dim(Hiso)[1]) {
    age.it[t] <- rtruncnorm(1, a=0, b=Inf, Hiso$age[t], Hiso$ageSD[t]/2)
    val.it[t] <- rnorm(1, Hiso$DHprecip[t], Hiso$DHprecipSD[t])
  } # end t
  Hiso.it <- data.frame(age.it, val.it)
  
  Hiso.it.approx <- approx(Hiso.it$age.it, Hiso.it$val.it, xout = agest)
  Hiso.it.sta <- data.frame(Hiso.it.approx$x, Hiso.it.approx$y)
  colnames(Hiso.it.sta) <- c("age", "DHprecip")
  
  # % tree
  age.it <- val.it <-  rep(NA,dim(tree)[1])
  for (t in 1:dim(tree)[1]) {
    age.it[t] <- rtruncnorm(1, a=0, b=Inf, tree$age[t], tree$ageSD[t]/2)
    val.it[t] <- rnorm(1, tree$tree[t], tree$treeSD[t])
  } # end t
  tree.it <- data.frame(age.it, val.it)
  
  tree.it.approx <- approx(tree.it$age.it, tree.it$val.it, xout = agest)
  tree.it.sta <- data.frame(tree.it.approx$x, tree.it.approx$y)
  colnames(tree.it.sta) <- c("age", "tree")
  
  # full correlation matrix
  all.dat <- (cbind(noerr.dat, chspeleo.it.sta$MaO18, dole.it.sta$dDE, prcpL.it.sta$prcpL, prcpH.it.sta$prcpH,
                    toc.it.sta$toc, Hiso.it.sta$DHprecip, tree.it.sta$tree))
  colnames(all.dat)[7:13] <- c("MaO18","dDE","prcpL","prcpH","toc","DHprecip","tree")
  
  LOVECLIM.NA.SA.dat <- cbind(prcpL.it.sta$prcpL, prcpSA.it.sta$prcpSA)
  colnames(LOVECLIM.NA.SA.dat) <- c("prcpL", "prcpSA")
  
  if (Heinrich.rem == 1) {
    for(h in 1:dim(Heinrich.dates)[1]) {
      all.dat[((Heinrich.dates[h,1]*2)+1):((Heinrich.dates[h,2]*2)+1), ] <- NA
      LOVECLIM.NA.SA.dat[((LOVECLIM.NA.SA.dat[h,1]*2)+1):((LOVECLIM.NA.SA.dat[h,2]*2)+1), ] <- NA
    } # end h
  } # end if
  
  alldat.arr[,,i] <- as.matrix(all.dat)
  LOVECLIM.NA.SAdat.arr[,,i] <- as.matrix(LOVECLIM.NA.SA.dat)
  
  if (Heinrich.rem == 1) {
    for(h in 1:dim(Heinrich.dates)[1]) {
      alldat.arr[((Heinrich.dates[h,1]*2)+1):((Heinrich.dates[h,2]*2)+1), , i] <- NA
      LOVECLIM.NA.SAdat.arr[((Heinrich.dates[h,1]*2)+1):((Heinrich.dates[h,2]*2)+1), , i] <- NA
    } # end h
  } # end if
  
  # correlation matrix
  allcor.arr[,,i] <- (cor(na.omit(all.dat), method="spearman"))
  LOVECLIM.NA.SAcor.arr[,,i] <- (cor(na.omit(LOVECLIM.NA.SA.dat), method="spearman"))
  
  # scale
  all.dat.sc <- scale(all.dat, scale=T, center=T)
  all.dat.sc.df <- as.data.frame(all.dat.sc)
  
  LOVECLIM.NA.SA.dat.sc <- scale(LOVECLIM.NA.SA.dat, scale=T, center=T)
  LOVECLIM.NA.SA.dat.sc.df <- as.data.frame(LOVECLIM.NA.SA.dat.sc)
  
  # linear fits
  # chspeleo vs. insol
  MaO18insol <- lm(all.dat.sc.df$MaO18 ~ all.dat.sc.df$insol)
  MaO18insol.orc <- cochrane.orcutt(MaO18insol, convergence = 5, max.iter=1000)
  MaO18insol.hl <- hildreth.lu.order.func(MaO18insol.orc$rho, MaO18insol, order=1)
  MaO18insol.slope[i] <- as.numeric(MaO18insol.hl$coefficients[2])
  MaO18insol.summ <- summary(MaO18insol.hl)
  MaO18insol.p[i] <- MaO18insol.summ$coefficients[8]
  MaO18insol.R2[i] <- MaO18insol.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    length(resample.ind)
    MaO18insol.rs <- lm(all.dat.sc.df$MaO18[resample.ind] ~ all.dat.sc.df$insol[resample.ind])
    MaO18insol.slope.rs[i] <- as.numeric(MaO18insol.rs$coefficients[2])
    MaO18insol.summ.rs <- summary(MaO18insol.rs)
    MaO18insol.p.rs[i] <- MaO18insol.summ.rs$coefficients[8]
    MaO18insol.R2.rs[i] <- MaO18insol.summ.rs$adj.r.squared
    MaO18insol.DWp.rs[i] <- as.numeric(dwtest(MaO18insol.rs)[4])
  
  #  chspeleo vs. prcpGcor.L
  MaO18prcpGcorL <- lm(all.dat.sc.df$MaO18 ~ all.dat.sc.df$prcpGcor.L)
  MaO18prcpGcorL.orc <- cochrane.orcutt(MaO18prcpGcorL, convergence = 5, max.iter=1000)
  MaO18prcpGcorL.hl <- hildreth.lu.order.func(MaO18prcpGcorL.orc$rho, MaO18prcpGcorL, order=1)
  MaO18prcpGcorL.slope[i] <- as.numeric(MaO18prcpGcorL.hl$coefficients[2])
  MaO18prcpGcorL.summ <- summary(MaO18prcpGcorL.hl)
  MaO18prcpGcorL.p[i] <- MaO18prcpGcorL.summ$coefficients[8]
  MaO18prcpGcorL.R2[i] <- MaO18prcpGcorL.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    MaO18prcpGcorL.rs <- lm(all.dat.sc.df$MaO18[resample.ind] ~ all.dat.sc.df$prcpGcor.L[resample.ind])
    MaO18prcpGcorL.slope.rs[i] <- as.numeric(MaO18prcpGcorL.rs$coefficients[2])
    MaO18prcpGcorL.summ.rs <- summary(MaO18prcpGcorL.rs)
    MaO18prcpGcorL.p.rs[i] <- MaO18prcpGcorL.summ.rs$coefficients[8]
    MaO18prcpGcorL.R2.rs[i] <- MaO18prcpGcorL.summ.rs$adj.r.squared
    MaO18prcpGcorL.DWp.rs[i] <- as.numeric(dwtest(MaO18prcpGcorL.rs)[4])
    
  # chspeleo vs. prcpG.L
  MaO18prcpGL <- lm(all.dat.sc.df$MaO18 ~ all.dat.sc.df$prcpG.L)
  MaO18prcpGL.orc <- cochrane.orcutt(MaO18prcpGL, convergence = 5, max.iter=1000)
  MaO18prcpGL.hl <- hildreth.lu.order.func(MaO18prcpGL.orc$rho, MaO18prcpGL, order=1)
  MaO18prcpGL.slope[i] <- as.numeric(MaO18prcpGL.hl$coefficients[2])
  MaO18prcpGL.summ <- summary(MaO18prcpGL.hl)
  MaO18prcpGL.p[i] <- MaO18prcpGL.summ$coefficients[8]
  MaO18prcpGL.R2[i] <- MaO18prcpGL.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    MaO18prcpGL.rs <- lm(all.dat.sc.df$MaO18[resample.ind] ~ all.dat.sc.df$prcpG.L[resample.ind])
    MaO18prcpGL.slope.rs[i] <- as.numeric(MaO18prcpGL.rs$coefficients[2])
    MaO18prcpGL.summ.rs <- summary(MaO18prcpGL.rs)
    MaO18prcpGL.p.rs[i] <- MaO18prcpGL.summ.rs$coefficients[8]
    MaO18prcpGL.R2.rs[i] <- MaO18prcpGL.summ.rs$adj.r.squared
    MaO18prcpGL.DWp.rs[i] <- as.numeric(dwtest(MaO18prcpGL.rs)[4])
  
  # chspeleo vs. prcp.L
  MaO18prcpL <- lm(all.dat.sc.df$MaO18 ~ all.dat.sc.df$prcpL)
  MaO18prcpL.orc <- cochrane.orcutt(MaO18prcpL, convergence = 5, max.iter=1000)
  MaO18prcpL.hl <- hildreth.lu.order.func(MaO18prcpL.orc$rho, MaO18prcpL, order=1)
  MaO18prcpL.slope[i] <- as.numeric(MaO18prcpL.hl$coefficients[2])
  MaO18prcpL.summ <- summary(MaO18prcpL.hl)
  MaO18prcpL.p[i] <- MaO18prcpL.summ$coefficients[8]
  MaO18prcpL.R2[i] <- MaO18prcpL.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    MaO18prcpL.rs <- lm(all.dat.sc.df$MaO18[resample.ind] ~ all.dat.sc.df$prcpL[resample.ind])
    MaO18prcpL.slope.rs[i] <- as.numeric(MaO18prcpL.rs$coefficients[2])
    MaO18prcpL.summ.rs <- summary(MaO18prcpL.rs)
    MaO18prcpL.p.rs[i] <- MaO18prcpL.summ.rs$coefficients[8]
    MaO18prcpL.R2.rs[i] <- MaO18prcpL.summ.rs$adj.r.squared
    MaO18prcpL.DWp.rs[i] <- as.numeric(dwtest(MaO18prcpL.rs)[4])
  
  #  chspeleo vs. prcpGcor.H
  MaO18prcpGcorH <- lm(all.dat.sc.df$MaO18 ~ all.dat.sc.df$prcpGcor.H)
  MaO18prcpGcorH.orc <- cochrane.orcutt(MaO18prcpGcorH, convergence = 5, max.iter=1000)
  MaO18prcpGcorH.hl <- hildreth.lu.order.func(MaO18prcpGcorH.orc$rho, MaO18prcpGcorH, order=1)
  MaO18prcpGcorH.slope[i] <- as.numeric(MaO18prcpGcorH.hl$coefficients[2])
  MaO18prcpGcorH.summ <- summary(MaO18prcpGcorH.hl)
  MaO18prcpGcorH.p[i] <- MaO18prcpGcorH.summ$coefficients[8]
  MaO18prcpGcorH.R2[i] <- MaO18prcpGcorH.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    MaO18prcpGcorH.rs <- lm(all.dat.sc.df$MaO18[resample.ind] ~ all.dat.sc.df$prcpGcor.H[resample.ind])
    MaO18prcpGcorH.slope.rs[i] <- as.numeric(MaO18prcpGcorH.rs$coefficients[2])
    MaO18prcpGcorH.summ.rs <- summary(MaO18prcpGcorH.rs)
    MaO18prcpGcorH.p.rs[i] <- MaO18prcpGcorH.summ.rs$coefficients[8]
    MaO18prcpGcorH.R2.rs[i] <- MaO18prcpGcorH.summ.rs$adj.r.squared
    MaO18prcpGcorH.DWp.rs[i] <- as.numeric(dwtest(MaO18prcpGcorH.rs)[4])
    
  # chspeleo vs. prcpG.H
  MaO18prcpGH <- lm(all.dat.sc.df$MaO18 ~ all.dat.sc.df$prcpG.H)
  MaO18prcpGH.orc <- cochrane.orcutt(MaO18prcpGH, convergence = 5, max.iter=1000)
  MaO18prcpGH.hl <- hildreth.lu.order.func(MaO18prcpGH.orc$rho, MaO18prcpGH, order=1)
  MaO18prcpGH.slope[i] <- as.numeric(MaO18prcpGH.hl$coefficients[2])
  MaO18prcpGH.summ <- summary(MaO18prcpGH.hl)
  MaO18prcpGH.p[i] <- MaO18prcpGH.summ$coefficients[8]
  MaO18prcpGH.R2[i] <- MaO18prcpGH.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    MaO18prcpGH.rs <- lm(all.dat.sc.df$MaO18[resample.ind] ~ all.dat.sc.df$prcpG.H[resample.ind])
    MaO18prcpGH.slope.rs[i] <- as.numeric(MaO18prcpGH.rs$coefficients[2])
    MaO18prcpGH.summ.rs <- summary(MaO18prcpGH.rs)
    MaO18prcpGH.p.rs[i] <- MaO18prcpGH.summ.rs$coefficients[8]
    MaO18prcpGH.R2.rs[i] <- MaO18prcpGH.summ.rs$adj.r.squared
    MaO18prcpGH.DWp.rs[i] <- as.numeric(dwtest(MaO18prcpGH.rs)[4])
  
  # chspeleo vs. prcp.H
  MaO18prcpH <- lm(all.dat.sc.df$MaO18 ~ all.dat.sc.df$prcpH)
  MaO18prcpH.orc <- cochrane.orcutt(MaO18prcpH, convergence = 5, max.iter=1000)
  MaO18prcpH.hl <- hildreth.lu.order.func(MaO18prcpH.orc$rho, MaO18prcpH, order=1)
  MaO18prcpH.slope[i] <- as.numeric(MaO18prcpH.hl$coefficients[2])
  MaO18prcpH.summ <- summary(MaO18prcpH.hl)
  MaO18prcpH.p[i] <- MaO18prcpH.summ$coefficients[8]
  MaO18prcpH.R2[i] <- MaO18prcpH.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    MaO18prcpH.rs <- lm(all.dat.sc.df$MaO18[resample.ind] ~ all.dat.sc.df$prcpH[resample.ind])
    MaO18prcpH.slope.rs[i] <- as.numeric(MaO18prcpH.rs$coefficients[2])
    MaO18prcpH.summ.rs <- summary(MaO18prcpH.rs)
    MaO18prcpH.p.rs[i] <- MaO18prcpH.summ.rs$coefficients[8]
    MaO18prcpH.R2.rs[i] <- MaO18prcpH.summ.rs$adj.r.squared
    MaO18prcpH.DWp.rs[i] <- as.numeric(dwtest(MaO18prcpH.rs)[4])
  
  # chspeleo vs. dole
  MaO18dole <- lm(all.dat.sc.df$MaO18 ~ all.dat.sc.df$dDE)
  MaO18dole.orc <- cochrane.orcutt(MaO18dole, convergence = 5, max.iter=1000)
  MaO18dole.hl <- hildreth.lu.order.func(MaO18dole.orc$rho, MaO18dole, order=1)
  MaO18dole.slope[i] <- as.numeric(MaO18dole.hl$coefficients[2])
  MaO18dole.summ <- summary(MaO18dole.hl)
  MaO18dole.p[i] <- MaO18dole.summ$coefficients[8]
  MaO18dole.R2[i] <- MaO18dole.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    MaO18dole.rs <- lm(all.dat.sc.df$MaO18[resample.ind] ~ all.dat.sc.df$dDE[resample.ind])
    MaO18dole.slope.rs[i] <- as.numeric(MaO18dole.rs$coefficients[2])
    MaO18dole.summ.rs <- summary(MaO18dole.rs)
    MaO18dole.p.rs[i] <- MaO18dole.summ.rs$coefficients[8]
    MaO18dole.R2.rs[i] <- MaO18dole.summ.rs$adj.r.squared
    MaO18dole.DWp.rs[i] <- as.numeric(dwtest(MaO18dole.rs)[4])
  
  # chspeleo vs. tree
  MaO18tree <- lm(all.dat.sc.df$MaO18 ~ all.dat.sc.df$tree)
  MaO18tree.orc <- cochrane.orcutt(MaO18tree, convergence = 5, max.iter=1000)
  MaO18tree.hl <- hildreth.lu.order.func(MaO18tree.orc$rho, MaO18tree, order=1)
  MaO18tree.slope[i] <- as.numeric(MaO18tree.hl$coefficients[2])
  MaO18tree.summ <- summary(MaO18tree.hl)
  MaO18tree.p[i] <- MaO18tree.summ$coefficients[8]
  MaO18tree.R2[i] <- MaO18tree.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    MaO18tree.rs <- lm(all.dat.sc.df$MaO18[resample.ind] ~ all.dat.sc.df$tree[resample.ind])
    MaO18tree.slope.rs[i] <- as.numeric(MaO18tree.rs$coefficients[2])
    MaO18tree.summ.rs <- summary(MaO18tree.rs)
    MaO18tree.p.rs[i] <- MaO18tree.summ.rs$coefficients[8]
    MaO18tree.R2.rs[i] <- MaO18tree.summ.rs$adj.r.squared
    MaO18tree.DWp.rs[i] <- as.numeric(dwtest(MaO18tree.rs)[4])
  
  # chspeleo vs. Hiso
  MaO18Hiso <- lm(all.dat.sc.df$MaO18 ~ all.dat.sc.df$DHprecip)
  MaO18Hiso.orc <- cochrane.orcutt(MaO18Hiso, convergence = 5, max.iter=1000)
  MaO18Hiso.hl <- hildreth.lu.order.func(MaO18Hiso.orc$rho, MaO18Hiso, order=1)
  MaO18Hiso.slope[i] <- as.numeric(MaO18Hiso.hl$coefficients[2])
  MaO18Hiso.summ <- summary(MaO18Hiso.hl)
  MaO18Hiso.p[i] <- MaO18Hiso.summ$coefficients[8]
  MaO18Hiso.R2[i] <- MaO18Hiso.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    MaO18Hiso.rs <- lm(all.dat.sc.df$MaO18[resample.ind] ~ all.dat.sc.df$DHprecip[resample.ind])
    MaO18Hiso.slope.rs[i] <- as.numeric(MaO18Hiso.rs$coefficients[2])
    MaO18Hiso.summ.rs <- summary(MaO18Hiso.rs)
    MaO18Hiso.p.rs[i] <- MaO18Hiso.summ.rs$coefficients[8]
    MaO18Hiso.R2.rs[i] <- MaO18Hiso.summ.rs$adj.r.squared
    MaO18Hiso.DWp.rs[i] <- as.numeric(dwtest(MaO18Hiso.rs)[4])
  
  # dole vs. insol
  doleinsol <- lm(all.dat.sc.df$dDE ~ all.dat.sc.df$insol)
  doleinsol.orc <- cochrane.orcutt(doleinsol, convergence = 5, max.iter=1000)
  doleinsol.hl <- hildreth.lu.order.func(doleinsol.orc$rho, doleinsol, order=1)
  doleinsol.slope[i] <- as.numeric(doleinsol.hl$coefficients[2])
  doleinsol.summ <- summary(doleinsol.hl)
  doleinsol.p[i] <- doleinsol.summ$coefficients[8]
  doleinsol.R2[i] <- doleinsol.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    doleinsol.rs <- lm(all.dat.sc.df$dDE[resample.ind] ~ all.dat.sc.df$insol[resample.ind])
    doleinsol.slope.rs[i] <- as.numeric(doleinsol.rs$coefficients[2])
    doleinsol.summ.rs <- summary(doleinsol.rs)
    doleinsol.p.rs[i] <- doleinsol.summ.rs$coefficients[8]
    doleinsol.R2.rs[i] <- doleinsol.summ.rs$adj.r.squared
    doleinsol.DWp.rs[i] <- as.numeric(dwtest(doleinsol.rs)[4])
  
  #  dole vs. prcpGcor.L
  doleprcpGcorL <- lm(all.dat.sc.df$dDE ~ all.dat.sc.df$prcpGcor.L)
  doleprcpGcorL.orc <- cochrane.orcutt(doleprcpGcorL, convergence = 5, max.iter=1000)
  doleprcpGcorL.hl <- hildreth.lu.order.func(doleprcpGcorL.orc$rho, doleprcpGcorL, order=1)
  doleprcpGcorL.slope[i] <- as.numeric(doleprcpGcorL.hl$coefficients[2])
  doleprcpGcorL.summ <- summary(doleprcpGcorL.hl)
  doleprcpGcorL.p[i] <- doleprcpGcorL.summ$coefficients[8]
  doleprcpGcorL.R2[i] <- doleprcpGcorL.summ$adj.r.squared
  
  # resampled
  resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
  doleprcpGcorL.rs <- lm(all.dat.sc.df$dDE[resample.ind] ~ all.dat.sc.df$prcpGcor.L[resample.ind])
  doleprcpGcorL.slope.rs[i] <- as.numeric(doleprcpGcorL.rs$coefficients[2])
  doleprcpGcorL.summ.rs <- summary(doleprcpGcorL.rs)
  doleprcpGcorL.p.rs[i] <- doleprcpGcorL.summ.rs$coefficients[8]
  doleprcpGcorL.R2.rs[i] <- doleprcpGcorL.summ.rs$adj.r.squared
  doleprcpGcorL.DWp.rs[i] <- as.numeric(dwtest(doleprcpGcorL.rs)[4])
  
  # dole vs. prcpG.L
  doleprcpGL <- lm(all.dat.sc.df$dDE ~ all.dat.sc.df$prcpG.L)
  doleprcpGL.orc <- cochrane.orcutt(doleprcpGL, convergence = 5, max.iter=1000)
  doleprcpGL.hl <- hildreth.lu.order.func(doleprcpGL.orc$rho, doleprcpGL, order=1)
  doleprcpGL.slope[i] <- as.numeric(doleprcpGL.hl$coefficients[2])
  doleprcpGL.summ <- summary(doleprcpGL.hl)
  doleprcpGL.p[i] <- doleprcpGL.summ$coefficients[8]
  doleprcpGL.R2[i] <- doleprcpGL.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    doleprcpGL.rs <- lm(all.dat.sc.df$dDE[resample.ind] ~ all.dat.sc.df$prcpG.L[resample.ind])
    doleprcpGL.slope.rs[i] <- as.numeric(doleprcpGL.rs$coefficients[2])
    doleprcpGL.summ.rs <- summary(doleprcpGL.rs)
    doleprcpGL.p.rs[i] <- doleprcpGL.summ.rs$coefficients[8]
    doleprcpGL.R2.rs[i] <- doleprcpGL.summ.rs$adj.r.squared
    doleprcpGL.DWp.rs[i] <- as.numeric(dwtest(doleprcpGL.rs)[4])
  
  # dole vs. prcp.L
  doleprcpL <- lm(all.dat.sc.df$dDE ~ all.dat.sc.df$prcpL)
  doleprcpL.orc <- cochrane.orcutt(doleprcpL, convergence = 5, max.iter=1000)
  doleprcpL.hl <- hildreth.lu.order.func(doleprcpL.orc$rho, doleprcpL, order=1)
  doleprcpL.slope[i] <- as.numeric(doleprcpL.hl$coefficients[2])
  doleprcpL.summ <- summary(doleprcpL.hl)
  doleprcpL.p[i] <- doleprcpL.summ$coefficients[8]
  doleprcpL.R2[i] <- doleprcpL.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    doleprcpL.rs <- lm(all.dat.sc.df$dDE[resample.ind] ~ all.dat.sc.df$prcpL[resample.ind])
    doleprcpL.slope.rs[i] <- as.numeric(doleprcpL.rs$coefficients[2])
    doleprcpL.summ.rs <- summary(doleprcpL.rs)
    doleprcpL.p.rs[i] <- doleprcpL.summ.rs$coefficients[8]
    doleprcpL.R2.rs[i] <- doleprcpL.summ.rs$adj.r.squared
    doleprcpL.DWp.rs[i] <- as.numeric(dwtest(doleprcpL.rs)[4])
  
  #  dole vs. prcpGcor.H
  doleprcpGcorH <- lm(all.dat.sc.df$dDE ~ all.dat.sc.df$prcpGcor.H)
  doleprcpGcorH.orc <- cochrane.orcutt(doleprcpGcorH, convergence = 5, max.iter=1000)
  doleprcpGcorH.hl <- hildreth.lu.order.func(doleprcpGcorH.orc$rho, doleprcpGcorH, order=1)
  doleprcpGcorH.slope[i] <- as.numeric(doleprcpGcorH.hl$coefficients[2])
  doleprcpGcorH.summ <- summary(doleprcpGcorH.hl)
  doleprcpGcorH.p[i] <- doleprcpGcorH.summ$coefficients[8]
  doleprcpGcorH.R2[i] <- doleprcpGcorH.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    doleprcpGcorH.rs <- lm(all.dat.sc.df$dDE[resample.ind] ~ all.dat.sc.df$prcpGcor.H[resample.ind])
    doleprcpGcorH.slope.rs[i] <- as.numeric(doleprcpGcorH.rs$coefficients[2])
    doleprcpGcorH.summ.rs <- summary(doleprcpGcorH.rs)
    doleprcpGcorH.p.rs[i] <- doleprcpGcorH.summ.rs$coefficients[8]
    doleprcpGcorH.R2.rs[i] <- doleprcpGcorH.summ.rs$adj.r.squared
    doleprcpGcorH.DWp.rs[i] <- as.numeric(dwtest(doleprcpGcorH.rs)[4])
  
  # dole vs. prcpG.H
  doleprcpGH <- lm(all.dat.sc.df$dDE ~ all.dat.sc.df$prcpG.H)
  doleprcpGH.orc <- cochrane.orcutt(doleprcpGH, convergence = 5, max.iter=1000)
  doleprcpGH.hl <- hildreth.lu.order.func(doleprcpGH.orc$rho, doleprcpGH, order=1)
  doleprcpGH.slope[i] <- as.numeric(doleprcpGH.hl$coefficients[2])
  doleprcpGH.summ <- summary(doleprcpGH.hl)
  doleprcpGH.p[i] <- doleprcpGH.summ$coefficients[8]
  doleprcpGH.R2[i] <- doleprcpGH.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    doleprcpGH.rs <- lm(all.dat.sc.df$dDE[resample.ind] ~ all.dat.sc.df$prcpG.H[resample.ind])
    doleprcpGH.slope.rs[i] <- as.numeric(doleprcpGH.rs$coefficients[2])
    doleprcpGH.summ.rs <- summary(doleprcpGH.rs)
    doleprcpGH.p.rs[i] <- doleprcpGH.summ.rs$coefficients[8]
    doleprcpGH.R2.rs[i] <- doleprcpGH.summ.rs$adj.r.squared
    doleprcpGH.DWp.rs[i] <- as.numeric(dwtest(doleprcpGH.rs)[4])
  
  # dole vs. prcp.H
  doleprcpH <- lm(all.dat.sc.df$dDE ~ all.dat.sc.df$prcpH)
  doleprcpH.orc <- cochrane.orcutt(doleprcpH, convergence = 5, max.iter=1000)
  doleprcpH.hl <- hildreth.lu.order.func(doleprcpH.orc$rho, doleprcpH, order=1)
  doleprcpH.slope[i] <- as.numeric(doleprcpH.hl$coefficients[2])
  doleprcpH.summ <- summary(doleprcpH.hl)
  doleprcpH.p[i] <- doleprcpH.summ$coefficients[8]
  doleprcpH.R2[i] <- doleprcpH.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    doleprcpH.rs <- lm(all.dat.sc.df$dDE[resample.ind] ~ all.dat.sc.df$prcpH[resample.ind])
    doleprcpH.slope.rs[i] <- as.numeric(doleprcpH.rs$coefficients[2])
    doleprcpH.summ.rs <- summary(doleprcpH.rs)
    doleprcpH.p.rs[i] <- doleprcpH.summ.rs$coefficients[8]
    doleprcpH.R2.rs[i] <- doleprcpH.summ.rs$adj.r.squared
    doleprcpH.DWp.rs[i] <- as.numeric(dwtest(doleprcpH.rs)[4])
  
  # dole vs. tree
  doletree <- lm(all.dat.sc.df$dDE ~ all.dat.sc.df$tree)
  doletree.orc <- cochrane.orcutt(doletree, convergence = 5, max.iter=1000)
  doletree.hl <- hildreth.lu.order.func(doletree.orc$rho, doletree, order=1)
  doletree.slope[i] <- as.numeric(doletree.hl$coefficients[2])
  doletree.summ <- summary(doletree.hl)
  doletree.p[i] <- doletree.summ$coefficients[8]
  doletree.R2[i] <- doletree.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    doletree.rs <- lm(all.dat.sc.df$dDE[resample.ind] ~ all.dat.sc.df$tree[resample.ind])
    doletree.slope.rs[i] <- as.numeric(doletree.rs$coefficients[2])
    doletree.summ.rs <- summary(doletree.rs)
    doletree.p.rs[i] <- doletree.summ.rs$coefficients[8]
    doletree.R2.rs[i] <- doletree.summ.rs$adj.r.squared
    doletree.DWp.rs[i] <- as.numeric(dwtest(doletree.rs)[4])
    
  # Hiso vs. insol
  Hisoinsol <- lm(all.dat.sc.df$DHprecip ~ all.dat.sc.df$insol)
  Hisoinsol.orc <- cochrane.orcutt(Hisoinsol, convergence = 5, max.iter=1000)
  Hisoinsol.hl <- hildreth.lu.order.func(Hisoinsol.orc$rho, Hisoinsol, order=1)
  Hisoinsol.slope[i] <- as.numeric(Hisoinsol.hl$coefficients[2])
  Hisoinsol.summ <- summary(Hisoinsol.hl)
  Hisoinsol.p[i] <- Hisoinsol.summ$coefficients[8]
  Hisoinsol.R2[i] <- Hisoinsol.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    Hisoinsol.rs <- lm(all.dat.sc.df$DHprecip[resample.ind] ~ all.dat.sc.df$insol[resample.ind])
    Hisoinsol.slope.rs[i] <- as.numeric(Hisoinsol.rs$coefficients[2])
    Hisoinsol.summ.rs <- summary(Hisoinsol.rs)
    Hisoinsol.p.rs[i] <- Hisoinsol.summ.rs$coefficients[8]
    Hisoinsol.R2.rs[i] <- Hisoinsol.summ.rs$adj.r.squared
    Hisoinsol.DWp.rs[i] <- as.numeric(dwtest(Hisoinsol.rs)[4])
  
  # Hiso vs. dole
  Hisodole <- lm(all.dat.sc.df$DHprecip ~ all.dat.sc.df$dDE)
  Hisodole.orc <- cochrane.orcutt(Hisodole, convergence = 5, max.iter=1000)
  Hisodole.hl <- hildreth.lu.order.func(Hisodole.orc$rho, Hisodole, order=1)
  Hisodole.slope[i] <- as.numeric(Hisodole.hl$coefficients[2])
  Hisodole.summ <- summary(Hisodole.hl)
  Hisodole.p[i] <- Hisodole.summ$coefficients[8]
  Hisodole.R2[i] <- Hisodole.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    Hisodole.rs <- lm(all.dat.sc.df$DHprecip[resample.ind] ~ all.dat.sc.df$dDE[resample.ind])
    Hisodole.slope.rs[i] <- as.numeric(Hisodole.rs$coefficients[2])
    Hisodole.summ.rs <- summary(Hisodole.rs)
    Hisodole.p.rs[i] <- Hisodole.summ.rs$coefficients[8]
    Hisodole.R2.rs[i] <- Hisodole.summ.rs$adj.r.squared
    Hisodole.DWp.rs[i] <- as.numeric(dwtest(Hisodole.rs)[4])
  
  #  Hiso vs. prcpGcor.L
  HisoprcpGcorL <- lm(all.dat.sc.df$DHprecip ~ all.dat.sc.df$prcpGcor.L)
  HisoprcpGcorL.orc <- cochrane.orcutt(HisoprcpGcorL, convergence = 5, max.iter=1000)
  HisoprcpGcorL.hl <- hildreth.lu.order.func(HisoprcpGcorL.orc$rho, HisoprcpGcorL, order=1)
  HisoprcpGcorL.slope[i] <- as.numeric(HisoprcpGcorL.hl$coefficients[2])
  HisoprcpGcorL.summ <- summary(HisoprcpGcorL.hl)
  HisoprcpGcorL.p[i] <- HisoprcpGcorL.summ$coefficients[8]
  HisoprcpGcorL.R2[i] <- HisoprcpGcorL.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    HisoprcpGcorL.rs <- lm(all.dat.sc.df$DHprecip[resample.ind] ~ all.dat.sc.df$prcpGcor.L[resample.ind])
    HisoprcpGcorL.slope.rs[i] <- as.numeric(HisoprcpGcorL.rs$coefficients[2])
    HisoprcpGcorL.summ.rs <- summary(HisoprcpGcorL.rs)
    HisoprcpGcorL.p.rs[i] <- HisoprcpGcorL.summ.rs$coefficients[8]
    HisoprcpGcorL.R2.rs[i] <- HisoprcpGcorL.summ.rs$adj.r.squared
    HisoprcpGcorL.DWp.rs[i] <- as.numeric(dwtest(HisoprcpGcorL.rs)[4])
  
  # Hiso vs. prcpG.L
  HisoprcpGL <- lm(all.dat.sc.df$DHprecip ~ all.dat.sc.df$prcpG.L)
  HisoprcpGL.orc <- cochrane.orcutt(HisoprcpGL, convergence = 5, max.iter=1000)
  HisoprcpGL.hl <- hildreth.lu.order.func(HisoprcpGL.orc$rho, HisoprcpGL, order=1)
  HisoprcpGL.slope[i] <- as.numeric(HisoprcpGL.hl$coefficients[2])
  HisoprcpGL.summ <- summary(HisoprcpGL.hl)
  HisoprcpGL.p[i] <- HisoprcpGL.summ$coefficients[8]
  HisoprcpGL.R2[i] <- HisoprcpGL.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    HisoprcpGL.rs <- lm(all.dat.sc.df$DHprecip[resample.ind] ~ all.dat.sc.df$prcpG.L[resample.ind])
    HisoprcpGL.slope.rs[i] <- as.numeric(HisoprcpGL.rs$coefficients[2])
    HisoprcpGL.summ.rs <- summary(HisoprcpGL.rs)
    HisoprcpGL.p.rs[i] <- HisoprcpGL.summ.rs$coefficients[8]
    HisoprcpGL.R2.rs[i] <- HisoprcpGL.summ.rs$adj.r.squared
    HisoprcpGL.DWp.rs[i] <- as.numeric(dwtest(HisoprcpGL.rs)[4])
  
  # Hiso vs. prcp.L
  HisoprcpL <- lm(all.dat.sc.df$DHprecip ~ all.dat.sc.df$prcpL)
  HisoprcpL.orc <- cochrane.orcutt(HisoprcpL, convergence = 5, max.iter=1000)
  HisoprcpL.hl <- hildreth.lu.order.func(HisoprcpL.orc$rho, HisoprcpL, order=1)
  HisoprcpL.slope[i] <- as.numeric(HisoprcpL.hl$coefficients[2])
  HisoprcpL.summ <- summary(HisoprcpL.hl)
  HisoprcpL.p[i] <- HisoprcpL.summ$coefficients[8]
  HisoprcpL.R2[i] <- HisoprcpL.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    HisoprcpL.rs <- lm(all.dat.sc.df$DHprecip[resample.ind] ~ all.dat.sc.df$prcpL[resample.ind])
    HisoprcpL.slope.rs[i] <- as.numeric(HisoprcpL.rs$coefficients[2])
    HisoprcpL.summ.rs <- summary(HisoprcpL.rs)
    HisoprcpL.p.rs[i] <- HisoprcpL.summ.rs$coefficients[8]
    HisoprcpL.R2.rs[i] <- HisoprcpL.summ.rs$adj.r.squared
    HisoprcpL.DWp.rs[i] <- as.numeric(dwtest(HisoprcpL.rs)[4])
  
  #  Hiso vs. prcpGcor.H
  HisoprcpGcorH <- lm(all.dat.sc.df$DHprecip ~ all.dat.sc.df$prcpGcor.H)
  HisoprcpGcorH.orc <- cochrane.orcutt(HisoprcpGcorH, convergence = 5, max.iter=1000)
  HisoprcpGcorH.hl <- hildreth.lu.order.func(HisoprcpGcorH.orc$rho, HisoprcpGcorH, order=1)
  HisoprcpGcorH.slope[i] <- as.numeric(HisoprcpGcorH.hl$coefficients[2])
  HisoprcpGcorH.summ <- summary(HisoprcpGcorH.hl)
  HisoprcpGcorH.p[i] <- HisoprcpGcorH.summ$coefficients[8]
  HisoprcpGcorH.R2[i] <- HisoprcpGcorH.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    HisoprcpGcorH.rs <- lm(all.dat.sc.df$DHprecip[resample.ind] ~ all.dat.sc.df$prcpGcor.H[resample.ind])
    HisoprcpGcorH.slope.rs[i] <- as.numeric(HisoprcpGcorH.rs$coefficients[2])
    HisoprcpGcorH.summ.rs <- summary(HisoprcpGcorH.rs)
    HisoprcpGcorH.p.rs[i] <- HisoprcpGcorH.summ.rs$coefficients[8]
    HisoprcpGcorH.R2.rs[i] <- HisoprcpGcorH.summ.rs$adj.r.squared
    HisoprcpGcorH.DWp.rs[i] <- as.numeric(dwtest(HisoprcpGcorH.rs)[4])
    
  # Hiso vs. prcpG.H
  HisoprcpGH <- lm(all.dat.sc.df$DHprecip ~ all.dat.sc.df$prcpG.H)
  HisoprcpGH.orc <- cochrane.orcutt(HisoprcpGH, convergence = 5, max.iter=1000)
  HisoprcpGH.hl <- hildreth.lu.order.func(HisoprcpGH.orc$rho, HisoprcpGH, order=1)
  HisoprcpGH.slope[i] <- as.numeric(HisoprcpGH.hl$coefficients[2])
  HisoprcpGH.summ <- summary(HisoprcpGH.hl)
  HisoprcpGH.p[i] <- HisoprcpGH.summ$coefficients[8]
  HisoprcpGH.R2[i] <- HisoprcpGH.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    HisoprcpGH.rs <- lm(all.dat.sc.df$DHprecip[resample.ind] ~ all.dat.sc.df$prcpG.H[resample.ind])
    HisoprcpGH.slope.rs[i] <- as.numeric(HisoprcpGH.rs$coefficients[2])
    HisoprcpGH.summ.rs <- summary(HisoprcpGH.rs)
    HisoprcpGH.p.rs[i] <- HisoprcpGH.summ.rs$coefficients[8]
    HisoprcpGH.R2.rs[i] <- HisoprcpGH.summ.rs$adj.r.squared
    HisoprcpGH.DWp.rs[i] <- as.numeric(dwtest(HisoprcpGH.rs)[4])
    
  # Hiso vs. prcp.H
  HisoprcpH <- lm(all.dat.sc.df$DHprecip ~ all.dat.sc.df$prcpH)
  HisoprcpH.orc <- cochrane.orcutt(HisoprcpH, convergence = 5, max.iter=1000)
  HisoprcpH.hl <- hildreth.lu.order.func(HisoprcpH.orc$rho, HisoprcpH, order=1)
  HisoprcpH.slope[i] <- as.numeric(HisoprcpH.hl$coefficients[2])
  HisoprcpH.summ <- summary(HisoprcpH.hl)
  HisoprcpH.p[i] <- HisoprcpH.summ$coefficients[8]
  HisoprcpH.R2[i] <- HisoprcpH.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    HisoprcpH.rs <- lm(all.dat.sc.df$DHprecip[resample.ind] ~ all.dat.sc.df$prcpH[resample.ind])
    HisoprcpH.slope.rs[i] <- as.numeric(HisoprcpH.rs$coefficients[2])
    HisoprcpH.summ.rs <- summary(HisoprcpH.rs)
    HisoprcpH.p.rs[i] <- HisoprcpH.summ.rs$coefficients[8]
    HisoprcpH.R2.rs[i] <- HisoprcpH.summ.rs$adj.r.squared
    HisoprcpH.DWp.rs[i] <- as.numeric(dwtest(HisoprcpH.rs)[4])
  
  # Hiso vs. tree
  Hisotree <- lm(all.dat.sc.df$DHprecip ~ all.dat.sc.df$tree)
  Hisotree.orc <- cochrane.orcutt(Hisotree, convergence = 5, max.iter=1000)
  Hisotree.hl <- hildreth.lu.order.func(Hisotree.orc$rho, Hisotree, order=1)
  Hisotree.slope[i] <- as.numeric(Hisotree.hl$coefficients[2])
  Hisotree.summ <- summary(Hisotree.hl)
  Hisotree.p[i] <- Hisotree.summ$coefficients[8]
  Hisotree.R2[i] <- Hisotree.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    Hisotree.rs <- lm(all.dat.sc.df$DHprecip[resample.ind] ~ all.dat.sc.df$tree[resample.ind])
    Hisotree.slope.rs[i] <- as.numeric(Hisotree.rs$coefficients[2])
    Hisotree.summ.rs <- summary(Hisotree.rs)
    Hisotree.p.rs[i] <- Hisotree.summ.rs$coefficients[8]
    Hisotree.R2.rs[i] <- Hisotree.summ.rs$adj.r.squared
    Hisotree.DWp.rs[i] <- as.numeric(dwtest(Hisotree.rs)[4])
  
  #  insol vs. prcpGcor.L
  insolprcpGcorL <- lm(all.dat.sc.df$insol ~ all.dat.sc.df$prcpGcor.L)
  insolprcpGcorL.orc <- cochrane.orcutt(insolprcpGcorL, convergence = 5, max.iter=1000)
  insolprcpGcorL.hl <- hildreth.lu.order.func(insolprcpGcorL.orc$rho, insolprcpGcorL, order=1)
  insolprcpGcorL.slope[i] <- as.numeric(insolprcpGcorL.hl$coefficients[2])
  insolprcpGcorL.summ <- summary(insolprcpGcorL.hl)
  insolprcpGcorL.p[i] <- insolprcpGcorL.summ$coefficients[8]
  insolprcpGcorL.R2[i] <- insolprcpGcorL.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    insolprcpGcorL.rs <- lm(all.dat.sc.df$insol[resample.ind] ~ all.dat.sc.df$prcpGcor.L[resample.ind])
    insolprcpGcorL.slope.rs[i] <- as.numeric(insolprcpGcorL.rs$coefficients[2])
    insolprcpGcorL.summ.rs <- summary(insolprcpGcorL.rs)
    insolprcpGcorL.p.rs[i] <- insolprcpGcorL.summ.rs$coefficients[8]
    insolprcpGcorL.R2.rs[i] <- insolprcpGcorL.summ.rs$adj.r.squared
    insolprcpGcorL.DWp.rs[i] <- as.numeric(dwtest(insolprcpGcorL.rs)[4])
  
  # insol vs. prcpG.L
  insolprcpGL <- lm(all.dat.sc.df$insol ~ all.dat.sc.df$prcpG.L)
  insolprcpGL.orc <- cochrane.orcutt(insolprcpGL, convergence = 5, max.iter=1000)
  insolprcpGL.hl <- hildreth.lu.order.func(insolprcpGL.orc$rho, insolprcpGL, order=1)
  insolprcpGL.slope[i] <- as.numeric(insolprcpGL.hl$coefficients[2])
  insolprcpGL.summ <- summary(insolprcpGL.hl)
  insolprcpGL.p[i] <- insolprcpGL.summ$coefficients[8]
  insolprcpGL.R2[i] <- insolprcpGL.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    insolprcpGL.rs <- lm(all.dat.sc.df$insol[resample.ind] ~ all.dat.sc.df$prcpG.L[resample.ind])
    insolprcpGL.slope.rs[i] <- as.numeric(insolprcpGL.rs$coefficients[2])
    insolprcpGL.summ.rs <- summary(insolprcpGL.rs)
    insolprcpGL.p.rs[i] <- insolprcpGL.summ.rs$coefficients[8]
    insolprcpGL.R2.rs[i] <- insolprcpGL.summ.rs$adj.r.squared
    insolprcpGL.DWp.rs[i] <- as.numeric(dwtest(insolprcpGL.rs)[4])
  
  # insol vs. prcpL
  insolprcpL <- lm(all.dat.sc.df$insol ~ all.dat.sc.df$prcpL)
  insolprcpL.orc <- cochrane.orcutt(insolprcpL, convergence = 5, max.iter=1000)
  insolprcpL.hl <- hildreth.lu.order.func(insolprcpL.orc$rho, insolprcpL, order=1)
  insolprcpL.slope[i] <- as.numeric(insolprcpL.hl$coefficients[2])
  insolprcpL.summ <- summary(insolprcpL.hl)
  insolprcpL.p[i] <- insolprcpL.summ$coefficients[8]
  insolprcpL.R2[i] <- insolprcpL.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    insolprcpL.rs <- lm(all.dat.sc.df$insol[resample.ind] ~ all.dat.sc.df$prcpL[resample.ind])
    insolprcpL.slope.rs[i] <- as.numeric(insolprcpL.rs$coefficients[2])
    insolprcpL.summ.rs <- summary(insolprcpL.rs)
    insolprcpL.p.rs[i] <- insolprcpL.summ.rs$coefficients[8]
    insolprcpL.R2.rs[i] <- insolprcpL.summ.rs$adj.r.squared
    insolprcpL.DWp.rs[i] <- as.numeric(dwtest(insolprcpL.rs)[4])
  
  #  insol vs. prcpGcor.H
  insolprcpGcorH <- lm(all.dat.sc.df$insol ~ all.dat.sc.df$prcpGcor.H)
  insolprcpGcorH.orc <- cochrane.orcutt(insolprcpGcorH, convergence = 5, max.iter=1000)
  insolprcpGcorH.hl <- hildreth.lu.order.func(insolprcpGcorH.orc$rho, insolprcpGcorH, order=1)
  insolprcpGcorH.slope[i] <- as.numeric(insolprcpGcorH.hl$coefficients[2])
  insolprcpGcorH.summ <- summary(insolprcpGcorH.hl)
  insolprcpGcorH.p[i] <- insolprcpGcorH.summ$coefficients[8]
  insolprcpGcorH.R2[i] <- insolprcpGcorH.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    insolprcpGcorH.rs <- lm(all.dat.sc.df$insol[resample.ind] ~ all.dat.sc.df$prcpGcor.H[resample.ind])
    insolprcpGcorH.slope.rs[i] <- as.numeric(insolprcpGcorH.rs$coefficients[2])
    insolprcpGcorH.summ.rs <- summary(insolprcpGcorH.rs)
    insolprcpGcorH.p.rs[i] <- insolprcpGcorH.summ.rs$coefficients[8]
    insolprcpGcorH.R2.rs[i] <- insolprcpGcorH.summ.rs$adj.r.squared
    insolprcpGcorH.DWp.rs[i] <- as.numeric(dwtest(insolprcpGcorH.rs)[4])
  
  # insol vs. prcpG.H
  insolprcpGH <- lm(all.dat.sc.df$insol ~ all.dat.sc.df$prcpG.H)
  insolprcpGH.orc <- cochrane.orcutt(insolprcpGH, convergence = 5, max.iter=1000)
  insolprcpGH.hl <- hildreth.lu.order.func(insolprcpGH.orc$rho, insolprcpGH, order=1)
  insolprcpGH.slope[i] <- as.numeric(insolprcpGH.hl$coefficients[2])
  insolprcpGH.summ <- summary(insolprcpGH.hl)
  insolprcpGH.p[i] <- insolprcpGH.summ$coefficients[8]
  insolprcpGH.R2[i] <- insolprcpGH.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    insolprcpGH.rs <- lm(all.dat.sc.df$insol[resample.ind] ~ all.dat.sc.df$prcpG.H[resample.ind])
    insolprcpGH.slope.rs[i] <- as.numeric(insolprcpGH.rs$coefficients[2])
    insolprcpGH.summ.rs <- summary(insolprcpGH.rs)
    insolprcpGH.p.rs[i] <- insolprcpGH.summ.rs$coefficients[8]
    insolprcpGH.R2.rs[i] <- insolprcpGH.summ.rs$adj.r.squared
    insolprcpGH.DWp.rs[i] <- as.numeric(dwtest(insolprcpGH.rs)[4])
    
  # insol vs. prcp.H
  insolprcpH <- lm(all.dat.sc.df$insol ~ all.dat.sc.df$prcpH)
  insolprcpH.orc <- cochrane.orcutt(insolprcpH, convergence = 5, max.iter=1000)
  insolprcpH.hl <- hildreth.lu.order.func(insolprcpH.orc$rho, insolprcpH, order=1)
  insolprcpH.slope[i] <- as.numeric(insolprcpH.hl$coefficients[2])
  insolprcpH.summ <- summary(insolprcpH.hl)
  insolprcpH.p[i] <- insolprcpH.summ$coefficients[8]
  insolprcpH.R2[i] <- insolprcpH.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    insolprcpH.rs <- lm(all.dat.sc.df$insol[resample.ind] ~ all.dat.sc.df$prcpH[resample.ind])
    insolprcpH.slope.rs[i] <- as.numeric(insolprcpH.rs$coefficients[2])
    insolprcpH.summ.rs <- summary(insolprcpH.rs)
    insolprcpH.p.rs[i] <- insolprcpH.summ.rs$coefficients[8]
    insolprcpH.R2.rs[i] <- insolprcpH.summ.rs$adj.r.squared
    insolprcpH.DWp.rs[i] <- as.numeric(dwtest(insolprcpH.rs)[4])
    
    # insol vs. tree
  insoltree <- lm(all.dat.sc.df$insol ~ all.dat.sc.df$tree)
  insoltree.orc <- cochrane.orcutt(insoltree, convergence = 5, max.iter=1000)
  insoltree.hl <- hildreth.lu.order.func(insoltree.orc$rho, insoltree, order=1)
  insoltree.slope[i] <- as.numeric(insoltree.hl$coefficients[2])
  insoltree.summ <- summary(insoltree.hl)
  insoltree.p[i] <- insoltree.summ$coefficients[8]
  insoltree.R2[i] <- insoltree.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    insoltree.rs <- lm(all.dat.sc.df$insol[resample.ind] ~ all.dat.sc.df$tree[resample.ind])
    insoltree.slope.rs[i] <- as.numeric(insoltree.rs$coefficients[2])
    insoltree.summ.rs <- summary(insoltree.rs)
    insoltree.p.rs[i] <- insoltree.summ.rs$coefficients[8]
    insoltree.R2.rs[i] <- insoltree.summ.rs$adj.r.squared
    insoltree.DWp.rs[i] <- as.numeric(dwtest(insoltree.rs)[4])
  
  # tree vs. prcpG.L
  treeprcpGL <- lm(all.dat.sc.df$tree ~ all.dat.sc.df$prcpG.L)
  treeprcpGL.orc <- cochrane.orcutt(treeprcpGL, convergence = 5, max.iter=1000)
  treeprcpGL.hl <- hildreth.lu.order.func(treeprcpGL.orc$rho, treeprcpGL, order=1)
  treeprcpGL.slope[i] <- as.numeric(treeprcpGL.hl$coefficients[2])
  treeprcpGL.summ <- summary(treeprcpGL.hl)
  treeprcpGL.p[i] <- treeprcpGL.summ$coefficients[8]
  treeprcpGL.R2[i] <- treeprcpGL.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    treeprcpGL.rs <- lm(all.dat.sc.df$tree[resample.ind] ~ all.dat.sc.df$prcpG.L[resample.ind])
    treeprcpGL.slope.rs[i] <- as.numeric(treeprcpGL.rs$coefficients[2])
    treeprcpGL.summ.rs <- summary(treeprcpGL.rs)
    treeprcpGL.p.rs[i] <- treeprcpGL.summ.rs$coefficients[8]
    treeprcpGL.R2.rs[i] <- treeprcpGL.summ.rs$adj.r.squared
    treeprcpGL.DWp.rs[i] <- as.numeric(dwtest(treeprcpGL.rs)[4])
  
  # tree vs. prcpL
  treeprcpL <- lm(all.dat.sc.df$tree ~ all.dat.sc.df$prcpL)
  treeprcpL.orc <- cochrane.orcutt(treeprcpL, convergence = 5, max.iter=1000)
  treeprcpL.hl <- hildreth.lu.order.func(treeprcpL.orc$rho, treeprcpL, order=1)
  treeprcpL.slope[i] <- as.numeric(treeprcpL.hl$coefficients[2])
  treeprcpL.summ <- summary(treeprcpL.hl)
  treeprcpL.p[i] <- treeprcpL.summ$coefficients[8]
  treeprcpL.R2[i] <- treeprcpL.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    treeprcpL.rs <- lm(all.dat.sc.df$tree[resample.ind] ~ all.dat.sc.df$prcpL[resample.ind])
    treeprcpL.slope.rs[i] <- as.numeric(treeprcpL.rs$coefficients[2])
    treeprcpL.summ.rs <- summary(treeprcpL.rs)
    treeprcpL.p.rs[i] <- treeprcpL.summ.rs$coefficients[8]
    treeprcpL.R2.rs[i] <- treeprcpL.summ.rs$adj.r.squared
    treeprcpL.DWp.rs[i] <- as.numeric(dwtest(treeprcpL.rs)[4])
    
  #  tree vs. prcpGcor.H
  treeprcpGcorH <- lm(all.dat.sc.df$tree ~ all.dat.sc.df$prcpGcor.H)
  treeprcpGcorH.orc <- cochrane.orcutt(treeprcpGcorH, convergence = 5, max.iter=1000)
  treeprcpGcorH.hl <- hildreth.lu.order.func(treeprcpGcorH.orc$rho, treeprcpGcorH, order=1)
  treeprcpGcorH.slope[i] <- as.numeric(treeprcpGcorH.hl$coefficients[2])
  treeprcpGcorH.summ <- summary(treeprcpGcorH.hl)
  treeprcpGcorH.p[i] <- treeprcpGcorH.summ$coefficients[8]
  treeprcpGcorH.R2[i] <- treeprcpGcorH.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    treeprcpGcorH.rs <- lm(all.dat.sc.df$tree[resample.ind] ~ all.dat.sc.df$prcpGcor.H[resample.ind])
    treeprcpGcorH.slope.rs[i] <- as.numeric(treeprcpGcorH.rs$coefficients[2])
    treeprcpGcorH.summ.rs <- summary(treeprcpGcorH.rs)
    treeprcpGcorH.p.rs[i] <- treeprcpGcorH.summ.rs$coefficients[8]
    treeprcpGcorH.R2.rs[i] <- treeprcpGcorH.summ.rs$adj.r.squared
    treeprcpGcorH.DWp.rs[i] <- as.numeric(dwtest(treeprcpGcorH.rs)[4])
  
  #  tree vs. prcpGcor.L
  treeprcpGcorL <- lm(all.dat.sc.df$tree ~ all.dat.sc.df$prcpGcor.L)
  treeprcpGcorL.orc <- cochrane.orcutt(treeprcpGcorL, convergence = 5, max.iter=1000)
  treeprcpGcorL.hl <- hildreth.lu.order.func(treeprcpGcorL.orc$rho, treeprcpGcorL, order=1)
  treeprcpGcorL.slope[i] <- as.numeric(treeprcpGcorL.hl$coefficients[2])
  treeprcpGcorL.summ <- summary(treeprcpGcorL.hl)
  treeprcpGcorL.p[i] <- treeprcpGcorL.summ$coefficients[8]
  treeprcpGcorL.R2[i] <- treeprcpGcorL.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    treeprcpGcorL.rs <- lm(all.dat.sc.df$tree[resample.ind] ~ all.dat.sc.df$prcpGcor.L[resample.ind])
    treeprcpGcorL.slope.rs[i] <- as.numeric(treeprcpGcorL.rs$coefficients[2])
    treeprcpGcorL.summ.rs <- summary(treeprcpGcorL.rs)
    treeprcpGcorL.p.rs[i] <- treeprcpGcorL.summ.rs$coefficients[8]
    treeprcpGcorL.R2.rs[i] <- treeprcpGcorL.summ.rs$adj.r.squared
    treeprcpGcorL.DWp.rs[i] <- as.numeric(dwtest(treeprcpGcorL.rs)[4])
    
  # tree vs. prcpG.H
  treeprcpGH <- lm(all.dat.sc.df$tree ~ all.dat.sc.df$prcpG.H)
  treeprcpGH.orc <- cochrane.orcutt(treeprcpGH, convergence = 5, max.iter=1000)
  treeprcpGH.hl <- hildreth.lu.order.func(treeprcpGH.orc$rho, treeprcpGH, order=1)
  treeprcpGH.slope[i] <- as.numeric(treeprcpGH.hl$coefficients[2])
  treeprcpGH.summ <- summary(treeprcpGH.hl)
  treeprcpGH.p[i] <- treeprcpGH.summ$coefficients[8]
  treeprcpGH.R2[i] <- treeprcpGH.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    treeprcpGH.rs <- lm(all.dat.sc.df$tree[resample.ind] ~ all.dat.sc.df$prcpG.H[resample.ind])
    treeprcpGH.slope.rs[i] <- as.numeric(treeprcpGH.rs$coefficients[2])
    treeprcpGH.summ.rs <- summary(treeprcpGH.rs)
    treeprcpGH.p.rs[i] <- treeprcpGH.summ.rs$coefficients[8]
    treeprcpGH.R2.rs[i] <- treeprcpGH.summ.rs$adj.r.squared
    treeprcpGH.DWp.rs[i] <- as.numeric(dwtest(treeprcpGH.rs)[4])
  
  # tree vs. prcp.H
  treeprcpH <- lm(all.dat.sc.df$tree ~ all.dat.sc.df$prcpH)
  treeprcpH.orc <- cochrane.orcutt(treeprcpH, convergence = 5, max.iter=1000)
  treeprcpH.hl <- hildreth.lu.order.func(treeprcpH.orc$rho, treeprcpH, order=1)
  treeprcpH.slope[i] <- as.numeric(treeprcpH.hl$coefficients[2])
  treeprcpH.summ <- summary(treeprcpH.hl)
  treeprcpH.p[i] <- treeprcpH.summ$coefficients[8]
  treeprcpH.R2[i] <- treeprcpH.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    treeprcpH.rs <- lm(all.dat.sc.df$tree[resample.ind] ~ all.dat.sc.df$prcpH[resample.ind])
    treeprcpH.slope.rs[i] <- as.numeric(treeprcpH.rs$coefficients[2])
    treeprcpH.summ.rs <- summary(treeprcpH.rs)
    treeprcpH.p.rs[i] <- treeprcpH.summ.rs$coefficients[8]
    treeprcpH.R2.rs[i] <- treeprcpH.summ.rs$adj.r.squared
    treeprcpH.DWp.rs[i] <- as.numeric(dwtest(treeprcpH.rs)[4])
  
  # tree vs. prcpL
  prcpLprcpSA <- lm(LOVECLIM.NA.SA.dat.sc.df$prcpL ~ LOVECLIM.NA.SA.dat.sc.df$prcpSA)
  prcpLprcpSA.orc <- cochrane.orcutt(prcpLprcpSA, convergence = 5, max.iter=1000)
  prcpLprcpSA.hl <- hildreth.lu.order.func(prcpLprcpSA.orc$rho, prcpLprcpSA, order=1)
  prcpLprcpSA.slope[i] <- as.numeric(prcpLprcpSA.hl$coefficients[2])
  prcpLprcpSA.summ <- summary(prcpLprcpSA.hl)
  prcpLprcpSA.p[i] <- prcpLprcpSA.summ$coefficients[8]
  prcpLprcpSA.R2[i] <- prcpLprcpSA.summ$adj.r.squared
  
    # resampled
    resample.ind <- sort(NonSeqSample(x=1:(dim(all.dat.sc.df)[1]), size=resample.n, replace=F))
    prcpLprcpSA.rs <- lm(LOVECLIM.NA.SA.dat.sc.df$prcpL[resample.ind] ~ LOVECLIM.NA.SA.dat.sc.df$prcpSA[resample.ind])
    prcpLprcpSA.slope.rs[i] <- as.numeric(prcpLprcpSA.rs$coefficients[2])
    prcpLprcpSA.summ.rs <- summary(prcpLprcpSA.rs)
    prcpLprcpSA.p.rs[i] <- prcpLprcpSA.summ.rs$coefficients[8]
    prcpLprcpSA.R2.rs[i] <- prcpLprcpSA.summ.rs$adj.r.squared
    prcpLprcpSA.DWp.rs[i] <- as.numeric(dwtest(prcpLprcpSA.rs)[4])
  
  if (i %% itdiv==0) print(i) 
  
} # end i

cor.mean <- apply(allcor.arr, MARGIN=c(1,2), mean, na.rm=T)
colnames(cor.mean) <- colnames(all.dat)
rownames(cor.mean) <- colnames(all.dat)
diag(cor.mean) <- NA
cor.mean[upper.tri(cor.mean)] <- NA
cor.mean <- cor.mean[-1,-dim(cor.mean)[2]]
cor.mean

cor.lo <- apply(allcor.arr, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
colnames(cor.lo) <- colnames(all.dat)
rownames(cor.lo) <- colnames(all.dat)
diag(cor.lo) <- NA
cor.lo[upper.tri(cor.lo)] <- NA
cor.lo <- cor.lo[-1,-dim(cor.lo)[2]]
round(cor.lo, 3)

cor.up <- apply(allcor.arr, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
colnames(cor.up) <- colnames(all.dat)
rownames(cor.up) <- colnames(all.dat)
diag(cor.up) <- NA
cor.up[upper.tri(cor.up)] <- NA
cor.up <- cor.up[-1,-dim(cor.up)[2]]
round(cor.up, 3)


# autocorrelation-corrected outputs
MaO18insol.slope.mn <- mean(MaO18insol.slope, na.rm=T)
MaO18insol.slope.lo <- quantile(MaO18insol.slope, probs=0.025, na.rm=T)
MaO18insol.slope.up <- quantile(MaO18insol.slope, probs=0.975, na.rm=T)
MaO18insol.p.mn <- mean(MaO18insol.p, na.rm=T)
MaO18insol.p.lo <- quantile(MaO18insol.p, probs=0.025, na.rm=T)
MaO18insol.p.up <- quantile(MaO18insol.p, probs=0.975, na.rm=T)
MaO18insol.R2.mn <- mean(MaO18insol.R2, na.rm=T)
MaO18insol.R2.lo <- quantile(MaO18insol.R2, probs=0.025, na.rm=T)
MaO18insol.R2.up <- quantile(MaO18insol.R2, probs=0.975, na.rm=T)

MaO18prcpGcorL.slope.mn <- mean(MaO18prcpGcorL.slope, na.rm=T)
MaO18prcpGcorL.slope.lo <- quantile(MaO18prcpGcorL.slope, probs=0.025, na.rm=T)
MaO18prcpGcorL.slope.up <- quantile(MaO18prcpGcorL.slope, probs=0.975, na.rm=T)
MaO18prcpGcorL.p.mn <- mean(MaO18prcpGcorL.p, na.rm=T)
MaO18prcpGcorL.p.lo <- quantile(MaO18prcpGcorL.p, probs=0.025, na.rm=T)
MaO18prcpGcorL.p.up <- quantile(MaO18prcpGcorL.p, probs=0.975, na.rm=T)
MaO18prcpGcorL.R2.mn <- mean(MaO18prcpGcorL.R2, na.rm=T)
MaO18prcpGcorL.R2.lo <- quantile(MaO18prcpGcorL.R2, probs=0.025, na.rm=T)
MaO18prcpGcorL.R2.up <- quantile(MaO18prcpGcorL.R2, probs=0.975, na.rm=T)

MaO18prcpGL.slope.mn <- mean(MaO18prcpGL.slope, na.rm=T)
MaO18prcpGL.slope.lo <- quantile(MaO18prcpGL.slope, probs=0.025, na.rm=T)
MaO18prcpGL.slope.up <- quantile(MaO18prcpGL.slope, probs=0.975, na.rm=T)
MaO18prcpGL.p.mn <- mean(MaO18prcpGL.p, na.rm=T)
MaO18prcpGL.p.lo <- quantile(MaO18prcpGL.p, probs=0.025, na.rm=T)
MaO18prcpGL.p.up <- quantile(MaO18prcpGL.p, probs=0.975, na.rm=T)
MaO18prcpGL.R2.mn <- mean(MaO18prcpGL.R2, na.rm=T)
MaO18prcpGL.R2.lo <- quantile(MaO18prcpGL.R2, probs=0.025, na.rm=T)
MaO18prcpGL.R2.up <- quantile(MaO18prcpGL.R2, probs=0.975, na.rm=T)

MaO18prcpL.slope.mn <- mean(MaO18prcpL.slope, na.rm=T)
MaO18prcpL.slope.lo <- quantile(MaO18prcpL.slope, probs=0.025, na.rm=T)
MaO18prcpL.slope.up <- quantile(MaO18prcpL.slope, probs=0.975, na.rm=T)
MaO18prcpL.p.mn <- mean(MaO18prcpL.p, na.rm=T)
MaO18prcpL.p.lo <- quantile(MaO18prcpL.p, probs=0.025, na.rm=T)
MaO18prcpL.p.up <- quantile(MaO18prcpL.p, probs=0.975, na.rm=T)
MaO18prcpL.R2.mn <- mean(MaO18prcpL.R2, na.rm=T)
MaO18prcpL.R2.lo <- quantile(MaO18prcpL.R2, probs=0.025, na.rm=T)
MaO18prcpL.R2.up <- quantile(MaO18prcpL.R2, probs=0.975, na.rm=T)

MaO18prcpGcorH.slope.mn <- mean(MaO18prcpGcorH.slope, na.rm=T)
MaO18prcpGcorH.slope.lo <- quantile(MaO18prcpGcorH.slope, probs=0.025, na.rm=T)
MaO18prcpGcorH.slope.up <- quantile(MaO18prcpGcorH.slope, probs=0.975, na.rm=T)
MaO18prcpGcorH.p.mn <- mean(MaO18prcpGcorH.p, na.rm=T)
MaO18prcpGcorH.p.lo <- quantile(MaO18prcpGcorH.p, probs=0.025, na.rm=T)
MaO18prcpGcorH.p.up <- quantile(MaO18prcpGcorH.p, probs=0.975, na.rm=T)
MaO18prcpGcorH.R2.mn <- mean(MaO18prcpGcorH.R2, na.rm=T)
MaO18prcpGcorH.R2.lo <- quantile(MaO18prcpGcorH.R2, probs=0.025, na.rm=T)
MaO18prcpGcorH.R2.up <- quantile(MaO18prcpGcorH.R2, probs=0.975, na.rm=T)

MaO18prcpGH.slope.mn <- mean(MaO18prcpGH.slope, na.rm=T)
MaO18prcpGH.slope.lo <- quantile(MaO18prcpGH.slope, probs=0.025, na.rm=T)
MaO18prcpGH.slope.up <- quantile(MaO18prcpGH.slope, probs=0.975, na.rm=T)
MaO18prcpGH.p.mn <- mean(MaO18prcpGH.p, na.rm=T)
MaO18prcpGH.p.lo <- quantile(MaO18prcpGH.p, probs=0.025, na.rm=T)
MaO18prcpGH.p.up <- quantile(MaO18prcpGH.p, probs=0.975, na.rm=T)
MaO18prcpGH.R2.mn <- mean(MaO18prcpGH.R2, na.rm=T)
MaO18prcpGH.R2.lo <- quantile(MaO18prcpGH.R2, probs=0.025, na.rm=T)
MaO18prcpGH.R2.up <- quantile(MaO18prcpGH.R2, probs=0.975, na.rm=T)

MaO18prcpH.slope.mn <- mean(MaO18prcpH.slope, na.rm=T)
MaO18prcpH.slope.lo <- quantile(MaO18prcpH.slope, probs=0.025, na.rm=T)
MaO18prcpH.slope.up <- quantile(MaO18prcpH.slope, probs=0.975, na.rm=T)
MaO18prcpH.p.mn <- mean(MaO18prcpH.p, na.rm=T)
MaO18prcpH.p.lo <- quantile(MaO18prcpH.p, probs=0.025, na.rm=T)
MaO18prcpH.p.up <- quantile(MaO18prcpH.p, probs=0.975, na.rm=T)
MaO18prcpH.R2.mn <- mean(MaO18prcpH.R2, na.rm=T)
MaO18prcpH.R2.lo <- quantile(MaO18prcpH.R2, probs=0.025, na.rm=T)
MaO18prcpH.R2.up <- quantile(MaO18prcpH.R2, probs=0.975, na.rm=T)

MaO18dole.slope.mn <- mean(MaO18dole.slope, na.rm=T)
MaO18dole.slope.lo <- quantile(MaO18dole.slope, probs=0.025, na.rm=T)
MaO18dole.slope.up <- quantile(MaO18dole.slope, probs=0.975, na.rm=T)
MaO18dole.p.mn <- mean(MaO18dole.p, na.rm=T)
MaO18dole.p.lo <- quantile(MaO18dole.p, probs=0.025, na.rm=T)
MaO18dole.p.up <- quantile(MaO18dole.p, probs=0.975, na.rm=T)
MaO18dole.R2.mn <- mean(MaO18dole.R2, na.rm=T)
MaO18dole.R2.lo <- quantile(MaO18dole.R2, probs=0.025, na.rm=T)
MaO18dole.R2.up <- quantile(MaO18dole.R2, probs=0.975, na.rm=T)

MaO18tree.slope.mn <- mean(MaO18tree.slope, na.rm=T)
MaO18tree.slope.lo <- quantile(MaO18tree.slope, probs=0.025, na.rm=T)
MaO18tree.slope.up <- quantile(MaO18tree.slope, probs=0.975, na.rm=T)
MaO18tree.p.mn <- mean(MaO18tree.p, na.rm=T)
MaO18tree.p.lo <- quantile(MaO18tree.p, probs=0.025, na.rm=T)
MaO18tree.p.up <- quantile(MaO18tree.p, probs=0.975, na.rm=T)
MaO18tree.R2.mn <- mean(MaO18tree.R2, na.rm=T)
MaO18tree.R2.lo <- quantile(MaO18tree.R2, probs=0.025, na.rm=T)
MaO18tree.R2.up <- quantile(MaO18tree.R2, probs=0.975, na.rm=T)

MaO18Hiso.slope.mn <- mean(MaO18Hiso.slope, na.rm=T)
MaO18Hiso.slope.lo <- quantile(MaO18Hiso.slope, probs=0.025, na.rm=T)
MaO18Hiso.slope.up <- quantile(MaO18Hiso.slope, probs=0.975, na.rm=T)
MaO18Hiso.p.mn <- mean(MaO18Hiso.p, na.rm=T)
MaO18Hiso.p.lo <- quantile(MaO18Hiso.p, probs=0.025, na.rm=T)
MaO18Hiso.p.up <- quantile(MaO18Hiso.p, probs=0.975, na.rm=T)
MaO18Hiso.R2.mn <- mean(MaO18Hiso.R2, na.rm=T)
MaO18Hiso.R2.lo <- quantile(MaO18Hiso.R2, probs=0.025, na.rm=T)
MaO18Hiso.R2.up <- quantile(MaO18Hiso.R2, probs=0.975, na.rm=T)


doleinsol.slope.mn <- mean(doleinsol.slope, na.rm=T)
doleinsol.slope.lo <- quantile(doleinsol.slope, probs=0.025, na.rm=T)
doleinsol.slope.up <- quantile(doleinsol.slope, probs=0.975, na.rm=T)
doleinsol.p.mn <- mean(doleinsol.p, na.rm=T)
doleinsol.p.lo <- quantile(doleinsol.p, probs=0.025, na.rm=T)
doleinsol.p.up <- quantile(doleinsol.p, probs=0.975, na.rm=T)
doleinsol.R2.mn <- mean(doleinsol.R2, na.rm=T)
doleinsol.R2.lo <- quantile(doleinsol.R2, probs=0.025, na.rm=T)
doleinsol.R2.up <- quantile(doleinsol.R2, probs=0.975, na.rm=T)

doleprcpGcorL.slope.mn <- mean(doleprcpGcorL.slope, na.rm=T)
doleprcpGcorL.slope.lo <- quantile(doleprcpGcorL.slope, probs=0.025, na.rm=T)
doleprcpGcorL.slope.up <- quantile(doleprcpGcorL.slope, probs=0.975, na.rm=T)
doleprcpGcorL.p.mn <- mean(doleprcpGcorL.p, na.rm=T)
doleprcpGcorL.p.lo <- quantile(doleprcpGcorL.p, probs=0.025, na.rm=T)
doleprcpGcorL.p.up <- quantile(doleprcpGcorL.p, probs=0.975, na.rm=T)
doleprcpGcorL.R2.mn <- mean(doleprcpGcorL.R2, na.rm=T)
doleprcpGcorL.R2.lo <- quantile(doleprcpGcorL.R2, probs=0.025, na.rm=T)
doleprcpGcorL.R2.up <- quantile(doleprcpGcorL.R2, probs=0.975, na.rm=T)

doleprcpGL.slope.mn <- mean(doleprcpGL.slope, na.rm=T)
doleprcpGL.slope.lo <- quantile(doleprcpGL.slope, probs=0.025, na.rm=T)
doleprcpGL.slope.up <- quantile(doleprcpGL.slope, probs=0.975, na.rm=T)
doleprcpGL.p.mn <- mean(doleprcpGL.p, na.rm=T)
doleprcpGL.p.lo <- quantile(doleprcpGL.p, probs=0.025, na.rm=T)
doleprcpGL.p.up <- quantile(doleprcpGL.p, probs=0.975, na.rm=T)
doleprcpGL.R2.mn <- mean(doleprcpGL.R2, na.rm=T)
doleprcpGL.R2.lo <- quantile(doleprcpGL.R2, probs=0.025, na.rm=T)
doleprcpGL.R2.up <- quantile(doleprcpGL.R2, probs=0.975, na.rm=T)

doleprcpL.slope.mn <- mean(doleprcpL.slope, na.rm=T)
doleprcpL.slope.lo <- quantile(doleprcpL.slope, probs=0.025, na.rm=T)
doleprcpL.slope.up <- quantile(doleprcpL.slope, probs=0.975, na.rm=T)
doleprcpL.p.mn <- mean(doleprcpL.p, na.rm=T)
doleprcpL.p.lo <- quantile(doleprcpL.p, probs=0.025, na.rm=T)
doleprcpL.p.up <- quantile(doleprcpL.p, probs=0.975, na.rm=T)
doleprcpL.R2.mn <- mean(doleprcpL.R2, na.rm=T)
doleprcpL.R2.lo <- quantile(doleprcpL.R2, probs=0.025, na.rm=T)
doleprcpL.R2.up <- quantile(doleprcpL.R2, probs=0.975, na.rm=T)

doleprcpGcorH.slope.mn <- mean(doleprcpGcorH.slope, na.rm=T)
doleprcpGcorH.slope.lo <- quantile(doleprcpGcorH.slope, probs=0.025, na.rm=T)
doleprcpGcorH.slope.up <- quantile(doleprcpGcorH.slope, probs=0.975, na.rm=T)
doleprcpGcorH.p.mn <- mean(doleprcpGcorH.p, na.rm=T)
doleprcpGcorH.p.lo <- quantile(doleprcpGcorH.p, probs=0.025, na.rm=T)
doleprcpGcorH.p.up <- quantile(doleprcpGcorH.p, probs=0.975, na.rm=T)
doleprcpGcorH.R2.mn <- mean(doleprcpGcorH.R2, na.rm=T)
doleprcpGcorH.R2.lo <- quantile(doleprcpGcorH.R2, probs=0.025, na.rm=T)
doleprcpGcorH.R2.up <- quantile(doleprcpGcorH.R2, probs=0.975, na.rm=T)

doleprcpGH.slope.mn <- mean(doleprcpGH.slope, na.rm=T)
doleprcpGH.slope.lo <- quantile(doleprcpGH.slope, probs=0.025, na.rm=T)
doleprcpGH.slope.up <- quantile(doleprcpGH.slope, probs=0.975, na.rm=T)
doleprcpGH.p.mn <- mean(doleprcpGH.p, na.rm=T)
doleprcpGH.p.lo <- quantile(doleprcpGH.p, probs=0.025, na.rm=T)
doleprcpGH.p.up <- quantile(doleprcpGH.p, probs=0.975, na.rm=T)
doleprcpGH.R2.mn <- mean(doleprcpGH.R2, na.rm=T)
doleprcpGH.R2.lo <- quantile(doleprcpGH.R2, probs=0.025, na.rm=T)
doleprcpGH.R2.up <- quantile(doleprcpGH.R2, probs=0.975, na.rm=T)

doleprcpH.slope.mn <- mean(doleprcpH.slope, na.rm=T)
doleprcpH.slope.lo <- quantile(doleprcpH.slope, probs=0.025, na.rm=T)
doleprcpH.slope.up <- quantile(doleprcpH.slope, probs=0.975, na.rm=T)
doleprcpH.p.mn <- mean(doleprcpH.p, na.rm=T)
doleprcpH.p.lo <- quantile(doleprcpH.p, probs=0.025, na.rm=T)
doleprcpH.p.up <- quantile(doleprcpH.p, probs=0.975, na.rm=T)
doleprcpH.R2.mn <- mean(doleprcpH.R2, na.rm=T)
doleprcpH.R2.lo <- quantile(doleprcpH.R2, probs=0.025, na.rm=T)
doleprcpH.R2.up <- quantile(doleprcpH.R2, probs=0.975, na.rm=T)

doletree.slope.mn <- mean(doletree.slope, na.rm=T)
doletree.slope.lo <- quantile(doletree.slope, probs=0.025, na.rm=T)
doletree.slope.up <- quantile(doletree.slope, probs=0.975, na.rm=T)
doletree.p.mn <- mean(doletree.p, na.rm=T)
doletree.p.lo <- quantile(doletree.p, probs=0.025, na.rm=T)
doletree.p.up <- quantile(doletree.p, probs=0.975, na.rm=T)
doletree.R2.mn <- mean(doletree.R2, na.rm=T)
doletree.R2.lo <- quantile(doletree.R2, probs=0.025, na.rm=T)
doletree.R2.up <- quantile(doletree.R2, probs=0.975, na.rm=T)

Hisoinsol.slope.mn <- mean(Hisoinsol.slope, na.rm=T)
Hisoinsol.slope.lo <- quantile(Hisoinsol.slope, probs=0.025, na.rm=T)
Hisoinsol.slope.up <- quantile(Hisoinsol.slope, probs=0.975, na.rm=T)
Hisoinsol.p.mn <- mean(Hisoinsol.p, na.rm=T)
Hisoinsol.p.lo <- quantile(Hisoinsol.p, probs=0.025, na.rm=T)
Hisoinsol.p.up <- quantile(Hisoinsol.p, probs=0.975, na.rm=T)
Hisoinsol.R2.mn <- mean(Hisoinsol.R2, na.rm=T)
Hisoinsol.R2.lo <- quantile(Hisoinsol.R2, probs=0.025, na.rm=T)
Hisoinsol.R2.up <- quantile(Hisoinsol.R2, probs=0.975, na.rm=T)

Hisodole.slope.mn <- mean(Hisodole.slope, na.rm=T)
Hisodole.slope.lo <- quantile(Hisodole.slope, probs=0.025, na.rm=T)
Hisodole.slope.up <- quantile(Hisodole.slope, probs=0.975, na.rm=T)
Hisodole.p.mn <- mean(Hisodole.p, na.rm=T)
Hisodole.p.lo <- quantile(Hisodole.p, probs=0.025, na.rm=T)
Hisodole.p.up <- quantile(Hisodole.p, probs=0.975, na.rm=T)
Hisodole.R2.mn <- mean(Hisodole.R2, na.rm=T)
Hisodole.R2.lo <- quantile(Hisodole.R2, probs=0.025, na.rm=T)
Hisodole.R2.up <- quantile(Hisodole.R2, probs=0.975, na.rm=T)

HisoprcpGcorL.slope.mn <- mean(HisoprcpGcorL.slope, na.rm=T)
HisoprcpGcorL.slope.lo <- quantile(HisoprcpGcorL.slope, probs=0.025, na.rm=T)
HisoprcpGcorL.slope.up <- quantile(HisoprcpGcorL.slope, probs=0.975, na.rm=T)
HisoprcpGcorL.p.mn <- mean(HisoprcpGcorL.p, na.rm=T)
HisoprcpGcorL.p.lo <- quantile(HisoprcpGcorL.p, probs=0.025, na.rm=T)
HisoprcpGcorL.p.up <- quantile(HisoprcpGcorL.p, probs=0.975, na.rm=T)
HisoprcpGcorL.R2.mn <- mean(HisoprcpGcorL.R2, na.rm=T)
HisoprcpGcorL.R2.lo <- quantile(HisoprcpGcorL.R2, probs=0.025, na.rm=T)
HisoprcpGcorL.R2.up <- quantile(HisoprcpGcorL.R2, probs=0.975, na.rm=T)

HisoprcpGL.slope.mn <- mean(HisoprcpGL.slope, na.rm=T)
HisoprcpGL.slope.lo <- quantile(HisoprcpGL.slope, probs=0.025, na.rm=T)
HisoprcpGL.slope.up <- quantile(HisoprcpGL.slope, probs=0.975, na.rm=T)
HisoprcpGL.p.mn <- mean(HisoprcpGL.p, na.rm=T)
HisoprcpGL.p.lo <- quantile(HisoprcpGL.p, probs=0.025, na.rm=T)
HisoprcpGL.p.up <- quantile(HisoprcpGL.p, probs=0.975, na.rm=T)
HisoprcpGL.R2.mn <- mean(HisoprcpGL.R2, na.rm=T)
HisoprcpGL.R2.lo <- quantile(HisoprcpGL.R2, probs=0.025, na.rm=T)
HisoprcpGL.R2.up <- quantile(HisoprcpGL.R2, probs=0.975, na.rm=T)

HisoprcpL.slope.mn <- mean(HisoprcpL.slope, na.rm=T)
HisoprcpL.slope.lo <- quantile(HisoprcpL.slope, probs=0.025, na.rm=T)
HisoprcpL.slope.up <- quantile(HisoprcpL.slope, probs=0.975, na.rm=T)
HisoprcpL.p.mn <- mean(HisoprcpL.p, na.rm=T)
HisoprcpL.p.lo <- quantile(HisoprcpL.p, probs=0.025, na.rm=T)
HisoprcpL.p.up <- quantile(HisoprcpL.p, probs=0.975, na.rm=T)
HisoprcpL.R2.mn <- mean(HisoprcpL.R2, na.rm=T)
HisoprcpL.R2.lo <- quantile(HisoprcpL.R2, probs=0.025, na.rm=T)
HisoprcpL.R2.up <- quantile(HisoprcpL.R2, probs=0.975, na.rm=T)

HisoprcpGcorH.slope.mn <- mean(HisoprcpGcorH.slope, na.rm=T)
HisoprcpGcorH.slope.lo <- quantile(HisoprcpGcorH.slope, probs=0.025, na.rm=T)
HisoprcpGcorH.slope.up <- quantile(HisoprcpGcorH.slope, probs=0.975, na.rm=T)
HisoprcpGcorH.p.mn <- mean(HisoprcpGcorH.p, na.rm=T)
HisoprcpGcorH.p.lo <- quantile(HisoprcpGcorH.p, probs=0.025, na.rm=T)
HisoprcpGcorH.p.up <- quantile(HisoprcpGcorH.p, probs=0.975, na.rm=T)
HisoprcpGcorH.R2.mn <- mean(HisoprcpGcorH.R2, na.rm=T)
HisoprcpGcorH.R2.lo <- quantile(HisoprcpGcorH.R2, probs=0.025, na.rm=T)
HisoprcpGcorH.R2.up <- quantile(HisoprcpGcorH.R2, probs=0.975, na.rm=T)

HisoprcpGH.slope.mn <- mean(HisoprcpGH.slope, na.rm=T)
HisoprcpGH.slope.lo <- quantile(HisoprcpGH.slope, probs=0.025, na.rm=T)
HisoprcpGH.slope.up <- quantile(HisoprcpGH.slope, probs=0.975, na.rm=T)
HisoprcpGH.p.mn <- mean(HisoprcpGH.p, na.rm=T)
HisoprcpGH.p.lo <- quantile(HisoprcpGH.p, probs=0.025, na.rm=T)
HisoprcpGH.p.up <- quantile(HisoprcpGH.p, probs=0.975, na.rm=T)
HisoprcpGH.R2.mn <- mean(HisoprcpGH.R2, na.rm=T)
HisoprcpGH.R2.lo <- quantile(HisoprcpGH.R2, probs=0.025, na.rm=T)
HisoprcpGH.R2.up <- quantile(HisoprcpGH.R2, probs=0.975, na.rm=T)

HisoprcpH.slope.mn <- mean(HisoprcpH.slope, na.rm=T)
HisoprcpH.slope.lo <- quantile(HisoprcpH.slope, probs=0.025, na.rm=T)
HisoprcpH.slope.up <- quantile(HisoprcpH.slope, probs=0.975, na.rm=T)
HisoprcpH.p.mn <- mean(HisoprcpH.p, na.rm=T)
HisoprcpH.p.lo <- quantile(HisoprcpH.p, probs=0.025, na.rm=T)
HisoprcpH.p.up <- quantile(HisoprcpH.p, probs=0.975, na.rm=T)
HisoprcpH.R2.mn <- mean(HisoprcpH.R2, na.rm=T)
HisoprcpH.R2.lo <- quantile(HisoprcpH.R2, probs=0.025, na.rm=T)
HisoprcpH.R2.up <- quantile(HisoprcpH.R2, probs=0.975, na.rm=T)

Hisotree.slope.mn <- mean(Hisotree.slope, na.rm=T)
Hisotree.slope.lo <- quantile(Hisotree.slope, probs=0.025, na.rm=T)
Hisotree.slope.up <- quantile(Hisotree.slope, probs=0.975, na.rm=T)
Hisotree.p.mn <- mean(Hisotree.p, na.rm=T)
Hisotree.p.lo <- quantile(Hisotree.p, probs=0.025, na.rm=T)
Hisotree.p.up <- quantile(Hisotree.p, probs=0.975, na.rm=T)
Hisotree.R2.mn <- mean(Hisotree.R2, na.rm=T)
Hisotree.R2.lo <- quantile(Hisotree.R2, probs=0.025, na.rm=T)
Hisotree.R2.up <- quantile(Hisotree.R2, probs=0.975, na.rm=T)

insolprcpGcorL.slope.mn <- mean(insolprcpGcorL.slope, na.rm=T)
insolprcpGcorL.slope.lo <- quantile(insolprcpGcorL.slope, probs=0.025, na.rm=T)
insolprcpGcorL.slope.up <- quantile(insolprcpGcorL.slope, probs=0.975, na.rm=T)
insolprcpGcorL.p.mn <- mean(insolprcpGcorL.p, na.rm=T)
insolprcpGcorL.p.lo <- quantile(insolprcpGcorL.p, probs=0.025, na.rm=T)
insolprcpGcorL.p.up <- quantile(insolprcpGcorL.p, probs=0.975, na.rm=T)
insolprcpGcorL.R2.mn <- mean(insolprcpGcorL.R2, na.rm=T)
insolprcpGcorL.R2.lo <- quantile(insolprcpGcorL.R2, probs=0.025, na.rm=T)
insolprcpGcorL.R2.up <- quantile(insolprcpGcorL.R2, probs=0.975, na.rm=T)

insolprcpGL.slope.mn <- mean(insolprcpGL.slope, na.rm=T)
insolprcpGL.slope.lo <- quantile(insolprcpGL.slope, probs=0.025, na.rm=T)
insolprcpGL.slope.up <- quantile(insolprcpGL.slope, probs=0.975, na.rm=T)
insolprcpGL.p.mn <- mean(insolprcpGL.p, na.rm=T)
insolprcpGL.p.lo <- quantile(insolprcpGL.p, probs=0.025, na.rm=T)
insolprcpGL.p.up <- quantile(insolprcpGL.p, probs=0.975, na.rm=T)
insolprcpGL.R2.mn <- mean(insolprcpGL.R2, na.rm=T)
insolprcpGL.R2.lo <- quantile(insolprcpGL.R2, probs=0.025, na.rm=T)
insolprcpGL.R2.up <- quantile(insolprcpGL.R2, probs=0.975, na.rm=T)

insolprcpL.slope.mn <- mean(insolprcpL.slope, na.rm=T)
insolprcpL.slope.lo <- quantile(insolprcpL.slope, probs=0.025, na.rm=T)
insolprcpL.slope.up <- quantile(insolprcpL.slope, probs=0.975, na.rm=T)
insolprcpL.p.mn <- mean(insolprcpL.p, na.rm=T)
insolprcpL.p.lo <- quantile(insolprcpL.p, probs=0.025, na.rm=T)
insolprcpL.p.up <- quantile(insolprcpL.p, probs=0.975, na.rm=T)
insolprcpL.R2.mn <- mean(insolprcpL.R2, na.rm=T)
insolprcpL.R2.lo <- quantile(insolprcpL.R2, probs=0.025, na.rm=T)
insolprcpL.R2.up <- quantile(insolprcpL.R2, probs=0.975, na.rm=T)

insolprcpGcorH.slope.mn <- mean(insolprcpGcorH.slope, na.rm=T)
insolprcpGcorH.slope.lo <- quantile(insolprcpGcorH.slope, probs=0.025, na.rm=T)
insolprcpGcorH.slope.up <- quantile(insolprcpGcorH.slope, probs=0.975, na.rm=T)
insolprcpGcorH.p.mn <- mean(insolprcpGcorH.p, na.rm=T)
insolprcpGcorH.p.lo <- quantile(insolprcpGcorH.p, probs=0.025, na.rm=T)
insolprcpGcorH.p.up <- quantile(insolprcpGcorH.p, probs=0.975, na.rm=T)
insolprcpGcorH.R2.mn <- mean(insolprcpGcorH.R2, na.rm=T)
insolprcpGcorH.R2.lo <- quantile(insolprcpGcorH.R2, probs=0.025, na.rm=T)
insolprcpGcorH.R2.up <- quantile(insolprcpGcorH.R2, probs=0.975, na.rm=T)

insolprcpGH.slope.mn <- mean(insolprcpGH.slope, na.rm=T)
insolprcpGH.slope.lo <- quantile(insolprcpGH.slope, probs=0.025, na.rm=T)
insolprcpGH.slope.up <- quantile(insolprcpGH.slope, probs=0.975, na.rm=T)
insolprcpGH.p.mn <- mean(insolprcpGH.p, na.rm=T)
insolprcpGH.p.lo <- quantile(insolprcpGH.p, probs=0.025, na.rm=T)
insolprcpGH.p.up <- quantile(insolprcpGH.p, probs=0.975, na.rm=T)
insolprcpGH.R2.mn <- mean(insolprcpGH.R2, na.rm=T)
insolprcpGH.R2.lo <- quantile(insolprcpGH.R2, probs=0.025, na.rm=T)
insolprcpGH.R2.up <- quantile(insolprcpGH.R2, probs=0.975, na.rm=T)

insolprcpH.slope.mn <- mean(insolprcpH.slope, na.rm=T)
insolprcpH.slope.lo <- quantile(insolprcpH.slope, probs=0.025, na.rm=T)
insolprcpH.slope.up <- quantile(insolprcpH.slope, probs=0.975, na.rm=T)
insolprcpH.p.mn <- mean(insolprcpH.p, na.rm=T)
insolprcpH.p.lo <- quantile(insolprcpH.p, probs=0.025, na.rm=T)
insolprcpH.p.up <- quantile(insolprcpH.p, probs=0.975, na.rm=T)
insolprcpH.R2.mn <- mean(insolprcpH.R2, na.rm=T)
insolprcpH.R2.lo <- quantile(insolprcpH.R2, probs=0.025, na.rm=T)
insolprcpH.R2.up <- quantile(insolprcpH.R2, probs=0.975, na.rm=T)

insoltree.slope.mn <- mean(insoltree.slope, na.rm=T)
insoltree.slope.lo <- quantile(insoltree.slope, probs=0.025, na.rm=T)
insoltree.slope.up <- quantile(insoltree.slope, probs=0.975, na.rm=T)
insoltree.p.mn <- mean(insoltree.p, na.rm=T)
insoltree.p.lo <- quantile(insoltree.p, probs=0.025, na.rm=T)
insoltree.p.up <- quantile(insoltree.p, probs=0.975, na.rm=T)
insoltree.R2.mn <- mean(insoltree.R2, na.rm=T)
insoltree.R2.lo <- quantile(insoltree.R2, probs=0.025, na.rm=T)
insoltree.R2.up <- quantile(insoltree.R2, probs=0.975, na.rm=T)

treeprcpGL.slope.mn <- mean(treeprcpGL.slope, na.rm=T)
treeprcpGL.slope.lo <- quantile(treeprcpGL.slope, probs=0.025, na.rm=T)
treeprcpGL.slope.up <- quantile(treeprcpGL.slope, probs=0.975, na.rm=T)
treeprcpGL.p.mn <- mean(treeprcpGL.p, na.rm=T)
treeprcpGL.p.lo <- quantile(treeprcpGL.p, probs=0.025, na.rm=T)
treeprcpGL.p.up <- quantile(treeprcpGL.p, probs=0.975, na.rm=T)
treeprcpGL.R2.mn <- mean(treeprcpGL.R2, na.rm=T)
treeprcpGL.R2.lo <- quantile(treeprcpGL.R2, probs=0.025, na.rm=T)
treeprcpGL.R2.up <- quantile(treeprcpGL.R2, probs=0.975, na.rm=T)

treeprcpL.slope.mn <- mean(treeprcpL.slope, na.rm=T)
treeprcpL.slope.lo <- quantile(treeprcpL.slope, probs=0.025, na.rm=T)
treeprcpL.slope.up <- quantile(treeprcpL.slope, probs=0.975, na.rm=T)
treeprcpL.p.mn <- mean(treeprcpL.p, na.rm=T)
treeprcpL.p.lo <- quantile(treeprcpL.p, probs=0.025, na.rm=T)
treeprcpL.p.up <- quantile(treeprcpL.p, probs=0.975, na.rm=T)
treeprcpL.R2.mn <- mean(treeprcpL.R2, na.rm=T)
treeprcpL.R2.lo <- quantile(treeprcpL.R2, probs=0.025, na.rm=T)
treeprcpL.R2.up <- quantile(treeprcpL.R2, probs=0.975, na.rm=T)

treeprcpGcorH.slope.mn <- mean(treeprcpGcorH.slope, na.rm=T)
treeprcpGcorH.slope.lo <- quantile(treeprcpGcorH.slope, probs=0.025, na.rm=T)
treeprcpGcorH.slope.up <- quantile(treeprcpGcorH.slope, probs=0.975, na.rm=T)
treeprcpGcorH.p.mn <- mean(treeprcpGcorH.p, na.rm=T)
treeprcpGcorH.p.lo <- quantile(treeprcpGcorH.p, probs=0.025, na.rm=T)
treeprcpGcorH.p.up <- quantile(treeprcpGcorH.p, probs=0.975, na.rm=T)
treeprcpGcorH.R2.mn <- mean(treeprcpGcorH.R2, na.rm=T)
treeprcpGcorH.R2.lo <- quantile(treeprcpGcorH.R2, probs=0.025, na.rm=T)
treeprcpGcorH.R2.up <- quantile(treeprcpGcorH.R2, probs=0.975, na.rm=T)

treeprcpGcorL.slope.mn <- mean(treeprcpGcorL.slope, na.rm=T)
treeprcpGcorL.slope.lo <- quantile(treeprcpGcorL.slope, probs=0.025, na.rm=T)
treeprcpGcorL.slope.up <- quantile(treeprcpGcorL.slope, probs=0.975, na.rm=T)
treeprcpGcorL.p.mn <- mean(treeprcpGcorL.p, na.rm=T)
treeprcpGcorL.p.lo <- quantile(treeprcpGcorL.p, probs=0.025, na.rm=T)
treeprcpGcorL.p.up <- quantile(treeprcpGcorL.p, probs=0.975, na.rm=T)
treeprcpGcorL.R2.mn <- mean(treeprcpGcorL.R2, na.rm=T)
treeprcpGcorL.R2.lo <- quantile(treeprcpGcorL.R2, probs=0.025, na.rm=T)
treeprcpGcorL.R2.up <- quantile(treeprcpGcorL.R2, probs=0.975, na.rm=T)

treeprcpGH.slope.mn <- mean(treeprcpGH.slope, na.rm=T)
treeprcpGH.slope.lo <- quantile(treeprcpGH.slope, probs=0.025, na.rm=T)
treeprcpGH.slope.up <- quantile(treeprcpGH.slope, probs=0.975, na.rm=T)
treeprcpGH.p.mn <- mean(treeprcpGH.p, na.rm=T)
treeprcpGH.p.lo <- quantile(treeprcpGH.p, probs=0.025, na.rm=T)
treeprcpGH.p.up <- quantile(treeprcpGH.p, probs=0.975, na.rm=T)
treeprcpGH.R2.mn <- mean(treeprcpGH.R2, na.rm=T)
treeprcpGH.R2.lo <- quantile(treeprcpGH.R2, probs=0.025, na.rm=T)
treeprcpGH.R2.up <- quantile(treeprcpGH.R2, probs=0.975, na.rm=T)

treeprcpH.slope.mn <- mean(treeprcpH.slope, na.rm=T)
treeprcpH.slope.lo <- quantile(treeprcpH.slope, probs=0.025, na.rm=T)
treeprcpH.slope.up <- quantile(treeprcpH.slope, probs=0.975, na.rm=T)
treeprcpH.p.mn <- mean(treeprcpH.p, na.rm=T)
treeprcpH.p.lo <- quantile(treeprcpH.p, probs=0.025, na.rm=T)
treeprcpH.p.up <- quantile(treeprcpH.p, probs=0.975, na.rm=T)
treeprcpH.R2.mn <- mean(treeprcpH.R2, na.rm=T)
treeprcpH.R2.lo <- quantile(treeprcpH.R2, probs=0.025, na.rm=T)
treeprcpH.R2.up <- quantile(treeprcpH.R2, probs=0.975, na.rm=T)

# place lo/up into matrices following cor.mean format
p.lo <- cor.mean
p.lo[,] <- NA
p.lo[upper.tri(p.lo)] <- 0
p.lo[2,1] <- insolprcpGcorL.p.lo
p.lo[3,1] <- insolprcpGL.p.lo
p.lo[4,1] <- insolprcpGcorH.p.lo
p.lo[5,1] <- insolprcpGH.p.lo
p.lo[6,1] <- MaO18insol.p.lo
p.lo[7,1] <- doleinsol.p.lo
p.lo[8,1] <- insolprcpL.p.lo
p.lo[9,1] <- insolprcpH.p.lo
p.lo[11,1] <- Hisoinsol.p.lo
p.lo[12,1] <- insoltree.p.lo
p.lo[6,3] <- MaO18prcpGcorL.p.lo
p.lo[7,3] <- doleprcpGcorL.p.lo
p.lo[11,3] <- HisoprcpGcorL.p.lo
p.lo[12,3] <- treeprcpGcorL.p.lo;
p.lo[6,4] <- MaO18prcpGL.p.lo
p.lo[7,4] <- doleprcpGL.p.lo
p.lo[11,4] <- HisoprcpGL.p.lo
p.lo[12,4] <- treeprcpGL.p.lo
p.lo[6,5] <- MaO18prcpGcorH.p.lo
p.lo[7,5] <- doleprcpGcorH.p.lo
p.lo[11,5] <- HisoprcpGcorH.p.lo
p.lo[12,5] <- treeprcpGcorH.p.lo 
p.lo[6,6] <- MaO18prcpGH.p.lo
p.lo[7,6] <- doleprcpGH.p.lo
p.lo[11,6] <- HisoprcpGH.p.lo
p.lo[12,6] <- treeprcpGH.p.lo;
p.lo[7,7] <- MaO18dole.p.lo
p.lo[8,7] <- MaO18prcpL.p.lo
p.lo[9,7] <- MaO18prcpH.p.lo
p.lo[11,7] <- MaO18Hiso.p.lo
p.lo[12,7] <- MaO18tree.p.lo
p.lo[8,8] <- doleprcpL.p.lo
p.lo[9,8] <- doleprcpH.p.lo
p.lo[11,8] <- Hisodole.p.lo
p.lo[12,8] <- doletree.p.lo
p.lo[11,9] <- HisoprcpL.p.lo
p.lo[12,9] <- treeprcpL.p.lo
p.lo[11,10] <- HisoprcpH.p.lo
p.lo[12,10] <- treeprcpH.p.lo
p.lo[12,12] <- Hisotree.p.lo
round(p.lo, 3)

p.up <- cor.mean
p.up[,] <- NA
p.up[upper.tri(p.up)] <- 0
p.up[2,1] <- insolprcpGcorL.p.up
p.up[3,1] <- insolprcpGL.p.up
p.up[4,1] <- insolprcpGcorH.p.up
p.up[5,1] <- insolprcpGH.p.up
p.up[6,1] <- MaO18insol.p.up
p.up[7,1] <- doleinsol.p.up
p.up[8,1] <- insolprcpL.p.up
p.up[9,1] <- insolprcpH.p.up
p.up[11,1] <- Hisoinsol.p.up
p.up[12,1] <- insoltree.p.up
p.up[6,3] <- MaO18prcpGcorL.p.up
p.up[7,3] <- doleprcpGcorL.p.up
p.up[11,3] <- HisoprcpGcorL.p.up
p.up[12,3] <- treeprcpGcorL.p.up;
p.up[6,4] <- MaO18prcpGL.p.up
p.up[7,4] <- doleprcpGL.p.up
p.up[11,4] <- HisoprcpGL.p.up
p.up[12,4] <- treeprcpGL.p.up
p.up[6,5] <- MaO18prcpGcorH.p.up
p.up[7,5] <- doleprcpGcorH.p.up
p.up[11,5] <- HisoprcpGcorH.p.up
p.up[12,5] <- treeprcpGcorH.p.up 
p.up[6,6] <- MaO18prcpGH.p.up
p.up[7,6] <- doleprcpGH.p.up
p.up[11,6] <- HisoprcpGH.p.up
p.up[12,6] <- treeprcpGH.p.up;
p.up[7,7] <- MaO18dole.p.up
p.up[8,7] <- MaO18prcpL.p.up
p.up[9,7] <- MaO18prcpH.p.up
p.up[11,7] <- MaO18Hiso.p.up
p.up[12,7] <- MaO18tree.p.up
p.up[8,8] <- doleprcpL.p.up
p.up[9,8] <- doleprcpH.p.up
p.up[11,8] <- Hisodole.p.up
p.up[12,8] <- doletree.p.up
p.up[11,9] <- HisoprcpL.p.up
p.up[12,9] <- treeprcpL.p.up
p.up[11,10] <- HisoprcpH.p.up
p.up[12,10] <- treeprcpH.p.up
p.up[12,12] <- Hisotree.p.up
round(p.up, 3)

p.mn <- cor.mean
p.mn[,] <- NA
p.mn[upper.tri(p.mn)] <- 0
p.mn[2,1] <- insolprcpGcorL.p.mn
p.mn[3,1] <- insolprcpGL.p.mn
p.mn[4,1] <- insolprcpGcorH.p.mn
p.mn[5,1] <- insolprcpGH.p.mn
p.mn[6,1] <- MaO18insol.p.mn
p.mn[7,1] <- doleinsol.p.mn
p.mn[8,1] <- insolprcpL.p.mn
p.mn[9,1] <- insolprcpH.p.mn
p.mn[11,1] <- Hisoinsol.p.mn
p.mn[12,1] <- insoltree.p.mn
p.mn[6,3] <- MaO18prcpGcorL.p.mn
p.mn[7,3] <- doleprcpGcorL.p.mn
p.mn[11,3] <- HisoprcpGcorL.p.mn
p.mn[12,3] <- treeprcpGcorL.p.mn;
p.mn[6,4] <- MaO18prcpGL.p.mn
p.mn[7,4] <- doleprcpGL.p.mn
p.mn[11,4] <- HisoprcpGL.p.mn
p.mn[12,4] <- treeprcpGL.p.mn
p.mn[6,5] <- MaO18prcpGcorH.p.mn
p.mn[7,5] <- doleprcpGcorH.p.mn
p.mn[11,5] <- HisoprcpGcorH.p.mn
p.mn[12,5] <- treeprcpGcorH.p.mn 
p.mn[6,6] <- MaO18prcpGH.p.mn
p.mn[7,6] <- doleprcpGH.p.mn
p.mn[11,6] <- HisoprcpGH.p.mn
p.mn[12,6] <- treeprcpGH.p.mn;
p.mn[7,7] <- MaO18dole.p.mn
p.mn[8,7] <- MaO18prcpL.p.mn
p.mn[9,7] <- MaO18prcpH.p.mn
p.mn[11,7] <- MaO18Hiso.p.mn
p.mn[12,7] <- MaO18tree.p.mn
p.mn[8,8] <- doleprcpL.p.mn
p.mn[9,8] <- doleprcpH.p.mn
p.mn[11,8] <- Hisodole.p.mn
p.mn[12,8] <- doletree.p.mn
p.mn[11,9] <- HisoprcpL.p.mn
p.mn[12,9] <- treeprcpL.p.mn
p.mn[11,10] <- HisoprcpH.p.mn
p.mn[12,10] <- treeprcpH.p.mn
p.mn[12,12] <- Hisotree.p.mn
round(p.mn, 3)

R2.lo <- cor.mean
R2.lo[,] <- NA
R2.lo[upper.tri(R2.lo)] <- 0
R2.lo[2,1] <- insolprcpGcorL.R2.lo
R2.lo[3,1] <- insolprcpGL.R2.lo
R2.lo[4,1] <- insolprcpGcorH.R2.lo
R2.lo[5,1] <- insolprcpGH.R2.lo
R2.lo[6,1] <- MaO18insol.R2.lo
R2.lo[7,1] <- doleinsol.R2.lo
R2.lo[8,1] <- insolprcpL.R2.lo
R2.lo[9,1] <- insolprcpH.R2.lo
R2.lo[11,1] <- Hisoinsol.R2.lo
R2.lo[12,1] <- insoltree.R2.lo
R2.lo[6,3] <- MaO18prcpGcorL.R2.lo
R2.lo[7,3] <- doleprcpGcorL.R2.lo
R2.lo[11,3] <- HisoprcpGcorL.R2.lo
R2.lo[12,3] <- treeprcpGcorL.R2.lo;
R2.lo[6,4] <- MaO18prcpGL.R2.lo
R2.lo[7,4] <- doleprcpGL.R2.lo
R2.lo[11,4] <- HisoprcpGL.R2.lo
R2.lo[12,4] <- treeprcpGL.R2.lo
R2.lo[6,5] <- MaO18prcpGcorH.R2.lo
R2.lo[7,5] <- doleprcpGcorH.R2.lo
R2.lo[11,5] <- HisoprcpGcorH.R2.lo
R2.lo[12,5] <- treeprcpGcorH.R2.lo 
R2.lo[6,6] <- MaO18prcpGH.R2.lo
R2.lo[7,6] <- doleprcpGH.R2.lo
R2.lo[11,6] <- HisoprcpGH.R2.lo
R2.lo[12,6] <- treeprcpGH.R2.lo;
R2.lo[7,7] <- MaO18dole.R2.lo
R2.lo[8,7] <- MaO18prcpL.R2.lo
R2.lo[9,7] <- MaO18prcpH.R2.lo
R2.lo[11,7] <- MaO18Hiso.R2.lo
R2.lo[12,7] <- MaO18tree.R2.lo
R2.lo[8,8] <- doleprcpL.R2.lo
R2.lo[9,8] <- doleprcpH.R2.lo
R2.lo[11,8] <- Hisodole.R2.lo
R2.lo[12,8] <- doletree.R2.lo
R2.lo[11,9] <- HisoprcpL.R2.lo
R2.lo[12,9] <- treeprcpL.R2.lo
R2.lo[11,10] <- HisoprcpH.R2.lo
R2.lo[12,10] <- treeprcpH.R2.lo
R2.lo[12,12] <- Hisotree.R2.lo
round(R2.lo, 3)
# R2.lo.lowp <- R2.lo
# R2.lo.lowp[highp.ind] <- NA
# R2.lo.lowp

R2.up <- cor.mean
R2.up[,] <- NA
R2.up[upper.tri(R2.up)] <- 0
R2.up[2,1] <- insolprcpGcorL.R2.up
R2.up[3,1] <- insolprcpGL.R2.up
R2.up[4,1] <- insolprcpGcorH.R2.up
R2.up[5,1] <- insolprcpGH.R2.up
R2.up[6,1] <- MaO18insol.R2.up
R2.up[7,1] <- doleinsol.R2.up
R2.up[8,1] <- insolprcpL.R2.up
R2.up[9,1] <- insolprcpH.R2.up
R2.up[11,1] <- Hisoinsol.R2.up
R2.up[12,1] <- insoltree.R2.up
R2.up[6,3] <- MaO18prcpGcorL.R2.up
R2.up[7,3] <- doleprcpGcorL.R2.up
R2.up[11,3] <- HisoprcpGcorL.R2.up
R2.up[12,3] <- treeprcpGcorL.R2.up;
R2.up[6,4] <- MaO18prcpGL.R2.up
R2.up[7,4] <- doleprcpGL.R2.up
R2.up[11,4] <- HisoprcpGL.R2.up
R2.up[12,4] <- treeprcpGL.R2.up
R2.up[6,5] <- MaO18prcpGcorH.R2.up
R2.up[7,5] <- doleprcpGcorH.R2.up
R2.up[11,5] <- HisoprcpGcorH.R2.up
R2.up[12,5] <- treeprcpGcorH.R2.up 
R2.up[6,6] <- MaO18prcpGH.R2.up
R2.up[7,6] <- doleprcpGH.R2.up
R2.up[11,6] <- HisoprcpGH.R2.up
R2.up[12,6] <- treeprcpGH.R2.up;
R2.up[7,7] <- MaO18dole.R2.up
R2.up[8,7] <- MaO18prcpL.R2.up
R2.up[9,7] <- MaO18prcpH.R2.up
R2.up[11,7] <- MaO18Hiso.R2.up
R2.up[12,7] <- MaO18tree.R2.up
R2.up[8,8] <- doleprcpL.R2.up
R2.up[9,8] <- doleprcpH.R2.up
R2.up[11,8] <- Hisodole.R2.up
R2.up[12,8] <- doletree.R2.up
R2.up[11,9] <- HisoprcpL.R2.up
R2.up[12,9] <- treeprcpL.R2.up
R2.up[11,10] <- HisoprcpH.R2.up
R2.up[12,10] <- treeprcpH.R2.up
R2.up[12,12] <- Hisotree.p.lo
round(R2.up, 3)

R2.mn <- cor.mean
R2.mn[,] <- NA
R2.mn[upper.tri(R2.rs.mn)] <- 0
R2.mn[2,1] <- insolprcpGcorL.R2.mn
R2.mn[3,1] <- insolprcpGL.R2.mn
R2.mn[4,1] <- insolprcpGcorH.R2.mn
R2.mn[5,1] <- insolprcpGH.R2.mn
R2.mn[6,1] <- MaO18insol.R2.mn
R2.mn[7,1] <- doleinsol.R2.mn
R2.mn[8,1] <- insolprcpL.R2.mn
R2.mn[9,1] <- insolprcpH.R2.mn
R2.mn[11,1] <- Hisoinsol.R2.mn
R2.mn[12,1] <- insoltree.R2.mn
R2.mn[6,3] <- MaO18prcpGcorL.R2.mn
R2.mn[7,3] <- doleprcpGcorL.R2.mn
R2.mn[11,3] <- HisoprcpGcorL.R2.mn
R2.mn[12,3] <- treeprcpGcorL.R2.mn;
R2.mn[6,4] <- MaO18prcpGL.R2.mn
R2.mn[7,4] <- doleprcpGL.R2.mn
R2.mn[11,4] <- HisoprcpGL.R2.mn
R2.mn[12,4] <- treeprcpGL.R2.mn
R2.mn[6,5] <- MaO18prcpGcorH.R2.mn
R2.mn[7,5] <- doleprcpGcorH.R2.mn
R2.mn[11,5] <- HisoprcpGcorH.R2.mn
R2.mn[12,5] <- treeprcpGcorH.R2.mn 
R2.mn[6,6] <- MaO18prcpGH.R2.mn
R2.mn[7,6] <- doleprcpGH.R2.mn
R2.mn[11,6] <- HisoprcpGH.R2.mn
R2.mn[12,6] <- treeprcpGH.R2.mn;
R2.mn[7,7] <- MaO18dole.R2.mn
R2.mn[8,7] <- MaO18prcpL.R2.mn
R2.mn[9,7] <- MaO18prcpH.R2.mn
R2.mn[11,7] <- MaO18Hiso.R2.mn
R2.mn[12,7] <- MaO18tree.R2.mn
R2.mn[8,8] <- doleprcpL.R2.mn
R2.mn[9,8] <- doleprcpH.R2.mn
R2.mn[11,8] <- Hisodole.R2.mn
R2.mn[12,8] <- doletree.R2.mn
R2.mn[11,9] <- HisoprcpL.R2.mn
R2.mn[12,9] <- treeprcpL.R2.mn
R2.mn[11,10] <- HisoprcpH.R2.mn
R2.mn[12,10] <- treeprcpH.R2.mn
R2.mn[12,12] <- Hisotree.p.lo
round(R2.mn, 3)

slope.lo <- cor.mean
slope.lo[,] <- NA
slope.lo[upper.tri(slope.lo)] <- 0
slope.lo[2,1] <- insolprcpGcorL.slope.lo
slope.lo[3,1] <- insolprcpGL.slope.lo
slope.lo[4,1] <- insolprcpGcorH.slope.lo
slope.lo[5,1] <- insolprcpGH.slope.lo
slope.lo[6,1] <- MaO18insol.slope.lo
slope.lo[7,1] <- doleinsol.slope.lo
slope.lo[8,1] <- insolprcpL.slope.lo
slope.lo[9,1] <- insolprcpH.slope.lo
slope.lo[11,1] <- Hisoinsol.slope.lo
slope.lo[12,1] <- insoltree.slope.lo
slope.lo[6,3] <- MaO18prcpGcorL.slope.lo
slope.lo[7,3] <- doleprcpGcorL.slope.lo
slope.lo[11,3] <- HisoprcpGcorL.slope.lo
slope.lo[12,3] <- treeprcpGcorL.slope.lo;
slope.lo[6,4] <- MaO18prcpGL.slope.lo
slope.lo[7,4] <- doleprcpGL.slope.lo
slope.lo[11,4] <- HisoprcpGL.slope.lo
slope.lo[12,4] <- treeprcpGL.slope.lo
slope.lo[6,5] <- MaO18prcpGcorH.slope.lo
slope.lo[7,5] <- doleprcpGcorH.slope.lo
slope.lo[11,5] <- HisoprcpGcorH.slope.lo
slope.lo[12,5] <- treeprcpGcorH.slope.lo 
slope.lo[6,6] <- MaO18prcpGH.slope.lo
slope.lo[7,6] <- doleprcpGH.slope.lo
slope.lo[11,6] <- HisoprcpGH.slope.lo
slope.lo[12,6] <- treeprcpGH.slope.lo;
slope.lo[7,7] <- MaO18dole.slope.lo
slope.lo[8,7] <- MaO18prcpL.slope.lo
slope.lo[9,7] <- MaO18prcpH.slope.lo
slope.lo[11,7] <- MaO18Hiso.slope.lo
slope.lo[12,7] <- MaO18tree.slope.lo
slope.lo[8,8] <- doleprcpL.slope.lo
slope.lo[9,8] <- doleprcpH.slope.lo
slope.lo[11,8] <- Hisodole.slope.lo
slope.lo[12,8] <- doletree.slope.lo
slope.lo[11,9] <- HisoprcpL.slope.lo
slope.lo[12,9] <- treeprcpL.slope.lo
slope.lo[11,10] <- HisoprcpH.slope.lo
slope.lo[12,10] <- treeprcpH.slope.lo
slope.lo[12,12] <- Hisotree.slope.lo
round(slope.lo, 3)

slope.up <- cor.mean
slope.up[,] <- NA
slope.up[upper.tri(slope.up)] <- 0
slope.up[2,1] <- insolprcpGcorL.slope.up
slope.up[3,1] <- insolprcpGL.slope.up
slope.up[4,1] <- insolprcpGcorH.slope.up
slope.up[5,1] <- insolprcpGH.slope.up
slope.up[6,1] <- MaO18insol.slope.up
slope.up[7,1] <- doleinsol.slope.up
slope.up[8,1] <- insolprcpL.slope.up
slope.up[9,1] <- insolprcpH.slope.up
slope.up[11,1] <- Hisoinsol.slope.up
slope.up[12,1] <- insoltree.slope.up
slope.up[6,3] <- MaO18prcpGcorL.slope.up
slope.up[7,3] <- doleprcpGcorL.slope.up
slope.up[11,3] <- HisoprcpGcorL.slope.up
slope.up[12,3] <- treeprcpGcorL.slope.up;
slope.up[6,4] <- MaO18prcpGL.slope.up
slope.up[7,4] <- doleprcpGL.slope.up
slope.up[11,4] <- HisoprcpGL.slope.up
slope.up[12,4] <- treeprcpGL.slope.up
slope.up[6,5] <- MaO18prcpGcorH.slope.up
slope.up[7,5] <- doleprcpGcorH.slope.up
slope.up[11,5] <- HisoprcpGcorH.slope.up
slope.up[12,5] <- treeprcpGcorH.slope.up 
slope.up[6,6] <- MaO18prcpGH.slope.up
slope.up[7,6] <- doleprcpGH.slope.up
slope.up[11,6] <- HisoprcpGH.slope.up
slope.up[12,6] <- treeprcpGH.slope.up;
slope.up[7,7] <- MaO18dole.slope.up
slope.up[8,7] <- MaO18prcpL.slope.up
slope.up[9,7] <- MaO18prcpH.slope.up
slope.up[11,7] <- MaO18Hiso.slope.up
slope.up[12,7] <- MaO18tree.slope.up
slope.up[8,8] <- doleprcpL.slope.up
slope.up[9,8] <- doleprcpH.slope.up
slope.up[11,8] <- Hisodole.slope.up
slope.up[12,8] <- doletree.slope.up
slope.up[11,9] <- HisoprcpL.slope.up
slope.up[12,9] <- treeprcpL.slope.up
slope.up[11,10] <- HisoprcpH.slope.up
slope.up[12,10] <- treeprcpH.slope.up
slope.up[12,12] <- Hisotree.p.lo
round(slope.up, 3)



# RESAMPLED CORRELATION OUTPUT SUMMARIES
MaO18insol.slope.rs.mn <- mean(MaO18insol.slope.rs, na.rm=T)
MaO18insol.slope.rs.lo <- quantile(MaO18insol.slope.rs, probs=0.025, na.rm=T)
MaO18insol.slope.rs.up <- quantile(MaO18insol.slope.rs, probs=0.975, na.rm=T)
MaO18insol.p.rs.mn <- mean(MaO18insol.p.rs, na.rm=T)
MaO18insol.p.rs.lo <- quantile(MaO18insol.p.rs, probs=0.025, na.rm=T)
MaO18insol.p.rs.up <- quantile(MaO18insol.p.rs, probs=0.975, na.rm=T)
MaO18insol.R2.rs.mn <- mean(MaO18insol.R2.rs, na.rm=T)
MaO18insol.R2.rs.lo <- quantile(MaO18insol.R2.rs, probs=0.025, na.rm=T)
MaO18insol.R2.rs.up <- quantile(MaO18insol.R2.rs, probs=0.975, na.rm=T)
MaO18insol.DWp.rs.mn <- mean(MaO18insol.DWp.rs, na.rm=T)
MaO18insol.DWp.rs.lo <- quantile(MaO18insol.DWp.rs, probs=0.025, na.rm=T)
MaO18insol.DWp.rs.up <- quantile(MaO18insol.DWp.rs, probs=0.975, na.rm=T)

MaO18prcpGcorL.slope.rs.mn <- mean(MaO18prcpGcorL.slope.rs, na.rm=T)
MaO18prcpGcorL.slope.rs.lo <- quantile(MaO18prcpGcorL.slope.rs, probs=0.025, na.rm=T)
MaO18prcpGcorL.slope.rs.up <- quantile(MaO18prcpGcorL.slope.rs, probs=0.975, na.rm=T)
MaO18prcpGcorL.p.rs.mn <- mean(MaO18prcpGcorL.p.rs, na.rm=T)
MaO18prcpGcorL.p.rs.lo <- quantile(MaO18prcpGcorL.p.rs, probs=0.025, na.rm=T)
MaO18prcpGcorL.p.rs.up <- quantile(MaO18prcpGcorL.p.rs, probs=0.975, na.rm=T)
MaO18prcpGcorL.R2.rs.mn <- mean(MaO18prcpGcorL.R2.rs, na.rm=T)
MaO18prcpGcorL.R2.rs.lo <- quantile(MaO18prcpGcorL.R2.rs, probs=0.025, na.rm=T)
MaO18prcpGcorL.R2.rs.up <- quantile(MaO18prcpGcorL.R2.rs, probs=0.975, na.rm=T)
MaO18prcpGcorL.DWp.rs.mn <- mean(MaO18prcpGcorL.DWp.rs, na.rm=T)
MaO18prcpGcorL.DWp.rs.lo <- quantile(MaO18prcpGcorL.DWp.rs, probs=0.025, na.rm=T)
MaO18prcpGcorL.DWp.rs.up <- quantile(MaO18prcpGcorL.DWp.rs, probs=0.975, na.rm=T)

MaO18prcpGL.slope.rs.mn <- mean(MaO18prcpGL.slope.rs, na.rm=T)
MaO18prcpGL.slope.rs.lo <- quantile(MaO18prcpGL.slope.rs, probs=0.025, na.rm=T)
MaO18prcpGL.slope.rs.up <- quantile(MaO18prcpGL.slope.rs, probs=0.975, na.rm=T)
MaO18prcpGL.p.rs.mn <- mean(MaO18prcpGL.p.rs, na.rm=T)
MaO18prcpGL.p.rs.lo <- quantile(MaO18prcpGL.p.rs, probs=0.025, na.rm=T)
MaO18prcpGL.p.rs.up <- quantile(MaO18prcpGL.p.rs, probs=0.975, na.rm=T)
MaO18prcpGL.R2.rs.mn <- mean(MaO18prcpGL.R2.rs, na.rm=T)
MaO18prcpGL.R2.rs.lo <- quantile(MaO18prcpGL.R2.rs, probs=0.025, na.rm=T)
MaO18prcpGL.R2.rs.up <- quantile(MaO18prcpGL.R2.rs, probs=0.975, na.rm=T)
MaO18prcpGL.DWp.rs.mn <- mean(MaO18prcpGL.DWp.rs, na.rm=T)
MaO18prcpGL.DWp.rs.lo <- quantile(MaO18prcpGL.DWp.rs, probs=0.025, na.rm=T)
MaO18prcpGL.DWp.rs.up <- quantile(MaO18prcpGL.DWp.rs, probs=0.975, na.rm=T)

MaO18prcpL.slope.rs.mn <- mean(MaO18prcpL.slope.rs, na.rm=T)
MaO18prcpL.slope.rs.lo <- quantile(MaO18prcpL.slope.rs, probs=0.025, na.rm=T)
MaO18prcpL.slope.rs.up <- quantile(MaO18prcpL.slope.rs, probs=0.975, na.rm=T)
MaO18prcpL.p.rs.mn <- mean(MaO18prcpL.p.rs, na.rm=T)
MaO18prcpL.p.rs.lo <- quantile(MaO18prcpL.p.rs, probs=0.025, na.rm=T)
MaO18prcpL.p.rs.up <- quantile(MaO18prcpL.p.rs, probs=0.975, na.rm=T)
MaO18prcpL.R2.rs.mn <- mean(MaO18prcpL.R2.rs, na.rm=T)
MaO18prcpL.R2.rs.lo <- quantile(MaO18prcpL.R2.rs, probs=0.025, na.rm=T)
MaO18prcpL.R2.rs.up <- quantile(MaO18prcpL.R2.rs, probs=0.975, na.rm=T)
MaO18prcpL.DWp.rs.mn <- mean(MaO18prcpL.DWp.rs, na.rm=T)
MaO18prcpL.DWp.rs.lo <- quantile(MaO18prcpL.DWp.rs, probs=0.025, na.rm=T)
MaO18prcpL.DWp.rs.up <- quantile(MaO18prcpL.DWp.rs, probs=0.975, na.rm=T)

MaO18prcpGcorH.slope.rs.mn <- mean(MaO18prcpGcorH.slope.rs, na.rm=T)
MaO18prcpGcorH.slope.rs.lo <- quantile(MaO18prcpGcorH.slope.rs, probs=0.025, na.rm=T)
MaO18prcpGcorH.slope.rs.up <- quantile(MaO18prcpGcorH.slope.rs, probs=0.975, na.rm=T)
MaO18prcpGcorH.p.rs.mn <- mean(MaO18prcpGcorH.p.rs, na.rm=T)
MaO18prcpGcorH.p.rs.lo <- quantile(MaO18prcpGcorH.p.rs, probs=0.025, na.rm=T)
MaO18prcpGcorH.p.rs.up <- quantile(MaO18prcpGcorH.p.rs, probs=0.975, na.rm=T)
MaO18prcpGcorH.R2.rs.mn <- mean(MaO18prcpGcorH.R2.rs, na.rm=T)
MaO18prcpGcorH.R2.rs.lo <- quantile(MaO18prcpGcorH.R2.rs, probs=0.025, na.rm=T)
MaO18prcpGcorH.R2.rs.up <- quantile(MaO18prcpGcorH.R2.rs, probs=0.975, na.rm=T)
MaO18prcpGcorH.DWp.rs.mn <- mean(MaO18prcpGcorH.DWp.rs, na.rm=T)
MaO18prcpGcorH.DWp.rs.lo <- quantile(MaO18prcpGcorH.DWp.rs, probs=0.025, na.rm=T)
MaO18prcpGcorH.DWp.rs.up <- quantile(MaO18prcpGcorH.DWp.rs, probs=0.975, na.rm=T)

MaO18prcpGH.slope.rs.mn <- mean(MaO18prcpGH.slope.rs, na.rm=T)
MaO18prcpGH.slope.rs.lo <- quantile(MaO18prcpGH.slope.rs, probs=0.025, na.rm=T)
MaO18prcpGH.slope.rs.up <- quantile(MaO18prcpGH.slope.rs, probs=0.975, na.rm=T)
MaO18prcpGH.p.rs.mn <- mean(MaO18prcpGH.p.rs, na.rm=T)
MaO18prcpGH.p.rs.lo <- quantile(MaO18prcpGH.p.rs, probs=0.025, na.rm=T)
MaO18prcpGH.p.rs.up <- quantile(MaO18prcpGH.p.rs, probs=0.975, na.rm=T)
MaO18prcpGH.R2.rs.mn <- mean(MaO18prcpGH.R2.rs, na.rm=T)
MaO18prcpGH.R2.rs.lo <- quantile(MaO18prcpGH.R2.rs, probs=0.025, na.rm=T)
MaO18prcpGH.R2.rs.up <- quantile(MaO18prcpGH.R2.rs, probs=0.975, na.rm=T)
MaO18prcpGH.DWp.rs.mn <- mean(MaO18prcpGH.DWp.rs, na.rm=T)
MaO18prcpGH.DWp.rs.lo <- quantile(MaO18prcpGH.DWp.rs, probs=0.025, na.rm=T)
MaO18prcpGH.DWp.rs.up <- quantile(MaO18prcpGH.DWp.rs, probs=0.975, na.rm=T)

MaO18prcpH.slope.rs.mn <- mean(MaO18prcpH.slope.rs, na.rm=T)
MaO18prcpH.slope.rs.lo <- quantile(MaO18prcpH.slope.rs, probs=0.025, na.rm=T)
MaO18prcpH.slope.rs.up <- quantile(MaO18prcpH.slope.rs, probs=0.975, na.rm=T)
MaO18prcpH.p.rs.mn <- mean(MaO18prcpH.p.rs, na.rm=T)
MaO18prcpH.p.rs.lo <- quantile(MaO18prcpH.p.rs, probs=0.025, na.rm=T)
MaO18prcpH.p.rs.up <- quantile(MaO18prcpH.p.rs, probs=0.975, na.rm=T)
MaO18prcpH.R2.rs.mn <- mean(MaO18prcpH.R2.rs, na.rm=T)
MaO18prcpH.R2.rs.lo <- quantile(MaO18prcpH.R2.rs, probs=0.025, na.rm=T)
MaO18prcpH.R2.rs.up <- quantile(MaO18prcpH.R2.rs, probs=0.975, na.rm=T)
MaO18prcpH.DWp.rs.mn <- mean(MaO18prcpH.DWp.rs, na.rm=T)
MaO18prcpH.DWp.rs.lo <- quantile(MaO18prcpH.DWp.rs, probs=0.025, na.rm=T)
MaO18prcpH.DWp.rs.up <- quantile(MaO18prcpH.DWp.rs, probs=0.975, na.rm=T)

MaO18dole.slope.rs.mn <- mean(MaO18dole.slope.rs, na.rm=T)
MaO18dole.slope.rs.lo <- quantile(MaO18dole.slope.rs, probs=0.025, na.rm=T)
MaO18dole.slope.rs.up <- quantile(MaO18dole.slope.rs, probs=0.975, na.rm=T)
MaO18dole.p.rs.mn <- mean(MaO18dole.p.rs, na.rm=T)
MaO18dole.p.rs.lo <- quantile(MaO18dole.p.rs, probs=0.025, na.rm=T)
MaO18dole.p.rs.up <- quantile(MaO18dole.p.rs, probs=0.975, na.rm=T)
MaO18dole.R2.rs.mn <- mean(MaO18dole.R2.rs, na.rm=T)
MaO18dole.R2.rs.lo <- quantile(MaO18dole.R2.rs, probs=0.025, na.rm=T)
MaO18dole.R2.rs.up <- quantile(MaO18dole.R2.rs, probs=0.975, na.rm=T)
MaO18dole.DWp.rs.mn <- mean(MaO18dole.DWp.rs, na.rm=T)
MaO18dole.DWp.rs.lo <- quantile(MaO18dole.DWp.rs, probs=0.025, na.rm=T)
MaO18dole.DWp.rs.up <- quantile(MaO18dole.DWp.rs, probs=0.975, na.rm=T)

MaO18tree.slope.rs.mn <- mean(MaO18tree.slope.rs, na.rm=T)
MaO18tree.slope.rs.lo <- quantile(MaO18tree.slope.rs, probs=0.025, na.rm=T)
MaO18tree.slope.rs.up <- quantile(MaO18tree.slope.rs, probs=0.975, na.rm=T)
MaO18tree.p.rs.mn <- mean(MaO18tree.p.rs, na.rm=T)
MaO18tree.p.rs.lo <- quantile(MaO18tree.p.rs, probs=0.025, na.rm=T)
MaO18tree.p.rs.up <- quantile(MaO18tree.p.rs, probs=0.975, na.rm=T)
MaO18tree.R2.rs.mn <- mean(MaO18tree.R2.rs, na.rm=T)
MaO18tree.R2.rs.lo <- quantile(MaO18tree.R2.rs, probs=0.025, na.rm=T)
MaO18tree.R2.rs.up <- quantile(MaO18tree.R2.rs, probs=0.975, na.rm=T)
MaO18tree.DWp.rs.mn <- mean(MaO18tree.DWp.rs, na.rm=T)
MaO18tree.DWp.rs.lo <- quantile(MaO18tree.DWp.rs, probs=0.025, na.rm=T)
MaO18tree.DWp.rs.up <- quantile(MaO18tree.DWp.rs, probs=0.975, na.rm=T)

MaO18Hiso.slope.rs.mn <- mean(MaO18Hiso.slope.rs, na.rm=T)
MaO18Hiso.slope.rs.lo <- quantile(MaO18Hiso.slope.rs, probs=0.025, na.rm=T)
MaO18Hiso.slope.rs.up <- quantile(MaO18Hiso.slope.rs, probs=0.975, na.rm=T)
MaO18Hiso.p.rs.mn <- mean(MaO18Hiso.p.rs, na.rm=T)
MaO18Hiso.p.rs.lo <- quantile(MaO18Hiso.p.rs, probs=0.025, na.rm=T)
MaO18Hiso.p.rs.up <- quantile(MaO18Hiso.p.rs, probs=0.975, na.rm=T)
MaO18Hiso.R2.rs.mn <- mean(MaO18Hiso.R2.rs, na.rm=T)
MaO18Hiso.R2.rs.lo <- quantile(MaO18Hiso.R2.rs, probs=0.025, na.rm=T)
MaO18Hiso.R2.rs.up <- quantile(MaO18Hiso.R2.rs, probs=0.975, na.rm=T)
MaO18Hiso.DWp.rs.mn <- mean(MaO18Hiso.DWp.rs, na.rm=T)
MaO18Hiso.DWp.rs.lo <- quantile(MaO18Hiso.DWp.rs, probs=0.025, na.rm=T)
MaO18Hiso.DWp.rs.up <- quantile(MaO18Hiso.DWp.rs, probs=0.975, na.rm=T)


doleinsol.slope.rs.mn <- mean(doleinsol.slope.rs, na.rm=T)
doleinsol.slope.rs.lo <- quantile(doleinsol.slope.rs, probs=0.025, na.rm=T)
doleinsol.slope.rs.up <- quantile(doleinsol.slope.rs, probs=0.975, na.rm=T)
doleinsol.p.rs.mn <- mean(doleinsol.p.rs, na.rm=T)
doleinsol.p.rs.lo <- quantile(doleinsol.p.rs, probs=0.025, na.rm=T)
doleinsol.p.rs.up <- quantile(doleinsol.p.rs, probs=0.975, na.rm=T)
doleinsol.R2.rs.mn <- mean(doleinsol.R2.rs, na.rm=T)
doleinsol.R2.rs.lo <- quantile(doleinsol.R2.rs, probs=0.025, na.rm=T)
doleinsol.R2.rs.up <- quantile(doleinsol.R2.rs, probs=0.975, na.rm=T)
doleinsol.DWp.rs.mn <- mean(doleinsol.DWp.rs, na.rm=T)
doleinsol.DWp.rs.lo <- quantile(doleinsol.DWp.rs, probs=0.025, na.rm=T)
doleinsol.DWp.rs.up <- quantile(doleinsol.DWp.rs, probs=0.975, na.rm=T)

doleprcpGcorL.slope.rs.mn <- mean(doleprcpGcorL.slope.rs, na.rm=T)
doleprcpGcorL.slope.rs.lo <- quantile(doleprcpGcorL.slope.rs, probs=0.025, na.rm=T)
doleprcpGcorL.slope.rs.up <- quantile(doleprcpGcorL.slope.rs, probs=0.975, na.rm=T)
doleprcpGcorL.p.rs.mn <- mean(doleprcpGcorL.p.rs, na.rm=T)
doleprcpGcorL.p.rs.lo <- quantile(doleprcpGcorL.p.rs, probs=0.025, na.rm=T)
doleprcpGcorL.p.rs.up <- quantile(doleprcpGcorL.p.rs, probs=0.975, na.rm=T)
doleprcpGcorL.R2.rs.mn <- mean(doleprcpGcorL.R2.rs, na.rm=T)
doleprcpGcorL.R2.rs.lo <- quantile(doleprcpGcorL.R2.rs, probs=0.025, na.rm=T)
doleprcpGcorL.R2.rs.up <- quantile(doleprcpGcorL.R2.rs, probs=0.975, na.rm=T)
doleprcpGcorL.DWp.rs.mn <- mean(doleprcpGcorL.DWp.rs, na.rm=T)
doleprcpGcorL.DWp.rs.lo <- quantile(doleprcpGcorL.DWp.rs, probs=0.025, na.rm=T)
doleprcpGcorL.DWp.rs.up <- quantile(doleprcpGcorL.DWp.rs, probs=0.975, na.rm=T)

doleprcpGL.slope.rs.mn <- mean(doleprcpGL.slope.rs, na.rm=T)
doleprcpGL.slope.rs.lo <- quantile(doleprcpGL.slope.rs, probs=0.025, na.rm=T)
doleprcpGL.slope.rs.up <- quantile(doleprcpGL.slope.rs, probs=0.975, na.rm=T)
doleprcpGL.p.rs.mn <- mean(doleprcpGL.p.rs, na.rm=T)
doleprcpGL.p.rs.lo <- quantile(doleprcpGL.p.rs, probs=0.025, na.rm=T)
doleprcpGL.p.rs.up <- quantile(doleprcpGL.p.rs, probs=0.975, na.rm=T)
doleprcpGL.R2.rs.mn <- mean(doleprcpGL.R2.rs, na.rm=T)
doleprcpGL.R2.rs.lo <- quantile(doleprcpGL.R2.rs, probs=0.025, na.rm=T)
doleprcpGL.R2.rs.up <- quantile(doleprcpGL.R2.rs, probs=0.975, na.rm=T)
doleprcpGL.DWp.rs.mn <- mean(doleprcpGL.DWp.rs, na.rm=T)
doleprcpGL.DWp.rs.lo <- quantile(doleprcpGL.DWp.rs, probs=0.025, na.rm=T)
doleprcpGL.DWp.rs.up <- quantile(doleprcpGL.DWp.rs, probs=0.975, na.rm=T)

doleprcpL.slope.rs.mn <- mean(doleprcpL.slope.rs, na.rm=T)
doleprcpL.slope.rs.lo <- quantile(doleprcpL.slope.rs, probs=0.025, na.rm=T)
doleprcpL.slope.rs.up <- quantile(doleprcpL.slope.rs, probs=0.975, na.rm=T)
doleprcpL.p.rs.mn <- mean(doleprcpL.p.rs, na.rm=T)
doleprcpL.p.rs.lo <- quantile(doleprcpL.p.rs, probs=0.025, na.rm=T)
doleprcpL.p.rs.up <- quantile(doleprcpL.p.rs, probs=0.975, na.rm=T)
doleprcpL.R2.rs.mn <- mean(doleprcpL.R2.rs, na.rm=T)
doleprcpL.R2.rs.lo <- quantile(doleprcpL.R2.rs, probs=0.025, na.rm=T)
doleprcpL.R2.rs.up <- quantile(doleprcpL.R2.rs, probs=0.975, na.rm=T)
doleprcpL.DWp.rs.mn <- mean(doleprcpL.DWp.rs, na.rm=T)
doleprcpL.DWp.rs.lo <- quantile(doleprcpL.DWp.rs, probs=0.025, na.rm=T)
doleprcpL.DWp.rs.up <- quantile(doleprcpL.DWp.rs, probs=0.975, na.rm=T)

doleprcpGcorH.slope.rs.mn <- mean(doleprcpGcorH.slope.rs, na.rm=T)
doleprcpGcorH.slope.rs.lo <- quantile(doleprcpGcorH.slope.rs, probs=0.025, na.rm=T)
doleprcpGcorH.slope.rs.up <- quantile(doleprcpGcorH.slope.rs, probs=0.975, na.rm=T)
doleprcpGcorH.p.rs.mn <- mean(doleprcpGcorH.p.rs, na.rm=T)
doleprcpGcorH.p.rs.lo <- quantile(doleprcpGcorH.p.rs, probs=0.025, na.rm=T)
doleprcpGcorH.p.rs.up <- quantile(doleprcpGcorH.p.rs, probs=0.975, na.rm=T)
doleprcpGcorH.R2.rs.mn <- mean(doleprcpGcorH.R2.rs, na.rm=T)
doleprcpGcorH.R2.rs.lo <- quantile(doleprcpGcorH.R2.rs, probs=0.025, na.rm=T)
doleprcpGcorH.R2.rs.up <- quantile(doleprcpGcorH.R2.rs, probs=0.975, na.rm=T)
doleprcpGcorH.DWp.rs.mn <- mean(doleprcpGcorH.DWp.rs, na.rm=T)
doleprcpGcorH.DWp.rs.lo <- quantile(doleprcpGcorH.DWp.rs, probs=0.025, na.rm=T)
doleprcpGcorH.DWp.rs.up <- quantile(doleprcpGcorH.DWp.rs, probs=0.975, na.rm=T)

doleprcpGH.slope.rs.mn <- mean(doleprcpGH.slope.rs, na.rm=T)
doleprcpGH.slope.rs.lo <- quantile(doleprcpGH.slope.rs, probs=0.025, na.rm=T)
doleprcpGH.slope.rs.up <- quantile(doleprcpGH.slope.rs, probs=0.975, na.rm=T)
doleprcpGH.p.rs.mn <- mean(doleprcpGH.p.rs, na.rm=T)
doleprcpGH.p.rs.lo <- quantile(doleprcpGH.p.rs, probs=0.025, na.rm=T)
doleprcpGH.p.rs.up <- quantile(doleprcpGH.p.rs, probs=0.975, na.rm=T)
doleprcpGH.R2.rs.mn <- mean(doleprcpGH.R2.rs, na.rm=T)
doleprcpGH.R2.rs.lo <- quantile(doleprcpGH.R2.rs, probs=0.025, na.rm=T)
doleprcpGH.R2.rs.up <- quantile(doleprcpGH.R2.rs, probs=0.975, na.rm=T)
doleprcpGH.DWp.rs.mn <- mean(doleprcpGH.DWp.rs, na.rm=T)
doleprcpGH.DWp.rs.lo <- quantile(doleprcpGH.DWp.rs, probs=0.025, na.rm=T)
doleprcpGH.DWp.rs.up <- quantile(doleprcpGH.DWp.rs, probs=0.975, na.rm=T)

doleprcpH.slope.rs.mn <- mean(doleprcpH.slope.rs, na.rm=T)
doleprcpH.slope.rs.lo <- quantile(doleprcpH.slope.rs, probs=0.025, na.rm=T)
doleprcpH.slope.rs.up <- quantile(doleprcpH.slope.rs, probs=0.975, na.rm=T)
doleprcpH.p.rs.mn <- mean(doleprcpH.p.rs, na.rm=T)
doleprcpH.p.rs.lo <- quantile(doleprcpH.p.rs, probs=0.025, na.rm=T)
doleprcpH.p.rs.up <- quantile(doleprcpH.p.rs, probs=0.975, na.rm=T)
doleprcpH.R2.rs.mn <- mean(doleprcpH.R2.rs, na.rm=T)
doleprcpH.R2.rs.lo <- quantile(doleprcpH.R2.rs, probs=0.025, na.rm=T)
doleprcpH.R2.rs.up <- quantile(doleprcpH.R2.rs, probs=0.975, na.rm=T)
doleprcpH.DWp.rs.mn <- mean(doleprcpH.DWp.rs, na.rm=T)
doleprcpH.DWp.rs.lo <- quantile(doleprcpH.DWp.rs, probs=0.025, na.rm=T)
doleprcpH.DWp.rs.up <- quantile(doleprcpH.DWp.rs, probs=0.975, na.rm=T)

doletree.slope.rs.mn <- mean(doletree.slope.rs, na.rm=T)
doletree.slope.rs.lo <- quantile(doletree.slope.rs, probs=0.025, na.rm=T)
doletree.slope.rs.up <- quantile(doletree.slope.rs, probs=0.975, na.rm=T)
doletree.p.rs.mn <- mean(doletree.p.rs, na.rm=T)
doletree.p.rs.lo <- quantile(doletree.p.rs, probs=0.025, na.rm=T)
doletree.p.rs.up <- quantile(doletree.p.rs, probs=0.975, na.rm=T)
doletree.R2.rs.mn <- mean(doletree.R2.rs, na.rm=T)
doletree.R2.rs.lo <- quantile(doletree.R2.rs, probs=0.025, na.rm=T)
doletree.R2.rs.up <- quantile(doletree.R2.rs, probs=0.975, na.rm=T)
doletree.DWp.rs.mn <- mean(doletree.DWp.rs, na.rm=T)
doletree.DWp.rs.lo <- quantile(doletree.DWp.rs, probs=0.025, na.rm=T)
doletree.DWp.rs.up <- quantile(doletree.DWp.rs, probs=0.975, na.rm=T)

Hisoinsol.slope.rs.mn <- mean(Hisoinsol.slope.rs, na.rm=T)
Hisoinsol.slope.rs.lo <- quantile(Hisoinsol.slope.rs, probs=0.025, na.rm=T)
Hisoinsol.slope.rs.up <- quantile(Hisoinsol.slope.rs, probs=0.975, na.rm=T)
Hisoinsol.p.rs.mn <- mean(Hisoinsol.p.rs, na.rm=T)
Hisoinsol.p.rs.lo <- quantile(Hisoinsol.p.rs, probs=0.025, na.rm=T)
Hisoinsol.p.rs.up <- quantile(Hisoinsol.p.rs, probs=0.975, na.rm=T)
Hisoinsol.R2.rs.mn <- mean(Hisoinsol.R2.rs, na.rm=T)
Hisoinsol.R2.rs.lo <- quantile(Hisoinsol.R2.rs, probs=0.025, na.rm=T)
Hisoinsol.R2.rs.up <- quantile(Hisoinsol.R2.rs, probs=0.975, na.rm=T)
Hisoinsol.DWp.rs.mn <- mean(Hisoinsol.DWp.rs, na.rm=T)
Hisoinsol.DWp.rs.lo <- quantile(Hisoinsol.DWp.rs, probs=0.025, na.rm=T)
Hisoinsol.DWp.rs.up <- quantile(Hisoinsol.DWp.rs, probs=0.975, na.rm=T)

Hisodole.slope.rs.mn <- mean(Hisodole.slope.rs, na.rm=T)
Hisodole.slope.rs.lo <- quantile(Hisodole.slope.rs, probs=0.025, na.rm=T)
Hisodole.slope.rs.up <- quantile(Hisodole.slope.rs, probs=0.975, na.rm=T)
Hisodole.p.rs.mn <- mean(Hisodole.p.rs, na.rm=T)
Hisodole.p.rs.lo <- quantile(Hisodole.p.rs, probs=0.025, na.rm=T)
Hisodole.p.rs.up <- quantile(Hisodole.p.rs, probs=0.975, na.rm=T)
Hisodole.R2.rs.mn <- mean(Hisodole.R2.rs, na.rm=T)
Hisodole.R2.rs.lo <- quantile(Hisodole.R2.rs, probs=0.025, na.rm=T)
Hisodole.R2.rs.up <- quantile(Hisodole.R2.rs, probs=0.975, na.rm=T)
Hisodole.DWp.rs.mn <- mean(Hisodole.DWp.rs, na.rm=T)
Hisodole.DWp.rs.lo <- quantile(Hisodole.DWp.rs, probs=0.025, na.rm=T)
Hisodole.DWp.rs.up <- quantile(Hisodole.DWp.rs, probs=0.975, na.rm=T)

HisoprcpGcorL.slope.rs.mn <- mean(HisoprcpGcorL.slope.rs, na.rm=T)
HisoprcpGcorL.slope.rs.lo <- quantile(HisoprcpGcorL.slope.rs, probs=0.025, na.rm=T)
HisoprcpGcorL.slope.rs.up <- quantile(HisoprcpGcorL.slope.rs, probs=0.975, na.rm=T)
HisoprcpGcorL.p.rs.mn <- mean(HisoprcpGcorL.p.rs, na.rm=T)
HisoprcpGcorL.p.rs.lo <- quantile(HisoprcpGcorL.p.rs, probs=0.025, na.rm=T)
HisoprcpGcorL.p.rs.up <- quantile(HisoprcpGcorL.p.rs, probs=0.975, na.rm=T)
HisoprcpGcorL.R2.rs.mn <- mean(HisoprcpGcorL.R2.rs, na.rm=T)
HisoprcpGcorL.R2.rs.lo <- quantile(HisoprcpGcorL.R2.rs, probs=0.025, na.rm=T)
HisoprcpGcorL.R2.rs.up <- quantile(HisoprcpGcorL.R2.rs, probs=0.975, na.rm=T)
HisoprcpGcorL.DWp.rs.mn <- mean(HisoprcpGcorL.DWp.rs, na.rm=T)
HisoprcpGcorL.DWp.rs.lo <- quantile(HisoprcpGcorL.DWp.rs, probs=0.025, na.rm=T)
HisoprcpGcorL.DWp.rs.up <- quantile(HisoprcpGcorL.DWp.rs, probs=0.975, na.rm=T)

HisoprcpGL.slope.rs.mn <- mean(HisoprcpGL.slope.rs, na.rm=T)
HisoprcpGL.slope.rs.lo <- quantile(HisoprcpGL.slope.rs, probs=0.025, na.rm=T)
HisoprcpGL.slope.rs.up <- quantile(HisoprcpGL.slope.rs, probs=0.975, na.rm=T)
HisoprcpGL.p.rs.mn <- mean(HisoprcpGL.p.rs, na.rm=T)
HisoprcpGL.p.rs.lo <- quantile(HisoprcpGL.p.rs, probs=0.025, na.rm=T)
HisoprcpGL.p.rs.up <- quantile(HisoprcpGL.p.rs, probs=0.975, na.rm=T)
HisoprcpGL.R2.rs.mn <- mean(HisoprcpGL.R2.rs, na.rm=T)
HisoprcpGL.R2.rs.lo <- quantile(HisoprcpGL.R2.rs, probs=0.025, na.rm=T)
HisoprcpGL.R2.rs.up <- quantile(HisoprcpGL.R2.rs, probs=0.975, na.rm=T)
HisoprcpGL.DWp.rs.mn <- mean(HisoprcpGL.DWp.rs, na.rm=T)
HisoprcpGL.DWp.rs.lo <- quantile(HisoprcpGL.DWp.rs, probs=0.025, na.rm=T)
HisoprcpGL.DWp.rs.up <- quantile(HisoprcpGL.DWp.rs, probs=0.975, na.rm=T)

HisoprcpL.slope.rs.mn <- mean(HisoprcpL.slope.rs, na.rm=T)
HisoprcpL.slope.rs.lo <- quantile(HisoprcpL.slope.rs, probs=0.025, na.rm=T)
HisoprcpL.slope.rs.up <- quantile(HisoprcpL.slope.rs, probs=0.975, na.rm=T)
HisoprcpL.p.rs.mn <- mean(HisoprcpL.p.rs, na.rm=T)
HisoprcpL.p.rs.lo <- quantile(HisoprcpL.p.rs, probs=0.025, na.rm=T)
HisoprcpL.p.rs.up <- quantile(HisoprcpL.p.rs, probs=0.975, na.rm=T)
HisoprcpL.R2.rs.mn <- mean(HisoprcpL.R2.rs, na.rm=T)
HisoprcpL.R2.rs.lo <- quantile(HisoprcpL.R2.rs, probs=0.025, na.rm=T)
HisoprcpL.R2.rs.up <- quantile(HisoprcpL.R2.rs, probs=0.975, na.rm=T)
HisoprcpL.DWp.rs.mn <- mean(HisoprcpL.DWp.rs, na.rm=T)
HisoprcpL.DWp.rs.lo <- quantile(HisoprcpL.DWp.rs, probs=0.025, na.rm=T)
HisoprcpL.DWp.rs.up <- quantile(HisoprcpL.DWp.rs, probs=0.975, na.rm=T)

HisoprcpGcorH.slope.rs.mn <- mean(HisoprcpGcorH.slope.rs, na.rm=T)
HisoprcpGcorH.slope.rs.lo <- quantile(HisoprcpGcorH.slope.rs, probs=0.025, na.rm=T)
HisoprcpGcorH.slope.rs.up <- quantile(HisoprcpGcorH.slope.rs, probs=0.975, na.rm=T)
HisoprcpGcorH.p.rs.mn <- mean(HisoprcpGcorH.p.rs, na.rm=T)
HisoprcpGcorH.p.rs.lo <- quantile(HisoprcpGcorH.p.rs, probs=0.025, na.rm=T)
HisoprcpGcorH.p.rs.up <- quantile(HisoprcpGcorH.p.rs, probs=0.975, na.rm=T)
HisoprcpGcorH.R2.rs.mn <- mean(HisoprcpGcorH.R2.rs, na.rm=T)
HisoprcpGcorH.R2.rs.lo <- quantile(HisoprcpGcorH.R2.rs, probs=0.025, na.rm=T)
HisoprcpGcorH.R2.rs.up <- quantile(HisoprcpGcorH.R2.rs, probs=0.975, na.rm=T)
HisoprcpGcorH.DWp.rs.mn <- mean(HisoprcpGcorH.DWp.rs, na.rm=T)
HisoprcpGcorH.DWp.rs.lo <- quantile(HisoprcpGcorH.DWp.rs, probs=0.025, na.rm=T)
HisoprcpGcorH.DWp.rs.up <- quantile(HisoprcpGcorH.DWp.rs, probs=0.975, na.rm=T)

HisoprcpGH.slope.rs.mn <- mean(HisoprcpGH.slope.rs, na.rm=T)
HisoprcpGH.slope.rs.lo <- quantile(HisoprcpGH.slope.rs, probs=0.025, na.rm=T)
HisoprcpGH.slope.rs.up <- quantile(HisoprcpGH.slope.rs, probs=0.975, na.rm=T)
HisoprcpGH.p.rs.mn <- mean(HisoprcpGH.p.rs, na.rm=T)
HisoprcpGH.p.rs.lo <- quantile(HisoprcpGH.p.rs, probs=0.025, na.rm=T)
HisoprcpGH.p.rs.up <- quantile(HisoprcpGH.p.rs, probs=0.975, na.rm=T)
HisoprcpGH.R2.rs.mn <- mean(HisoprcpGH.R2.rs, na.rm=T)
HisoprcpGH.R2.rs.lo <- quantile(HisoprcpGH.R2.rs, probs=0.025, na.rm=T)
HisoprcpGH.R2.rs.up <- quantile(HisoprcpGH.R2.rs, probs=0.975, na.rm=T)
HisoprcpGH.DWp.rs.mn <- mean(HisoprcpGH.DWp.rs, na.rm=T)
HisoprcpGH.DWp.rs.lo <- quantile(HisoprcpGH.DWp.rs, probs=0.025, na.rm=T)
HisoprcpGH.DWp.rs.up <- quantile(HisoprcpGH.DWp.rs, probs=0.975, na.rm=T)

HisoprcpH.slope.rs.mn <- mean(HisoprcpH.slope.rs, na.rm=T)
HisoprcpH.slope.rs.lo <- quantile(HisoprcpH.slope.rs, probs=0.025, na.rm=T)
HisoprcpH.slope.rs.up <- quantile(HisoprcpH.slope.rs, probs=0.975, na.rm=T)
HisoprcpH.p.rs.mn <- mean(HisoprcpH.p.rs, na.rm=T)
HisoprcpH.p.rs.lo <- quantile(HisoprcpH.p.rs, probs=0.025, na.rm=T)
HisoprcpH.p.rs.up <- quantile(HisoprcpH.p.rs, probs=0.975, na.rm=T)
HisoprcpH.R2.rs.mn <- mean(HisoprcpH.R2.rs, na.rm=T)
HisoprcpH.R2.rs.lo <- quantile(HisoprcpH.R2.rs, probs=0.025, na.rm=T)
HisoprcpH.R2.rs.up <- quantile(HisoprcpH.R2.rs, probs=0.975, na.rm=T)
HisoprcpH.DWp.rs.mn <- mean(HisoprcpH.DWp.rs, na.rm=T)
HisoprcpH.DWp.rs.lo <- quantile(HisoprcpH.DWp.rs, probs=0.025, na.rm=T)
HisoprcpH.DWp.rs.up <- quantile(HisoprcpH.DWp.rs, probs=0.975, na.rm=T)

Hisotree.slope.rs.mn <- mean(Hisotree.slope.rs, na.rm=T)
Hisotree.slope.rs.lo <- quantile(Hisotree.slope.rs, probs=0.025, na.rm=T)
Hisotree.slope.rs.up <- quantile(Hisotree.slope.rs, probs=0.975, na.rm=T)
Hisotree.p.rs.mn <- mean(Hisotree.p.rs, na.rm=T)
Hisotree.p.rs.lo <- quantile(Hisotree.p.rs, probs=0.025, na.rm=T)
Hisotree.p.rs.up <- quantile(Hisotree.p.rs, probs=0.975, na.rm=T)
Hisotree.R2.rs.mn <- mean(Hisotree.R2.rs, na.rm=T)
Hisotree.R2.rs.lo <- quantile(Hisotree.R2.rs, probs=0.025, na.rm=T)
Hisotree.R2.rs.up <- quantile(Hisotree.R2.rs, probs=0.975, na.rm=T)
Hisotree.DWp.rs.mn <- mean(Hisotree.DWp.rs, na.rm=T)
Hisotree.DWp.rs.lo <- quantile(Hisotree.DWp.rs, probs=0.025, na.rm=T)
Hisotree.DWp.rs.up <- quantile(Hisotree.DWp.rs, probs=0.975, na.rm=T)

insolprcpGcorL.slope.rs.mn <- mean(insolprcpGcorL.slope.rs, na.rm=T)
insolprcpGcorL.slope.rs.lo <- quantile(insolprcpGcorL.slope.rs, probs=0.025, na.rm=T)
insolprcpGcorL.slope.rs.up <- quantile(insolprcpGcorL.slope.rs, probs=0.975, na.rm=T)
insolprcpGcorL.p.rs.mn <- mean(insolprcpGcorL.p.rs, na.rm=T)
insolprcpGcorL.p.rs.lo <- quantile(insolprcpGcorL.p.rs, probs=0.025, na.rm=T)
insolprcpGcorL.p.rs.up <- quantile(insolprcpGcorL.p.rs, probs=0.975, na.rm=T)
insolprcpGcorL.R2.rs.mn <- mean(insolprcpGcorL.R2.rs, na.rm=T)
insolprcpGcorL.R2.rs.lo <- quantile(insolprcpGcorL.R2.rs, probs=0.025, na.rm=T)
insolprcpGcorL.R2.rs.up <- quantile(insolprcpGcorL.R2.rs, probs=0.975, na.rm=T)
insolprcpGcorL.DWp.rs.mn <- mean(insolprcpGcorL.DWp.rs, na.rm=T)
insolprcpGcorL.DWp.rs.lo <- quantile(insolprcpGcorL.DWp.rs, probs=0.025, na.rm=T)
insolprcpGcorL.DWp.rs.up <- quantile(insolprcpGcorL.DWp.rs, probs=0.975, na.rm=T)

insolprcpGL.slope.rs.mn <- mean(insolprcpGL.slope.rs, na.rm=T)
insolprcpGL.slope.rs.lo <- quantile(insolprcpGL.slope.rs, probs=0.025, na.rm=T)
insolprcpGL.slope.rs.up <- quantile(insolprcpGL.slope.rs, probs=0.975, na.rm=T)
insolprcpGL.p.rs.mn <- mean(insolprcpGL.p.rs, na.rm=T)
insolprcpGL.p.rs.lo <- quantile(insolprcpGL.p.rs, probs=0.025, na.rm=T)
insolprcpGL.p.rs.up <- quantile(insolprcpGL.p.rs, probs=0.975, na.rm=T)
insolprcpGL.R2.rs.mn <- mean(insolprcpGL.R2.rs, na.rm=T)
insolprcpGL.R2.rs.lo <- quantile(insolprcpGL.R2.rs, probs=0.025, na.rm=T)
insolprcpGL.R2.rs.up <- quantile(insolprcpGL.R2.rs, probs=0.975, na.rm=T)
insolprcpGL.DWp.rs.mn <- mean(insolprcpGL.DWp.rs, na.rm=T)
insolprcpGL.DWp.rs.lo <- quantile(insolprcpGL.DWp.rs, probs=0.025, na.rm=T)
insolprcpGL.DWp.rs.up <- quantile(insolprcpGL.DWp.rs, probs=0.975, na.rm=T)

insolprcpL.slope.rs.mn <- mean(insolprcpL.slope.rs, na.rm=T)
insolprcpL.slope.rs.lo <- quantile(insolprcpL.slope.rs, probs=0.025, na.rm=T)
insolprcpL.slope.rs.up <- quantile(insolprcpL.slope.rs, probs=0.975, na.rm=T)
insolprcpL.p.rs.mn <- mean(insolprcpL.p.rs, na.rm=T)
insolprcpL.p.rs.lo <- quantile(insolprcpL.p.rs, probs=0.025, na.rm=T)
insolprcpL.p.rs.up <- quantile(insolprcpL.p.rs, probs=0.975, na.rm=T)
insolprcpL.R2.rs.mn <- mean(insolprcpL.R2.rs, na.rm=T)
insolprcpL.R2.rs.lo <- quantile(insolprcpL.R2.rs, probs=0.025, na.rm=T)
insolprcpL.R2.rs.up <- quantile(insolprcpL.R2.rs, probs=0.975, na.rm=T)
insolprcpL.DWp.rs.mn <- mean(insolprcpL.DWp.rs, na.rm=T)
insolprcpL.DWp.rs.lo <- quantile(insolprcpL.DWp.rs, probs=0.025, na.rm=T)
insolprcpL.DWp.rs.up <- quantile(insolprcpL.DWp.rs, probs=0.975, na.rm=T)

insolprcpGcorH.slope.rs.mn <- mean(insolprcpGcorH.slope.rs, na.rm=T)
insolprcpGcorH.slope.rs.lo <- quantile(insolprcpGcorH.slope.rs, probs=0.025, na.rm=T)
insolprcpGcorH.slope.rs.up <- quantile(insolprcpGcorH.slope.rs, probs=0.975, na.rm=T)
insolprcpGcorH.p.rs.mn <- mean(insolprcpGcorH.p.rs, na.rm=T)
insolprcpGcorH.p.rs.lo <- quantile(insolprcpGcorH.p.rs, probs=0.025, na.rm=T)
insolprcpGcorH.p.rs.up <- quantile(insolprcpGcorH.p.rs, probs=0.975, na.rm=T)
insolprcpGcorH.R2.rs.mn <- mean(insolprcpGcorH.R2.rs, na.rm=T)
insolprcpGcorH.R2.rs.lo <- quantile(insolprcpGcorH.R2.rs, probs=0.025, na.rm=T)
insolprcpGcorH.R2.rs.up <- quantile(insolprcpGcorH.R2.rs, probs=0.975, na.rm=T)
insolprcpGcorH.DWp.rs.mn <- mean(insolprcpGcorH.DWp.rs, na.rm=T)
insolprcpGcorH.DWp.rs.lo <- quantile(insolprcpGcorH.DWp.rs, probs=0.025, na.rm=T)
insolprcpGcorH.DWp.rs.up <- quantile(insolprcpGcorH.DWp.rs, probs=0.975, na.rm=T)

insolprcpGH.slope.rs.mn <- mean(insolprcpGH.slope.rs, na.rm=T)
insolprcpGH.slope.rs.lo <- quantile(insolprcpGH.slope.rs, probs=0.025, na.rm=T)
insolprcpGH.slope.rs.up <- quantile(insolprcpGH.slope.rs, probs=0.975, na.rm=T)
insolprcpGH.p.rs.mn <- mean(insolprcpGH.p.rs, na.rm=T)
insolprcpGH.p.rs.lo <- quantile(insolprcpGH.p.rs, probs=0.025, na.rm=T)
insolprcpGH.p.rs.up <- quantile(insolprcpGH.p.rs, probs=0.975, na.rm=T)
insolprcpGH.R2.rs.mn <- mean(insolprcpGH.R2.rs, na.rm=T)
insolprcpGH.R2.rs.lo <- quantile(insolprcpGH.R2.rs, probs=0.025, na.rm=T)
insolprcpGH.R2.rs.up <- quantile(insolprcpGH.R2.rs, probs=0.975, na.rm=T)
insolprcpGH.DWp.rs.mn <- mean(insolprcpGH.DWp.rs, na.rm=T)
insolprcpGH.DWp.rs.lo <- quantile(insolprcpGH.DWp.rs, probs=0.025, na.rm=T)
insolprcpGH.DWp.rs.up <- quantile(insolprcpGH.DWp.rs, probs=0.975, na.rm=T)

insolprcpH.slope.rs.mn <- mean(insolprcpH.slope.rs, na.rm=T)
insolprcpH.slope.rs.lo <- quantile(insolprcpH.slope.rs, probs=0.025, na.rm=T)
insolprcpH.slope.rs.up <- quantile(insolprcpH.slope.rs, probs=0.975, na.rm=T)
insolprcpH.p.rs.mn <- mean(insolprcpH.p.rs, na.rm=T)
insolprcpH.p.rs.lo <- quantile(insolprcpH.p.rs, probs=0.025, na.rm=T)
insolprcpH.p.rs.up <- quantile(insolprcpH.p.rs, probs=0.975, na.rm=T)
insolprcpH.R2.rs.mn <- mean(insolprcpH.R2.rs, na.rm=T)
insolprcpH.R2.rs.lo <- quantile(insolprcpH.R2.rs, probs=0.025, na.rm=T)
insolprcpH.R2.rs.up <- quantile(insolprcpH.R2.rs, probs=0.975, na.rm=T)
insolprcpH.DWp.rs.mn <- mean(insolprcpH.DWp.rs, na.rm=T)
insolprcpH.DWp.rs.lo <- quantile(insolprcpH.DWp.rs, probs=0.025, na.rm=T)
insolprcpH.DWp.rs.up <- quantile(insolprcpH.DWp.rs, probs=0.975, na.rm=T)

insoltree.slope.rs.mn <- mean(insoltree.slope.rs, na.rm=T)
insoltree.slope.rs.lo <- quantile(insoltree.slope.rs, probs=0.025, na.rm=T)
insoltree.slope.rs.up <- quantile(insoltree.slope.rs, probs=0.975, na.rm=T)
insoltree.p.rs.mn <- mean(insoltree.p.rs, na.rm=T)
insoltree.p.rs.lo <- quantile(insoltree.p.rs, probs=0.025, na.rm=T)
insoltree.p.rs.up <- quantile(insoltree.p.rs, probs=0.975, na.rm=T)
insoltree.R2.rs.mn <- mean(insoltree.R2.rs, na.rm=T)
insoltree.R2.rs.lo <- quantile(insoltree.R2.rs, probs=0.025, na.rm=T)
insoltree.R2.rs.up <- quantile(insoltree.R2.rs, probs=0.975, na.rm=T)
insoltree.DWp.rs.mn <- mean(insoltree.DWp.rs, na.rm=T)
insoltree.DWp.rs.lo <- quantile(insoltree.DWp.rs, probs=0.025, na.rm=T)
insoltree.DWp.rs.up <- quantile(insoltree.DWp.rs, probs=0.975, na.rm=T)

treeprcpGL.slope.rs.mn <- mean(treeprcpGL.slope.rs, na.rm=T)
treeprcpGL.slope.rs.lo <- quantile(treeprcpGL.slope.rs, probs=0.025, na.rm=T)
treeprcpGL.slope.rs.up <- quantile(treeprcpGL.slope.rs, probs=0.975, na.rm=T)
treeprcpGL.p.rs.mn <- mean(treeprcpGL.p.rs, na.rm=T)
treeprcpGL.p.rs.lo <- quantile(treeprcpGL.p.rs, probs=0.025, na.rm=T)
treeprcpGL.p.rs.up <- quantile(treeprcpGL.p.rs, probs=0.975, na.rm=T)
treeprcpGL.R2.rs.mn <- mean(treeprcpGL.R2.rs, na.rm=T)
treeprcpGL.R2.rs.lo <- quantile(treeprcpGL.R2.rs, probs=0.025, na.rm=T)
treeprcpGL.R2.rs.up <- quantile(treeprcpGL.R2.rs, probs=0.975, na.rm=T)
treeprcpGL.DWp.rs.mn <- mean(treeprcpGL.DWp.rs, na.rm=T)
treeprcpGL.DWp.rs.lo <- quantile(treeprcpGL.DWp.rs, probs=0.025, na.rm=T)
treeprcpGL.DWp.rs.up <- quantile(treeprcpGL.DWp.rs, probs=0.975, na.rm=T)

treeprcpL.slope.rs.mn <- mean(treeprcpL.slope.rs, na.rm=T)
treeprcpL.slope.rs.lo <- quantile(treeprcpL.slope.rs, probs=0.025, na.rm=T)
treeprcpL.slope.rs.up <- quantile(treeprcpL.slope.rs, probs=0.975, na.rm=T)
treeprcpL.p.rs.mn <- mean(treeprcpL.p.rs, na.rm=T)
treeprcpL.p.rs.lo <- quantile(treeprcpL.p.rs, probs=0.025, na.rm=T)
treeprcpL.p.rs.up <- quantile(treeprcpL.p.rs, probs=0.975, na.rm=T)
treeprcpL.R2.rs.mn <- mean(treeprcpL.R2.rs, na.rm=T)
treeprcpL.R2.rs.lo <- quantile(treeprcpL.R2.rs, probs=0.025, na.rm=T)
treeprcpL.R2.rs.up <- quantile(treeprcpL.R2.rs, probs=0.975, na.rm=T)
treeprcpL.DWp.rs.mn <- mean(treeprcpL.DWp.rs, na.rm=T)
treeprcpL.DWp.rs.lo <- quantile(treeprcpL.DWp.rs, probs=0.025, na.rm=T)
treeprcpL.DWp.rs.up <- quantile(treeprcpL.DWp.rs, probs=0.975, na.rm=T)

treeprcpGcorH.slope.rs.mn <- mean(treeprcpGcorH.slope.rs, na.rm=T)
treeprcpGcorH.slope.rs.lo <- quantile(treeprcpGcorH.slope.rs, probs=0.025, na.rm=T)
treeprcpGcorH.slope.rs.up <- quantile(treeprcpGcorH.slope.rs, probs=0.975, na.rm=T)
treeprcpGcorH.p.rs.mn <- mean(treeprcpGcorH.p.rs, na.rm=T)
treeprcpGcorH.p.rs.lo <- quantile(treeprcpGcorH.p.rs, probs=0.025, na.rm=T)
treeprcpGcorH.p.rs.up <- quantile(treeprcpGcorH.p.rs, probs=0.975, na.rm=T)
treeprcpGcorH.R2.rs.mn <- mean(treeprcpGcorH.R2.rs, na.rm=T)
treeprcpGcorH.R2.rs.lo <- quantile(treeprcpGcorH.R2.rs, probs=0.025, na.rm=T)
treeprcpGcorH.R2.rs.up <- quantile(treeprcpGcorH.R2.rs, probs=0.975, na.rm=T)
treeprcpGcorH.DWp.rs.mn <- mean(treeprcpGcorH.DWp.rs, na.rm=T)
treeprcpGcorH.DWp.rs.lo <- quantile(treeprcpGcorH.DWp.rs, probs=0.025, na.rm=T)
treeprcpGcorH.DWp.rs.up <- quantile(treeprcpGcorH.DWp.rs, probs=0.975, na.rm=T)

treeprcpGcorL.slope.rs.mn <- mean(treeprcpGcorL.slope.rs, na.rm=T)
treeprcpGcorL.slope.rs.lo <- quantile(treeprcpGcorL.slope.rs, probs=0.025, na.rm=T)
treeprcpGcorL.slope.rs.up <- quantile(treeprcpGcorL.slope.rs, probs=0.975, na.rm=T)
treeprcpGcorL.p.rs.mn <- mean(treeprcpGcorL.p.rs, na.rm=T)
treeprcpGcorL.p.rs.lo <- quantile(treeprcpGcorL.p.rs, probs=0.025, na.rm=T)
treeprcpGcorL.p.rs.up <- quantile(treeprcpGcorL.p.rs, probs=0.975, na.rm=T)
treeprcpGcorL.R2.rs.mn <- mean(treeprcpGcorL.R2.rs, na.rm=T)
treeprcpGcorL.R2.rs.lo <- quantile(treeprcpGcorL.R2.rs, probs=0.025, na.rm=T)
treeprcpGcorL.R2.rs.up <- quantile(treeprcpGcorL.R2.rs, probs=0.975, na.rm=T)
treeprcpGcorL.DWp.rs.mn <- mean(treeprcpGcorL.DWp.rs, na.rm=T)
treeprcpGcorL.DWp.rs.lo <- quantile(treeprcpGcorL.DWp.rs, probs=0.025, na.rm=T)
treeprcpGcorL.DWp.rs.up <- quantile(treeprcpGcorL.DWp.rs, probs=0.975, na.rm=T)

treeprcpGH.slope.rs.mn <- mean(treeprcpGH.slope.rs, na.rm=T)
treeprcpGH.slope.rs.lo <- quantile(treeprcpGH.slope.rs, probs=0.025, na.rm=T)
treeprcpGH.slope.rs.up <- quantile(treeprcpGH.slope.rs, probs=0.975, na.rm=T)
treeprcpGH.p.rs.mn <- mean(treeprcpGH.p.rs, na.rm=T)
treeprcpGH.p.rs.lo <- quantile(treeprcpGH.p.rs, probs=0.025, na.rm=T)
treeprcpGH.p.rs.up <- quantile(treeprcpGH.p.rs, probs=0.975, na.rm=T)
treeprcpGH.R2.rs.mn <- mean(treeprcpGH.R2.rs, na.rm=T)
treeprcpGH.R2.rs.lo <- quantile(treeprcpGH.R2.rs, probs=0.025, na.rm=T)
treeprcpGH.R2.rs.up <- quantile(treeprcpGH.R2.rs, probs=0.975, na.rm=T)
treeprcpGH.DWp.rs.mn <- mean(treeprcpGH.DWp.rs, na.rm=T)
treeprcpGH.DWp.rs.lo <- quantile(treeprcpGH.DWp.rs, probs=0.025, na.rm=T)
treeprcpGH.DWp.rs.up <- quantile(treeprcpGH.DWp.rs, probs=0.975, na.rm=T)

treeprcpH.slope.rs.mn <- mean(treeprcpH.slope.rs, na.rm=T)
treeprcpH.slope.rs.lo <- quantile(treeprcpH.slope.rs, probs=0.025, na.rm=T)
treeprcpH.slope.rs.up <- quantile(treeprcpH.slope.rs, probs=0.975, na.rm=T)
treeprcpH.p.rs.mn <- mean(treeprcpH.p.rs, na.rm=T)
treeprcpH.p.rs.lo <- quantile(treeprcpH.p.rs, probs=0.025, na.rm=T)
treeprcpH.p.rs.up <- quantile(treeprcpH.p.rs, probs=0.975, na.rm=T)
treeprcpH.R2.rs.mn <- mean(treeprcpH.R2.rs, na.rm=T)
treeprcpH.R2.rs.lo <- quantile(treeprcpH.R2.rs, probs=0.025, na.rm=T)
treeprcpH.R2.rs.up <- quantile(treeprcpH.R2.rs, probs=0.975, na.rm=T)
treeprcpH.DWp.rs.mn <- mean(treeprcpH.DWp.rs, na.rm=T)
treeprcpH.DWp.rs.lo <- quantile(treeprcpH.DWp.rs, probs=0.025, na.rm=T)
treeprcpH.DWp.rs.up <- quantile(treeprcpH.DWp.rs, probs=0.975, na.rm=T)

# place lo/up into matrices following cor.mean format
p.rs.lo <- cor.mean
p.rs.lo[,] <- NA
p.rs.lo[upper.tri(p.rs.lo)] <- 0
p.rs.lo[2,1] <- insolprcpGcorL.p.rs.lo
p.rs.lo[3,1] <- insolprcpGL.p.rs.lo
p.rs.lo[4,1] <- insolprcpGcorH.p.rs.lo
p.rs.lo[5,1] <- insolprcpGH.p.rs.lo
p.rs.lo[6,1] <- MaO18insol.p.rs.lo
p.rs.lo[7,1] <- doleinsol.p.rs.lo
p.rs.lo[8,1] <- insolprcpL.p.rs.lo
p.rs.lo[9,1] <- insolprcpH.p.rs.lo
p.rs.lo[11,1] <- Hisoinsol.p.rs.lo
p.rs.lo[12,1] <- insoltree.p.rs.lo
p.rs.lo[6,3] <- MaO18prcpGcorL.p.rs.lo
p.rs.lo[7,3] <- doleprcpGcorL.p.rs.lo
p.rs.lo[11,3] <- HisoprcpGcorL.p.rs.lo
p.rs.lo[12,3] <- treeprcpGcorL.p.rs.lo;
p.rs.lo[6,4] <- MaO18prcpGL.p.rs.lo
p.rs.lo[7,4] <- doleprcpGL.p.rs.lo
p.rs.lo[11,4] <- HisoprcpGL.p.rs.lo
p.rs.lo[12,4] <- treeprcpGL.p.rs.lo
p.rs.lo[6,5] <- MaO18prcpGcorH.p.rs.lo
p.rs.lo[7,5] <- doleprcpGcorH.p.rs.lo
p.rs.lo[11,5] <- HisoprcpGcorH.p.rs.lo
p.rs.lo[12,5] <- treeprcpGcorH.p.rs.lo 
p.rs.lo[6,6] <- MaO18prcpGH.p.rs.lo
p.rs.lo[7,6] <- doleprcpGH.p.rs.lo
p.rs.lo[11,6] <- HisoprcpGH.p.rs.lo
p.rs.lo[12,6] <- treeprcpGH.p.rs.lo;
p.rs.lo[7,7] <- MaO18dole.p.rs.lo
p.rs.lo[8,7] <- MaO18prcpL.p.rs.lo
p.rs.lo[9,7] <- MaO18prcpH.p.rs.lo
p.rs.lo[11,7] <- MaO18Hiso.p.rs.lo
p.rs.lo[12,7] <- MaO18tree.p.rs.lo
p.rs.lo[8,8] <- doleprcpL.p.rs.lo
p.rs.lo[9,8] <- doleprcpH.p.rs.lo
p.rs.lo[11,8] <- Hisodole.p.rs.lo
p.rs.lo[12,8] <- doletree.p.rs.lo
p.rs.lo[11,9] <- HisoprcpL.p.rs.lo
p.rs.lo[12,9] <- treeprcpL.p.rs.lo
p.rs.lo[11,10] <- HisoprcpH.p.rs.lo
p.rs.lo[12,10] <- treeprcpH.p.rs.lo
p.rs.lo[12,12] <- Hisotree.p.rs.lo
round(p.rs.lo, 3)

p.rs.up <- cor.mean
p.rs.up[,] <- NA
p.rs.up[upper.tri(p.rs.up)] <- 0
p.rs.up[2,1] <- insolprcpGcorL.p.rs.up
p.rs.up[3,1] <- insolprcpGL.p.rs.up
p.rs.up[4,1] <- insolprcpGcorH.p.rs.up
p.rs.up[5,1] <- insolprcpGH.p.rs.up
p.rs.up[6,1] <- MaO18insol.p.rs.up
p.rs.up[7,1] <- doleinsol.p.rs.up
p.rs.up[8,1] <- insolprcpL.p.rs.up
p.rs.up[9,1] <- insolprcpH.p.rs.up
p.rs.up[11,1] <- Hisoinsol.p.rs.up
p.rs.up[12,1] <- insoltree.p.rs.up
p.rs.up[6,3] <- MaO18prcpGcorL.p.rs.up
p.rs.up[7,3] <- doleprcpGcorL.p.rs.up
p.rs.up[11,3] <- HisoprcpGcorL.p.rs.up
p.rs.up[12,3] <- treeprcpGcorL.p.rs.up;
p.rs.up[6,4] <- MaO18prcpGL.p.rs.up
p.rs.up[7,4] <- doleprcpGL.p.rs.up
p.rs.up[11,4] <- HisoprcpGL.p.rs.up
p.rs.up[12,4] <- treeprcpGL.p.rs.up
p.rs.up[6,5] <- MaO18prcpGcorH.p.rs.up
p.rs.up[7,5] <- doleprcpGcorH.p.rs.up
p.rs.up[11,5] <- HisoprcpGcorH.p.rs.up
p.rs.up[12,5] <- treeprcpGcorH.p.rs.up 
p.rs.up[6,6] <- MaO18prcpGH.p.rs.up
p.rs.up[7,6] <- doleprcpGH.p.rs.up
p.rs.up[11,6] <- HisoprcpGH.p.rs.up
p.rs.up[12,6] <- treeprcpGH.p.rs.up;
p.rs.up[7,7] <- MaO18dole.p.rs.up
p.rs.up[8,7] <- MaO18prcpL.p.rs.up
p.rs.up[9,7] <- MaO18prcpH.p.rs.up
p.rs.up[11,7] <- MaO18Hiso.p.rs.up
p.rs.up[12,7] <- MaO18tree.p.rs.up
p.rs.up[8,8] <- doleprcpL.p.rs.up
p.rs.up[9,8] <- doleprcpH.p.rs.up
p.rs.up[11,8] <- Hisodole.p.rs.up
p.rs.up[12,8] <- doletree.p.rs.up
p.rs.up[11,9] <- HisoprcpL.p.rs.up
p.rs.up[12,9] <- treeprcpL.p.rs.up
p.rs.up[11,10] <- HisoprcpH.p.rs.up
p.rs.up[12,10] <- treeprcpH.p.rs.up
p.rs.up[12,12] <- Hisotree.p.rs.up
round(p.rs.up, 3)

p.rs.mn <- cor.mean
p.rs.mn[,] <- NA
p.rs.mn[upper.tri(p.rs.mn)] <- 0
p.rs.mn[2,1] <- insolprcpGcorL.p.rs.mn
p.rs.mn[3,1] <- insolprcpGL.p.rs.mn
p.rs.mn[4,1] <- insolprcpGcorH.p.rs.mn
p.rs.mn[5,1] <- insolprcpGH.p.rs.mn
p.rs.mn[6,1] <- MaO18insol.p.rs.mn
p.rs.mn[7,1] <- doleinsol.p.rs.mn
p.rs.mn[8,1] <- insolprcpL.p.rs.mn
p.rs.mn[9,1] <- insolprcpH.p.rs.mn
p.rs.mn[11,1] <- Hisoinsol.p.rs.mn
p.rs.mn[12,1] <- insoltree.p.rs.mn
p.rs.mn[6,3] <- MaO18prcpGcorL.p.rs.mn
p.rs.mn[7,3] <- doleprcpGcorL.p.rs.mn
p.rs.mn[11,3] <- HisoprcpGcorL.p.rs.mn
p.rs.mn[12,3] <- treeprcpGcorL.p.rs.mn;
p.rs.mn[6,4] <- MaO18prcpGL.p.rs.mn
p.rs.mn[7,4] <- doleprcpGL.p.rs.mn
p.rs.mn[11,4] <- HisoprcpGL.p.rs.mn
p.rs.mn[12,4] <- treeprcpGL.p.rs.mn
p.rs.mn[6,5] <- MaO18prcpGcorH.p.rs.mn
p.rs.mn[7,5] <- doleprcpGcorH.p.rs.mn
p.rs.mn[11,5] <- HisoprcpGcorH.p.rs.mn
p.rs.mn[12,5] <- treeprcpGcorH.p.rs.mn 
p.rs.mn[6,6] <- MaO18prcpGH.p.rs.mn
p.rs.mn[7,6] <- doleprcpGH.p.rs.mn
p.rs.mn[11,6] <- HisoprcpGH.p.rs.mn
p.rs.mn[12,6] <- treeprcpGH.p.rs.mn;
p.rs.mn[7,7] <- MaO18dole.p.rs.mn
p.rs.mn[8,7] <- MaO18prcpL.p.rs.mn
p.rs.mn[9,7] <- MaO18prcpH.p.rs.mn
p.rs.mn[11,7] <- MaO18Hiso.p.rs.mn
p.rs.mn[12,7] <- MaO18tree.p.rs.mn
p.rs.mn[8,8] <- doleprcpL.p.rs.mn
p.rs.mn[9,8] <- doleprcpH.p.rs.mn
p.rs.mn[11,8] <- Hisodole.p.rs.mn
p.rs.mn[12,8] <- doletree.p.rs.mn
p.rs.mn[11,9] <- HisoprcpL.p.rs.mn
p.rs.mn[12,9] <- treeprcpL.p.rs.mn
p.rs.mn[11,10] <- HisoprcpH.p.rs.mn
p.rs.mn[12,10] <- treeprcpH.p.rs.mn
p.rs.mn[12,12] <- Hisotree.p.rs.mn
round(p.rs.mn, 3)

R2.rs.lo <- cor.mean
R2.rs.lo[,] <- NA
R2.rs.lo[upper.tri(R2.rs.lo)] <- 0
R2.rs.lo[2,1] <- insolprcpGcorL.R2.rs.lo
R2.rs.lo[3,1] <- insolprcpGL.R2.rs.lo
R2.rs.lo[4,1] <- insolprcpGcorH.R2.rs.lo
R2.rs.lo[5,1] <- insolprcpGH.R2.rs.lo
R2.rs.lo[6,1] <- MaO18insol.R2.rs.lo
R2.rs.lo[7,1] <- doleinsol.R2.rs.lo
R2.rs.lo[8,1] <- insolprcpL.R2.rs.lo
R2.rs.lo[9,1] <- insolprcpH.R2.rs.lo
R2.rs.lo[11,1] <- Hisoinsol.R2.rs.lo
R2.rs.lo[12,1] <- insoltree.R2.rs.lo
R2.rs.lo[6,3] <- MaO18prcpGcorL.R2.rs.lo
R2.rs.lo[7,3] <- doleprcpGcorL.R2.rs.lo
R2.rs.lo[11,3] <- HisoprcpGcorL.R2.rs.lo
R2.rs.lo[12,3] <- treeprcpGcorL.R2.rs.lo;
R2.rs.lo[6,4] <- MaO18prcpGL.R2.rs.lo
R2.rs.lo[7,4] <- doleprcpGL.R2.rs.lo
R2.rs.lo[11,4] <- HisoprcpGL.R2.rs.lo
R2.rs.lo[12,4] <- treeprcpGL.R2.rs.lo
R2.rs.lo[6,5] <- MaO18prcpGcorH.R2.rs.lo
R2.rs.lo[7,5] <- doleprcpGcorH.R2.rs.lo
R2.rs.lo[11,5] <- HisoprcpGcorH.R2.rs.lo
R2.rs.lo[12,5] <- treeprcpGcorH.R2.rs.lo 
R2.rs.lo[6,6] <- MaO18prcpGH.R2.rs.lo
R2.rs.lo[7,6] <- doleprcpGH.R2.rs.lo
R2.rs.lo[11,6] <- HisoprcpGH.R2.rs.lo
R2.rs.lo[12,6] <- treeprcpGH.R2.rs.lo;
R2.rs.lo[7,7] <- MaO18dole.R2.rs.lo
R2.rs.lo[8,7] <- MaO18prcpL.R2.rs.lo
R2.rs.lo[9,7] <- MaO18prcpH.R2.rs.lo
R2.rs.lo[11,7] <- MaO18Hiso.R2.rs.lo
R2.rs.lo[12,7] <- MaO18tree.R2.rs.lo
R2.rs.lo[8,8] <- doleprcpL.R2.rs.lo
R2.rs.lo[9,8] <- doleprcpH.R2.rs.lo
R2.rs.lo[11,8] <- Hisodole.R2.rs.lo
R2.rs.lo[12,8] <- doletree.R2.rs.lo
R2.rs.lo[11,9] <- HisoprcpL.R2.rs.lo
R2.rs.lo[12,9] <- treeprcpL.R2.rs.lo
R2.rs.lo[11,10] <- HisoprcpH.R2.rs.lo
R2.rs.lo[12,10] <- treeprcpH.R2.rs.lo
R2.rs.lo[12,12] <- Hisotree.R2.rs.lo
round(R2.rs.lo, 3)

R2.rs.up <- cor.mean
R2.rs.up[,] <- NA
R2.rs.up[upper.tri(R2.rs.up)] <- 0
R2.rs.up[2,1] <- insolprcpGcorL.R2.rs.up
R2.rs.up[3,1] <- insolprcpGL.R2.rs.up
R2.rs.up[4,1] <- insolprcpGcorH.R2.rs.up
R2.rs.up[5,1] <- insolprcpGH.R2.rs.up
R2.rs.up[6,1] <- MaO18insol.R2.rs.up
R2.rs.up[7,1] <- doleinsol.R2.rs.up
R2.rs.up[8,1] <- insolprcpL.R2.rs.up
R2.rs.up[9,1] <- insolprcpH.R2.rs.up
R2.rs.up[11,1] <- Hisoinsol.R2.rs.up
R2.rs.up[12,1] <- insoltree.R2.rs.up
R2.rs.up[6,3] <- MaO18prcpGcorL.R2.rs.up
R2.rs.up[7,3] <- doleprcpGcorL.R2.rs.up
R2.rs.up[11,3] <- HisoprcpGcorL.R2.rs.up
R2.rs.up[12,3] <- treeprcpGcorL.R2.rs.up;
R2.rs.up[6,4] <- MaO18prcpGL.R2.rs.up
R2.rs.up[7,4] <- doleprcpGL.R2.rs.up
R2.rs.up[11,4] <- HisoprcpGL.R2.rs.up
R2.rs.up[12,4] <- treeprcpGL.R2.rs.up
R2.rs.up[6,5] <- MaO18prcpGcorH.R2.rs.up
R2.rs.up[7,5] <- doleprcpGcorH.R2.rs.up
R2.rs.up[11,5] <- HisoprcpGcorH.R2.rs.up
R2.rs.up[12,5] <- treeprcpGcorH.R2.rs.up 
R2.rs.up[6,6] <- MaO18prcpGH.R2.rs.up
R2.rs.up[7,6] <- doleprcpGH.R2.rs.up
R2.rs.up[11,6] <- HisoprcpGH.R2.rs.up
R2.rs.up[12,6] <- treeprcpGH.R2.rs.up;
R2.rs.up[7,7] <- MaO18dole.R2.rs.up
R2.rs.up[8,7] <- MaO18prcpL.R2.rs.up
R2.rs.up[9,7] <- MaO18prcpH.R2.rs.up
R2.rs.up[11,7] <- MaO18Hiso.R2.rs.up
R2.rs.up[12,7] <- MaO18tree.R2.rs.up
R2.rs.up[8,8] <- doleprcpL.R2.rs.up
R2.rs.up[9,8] <- doleprcpH.R2.rs.up
R2.rs.up[11,8] <- Hisodole.R2.rs.up
R2.rs.up[12,8] <- doletree.R2.rs.up
R2.rs.up[11,9] <- HisoprcpL.R2.rs.up
R2.rs.up[12,9] <- treeprcpL.R2.rs.up
R2.rs.up[11,10] <- HisoprcpH.R2.rs.up
R2.rs.up[12,10] <- treeprcpH.R2.rs.up
R2.rs.up[12,12] <- Hisotree.p.lo
round(R2.rs.up, 3)

R2.rs.mn <- cor.mean
R2.rs.mn[,] <- NA
R2.rs.mn[upper.tri(R2.rs.mn)] <- 0
R2.rs.mn[2,1] <- insolprcpGcorL.R2.rs.mn
R2.rs.mn[3,1] <- insolprcpGL.R2.rs.mn
R2.rs.mn[4,1] <- insolprcpGcorH.R2.rs.mn
R2.rs.mn[5,1] <- insolprcpGH.R2.rs.mn
R2.rs.mn[6,1] <- MaO18insol.R2.rs.mn
R2.rs.mn[7,1] <- doleinsol.R2.rs.mn
R2.rs.mn[8,1] <- insolprcpL.R2.rs.mn
R2.rs.mn[9,1] <- insolprcpH.R2.rs.mn
R2.rs.mn[11,1] <- Hisoinsol.R2.rs.mn
R2.rs.mn[12,1] <- insoltree.R2.rs.mn
R2.rs.mn[6,3] <- MaO18prcpGcorL.R2.rs.mn
R2.rs.mn[7,3] <- doleprcpGcorL.R2.rs.mn
R2.rs.mn[11,3] <- HisoprcpGcorL.R2.rs.mn
R2.rs.mn[12,3] <- treeprcpGcorL.R2.rs.mn;
R2.rs.mn[6,4] <- MaO18prcpGL.R2.rs.mn
R2.rs.mn[7,4] <- doleprcpGL.R2.rs.mn
R2.rs.mn[11,4] <- HisoprcpGL.R2.rs.mn
R2.rs.mn[12,4] <- treeprcpGL.R2.rs.mn
R2.rs.mn[6,5] <- MaO18prcpGcorH.R2.rs.mn
R2.rs.mn[7,5] <- doleprcpGcorH.R2.rs.mn
R2.rs.mn[11,5] <- HisoprcpGcorH.R2.rs.mn
R2.rs.mn[12,5] <- treeprcpGcorH.R2.rs.mn 
R2.rs.mn[6,6] <- MaO18prcpGH.R2.rs.mn
R2.rs.mn[7,6] <- doleprcpGH.R2.rs.mn
R2.rs.mn[11,6] <- HisoprcpGH.R2.rs.mn
R2.rs.mn[12,6] <- treeprcpGH.R2.rs.mn;
R2.rs.mn[7,7] <- MaO18dole.R2.rs.mn
R2.rs.mn[8,7] <- MaO18prcpL.R2.rs.mn
R2.rs.mn[9,7] <- MaO18prcpH.R2.rs.mn
R2.rs.mn[11,7] <- MaO18Hiso.R2.rs.mn
R2.rs.mn[12,7] <- MaO18tree.R2.rs.mn
R2.rs.mn[8,8] <- doleprcpL.R2.rs.mn
R2.rs.mn[9,8] <- doleprcpH.R2.rs.mn
R2.rs.mn[11,8] <- Hisodole.R2.rs.mn
R2.rs.mn[12,8] <- doletree.R2.rs.mn
R2.rs.mn[11,9] <- HisoprcpL.R2.rs.mn
R2.rs.mn[12,9] <- treeprcpL.R2.rs.mn
R2.rs.mn[11,10] <- HisoprcpH.R2.rs.mn
R2.rs.mn[12,10] <- treeprcpH.R2.rs.mn
R2.rs.mn[12,12] <- Hisotree.p.lo
round(R2.rs.mn, 3)

slope.rs.lo <- cor.mean
slope.rs.lo[,] <- NA
slope.rs.lo[upper.tri(slope.rs.lo)] <- 0
slope.rs.lo[2,1] <- insolprcpGcorL.slope.rs.lo
slope.rs.lo[3,1] <- insolprcpGL.slope.rs.lo
slope.rs.lo[4,1] <- insolprcpGcorH.slope.rs.lo
slope.rs.lo[5,1] <- insolprcpGH.slope.rs.lo
slope.rs.lo[6,1] <- MaO18insol.slope.rs.lo
slope.rs.lo[7,1] <- doleinsol.slope.rs.lo
slope.rs.lo[8,1] <- insolprcpL.slope.rs.lo
slope.rs.lo[9,1] <- insolprcpH.slope.rs.lo
slope.rs.lo[11,1] <- Hisoinsol.slope.rs.lo
slope.rs.lo[12,1] <- insoltree.slope.rs.lo
slope.rs.lo[6,3] <- MaO18prcpGcorL.slope.rs.lo
slope.rs.lo[7,3] <- doleprcpGcorL.slope.rs.lo
slope.rs.lo[11,3] <- HisoprcpGcorL.slope.rs.lo
slope.rs.lo[12,3] <- treeprcpGcorL.slope.rs.lo;
slope.rs.lo[6,4] <- MaO18prcpGL.slope.rs.lo
slope.rs.lo[7,4] <- doleprcpGL.slope.rs.lo
slope.rs.lo[11,4] <- HisoprcpGL.slope.rs.lo
slope.rs.lo[12,4] <- treeprcpGL.slope.rs.lo
slope.rs.lo[6,5] <- MaO18prcpGcorH.slope.rs.lo
slope.rs.lo[7,5] <- doleprcpGcorH.slope.rs.lo
slope.rs.lo[11,5] <- HisoprcpGcorH.slope.rs.lo
slope.rs.lo[12,5] <- treeprcpGcorH.slope.rs.lo 
slope.rs.lo[6,6] <- MaO18prcpGH.slope.rs.lo
slope.rs.lo[7,6] <- doleprcpGH.slope.rs.lo
slope.rs.lo[11,6] <- HisoprcpGH.slope.rs.lo
slope.rs.lo[12,6] <- treeprcpGH.slope.rs.lo;
slope.rs.lo[7,7] <- MaO18dole.slope.rs.lo
slope.rs.lo[8,7] <- MaO18prcpL.slope.rs.lo
slope.rs.lo[9,7] <- MaO18prcpH.slope.rs.lo
slope.rs.lo[11,7] <- MaO18Hiso.slope.rs.lo
slope.rs.lo[12,7] <- MaO18tree.slope.rs.lo
slope.rs.lo[8,8] <- doleprcpL.slope.rs.lo
slope.rs.lo[9,8] <- doleprcpH.slope.rs.lo
slope.rs.lo[11,8] <- Hisodole.slope.rs.lo
slope.rs.lo[12,8] <- doletree.slope.rs.lo
slope.rs.lo[11,9] <- HisoprcpL.slope.rs.lo
slope.rs.lo[12,9] <- treeprcpL.slope.rs.lo
slope.rs.lo[11,10] <- HisoprcpH.slope.rs.lo
slope.rs.lo[12,10] <- treeprcpH.slope.rs.lo
slope.rs.lo[12,12] <- Hisotree.slope.rs.lo
round(slope.rs.lo, 3)

slope.rs.up <- cor.mean
slope.rs.up[,] <- NA
slope.rs.up[upper.tri(slope.rs.up)] <- 0
slope.rs.up[2,1] <- insolprcpGcorL.slope.rs.up
slope.rs.up[3,1] <- insolprcpGL.slope.rs.up
slope.rs.up[4,1] <- insolprcpGcorH.slope.rs.up
slope.rs.up[5,1] <- insolprcpGH.slope.rs.up
slope.rs.up[6,1] <- MaO18insol.slope.rs.up
slope.rs.up[7,1] <- doleinsol.slope.rs.up
slope.rs.up[8,1] <- insolprcpL.slope.rs.up
slope.rs.up[9,1] <- insolprcpH.slope.rs.up
slope.rs.up[11,1] <- Hisoinsol.slope.rs.up
slope.rs.up[12,1] <- insoltree.slope.rs.up
slope.rs.up[6,3] <- MaO18prcpGcorL.slope.rs.up
slope.rs.up[7,3] <- doleprcpGcorL.slope.rs.up
slope.rs.up[11,3] <- HisoprcpGcorL.slope.rs.up
slope.rs.up[12,3] <- treeprcpGcorL.slope.rs.up;
slope.rs.up[6,4] <- MaO18prcpGL.slope.rs.up
slope.rs.up[7,4] <- doleprcpGL.slope.rs.up
slope.rs.up[11,4] <- HisoprcpGL.slope.rs.up
slope.rs.up[12,4] <- treeprcpGL.slope.rs.up
slope.rs.up[6,5] <- MaO18prcpGcorH.slope.rs.up
slope.rs.up[7,5] <- doleprcpGcorH.slope.rs.up
slope.rs.up[11,5] <- HisoprcpGcorH.slope.rs.up
slope.rs.up[12,5] <- treeprcpGcorH.slope.rs.up 
slope.rs.up[6,6] <- MaO18prcpGH.slope.rs.up
slope.rs.up[7,6] <- doleprcpGH.slope.rs.up
slope.rs.up[11,6] <- HisoprcpGH.slope.rs.up
slope.rs.up[12,6] <- treeprcpGH.slope.rs.up;
slope.rs.up[7,7] <- MaO18dole.slope.rs.up
slope.rs.up[8,7] <- MaO18prcpL.slope.rs.up
slope.rs.up[9,7] <- MaO18prcpH.slope.rs.up
slope.rs.up[11,7] <- MaO18Hiso.slope.rs.up
slope.rs.up[12,7] <- MaO18tree.slope.rs.up
slope.rs.up[8,8] <- doleprcpL.slope.rs.up
slope.rs.up[9,8] <- doleprcpH.slope.rs.up
slope.rs.up[11,8] <- Hisodole.slope.rs.up
slope.rs.up[12,8] <- doletree.slope.rs.up
slope.rs.up[11,9] <- HisoprcpL.slope.rs.up
slope.rs.up[12,9] <- treeprcpL.slope.rs.up
slope.rs.up[11,10] <- HisoprcpH.slope.rs.up
slope.rs.up[12,10] <- treeprcpH.slope.rs.up
slope.rs.up[12,12] <- Hisotree.p.lo
round(slope.rs.up, 3)

DWp.rs.mn <- cor.mean
DWp.rs.mn[,] <- NA
DWp.rs.mn[upper.tri(DWp.rs.mn)] <- 0
DWp.rs.mn[2,1] <- insolprcpGcorL.DWp.rs.mn
DWp.rs.mn[3,1] <- insolprcpGL.DWp.rs.mn
DWp.rs.mn[4,1] <- insolprcpGcorH.DWp.rs.mn
DWp.rs.mn[5,1] <- insolprcpGH.DWp.rs.mn
DWp.rs.mn[6,1] <- MaO18insol.DWp.rs.mn
DWp.rs.mn[7,1] <- doleinsol.DWp.rs.mn
DWp.rs.mn[8,1] <- insolprcpL.DWp.rs.mn
DWp.rs.mn[9,1] <- insolprcpH.DWp.rs.mn
DWp.rs.mn[11,1] <- Hisoinsol.DWp.rs.mn
DWp.rs.mn[12,1] <- insoltree.DWp.rs.mn
DWp.rs.mn[6,3] <- MaO18prcpGcorL.DWp.rs.mn
DWp.rs.mn[7,3] <- doleprcpGcorL.DWp.rs.mn
DWp.rs.mn[11,3] <- HisoprcpGcorL.DWp.rs.mn
DWp.rs.mn[12,3] <- treeprcpGcorL.DWp.rs.mn;
DWp.rs.mn[6,4] <- MaO18prcpGL.DWp.rs.mn
DWp.rs.mn[7,4] <- doleprcpGL.DWp.rs.mn
DWp.rs.mn[11,4] <- HisoprcpGL.DWp.rs.mn
DWp.rs.mn[12,4] <- treeprcpGL.DWp.rs.mn
DWp.rs.mn[6,5] <- MaO18prcpGcorH.DWp.rs.mn
DWp.rs.mn[7,5] <- doleprcpGcorH.DWp.rs.mn
DWp.rs.mn[11,5] <- HisoprcpGcorH.DWp.rs.mn
DWp.rs.mn[12,5] <- treeprcpGcorH.DWp.rs.mn 
DWp.rs.mn[6,6] <- MaO18prcpGH.DWp.rs.mn
DWp.rs.mn[7,6] <- doleprcpGH.DWp.rs.mn
DWp.rs.mn[11,6] <- HisoprcpGH.DWp.rs.mn
DWp.rs.mn[12,6] <- treeprcpGH.DWp.rs.mn;
DWp.rs.mn[7,7] <- MaO18dole.DWp.rs.mn
DWp.rs.mn[8,7] <- MaO18prcpL.DWp.rs.mn
DWp.rs.mn[9,7] <- MaO18prcpH.DWp.rs.mn
DWp.rs.mn[11,7] <- MaO18Hiso.DWp.rs.mn
DWp.rs.mn[12,7] <- MaO18tree.DWp.rs.mn
DWp.rs.mn[8,8] <- doleprcpL.DWp.rs.mn
DWp.rs.mn[9,8] <- doleprcpH.DWp.rs.mn
DWp.rs.mn[11,8] <- Hisodole.DWp.rs.mn
DWp.rs.mn[12,8] <- doletree.DWp.rs.mn
DWp.rs.mn[11,9] <- HisoprcpL.DWp.rs.mn
DWp.rs.mn[12,9] <- treeprcpL.DWp.rs.mn
DWp.rs.mn[11,10] <- HisoprcpH.DWp.rs.mn
DWp.rs.mn[12,10] <- treeprcpH.DWp.rs.mn
DWp.rs.mn[12,12] <- Hisotree.p.lo
round(DWp.rs.mn, 3)



###########
## plots ##
###########

# rescale for plotting
  # do not centre for Hiso & toc
head(alldat.arr[,,1])
deplen <- dim(alldat.arr)[3]
alldat.arr.sc <- alldat.arr
for (d in 1:deplen) {
  alldat.arr.sc[,,d] <- as.matrix(scale(alldat.arr[,,d], scale=T, center=T))
} # end d loop
head(alldat.arr.sc[,,1])
colnames(all.dat)

# chspeleo
MaO18.up <- apply(alldat.arr.sc[,which(colnames(all.dat) == "MaO18"),], MARGIN=c(1), quantile, probs=0.975, na.rm=T)
MaO18.lo <- apply(alldat.arr.sc[,which(colnames(all.dat) == "MaO18"),], MARGIN=c(1), quantile, probs=0.025, na.rm=T)
plot(agest, MaO18.lo, type="l", lty=2, ylab="MaO18", col=NULL, xlab="", ylim=c(min(MaO18.lo, na.rm=T), max(MaO18.up, na.rm=T)))
lines(agest, MaO18.up, lty=2, col=NULL)
polygon(c(agest, rev(agest)), c(MaO18.up, rev(MaO18.lo)), col = "blue", density=50, border=NA)

# LOVECLIM regional precip
prcpL.up <- apply(alldat.arr.sc[,which(colnames(all.dat) == "prcpL"),], MARGIN=c(1), quantile, probs=0.975, na.rm=T)
prcpL.lo <- apply(alldat.arr.sc[,which(colnames(all.dat) == "prcpL"),], MARGIN=c(1), quantile, probs=0.025, na.rm=T)
plot(agest, prcpL.lo, type="l", lty=2, ylab="prcpL", col=NULL, xlab="", ylim=c(min(prcpL.lo, na.rm=T), max(prcpL.up, na.rm=T)))
lines(agest, prcpL.up, lty=2, col=NULL)
polygon(c(agest, rev(agest)), c(prcpL.up, rev(prcpL.lo)), col = "lightblue", density=50, border=NA)

# HadCM3 regional precip
prcpH.up <- apply(alldat.arr.sc[,which(colnames(all.dat) == "prcpH"),], MARGIN=c(1), quantile, probs=0.975, na.rm=T)
prcpH.lo <- apply(alldat.arr.sc[,which(colnames(all.dat) == "prcpH"),], MARGIN=c(1), quantile, probs=0.025, na.rm=T)
plot(agest, prcpH.lo, type="l", lty=2, ylab="prcpH", col=NULL, xlab="", ylim=c(min(prcpH.lo, na.rm=T), max(prcpH.up, na.rm=T)))
lines(agest, prcpH.up, lty=2, col=NULL)
polygon(c(agest, rev(agest)), c(prcpH.up, rev(prcpH.lo)), col = "lightblue", density=50, border=NA)

# dole
dole.up <- apply(alldat.arr.sc[,which(colnames(all.dat) == "dDE"),], MARGIN=c(1), quantile, probs=0.975, na.rm=T)
dole.lo <- apply(alldat.arr.sc[,which(colnames(all.dat) == "dDE"),], MARGIN=c(1), quantile, probs=0.025, na.rm=T)
plot(agest, dole.lo, type="l", lty=2, ylab="dDE", col=NULL, xlab="", ylim=c(min(dole.lo, na.rm=T), max(dole.up, na.rm=T)))
lines(agest, dole.up, lty=2, col=NULL)
polygon(c(agest, rev(agest)), c(dole.up, rev(dole.lo)), col = "green", density=50, border=NA)


# tree
tree.up <- apply(alldat.arr.sc[,which(colnames(all.dat) == "tree"),], MARGIN=c(1), quantile, probs=0.975, na.rm=T)
tree.lo <- apply(alldat.arr.sc[,which(colnames(all.dat) == "tree"),], MARGIN=c(1), quantile, probs=0.025, na.rm=T)
plot(agest, tree.lo, type="l", lty=2, ylab="tree", col=NULL, xlab="", ylim=c(min(tree.lo, na.rm=T), max(tree.up, na.rm=T)))
lines(agest, tree.up, lty=2, col=NULL)
polygon(c(agest, rev(agest)), c(tree.up, rev(tree.lo)), col = "darkgreen", density=50, border=NA)

# Hiso
Hiso.up <- apply(alldat.arr.sc[,which(colnames(all.dat) == "DHprecip"),], MARGIN=c(1), quantile, probs=0.975, na.rm=T)
Hiso.lo <- apply(alldat.arr.sc[,which(colnames(all.dat) == "DHprecip"),], MARGIN=c(1), quantile, probs=0.025, na.rm=T)
plot(agest, Hiso.lo, type="l", lty=2, ylab="DHprecip", col=NULL, xlab="", ylim=c(min(Hiso.lo, na.rm=T), max(Hiso.up, na.rm=T)))
lines(agest, Hiso.up, lty=2, col=NULL)
polygon(c(agest, rev(agest)), c(Hiso.up, rev(Hiso.lo)), col = "darkblue", density=50, border=NA)

# toc
toc.up <- apply(alldat.arr.sc[,which(colnames(all.dat) == "toc"),], MARGIN=c(1), quantile, probs=0.975, na.rm=T)
toc.lo <- apply(alldat.arr.sc[,which(colnames(all.dat) == "toc"),], MARGIN=c(1), quantile, probs=0.025, na.rm=T)
plot(agest, toc.lo, type="l", lty=2, ylab="toc", col=NULL, xlab="", ylim=c(min(toc.lo, na.rm=T), max(toc.up, na.rm=T)))
lines(agest, toc.up, lty=2, col=NULL)
polygon(c(agest, rev(agest)), c(toc.up, rev(toc.lo)), col = "brown", density=50, border=NA)


###############################
# scale original chspeleo data
chspeleo.orig.arr <- array(data=NA, dim=c(dim(chspeleo)[1], 2, iter))
for (i in 1:iter) {
  age.it <- val.it <-  rep(NA,dim(chspeleo)[1])
  for (t in 1:dim(chspeleo)[1]) {
    age.it[t] <- rtruncnorm(1, a=0, b=Inf, chspeleo$age[t], chspeleo$ageSD[t])
    val.it[t] <- rnorm(1, chspeleo$MaO18[t], chspeleo$MaO18SD[t])
  } # end t
  chspeleo.orig.arr[,,i] <- cbind(age.it, val.it)
  if (i %% itdiv==0) print(i)
} # end i

# rescale for plotting
head(chspeleo.orig.arr[,,1])
deplen <- dim(chspeleo.orig.arr)[3]
chspeleo.orig.arr.sc <- chspeleo.orig.arr
for (d in 1:deplen) {
  chspeleo.orig.arr.sc[,2,d] <- as.matrix(scale(chspeleo.orig.arr[,2,d], center=T, scale=T))
}
head(chspeleo.orig.arr.sc[,,1])
dim(chspeleo.orig.arr.sc)

MaO18.orig.mn <- apply(chspeleo.orig.arr.sc[,2,], MARGIN=c(1), mean, na.rm=T)
MaO18.orig.lo <- apply(chspeleo.orig.arr.sc[,2,], MARGIN=c(1), quantile, probs=0.025, na.rm=T)
MaO18.orig.up <- apply(chspeleo.orig.arr.sc[,2,], MARGIN=c(1), quantile, probs=0.975, na.rm=T)
MaO18.orig.age.mn <- apply(chspeleo.orig.arr.sc[,1,], MARGIN=c(1), mean, na.rm=T)
plot(MaO18.orig.age.mn, MaO18.orig.mn, pch=19, type="b", cex=0.5, lty=2, ylab="MaO18", col="blue", xlab="", ylim=c(min(MaO18.orig.lo, na.rm=T), max(MaO18.orig.up, na.rm=T)))
polygon(c(MaO18.orig.age.mn, rev(MaO18.orig.age.mn)), c(MaO18.orig.up, rev(MaO18.orig.lo)), col = "black", density=70, border=NA)


#######################################
# combine chspeleo, prcpL, prcpGcorL ##
#######################################
plot(agest, -MaO18.up, type="l", lty=2, ylab="value", col=NULL, xlab="", ylim=c(min(-MaO18.orig.up, na.rm=T), max(prcpL.up, na.rm=T)))
polygon(c(agest, rev(agest)), c(-MaO18.lo, rev(-MaO18.up)), col = "blue", density=50, border=NA)
polygon(c(agest, rev(agest)), c(prcpL.up, rev(prcpL.lo)), col = "grey9", density=50, border=NA)
lines(agest, all.dat.sc.df$prcpGcor.L, lwd=2)
lines(agest, all.dat.sc.df$insol, col="orange", lwd=2)
polygon(c(MaO18.orig.age.mn, rev(MaO18.orig.age.mn)), c(-MaO18.orig.lo, rev(-MaO18.orig.up)), col = "black", density=70, border=NA)

spelprcpLprcpGcorL.dat <- data.frame(agest, -MaO18.up, -MaO18.lo, prcpL.lo, prcpL.up, all.dat.sc.df$prcpGcor.L, all.dat.sc.df$insol)
colnames(spelprcpLprcpGcorL.dat) <- c("age","speleo.lo","speleo.up","prcpL.lo","prcpL.up", "prcpGcorL", "insol")
head(spelprcpLprcpGcorL.dat)

MaO18orig.dat <- data.frame(MaO18.orig.age.mn, -MaO18.orig.up, -MaO18.orig.lo)
colnames(MaO18orig.dat) <- c("age","speleoOrig.lo","speleoOrig.up")
head(MaO18orig.dat)


#######################################
# combine chspeleo, prcpH, prcpGcorH ##
#######################################
plot(agest, -MaO18.up, type="l", lty=2, ylab="value", col=NULL, xlab="", ylim=c(min(-MaO18.orig.up, na.rm=T), max(prcpH.up, na.rm=T)))
polygon(c(agest, rev(agest)), c(-MaO18.lo, rev(-MaO18.up)), col = "blue", density=50, border=NA)
polygon(c(agest, rev(agest)), c(prcpH.up, rev(prcpH.lo)), col = "grey9", density=50, border=NA)
lines(agest, all.dat.sc.df$prcpGcor.H, lwd=2)
polygon(c(MaO18.orig.age.mn, rev(MaO18.orig.age.mn)), c(-MaO18.orig.lo, rev(-MaO18.orig.up)), col = "black", density=70, border=NA)

spelprcpHprcpGcorH.dat <- data.frame(agest, -MaO18.up, -MaO18.lo, prcpH.up, prcpH.lo, all.dat.sc.df$prcpGcor.H, all.dat.sc.df$insol)
colnames(spelprcpHprcpGcorH.dat) <- c("age","speleo.lo","speleo.up","prcpH.up","prcpH.lo", "prcpGcorH", "insol")
head(spelprcpHprcpGcorH.dat)


# scale original dole data
head(dole)
dole.orig.arr <- array(data=NA, dim=c(dim(dole)[1], 2, iter))
for (i in 1:iter) {
  age.it <- val.it <-  rep(NA,dim(dole)[1])
  for (t in 1:dim(dole)[1]) {
    age.it[t] <- rtruncnorm(1, a=0, b=Inf, dole$age[t], dole$ageSD[t])
    val.it[t] <- rnorm(1, dole$dDE[t], dole$dDESD[t])
  } # end t
  dole.orig.arr[,,i] <- cbind(age.it, val.it)
  if (i %% itdiv==0) print(i)
} # end i

# rescale for plotting
head(dole.orig.arr[,,1])
deplen <- dim(dole.orig.arr)[3]
dole.orig.arr.sc <- dole.orig.arr
for (d in 1:deplen) {
  dole.orig.arr.sc[,2,d] <- as.matrix(scale(dole.orig.arr[,2,d], center=T, scale=T))
}
head(dole.orig.arr.sc[,,1])
dim(dole.orig.arr.sc)

dDE.orig.mn <- apply(dole.orig.arr.sc[,2,], MARGIN=c(1), mean, na.rm=T)
dDE.orig.lo <- apply(dole.orig.arr.sc[,2,], MARGIN=c(1), quantile, probs=0.025, na.rm=T)
dDE.orig.up <- apply(dole.orig.arr.sc[,2,], MARGIN=c(1), quantile, probs=0.975, na.rm=T)
dDE.orig.age.mn <- apply(dole.orig.arr.sc[,1,], MARGIN=c(1), mean, na.rm=T)
plot(dDE.orig.age.mn, dDE.orig.mn, pch=19, type="b", cex=0.5, lty=2, ylab="MaO18", col="blue", xlab="", ylim=c(min(dDE.orig.lo, na.rm=T), max(dDE.orig.up, na.rm=T)))
polygon(c(dDE.orig.age.mn, rev(dDE.orig.age.mn)), c(dDE.orig.up, rev(dDE.orig.lo)), col = "black", density=70, border=NA)


###########################
# combine chspeleo, dole ##
###########################
plot(agest, -MaO18.up, type="l", lty=2, ylab="value", col=NULL, xlab="", ylim=c(min(dole.lo, na.rm=T), max(dole.up, na.rm=T)))
polygon(c(agest, rev(agest)), c(-MaO18.lo, rev(-MaO18.up)), col = "blue", density=50, border=NA)
polygon(c(MaO18.orig.age.mn, rev(MaO18.orig.age.mn)), c(-MaO18.orig.lo, rev(-MaO18.orig.up)), col = "black", density=70, border=NA)
polygon(c(agest, rev(agest)), c(dole.up, rev(dole.lo)), col = "green", density=50, border=NA)
lines(dDE.orig.age.mn, dDE.orig.mn, pch=19, type="b", cex=0.5, col="grey2", lty=2, lwd=2)
polygon(c(dDE.orig.age.mn, rev(dDE.orig.age.mn)), c(dDE.orig.up, rev(dDE.orig.lo)), col = "darkgreen", density=70, border=NA)

speldole.dat <- data.frame(agest, -MaO18.lo, -MaO18.up, dole.up, dole.lo)
colnames(speldole.dat) <- c("age","speleo.up","speleo.lo","dole.up","dole.lo")
head(speldole.dat)

doleorig.dat <- data.frame(dDE.orig.age.mn, dDE.orig.mn, dDE.orig.up, dDE.orig.lo)
colnames(doleorig.dat) <- c("age","doleOrig.mn","doleOrig.up","doleOrig.lo")
head(doleorig.dat)


# scale original tree data
head(tree)
tree.orig.arr <- array(data=NA, dim=c(dim(tree)[1], 2, iter))
for (i in 1:iter) {
  age.it <- val.it <-  rep(NA,dim(tree)[1])
  for (t in 1:dim(tree)[1]) {
    age.it[t] <- rtruncnorm(1, a=0, b=Inf, tree$age[t], tree$ageSD[t]/2)
    val.it[t] <- rnorm(1, tree$tree[t], tree$treeSD[t])
  } # end t
  tree.orig.arr[,,i] <- cbind(age.it, val.it)
  if (i %% itdiv==0) print(i)
} # end i

# rescale for plotting
head(tree.orig.arr[,,1])
deplen <- dim(tree.orig.arr)[3]
tree.orig.arr.sc <- tree.orig.arr
for (d in 1:deplen) {
  tree.orig.arr.sc[,2,d] <- as.matrix(scale(tree.orig.arr[,2,d], center=T, scale=T))
}
head(tree.orig.arr.sc[,,1])
dim(tree.orig.arr.sc)

tree.orig.mn <- apply(tree.orig.arr.sc[,2,], MARGIN=c(1), mean, na.rm=T)
tree.orig.lo <- apply(tree.orig.arr.sc[,2,], MARGIN=c(1), quantile, probs=0.025, na.rm=T)
tree.orig.up <- apply(tree.orig.arr.sc[,2,], MARGIN=c(1), quantile, probs=0.975, na.rm=T)
tree.orig.age.mn <- apply(tree.orig.arr.sc[,1,], MARGIN=c(1), mean, na.rm=T)
plot(tree.orig.age.mn,tree.orig.mn, pch=19, type="b", cex=0.5, lty=2, ylab="tree", col="green", xlab="", ylim=c(min(tree.orig.lo, na.rm=T), max(tree.orig.up, na.rm=T)))
polygon(c(tree.orig.age.mn, rev(tree.orig.age.mn)), c(tree.orig.up, rev(tree.orig.lo)), col = "darkgreen", density=70, border=NA)

# scale original Hiso data
head(Hiso)
Hiso.orig.arr <- array(data=NA, dim=c(dim(Hiso)[1], 2, iter))
for (i in 1:iter) {
  age.it <- val.it <-  rep(NA,dim(Hiso)[1])
  for (t in 1:dim(Hiso)[1]) {
    age.it[t] <- rtruncnorm(1, a=0, b=Inf, Hiso$age[t], Hiso$ageSD[t]/2)
    val.it[t] <- rnorm(1, Hiso$DHprecip[t], Hiso$DHprecipSD[t])
  } # end t
  Hiso.orig.arr[,,i] <- cbind(age.it, val.it)
  if (i %% itdiv==0) print(i)
} # end i

# rescale for plotting
head(Hiso.orig.arr[,,1])
deplen <- dim(Hiso.orig.arr)[3]
Hiso.orig.arr.sc <- Hiso.orig.arr
for (d in 1:deplen) {
  Hiso.orig.arr.sc[,2,d] <- as.matrix(scale(Hiso.orig.arr[,2,d], center=T, scale=T))
}
head(Hiso.orig.arr.sc[,,1])
dim(Hiso.orig.arr.sc)

Hiso.orig.mn <- apply(Hiso.orig.arr.sc[,2,], MARGIN=c(1), mean, na.rm=T)
Hiso.orig.lo <- apply(Hiso.orig.arr.sc[,2,], MARGIN=c(1), quantile, probs=0.025, na.rm=T)
Hiso.orig.up <- apply(Hiso.orig.arr.sc[,2,], MARGIN=c(1), quantile, probs=0.975, na.rm=T)
Hiso.orig.age.mn <- apply(Hiso.orig.arr.sc[,1,], MARGIN=c(1), mean, na.rm=T)
plot(Hiso.orig.age.mn,Hiso.orig.mn, pch=19, type="b", cex=0.5, lty=2, ylab="Hiso", col="green", xlab="", ylim=c(min(Hiso.orig.lo, na.rm=T), max(Hiso.orig.up, na.rm=T)))
polygon(c(Hiso.orig.age.mn, rev(Hiso.orig.age.mn)), c(Hiso.orig.up, rev(Hiso.orig.lo)), col = "darkgreen", density=70, border=NA)


#######################
# combine tree, Hiso ##
#######################
plot(agest, tree.lo, type="b", lty=2, ylab="value", col=NULL, xlab="", ylim=c(min(-Hiso.orig.up, na.rm=T), max(toc.up, na.rm=T)))
polygon(c(agest, rev(agest)), c(tree.up, rev(tree.lo)), col = "darkgreen", density=50, border=NA)
polygon(c(tree.orig.age.mn, rev(tree.orig.age.mn)), c(tree.orig.up, rev(tree.orig.lo)), col = "black", density=90, border=NA)
#polygon(c(agest, rev(agest)), c(-Hiso.lo, rev(-Hiso.up)), col = "red", density=40, border=NA)
lines(Hiso.orig.age.mn,-Hiso.orig.mn, pch=19, type="b", cex=0.5, lty=1, lwd=2, ylab="Hiso", col="black")
polygon(c(Hiso.orig.age.mn, rev(Hiso.orig.age.mn)), c(-Hiso.orig.lo, rev(-Hiso.orig.up)), col = "purple4", density=90, border=NA)

treeHiso.dat <- data.frame(agest, tree.up, tree.lo, -Hiso.lo, -Hiso.up)
colnames(treeHiso.dat) <- c("age","tree.up","tree.lo","Hiso.up","Hiso.lo")
head(treeHiso.dat)

treeorig.dat <- data.frame(tree.orig.age.mn, tree.orig.mn, tree.orig.up, tree.orig.lo)
colnames(treeorig.dat) <- c("age","treeOrig.mn","treeOrig.up","treeOrig.lo")
head(treeorig.dat)

Hisoorig.dat <- data.frame(Hiso.orig.age.mn, -Hiso.orig.mn, -Hiso.orig.lo, -Hiso.orig.up)
colnames(Hisoorig.dat) <- c("age","HisoOrig.mn","HisoOrig.up","HisoOrig.lo")
head(Hisoorig.dat)


#########################
# combine speleo, tree ##
#########################
plot(agest, tree.lo, type="b", lty=2, ylab="value", col=NULL, xlab="", ylim=c(min(-MaO18.up, na.rm=T), max(tree.up, na.rm=T)))
polygon(c(agest, rev(agest)), c(tree.up, rev(tree.lo)), col = "darkgreen", density=50, border=NA)
polygon(c(tree.orig.age.mn, rev(tree.orig.age.mn)), c(tree.orig.up, rev(tree.orig.lo)), col = "black", density=90, border=NA)
polygon(c(agest, rev(agest)), c(-MaO18.lo, rev(-MaO18.up)), col = "blue", density=50, border=NA)
polygon(c(MaO18.orig.age.mn, rev(MaO18.orig.age.mn)), c(-MaO18.orig.lo, rev(-MaO18.orig.up)), col = "black", density=70, border=NA)


#########################
# combine speleo, Hiso ##
#########################
plot(agest, -Hiso.up, type="b", lty=2, ylab="value", col=NULL, xlab="", ylim=c(min(-Hiso.orig.up, na.rm=T), max(-Hiso.orig.lo, na.rm=T)))
#polygon(c(agest, rev(agest)), c(-Hiso.lo, rev(-Hiso.up)), col = "red", density=40, border=NA)
polygon(c(agest, rev(agest)), c(-MaO18.lo, rev(-MaO18.up)), col = "blue", density=50, border=NA)
polygon(c(MaO18.orig.age.mn, rev(MaO18.orig.age.mn)), c(-MaO18.orig.lo, rev(-MaO18.orig.up)), col = "black", density=70, border=NA)
lines(Hiso.orig.age.mn,-Hiso.orig.mn, pch=19, type="b", cex=0.5, lty=1, lwd=2, ylab="Hiso", col="black")
polygon(c(Hiso.orig.age.mn, rev(Hiso.orig.age.mn)), c(-Hiso.orig.lo, rev(-Hiso.orig.up)), col = "purple4", density=90, border=NA)

speleotreeHiso.dat <- data.frame(agest, apply(data.frame(-MaO18.up, -MaO18.lo), MARGIN=1, mean, na.rm=T),
                                 -MaO18.lo, -MaO18.up, apply(data.frame(tree.up, tree.lo), MARGIN=1, mean, na.rm=T), tree.up, tree.lo,
                                 apply(data.frame(-Hiso.up, -Hiso.lo), MARGIN=1, mean, na.rm=T), -Hiso.lo, -Hiso.up)
colnames(speleotreeHiso.dat) <- c("age","MaO18.mn","MaO18.up","MaO18.lo","tree.mn","tree.up","tree.lo","Hiso.mn","Hiso.up","Hiso.lo")
head(speleotreeHiso.dat)


# scale original toc data
head(toc)
toc.orig.arr <- array(data=NA, dim=c(dim(toc)[1], 2, iter))
for (i in 1:iter) {
  age.it <- val.it <-  rep(NA,dim(toc)[1])
  for (t in 1:dim(toc)[1]) {
    age.it[t] <- rtruncnorm(1, a=0, b=Inf, toc$age[t], toc$ageSD[t]/2)
    val.it[t] <- rnorm(1, toc$toc[t], toc$tocSD[t])
  } # end t
   toc.orig.arr[,,i] <- cbind(age.it, val.it)
  if (i %% itdiv==0) print(i)
} # end i

# rescale for plotting
head(toc.orig.arr[,,1])
deplen <- dim(toc.orig.arr)[3]
toc.orig.arr.sc <- toc.orig.arr
for (d in 1:deplen) {
  toc.orig.arr.sc[,2,d] <- as.matrix(scale(toc.orig.arr[,2,d], center=F, scale=T))
}
head(toc.orig.arr.sc[,,1])
dim(toc.orig.arr.sc)

toc.orig.mn <- apply(toc.orig.arr[,2,], MARGIN=c(1), mean, na.rm=T)
toc.orig.lo <- apply(toc.orig.arr[,2,], MARGIN=c(1), quantile, probs=0.025, na.rm=T)
toc.orig.up <- apply(toc.orig.arr[,2,], MARGIN=c(1), quantile, probs=0.975, na.rm=T)
toc.orig.age.mn <- apply(toc.orig.arr[,1,], MARGIN=c(1), mean, na.rm=T)
plot(toc.orig.age.mn,toc.orig.mn, pch=19, type="b", cex=0.5, lty=2, ylab="toc", col="black", xlab="", ylim=c(min(toc.orig.lo, na.rm=T), max(toc.orig.up, na.rm=T)))
polygon(c(toc.orig.age.mn, rev(toc.orig.age.mn)), c(toc.orig.up, rev(toc.orig.lo)), col = "grey", density=70, border=NA)

# toc not scaled
tocNS.up <- apply(alldat.arr[,which(colnames(all.dat) == "toc"),], MARGIN=c(1), quantile, probs=0.975, na.rm=T)
tocNS.lo <- apply(alldat.arr[,which(colnames(all.dat) == "toc"),], MARGIN=c(1), quantile, probs=0.025, na.rm=T)


#############
# toc plot ##
#############
plot(agest, tocNS.lo, type="l", lty=2, ylab="toc", col=NULL, xlab="", ylim=c(min(tocNS.lo, na.rm=T), max(toc.orig.up, na.rm=T)))
lines(agest, tocNS.up, lty=2, col=NULL)
polygon(c(agest, rev(agest)), c(tocNS.up, rev(tocNS.lo)), col = "brown", density=50, border=NA)
lines(toc.orig.age.mn,toc.orig.mn, pch=19, type="b", lwd=2, cex=0.5, lty=2)
polygon(c(toc.orig.age.mn, rev(toc.orig.age.mn)), c(toc.orig.up, rev(toc.orig.lo)), col = "grey", density=70, border=NA)

toc.dat <- data.frame(agest, apply(data.frame(tocNS.up, tocNS.lo), MARGIN=1, mean, na.rm=T), tocNS.up, tocNS.lo)
colnames(toc.dat) <- c("age","toc.mn","toc.up","toc.lo")
head(toc.dat)

tocOrig.dat <- data.frame(toc.orig.age.mn, apply(data.frame(toc.orig.up, toc.orig.lo), MARGIN=1, mean, na.rm=T), toc.orig.up, toc.orig.lo)
colnames(tocOrig.dat) <- c("age","tocOrig.mn","tocOrig.up","tocOrig.lo")
head(tocOrig.dat)


#####################################################
# compare LOVECLIM northern Australia & South Asia ##
#####################################################

#######################################
# combine prcpL, prcpSA ##
#######################################
LOVECLIM.NA.SAcor.mean <- apply(LOVECLIM.NA.SAcor.arr, MARGIN=c(1,2), mean, na.rm=T)
round(LOVECLIM.NA.SAcor.mean, 3)

LOVECLIM.NA.SAcor.lo <- apply(LOVECLIM.NA.SAcor.arr, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
round(LOVECLIM.NA.SAcor.lo, 3)

LOVECLIM.NA.SAcor.up <- apply(LOVECLIM.NA.SAcor.arr, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
round(LOVECLIM.NA.SAcor.up, 3)

print(round(c(prcpLprcpSA.slope.lo, prcpLprcpSA.slope.mn, prcpLprcpSA.slope.up), 3))
print(round(c(prcpLprcpSA.R2.lo, prcpLprcpSA.R2.mn, prcpLprcpSA.R2.up), 3))
print(round(c(prcpLprcpSA.p.lo, prcpLprcpSA.p.mn, prcpLprcpSA.p.up), 3))

print(round(c(prcpLprcpSA.slope.rs.lo, prcpLprcpSA.slope.rs.mn, prcpLprcpSA.slope.rs.up), 3))
print(round(c(prcpLprcpSA.R2.rs.lo, prcpLprcpSA.R2.rs.mn, prcpLprcpSA.R2.rs.up), 3))
print(round(c(prcpLprcpSA.p.rs.lo, prcpLprcpSA.p.rs.mn, prcpLprcpSA.p.rs.up), 3))
print(round(c(prcpLprcpSA.DWp.rs.lo, prcpLprcpSA.DWp.rs.mn, prcpLprcpSA.DWp.rs.up), 3))

# rescale for plotting
head(LOVECLIM.NA.SAdat.arr[,,1])
deplen <- dim(LOVECLIM.NA.SAdat.arr)[3]
LOVECLIM.NA.SAdat.arr.sc <- LOVECLIM.NA.SAdat.arr
for (d in 1:deplen) {
  LOVECLIM.NA.SAdat.arr.sc[,,d] <- as.matrix(scale(LOVECLIM.NA.SAdat.arr[,,d], scale=T, center=T))
} # end d loop
head(LOVECLIM.NA.SAdat.arr[,,1])
colnames(LOVECLIM.NA.SA.dat)

# LOVECLIM South Asia
prcpSA.up <- apply(LOVECLIM.NA.SAdat.arr.sc[,which(colnames(LOVECLIM.NA.SA.dat) == "prcpSA"),], MARGIN=c(1), quantile, probs=0.975, na.rm=T)
prcpSA.lo <- apply(LOVECLIM.NA.SAdat.arr.sc[,which(colnames(LOVECLIM.NA.SA.dat) == "prcpL"),], MARGIN=c(1), quantile, probs=0.025, na.rm=T)
plot(agest, prcpSA.lo, type="l", lty=2, ylab="prcpL", col=NULL, xlab="", ylim=c(min(prcpSA.lo, na.rm=T), max(prcpSA.up, na.rm=T)))
lines(agest, prcpSA.up, lty=2, col=NULL)
polygon(c(agest, rev(agest)), c(prcpSA.up, rev(prcpSA.lo)), col = "lightblue", density=50, border=NA)

LOVECLIMprpNASA.dat <- data.frame(agest, apply(data.frame(prcpSA.up, prcpSA.lo), MARGIN=1, mean, na.rm=T), prcpSA.up, prcpSA.lo,
                                  apply(data.frame(prcpL.up, prcpL.lo), MARGIN=1, mean, na.rm=T), prcpL.up, prcpL.lo)
colnames(LOVECLIMprpNASA.dat) <- c("age","prcpSA.mn","prcpSA.up","prcpSA.lo","prcpL.mn","prcpL.up","prcpL.lo")
head(LOVECLIMprpNASA.dat)

plot(agest, prcpSA.lo, type="l", lty=2, ylab="value", col=NULL, xlab="", ylim=c(min(prcpSA.lo, na.rm=T), max(prcpSA.up, na.rm=T)))
polygon(c(agest, rev(agest)), c(prcpSA.up, rev(prcpSA.lo)), col = "blue", density=50, border=NA)
polygon(c(agest, rev(agest)), c(prcpL.up, rev(prcpL.lo)), col = "pink", density=50, border=NA)
lines(LOVECLIMprpNASA.dat$age, LOVECLIMprpNASA.dat$prcpSA.mn, lty=2, lwd=2, col="blue")
lines(LOVECLIMprpNASA.dat$age, LOVECLIMprpNASA.dat$prcpL.mn, lty=2, lwd=2, col="red")

plot(LOVECLIMprpNASA.dat$prcpSA.mn, LOVECLIMprpNASA.dat$prcpL.mn, pch=19, cex=0.7, xlab="South Asia", ylab="northern Australia")
LOVECLIMprpNASA.fit <- lm( LOVECLIMprpNASA.dat$prcpL.mn ~  LOVECLIMprpNASA.dat$prcpSA.mn)
abline(LOVECLIMprpNASA.fit, lty=2, lwd=2, col="red")
cor(na.omit(LOVECLIMprpNASA.dat[,c(2,5)]), method="spearman")
