## simulation testing whether Girraween tree pollen counts higher during Heinrich events
## Corey Bradshaw
## November 2024

# libraries
library(ggplot2)
library(truncnorm)

##########
## DATA ##
##########

# % tree
# % of tree pollen in total pollen (so trees:grass)
# with a flat assumption of a 2% counting error (interpretation is more trees = wetter climate),
# while age uncertainty is 2 standard deviations provided by rbacon model
tree <- read.csv("tree.csv", header=T)
tree$treeSD <- sqrt(((tree$tree/100)*(1 - (tree$tree/100)))/tree$pollenN)
head(tree)

## Heinrich event dates
Heinrich0 <- c(12,12.5); Heinrich1 <- c(15.5,18); Heinrich2 <- c(24,26); Heinrich3 <- c(30,32)
Heinrich4 <- c(38,40); Heinrich5 <- c(47,49); Heinrich5a <- c(54,56); Heinrich6 <- c(60,63)
Heinrich7 <- c(66,68); Heinrich8 <- c(85.5,88); Heinrich9 <- c(103,105); Heinrich10 <- c(109,110)
Heinrich11 <- c(129,136); Heinrich12 <- c(139.5,140.5)

Heinrich.dates <- rbind(Heinrich0, Heinrich1, Heinrich2, Heinrich3, Heinrich4, Heinrich5, Heinrich5a,
                        Heinrich6, Heinrich7, Heinrich8, Heinrich9, Heinrich10, Heinrich11, Heinrich12)
colnames(Heinrich.dates) <- c("en","st")
rownames(Heinrich.dates) <- c(paste("H",seq(0,5,1),sep=""), "H5a", paste("H",seq(6,12,1),sep=""))
Heinrich.dates <- as.data.frame(Heinrich.dates)
Heinrich.dates$dur <- Heinrich.dates$st - Heinrich.dates$en
Heinrich.dates

#########################################
## Heinrich event sensitivity analysis ##
#########################################

## create randomly sampled tree ages & pollen scores
reps <- 2000
tree.dat.rnd <- data.frame(age=rep(NA,1),poll=rep(NA,1))
for (t in 1:dim(tree)[1]) {
  tr.age.rnd <- rtruncnorm(reps, a=0, b=Inf, mean=tree$age[t], sd=tree$ageSD[t])
  tr.poll.rnd <- rtruncnorm(reps, a=0, b=Inf, mean=tree$tree[t], sd=tree$treeSD[t])
  dat.it <- data.frame(age=tr.age.rnd,poll=tr.poll.rnd)
  tree.dat.rnd <- rbind(tree.dat.rnd, dat.it)
} # end t
tree.dat.rnd <- tree.dat.rnd[-1,]
dim(tree.dat.rnd)
head(tree.dat.rnd)

Heinrich.tree <- Heinrich.dates
nonHeinrich.dates <-
  data.frame(lab=c("NH0","NH1","NH2","NH3","NH4","NH5","NH5a","NH6","NH7","NH8","NH9","NH10","NH11","NH12","NH13"),
             en=c(11.5, 12.5, 18, 26, 32, 40, 49, 56, 63, 68, 88, 105, 110, 136, 140.5),
             st=c(12, 15.5, 24, 30, 38, 47, 54, 60, 66, 85.5, 103, 109, 129, 139.5, 141.5))
nonHeinrich.tree <- nonHeinrich.dates

## select randomised tree pollen values inside/outside Heinrich events
Heinrich.tree$pollmed <- rep(NA,dim(Heinrich.dates)[1])
for (h in 1:dim(Heinrich.dates)[1]) {
  h.subs <- which(tree.dat.rnd$age >= Heinrich.dates[h,1] & tree.dat.rnd$age <= Heinrich.dates[h,2])
  Heinrich.tree$pollmed[h] <- median(tree.dat.rnd[h.subs,]$poll)
  
} # end h

nonHeinrich.tree$pollmed <- rep(NA,dim(nonHeinrich.dates)[1])
for (n in 1:dim(nonHeinrich.dates)[1]) {
  nh.subs <- which(tree.dat.rnd$age >= nonHeinrich.dates[n,2] & tree.dat.rnd$age <= nonHeinrich.dates[n,3])
  nonHeinrich.tree$pollmed[n] <- median(tree.dat.rnd[nh.subs,]$poll)
} # end n

par(mfrow=c(1,2))
plot(apply(Heinrich.tree[,2:3],1,median), Heinrich.tree$pollmed, pch=19, col="red",
     xlab="Heinrich event mid age (ka)", ylab="%tree pollen", ylim=c(0,60))

plot(apply(nonHeinrich.tree[,3:4],1,median), nonHeinrich.tree$pollmed, pch=19, col="black",
     xlab="non Heinrich event mid age (ka)", ylab="%tree pollen", ylim=c(0,60))
par(mfrow=c(1,1))


H2 <- data.frame(lab=row.names(Heinrich.tree), Heinrich.tree)
NH2 <- data.frame(nonHeinrich.tree[,1:3], dur=(nonHeinrich.tree$st - nonHeinrich.tree$en),
                  pollmed=nonHeinrich.tree$pollmed)
row.names(NH2) <- nonHeinrich.tree$lab

all.tree <- rbind(H2,NH2)
all.tree.sort <- all.tree[order(all.tree$en,decreasing=F),]
all.tree.sort

## calculate sampling windows either side of Heinrich event
yng.side <- old.side <- rep(NA,dim(H2)[1])
for (h in 1:dim(H2)[1]) {
  yng.side[h] <- H2[h,2] - H2[h,4]
  old.side[h] <- H2[h,3] + H2[h,4]
} # end h

event.comb <- rbind(data.frame(lab=paste(H2$lab,"yng",sep=""), en=yng.side, st=H2$en, type="NH"),
                    data.frame(lab=paste(H2$lab,"old",sep=""), en=H2$st,st=old.side, type="NH"),
                    data.frame(lab=H2$lab, en=H2$en, st=H2$st, type="H"))
event.sort <- event.comb[order(event.comb$en,decreasing=F),]

# correct for H5a & H6 overlap
event.sort[which(event.sort$lab=="H5aold"), 3] <- 57.5
event.sort[which(event.sort$lab=="H6yng"), 2] <- 57.5

# correct for H6 & H7 overlap
event.sort[which(event.sort$lab=="H6old"), 3] <- 65
event.sort[which(event.sort$lab=="H7yng"), 2] <- 65

# correct for H11 & H12 overlap
event.sort[which(event.sort$lab=="H11old"), 3] <- 138.5
event.sort

##########################################################
## assume uncertainty in onset & end of Heinrich events ##
##########################################################

uncertp <- 0.01 # proportion
Heinrich.uncert <- Heinrich.dates
Heinrich.uncert$enUncrt <-  Heinrich.uncert$en * uncertp
Heinrich.uncert$stUncrt <-  Heinrich.uncert$st * uncertp
H.uncrt <- data.frame(lab=row.names(Heinrich.uncert), Heinrich.uncert)

# select randomised tree pollen values inside/outside Heinrich events with new yng/old overlaps
rsamp <- 100 # number of random samples within the interval

iter <- 10000
itdiv <- iter/100

HgtNHsums <- rep(NA,iter)
for (s in 1:iter) {
  
  # generate resampled Heinrich onset times
  st.it <- rnorm(length(H.uncrt$st), H.uncrt$st, H.uncrt$stUncrt)
  en.it <- st.it - H.uncrt$dur
  H.it <- data.frame(lab=row.names(H.uncrt), en=en.it, st=st.it, dur=(st.it-en.it))
  
  # establish sampling windows before and after Heinrichs
  yng.side <- old.side <- rep(NA,dim(H.it)[1])
  for (h in 1:dim(H.it)[1]) {
    yng.side[h] <- H.it[h,2] - H.it[h,4]
    old.side[h] <- H.it[h,3] + H.it[h,4]
  } # end h
  
  # combine Heinrich and non-Heinrich windows
  event.comb <- rbind(data.frame(lab=paste(H.it$lab,"yng",sep=""), en=yng.side, st=H.it$en, type="NH"),
                      data.frame(lab=paste(H.it$lab,"old",sep=""), en=H.it$st,st=old.side, type="NH"),
                      data.frame(lab=H.it$lab, en=H.it$en, st=H.it$st, type="H"))
  event.it <- event.comb[order(event.comb$en,decreasing=F),]
  
  # event overlap correction
  for (o in 2:(dim(event.it)[1])-1) {
    event.it[o,3] <- ifelse(event.it[o,3] > event.it[o+1,2],
                            median(c(event.it[o,3],event.it[o+1,2])),
                            event.it[o,3])
    event.it[o+1,2] <- ifelse(event.it[o,3] < event.it[o+1,2],
                              event.it[o+1,2], event.it[o,3])
  } # end s
  
  # remove negative durations
  negdur.sub <- which(event.it$st - event.it$en < 0)
  if (length(negdur.sub) > 0) {
    event.it <- event.it[-negdur.sub,]
  } # end if
  
  # redo overlap correction
  overlap.sub <- which(event.it[2:(length(event.it[,3])),2] - event.it[1:(length(event.it[,3])-1),3] < 0)
  
  # remove overlaps while loop
  while (length(overlap.sub) > 0) {
    event.it <- event.it[-(overlap.sub+1),]
    overlap.sub <- which(event.it[2:(length(event.it[,3])),2] - event.it[1:(length(event.it[,3])-1),3] < 0)
  } # end while
  
  # determine %tree values per event
  pollMd <- HgtNH <- rep(NA,dim(event.it)[1])
  for (e in 1:dim(event.it)[1]) {
    e.subs <- which(tree.dat.rnd$age >= event.it[e,2] & tree.dat.rnd$age <= event.it[e,3])
    if (length(e.subs) > 0) {
      e.subs.rsamp <- sample(e.subs, rsamp, replace=T)
      pollMd[e] <- median(tree.dat.rnd[e.subs.rsamp,]$poll)} # end if
    if (length(e.subs) == 0) {
      pollMd[e] <- NA
    }
    
    # determine which event.it rows are H & NH
    Hsubs <- which(event.it$type == "H")
    NHsubs <- which(event.it$type == "NH")
    i <- 1
    NHsubs2 <- NHsubs[1:(i+1)==(i+1)]
    
    if ((e %in% NHsubs2) && (e > 2)) {
      HgtNH[e] <- ifelse((pollMd[e-1] - pollMd[e] > 0) | (pollMd[e-1] - pollMd[e-2] > 0), 1, 0)
    } # end if
    
  } # end e
  
  HgtNHsums[s] <- ifelse(sum(HgtNH, na.rm=T) > rbinom(1,length(NHsubs2),0.5), 0, 1)
  
  if (s %% itdiv==0) print(paste("iter = ", s, ";", " p = ",round(sum(HgtNHsums, na.rm=T)/s, 4), sep="")) # running probability
} # end s

prNotRnd <- sum(HgtNHsums)/iter
prNotRnd # probability not random that higher %tree in H vs. NH
