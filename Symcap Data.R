setwd("~/Symcap")
library(data.table)
library(lsmeans)
library(devtools)
library(plyr)
library(reshape2)
library(popbio)
library(RgoogleMaps)
library(plotrix)
library(zoo)
library(rgdal)
library(car)
library(scales)
library(png)
library(pixmap)
library(ecodist)
library(cluster)
library(fpc)
library(clustsig)
library(foreign)
library(nnet)
library(ggplot2)
library(mlogit)

#import collection data
Coral_Data <- read.csv("Coral_Collection.csv")
Coral_Data$Depth..m. <- as.numeric(as.character(Coral_Data$Depth..m.))

# Importing qPCR data
source_url("https://raw.githubusercontent.com/jrcunning/steponeR/master/steponeR.R")
Mcap.plates <- list.files(path = "qPCR_data", pattern = "txt$", full.names = T)
Mcap.plates
Mcap <- steponeR(files = Mcap.plates, delim="\t",
                 target.ratios=c("C.D"), 
                 fluor.norm = list(C=2.26827, D=0), 
                 copy.number=list(C=33, D=3))
Mcap <- Mcap$result

#Remove positive controls
Mcap <- Mcap[grep("+", Mcap$Sample.Name, fixed=T, invert = T), ]

#Remove NTCs
Mcap <- Mcap[grep("NTC", Mcap$Sample.Name, fixed = T, invert = T), ]

#Remove PCTs
Mcap <- Mcap[grep("PCT", Mcap$Sample.Name, fixed = T, invert = T), ]

#Replace "Sample.Name" column with "Colony"
colnames(Mcap)[which(colnames(Mcap)=="Sample.Name")] <- "Colony"

#Remove failed samples
Mcap$fail <- ifelse(Mcap$C.reps < 2 & Mcap$D.reps < 2, TRUE, FALSE)
fails <- Mcap[Mcap$fail==TRUE, ]
Mcap <- Mcap[which(Mcap$fail==FALSE),]

#When C or D detected in only 1 sample, set ratio to +/- Inf
Mcap$C.D[which(Mcap$C.reps<2)] <- -Inf
Mcap$C.D[which(Mcap$D.reps<2)] <- Inf

#Order Mcap by Colony
Mcap <- Mcap[with(Mcap, order(Colony)), ]

#Proportion C
Mcap$propC <- Mcap$C.D / (Mcap$C.D+1)

#Proportion D
Mcap$propD <- 1-Mcap$propC

#Change NA to +/- Inf
Mcap$propD[which(Mcap$C.D==-Inf)] <-1
Mcap$propC[which(Mcap$C.D==-Inf)] <-0
Mcap$propD[which(Mcap$C.D==Inf)] <-0
Mcap$propC[which(Mcap$C.D==Inf)] <-1

#Dominant Symbiont type
Mcap$Dom <- ifelse(Mcap$propC>Mcap$propD, "C", "D")

#Merge datasets
Symcap<-merge(Coral_Data, Mcap, by="Colony", all=T)

#Order Symcap by Colony
Symcap <- Symcap[with(Symcap, order(Colony)), ]

#Identify Symbiont clades present
Symcap$Mix <- factor(ifelse(Symcap$propC>Symcap$propD, ifelse(Symcap$propD!=0, "CD", "C"), ifelse(Symcap$propD>Symcap$propC, ifelse(Symcap$propC!=0, "DC", "D"), NA)), levels = c("C", "CD", "DC", "D"))

#Leeward and Windward to Slope
Symcap$Reef.Area <- ifelse(Symcap$Reef.Area!="Top", yes = "Slope", no = "Top")

#Adjust depth by sea level
JuneTide=read.csv("Station_1612480_tide_ht_20160601-20160630.csv")
JulyTide=read.csv("Station_1612480_tide_ht_20160701-20160731.csv")
AugustTide=read.csv("Station_1612480_tide_ht_20160801-20160812.csv")
Tide<-rbind(JuneTide, JulyTide, AugustTide)
Tide$Time <- as.POSIXct(Tide$TimeUTC, format="%Y-%m-%d %H:%M:%S", tz="UTC")
attributes(Tide$Time)$tzone <- "Pacific/Honolulu"
plot(Tide$Time, Tide$TideHT, type = "l", xlab = "Date", ylab = "Tide Height, meters")
with(Tide[Tide$Time > "2016-06-07 00:00:00" & Tide$Time < "2016-06-09 00:00:00", ], 
     plot(Time, TideHT, type="l"))
Symcap$Time2 <- as.POSIXct(paste(as.character(Symcap$Date), as.character(Symcap$Time)),
                                format="%m/%d/%y %H:%M", tz="Pacific/Honolulu")
Symcap$Time=Symcap$Time2

# Add estimated times for missing values
which(is.na(Symcap$Time))
Symcap$Time[170] <- as.POSIXct("2016-06-14 12:07:00")
Symcap$Time[177] <- as.POSIXct("2016-06-14 12:20:00")
Symcap$Time[178] <- as.POSIXct("2016-06-14 12:22:00")
Symcap$Time[180] <- as.POSIXct("2016-06-14 13:08:00")
Symcap$Time[187] <- as.POSIXct("2016-06-14 12:42:00")
Symcap$Time[188] <- as.POSIXct("2016-06-14 12:44:00")
Symcap$Time[206] <- as.POSIXct("2016-06-16 13:10:00")
Symcap$Time[208] <- as.POSIXct("2016-06-16 13:24:00")
Symcap$Time[211] <- as.POSIXct("2016-06-16 12:37:00")
Symcap$Time[218] <- as.POSIXct("2016-06-16 12:27:00")
Symcap$Time[448] <- as.POSIXct("2016-07-16 13:32:00")

Round6 <- function (times)  {
  x <- as.POSIXlt(times)
  mins <- x$min
  mult <- mins %/% 6
  remain <- mins %% 6
  if(remain > 3L) {
    mult <- mult + 1
  } else {
    x$min <- 6 * mult
  }
  x <- as.POSIXct(x)
  x
}

Symcap$Time.r <- Round6(Symcap$Time)
Tide$Time.r <- Tide$Time
merged<-merge(Symcap, Tide, by="Time.r", all.x=T)
merged$newDepth <- merged$Depth..m.- merged$TideHT

#Chi Squared test for independence
merged$Reef.Area <- ifelse(merged$Reef.Area!="Top", yes = "Slope", no = "Top")
merged$DepthInt <- cut(merged$Depth..m., breaks = 0:13)
merged$Dominant <- ifelse(merged$Dom=="C", 0, 1)
merged$Dominant2 <- ifelse(merged$Dom=="C", 1, 0)
results=table(merged$Dominant2, merged$DepthInt)
results
props <- prop.table(results, margin = 2)
par(mar=c(4, 4, 2, 6), lwd = 0.25)
barplot(props[,1:11], col = c(alpha("red", 0.25), alpha("blue", 0.25)), 
        xlab = "", ylab = "",
        space = 0, xaxs="i", yaxs="i", axisnames = FALSE)
par(lwd=1)
legend("topright", legend=c("D", "C"), fill=c("red", "blue"), inset = c(-.2, 0), xpd = NA)
par(new = T)
par(mar=c(4, 4, 2, 6))
results=glm(Dominant~Depth..m., family = "binomial", data = merged)
fitted <- predict(results, newdata = list(Depth..m.=seq(0,11,0.1)), type = "response")
plot(fitted~seq(0,11,0.1), xaxs="i", yaxs="i", xlim=c(0,11), ylim=c(0,1), type="l", lwd = 3)

merged$Color <- ifelse(merged$Color.Morph=="Orange", 0, 1)
results=table(merged$Color, merged$DepthInt)
results
props <- prop.table(results, margin = 2)
par(mar=c(4, 4, 2, 6), lwd = 0.25)
barplot(props[,1:11], col = c(alpha("orange", 0.25), alpha("sienna", 0.25)), 
        xlab = "Depth (m)", ylab = "Color Morph Proportion",
        space = 0, xaxs="i", yaxs="i")
par(lwd=1)
legend("topright", legend=c("Brown", "Orange"), fill=c("sienna", "orange"), inset = c(-.2, 0), xpd = NA)
par(new = T)
par(mar=c(4.2, 4, 2, 6))
merged$Color2 <- ifelse(merged$Color=="0", 1, 0)
results=glm(Color2~Depth..m., family = "binomial", data = merged)
fitted <- predict(results, newdata = list(Depth..m.=seq(0,11,0.1)), type = "response")
plot(fitted~seq(0,11,0.1), xaxs="i", yaxs="i", xlim=c(0,11), ylim=c(0,1), type="l", lwd = 3)

Type=table(merged$ColDom, merged$Bay.Area)
Type
chisq.test(Type)

B <- subset(XY2, !Reef.ID=="13")
reef.dists <- dist(cbind(B$Longitude, B$Latitude))
dom.dists <- bcdist(B$c.adj)
set.seed(12456)
mantel(dom.dists~reef.dists)

merged$Mixture <- ifelse(!merged$Mix=="C" & !merged$Mix=="D", 1, 0)
results=glm(Mixture~newDepth, family = "binomial", data = merged)
anova(results, test = "Chisq")

#D-only colonies
D <- subset(merged, Mix=="D")

#Mosaic Plot
mosaicplot(total, ylab = "Reef Area", xlab = "Color Morph", main = "")

#Logistic Regression
Symcap$Dom <- as.factor(Symcap$Dom)
results=glm(Color.Morph~Depth..m., family = "binomial", data = Symcap)
anova(results, test = "Chisq")
summary(results)
plot(results)
plot(Symcap$Dom~Symcap$Depth..m.)

#Plot Color Morph and Depth
merged$Color <- ifelse(merged$Color.Morph=="Orange", 1, 0)
plot(merged$Color~merged$newDepth, xlab="Depth (m)", ylab = "Proportion of Color Morph")
results=glm(Color~newDepth, family = "binomial", data = merged)
anova(results, test = "Chisq")
summary(results)
fitted <- predict(results, newdata = list(newDepth=seq(0,12,0.1)), type = "response")
lines(fitted ~ seq(0,12,0.1))

#Plot Dominant Symbiont and Depth
merged$Dominant <- ifelse(merged$Dom=="C", 0, 1)
plot(merged$propC~merged$newDepth, xlab="Depth (m)", ylab = "Proportion of Clade C Symbiont")
results=glm(propC~newDepth, family = "binomial", data = merged)
anova(results, test = "Chisq")
summary(results)
fitted <- predict(results, newdata = list(newDepth=seq(0,12,0.1)), type = "response")
lines(fitted ~ seq(0,12,0.1))

#Logistic Regression/Histogram Plot
merged$Color <- ifelse(merged$Color.Morph=="Orange", 1, 0)
Color <- subset(merged, !is.na(newDepth) & !is.na(Color))
logi.hist.plot(independ = Color$newDepth, depend = Color$Color, type = "hist", boxp = FALSE, ylabel = "", col="gray", ylabel2 = "", xlabel = "Depth (m)")
mtext(side = 4, text = "Frequency", line = 3, cex=1)
mtext(side = 4, text = "Brown                                Orange", line = 2, cex = 0.75)
mtext(side = 2, text = "Probability of Orange Color Morph", line = 3, cex = 1)

Dom1 <- subset(merged, !is.na(newDepth) & !is.na(Dominant))
logi.hist.plot(Dom1$newDepth, Dom1$Dominant, boxp = FALSE, type = "hist", col="gray", xlabel = "Depth (m)", ylabel = "", ylabel2 = "")
mtext(side = 4, text = "Frequency", line = 3, cex=1)
mtext(side = 4, text = "C                                       D", line = 2, cex = 0.75)
mtext(side = 2, text = "Probability of clade C Symbiont", line = 3, cex = 1)

#Export Image
pdf(file="ColDom~Depth")
par(mfrow=c(3,1))

merged$DepthInt <- cut(merged$newDepth, breaks = 0:13)
merged$Dominant <- ifelse(merged$Dom=="C", 0, 1)
merged$Dominant2 <- ifelse(merged$Dom=="C", 1, 0)
results=table(merged$Dominant2, merged$DepthInt)
results
props <- prop.table(results, margin = 2)
par(mar=c(2, 4, 2, 6), lwd = 0.25)
barplot(props[,1:11], col = c(alpha("red", 0.25), alpha("blue", 0.25)), 
        xlab = "", ylab = "",
        space = 0, xaxs="i", yaxs="i", axisnames = FALSE)
par(lwd=1)
legend("topright", legend=c("C", "D"), fill=c(alpha("blue", 0.25), alpha("red", 0.25)), inset = c(0, 0), xpd = NA)
par(new = T)
par(mar=c(2.1, 4, 2, 6))
results=glm(Dominant~newDepth, family = "binomial", data = merged)
fitted <- predict(results, newdata = list(newDepth=seq(0,11,0.1)), type = "response")
plot(fitted~seq(0,11,0.1), xaxs="i", yaxs="i", xlim=c(0,11), ylim=c(0,1), type="l", lwd = 3, xlab="", ylab="Probability of D-Dominance")

merged$Color <- ifelse(merged$Color.Morph=="Orange", 0, 1)
results=table(merged$Color, merged$DepthInt)
results
props <- prop.table(results, margin = 2)
par(mar=c(3, 4, 1, 6), lwd = 0.25)
barplot(props[,1:11], col = c(alpha("orange", 0.25), alpha("sienna", 0.25)), 
        xlab = "", ylab = "Probability of Orange-Dominance",
        space = 0, xaxs="i", yaxs="i", axisnames = FALSE)
par(lwd=1)
legend("topright", legend=c("Brown", "Orange"), fill=c(alpha("sienna", 0.25), alpha("orange", 0.25)), inset = c(0, 0), xpd = NA)
par(new = T)
par(mar=c(3.1, 4, 1, 6))
merged$Color2 <- ifelse(merged$Color=="0", 1, 0)
results=glm(Color2~newDepth, family = "binomial", data = merged)
fitted <- predict(results, newdata = list(newDepth=seq(0,11,0.1)), type = "response")
plot(fitted~seq(0,11,0.1), xaxs="i", yaxs="i", xlim=c(0,11), ylim=c(0,1), type="l", lwd = 3, xlab="", ylab="")

df <- subset(merged, Color.Morph=="Orange")
results=glm(Dominant~newDepth, family = "binomial", data = df)
newdata <- list(newDepth=seq(0,11,0.01))
par(mar=c(4, 4, 0, 6))
fitted <- predict(results, newdata = newdata, type = "response")
plot(fitted ~ seq(0,11,0.01), ylim = c(0,1), type="l", col="orange", lwd=3, xlab="Depth (m)", ylab="Probabilty of D-Dominance", axisnames=FALSE, xaxs = "i", yaxs = "i")
abline(h = 0.5, lty=2)
df <- subset(merged, Color.Morph=="Brown")
results=glm(Dominant~newDepth, family = "binomial", data = df)
newdata <- list(newDepth=seq(0,11,0.01))
fitted <- predict(results, newdata = newdata, type = "response")
lines(fitted~seq(0, 11, 0.01), col="sienna", lwd=3)
legend("topright", legend=c("Brown", "Orange"), fill=c("sienna", "orange"), inset = c(0, 0), xpd = NA)
dev.off()

#RGoogleMaps
KB <- c(21.46087401, -157.809907) 
KBMap <- GetMap(center = KB, zoom = 13, maptype = "satellite", SCALE = 2, GRAYSCALE = FALSE)
Latitude=aggregate(Latitude~Reef.ID, data=Symcap, FUN = mean)
Longitude=aggregate(Longitude~Reef.ID, data = Symcap, FUN=mean)
XY<-merge(Latitude, Longitude, by="Reef.ID", all=T)
propDom=table(Symcap$Dom, Symcap$Reef.ID)
propDom=prop.table(propDom, margin = 2)
propDom <- as.data.frame.matrix(propDom)
props <- data.frame(t(propDom))
props$Reef.ID <- rownames(props)
XY<-merge(XY, props, by="Reef.ID", all=T)
newcoords <- LatLon2XY.centered(KBMap, XY$Latitude, XY$Longitude, zoom=13)
XY$X <- newcoords$newX
XY$Y <- newcoords$newY
par(oma=c(3,3,0,0))
PlotOnStaticMap(KBMap, XY$Latitude, XY$Longitude, col=153, pch=21, bg="#FF7F50", lwd=2)
axis(1, at = LatLon2XY.centered(KBMap, NA, c(-157.85, -157.81, -157.77))$newX, 
     tcl=0.5, line = 0.5, col = "ghostwhite", col.ticks = "black", lwd = 1, outer = TRUE,
     labels = c("157.85°W", "157.81°W", "157.77°W"), padj = -2.5, cex.axis = 0.75)
axis(2, at = LatLon2XY.centered(KBMap, c(21.42, 21.46, 21.50), NA)$newY, 
     tcl=0.5, line = 0.5, col = "ghostwhite", col.ticks = "black", lwd = 1, outer = TRUE,
     labels = c("21.42°N", "21.46°N", "21.50°N"), padj = 0.5, hadj = 0.60, las = 1, 
     cex.axis = 0.75)
par(new=T, mar=c(9,17,0,0))
HI <- readOGR("coast_n83.shp", "coast_n83") 
HI <- spTransform(HI, CRS("+proj=longlat +datum=NAD83")) 
plot(HI, xlim=c(-158.3, -157.6), ylim=c(21.35, 21.6), lwd=0.4, col="gray", bg="white")
rect(-157.9, 21.41, -157.75, 21.53)
box()

rownames(XY) <- XY$Reef.ID
XY <- XY[, -1]
XY <- na.omit(XY)
apply(XY, MARGIN=1, FUN=function(reef) {
  floating.pie(xpos = reef["X"], ypos = reef["Y"], 
               x=c(reef["Orange"], reef["Brown"]), radius = 10, col = c("#7FFFD4"))
})

Latitude=aggregate(Latitude~Reef.ID, data=Symcap, FUN = mean)
Longitude=aggregate(Longitude~Reef.ID, data = Symcap, FUN=mean)
XY<-merge(Latitude, Longitude, by="Reef.ID", all=T)
propMix=table(Symcap$Mix, Symcap$Reef.ID)
propMix=prop.table(propMix, margin = 2)
propMix[which(propMix==0)] <- 0.0000001
propMix <- as.data.frame.matrix(propMix)
props <- data.frame(t(propMix))
props$Reef.ID <- rownames(props)
XY<-merge(XY, props, by="Reef.ID", all=T)
newcoords <- LatLon2XY.centered(KBMap, XY$Latitude, XY$Longitude, zoom=13)
XY$X <- newcoords$newX
XY$Y <- newcoords$newY

# convert to matrix and get rid of reef.id column
rownames(XY) <- XY$Reef.ID
XY <- XY[, -1]
XY <- na.omit(XY)
apply(XY, MARGIN=1, FUN=function(reef) {
  floating.pie(xpos = reef["X"], ypos = reef["Y"], 
               x=c(reef["C"], reef["CD"], reef["DC"], reef["D"]), 
               radius = 7, col = c("#0571b0","#92c5de","#f4a582","#ca0020"))
})

#Plot Color Morph Proportions per Reef
KB <- c(21.46087401, -157.809907) 
KBMap <- GetMap(center = KB, zoom = 13, maptype = "satellite", SCALE = 2, GRAYSCALE = FALSE)
Latitude=aggregate(Latitude~Reef.ID, data=Symcap, FUN = mean)
Longitude=aggregate(Longitude~Reef.ID, data = Symcap, FUN=mean)
XY<-merge(Latitude, Longitude, by="Reef.ID", all=T)
propCol=table(Symcap$Color.Morph, Symcap$Reef.ID)
propCol=prop.table(propCol, margin = 2)
propCol <- as.data.frame.matrix(propCol)
props <- data.frame(t(propCol))
props$Reef.ID <- rownames(props)
XY<-merge(XY, props, by="Reef.ID", all=T)
newcoords <- LatLon2XY.centered(KBMap, XY$Latitude, XY$Longitude, zoom=13)
XY$X <- newcoords$newX
XY$Y <- newcoords$newY

rownames(XY) <- XY$Reef.ID
XY <- XY[, -1]
XY <- na.omit(XY)
apply(XY, MARGIN=1, FUN=function(reef) {
  floating.pie(xpos = reef["X"], ypos = reef["Y"], 
               x=c(reef["Orange"], reef["Brown"]), radius = 7, col = c("orange","sienna"))
})

#Plot Dominant Symbiont vs. Depth and find threshold depth of D to C dominance
threshdepth <- function(reef) {
  df <- subset(merged, Reef.ID==reef)
  plot(df$Dominant~df$newDepth, xlab="Depth (m)", ylab = "Proportion of Dominant Symbiont",
       main=reef)
  abline(h = 0.5, lty=2)
  results=glm(Dominant~newDepth, family = "binomial", data = df)
  pval <- data.frame(coef(summary(results)))$`Pr...z..`[2]
  mtext(side=3, text=round(pval, 3))
  newdata <- list(newDepth=seq(0,12,0.01))
  fitted <- predict(results, newdata = newdata, type = "response")
  lines(fitted ~ seq(0,12,0.01))
  thresh <- ifelse(pval < 0.05,
                   newdata$newDepth[which(diff(sign(fitted - 0.5))!=0)], NA)
  return(thresh)
}

sapply(levels(merged$Reef.ID), FUN=threshdepth)
levels(merged$Reef.ID)

#Depth at which Orange switches from D to C dominance
threshdepth <- function(color, reef) {
  df <- subset(merged, Color.Morph==color & Reef.ID==reef)
  plot(df$Dominant~df$newDepth, xlab="Depth (m)", ylab = "Probability of Clade C Symbiont",
       main=c(color, reef))
  abline(h = 0.5, lty=2)
  results=glm(Dominant~newDepth, family = "binomial", data = df)
  pval <- data.frame(coef(summary(results)))$`Pr...z..`[2]
  mtext(side=3, text=pval)
  newdata <- list(newDepth=seq(0,12,0.01))
  fitted <- predict(results, newdata = newdata, type = "response")
  lines(fitted ~ seq(0,12,0.01))
  thresh <- ifelse(pval < 0.05,
                   newdata$newDepth[which(diff(sign(fitted - 0.5))!=0)], NA)
  return(c(thresh, pval))
}

threshdepth(color = "Orange", reef = "42")
sapply(levels(merged$Reef.ID), FUN=function(x) threshdepth(color="Orange", reef=x))

threshdepth <- function(color) {
  df <- subset(merged, Color.Morph==color)
  plot(df$Dominant~df$newDepth, xlab="Depth (m)", ylab = "Probability of Clade C Symbiont",
       main=color)
  abline(h = 0.5, lty=2)
  results=glm(Dominant~newDepth, family = "binomial", data = df)
  pval <- data.frame(coef(summary(results)))$`Pr...z..`[2]
  mtext(side=3, text=pval)
  newdata <- list(newDepth=seq(0,12,0.01))
  fitted <- predict(results, newdata = newdata, type = "response")
  lines(fitted ~ seq(0,12,0.01))
  thresh <- ifelse(pval < 0.05,
                   newdata$newDepth[which(diff(sign(fitted - 0.5))!=0)], NA)
  return(thresh)
}

#Plot
df <- subset(merged, Color.Morph=="Orange")
results=glm(Dominant~newDepth, family = "binomial", data = df)
newdata <- list(newDepth=seq(0,12,0.01))
fitted <- predict(results, newdata = newdata, type = "response")
plot(fitted ~ seq(0,12,0.01), ylim = c(0,1), type="l", col="orange", lwd=3)
abline(h = 0.5, lty=2)
df <- subset(merged, Color.Morph=="Brown")
results=glm(Dominant~newDepth, family = "binomial", data = df)
newdata <- list(newDepth=seq(0,12,0.01))
fitted <- predict(results, newdata = newdata, type = "response")
lines(fitted~seq(0,12,0.01), col="sienna", lwd=3)

# 3 Variables 
merged$Reef.Area <- ifelse(merged$Reef.Area!="Top", yes = "Slope", no = "Top")
table(merged$Dom, merged$Color.Morph, merged$Bay.Area)

#Patch vs. Fringe
Patch <- subset(merged, Reef.Type=="Patch")
Fringe <- subset(merged, Reef.Type=="Fringe")

merged$Color <- ifelse(merged$Color.Morph=="Orange", 1, 0)
results=glm(Color~newDepth, family = "binomial", data = Fringe)
anova(results, test = "Chisq")
Color <- subset(Fringe, !is.na(newDepth) & !is.na(Color))
logi.hist.plot(independ = Color$newDepth, depend = Color$Color, type = "hist", boxp = FALSE, ylabel = "", col="gray", ylabel2 = "", xlabel = "Depth (m)")
mtext(side = 4, text = "Frequency", line = 3, cex=1)
mtext(side = 4, text = "Brown                                Orange", line = 2, cex = 0.75)
mtext(side = 2, text = "Probability of Orange Color Morph", line = 3, cex = 1)

df <- subset(merged, Reef.Type=="Patch")
results=glm(Color2~newDepth, family = "binomial", data = df)
newdata <- list(newDepth=seq(0,12,0.01))
fitted <- predict(results, newdata = newdata, type = "response")
plot(fitted ~ seq(0,12,0.01), ylim = c(0,1), type="l", col="black", lwd=3)
abline(h = 0.5, lty=2)
df <- subset(merged, Reef.Type=="Fringe")
results=glm(Color2~newDepth, family = "binomial", data = df)
newdata <- list(newDepth=seq(0,12,0.01))
fitted <- predict(results, newdata = newdata, type = "response")
lines(fitted~seq(0,12,0.01), col="red", lwd=3)

#2-Way ANOVA - Interactive Effects 
merged$Dominant <- ifelse(merged$Dom=="C", 1, 0)
model3=aov(Dominant~Reef.Area*Color, data = merged)
Anova(model3, type=2)

model1=lm(Dominant~Reef.Type*newDepth, data = merged)
anova(model1)

merged$Dominant <- ifelse(merged$Dom=="C", 1, 0)
model4=aov(Dominant~Depth..m.*Reef.ID, data = merged)
Anova(model4, type = 2)

model3=aov(Dominant2~Depth..m.*Reef.Type, data = merged)
Anova(model3, type = 2)

#Dominant Symbiont per Depth and Reef ID
merged$Dominant2 <- ifelse(merged$Dom=="C", 0, 1)

domreef <- function(id) {
  df <- subset(merged, Reef.ID==id)
  results=glm(Dominant2~Depth..m., family = "binomial", data = df)
  newdata <- list(Depth..m.=seq(0,12,0.01))
  par(mar=c(4, 4, 2, 6))
  fitted <- predict(results, newdata = newdata, type = "response")
  plot(fitted ~ seq(0,12,0.01), ylim = c(0,1), type="l", 
       col="dodgerblue3", lwd=3, xlab="", ylab="")
  mtext(side = 3, text = id)
  mtext(side = 2, text = "Probability of Clade D Symbiont", line = 3, cex = 1)
  mtext(side = 1, text = "Depth (m)", line = 3, cex = 1)
}

sapply(levels(merged$Reef.ID), FUN=domreef)

#MANOVA
merged$Dominant <- ifelse(merged$Dom=="C", 0, 1)
merged$Color <- ifelse(merged$Color.Morph=="Orange", 1, 0)
DomCol <- cbind(merged$Dominant, merged$Color)
fit <- manova(DomCol~merged$Depth..m.)
manova(fit)
summary(fit)

merged$newDepth <- as.factor(merged$newDepth)
fit2 <- manova(DomCol~merged$newDepth)
manova(fit2)
summary(fit2)

# 3-panel figure
attach(mtcars)
par(mfrow=c(3,1))
merged$DepthInt <- cut(merged$Depth..m., breaks = 0:13)
merged$Dominant <- ifelse(merged$Dom=="C", 0, 1)
merged$Dominant2 <- ifelse(merged$Dom=="C", 1, 0)
results=table(merged$Dominant2, merged$DepthInt)
results
props <- prop.table(results, margin = 2)
par(mar=c(4, 4, 2, 6), lwd = 0.25)
barplot(props[,1:11], col = c(alpha("red", 0.25), alpha("blue", 0.25)), 
        xlab = "", ylab = "",
        space = 0, xaxs="i", yaxs="i", axisnames = FALSE)
par(lwd=1)
legend("topright", legend=c("C", "D"), fill=c(alpha("blue", 0.25), alpha("red", 0.25)), inset = c(-.1, 0), xpd = NA)
par(new = T)
par(mar=c(4.2, 4, 2, 6))
results=glm(Dominant~Depth..m., family = "binomial", data = merged)
fitted <- predict(results, newdata = list(Depth..m.=seq(0,11,0.1)), type = "response")
plot(fitted~seq(0,11,0.1), xaxs="i", yaxs="i", xlim=c(0,11), ylim=c(0,1), type="l", lwd = 3, xlab="Depth (m)", ylab="Dominant Symbiont Proportion")
merged$Color <- ifelse(merged$Color.Morph=="Orange", 0, 1)
results=table(merged$Color, merged$DepthInt)
results
props <- prop.table(results, margin = 2)
par(mar=c(4, 4, 2, 6), lwd = 0.25)
barplot(props[,1:11], col = c(alpha("orange", 0.25), alpha("sienna", 0.25)), 
        xlab = "", ylab = "",
        space = 0, xaxs="i", yaxs="i", axisnames = FALSE)
par(lwd=1)
legend("topright", legend=c("Brown", "Orange"), fill=c(alpha("sienna", 0.25), alpha("orange", 0.25)), inset = c(-0.16, 0), xpd = NA)
par(new = T)
par(mar=c(4.2, 4, 2, 6))
merged$Color2 <- ifelse(merged$Color=="0", 1, 0)
results=glm(Color2~Depth..m., family = "binomial", data = merged)
fitted <- predict(results, newdata = list(Depth..m.=seq(0,11,0.1)), type = "response")
plot(fitted~seq(0,11,0.1), xaxs="i", yaxs="i", xlim=c(0,11), ylim=c(0,1), type="l", lwd = 3, xlab="Depth (m)", ylab="Color Morph Proportion")
df <- subset(merged, Color.Morph=="Orange")
results=glm(Dominant~newDepth, family = "binomial", data = df)
newdata <- list(newDepth=seq(0,12,0.01))
par(mar=c(4, 4, 2, 6))
fitted <- predict(results, newdata = newdata, type = "response")
plot(fitted ~ seq(0,12,0.01), ylim = c(0,1), type="l", col="orange", lwd=3, xlab="", ylab="", axisnames=FALSE)
abline(h = 0.5, lty=2)
df <- subset(merged, Color.Morph=="Brown")
results=glm(Dominant~newDepth, family = "binomial", data = df)
newdata <- list(newDepth=seq(0,12,0.01))
fitted <- predict(results, newdata = newdata, type = "response")
lines(fitted~seq(0,12,0.01), col="sienna", lwd=3)
mtext(side = 1, text = "Depth (m)", line = 3, cex = 1)
mtext(side = 2, text = "Probability of Clade C Symbiont", line = 3, cex = 1)
legend("topright", legend=c("Brown", "Orange"), fill=c("sienna", "orange"), inset = c(-.16, 0), xpd = NA)

# Proportion of D 
propD <- merged$propD[which(merged$propD > 0 & merged$propD < 1.1)]
hist(propD, xlab = "Proportion of Clade D", ylab = "Number of Samples", main = "", col = "gray75")

# C Mixture by Depth
C <- subset(merged, Dom=="C")
C$Mixture <- ifelse(C$Mix=="C", 1, 0)
C$Mixture2 <- ifelse(C$Mix=="C", 0, 1)
results=glm(Mixture~newDepth, family = "binomial", data = C)
anova(results, test = "Chisq")
results=table(C$Mixture2, C$DepthInt)
results
props <- prop.table(results, margin = 2)
par(mar=c(4, 4, 2, 6), lwd = 0.25)
barplot(props[,1:11], col = c(alpha("red", 0.25), alpha("blue", 0.25)), 
        xlab = "", ylab = "",
        space = 0, xaxs="i", yaxs="i", axisnames = FALSE)
par(lwd=1)
legend("topright", legend=c("CD", "C"), fill=c(alpha("blue", 0.25), alpha("red", 0.25)), inset = c(-.23, 0), xpd = NA)
par(new = T)
par(mar=c(4.2, 4, 2, 6))
results=glm(Mixture~Depth..m., family = "binomial", data = C)
fitted <- predict(results, newdata = list(Depth..m.=seq(0,11,0.1)), type = "response")
plot(fitted~seq(0,11,0.1), xaxs="i", yaxs="i", xlim=c(0,11), ylim=c(0,1), type="l", lwd = 3, xlab="Depth (m)", ylab="Mixture Proportion")

results=table(merged$Dominant, merged$Color.Morph)
chisq.test(results)
prop.table(results, margin = 2)
par(mar=c(4, 4, 2, 6))
barplot(prop.table(results, margin = 2), col = c("gray10", "gray100"), xlab = "Color Morph", ylab = "Symbiont Proportion")
legend("topright", legend=c("C", "D"), fill=c("gray10", "gray100"), inset = c(-.2, 0), xpd = NA)

# Dominant Symbiont per Color Morph by Depth and Reef Type
df <- subset(merged, Reef.Type=="Patch")
dfo <- subset(df, Color.Morph=="Orange")
results=glm(Dominant~Depth..m., family = "binomial", data = dfo)
newdata <- list(Depth..m.=seq(0,12,0.01))
par(mar=c(4, 4, 2, 6))
fitted <- predict(results, newdata = newdata, type = "response")
plot(fitted ~ seq(0,12,0.01), ylim = c(0,1), type="l", col="dodgerblue3", lwd=3, xlab="", ylab="", axisnames=FALSE)
abline(h = 0.5, lty=2)

df <- subset(merged, Reef.Type=="Patch")
dfb <- subset(df, Color.Morph=="Brown")
results=glm(Dominant~Depth..m., family = "binomial", data = dfb)
newdata <- list(Depth..m.=seq(0,12,0.01))
par(new = T, mar=c(4, 4, 2, 6))
fitted <- predict(results, newdata = newdata, type = "response")
plot(fitted ~ seq(0,12,0.01), ylim = c(0,1), type="l", col="black", lwd=3, xlab="", ylab="", axisnames=FALSE)
abline(h = 0.5, lty=2)

df <- subset(merged, Reef.Type=="Fringe")
dfo <- subset(df, Color.Morph=="Orange")
results=glm(Dominant~Depth..m., family = "binomial", data = dfo)
newdata <- list(Depth..m.=seq(0,12,0.01))
par(new = T, mar=c(4, 4, 2, 6))
fitted <- predict(results, newdata = newdata, type = "response")
plot(fitted ~ seq(0,12,0.01), ylim = c(0,1), type="l", col="red", lwd=3, xlab="", ylab="", axisnames=FALSE)
abline(h = 0.5, lty=2)

df <- subset(merged, Reef.Type=="Fringe")
dfb <- subset(df, Color.Morph=="Brown")
results=glm(Dominant~Depth..m., family = "binomial", data = dfb)
newdata <- list(Depth..m.=seq(0,12,0.01))
par(new = T, mar=c(4, 4, 2, 6))
fitted <- predict(results, newdata = newdata, type = "response")
plot(fitted ~ seq(0,12,0.01), ylim = c(0,1), type="l", col="orange", lwd=3, xlab="", ylab="", axisnames=FALSE)
abline(h = 0.5, lty=2)
legend("topright", legend=c("PO", "PB", "FO", "FB"), fill=c("dodgerblue3", "black", "red", "orange"), inset = c(-.2, 0), xpd = NA)

#Interaction of Color and Symbiont at Depth
merged$DomCol <- interaction(merged$Dom, merged$Color.Morph)
merged$DomCol <- factor(merged$DomCol, levels=rev(levels(merged$DomCol)))
results=table(merged$DomCol, merged$DepthInt)
results
props <- prop.table(results, margin = 2)
par(mar=c(4, 4, 2, 6), lwd = 0.25)
barplot(props[,1:11], col = c(alpha("orange", 0.75), alpha("orange", 0.25), alpha("sienna", 0.75), alpha("sienna", 0.25)), 
        xlab = "", ylab = "",
        space = 0, xaxs="i", yaxs="i", axisnames = FALSE)
par(lwd=1)
legend("topright", legend=c("OD", "OC", "BD", "BC"), fill=c(alpha("orange", 0.75), alpha("orange", 0.25), alpha("sienna", 0.75), alpha("sienna", 0.25)), inset = c(-.23, 0), xpd = NA)

merged$DomCol2 <- interaction(merged$Color.Morph, merged$Dom)
merged$DomCol2 <- factor(merged$DomCol2, levels=rev(levels(merged$DomCol2)))
results=table(merged$DomCol2, merged$DepthInt)
results
props <- prop.table(results, margin = 2)
par(mar=c(4, 4, 2, 6), lwd = 0.25)
barplot(props[,1:11], col = c(alpha("orange", 0.75), alpha("sienna", 0.75), alpha("orange", 0.25), alpha("sienna", 0.25)), 
        xlab = "", ylab = "",
        space = 0, xaxs="i", yaxs="i", axisnames = FALSE)
par(lwd=1)
legend("topright", legend=c("OD", "BD", "OC", "BC"), fill=c(alpha("orange", 0.75), alpha("sienna", 0.75), alpha("orange", 0.25), alpha("sienna", 0.25)), inset = c(-.23, 0), xpd = NA)

merged$Mixture <- ifelse(!merged$Mix=="C" & !merged$Mix=="D", 1, 0)
model4=aov(Mixture~newDepth*Bay.Area, data = merged)
Anova(model4, type = 2)

results=table(Symcap$Dom, Symcap$Reef.ID)
chisq.test(results)
prop.table(results, margin = 2)

KB <- c(21.46087401, -157.809907) 
KBMap <- GetMap(center = KB, zoom = 13, maptype = "satellite", SCALE = 2, GRAYSCALE = FALSE)
Latitude=aggregate(Latitude~Reef.ID, data=Symcap, FUN = mean)
Longitude=aggregate(Longitude~Reef.ID, data = Symcap, FUN=mean)
XY<-merge(Latitude, Longitude, by="Reef.ID", all=T)
propDom=table(Symcap$Dom, Symcap$Reef.ID)
propDom=prop.table(propDom, margin = 2)
propDom <- as.data.frame.matrix(propDom)
props <- data.frame(t(propDom))
props$Reef.ID <- rownames(props)
XY<-merge(XY, props, by="Reef.ID", all=T)
newcoords <- LatLon2XY.centered(KBMap, XY$Latitude, XY$Longitude, zoom=13)
XY$X <- newcoords$newX
XY$Y <- newcoords$newY
PlotOnStaticMap(KBMap, df$Latitude, df$Longitude)

rownames(XY) <- XY$Reef.ID
XY <- XY[, -1]
XY <- na.omit(XY)
apply(XY, MARGIN=1, FUN=function(reef) {
  floating.pie(xpos = reef["X"], ypos = reef["Y"], 
               x=c(reef["C"], reef["D"]), radius = 7, col = c("blue", "red"))
})

#Prop D Histogram
merged$propDHist <- cut(merged$propD, breaks = 10)
DCol=table(merged$Color.Morph, merged$propDHist)
barplot(DCol)

propD <- merged$propD[which(merged$propD > 0 & merged$propD < 1)]
hist(propD, xlab = "Proportion of Clade D", ylab = "Number of Samples", 
     main = "", col = "ghostwhite")

par(mfrow=c(2,1))
propDHist <- subset(merged, propD > 0 & propD < 1)
propDHist$propD <- cut(propDHist$propD, breaks = 10)
DCol=table(propDHist$Color.Morph, propDHist$propD)
par(mar=c(2, 4, 2, 6))
barplot(DCol, xlab = "Clade D Proportion", ylab = "Number of Samples", 
        main = "", col = c(alpha("sienna", 0.55), alpha("orange", 0.55)), 
        axisnames = FALSE, space = 0)

img <- readPNG("brown_orange.png")
img2 <- pixmapRGB(img)
plot(img2)

#Mantel Test for Symbiont and Color
manteltable = table(merged$Dom, merged$Color.Morph, merged$Reef.ID)
nc <- aggregate(interaction(merged$Color.Morph, merged$Dom), by=list(merged$Reef.ID), FUN=table)
nc <- data.frame(Reef.ID=as.character(nc[,1]), prop.table(nc[,2], margin=1))
nc
XY3 <- merge(XY2, nc, by="Reef.ID", all=T)

#Mantel Test for Dominant Symbiont
Latitude=aggregate(Latitude~Reef.ID, data=Symcap, FUN = mean)
Longitude=aggregate(Longitude~Reef.ID, data = Symcap, FUN=mean)
XY<-merge(Latitude, Longitude, by="Reef.ID", all=T)
propDom=table(Symcap$Dom, Symcap$Reef.ID)
propDom=prop.table(propDom, margin = 2)
propDom <- as.data.frame.matrix(propDom)
props <- data.frame(t(propDom))
props$Reef.ID <- rownames(props)
XY<-merge(XY, props, by="Reef.ID", all=T)

reef.dists <- dist(cbind(XY3$Longitude, XY3$Latitude))
dom.dists <- dist(cbind(XY3$Brown.C, XY3$Brown.D, XY3$Orange.C, XY3$Orange.D),
                  method="manhattan")
set.seed(12456)
mantel(dom.dists~reef.dists)

rownames(XY3) <- XY3$Reef.ID
XY3 <- XY3[, -1]
XY3 <- na.omit(XY3)
apply(XY3, MARGIN=1, FUN=function(reef) {
  floating.pie(xpos = reef["X"], ypos = reef["Y"], 
               x=c(reef["Brown.C"], reef["Brown.D"], reef["Orance.C"], reef["Orange.D"]), radius = 7, 
               col = c("blue", "purple", "green", "red"))
})

domsk <- pamk(data = dom.dists)
domsk

km <- kmeans(dom.dists, centers = 6)
km

df <- cbind(XY, km$cluster)
km$cluster

PlotOnStaticMap(KBMap, XY2$Latitude, XY2$Longitude, col="red", 
                pch=as.character(km$cluster), lwd = 2)

KB <- c(21.46087401, -157.809907) 
KBMap <- GetMap(center = KB, zoom = 13, maptype = "satellite", SCALE = 2, GRAYSCALE = FALSE)
Latitude=aggregate(Latitude~Reef.ID, data=Symcap, FUN = mean)
Longitude=aggregate(Longitude~Reef.ID, data = Symcap, FUN=mean)
XY<-merge(Latitude, Longitude, by="Reef.ID", all=T)
newcoords <- LatLon2XY.centered(KBMap, XY$Latitude, XY$Longitude, zoom=13)
XY$X <- newcoords$newX
XY$Y <- newcoords$newY
XY <- subset(XY, Reef.ID!="37")

# Prevalence of C vs. D dominance, adjusted for depth.
merged$Dominant <- ifelse(merged$Dom=="C", 0, 1)
dat <- aggregate(data.frame(prop=merged$Dominant),
                 by=list(Reef.ID=merged$Reef.ID), FUN=mean, na.rm=T) #prop D/reef
mod <- glm(Dominant ~ newDepth + Reef.ID, data=merged)
mod2 <- glm(Dominant ~ newDepth * Reef.ID, data=merged)
library(lsmeans)
lsm <- lsmeans(mod, specs="Reef.ID")
lsm2 <- lsmeans(mod2, specs="Reef.ID")
res <- merge(dat, data.frame(summary(lsm))[,c(1:2)], by="Reef.ID")
res <- merge(res, data.frame(summary(lsm2))[,c(1:2)], by="Reef.ID")
colnames(res) <- c("Reef.ID", "raw.c", "c.adj", "c.adj.int")

XY2 <- merge(XY, res, by = "Reef.ID")
XY2$d.adj <- 1-XY2$c.adj
XY2$d.adj.int <- 1-XY2$c.adj.int
XY2$raw.d <- 1-XY2$raw.c

KB <- c(21.46087401, -157.809907) 
KBMap <- GetMap(center = KB, zoom = 13, maptype = "satellite", SCALE = 2, GRAYSCALE = FALSE)
PlotOnStaticMap(KBMap, XY$Latitude, XY$Longitude)

rownames(XY2) <- XY2$Reef.ID
XY2 <- XY2[, -1]
XY2 <- na.omit(XY2)
apply(XY2, MARGIN=1, FUN=function(reef) {
  floating.pie(xpos = reef["X"], ypos = reef["Y"], 
               x=c(reef["c.adj"], reef["d.adj"]), radius = 7, 
               col = c("blue", "red"))
})

# Symbiont and Color Adjusted for Depth
merged$Dominant <- ifelse(merged$Dom=="C", 0, 1)
dat <- aggregate(data.frame(prop=merged$Dominant), by=list(Reef.ID=merged$Reef.ID), 
                 FUN=mean, na.rm=T) #prop D/reef
mod <- glm(Dominant ~ newDepth + Reef.ID, data=merged)
mod2 <- glm(Dominant ~ newDepth * Reef.ID, data=merged)
library(lsmeans)
lsm <- lsmeans(mod, specs="Reef.ID")
lsm2 <- lsmeans(mod2, specs="Reef.ID")
res <- merge(dat, data.frame(summary(lsm))[,c(1:2)], by="Reef.ID")
res <- merge(res, data.frame(summary(lsm2))[,c(1:2)], by="Reef.ID")
colnames(res) <- c("Reef.ID", "raw.c", "c.adj", "c.adj.int")

XY2 <- merge(XY, res, by = "Reef.ID")
XY2$d.adj <- 1-XY2$c.adj
XY2$d.adj.int <- 1-XY2$c.adj.int
XY2$raw.d <- 1-XY2$raw.c

manteltable = table(merged$Dom, merged$Color.Morph, merged$Reef.ID)
nc <- aggregate(interaction(merged$Color.Morph, merged$Dom), 
                by=list(merged$Reef.ID), FUN=table)
nc <- data.frame(Reef.ID=as.character(nc[,1]), prop.table(nc[,2], margin=1))
XY3 <- merge(XY2, nc, by="Reef.ID", all=T)

merged$ColDom <- interaction(merged$Color.Morph, merged$Dom)
model <- multinom(ColDom~newDepth+Reef.ID, merged)
means <- lsmeans(model, specs = "Reef.ID")
means
pp <- fitted(model)

newdat <- data.frame(Reef.ID = levels(merged$Reef.ID),
                   newDepth = mean(merged$newDepth, na.rm=T))
pred <- predict(model, newdata = newdat, "probs")
new <- data.frame(Reef.ID=as.character(newdat[,1]), pred)
XY4 <- merge(XY3, new, by="Reef.ID", all=T)

reef.dists <- dist(cbind(XY4$Longitude, XY4$Latitude))
dom.dists <- bcdist(cbind(XY4$Brown.C.y, XY4$Orange.C.y, XY4$Brown.D.y, XY4$Orange.D.y))
set.seed(12456)
mantel(dom.dists~reef.dists)

KB <- c(21.46087401, -157.809907) 
KBMap <- GetMap(center = KB, zoom = 13, maptype = "satellite", SCALE = 2, GRAYSCALE = FALSE)
PlotOnStaticMap(KBMap, XY4$Latitude, XY4$Longitude)

XY4 <- merge(XY3, new, by="Reef.ID", all=T)
XY4[XY4 == 0.000] <- 0.0000000001
rownames(XY4) <- XY4$Reef.ID
XY4 <- XY4[, -1]
XY4 <- na.omit(XY4)
apply(XY4, MARGIN=1, FUN=function(reef) {
  floating.pie(xpos = reef["X"], ypos = reef["Y"], 
               x=c(reef["Brown.C.x"], reef["Orange.C.x"], reef["Brown.D.x"], 
                   reef["Orange.D.x"]), radius = 7, 
               col = c("blue", "purple", "green", "red"))
})
legend("topright", legend=c("BC", "OC", "BD", "OD"), fill = c("blue", "purple", "green", "red"))


KB <- c(21.46087401, -157.809907) 
KBMap <- GetMap(center = KB, zoom = 13, maptype = "satellite", SCALE = 2, GRAYSCALE = FALSE)
PlotOnStaticMap(KBMap, XY4$Latitude, XY4$Longitude)

XY4 <- merge(XY3, new, by="Reef.ID", all=T)
XY4[XY4 == 0.000] <- 0.0000000001
rownames(XY4) <- XY4$Reef.ID
XY4 <- XY4[, -1]
XY4 <- na.omit(XY4)
apply(XY4, MARGIN=1, FUN=function(reef) {
  floating.pie(xpos = reef["X"], ypos = reef["Y"], 
               x=c(reef["Brown.C.y"], reef["Orange.C.y"], reef["Brown.D.y"], 
                   reef["Orange.D.y"]), radius = 7, 
               col = c("blue", "purple", "green", "red"))
})
legend("topright", legend=c("BC", "OC", "BD", "OD"), fill = c("blue", "purple", "green", "red"))


# plot color morph and symbiont bars
merged$DomCol <- interaction(merged$Dom, merged$Color.Morph)
merged$DomCol <- factor(merged$DomCol, levels=rev(levels(merged$DomCol)))
results=table(merged$DomCol, merged$DepthInt)
results
props <- prop.table(results, margin = 2)
par(mar=c(4, 4, 2, 6), lwd = 0.25)
barplot(props[,1:11], col = c(alpha("orange", 0.75), alpha("orange", 0.25), alpha("sienna", 0.75), alpha("sienna", 0.25)), 
        xlab = "", ylab = "",
        space = 0, xaxs="i", yaxs="i", axisnames = FALSE)
par(lwd=1)
legend("topright", legend=c("OD", "OC", "BD", "BC"), fill=c(alpha("orange", 0.75), alpha("orange", 0.25), alpha("sienna", 0.75), alpha("sienna", 0.25)), inset = c(-.23, 0), xpd = NA)
par(new=T, mar=c(4.2, 4, 2, 6))
box()

merged$DomCol2 <- interaction(merged$Color.Morph, merged$Dom)
merged$DomCol2 <- factor(merged$DomCol2, levels=rev(levels(merged$DomCol2)))
results=table(merged$DomCol2, merged$DepthInt)
results
props <- prop.table(results, margin = 2)
par(mar=c(4, 4, 2, 6), lwd = 0.25)
barplot(props[,1:11], col = c(alpha("orange", 0.75), alpha("sienna", 0.75), alpha("orange", 0.25), alpha("sienna", 0.25)), 
        xlab = "", ylab = "",
        space = 0, xaxs="i", yaxs="i", axisnames = FALSE)
par(lwd=1)
legend("topright", legend=c("OD", "BD", "OC", "BC"), fill=c(alpha("orange", 0.75), alpha("sienna", 0.75), alpha("orange", 0.25), alpha("sienna", 0.25)), inset = c(-.23, 0), xpd = NA)
par(new=T, mar=c(4.2, 4, 2, 6))
box()

merged$DepthInt <- cut(merged$newDepth, breaks = 0:13)
results=table(merged$Mix, merged$DepthInt)
results
props <- prop.table(results, margin = 2)
par(mar=c(4, 4, 2, 6), lwd = 0.25)
barplot(props[,1:11], col = c("red", "orange", "yellow", "green"), 
        xlab = "", ylab = "",
        space = 0, xaxs="i", yaxs="i", axisnames = FALSE)
par(lwd=1)
legend("topright", legend=c("C", "CD", "DC", "D"), 
       fill=c("red", "orange", "yellow", "green"), inset = c(-.2, 0), xpd = NA)
par(new = T)
par(mar=c(4.2, 4, 2, 6))
results=glm(Mix~newDepth, family = "binomial", data = merged)
fitted <- predict(results, newdata = list(newDepth=seq(0,11,0.1)), type = "response")
plot(fitted~seq(0,11,0.1), xaxs="i", yaxs="i", xlim=c(0,11), ylim=c(0,1), type="l", lwd = 3, xlab="Depth (m)", ylab="Dominant Symbiont Proportion")

merged$Color <- ifelse(merged$Color.Morph=="Brown", 0, 1)
dat <- aggregate(data.frame(prop=merged$Color), by=list(Reef.ID=merged$Reef.ID), FUN=mean, na.rm=T) #prop D/reef
mod <- glm(Color ~ newDepth + Reef.ID, data=merged)
mod2 <- glm(Dominant ~ newDepth * Reef.ID, data=merged)
library(lsmeans)
lsm <- lsmeans(mod, specs="Reef.ID")
lsm2 <- lsmeans(mod2, specs="Reef.ID")
res <- merge(dat, data.frame(summary(lsm))[,c(1:2)], by="Reef.ID")
res <- merge(res, data.frame(summary(lsm2))[,c(1:2)], by="Reef.ID")
colnames(res) <- c("Reef.ID", "raw", "Brown")

XY2 <- merge(XY, res, by = "Reef.ID")
XY2$Orange <- 1-XY2$Brown.y

reef.dists <- dist(cbind(XY2$Longitude, XY2$Latitude))
col.dists <- dist(XY2$Orange)
set.seed(12456)
mantel(col.dists~reef.dists)

model3=aov(Dominant~Bay.Area, data = c)
Anova(model3, type = 2)

merged$Bay.Area[which(merged$Reef.ID=="26")] <- "Central"
merged$Bay.Area[which(merged$Reef.ID=="F4-34")] <- "Central"

results=table(merged$Dom, merged$Bay.Area)
chisq.test(results)

c <- subset(merged, !Reef.ID=="Deep" & !Reef.ID=="13")
results=table(c$Dom, c$Bay.Area)
chisq.test(results)

b$Dominant <- ifelse(b$Dom=="C", 0, 1)
dat <- aggregate(data.frame(prop=b$Dominant), by=list(Bay.Area=b$Bay.Area), FUN=mean, na.rm=T) #prop D/reef
mod <- glm(Dominant ~ newDepth + Bay.Area, data=b)
library(lsmeans)
lsm <- lsmeans(mod, specs="Bay.Area")
res <- merge(dat, data.frame(summary(lsm))[,c(1:2)], by="Bay.Area")
colnames(res) <- c("Bay.Area", "raw.c", "c.adj")

results=table(res$c.adj, res$Bay.Area)
chisq.test(results)

b$ColDom <- interaction(b$Color.Morph, b$Dom)
model <- multinom(ColDom~newDepth+Bay.Area, b)
means <- lsmeans(model, specs = "Bay.Area")

Type=table(b$ColDom, b$Bay.Area)
prop.table(Type, margin = 2)
chisq.test(Type)

merged$Bay.Area[which(merged$Reef.ID=="26")] <- "Central"
merged$Bay.Area[which(merged$Reef.ID=="18")] <- "Southern"
merged$Bay.Area[which(merged$Reef.ID=="F7-18")] <- "Southern"
model3=aov(Color~newDepth*Bay.Area, data = merged)
Anova(model3, type = 2)

df <- subset(merged, Bay.Area=="Northern")
results=glm(Dominant~newDepth, family = "binomial", data = df)
newdata <- list(newDepth=seq(0,11,0.01))
par(mar=c(4, 4, 2, 6))
fitted <- predict(results, newdata = newdata, type = "response")
plot(fitted ~ seq(0,11,0.01), ylim = c(0,1), type="l", col="red", lwd=3, xlab="", ylab="", axisnames=FALSE)
abline(h = 0.5, lty=2)
df <- subset(merged, Bay.Area=="Central")
results=glm(Dominant~newDepth, family = "binomial", data = df)
newdata <- list(newDepth=seq(0,11,0.01))
fitted <- predict(results, newdata = newdata, type = "response")
lines(fitted~seq(0,11,0.01), col="gold", lwd=3)
df <- subset(merged, Bay.Area=="Southern")
results=glm(Dominant~newDepth, family = "binomial", data = df)
newdata <- list(newDepth=seq(0,11,0.01))
fitted <- predict(results, newdata = newdata, type = "response")
lines(fitted~seq(0,11,0.01), col="green", lwd=3)

z <- subset(merged, !Reef.ID=="Deep")
results=table(z$DomCol, z$Bay.Area)
prop.table(results, margin = 2)
chisq.test(results)
