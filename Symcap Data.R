library(data.table)
library(devtools)
library(plyr)
library(reshape2)
library(popbio)
library(RgoogleMaps)
library(plotrix)
library(zoo)

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

#Adjust depth by sea level
JuneTide=read.csv("Station_1612480_tide_ht_20160601-20160630.csv")
JulyTide=read.csv("Station_1612480_tide_ht_20160701-20160731.csv")
Tide<-rbind(JuneTide, JulyTide)
Tide$Time <- as.POSIXct(Tide$TimeUTC, format="%Y-%m-%d %H:%M:%S", tz="UTC")
attributes(Tide$Time)$tzone <- "Pacific/Honolulu"
plot(Tide$Time, Tide$TideHT, type = "l", xlab = "Date", ylab = "Tide Height, meters")
with(Tide[Tide$Time > "2016-06-07 00:00:00" & Tide$Time < "2016-06-09 00:00:00", ], 
     plot(Time, TideHT, type="l"))
Symcap$Time2 <- as.POSIXct(paste(as.character(Symcap$Date), as.character(Symcap$Time)),
                                format="%m/%d/%y %H:%M", tz="Pacific/Honolulu")
Symcap$Time=Symcap$Time2

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
Symcap$Reef.Area <- ifelse(Symcap$Reef.Area!="Top", yes = "Slope", no = "Top")
results=table(Symcap$Mix, Symcap$Reef.Area)
results
chisq.test(results)
prop.table(results, margin = 2)
par(mar=c(3, 4, 2, 6))
barplot(prop.table(results, margin = 2), col = c("gray10", "gray40", "gray70", "gray100"), xlab = "Reef Area", ylab = "Symbiont Community Composition")
legend("topright", legend=c("C", "CD", "DC", "D"), fill=c("gray10", "gray40", "gray70", "gray100"), inset = c(-.2, 0), xpd = NA)

Type=table(Symcap$Mix, Symcap$Reef.Type)
Type
chisq.test(Type)
prop.table(Type, margin = 2)

#Mosaic Plot
mosaicplot(total, ylab = "Reef Area", xlab = "Color Morph", main = "")

#Logistic Regression/ANOVA
Symcap$Dom <- as.factor(Symcap$Dom)
results=glm(Color.Morph~Depth..m., family = "binomial", data = Symcap)
anova(results, test = "Chisq")
summary(results)
plot(results)
plot(Symcap$Dom~Symcap$Depth..m.)

#Plot Color Morph and Depth
Symcap$Color <- ifelse(Symcap$Color.Morph=="Orange", 1, 0)
plot(Symcap$Color~Symcap$Depth..m., xlab="Depth (m)", ylab = "Proportion of Color Morph")
results=glm(Color~Depth..m., family = "binomial", data = Symcap)
anova(results, test = "Chisq")
summary(results)
fitted <- predict(results, newdata = list(Depth..m.=seq(0,12,0.1)), type = "response")
lines(fitted ~ seq(0,12,0.1))

#Plot Dominant Symbiont and Depth
Symcap$Dominant <- ifelse(Symcap$Dom=="C", 0, 1)
plot(Symcap$propC~Symcap$Depth..m., xlab="Depth (m)", ylab = "Proportion of Dominant Symbiont")
results=glm(propC~Depth..m., family = "binomial", data = Symcap)
anova(results, test = "Chisq")
summary(results)
fitted <- predict(results, newdata = list(Depth..m.=seq(0,12,0.1)), type = "response")
lines(fitted ~ seq(0,12,0.1))

#Logistic Regression/Histogram Plot
Symcap$Color <- ifelse(Symcap$Color.Morph=="Orange", 1, 0)
logi.hist.plot(independ = Symcap$Depth..m., depend = Symcap$Color, type = "hist", boxp = FALSE, ylabel = "", col="gray", ylabel2 = "", xlabel = "Depth (m)")
mtext(side = 4, text = "Frequency", line = 3, cex=1)
mtext(side = 4, text = "Brown                                Orange", line = 2, cex = 0.75)
mtext(side = 2, text = "Probability of Orange Color Morph", line = 3, cex = 1)

Dom1 <- subset(Symcap, !is.na(Depth..m.) & !is.na(Dominant))
logi.hist.plot(Dom1$Depth..m., Dom1$Dominant, boxp = FALSE, type = "hist", col="gray", xlabel = "Depth (m)", ylabel = "", ylabel2 = "")
mtext(side = 4, text = "Frequency", line = 3, cex=1)
mtext(side = 4, text = "C                                       D", line = 2, cex = 0.75)
mtext(side = 2, text = "Probability of clade D Symbiont", line = 3, cex = 1)

#Export Image
pdf(file="Mix~Color Morph", height = 4, width = 6)
Symcap$Reef.Area <- ifelse(Symcap$Reef.Area!="Top", yes = "Slope", no = "Top")
results=table(Symcap$Mix, Symcap$Color.Morph)
results
chisq.test(results)
prop.table(results, margin = 2)
par(mar=c(4, 4, 2, 6))
barplot(prop.table(results, margin = 2), col = c("gray10", "gray40", "gray70", "gray100"), xlab = "Color Morph", ylab = "Symbiont Community Composition")
legend("topright", legend=c("C", "CD", "DC", "D"), fill=c("gray10", "gray40", "gray70", "gray100"), inset = c(-.2, 0), xpd = NA)
dev.off()

#RGoogleMaps
KB <- c(21.46087401, -157.809907) 
KBMap <- GetMap(center = KB, zoom = 13, maptype = "roadmap", SCALE = 2, GRAYSCALE = TRUE)
PlotOnStaticMap(KBMap, XY$Latitude, XY$Longitude, col=153, pch=21, bg="#7FFFD4", lwd=2)
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

rownames(XY) <- XY$Reef.ID
XY <- XY[, -1]
XY <- na.omit(XY)
apply(XY, MARGIN=1, FUN=function(reef) {
  floating.pie(xpos = reef["X"], ypos = reef["Y"], 
               x=c(reef["C"]), radius = 5, col = c("#7FFFD4"))
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
               x=c(reef["C"], reef["CD"], reef["DC"], reef["D"]), radius = 7, col = c("#0571b0","#92c5de","#f4a582","#ca0020"))
})

#Plot Dominant Symbiont vs. Depth and find threshold depth of D to C dominance
threshdepth <- function(reef) {
  df <- subset(Symcap, Reef.ID==reef)
  plot(df$Dominant~df$Depth..m., xlab="Depth (m)", ylab = "Proportion of Dominant Symbiont",
       main=reef)
  abline(h = 0.5, lty=2)
  results=glm(Dominant~Depth..m., family = "binomial", data = df)
  pval <- data.frame(coef(summary(results)))$`Pr...z..`[2]
  mtext(side=3, text=round(pval, 3))
  newdata <- list(Depth..m.=seq(0,12,0.01))
  fitted <- predict(results, newdata = newdata, type = "response")
  lines(fitted ~ seq(0,12,0.01))
  thresh <- ifelse(pval < 0.05,
                   newdata$Depth..m.[which(diff(sign(fitted - 0.5))!=0)], NA)
  return(thresh)
}

sapply(levels(Symcap$Reef.ID), FUN=threshdepth)
levels(Symcap$Reef.ID)

threshdepth("Deep")
threshdepth(42)
threshdepth("HIMB")
threshdepth(21)
threshdepth(46)
threshdepth(18)
threshdepth("F9-5")
threshdepth("F8-10")