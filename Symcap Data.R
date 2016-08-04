library(data.table)
library(devtools)
library(plyr)
library(reshape2)
library(popbio)
library(RgoogleMaps)

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

#Order by Colony
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

#Identify Symbiont clades present
Symcap$Mix <- factor(ifelse(Symcap$propC>Symcap$propD, ifelse(Symcap$propD!=0, "CD", "C"), ifelse(Symcap$propD>Symcap$propC, ifelse(Symcap$propC!=0, "DC", "D"), NA)), levels = c("C", "CD", "DC", "D"))

#Chi Squared test for independence
Symcap$Reef.Area <- ifelse(Symcap$Reef.Area!="Top", yes = "Slope", no = "Top")
total=table(Symcap$Color.Morph, Symcap$Reef.Area)
total
chisq.test(total)
prop.table(total, margin = 2)
par(mar=c(3, 4, 2, 6))
barplot(prop.table(total, margin = 2), col = c("gray25", "gray100"), xlab = "Proportion of Colonies", ylab = "Dominant Symbiont")
legend("topright", legend=c("Brown", "Orange"), fill=c("gray25", "gray100"), inset = c(-.2, 0), xpd = NA)

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
fitted <- predict(results, newdata = list(Depth..m.=seq(0,8,0.1)), type = "response")
lines(fitted ~ seq(0,8,0.1))

#Plot Dominant Symbiont and Depth
Symcap$Dominant <- ifelse(Symcap$Dom=="C", 0, 1)
plot(Symcap$Dominant~Symcap$Depth..m., xlab="Depth (m)", ylab = "Proportion of Dominant Symbiont")
results=glm(Dominant~Depth..m., family = "binomial", data = Symcap)
anova(results, test = "Chisq")
summary(results)
fitted <- predict(results, newdata = list(Depth..m.=seq(0,8,0.1)), type = "response")
lines(fitted ~ seq(0,8,0.1))

#Logistic Regression/Histogram Plot
logi.hist.plot(independ = Symcap$Depth..m., depend = Symcap$Color, type = "hist", boxp = FALSE, ylabel = "", col="gray", ylabel2 = "", xlabel = "Depth (m)")
mtext(side = 4, text = "Frequency", line = 3, cex=1)
mtext(side = 4, text = "Brown                                Orange", line = 2, cex = 0.75)
mtext(side = 2, text = "Probability of Orange Color Morph", line = 3, cex = 1)

Dom1 <- subset(Symcap, Dominant!="NA")
logi.hist.plot(Dom1$Depth..m., Dom1$Dominant, boxp = FALSE, type = "hist", col="gray", xlabel = "Depth (m)", ylabel = "", ylabel2 = "")
mtext(side = 4, text = "Frequency", line = 3, cex=1)
mtext(side = 4, text = "C                                     D", line = 2, cex = 0.75)
mtext(side = 2, text = "Probability of clade D Symbiont", line = 3, cex = 1)

#Export Image
pdf(file="Color~Depth.pdf", height = 4, width = 5)
par(mar=c(4, 4, 4, 4))
logi.hist.plot(Symcap$Depth..m., Symcap$Color, boxp = FALSE, type = "hist", col="gray", xlabel = "Depth (m)", ylabel = "", ylabel2 = "")
mtext(side = 4, text = "Frequency", line = 3, cex=1)
mtext(side = 4, text = "Brown                                   Orange", line = 2, cex = 0.75)
mtext(side = 2, text = "Probability of Orange Color Morph", line = 3, cex = 1)
dev.off()

#RGoogleMaps
KB <- c(21.47285, -157.82936) 
KBMap <- GetMap(center = KB, zoom = 15, maptype = "satellite", SCALE = 2)
PlotOnStaticMap(KBMap, Symcap$Latitude, Symcap$Longitude, col=c("red", "blue")[Symcap$Color.Morph], pch=c(1, 2)[Symcap$Dom])
