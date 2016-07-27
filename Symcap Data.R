library(data.table)
library(devtools)
library(plyr)
library(reshape2)

#import collection data
Coral_Data <- read.csv("Coral_Collection.csv")
Coral_Data$Depth..m. <- as.numeric(as.character(Coral_Data$Depth..m.))

#analyzing depth distribution
hist(Coral_Data$Depth..m.)
Leeward=subset(Coral_Data, Reef.Area=="Leeward")
hist(Leeward$Depth..m.)
Windward=subset(Coral_Data, Reef.Area=="Windward")
hist(Windward$Depth..m.)
subset(Coral_Data, Reef.Area!="Top")
Slope <- subset(Coral_Data, Reef.Area!="Top")
hist(Slope$Depth..m.)

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

#Replace "Sample.Name" column with "Colony"
colnames(Mcap)[which(colnames(Mcap)=="Sample.Name")] <- "Colony"

#Remove failed samples
Mcap$fail <- ifelse(Mcap$C.reps < 2 & Mcap$D.reps < 2, TRUE, FALSE)
fails <- Mcap[Mcap$fail==TRUE, ]
Mcap <- Mcap[which(Mcap$fail==FALSE),]

#When C or D detected in only 1 sample, set ratio to +/- Inf
Mcap$C.D[which(Mcap$C.reps<2)] <--Inf
Mcap$C.D[which(Mcap$D.reps<2)] <-Inf

#Order by Colony
Mcap <- Mcap[with(Mcap, order(Colony)), ]

#Proportion C
Mcap$propC <- Mcap$C.D / (Mcap$C.D+1)

#Proportion D
Mcap$propD <- 1-Mcap$propC
Mcap$propD[which(Mcap$C.D==-Inf)] <-1
Mcap$propC[which(Mcap$C.D==-Inf)] <-0
Mcap$propD[which(Mcap$C.D==Inf)] <-0
Mcap$propC[which(Mcap$C.D==Inf)] <-1

#Dominant Symbiont type
Mcap$Dom <- ifelse(Mcap$propC>Mcap$propD, "C", "D")

#Merge datasets
Symcap<-merge(Coral_Data, Mcap, by="Colony", all=T)

#Chi Squared test for independence
total=table(Symcap$Color.Morph, Symcap$Reef.Area)
chisq.test(total)
total
Slope <- subset(Symcap, Reef.Area!="Slope")
Slope$Reef.Area <- droplevels(Slope$Reef.Area)
total=table(Slope$Reef.Area, Slope$Dom)
chisq.test(total)

#Mosaic Plot
mosaicplot(total, ylab = "Reef Area", xlab = "Color Morph", main = "")

#Logistic Regression/ANOVA
Symcap$Dom <- as.factor(Symcap$Dom)
results=glm(Dom~Depth..m., family = "binomial", data = Symcap)
anova(results, test = "Chisq")
summary(results)

plot(results)
 
exp(coef(results)) #how odds change
exp(confint.default(results)) #95% confidence interval
pi.hat=predict.glm(results, data.frame(Depth..m.=5.1), type = "response", se.fit = TRUE) #predict probability
pi.hat$fit
I.hat=predict.glm(results, data.frame(Depth..m.=5.1), se.fit = TRUE) #95% confidence interval for estimate
ci=c(I.hat$se.fit-1.96*I.hat$se.fit, I.hat$fit+1.96*I.hat$se.fit)
exp(ci)/(1+exp(ci)) #transform results to probabilities 
