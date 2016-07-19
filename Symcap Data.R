library(data.table)
library(devtools)
library(plyr)
library(reshape2)

#import collection data
Coral_Data <- read.csv("Coral_Collection.csv")
Coral_Data$Depth..m. <- as.numeric(as.character(Coral_Data$Depth..m.))
View(Coral_Data)

#analyzing depth distribution
hist(Coral_Data$Depth..m.)
Leeward=subset(Coral_Data, Reef.Area=="Leeward")
hist(Leeward$Depth..m.)
Windward=subset(Coral_Data, Reef.Area=="Windward")
hist(Windward$Depth..m.)
subset(Coral_Data, Reef.Area!="Top")
Slope <- subset(Coral_Data, Reef.Area!="Top")
hist(Slope$Depth..m.)

#import gps data
GPSData <- read.csv("GPS_Data.csv", skip = 27, stringsAsFactors = F)
GPSData$ID <- as.numeric(as.character(GPSData$ID))
GPSData <- GPSData[GPSData$ID>158,]
stopat <- which(GPSData$ID=="")[1]  # Find first row where ID is blank
GPSData <- GPSData[1:stopat-1, ]  # Remove all rows past the first blank row
View(GPSData)

# Convert time from character object to time object
GPSData$time <- as.POSIXct(GPSData$time, format="%Y-%m-%dT%H:%M:%SZ", tz="GMT")
GPSData$hst <- format(GPSData$time, tz="Pacific/Honolulu")
View(GPSData)

# Importing qPCR data
source_url("https://raw.githubusercontent.com/jrcunning/steponeR/master/steponeR.R")
Mcap.plates <- list.files(path = "qPCR_data", pattern = "txt$", full.names = T)
Mcap.plates
Mcap <- steponeR(files = Mcap.plates, delim="\t",
                 target.ratios=c("C.D"), 
                 fluor.norm = list(C=2.26827, D=0), 
                 copy.number=list(C=33, D=3))
Mcap <- Mcap$result
View(Mcap)

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
View(Mcap)

#When C or D detected in only 1 sample, set ratio to +/- Inf
Mcap$C.D[which(Mcap$C.reps<2)] <--Inf
Mcap$C.D[which(Mcap$D.reps<2)] <-Inf

#Order by Colony
Mcap <- Mcap[with(Mcap, order(Colony)), ]
View(Mcap)

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
total=table(Symcap$Color.Morph, Symcap$Dom)
total
chisq.test(total)

#Mosaic Plot
mosaicplot(total, ylab = "Color Morph", xlab = "Dominant Symbiont", main = "")
