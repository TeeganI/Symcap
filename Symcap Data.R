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

#import gps data
GPSData <- read.csv("GPS_Data.csv", skip = 27, stringsAsFactors = F)
GPSData$ID <- as.numeric(as.character(GPSData$ID))
GPSData <- GPSData[GPSData$ID>158,]
View(GPSData)
stopat <- which(GPSData$ID=="")[1]  # Find first row where ID is blank...
GPSData <- GPSData[1:stopat-1, ]  # Remove all rows past the first blank row...

# Convert time from character object to time object
GPSData$time <- as.POSIXct(GPSData$time, format="%Y-%m-%dT%H:%M:%SZ", tz="GMT")
GPSData$hst <- format(GPSData$time, tz="Pacific/Honolulu")

# Importing qPCR data
source_url("https://raw.githubusercontent.com/jrcunning/steponeR/master/steponeR.R")
Mcap.plates <- list.files(path = "qPCR_data", pattern = "txt$", full.names = T)
Mcap.plates
Mcap <- steponeR(files = Mcap.plates, delim="\t",
                 target.ratios=c("C.D"), 
                 fluor.norm = list(C=2.26827, D=0), 
                 copy.number=list(C=33, D=3))
Mcap <- Mcap$result
Mcap
