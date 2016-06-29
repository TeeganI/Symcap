library(data.table)

Coral_Data <- read.csv("Coral_Collection.csv")
Coral_Data$Depth..m. <- as.numeric(as.character(Coral_Data$Depth..m.))
hist(Coral_Data$Depth..m.)
Leeward=subset(Coral_Data, Reef.Area=="Leeward")
hist(Leeward$Depth..m.)
Windward=subset(Coral_Data, Reef.Area=="Windward")
hist(Windward$Depth..m.)
subset(Coral_Data, Reef.Area!="Top")
Slope <- subset(Coral_Data, Reef.Area!="Top")
hist(Slope$Depth..m.)

GPSData <- read.csv("GPS_Data.csv", skip = 27, stringsAsFactors = F)
GPSData$ID <- as.numeric(as.character(GPSData$ID))
GPSData <- GPSData[GPSData$ID>158,]
View(GPSData)

stopat <- which(GPSData$ID=="")[1]  # Find first row where ID is blank...
GPSData <- GPSData[1:stopat-1, ]  # Remove all rows past the first blank row...

# Convert time from character object to time object
GPSData$time <- as.POSIXct(GPSData$time, format="%Y-%m-%dT%H:%M:%SZ", tz="GMT")
GPSData$hst <- format(GPSData$time, tz="Pacific/Honolulu")
