KB <- c(21.46087401, -157.809907) 
KBMap <- GetMap(center = KB, zoom = 13, maptype = "satellite", 
                SCALE = 2, GRAYSCALE = FALSE)
save(KBMap, file = "KBMap.Rdata")
load("KBMap.Rdata")
North <- subset(Symcap, Bay.Area=="Northern")
Central <- subset(Symcap, Bay.Area=="Central")
South <- subset(Symcap, Bay.Area=="Southern")
LatitudeN=aggregate(Latitude~Reef.ID, data=North, FUN = mean)
LongitudeN=aggregate(Longitude~Reef.ID, data = North, FUN=mean)
XYN<-merge(LatitudeN, LongitudeN, by="Reef.ID", all=T)
newcoordsN <- LatLon2XY.centered(KBMap, XYN$Latitude, XYN$Longitude, zoom=13)
XYN$X <- newcoordsN$newX
XYN$Y <- newcoordsN$newY
XYN <- subset(XYN, Reef.ID!="37")

LatitudeC=aggregate(Latitude~Reef.ID, data=Central, FUN = mean)
LongitudeC=aggregate(Longitude~Reef.ID, data = Central, FUN=mean)
XYC<-merge(LatitudeC, LongitudeC, by="Reef.ID", all=T)
newcoordsC <- LatLon2XY.centered(KBMap, XYC$Latitude, XYC$Longitude, zoom=13)
XYC$X <- newcoordsC$newX
XYC$Y <- newcoordsC$newY
XYC <- subset(XYC, Reef.ID!="37")

LatitudeS=aggregate(Latitude~Reef.ID, data=South, FUN = mean)
LongitudeS=aggregate(Longitude~Reef.ID, data = South, FUN=mean)
XYS<-merge(LatitudeS, LongitudeS, by="Reef.ID", all=T)
newcoordsS <- LatLon2XY.centered(KBMap, XYS$Latitude, XYS$Longitude, zoom=13)
XYS$X <- newcoordsS$newX
XYS$Y <- newcoordsS$newY
XYS <- subset(XYS, Reef.ID!="37")

par(oma=c(3,3,0,0))

PlotOnStaticMap(KBMap, XYN$Latitude, XYN$Longitude, col=153, pch=21, 
                bg="#FF7F50", lwd=2, cex = 1.2)
PlotOnStaticMap(KBMap, XYC$Latitude, XYC$Longitude, col=153, pch=21, 
                bg="#CD9B1D", lwd=2, cex = 1.2)
PlotOnStaticMap(KBMap, XYS$Latitude, XYS$Longitude, col=153, pch=21, 
                bg="#CDB79E", lwd=2, cex = 1.2)

axis(1, at = LatLon2XY.centered(KBMap, NA, c(-157.85, -157.81, -157.77))$newX, tcl=0.5, line = 0.5, col = "ghostwhite", col.ticks = "black", lwd = 1, outer = TRUE, labels = c("157.85°W", "157.81°W", "157.77°W"), padj = -2.5, cex.axis = 0.75)
axis(2, at = LatLon2XY.centered(KBMap, c(21.42, 21.46, 21.50), NA)$newY, tcl=0.5, line = 0.5, col = "ghostwhite", col.ticks = "black", lwd = 1, outer = TRUE, labels = c("21.42°N", "21.46°N", "21.50°N"), padj = 0.5, hadj = 0.60, las = 1, cex.axis = 0.75)
par(new=T, mar=c(14,21,0,0))
HI <- readOGR("coast_n83.shp", "coast_n83") 
HI <- spTransform(HI, CRS("+proj=longlat +datum=NAD83")) 
plot(HI, xlim=c(-158.3, -157.6), ylim=c(21.35, 21.6), lwd=0.4, col="gray", bg="white")
rect(-157.87, 21.41, -157.75, 21.52)
box()
