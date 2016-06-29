data=read.csv(file="demo island regression .csv")
#check that both variables are continuous
str(data)

#check assumptions
plot(species~area, data=data)
model=lm(species~area, data=data)
plot(model)
#1st plot is residuals and homogeneity of variances
#second plot is Q-Q normality
#Assumptions not met

#transform data
#can transform either or both of the variables, whatever gives you the best line
data$log.species=log(data$species)
data$log.area=log(data$area)
modeltrans=lm(log.species~log.area, data=data)
plot(modeltrans)
#assumptions met

#get anova table 
anova(modeltrans)

#get the R2 value and coefficents
summary(modeltrans)

#plot figure
plot(log.species~log.area, data=data, ylab="Species Richness", xlab="Island Area")
abline(modeltrans)

#results
#Species richness significantly increased with island size (linear regression;
#log(species)=0.043 + 0.330*log(area); F=21.044, df=1,15, p<0.001; R2=0.584).
#Island size explained 58% of the variation in species diversity.