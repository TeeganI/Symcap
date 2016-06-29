data=read.csv(file="anova demo.csv")

plot(mussels~site,data=data) #initial check of differences
#fit anova model:
model=aov(mussels~site,data=data)
#visual check of assumptions
plot(model)
#check assumptions via tests
shapiro.test(model$residuals)
library(car)
leveneTest(model)  #neither assumption is met so need to transform...

#transform data:
#log transformation
data$log.mussels=log(data$mussels)
#refit model to transformed data
modeltrans=aov(log.mussels~site,data=data)
#visually check assumptions, see if the data looks "better"
plot(modeltrans)
#check assumptions again
shapiro.test(modeltrans$residuals)
leveneTest(modeltrans)  #both assumptions now met

anova(modeltrans)

#levenetest via anova on residuals
data$residuals=abs(model$residuals)
anova(aov(residuals~site,data=data))
leveneTest(modeltrans, center=mean) #equivalent results as anova on residuals

#sqrt transform
data$sqrt.mussels=sqrt(data$mussels)

#arcsine
data$a.mussels=asin(sqrt(data$mussels/100))

#pairwise posthoc to compare differences in means
#tukeyHSD posthoc test
(tukey=TukeyHSD(modeltrans))

#assign letters to the comparisons
library(multcompView)
multcompLetters(tukey$site[, "p adj"])

#grah results with letters 
library(sciplot)
bar=bargraph.CI(x.factor=site, response=mussels, data=data,
                ylab="Mean density mussels", xlab="Site", ylim=c(0,70))
box()
#add letters
means=aggregate(mussels~site,data=data, FUN=mean)
ses=aggregate(mussels~site,data=data, FUN=se)
text(x=bar$xvals, y=means$mussels+ses$mussels, pos=3, labels=c("a", "b", "b"))

#Results
#calculate effect size to describe biological significance:
The mean density of mussels in CT was significantly different than the mean density of mussels
in ME and MA (ANOVA, F= 5.094, df=4,20, p<0.05). 
The mean density in CT was ... greater than the density in ME (Tukey HSD, p<0.001) 
and ... greater than the density in MA (Tukey HSD, p<0.001).
The densities in ME and MA were not significantly different from each other
(Tukey HSD, p>0.05).