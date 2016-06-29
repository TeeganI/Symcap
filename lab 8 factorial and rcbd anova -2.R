#FACTORIAL ANOVA
data=read.csv(file="factorial anova demo data.csv")
str(data)
data$snails=as.factor(data$snails)
str(data)

#plot data to visualize
lineplot.CI(x.factor=snails, response=plant.height, group=distance, data=data, 
            x.leg="topright", lty=1, col=c("red", "blue"))

#fit anova model
model=aov(plant.height~snails*distance, data=data)

#check assumptions
shapiro.test(model$residuals) #met
leveneTest(model) #met

#run anova
anova(model)
#the interaction term is NOT significant, which means we need to check
#the assumptions for the main effects of snails and distance

leveneTest(plant.height~snails, data=data)
leveneTest(plant.height~distance, data=data)
shapiro.test(residuals(aov(plant.height~snails, data=data)))
shapiro.test(residuals(aov(plant.height~distance, data=data)))
#all assumptions are met

#run Tukey test to get any significant main effects
(tukey=TukeyHSD(model))
library(multcompView)
multcompLetters(tukey$snails[, "p adj"])

bargraph.CI(x.factor=snails, response=plant.height, group=distance, data=data, ylab="Mean Plant Height", xlab=
              "Number of Snails")
legend(x="topright", legend=c("10m", "3m"), fill=c("black", "grey"))
box()

#RCBD ANOVA
data=read.csv(file="RCBD demo data.csv")
str(data)
data$branch=as.factor(data$branch)
str(data)
#put treatments in logical order
data$treatment=factor(data$treatment, level=c("Low N", "Moderate N", "High N", "Very High N"))
#plot, visualize data
plot(spartina.biomass~treatment, data=data)
plot(spartina.biomass~branch, data=data) #see the variation
library(sciplot)
lineplot.CI(x.factor=treatment, response=spartina.biomass, group=branch, data=data, col=c("red", "blue", "orange", "green"))

#without block
mod.no.block=aov(spartina.biomass~treatment, data=data)
anova(mod.no.block)

#with block
mod.block=aov(spartina.biomass~treatment+branch, data=data)
anova(mod.block)

#check assumptions
leveneTest(mod.no.block) #cant run levene on mod.block because only one replicate per block..?
shapiro.test(mod.block$residuals)

#run tukey on mod.block to get significance
(tukey=TukeyHSD(mod.block))
multcompLetters(tukey$treatment[, "p adj"])

bargraph.CI(x.factor=treatment, response=spartina.biomass, data=data)
box()
bargraph.CI(x.factor=branch, response=spartina.biomass, data=data)
box()
