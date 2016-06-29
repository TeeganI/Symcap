#Do means vary across these 3 sites? Main question

#green crabs
plot(green.crabs~site, data=anova.demo)
#2 steps to doing ANOVA

#model = aov(green.crabs~site, data=anova.demo)
aov(green.crabs~site, data=anova.demo)
Terms:
  site Residuals
Sum of Squares  1228.2778  181.7222
Deg. of Freedom         2        22

Residual standard error: 2.874039
Estimated effects may be unbalanced

#anova(model)
anova(model)
Analysis of Variance Table

Response: green.crabs
Df  Sum Sq Mean Sq F value   Pr(>F)    
site       2 1228.28  614.14   74.35 1.63e-10 ***
  Residuals 22  181.72    8.26                     
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> #really small p value here = 1.63e -10
  
  #can do this but not recommended because does not store the model, which you need
  anova(aov(green.crabs~site, data=anova.demo))

#assumptions of ANOVA, visual test of assumptions
plot(model)  #to get residual plot and Q-Q plot
#first plot is the residual plot, to see if the variances are homogenous/equal
#they should all have roughly the same spread, good here

#QQ plot is to see if the data is normally distributed

#if we do this and the data doesnt fit the assumptions, you would need to transform your data

#normality --> shapiro test
shapiro.test
shapiro.test(model$residual)

Shapiro-Wilk normality test

data:  model$residual
W = 0.9846, p-value = 0.9588
#p-value greater than 0.05 so accept null, which is that data are normal

#equal variances --> levene test
install.packages("car")
library(car)

leveneTest(model)
Levene's Test for Homogeneity of Variance (center = median)
      Df F value Pr(>F)
group  2  1.4847 0.2484
      22               
#greater than 0.05 so accept null that variances are homogenous

#assumptions are met
#run anova test
anova(model)


#snails
plot(snails~site, data=anova.demo)
model=aov(snails~site, data=anova.demo)
plot(snails~site, data=anova.demo)
> model=aov(snails~site, data=anova.demo)
> aov(model)
Call:
aov(formula = model)

Terms:
site Residuals
Sum of Squares   17.06778 237.97222
Deg. of Freedom         2        22

Residual standard error: 3.288909
Estimated effects may be unbalanced
> anova(model)
Analysis of Variance Table

Response: snails
Df  Sum Sq Mean Sq F value Pr(>F)
site       2  17.068  8.5339  0.7889 0.4668
Residuals 22 237.972 10.8169               
> 
#snail data meets the assumptions (just know that)

power.anova.test(groups=3, n=8, between.var=8.5339, within.var=10.8169, sig.level=0.05)
#groups = number of treatments
#n= number of replicates in each treatment
#between.var = mean square, variance betwee ngroups, get the number from table
#within groups= from table, residual number, error
power.anova.test(groups=3, n=8, between.var=8.5339, within.var=10.8169, sig.level=0.05)
power.anova.test(groups=3, power=0.9, between.var=8.5339, within.var=10.8169, sig.level=0.05)
    