1,1,3,5,11

A=c(1,1,3,5,11)
length(A)
sum(A)
mean(A)
meanA=mean(A)

example=read.csv(file="demo_lab 1.csv")
plot(snails~site, data=example, ylab="Number of snails", col="blue")
example
head(example)
title(main="Snails at a site")

library(sciplot)
bargraph.CI(x.factor=tide.height, response=mussels, data=example, ylab="Mussels")
bargraph.CI(x.factor=tide.height, response=mussels, group=site,data=example, ylab="Mussels")
box()

means=aggregate(snails~tide.height, data=example, FUN=mean)

lowA=subset(example, tide.height=="low" & site=="A")

table(example$snails)

hist(example$mussels, breaks=20)
