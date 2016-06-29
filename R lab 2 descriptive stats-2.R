example=read.csv(file="demo_lab 1.csv")

#mean
mean(example$mussels)

#median = middle number in set of numbers
median(example$mussels)

#mode = the value that occurs most often
table(example$mussels)


#variance = sample variance
var(example$mussels)
var(example$snails)
#standard deviation = sq. root of variance = 
sqrt(var(example$mussels))
#or
sd(example$mussels)
#standard error (of mean)
sd(example$snails)/sqrt(length(example$snails))
#or
library(sciplot)
se(example$snails)

