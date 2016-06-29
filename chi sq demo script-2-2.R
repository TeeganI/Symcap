data=read.csv(file="chi sq demo data.csv")

#two continuous variables
#test for independence to see whether these variables are associated with each other

#2 steps:
#1 put your data in count form, make that 2-way contingency table
total=table(data$habitat, data$isopod.presence)
total
#2 run the test on the table, not your raw data
chisq.test(total)

mosaicplot(total, ylab="Isopod Presence", xlab="Habitat Type", main="")
box("outer")
#WHAT THIS MEANS: column widths are proportional to the column (or row) totals.
#Column heights are proportional to the cell frequency. 
#Thus the columns overall represent the relative frequency of that combination.
#width: 17/30 are eelgrass so slightly over half
#height: 12/17 are yes/eelgrass so higher

#A chi square test was performed to test whether isopod presence was `related to
#habitat type. Results showed that snail presence significanlty depended on habitat choice.
#(df=2, n=40, X=4.887, p<0.05), and was greater in eelgrass habitat.
