#calculate t critical value
qt(p=0.25, df=8,lower.tail=F) # default is TRUE, which is prob 
qt(p=0.25, df=8,lower.tail=T)

pt(2.39, df=8, lower.tail=F)

data=read.csv("demo rocky data.csv")
str(data)

par(mfcol=c(2,1), pin=c(4,1.5))
bargraph.CI(x.factor=Site, response=Mytilus, data=data)
box()
bargraph.CI(x.factor=Site, response=Mytilus, data=data)
box()
box("outer")

t.test(Mytilus~Site, data=data, var.equal=T)
t.test(Mytilus~Site, data=data)

power.t.test(n=16, delta=20.0562, sd=19.342, sig.level=0.05)
#.8 is high so thats a good power

#if you want higher power, need to know the number of replicates
power.t.test(n=16, delta=20.0562, sd=19.342, sig.level=0.05)
power.t.test(power=0.9, delta=20.0562, sd=19.342, sig.level=0.05)
#would need about 21 additional samples to get that increased power


#Results (t-test) example
#We found that the mean percent cover of Mytilus edulis at Cunner Ledge was significantly greater than the mean percent cover of mussels at Canoe Beach (standard two-sample t test, t=2.9328, df = 15.575, p-value = 0.00997).

