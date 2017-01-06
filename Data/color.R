# color calling

# import each person's data
rc <- read.csv("Data/Collection Data/RC_color.csv", header=F)
colnames(rc) <- c("Colony", "RC")
rrw <- read.csv("Data/Collection Data/RRW_color.csv", header=F)
colnames(rrw) <- c("Colony", "RRW")
ti <- read.csv("Data/Collection Data/Coral_Collection.csv")
ti <- ti[,c("Colony", "Color.Morph")]
colnames(ti) <- c("Colony", "TI")
ti$TI <- ifelse(ti$TI=="Orange", "o", "b")

# merge data
df <- merge(merge(rc, rrw, by="Colony"), ti, by="Colony")
head(df)

# calculate percent agreement, pairwise combinations
prop.table(table(df$RC==df$RRW))
prop.table(table(df$RC==df$TI))
prop.table(table(df$RRW==df$TI))

# calculate percent unanimous calls among all observers
df$all <- ifelse(df$RC==df$RRW & df$RC==df$TI, 1, 0)
prop.table(table(df$all))

