# color calling

# import each person's data
rc <- read.csv("Data/Collection Data/RC_color.csv", header=F)
colnames(rc) <- c("Colony", "RC")
cw <- read.csv("Data/Collection Data/CW_color.csv", header=F)
colnames(cw) <- c("Colony", "CW")
rrw <- read.csv("Data/Collection Data/RRW_color.csv", header=F)
colnames(rrw) <- c("Colony", "RRW")
ti <- read.csv("Data/Collection Data/Coral_Collection.csv")
ti <- ti[,c("Colony", "Color.Morph")]
colnames(ti) <- c("Colony", "TI")
ti$TI <- ifelse(ti$TI=="Orange", "o", "b")

# merge data
df <- merge(merge(merge(rc, rrw, by="Colony"), ti, by="Colony"), cw, by="Colony")
df$RC <- ifelse(df$RC=="o", "Orange", "Brown")
df$RRW <- ifelse(df$RRW=="o", "Orange", "Brown")
df$TI <- ifelse(df$TI=="o", "Orange", "Brown")
df$CW <- ifelse(df$CW=="o", "Orange", "Brown")

# calculate percent agreement, pairwise combinations
#  Vector source for column combinations
n <- colnames(df[,-1])
#  Make combinations
id <- expand.grid(n, n)
#  Get result
out <- matrix(colSums(df[,-1][,id[,1]]==df[,-1][,id[,2]]), ncol=length(n))
diag(out) <- NA
colnames(out) <- n
rownames(out) <- n
# Convert to proportion
out/nrow(df)

# calculate percent unanimous calls among all observers
agree <- ifelse(df$RC==df$RRW & df$RC==df$TI & df$RC==df$CW, "all agree", "disagree")
prop.table(table(agree))

# calculate number of agreements per colony
nag <- apply(df[,-1], 1, function(x) max(table(x)))
nag
table(nag)
