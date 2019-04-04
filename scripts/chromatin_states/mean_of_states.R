### This script tries to find the optimal states in chromHMM
### Starts from the CompareModels output using the biggest state as the primary
### Calculates the mean of each state and plots
### Should find where the states level off
### From paper "Systematic mapping of chromatin state landscapes during mouse development"

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/')

df <- read.table("states.txt", header=T, sep="\t")
rownames(df) <- df[,1]
df <- df[,2:(NCOL(df))]


### Calculate the mean
df$mean <- rowMeans(df)
test <- as.data.frame(colMeans(df))
test$state <- c(5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)
plot(test$state,test$`colMeans(df)`)

### Calculate the median
test2 <- as.data.frame(apply(df, 2, FUN = median))
test2$state <- c(5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)
plot(test2$state,test2$`apply(df, 2, FUN = median)`)

test13 <- kmeans(df,13,iter.max=100)
test9 <- kmeans(df,9,iter.max=100)
test10 <- kmeans(df,10,iter.max=100)