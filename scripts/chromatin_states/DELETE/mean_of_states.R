### This script tries to find the optimal states in chromHMM
### Starts from the CompareModels output using the biggest state as the primary
### Calculates the mean of each state and plots
### Should find where the states level off
### From paper "Systematic mapping of chromatin state landscapes during mouse development"

setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/chromatin_states/chromHMM/states_all_data/')

df <- read.table("allbam_states.txt", header=T, sep="\t")
rownames(df) <- df[,1]
df <- df[,2:(NCOL(df))]

### Calculate the median correlation
test2 <- as.data.frame(apply(df, 2, FUN = median))
test2$state <- c(5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)
plot(test2$state,test2$`apply(df, 2, FUN = median)`)

#### Use kmeans
df2 <- data.frame()
for (i in 5:19) {
  du <- kmeans(df,i,iter.max=100)
  ratio <- du$betweenss/du$totss
  list <- as.data.frame(matrix(c(i,ratio), ncol = 2))
  df2 <- rbind(df2, list)
}