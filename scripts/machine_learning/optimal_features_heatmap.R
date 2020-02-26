########## This script makes a heatmap of HIV insertions divided by H, M, L and
########## Lamin sub-comparments (A1, A2, B1, B2, B3, B4)

##### High: 455, Intermediate: 753, Low:351
### Divide HIV by H, M, L
setwd('/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results')
library(pheatmap)

### Load the data
HIV_high <- read.table('/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_high_Lamin.bed', sep = "\t")
HIV_inter <- read.table('/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_intermediate_Lamin.bed', sep = "\t")
HIV_low <- read.table('/project/BICF/BICF_Core/shared/Projects/Dorso/machine_learning/heatmaps_from_results/HIV_IS_low_Lamin.bed', sep = "\t")

### Add variable to df
HIV_high$exp <- "High"
HIV_inter$exp <- "Intermediate"
HIV_low$exp <- "Low"

### Combind df's and remove unneccessary celss
df <- rbind(HIV_high, HIV_inter)
df <- rbind(df, HIV_low)

testdf <- data.frame(table(df$exp,df$V10))
library(reshape2)
testdf2 <- data.frame(dcast(testdf, Var1~Var2))
testdf3 <- sweep(testdf2[,2:ncol(testdf2)], MARGIN = 1, FUN="/", STATS=rowSums(testdf2[,2:ncol(testdf2)]))
rownames(testdf3) <- c("High", "Intermediate", "Low")


pheatmap(testdf3, cluster_rows = F, cluster_cols = F, display_numbers = round(testdf3,2), fontsize = 20, cellwidth=50, cellheight=50, filename = "HIV_lamin_optimalfeatures_heatmap.pdf")
dev.off()

### Get heatmap of mean by group
df3r <- df[,c(11, 10, 5)]

newdf <- aggregate(x=df3r$V5, by = list(df3r$exp, df3r$V10), FUN = mean)
newdf2 <- data.frame(dcast(newdf, Group.1~Group.2))
rownames(newdf2) <- c("High", "Intermediate", "Low")
newdf3 <- newdf2[-c(1)]
pheatmap(newdf3, cluster_rows = F, cluster_cols = F, display_numbers = round(newdf3,2), fontsize = 20, cellwidth=50, cellheight=50,filename = "HIV_lamin_optimalfeatures_heatmap_meanexp.pdf")
dev.off()
