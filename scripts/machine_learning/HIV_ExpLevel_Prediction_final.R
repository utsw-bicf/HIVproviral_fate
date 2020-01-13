### Developing an HIV expression level prediction model ###

CASE <- 2 #select a case {1,2,3}
#============================ Case descriptions ==================================#
#===== Case I: numerical ChrommHMM + numerical Lamin              ===============#
#===== Case II: numerical ChrommHMM + categorical Lamin           ===============#
#===== Case III: previous numerical ChrommHMM + categorical Lamin ===============#
#=================================================================================#

#==== Data import and preprocessing
library(stringr)
setwd("~/R/IvanDorso_lab")
inputData <- read.table("HIVexp_His7Plus4_chrom_lamin_damID_4ML.tsv", sep="\t", header = TRUE) #damID data added

newData <- inputData[, c(6:285)]
newID <- paste(inputData$chr, inputData$start, sep = '_')

# remove duplicate rows
idx_dup <- which(duplicated(newID))
newData <- newData[-idx_dup,]
rownames(newData) <- newID[-idx_dup]

# log10 transform for damID data columns
newData[,c(261:280)] <- data.frame(lapply(newData[,c(261:280)], function(x) {log(x+1e-3,10)}))

# replace '-' sign with 'neg' for column names 
colnames(newData) <- gsub("\\.", "neg", colnames(newData))

# convert ChromHMM and Lamin values into numeric values using lookup tables
tb_ChromHMM <- data.frame(Status=c("U1","U4","U3","U6","U2","U7","U5","U10","U8","U11","U12","U14","U13","U9", "U15"), 
                          Score=c(14:0)/14)
tb_Lamin <- data.frame(Protein=c("A1","A2","B1","B2","B3","B4"), 
                       Score=(c(2,2,-2,-1,-1,0)+2)/4)
if(CASE==1 | CASE==2){
  newData[,c(1:20)]  <- data.frame(lapply(newData[,c(1:20)], 
                                          function(x) {tb_ChromHMM[match(x, tb_ChromHMM$Status), "Score"]}))
}else if(CASE == 3) {
  newData[,c(1:20)] <- data.frame(lapply(newData[,c(1:20)], function(x) {as.numeric(as.character(gsub("U", "", x)))}))
} else {
  stop("Wrong case number is selected")
}

if(CASE==1) {
  newData[,c(21:40)] <- data.frame(lapply(newData[,c(21:40)], function(x) {tb_Lamin[match(x, tb_Lamin$Protein), "Score"]}))
}

# check data dimension
n_obs <- dim(newData)[1]
n_fea <- dim(newData)[2]

# scale HIVexp values and add an HIVexp column to the data matrix
library(dplyr)
newData$HIVexp <- scale(inputData$HIVexp[-idx_dup])
hist(newData$HIVexp, breaks=30)

# code HIVexp values into (0, 0.5, 1) which represent (Low, Intermediate, High) expression levels
newData <- newData %>% mutate(ExpLevel = case_when(HIVexp >= 0.5 ~ 1, HIVexp >= -0.5 ~ 0.5, TRUE ~ 0))
table(newData$ExpLevel)

#==== Create Training Data
newData_High <- newData[which(newData$ExpLevel == 1), ]
newData_Intermediate <- newData[which(newData$ExpLevel == 0.5), ]
newData_Low <- newData[which(newData$ExpLevel == 0), ]

set.seed(1234)
splitRatio <- 0.75
High_train_index <- sample(1:nrow(newData_High), splitRatio*nrow(newData_High))
Intermediate_train_index <- sample(1:nrow(newData_Intermediate), splitRatio*nrow(newData_Intermediate))
Low_train_index <- sample(1:nrow(newData_Low), splitRatio*nrow(newData_Low))

train_High <- newData_High[High_train_index, ]  
train_Intermediate <- newData_Intermediate[Intermediate_train_index, ]
train_Low <- newData_Low[Low_train_index, ]

train_Data <- rbind(train_High, train_Low) 
#train_Data <- rbind(train_High, train_Intermediate, train_Low)

#==== Create Test Data
test_High <- newData_High[-High_train_index, ]
test_Intermediate <- newData_Intermediate[-Intermediate_train_index, ]
test_Low <- newData_Low[-Low_train_index, ]

test_Data <- rbind(test_High, test_Low) 
#test_Data <- rbind(test_High, test_Intermediate, test_Low)

#==== Feature selection
library(smbinning)
# segregate continuous and factor variables
if(CASE==1) {
  factor_vars <- vector()
  continuous_vars <- setdiff(colnames(newData), c("HIVexp", "ExpLevel"))
  iv_df <- data.frame(VARS=c(factor_vars, continuous_vars), IV=numeric(n_fea))  # init for IV results
} else if(CASE==2 | CASE==3) {
  factor_vars <- colnames(newData)[21:40]
  continuous_vars <- setdiff(colnames(newData)[-c(21:40)], c("HIVexp", "ExpLevel"))
  iv_df <- data.frame(VARS=c(factor_vars, continuous_vars), IV=numeric(n_fea))  # init for IV results
}

# compute IV for categorical vars
for(factor_var in factor_vars){
  smb <- smbinning.factor(train_Data, y="ExpLevel", x=factor_var)  # WOE table
  if(class(smb) != "character"){ # check if some error occured
    iv_df[iv_df$VARS == factor_var, "IV"] <- smb$iv
  }
}
  
# compute IV for continuous vars
for(continuous_var in continuous_vars){
  smb <- smbinning(train_Data, y="ExpLevel", x=continuous_var)  # WOE table
  if(class(smb) != "character"){  # any error while calculating scores.
    iv_df[iv_df$VARS == continuous_var, "IV"] <- smb$iv
  }
}
  
iv_df <- iv_df[order(-iv_df$IV), ]  # sort
iv_df <- within(iv_df, VARS <- factor(VARS, levels=iv_df$VARS))
iv_df
  
ggplot(data=iv_df, aes(x=VARS,y=IV)) + 
  geom_bar(position="dodge",stat="identity") + 
  coord_flip() + theme(axis.text.y=element_blank()) +
  ggtitle("Information Values") 


#==== Model training with selected features  
var_notImp <- iv_df$VARS[iv_df$IV<=0.2] # variables not important by Optimal Binning

train_Data2 <- train_Data[,!names(train_Data) %in% var_notImp]
test_Data2 <- test_Data[,!names(train_Data) %in% var_notImp]

pModel <- glm(ExpLevel ~ . -HIVexp, data=train_Data2, family=binomial(link="logit"))
summary(pModel)

#==== Evaluate the model with test dataset
predicted <- plogis(predict(pModel, test_Data2))  #or# predicted <- predict(pModel, testData, type="response")

# Decide on optimal prediction probability cutoff for the model
library(InformationValue)
optCutOff <- optimalCutoff(test_Data2$ExpLevel, predicted)[1] 
optCutOff

# Misclassification Error
misClassError(test_Data2$ExpLevel, predicted, threshold = optCutOff)

# ROC
plotROC(test_Data2$ExpLevel, predicted)

# Concordance
Concordance(test_Data2$ExpLevel, predicted)

# Sensitivity and specificity
InformationValue::sensitivity(test_Data2$ExpLevel, predicted, threshold = optCutOff)
InformationValue::specificity(test_Data2$ExpLevel, predicted, threshold = optCutOff)

# Confusion Matrix
InformationValue::confusionMatrix(test_Data2$ExpLevel, predicted, threshold = optCutOff)

# box plot for 'High' and 'Low' test data
library(ggplot2)
library(ggpubr)
df <- data.frame(cbind(test_Data2$ExpLevel, test_Data2$HIVexp, predicted))
colnames(df) <- c("ExpLevel", "ExpValue", "Predicted")
df$ExpLevel <- ifelse(df$ExpLevel==0, "Low", "High")
df$ExpLevel <- as.factor(df$ExpLevel)

ggplot(df, aes(x=ExpLevel, y=Predicted, fill=ExpLevel)) + 
  geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) +
  stat_compare_means(label.x = 1.3, label.y = 1.1)

# box plot including 'Intermediate' test data
predicted_Intermediate <- plogis(predict(pModel, test_Intermediate))  #or# predicted <- predict(pModel, testData, type="response")
df2 <- data.frame(rbind(cbind(test_Data2$ExpLevel, test_Data2$HIVexp, predicted), cbind(test_Intermediate$ExpLevel, test_Intermediate$HIVexp, predicted_Intermediate)))
colnames(df2) <- c("ExpLevel", "ExpValue", "Predicted")
df2 <- df2 %>% mutate(ExpLevel = case_when(ExpLevel == 1 ~ "High", ExpLevel == 0 ~ "Low", TRUE ~ "Intermediate"))
df2$ExpLevel <- factor(df2$ExpLevel, ordered = TRUE, levels = c("High", "Intermediate", "Low"))

my_comparisons=list(c("High","Low"), c("Intermediate", "Low"))
ggplot(df2, aes(x=ExpLevel, y=Predicted, fill=ExpLevel)) +
  geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) +
  geom_signif(comparisons=my_comparisons, method = "t.test", y_position=c(1.10, 1.05)) + 
  scale_fill_manual(values=c("#F8766D", "#00BA38", "#619CFF"))

print(sort(names(test_Data2)))


# scatter plot including 'Intermediate' test data
library("gridExtra")

df3 <- df2[df2$ExpLevel %in% c('High', 'Low'),]

scatterPlot <- ggplot(df2, aes(x=ExpValue, y=Predicted, color=ExpLevel)) + geom_point() +
  geom_smooth(method=lm, color='black') +
  scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF")) + 
  annotate("text", label="italic(R)==0.19", x=2.2, y=0.77,hjust = "left",parse=TRUE) +
  annotate("text", label="italic(p)==0.00018", x=2.2, y=0.73,hjust = "left",parse=TRUE) +
  theme(legend.position=c(0,1), legend.justification=c(0,1))  # +stat_cor(label.x = 1, label.y = 0.85) 
  

xdensity <- ggplot(df3, aes(x=ExpValue, fill=ExpLevel)) + 
  geom_density(alpha=.7) + 
  scale_fill_manual(values = c("#F8766D", "#619CFF")) + 
  theme(legend.position = "none")

ydensity <- ggplot(df3, aes(Predicted, fill=ExpLevel)) + 
  geom_density(alpha=.7) + 
  scale_fill_manual(values = c("#F8766D", "#619CFF")) + 
  theme(legend.position = "none")

blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
  )

#grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
scatterPlot
