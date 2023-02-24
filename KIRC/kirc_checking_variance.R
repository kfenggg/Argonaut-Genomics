library(matrixStats)
data = readRDS('kirc_824.RDS') # [1] 4022634     824
class(data) # this is matrix when we read it in 

# get the following for each feature
variance <- rowVars(data)
rowIndex <- c(1:nrow(data))

# compare variance here to python variance
data.var.filtered <- fread('kirc_824_varFiltered_cutoff_4.csv')
dim(data.var.filtered) # [1] 9693  826, 2 extra columns for featureNames and var
data.var.filtered[data.var.filtered$featureNames == '15939']$featureVar == variance[15939]

# they are not the same, we will proceed with R variance calculations
quantile(variance, probs = c(0.05, 0.5, 0.95))
max(variance) # this looks the same as the python output
sum(variance > 4) # this also gives the same number as python

# test on small subset to see if constructing a df works
test <- as.data.frame(data[1:5, 1:5])
test$idx <- rowIndex[1:5]
test$var <- variance[1:5]
test # works!

# construct whole df
data.df <- as.data.frame(data)
data.df$idx <- rowIndex
data.df$var <- variance
dim(data.df) # [1] 4022634     826
data.df[1:5, 820:826]

# take a look at the sample info 
sample.info = readRDS('sample_info.RDS')


# this should be done after splitting the dataset into norm v tumor 
mean <- rowMeans(data)
q5 <- rowQuantiles(data, probs = 0.05)
q95 <- rowQuantiles(data, probs = 0.95)
