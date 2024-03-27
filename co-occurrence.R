Cor2 <- data.frame(matrix(NA, ncol = length(prevalence_df$rowname)+1, nrow = length(prevalence_df$rowname)))
#Cor <- data.frame(matrix(NA, ncol = length(df$Contig)+1, nrow = length(df$Contig)))# create empty df
colnames(Cor2) <- c("Contig", prevalence_df$rowname)
Cor2$Contig <- prevalence_df$rowname

for (i in 1:length(Cor$Contig)){
  for (j in 1:length(Cor$Contig)){
    if (i < j){ 
      next
    }  
    Cor2[j,i+1] <- cor(as.numeric(prevalence_df[i,-1]), as.numeric(prevalence_df[j,-1]))
  }
}

# Parallelization
library(doParallel)
library(foreach)

# Set up a parallel backend with a specified number of cores
num_cores <- detectCores()-2
registerDoParallel(cores = num_cores)

# Create an empty dataframe Cor2
Cor2 <- data.frame(matrix(NA, ncol = length(prevalence_df$rowname) + 1, nrow = length(prevalence_df$rowname)))
colnames(Cor2) <- c("Contig", prevalence_df$rowname)
Cor2$Contig <- prevalence_df$rowname

# Define a function for calculating correlations
correlation_function <- function(i, j, prevalence_df) {
  cor(as.numeric(prevalence_df[i, -1]), as.numeric(prevalence_df[j, -1]))
}

# Parallelize the nested loop using foreach
results <- foreach(i = 1:length(Cor2$Contig), .combine = rbind) %:%
  foreach(j = 1:length(Cor2$Contig), .combine = c) %dopar% {
    if (i >= j) {
      correlation_function(i, j, prevalence_df)
    } else {
      NA
    }
  }

# Stop the parallel backend
stopImplicitCluster()

# Update Cor2 with the correlation values
Cor2[, 2:(length(prevalence_df$rowname) + 1)] <- results
