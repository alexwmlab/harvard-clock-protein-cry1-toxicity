
################################################################################
# Ha W.M.Alex   Dec 10,2024                                                    #  
#                                                                              #
# Harvard Capstone Project Two:                                                #
#                                                                              #
# Ensemble Model of Core Clock Protein Cryptochrome (CRY1) Toxicity Prediction #
#                                                                              #
################################################################################
options(warn=-1)

# Install Necessary Packages if required
#
if(!require(tidyverse)) install.packages("tidyverse", repos = "http://cran.us.r-project.org")
if(!require(tidyr)) install.packages("tidyr", repos = "http://cran.us.r-project.org")
if(!require(dplyr)) install.packages("dplyr", repos = "http://cran.us.r-project.org")
if(!require(lubridate)) install.packages("lubridate", repos = "http://cran.us.r-project.org" )
if(!require(stringr)) install.packages("stringr", repos = "http://cran.us.r-project.org")
if(!require(ggplot2)) install.packages("ggplot2", repos = "http://cran.us.r-project.org")
if(!require(gridExtra)) install.packages("gridExtra", repos = "http://cran.us.r-project.org")
if(!require(knitr)) install.packages("knitr", repos = "http://cran.us.r-project.org")
if(!require(rstudioapi)) install.packages("rstudioapi", repos = "http://cran.us.r-project.org")
if(!require(caret)) install.packages("caret", repos = "http://cran.us.r-project.org")
if(!require(tinytex)){install.packages("tinytex", repos = "http://cran.us.r-project.org")
  tinytex::install_tinytex()}
if(!require(reshape2)) install.packages("reshape2", repos = "http://cran.us.r-project.org")
if(!require(matrixStats)) install.packages("matrixStats", repos = "http://cran.us.r-project.org")
if(!require(pROC)) install.packages("pROC", repos = "http://cran.us.r-project.org")
if(!require(PRROC)) install.packages("PRROC", repos = "http://cran.us.r-project.org")
if(!require(gbm)) install.packages("gbm", repos = "http://cran.us.r-project.org")
if(!require(randomForest)) install.packages("randomForest", repos = "http://cran.us.r-project.org")
if(!require(xgboost)) install.packages("xgboost", repos = "http://cran.us.r-project.org")


# Loading Necessary Libraries
#
library(dslabs)
library(tidyverse)
library(dplyr)
library(tidyr)
library(caret)
library(knitr)
library(lubridate)
library(stringr)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(matrixStats)
library(pROC)
library(PRROC)
library(gbm)
library(randomForest)
library(xgboost)


# Set number of significant digits=4 globally
options(digits = 4)

# set maximum display output= 2000 globally
options(max.print=2000)

set.seed(755)


##############################################################################################
#                                                                                            #
#  Beginning - Ensemble Model of Core Clock Protein Cryptochrome (CRY1) Toxicity Prediction  #
#                                                                                            #
##############################################################################################


#####################################################################################
# Initialization:                                                                   #
#                                                                                   #
# 1) Download cryptochrome-d.csv and benchmark-13.csv Files                         #
#                                                                                   #
# 2) Read csv Files cryptochrome-d.csv and benchmark-13.csv from current directory. #
#                                                                                   #
#####################################################################################

#################################################
# Files Handling for R script current directory #
#################################################

# Get R script current directory
current_path_r_script <- rstudioapi::getActiveDocumentContext()$path
current_directory_r_script <- dirname(current_path_r_script)
print(current_directory_r_script)


# Join R script current directory path and file name "benchmark-13.csv" 
benchmark_file_path <- str_c(current_directory_r_script,"/benchmark-13.csv")
print(benchmark_file_path)

# Join R script current directory path and file name "cryptochrome-d.csv" 
cryptochrome_file_path <- str_c(current_directory_r_script,"/cryptochrome-d.csv")
print(cryptochrome_file_path)


# Download cryptochrome-d.csv File
download.file("https://github.com/alexwmlab/test-one/blob/main/cryptochrome-d.csv?raw=true",cryptochrome_file_path)

# Download benchmark-13.csv File 
download.file("https://github.com/alexwmlab/test-one/blob/main/benchmark-13.csv?raw=true",benchmark_file_path)


# Read csv File cryptochrome-d.csv and Create ccpc_data dataset
ccpc_data <- read_csv(col_names=TRUE, cryptochrome_file_path)

# Read csv File benchmark-13.csv and Create  benchmark_d dataset
benchmark_d <- read_csv(col_names=TRUE, benchmark_file_path)


#####################################
#  Data Exploration Analysis (EDA)  #
#####################################

# Dimension of ccpc_data (number of column and number of row)
ncol(ccpc_data)
nrow(ccpc_data)

# Dimension of ccpc_data - Number of rows and columns
dim_df <- data.frame(Rows = dim(ccpc_data)[1], Columns = dim(ccpc_data)[2])
dim_df %>% knitr::kable(caption="Dimension of ccpc_data set")                   


# Dimension of benchmark_d (number of column and number or row)
ncol(benchmark_d)
nrow(benchmark_d)

# Dimension of benchmark_d - Number of rows and columns
dim_df <- data.frame(Rows = dim(benchmark_d)[1], Columns = dim(benchmark_d)[2])
dim_df %>% knitr::kable(caption="Dimension of benchmark_d dataset")                   


#
# Description of Columns/Features and Data Type
#


# Display Summary Table of ccpc_data column and data type using kable function
#
data_type <- data.frame(Column="1 to 1023", Data_Type="Numeric", Value="")

data_type<- bind_rows(data_type,
                      data.frame(Column="Class", Data_Type="Character", Value="'Toxic' or 'NonToxic'"))

data_type  %>% knitr::kable(caption="ccpc_data set")


# Display Summary Table of benchmark_d column and data type using kable function
#
col_data_type <- data.frame(Data_Type=sapply(benchmark_d, class))
col_data_type %>% knitr::kable(caption="benchmark_d dataset")


#
#  Data Cleaning
#

# check and count number of missing value in ccpc_data set
na_results  <- ccpc_data[apply(is.na(ccpc_data),1,any),]
nrow(na_results)


#
# Handling Zero Values
#

# Create a logical matrix where TRUE represents zero value
zero_check <- ccpc_data == 0 

# Sum the TRUE values for each column to count zeros 
zero_count <- colSums(zero_check) 

# Display sum of zeros of ccpc_data set
sum(zero_count)


#
# Pre-Processing
#

# Display first 6 rows of column 99 descriptors within range [95:108]
head(ccpc_data[,99])

# Data Wrangling (e.g., MDEC-23 to MDEC_23)

# Data Wrangling of ccpc_data for Molecular Descriptors/Features/Columns
#   1) convert hyphens "-" to under_scores "_"
#
colnames(ccpc_data) <- str_replace_all(colnames(ccpc_data), "^'|'$","")
colnames(ccpc_data) <- str_replace_all(colnames(ccpc_data), "-","_")

# Data Wrangling of benchmark_d for Molecular Descriptors/Features/Columns
#   1) convert hyphens "-" to under_scores "_"
#
colnames(benchmark_d) <- str_replace_all(colnames(benchmark_d), "^'|'$","")
colnames(benchmark_d) <- str_replace_all(colnames(benchmark_d), "-","_")

# Display first 6 rows of column 99 descriptors within range [95:108]
head(ccpc_data[,99])


#####################
# Dataset Splitting #
#####################

#
# ccpc_data set:
#
# Column Name of Columns 1 to 1203: Contain Name of 1023 Molecular Descriptors
# Value of Columns 1 to 1023      : Contain Value of 1023 Molecular Descriptors
# value of Column "Class"         : Contain Toxicity Class "Toxic" or "NonToxic" of Molecules/Instances


#
# benchmark_d set:
#
# Column Name of Columns 1 to 13: Contain Name of 13 Redundant Reduced Molecular Descriptors
# Value of Columns 1 to 13      : Contain Value of 13 Redundant Reduced Molecular Descriptors
# Value of Column "Class"       : Contain Toxicity Class "Toxic" or "NonToxic" of Molecules/Instances


#
# df_ccpc_data set:
#
# Column Name of Columns 1 to 1023: Contain Name of 1203 Molecular Descriptors
# Value of Columns 1 to 1023      : Contain Value of 1023 Molecular Descriptors


#
# df_ccpc_tox set:
#
# value of Column "Class"         : Contain Toxicity Class "Toxic" or "NonToxic" of Molecules/Instances


#    
# df_benchmark_d set:
#
# Column Name of Columns 1 to 13: Contain Name of 13 Redundant Reduced Molecular Descriptors
# Value of Columns 1 to 13      : Contain Value of 13 Redundant Reduced Molecular Descriptors


#
# df_benchmark_d_tox set:
#
# Value of Column "Class"       : Contain Toxicity Class "Toxic" or "NonToxic" of Molecules/Instances


# Generate df_ccpc_data comprising 1203 Molecular Descriptors of 171 Molecules/Instances
df_ccpc_data <- ccpc_data  %>% select(-Class)

# Generate df_ccpc_tox comprising Toxicity Class ("Toxic" or "NonToxic") of 171 Molecules/Instances
df_ccpc_tox <- ccpc_data %>% select(Class)

# Generate df_benchmark_d comprising Molecular Descriptors for benchmarking
df_benchmark_d <- benchmark_d %>% select(-Class)

# Generate df_benchmark_d_tox comprising Toxicity Class ("Toxic" or "NonToxic") for benchmarking
df_benchmark_d_tox <- benchmark_d %>% select(Class)


# Dimension of df_ccpc_data (number of Column and number of Row)
ncol(df_ccpc_data)
nrow(df_ccpc_data)

# Table Display Dimension of df_ccpc_data using kable function of knitr packagae
dim_df <- data.frame(Rows = dim(df_ccpc_data)[1], Columns = dim(df_ccpc_data)[2])
dim_df %>% knitr::kable(caption="Dimension of df_ccpc_data set")                   


# Dimension of df_ccpc_tox (number of Column and number of Row)
ncol(df_ccpc_tox)
nrow(df_ccpc_tox)

# Table Display Dimension of df_ccpc_tox using kable function of knitr of package
dim_df <- data.frame(Rows = dim(df_ccpc_tox)[1], Columns = dim(df_ccpc_tox)[2])
dim_df %>% knitr::kable(caption="Dimension of df_ccpc_tox set")                   


# Dimension of df_benchmark_d (number of Column and number of Row)
ncol(df_benchmark_d)
nrow(df_benchmark_d)

# Table Display Dimension of df_benchmark_d using kable function of knitr package
dim_df <- data.frame(Rows = dim(df_benchmark_d)[1], Columns = dim(df_benchmark_d)[2])
dim_df %>% knitr::kable(caption="Dimension of df_benchmark_d set")                   


# Dimension of df_benchmark_d_tox (number of Column and number of Row)
ncol(df_benchmark_d_tox)
nrow(df_benchmark_d_tox)

# Table Display Dimension of df_benchmark_d_tox using kable function of knitr package
dim_df <- data.frame(Rows = dim(df_benchmark_d_tox)[1], Columns = dim(df_benchmark_d_tox)[2])
dim_df %>% knitr::kable(caption="Dimension of df_benchmark_d_tox set")                   


#
# Redundant Features Reducing - benchmark set
#

###############################################################################
#  Eliminate highly correlated features of df_benchmark_d (coefficient > 90%) #
###############################################################################

#
# Eliminate highly correlated Features of df_benchmark_d with coefficient over 90%.
#
#   * "findCorrelation" is a function of Caret package in R, used to detect and remove highly correlated variables
#     in  df_benchmark_d dataset.
#
#   * cutoff: A Parameter/Threshold for correlation. Variables with correlation higher than this value will be flagged.
#

# Compute Correlation Coefficient
benchmark_correlation <- cor(df_benchmark_d)

# Find highly correlated features (cutoff = 0.9)
highly_correlated <- findCorrelation(benchmark_correlation, cutoff = 0.9)

# Remove the highly correlated features
df_benchmark_d_reduced  <- df_benchmark_d[, -highly_correlated]

# Assign df_benchmark_d_tox to Reduced toxicity dataset 
df_benchmark_d_tox_reduced <- df_benchmark_d_tox

# Verify Dimensions
cat("Original df_benchmark_d dimensions:", dim(df_benchmark_d), "\n")
cat("Reduced df_benchmark_d dimensions:", dim(df_benchmark_d_reduced), "\n")

cat("Original df_benchmark_d_tox dimensions:", dim(df_benchmark_d_tox), "\n")
cat("Reduced df_benchmark_d_tox dimensions:", dim(df_benchmark_d_tox_reduced), "\n")


# Number of Molecular Descriptors/Features/Predictors in Redundant Features Reduced benchmarking set
ncol(df_benchmark_d_reduced)

# Column Name of Molecular Descriptors/Features/Predictors in Redundant Features Reduced benchmarking set
colnames(df_benchmark_d_reduced)

# Display Redundant Features Reduced benchmarking set - df_benchmark_d_reduced
tail(df_benchmark_d_reduced)


#
# Feature Selection - Benchmark set for Molecular Descriptors Prioritization 
#


######################################################################
#  Generate  benchmark_core_clock_protein_cry set for normalization  #
######################################################################

# Transform df_benchmark_d_reduced to moecular_descriptors matrix
molecular_descriptors <- df_benchmark_d_reduced %>% as.matrix()

# Transform df_benchmark_d_tox_reduced$Class to toxicity vector
toxicity <- factor(df_benchmark_d_tox_reduced$Class)

# Generate list - benchmark_core_clock_protein_cry
benchmark_core_clock_protein_cry <- list(molecular_descriptors=molecular_descriptors, 
                                         toxicity=toxicity )

# Number of Molecules/Instances/Column  are in the Matrix
nrow(benchmark_core_clock_protein_cry$molecular_descriptors)


# Number of Molecular Descriptors/Features/Predictors are in the Matrix
ncol(benchmark_core_clock_protein_cry$molecular_descriptors)


# Proportion of the Molecules/Instances are "Toxic"
mean(benchmark_core_clock_protein_cry$toxicity=="Toxic")


# Proportion of the Molecules/Instances are "NonToxic"
mean(benchmark_core_clock_protein_cry$toxicity=="NonToxic")


###################################################################
#  Scaling and Normalizing  benchmark_core_clock_protein_cry set  #
###################################################################

# Compute Mean of each Column
col_means <- colMeans(benchmark_core_clock_protein_cry$molecular_descriptors, na.rm=TRUE)

# Compute Standard Deviations of each Column
col_sds <- colSds(benchmark_core_clock_protein_cry$molecular_descriptors, na.rm=TRUE)

# Subtract Column Means
benchmark_molecular_descriptors_centered <- sweep(benchmark_core_clock_protein_cry$molecular_descriptors,
                                                  2, col_means, FUN="-")

# Divide by Column Standard Deviations
benchmark_molecular_descriptors_scaled <- sweep(benchmark_molecular_descriptors_centered,
                                                2 , col_sds, FUN="/")

# List 6 rows of scaled and normalized "benchmark_molecular_descriptors_scalad" set 
head(benchmark_molecular_descriptors_scaled)


#######################################################################################
#  Splits dataset into train_set_x, train_set_y (80%) & test_set_x, test_set_y (20%)  #
#######################################################################################

# Partitions the dataset into 80% training and 20% testing sets as follow:
#
# * Splits 80% of "benchmark_molecular_descriptors_scaled" Matrix into "train_set_x"
# * Splits 20% of "benchmark_molecular_descriptors_scaled" Matrix into "test_set_x"
#
# * Splits 80% of "benchmark_core_clock_protein_cry$toxicity" Class into "train_set_y"
# * Splits 20% of "benchmark_core_clock_protein_cry$toxicity" Class into "test_set_y"
#
#
# Given our dataset's small size , we allocate 20% for testing rather than the typical 10%.
# This ensures an adequate amount of data for a reliable and effective evaluation.


set.seed(1)

test_index <- createDataPartition(benchmark_core_clock_protein_cry$toxicity, times = 1, p = 0.2, list = FALSE)

test_set_x   <-  benchmark_molecular_descriptors_scaled[test_index,]
test_set_y   <-  benchmark_core_clock_protein_cry$toxicity[test_index]

train_set_x  <- benchmark_molecular_descriptors_scaled[-test_index,]
train_set_y  <- benchmark_core_clock_protein_cry$toxicity[-test_index]


#
#Feature Importance of Benchmark set
#

##############################################################################################################
# Select Most Relevant Molecular Descriptors from Vars Importance as Benchmark using Random Forest Algorithm #
##############################################################################################################

set.seed(9)


####################################################
# Optimized Parameters of Random Forest Algorithm  #
#    * Sequence "mtry" from 1 to 10                #     
#    * Set Number of "tree" equal 500  (ntree=500) #
#    * Set Number of "cross validation" equal to 5 #
####################################################

train_ccpc_benchmark <- caret::train(train_set_x, train_set_y , method="rf",
                                     tuneGrid=data.frame(mtry=seq(1:10)),
                                     ntree=500,
                                     trControl=trainControl(method="cv", number=5),
                                     importance= TRUE)



# Process and Generate Vars Importance of train_ccpc_benchmark
importance_rf_b <- varImp(train_ccpc_benchmark)
print(importance_rf_b)

# Vars Importance of benchmark set
benchmark_importance_rf <- rownames(importance_rf_b$importance)[order(-importance_rf_b$importance[,1])]
#benchmark_importance_rf[1:12]


# Display Top Pairwise Most Relevant Molecular Descriptors/Features as Benchmark.
paste(benchmark_importance_rf[1]," and ",benchmark_importance_rf[2],
      " are selected as Benchmark Descriptors for Features/Molecular Descriptors Prioritization")

# max vars importance
most_importance_rf <-rownames(importance_rf_b$importance)[which.max(importance_rf_b$importance[,1])]
#most_importance_rf


# Compute Overall Accuracy for Molecular Descriptors/Feature Prioritization purpose 
predict_ccpc_benchmark <- predict(train_ccpc_benchmark, test_set_x)
accuracy_ccpc_benchmark <- confusionMatrix(table(predict_ccpc_benchmark, test_set_y),mode="everything")
#accuracy_ccpc_benchmark
benchmark_accuracy <- accuracy_ccpc_benchmark$overall["Accuracy"]
benchmark_accuracy

# Display Overall Accuracy
paste(round(accuracy_ccpc_benchmark$overall["Accuracy"],4),
      " will be adpoted as Benchmark Accuracy for Features/Molecular Descriptors Prioritization")


#
# Principal Component Analysis (PCA) - "core_clock_protein_cry" set
#

#############################################################################
#   Generate core_clock_protein_cry set (1203 features) for normalization   #
#############################################################################

# Transform df_ccpc_data to molecular_descriptors Matrix
molecular_descriptors <- df_ccpc_data %>% as.matrix()

# Transform df_ccpc_tox$Class to toxicity vector as factor
toxicity <- factor(df_ccpc_tox$Class)

# Generate list - core_clock_protein_cry
core_clock_protein_cry <- list(molecular_descriptors=molecular_descriptors, 
                               toxicity=toxicity )
# core_clock_protein_cry

# Number of Molecules/Instances/Column  are in the Matrix
nrow(core_clock_protein_cry$molecular_descriptors)

# Number of Molecular Descriptors/Features/Predictors are in the Matrix
ncol(core_clock_protein_cry$molecular_descriptors)

# Proportion of the Molecules/Instances are "Toxic"
mean(core_clock_protein_cry$toxicity=="Toxic")

# Proportion of the Molecules/Instances are "NonToxic"
mean(core_clock_protein_cry$toxicity=="NonToxic")


#
# Scaling and Normalizing "core_clock_protein_cry" set
#

###########################################################
#  Scaling and Normalizing  "core_clock_protein_cry" set  #
###########################################################

# Compute Mean of each Column
col_means <- colMeans(core_clock_protein_cry$molecular_descriptors, na.rm=TRUE)

# Compute Standard Deviations of each Column
col_sds <- colSds(core_clock_protein_cry$molecular_descriptors, na.rm=TRUE)

# Subtract Column Means
molecular_descriptors_centered <- sweep(core_clock_protein_cry$molecular_descriptors,
                                        2, col_means, FUN="-")

# Divide by Column Standard Deviations
molecular_descriptors_scaled <- sweep(molecular_descriptors_centered,
                                      2 , col_sds, FUN="/")

# Examples 6 rows and Columns 1 to 10  of scaled and normalized "molecular_descriptors_scalad" set 
head(molecular_descriptors_scaled[,170:178])


######################################
#    PCA - Proportion of variance    #
######################################

# PCA is performed using the "prcomp" function
pca_result_core_clock_protein_cry <- prcomp(molecular_descriptors_scaled)

# Create Dataframe with column pca_result_core_clock_protein$x
pca_data_core_clock_protein_cry <- data.frame(pca_result_core_clock_protein_cry$x)

# Add toxicity Class to the Dataframe
pca_data_core_clock_protein_cry$toxicity <- core_clock_protein_cry$toxicity

# Compute Proportion of Variance (Eigenvalues and Eigenvectors)
proportions_var_ccpc <-pca_result_core_clock_protein_cry$sdev^2/sum(pca_result_core_clock_protein_cry$sdev^2)
#proportions_var_ccpc[1:42]

# Proportion of Variance explained by First Principal Component (PC1)
proportions_var_ccpc_first_pc <- proportions_var_ccpc[1]
proportions_var_ccpc_first_pc

# Number of Principal Components are required to explain at least 90% of Proportion of  variance
cumsum_proportions_var_ccpc   <- cumsum(proportions_var_ccpc)
no_pc_90percent_var_ccpc      <- which(cumsum_proportions_var_ccpc >= 0.9)[1]
paste(no_pc_90percent_var_ccpc, "Components are required to explain at least 90% of the variance")

# List of 42 Proportion Variance
paste("List of",no_pc_90percent_var_ccpc,"Proportion Variance")
proportions_var_ccpc[1:no_pc_90percent_var_ccpc]


###########################################
#   PCA - plotting Principal components   #
###########################################

# Create data frame for plotting 
plot_df_core_clock_protein_cry <-data.frame(PC1=pca_result_core_clock_protein_cry$x[,1],
                                            PC2=pca_result_core_clock_protein_cry$x[,2],
                                            Toxicity=core_clock_protein_cry$toxicity)


# ggplot First Two Principal Components (PC1 & PC2) with color representing Toxicity
ggplot(plot_df_core_clock_protein_cry,aes(x=PC1, y=PC2, color=Toxicity)) +
  geom_point() +
  labs(x="PC1", y="PC2", color="Toxicity") +
  ggtitle("Core Clock Protein Cryptochrome (CRY1) Toxicity: Scatter plot of PC1 & PC2")


################################
#  PCA - Plotting PC Box-plot  #
################################

plot_df_10 <- data.frame(pca_result_core_clock_protein_cry$x[,1:no_pc_90percent_var_ccpc])
plot_df_10$Toxicity <- core_clock_protein_cry$toxicity

melted_df_10 <- melt(plot_df_10, id.vars="Toxicity")

ggplot(melted_df_10, aes(x=variable, y=value, fill=Toxicity)) +
  geom_boxplot() +
  labs(x="Principal Component", y="Value", fill="Toxicity") +
  theme(axis.text.x=element_text(angle=90)) +
  ggtitle("Core Clock Protein Cryptochrome (CRY1) Toxicity: Principal Components")



#######################################################################################
#  Splits dataset into train_set_x, train_set_y (80%) & test_set_x, test_set_y (20%)  #
#######################################################################################

# Partitions the dataset into 80% training and 20% testing sets as follow:
#
# * Splits 80% of "molecular_descriptors_scaled" Matrix into "train_set_x"
# * Splits 20% of "molecular_descriptors_scaled" Matrix into "test_set_x"
#
# * Splits 80% of "core_clock_protein_cry$toxicity" Class into "train_set_y"
# * Splits 20% of "core_clock_protein_cry$toxicity" Class into "test_set_y"
#
#
# Given our dataset's small size but high dimensionality, we allocate 20% for testing rather than the typical 10%.
# This ensures an adequate amount of data for a reliable and effective evaluation.


set.seed(1)

test_index <- createDataPartition(core_clock_protein_cry$toxicity, times = 1, p = 0.2, list = FALSE)

test_set_x   <-  molecular_descriptors_scaled[test_index,]
test_set_y   <-  core_clock_protein_cry$toxicity[test_index]

train_set_x  <- molecular_descriptors_scaled[-test_index,]
train_set_y  <- core_clock_protein_cry$toxicity[-test_index]


# Proportion of the train_set_y is "Toxic"
#train_set_y
mean(train_set_y =="Toxic")

# Proportion of the test_set_y is "NonToxic"
#test_set_y
mean(test_set_y  =="NonToxic")


#
# Feature Importance - core_clock_protein_cry set
#

##############################################################################################################
# Select Most Relevant Molecular Descriptors from Vars Importance as Benchmark using Random Forest Algorithm #
##############################################################################################################

set.seed(9)


#######################################################################
#  Optimized Parameters of Random Forest Algorithm                    #
#   * Sequentially increment "mtry" from 1 to 10 with increment by 1  #     
#   * Set Number of Tree ("ntree") equal to 500                       #
#   * Set Number of "Cross Validation" equal to 5                     #
#######################################################################

train_ccpc_1203 <- caret::train(train_set_x, train_set_y , method="rf",
                                tuneGrid=data.frame(mtry=seq(1:10)),
                                ntree=500,
                                trControl=trainControl(method="cv", number=5),
                                importance= TRUE)


# Process and Generate Vars Importance of train_ccpc_rf
importance_rf_1203 <- varImp(train_ccpc_1203)
print(importance_rf_1203)


# Vars Importance of benchmark set
ccpc_importance_1203 <- rownames(importance_rf_1203$importance)[order(-importance_rf_1203$importance[,1])]
ccpc_importance_1203[1:4]


# Display Top 4 Most Relevant Molecular Descriptors/Features as Benchmark.
paste(ccpc_importance_1203[1],",",ccpc_importance_1203[2],"and",ccpc_importance_1203[3],
      "and",ccpc_importance_1203[4],
      "are selected as Benchmark for Features/Molecular Descriptors Prioritization")

# max Vars Importance
most_importance_1203 <-rownames(importance_rf_1203$importance)[which.max(importance_rf_1203$importance[,1])]
#most_importance_1203


# Compute Overall Accuracy for Molecular Descriptors/Feature Prioritization Method 
predict_ccpc_1203 <- predict(train_ccpc_1203, test_set_x)
accuracy_ccpc_1203 <- confusionMatrix(table(predict_ccpc_1203, test_set_y),mode="everything")
#accuracy_ccpc_1203
accuracy_ccpc_1203$overall["Accuracy"]

#Display Overall Accuracy
paste(round(accuracy_ccpc_1203$overall["Accuracy"],4),
      " will be adpoted as Benchmark Accuracy for Features/Molecular Descriptors Prioritization")


#########################################################
#                                                       #
#     Molecular Descriptors/Features Prioritization     #
#                                                       #
#########################################################

#  Formulating of MD1 and MD2 set are as follows:
#
#    * "MD1" set : Select Top pairwise Molecular Descriptors/Features generated from Vars Importance of Benchmark set
#                  using Random Forest Algorithm to formulate appropriate Molecular Descriptors "MD1" set.
#                  (Refer to Section 3.2.4 for Vars Importance of "MD1" set)
# 
#    * "MD2" set : Select Top ranked Molecular Descriptors/Features generated from Vars Importance of
#                  "1203 Molecular Descriptors" set using Random Forest Algorithm to formula appropriate
#                   Molecular Descriptors "MD2" set. (Refer to Section 3.3.5 for Vars Importance of "MD2" set)
#
#    * "RFE" set : Perform RFE function using Caret Package in R Environment. It recursively removes the least
#                  important Features until specified Number of Features is reached.
#
#    * "MDU" set : The MD1", "MD2" and "RFE" set are used to formulate different combination of Molecular Descriptors
#                  
#  Formulation of the "MDU" set will be conducted in the RStudio environment as follows:
# 
#    * Utilize a trial and error strategy as the methodology, offering flexibility in handling the testing procedure.
#
#    * The accuracy of "MDU" sets will be tested using Machine Learning Algorithms such as KNN, GBM, Random Forest, and
#      XGBoost.
#
#    * The Accuracy of various "MDU" set combinations will be compared with the Benchmark Accuracy (0.6571).
#      Sets exceeding the Benchmark will be selected and added to the "Candidate MDU" sets. Variances in Accuracy 
#      indicate interactions between Molecular Descriptors, which can explain Variability among them.
#
#    * Prioritized Molecular Descriptors (PMD) from the "Candidate MDU" sets with the highest Accuracy will be selected
#      as "PMD" Set for each Machine Learning Model.
#
#    * These Prioritized Molecular Descriptors for each Machine Learning Model will be applied for Final Model
#      training.
#
#  For different Model (e.g. KNN), we will adopt different Method/Technique such as Redundant Feature Elimination
#    for Molecular Descriptors/Feature Selection in order to achieve better Performance hence Enhance Robustness of
#    the Models.  Recursive Feature Elimination (RFE) is a Feature Selection Method in Machine Learning.
#    Perform RFE using the Caret package in R. Recursive Feature Elimination recursively removes the least
#    important Features until specified number of Features is reached.
#
#  For each ML Algorithm (KNN, GBM, Random Forest, XGBoost), Prioritized Molecular Descriptors/Features Selection
#  will be employed for Training the Model.


#
# Summary of Molecular Descriptors (MD1 Set , MD2 Set and RFE Set)
#

all_MD <- data.frame(SET="MD1 Set", Descriptors="MDEC_23, GATS8s")

all_MD <- bind_rows(all_MD,
                    data.frame(SET="MD2 Set",
                               Descriptors="MATS1e, SaasC, apol, GATS8i"))

all_MD <- bind_rows(all_MD,
                    data.frame(SET="RFE Set",
                               Descriptors="ZMIC1, ATSC8i, ATSC7p, MATS2s"))

all_MD %>% knitr::kable(caption="MD1 ,MD2 and RFE Set - Molecular Descriptors Prioritization")


#
# Formation of MDU Set: 1 or 2 Descriptors from MD1 Set + Any Combination of MD2 Set OR None of MD2 Set +
#                       Any Combination of RFE Set of OR None of RFE Set
# (e.g., MDU Set contains Molecular Descriptors MDEC_23, ZMIC1, ATSC8i )
#

#
# Summary of Prioritized Molecular Descriptors of each Model (PMD Set)
#

Model_MD <- data.frame(Model="Random Forest", Prioritized_Molecular_Descriptors="MDEC_23, GATS8s, MATS1e, SaasC")

Model_MD <- bind_rows(Model_MD,
                      data.frame(Model="K-Nearest Neighbor(KNN)",
                                 Prioritized_Molecular_Descriptors="MDEC_23, ZMIC1, ATSC8i"))

Model_MD <- bind_rows(Model_MD,
                      data.frame(Model="Gradient Boosted Machine (GBM)",
                                 Prioritized_Molecular_Descriptors="MDEC_23, GATS8s, MATS1e, SaasC"))

Model_MD <- bind_rows(Model_MD,
                      data.frame(Model="Extreme Gradient Boosting (XGBoost)",
                                 Prioritized_Molecular_Descriptors="MDEC_23, GATS8s, MATS1e, SaasC, apol, GATS8i"))

Model_MD %>% knitr::kable(caption="PMD Set - Molecular Descriptors Prioritization")


#######################################################################################
#                                                                                     #
# Apply Prioritized Molecular Descriptors on K-Nearest Neighbors Model (KNN) Training #
#                                                                                     #
#######################################################################################

# Generate df_ccpc_data dataset comprising Prioritized Molecular Descriptors of KNN Model
df_ccpc_data <- ccpc_data  %>% select(MDEC_23, ZMIC1, ATSC8i, -Class)


# Compute Correlation coefficient of df_ccpc_data
cor_ccpc_data <- cor(na.omit(df_ccpc_data[, unlist(lapply(df_ccpc_data, is.numeric))]))


######################################################################################
# Heat Map -  Prioritized Molecular Descriptors Correlation coefficient of KNN Model #
######################################################################################
df_cor_ccpc_data <- melt(cor_ccpc_data)
colnames(df_cor_ccpc_data) <- c("molecular_descriptors","toxicity","value")

ggplot(df_cor_ccpc_data , aes(x=molecular_descriptors,y=toxicity, fill=value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#006EBB", high = "#D883B7" , mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name="Correlation") +
  #theme(axis.text.x=element_text(angle=90)) +
  labs(x="", title="Prioritized Molecular Descriptors of Core Clock Protein Cryptochrome (CRY1) - KNN")


# Generate molecular_descriptors Matrix and toxicity Vector
molecular_descriptors <- df_ccpc_data %>% as.matrix()
toxicity <- factor(ccpc_data$Class)

# Generate "core_clock_protein_cry" set for Normalization
core_clock_protein_cry <- list(molecular_descriptors=molecular_descriptors, 
                               toxicity=toxicity )

# Display column names of df_ccpc_data - KNN
#colnames(df_ccpc_data)

# Display last 6 rows of df_ccpc_data- KNN
tail(df_ccpc_data)


######################################################################
#  Scaling and Normalizing data for K-Nearest Neighbors Model (KNN)  #
######################################################################

# Compute Mean of each Column
col_means <- colMeans(core_clock_protein_cry$molecular_descriptors, na.rm=TRUE)


# Compute Standard Deviations of each Column
col_sds <- colSds(core_clock_protein_cry$molecular_descriptors, na.rm=TRUE)


# Subtract Column Means
molecular_descriptors_centered <- sweep(core_clock_protein_cry$molecular_descriptors,
                                        2, col_means, FUN="-")

# Divide by Column Standard Deviations
molecular_descriptors_scaled <- sweep(molecular_descriptors_centered,
                                      2 , col_sds, FUN="/")

# Display first 6 rows of molecular_descriptors_scaled -KNN
head(molecular_descriptors_scaled)


##############################################################################
#  Generate 80% training set and 20% test set for K-Nearest Neighbors Model  #
##############################################################################
set.seed(1)

test_index <- createDataPartition(core_clock_protein_cry$toxicity, times = 1, p = 0.2, list = FALSE)

test_set_x   <-  molecular_descriptors_scaled[test_index,]
test_set_y   <-  core_clock_protein_cry$toxicity[test_index]

train_set_x  <- molecular_descriptors_scaled[-test_index,]
train_set_y  <- core_clock_protein_cry$toxicity[-test_index]


###############################################
#                                             #
#  K-Nearest Neighbors Model (KNN) Training   #
#                                             #  
###############################################

set.seed(7)


#####################################################
# Optimize Parameters of K-Nearest Neighbors Model  #
#  * Sequentially increment  "k" from 3 to 31       #     
#  * Set Number of "cross validation" equal to 10   #
#####################################################

control <- trainControl(method="cv", number=10)

train_ccpc_knn <- caret::train(train_set_x, train_set_y ,
                               method = "knn",
                               tuneGrid = data.frame(k = seq(3, 31, 1)),
                               trControl = control)


# ggplot number of Neighbors vs Accuracy
ggplot(train_ccpc_knn, highlight = TRUE) + 
  ggtitle("Core Clock Protein Cryptochrome (CRY1) Toxicity: K-Nearest Neighbors model (KNN)")

#Display the Best "k" Number of Neighbors for KNN Model 
train_ccpc_knn$bestTune

#Display Final Model of KNN
train_ccpc_knn$finalModel

# Makes Prediction on test set and Compute Overall Accuracy of KNN Model
predict_ccpc_knn <- predict(train_ccpc_knn, test_set_x)
accuracy_ccpc_knn <- confusionMatrix(table(predict_ccpc_knn, test_set_y),mode="everything")
#accuracy_ccpc_knn
accuracy_ccpc_knn$byClass["F1"]
accuracy_ccpc_knn$overall["Accuracy"]
#sum(predict_ccpc_knn == test_set_y) / length(test_set_y)


##############################################################################################################
#                                                                                                            #
# Apply Prioritized Molecular Descriptors on Random Forest and Gradient Boosted Machine Model (GBM) Training #
#                                                                                                            #
##############################################################################################################

# Generate df_ccpc_data dataset comprising Prioritized Molecular Descriptors of Random Forest and GBM Model
df_ccpc_data <- ccpc_data  %>% select(MDEC_23, GATS8s, MATS1e, SaasC, -Class)


# Compute Correlation coefficient of df_ccpc_data
cor_ccpc_data <- cor(na.omit(df_ccpc_data[, unlist(lapply(df_ccpc_data, is.numeric))]))


######################################################################################################
# Heat Map -  Prioritized Molecular Descriptors Correlation coefficient of Random Forest & GBM Model #
######################################################################################################
df_cor_ccpc_data <- melt(cor_ccpc_data)
colnames(df_cor_ccpc_data) <- c("molecular_descriptors","toxicity","value")

ggplot(df_cor_ccpc_data , aes(x=molecular_descriptors,y=toxicity, fill=value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#006EBB", high = "#D883B7" , mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name="Correlation") +
  #theme(axis.text.x=element_text(angle=90)) +
  labs(x="", title="Prioritized Features of Core Clock Protein Cryptochrome (CRY1) - Random Forest & GBM")


# Generate molecular_descriptors Matrix and toxicity Vector
molecular_descriptors <- df_ccpc_data %>% as.matrix()
toxicity <- factor(ccpc_data$Class)

# Generate core_clock_protein_cry for Normalization
core_clock_protein_cry <- list(molecular_descriptors=molecular_descriptors, 
                               toxicity=toxicity)

# Display column names of df_ccpc_data
#colnames(df_ccpc_data)

# Display last 6 rows of df_ccpc_data
tail(df_ccpc_data)


##################################################################
#  Scaling and Normalizing data for Random Forest and GBM Model  #
##################################################################

# Compute Mean of each Column
col_means <- colMeans(core_clock_protein_cry$molecular_descriptors, na.rm=TRUE)

# Compute Standard Deviations of each Column
col_sds <- colSds(core_clock_protein_cry$molecular_descriptors, na.rm=TRUE)

# Subtract Column Means
molecular_descriptors_centered <- sweep(core_clock_protein_cry$molecular_descriptors,
                                        2, col_means, FUN="-")

# Divide by Column Standard Deviations
molecular_descriptors_scaled <- sweep(molecular_descriptors_centered,
                                      2 , col_sds, FUN="/")

# Display first 6 rows of molecular_descriptors_scaled
head(molecular_descriptors_scaled)


################################################################################
#  Generate 80% training set and 20% test set for Random Forest and GBM  Model #
################################################################################
set.seed(1)

test_index <- createDataPartition(core_clock_protein_cry$toxicity, times = 1, p = 0.2, list = FALSE)

test_set_x   <-  molecular_descriptors_scaled[test_index,]
test_set_y   <-  core_clock_protein_cry$toxicity[test_index]

train_set_x  <- molecular_descriptors_scaled[-test_index,]
train_set_y  <- core_clock_protein_cry$toxicity[-test_index]


######################################
#                                    #
#    Random Forest Model Training    #
#                                    #
######################################


######################################################################
# Optimize Parameters of Random Forest Model                         #
#                                                                    #
#  * Sequentially increment "mtry" from 1 to 10 with increment by 1  #
#  * Set Number of Tree ("ntree") equal to 500                                          #
#  * Set Number of "cross validation" equal to 5                     #        
######################################################################

set.seed(9)

train_ccpc_rf <- caret::train(train_set_x, train_set_y , method="rf",
                              #tuneGrid=grid,
                              tuneGrid=data.frame(mtry=seq(1:10)),
                              ntree=500,
                              #ntree=30,
                              trControl=trainControl(method="cv", number=5),
                              importance= TRUE)


# Print Vars Importance of train_ccpc_rf
importance_rf <- varImp(train_ccpc_rf)
#print(importance_rf)

# print Most Vars Importance of train_ccpc_rf
most_importance_rf <-rownames(importance_rf$importance)[which.max(importance_rf$importance[,1])]
#most_importance_rf

# Plot "mtry" value vs Accuracy
ggplot(train_ccpc_rf, highlight = TRUE) + 
  ggtitle("Core Clock Protein Cryptochrome (CRY1) Toxicity: Random Forest Model")

#Display the best tine of Random Forest Model
train_ccpc_rf$bestTune


# Makes Prediction on test set and Compute Overall Accuracy of Random Forest Model
predict_ccpc_rf <- predict(train_ccpc_rf, test_set_x)
accuracy_ccpc_rf <- confusionMatrix(table(predict_ccpc_rf, test_set_y),mode="everything")
#accuracy_ccpc_rf
accuracy_ccpc_rf$byClass["F1"]
accuracy_ccpc_rf$overall["Accuracy"]

# Display accuracy 
#sum(predict_ccpc_rf == test_set_y)/length(test_set_y)

# Find and display best mtry value
#train_ccpc_rf$results$mtry[which.max(train_ccpc_rf$results$Accuracy)]


#######################################################
#                                                     #
#    Gradient Boosted Machine Model (GBM) Training    #  
#                                                     #
#######################################################


#######################################################################################################
# Optimize Parameters of Gradient Boosted Machine (GBM) Model                                         #
#                                                                                                     #
#  * Sequentially increment "n.trees" from 100 to 1000 with increment by 50                           #
#  * Set Number of splits ("interaction.depth") in each tree equal 1,3 and 5                          #
#  * Set Learning Rate ("shrinkage") equal to 0.01, 0.1 and 0.3                                       #   
#  * Set n.minobsinnode equal to 10 and 20                                                            #
#    (Minimum number of observations that must exist in a node for it to be considered for splitting) #
#  * Set Number of "cross validation" equal to 5                                                      #        
#######################################################################################################

set.seed(9)

control <- trainControl(method="cv",number=5)

grid <- expand.grid(n.trees = seq(100, 1000, by = 50),
                    interaction.depth=c(1,3,5),
                    shrinkage=c(0.01,0.1,0.3),
                    n.minobsinnode=c(10,20))

train_ccpc_GBM <- train(train_set_x,train_set_y, method="gbm",
                        trControl=control,
                        tuneGrid=grid,
                        verbose=FALSE)


# results contain Hyper-Parameters of GBM Model
results <- train_ccpc_GBM$results

# Convert the results to a long format for ggplot -GBM
results_long <- results %>% pivot_longer(cols = -c(Accuracy, Kappa, KappaSD),
                                         names_to = "Parameter",
                                         values_to = "Value")

# ggplot Hyper-Parameters of GBM Model
ggplot(results_long, aes(x = Value, y = Accuracy, color = Parameter)) +
  geom_point() +
  facet_wrap(~ Parameter, scales = "free_x") +
  labs(title = "Hyperparameter Tuning Results - GBM", x = "Parameter Value", y = "Accuracy")


# ggplot Accuracy vs Boosting Iterations - GBM Model
ggplot(train_ccpc_GBM, highlight = TRUE) + 
  ggtitle("Core Clock Protein Cryptochrome (CRY1) Toxicity: Gradient Boosting Machines(GBM)")

# Display the Best GBM Model
train_ccpc_GBM$bestTune


# Plot Boosting Iterations vs Accuracy
plot(train_ccpc_GBM)


# Makes Prediction on test set and Compute Overall Accuracy of GBM Model
predict_ccpc_GBM <- predict(train_ccpc_GBM, test_set_x)
accuracy_ccpc_GBM <- confusionMatrix(table(predict_ccpc_GBM, test_set_y),mode="everything")
#accuracy_ccpc_GBM 
accuracy_ccpc_GBM$byClass["F1"]
accuracy_ccpc_GBM$overall["Accuracy"]
#sum(predict_ccpc_GBM == test_set_y) / length(test_set_y)


#################################################################################################
#                                                                                               #
# Apply Prioritized Molecular Descriptors on Extreme Gradient Boosting Model (XGBoost) Training #
#                                                                                               #
#################################################################################################

# Generate df_ccpc_data dataset comprising Prioritized Molecular Descriptors of XGBoost Model
df_ccpc_data <- ccpc_data  %>% select(MDEC_23, GATS8s, MATS1e, SaasC, apol, GATS8i, -Class)

# Compute Correlation coefficient of df_ccpc_data
cor_ccpc_data <- cor(na.omit(df_ccpc_data[, unlist(lapply(df_ccpc_data, is.numeric))]))


##########################################################################################
# Heat Map -  Prioritized Molecular Descriptors Correlation coefficient of XGBoost Model #
##########################################################################################
df_cor_ccpc_data <- melt(cor_ccpc_data)
colnames(df_cor_ccpc_data) <- c("molecular_descriptors","toxicity","value")

ggplot(df_cor_ccpc_data , aes(x=molecular_descriptors,y=toxicity, fill=value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#006EBB", high = "#D883B7" , mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name="Correlation") +
  #theme(axis.text.x=element_text(angle=90)) +
  labs(x="", title="Prioritized Molecular Descriptors of Core Clock Protein Cryptochrome (CRY1) - XGBoost")


# Generate molecular_descriptors Matrix and toxicity Vector
molecular_descriptors <- df_ccpc_data %>% as.matrix()
toxicity <- factor(ccpc_data$Class)

# Generate core_clock_protein_cry for Normalization
core_clock_protein_cry <- list(molecular_descriptors=molecular_descriptors, 
                               toxicity=toxicity )

# Display column names of df_ccpc_data - XGBoost
#colnames(df_ccpc_data)

# Display last 6 rows of df_ccpc_data - XGBoost
tail(df_ccpc_data)


################################################################################
#  Scaling and Normalizing data for Extreme Gradient Boosting Model (XGBoost)  #
################################################################################

# Compute Mean of each Column
col_means <- colMeans(core_clock_protein_cry$molecular_descriptors, na.rm=TRUE)

# Compute Standard Deviations of each Column
col_sds <- colSds(core_clock_protein_cry$molecular_descriptors, na.rm=TRUE)

# Subtract Column Means
molecular_descriptors_centered <- sweep(core_clock_protein_cry$molecular_descriptors,
                                        2, col_means, FUN="-")

# Divide by Column Standard Deviations
molecular_descriptors_scaled <- sweep(molecular_descriptors_centered,
                                      2 , col_sds, FUN="/")

# Display first 6 rows of molecular_descriptors_scaled
head(molecular_descriptors_scaled)


###############################################################################################
#  Generate 80% training set and 20% test set for  Extreme Gradient Boosting Model (XGBoost)  #
###############################################################################################
set.seed(1)

test_index <- createDataPartition(core_clock_protein_cry$toxicity, times = 1, p = 0.2, list = FALSE)

test_set_x   <-  molecular_descriptors_scaled[test_index,]
test_set_y   <-  core_clock_protein_cry$toxicity[test_index]

train_set_x  <- molecular_descriptors_scaled[-test_index,]
train_set_y  <- core_clock_protein_cry$toxicity[-test_index]


#############################################################
#                                                           #
#     Extreme Gradient Boosting Model (XGBoost) Training    #
#                                                           #
#############################################################


#####################################################################################
# Optimize Parameters of Extreme Gradient Boosting (XGBoost) Model                  #
#                                                                                   #
#       * Set Number of boosting iterations ("nrounds") equal to 100                #
#       * Set Maximum depth of each tree ("max_depth") equal to 8                   #
#       * Set Learning rate ("eta") equal to 0.01, 0.1 and 0.3                      #
#       * Set "gamma" equal to 0, 0.2 and 0.4                                       #
#         (Minimum loss reduction required to split a node)                         #
#       * Set "colsample_bytree" equal to 0.5, 0.7 and 1                            #
#         (Fraction of Molecular Descriptors/Features used for training each tree)  #
#       * Set "min_child_weight" equal to 1, 3, and 5                               #
#         (Minimum sum of instance weight needed in a child )                       #
#       * Set "subsample" equal to 0.5, 0.7 and 1                                   #
#         (Fraction of samples used for training each Molecular Descriptors tree)   #
#       * Set Number of "Cross validation" equal to 5                               #
#####################################################################################

set.seed(123)

control <- trainControl(method="cv",number=5)

grid <- expand.grid( nrounds=100,
                     max_depth=8,
                     eta=c(0.01, 0.1, 0.3),
                     gamma=c(0, 0.2, 0.4),
                     colsample_bytree=c(0.5, 0.7, 1),
                     min_child_weight=c(1, 3, 5),
                     subsample=c(0.5, 0.7, 1))

train_ccpc_XGB <- train(train_set_x,train_set_y,
                        method="xgbTree",
                        trControl=control,
                        tuneGrid=grid)


# results contain Hyper-Parameters of XGBoost Model
results <- train_ccpc_XGB$results

# Convert the results to a long format for ggplot -XGBoost
results_long <- results %>% pivot_longer(cols = -c(Accuracy, Kappa, KappaSD, nrounds, max_depth),
                                         names_to = "Parameter",
                                         values_to = "Value")

# ggplot Hyper-Parameters of XGBoost Model
ggplot(results_long, aes(x = Value, y = Accuracy, color = Parameter)) +
  geom_point() +
  facet_wrap(~ Parameter, scales = "free_x") +
  labs(title = "Hyperparameter Tuning Results - XGBoost", x = "Parameter Value", y = "Accuracy")

# Display the Best Model of XGBoost
train_ccpc_XGB$bestTune

# Display the Final Model of XGBoost
train_ccpc_XGB$finalModel

# print Feature Importance of XGBoost Model
importance_matrix <- xgb.importance(model=train_ccpc_XGB$finalModel)
print(importance_matrix)

# plot Molecular Descriptors of Vars Importance - XGBoost
xgb.plot.importance(importance_matrix=importance_matrix)

# plot Hyper-Parameters - XGBoost
plot(train_ccpc_XGB)


# Makes Prediction on test set and Compute Overall Accuracy of XGBoost Model
predict_ccpc_XGB <- predict(train_ccpc_XGB, test_set_x)
accuracy_ccpc_XGB <- confusionMatrix(table(predict_ccpc_XGB, test_set_y),mode="everything")
#accuracy_ccpc_XGB
accuracy_ccpc_XGB$byClass["F1"]
accuracy_ccpc_XGB$overall["Accuracy"]
#sum(predict_ccpc_XGB == test_set_y) / length(test_set_y)


################################################################################################
#                                                                                              #
#   Generate an Ensemble Model of Core Clock Protein Cryptochrome (CRY1) Toxicity Prediction   #
#                                                                                              #
################################################################################################

set.seed(1)

# Create Dataframe contains predictions of GBM, KNN, Random Forest and XGBoost Model
combined_predict_ccpc <- data.frame( predict_ccpc_GBM, predict_ccpc_knn, predict_ccpc_rf,
                                     predict_ccpc_XGB)


# Generate a Majority Prediction of Toxicity Class using "combined_predict_ccpc"- My Version
ensemble_ccpc <- apply(combined_predict_ccpc, 1, function(row)
{ ifelse(mean(row =="Toxic") > 0.5, "Toxic","NonToxic") })

# Display Result of Majority Prediction of Toxicity Class
#ensemble_ccpc


# Compute Accuracy of an Ensemble Model- My Version
accuracy_ccpc_ensemble <- confusionMatrix(factor(ensemble_ccpc), test_set_y, mode="everything")
accuracy_ccpc_ensemble

# Display Overall Accuracy of Ensemble Model from confusionMatrix
accuracy_ccpc_ensemble$overall["Accuracy"]

# Create a table of Accuracy encompassing 4 Models + Ensemble Model (My version- Majority Prediction)
accuracies_table_ccpc <- data.frame(Model    = c( "K-Nearest Neighbors(KNN)",
                                                  "Gradient Boosting Machines(GBM)",
                                                  "Extreme Gradient Boost(XGBoost)",
                                                  "Random Forest",
                                                  "Ensemble"),
                                    Accuracy = c( accuracy_ccpc_knn$overall["Accuracy"],
                                                  accuracy_ccpc_GBM$overall["Accuracy"],
                                                  accuracy_ccpc_XGB$overall["Accuracy"],
                                                  accuracy_ccpc_rf$overall["Accuracy"],
                                                  accuracy_ccpc_ensemble$overall["Accuracy"]))


# Display Summary table of Accuracy for all Models using knitr function of Markdown (My version- Majority Prediction)
knitr::kable(accuracies_table_ccpc, caption="Accuracy of all Model")
best_ensemble_ccpc <- accuracies_table_ccpc$Model[which.max(accuracies_table_ccpc$Accuracy)]
best_ensemble_ccpc


##################################################
#  Ensemble Model - Plot Precision-Recall Curve  #
##################################################

# Convert factors to binary (0, 1) for ROC analysis
binary_labels <- ifelse(test_set_y == "NonToxic", 1, 0) 
binary_preds <- ifelse(ensemble_ccpc == "NonToxic", 1, 0) 

# Create data for precision-recall plot 
pr_data <- pr.curve(scores.class0 = binary_preds,
                    weights.class0 = binary_labels == 1,
                    curve = TRUE) 

# Extract precision-recall data 
pr_df <- data.frame(Recall = pr_data$curve[, 1], Precision = pr_data$curve[, 2]) 

# Plot Precision-Recall curve 
pr_plot <- ggplot(pr_df, aes(x = Recall, y = Precision)) + 
  geom_line(color = "green", size = 1) + 
  ggtitle(" Core Clock Protein Cryptochrome (CRY1) Toxicity Prediction - Precision-Recall Curve") + 
  xlab("Recall") +
  ylab("Precision") +
  theme_minimal() 
print(pr_plot)



#################################################################################
##     Core Clock Protein Cryptochrome (CRY1) Toxicity Prediction - Ending     ##
#################################################################################

