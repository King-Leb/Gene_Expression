# Install GEOquery from Bioconductor
BiocManager::install("GEOquery")
library(dplyr)
library(tidyverse)
library(GEOquery)
library(ggplot2)

# List files in the directory
list.files("/kaggle/input/gene-expression-omnibus-geo-dataset-gse68086")

# Load the data from CSV file
bio_data <- read.csv("/kaggle/input/gene-expression-omnibus-geo-dataset-gse68086/GSE68086_TEP_data_matrix.csv")

# get Metadata
gmd <- getGEO(GEO = "GSE68086", GSEMatrix = TRUE)

# get phenodata from the metadata
metaD <- pData(phenoData(gmd[[1]]))
head(metaD)

# select and format the structure of the metadata
metaD.modify <- metaD %>%
  select(1,8,10,46,47,48) %>%
  rename(Sample_ID = title) %>%
  rename(Sample_Description = source_name_ch1)%>%
  rename(Tissue_Type = characteristics_ch1)%>%
  rename(Cancer_Type = 'cancer type:ch1')%>%
  rename(Cell_Type = 'cell type:ch1')%>%
  rename(Mutational_Subclass = 'mutational subclass:ch1')%>%
  select(-Tissue_Type)
head(metaD.modify)

#reshaping data from wide to long
bio_data.long <- bio_data %>%
  rename(gene_ID = X)%>%
  gather(key = 'samples', value = 'FPKM', -gene_ID)

head(bio_data.long)

# time to combine the two data frames, that is, bio_data.long + metaD.modify

# first we need to remove the prefix X attached to the samples column of the bio_data.long dataframe

bio_data.long$samples <- gsub("^X", "", bio_data.long$samples)

# we replace periods with hyphens in the 'samples' column
bio_data.long$samples <- gsub("\\.", "-", bio_data.long$samples)

# we apply the left join now
combined <- bio_data.long%>%
  left_join(., metaD.modify, by = c('samples' = 'Sample_Description'))

head(combined)

# Boxplot comparing FPKM across samples
ggplot(combined, aes(x = samples, y = FPKM)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "FPKM across Samples", x = "Samples", y = "FPKM")

# Boxplot comparing FPKM between Cancer Types
ggplot(combined, aes(x = Cancer_Type, y = FPKM, fill = Cancer_Type)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Gene Expression by Cancer Type", x = "Cancer Type", y = "FPKM")

# Bar Plot for average FPKM by Mutational Subclass
df_avg <- combined %>%
  group_by(Mutational_Subclass, gene_ID) %>%
  summarise(avg_FPKM = mean(FPKM))

ggplot(df_avg, aes(x = gene_ID, y = avg_FPKM, fill = Mutational_Subclass)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Average FPKM by Mutational Subclass", x = "Gene ID", y = "Average FPKM")

# Randomly sample 5000 data points
set.seed(123)  # for reproducibility
sampled_data <- combined[sample(1:nrow(combined), 5000), ]

# Perform Shapiro-Wilk test on the sampled data
shapiro_test <- shapiro.test(sampled_data$FPKM)
shapiro_test

hist(combined$FPKM, breaks = 10, main = "Histogram of FPKM values", xlab = "FPKM")

qqnorm(combined$FPKM)
qqline(combined$FPKM, col = "red")

# Perform Kruskal-Wallis test
kruskal_test <- kruskal.test(FPKM ~ Cancer_Type, data = combined)

# View the results
kruskal_test

# Create a violin plot to show the density of the data
ggplot(combined, aes(x = Cancer_Type, y = FPKM)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +  # Add boxplot inside violin
  labs(title = "Violin Plot of FPKM by Cancer Type",
       x = "Cancer Type",
       y = "FPKM") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(dunn.test)
library(caret)

# Perform Dunn's test
dunn_test <- dunn.test(combined$FPKM, combined$Cancer_Type, method="bonferroni")
print(dunn_test)

# Fit a linear model
linear_model <- lm(FPKM ~ Cancer_Type, data = combined)

summary(linear_model)


# ML-----------------------------------------
set.seed(123)
trainIndex <- createDataPartition(combined$Cancer_Type, p = .7, 
                                  list = FALSE, 
                                  times = 1)
train_set <- combined[trainIndex, ]
test_set <- combined[-trainIndex, ]

library(nnet)
multinomial_model <- multinom(Cancer_Type ~ FPKM, data = train_set)

# For multinomial logistic regression
predicted_classes <- predict(multinomial_model, newdata = test_set)

confusion_matrix <- table(test_set$Cancer_Type, predicted_classes)
print(confusion_matrix)

precision <- diag(confusion_matrix) / rowSums(confusion_matrix)
recall <- diag(confusion_matrix) / colSums(confusion_matrix)
f1_score <- 2 * (precision * recall) / (precision + recall)

# Print precision, recall, and F1-score
print(paste("Precision: ", round(precision, 4)))
print(paste("Recall: ", round(recall, 4)))
print(paste("F1 Score: ", round(f1_score, 4)))

# Visualizing the confusion matrix using a heatmap
library(ggplot2)
library(reshape2)

# Melt the confusion matrix
conf_matrix_melted <- as.data.frame(as.table(confusion_matrix))

# Visualize the confusion matrix with a heatmap
ggplot(conf_matrix_melted, aes(x = Actual, y = Predicted, fill = Count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "Actual", y = "Predicted", fill = "Count") +
  ggtitle("Confusion Matrix Heatmap") +
  theme_minimal()


library(xgboost)
train_set$Cancer_Type <- as.numeric(as.factor(train_set$Cancer_Type)) - 1
numeric_features <- train_set[, sapply(train_set, is.numeric)]
feature_matrix <- as.matrix(numeric_features)

# Create an xgb.DMatrix
dtrain <- xgb.DMatrix(data = feature_matrix, label = train_set$Cancer_Type)

# Define parameters for XGBoost
param <- list(max_depth = 6, eta = 0.1, objective = "multi:softmax", num_class = length(unique(train_set$Cancer_Type)))

# Train the XGBoost model
xgb_model <- xgboost(params = param, data = dtrain, nrounds = 100)

# Prepare test set
test_set$Cancer_Type <- as.numeric(as.factor(test_set$Cancer_Type)) - 1
test_features <- as.matrix(test_set[, sapply(test_set, is.numeric)])

# Make predictions on test data
preds <- predict(xgb_model, test_features)

# Convert predictions to classes
predicted_classes <- as.factor(round(preds))

confusion_matrix <- confusionMatrix(predicted_classes, as.factor(test_set$Cancer_Type))
print(confusion_matrix)

train_set <- train_set[, c("FPKM", "Cancer_Type")]
train_set$Cancer_Type <- as.factor(train_set$Cancer_Type)

# Prepare the data for xgboost
dtrain <- xgb.DMatrix(data = as.matrix(train_set[, "FPKM", drop = FALSE]), label = as.numeric(train_set$Cancer_Type) - 1)


# Set parameters for the XGBoost model
params <- list(
  objective = "multi:softmax",  # Multiclass classification
  num_class = length(levels(train_set$Cancer_Type)),  # Number of classes
  eta = 0.1,                     # Learning rate
  max_depth = 6,                 # Maximum depth of trees
  eval_metric = "mlogloss"       # Evaluation metric
)

# Train the model with cross-validation
set.seed(123)  # For reproducibility
cv_model <- xgb.cv(
  params = params,
  data = dtrain,
  nrounds = 100,               # Number of boosting rounds
  nfold = 5,                   # Number of folds for cross-validation
  early_stopping_rounds = 10,  # Stop if no improvement
  verbose = 1
)

# Train final model on the full dataset
final_model <- xgboost(
  params = params,
  data = dtrain,
  nrounds = cv_model$best_iteration  # Use the best iteration found during CV
)

# Extract the trees from the model and print them
trees <- xgb.model.dt.tree(model = final_model)
print(trees)

test_X <- as.matrix(test_set[, "FPKM"])
test_y <- as.numeric(test_set$Cancer_Type) - 1

test_preds <- predict(final_model, newdata = xgb.DMatrix(test_X))

test_preds_factor <- as.factor(test_preds)
test_y_factor <- as.factor(test_y)

# Confusion matrix
confusion_matrix <- table(Predicted = test_preds_factor, Actual = test_y_factor)
print(confusion_matrix)

# Detailed confusion matrix and statistics
conf_matrix_stats <- confusionMatrix(test_preds_factor, test_y_factor)
print(conf_matrix_stats)




