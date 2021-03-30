# https://www.r-bloggers.com/2018/01/how-to-implement-random-forests-in-r/
# https://www.blopig.com/blog/2017/04/a-very-basic-introduction-to-random-forests-using-r/
library(randomForest)
library(ROCR)
library(ggplot2)
library(ggpubr)
library(pheatmap)

# Read in guide score data
dat_file <- "output_data/out_val/table_val_guideScores_whole.txt"
dat <- read.delim(dat_file, h = TRUE, as.is = TRUE)
dat <- dat[-which(is.na(dat$distance_tss_0_stop_1)),]

output_folder <- file.path(dirname(dat_file), "guidePosition")
if (!dir.exists(output_folder)) dir.create(output_folder)

# Define classifier for random forest model
#classifier <- "escapeNMD"
classifier <- "distance_tss_0_stop_1"
classifier_name <- "position"

# Prepare data for random forest classification
pred <- dat[,c("logFC", "Exon_Length", "vbc_score",
  "doench_score", "provean_score", "disorder_score",
  "guideEfficiency_pred", "predicted_in_frame", "predicted_oof",
  "sequence_score", "escapeNMD", "distance_tss_0_stop_1")]
pred <- na.omit(pred)

# Set classifier levels
quant1 <- quantile(pred[[classifier]], 0.25)
quant3 <- quantile(pred[[classifier]], 0.75)
pred[[classifier_name]] <- "middle"
pred[which(pred[[classifier]] <= quant1), classifier_name] <- "start"
pred[which(pred[[classifier]] >= quant3), classifier_name] <- "end"
pred[[classifier_name]] <- as.factor(pred[[classifier_name]])
pred <- pred[-which(colnames(pred) %in% classifier)]

# Simple correlation analysis
pred_corr <- cor(pred, use = "pairwise.complete.obs")
pdf(file.path(output_folder, "guide_feature_pcc.pdf"))
pheatmap(pred_corr, display_numbers = TRUE)
dev.off()

# Split into Train and Validation sets
# Training Set : Validation Set = 70 : 30 (random)
set.seed(42)
train <- sample(nrow(pred), 0.7 * nrow(pred), replace = FALSE)
TrainSet <- pred[train,]
ValidSet <- pred[-train,]

# Create a Random Forest model with default parameters
model1 <- randomForest(position ~ ., data = TrainSet, importance = TRUE)

# Validation set assessment #1: looking at confusion matrix
# Predicting on train set
predTrain <- predict(model1, TrainSet, type = "class")

# Checking classification accuracy
table(predTrain, TrainSet$position)

# Predicting on Validation set
predValid <- predict(model1, ValidSet, type = "class")

# Checking classification accuracy
mean(predValid == ValidSet$position)
table(predValid, ValidSet$position)

# To check important variables
imp <- importance(model1)
imp <- as.data.frame(imp)
imp <- cbind(feature = rownames(imp), imp)
imp <- imp[order(imp$`%IncMSE`),]
imp$feature <- factor(imp$feature, levels = imp$feature)
#varImpPlot(model1)

# Plot
p <- ggplot(imp, aes(x = `%IncMSE`, y = feature)) +
      geom_bar(stat = "identity") +
      labs(y = NULL, x = "Feature importance measure",
           title = "Features that best predict guide\ntarget position (random forest)") +
      theme_pubr()

# Draw out
plot_out <- file.path(output_folder, "rf_guide_feature_importance_position.pdf")
ggsave(plot_out, p,, width = 5, height = 5)

# Validation set assessment #2: ROC curves and AUC
# Calculate the probability of new observations belonging to each class
# prediction_for_roc_curve will be a matrix with dimensions data_set_size x number_of_classes
predROC <- predict(model1, ValidSet, type = "prob")

# Specify the different classes
classes <- levels(ValidSet$position)

# For each class
for (i in seq_along(classes)) {

 # Define which observations belong to class[i]
 true_values <- ifelse(ValidSet[,11] == classes[i],1,0)

 # Assess the performance of classifier for class[i]
 pred <- prediction(predROC[,i], true_values)
 perf <- performance(pred, "tpr", "fpr")

 if (i == 1) {
     plot(perf, main = "ROC Curve")
 } else {
     plot(perf, main = "ROC Curve", add = TRUE)
 }

 # Calculate the AUC and print it to screen
 auc.perf <- performance(pred, measure = "auc")
 print(auc.perf@y.values)
}
