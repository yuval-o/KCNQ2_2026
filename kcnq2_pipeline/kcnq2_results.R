# Packages
pkgs <- c("readxl","data.table","caret","gbm","pROC","PRROC","doParallel")
for (p in pkgs) {
  if (!require(p, character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

# Read Excel "kcnq2_tabletorun.xlsx" - EDIT PATH
file_path <- ""
df <- as.data.table(read_excel(file_path))

# Build Class from label (GOF=1, LOF=0)
stopifnot("label" %in% names(df))
df[, Class := ifelse(label == 1, "gof", "lof")]
df[, Class := factor(Class, levels = c("gof", "lof"))]  # gof = positive

cat("Class balance:\n")
print(table(df$Class)); cat("\n")

info_cols <- c("protid","USED_REF","alt","pos","refAA","altAA","label")
feature_cols <- setdiff(names(df), c(info_cols, "Class"))

non_numeric <- feature_cols[!sapply(df[, ..feature_cols], is.numeric)]
if (length(non_numeric) > 0) stop(paste("Non-numeric predictors:", paste(non_numeric, collapse=", ")))

# ensure no NA in predictors
na_in_feat <- colSums(is.na(df[, ..feature_cols]))
if (any(na_in_feat > 0)) {
  print(na_in_feat[na_in_feat > 0])
  stop("NA values exist in predictors.")
}

# modeling table
dat <- as.data.frame(df[, c(feature_cols, "Class"), with = FALSE])

cat("Rows:", nrow(dat), "\n")
cat("Predictors:", length(feature_cols), "\n\n")

fitControl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,  # caret will report ROC
  savePredictions = "final"
)

# Train GBM (with parallel)
cl <- makePSOCKcluster(5)
registerDoParallel(cl)

set.seed(42)
model <- train(
  Class ~ .,
  data = dat,
  method = "gbm",
  metric = "ROC",
  trControl = fitControl,
  verbose = FALSE
)

stopCluster(cl)

cat("=== caret model summary ===\n")
print(model)
cat("\nBest tuning:\n")
print(model$bestTune)
cat("\n")

# Compute metrics on UNIQUE variants
pred <- as.data.table(model$pred)

bt <- model$bestTune
for (nm in names(bt)) {
  pred <- pred[get(nm) == bt[[nm]]]
}

oof <- pred[, .(
  obs = obs[1],
  gof_prob = mean(gof)
), by = rowIndex]

oof[, obs := factor(obs, levels = c("gof","lof"))]

cat("N unique variants (should be 105):", nrow(oof), "\n\n")


# Metrics function (BA, Kappa, MCC, ROC, prAUC)

calc_metrics_oof <- function(obs, prob_gof, threshold = 0.5, positive = "gof") {
  neg <- setdiff(levels(obs), positive)[1]
  
  pred_class <- ifelse(prob_gof >= threshold, positive, neg)
  pred_class <- factor(pred_class, levels = levels(obs))
  
  cm <- caret::confusionMatrix(pred_class, obs, positive = positive)
  
  BA <- unname(cm$byClass["Balanced Accuracy"])
  Kappa <- unname(cm$overall["Kappa"])
  
  tab <- table(pred_class, obs)
  TP <- as.numeric(tab[positive, positive])
  TN <- as.numeric(tab[neg, neg])
  FP <- as.numeric(tab[positive, neg])
  FN <- as.numeric(tab[neg, positive])
  
  denom <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  MCC <- ifelse(denom == 0, NA, (TP*TN - FP*FN)/denom)
  
  roc_obj <- pROC::roc(response = obs, predictor = prob_gof,
                       levels = c("lof","gof"), direction = "<")
  ROC <- as.numeric(pROC::auc(roc_obj))
  
  scores_pos <- prob_gof[obs == positive]
  scores_neg <- prob_gof[obs != positive]
  pr <- PRROC::pr.curve(scores.class0 = scores_pos, scores.class1 = scores_neg, curve = FALSE)
  prAUC <- pr$auc.integral
  
  list(confusion = cm$table, BA = BA, Kappa = Kappa, MCC = MCC, ROC = ROC, prAUC = prAUC)
}

# Compute & print metrics at thresholds 0.3 / 0.5 / 0.7

m05 <- calc_metrics_oof(oof$obs, oof$gof_prob, threshold = 0.5)
m03 <- calc_metrics_oof(oof$obs, oof$gof_prob, threshold = 0.3)
m07 <- calc_metrics_oof(oof$obs, oof$gof_prob, threshold = 0.7)

cat("=== Metrics (UNIQUE OOF, threshold=0.5) ===\n")
print(m05$confusion)
cat(sprintf("BA=%.3f | Kappa=%.3f | MCC=%.3f | ROC=%.3f | prAUC=%.3f\n\n",
            m05$BA, m05$Kappa, m05$MCC, m05$ROC, m05$prAUC))

cat("=== Metrics (UNIQUE OOF, threshold=0.3) ===\n")
print(m03$confusion)
cat(sprintf("BA=%.3f | Kappa=%.3f | MCC=%.3f | ROC=%.3f | prAUC=%.3f\n\n",
            m03$BA, m03$Kappa, m03$MCC, m03$ROC, m03$prAUC))

cat("=== Metrics (UNIQUE OOF, threshold=0.7) ===\n")
print(m07$confusion)
cat(sprintf("BA=%.3f | Kappa=%.3f | MCC=%.3f | ROC=%.3f | prAUC=%.3f\n\n",
            m07$BA, m07$Kappa, m07$MCC, m07$ROC, m07$prAUC))


## Variable importance (gbm loaded -> no relative.influence error)
#cat("=== Top variable importance ===\n")
#imp <- varImp(model, scale = TRUE)
#imp_df <- imp$importance
#imp_df <- imp_df[order(-imp_df$Overall), , drop = FALSE]
#print(head(imp_df, 20))