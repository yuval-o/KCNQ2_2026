# Packages
pkgs <- c("data.table","caret","gbm","pROC","PRROC","doParallel")
for (p in pkgs) {
  if (!require(p, character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}
set.seed(1)

# Paths (EDIT THESE)
# supp_path = SupplementaryTable_S1_pathvariantsusedintraining_revision2.txt
# featuretable_path = featuretable4github_revision.txt
supp_path <- ""
featuretable_path <- ""


# Read labels (Heyne) - TXT
varall <- data.table::fread(
  supp_path,
  sep = "\t",
  data.table = FALSE
)

need_cols <- c("protid","gene","pos","altAA","prd_mech_revised","used_in_functional_prediction")
miss_lab <- setdiff(need_cols, colnames(varall))
if (length(miss_lab) > 0) stop(paste("Missing columns in SupplementaryTable:", paste(miss_lab, collapse=", ")))

# Filter like Heyne
varall <- varall[varall$used_in_functional_prediction == 1, ]
varall <- varall[varall$prd_mech_revised %in% c("gof","lof"), ]

# Deduplicate like Heyne: unique by (gene, pos, altAA)
varall <- varall[!duplicated(varall[, c("gene","altAA","pos")]), ]

# key for aggregation
varall$variant_key <- paste(varall$gene, varall$pos, varall$altAA, sep=":")

cat("\nAfter Heyne filtering:\n")
cat("n_rows:", nrow(varall), "\n")
print(table(varall$prd_mech_revised))


# Evaluation 3: Hybrid Classifier

mixed_genes <- c("SCN5A","SCN1A","SCN2A","CACNA1A","SCN4A")
single_lof_genes <- c("CACNA1S","CACNA1F")
single_gof_genes <- c("SCN8A","CACNA1C","SCN9A","CACNA1D","CACNA1E")

varall$gene_group <- ifelse(varall$gene %in% mixed_genes, "mixed",
                            ifelse(varall$gene %in% single_lof_genes, "single_lof",
                                   ifelse(varall$gene %in% single_gof_genes, "single_gof", "other")))

cat("\nGene-group counts (after filtering + dedup):\n")
print(table(varall$gene_group, useNA = "ifany"))

keep_groups <- c("mixed","single_lof","single_gof")
var_exp3 <- varall[varall$gene_group %in% keep_groups, ]

cat("\nEXP3 kept groups only:\n")
cat("n_rows:", nrow(var_exp3), "\n")
print(table(var_exp3$gene_group))
cat("\nLabel counts by group:\n")
print(table(var_exp3$prd_mech_revised, var_exp3$gene_group))


# Read features (Heyne) - BIG TSV

featuretable <- data.table::fread(
  featuretable_path,
  sep = "\t",
  data.table = FALSE,
  showProgress = TRUE
)

if (!("protid" %in% colnames(featuretable))) stop("featuretable must contain column 'protid'")

# Join by protid (like Heyne) on EXP3 subset
feat_all <- featuretable[match(var_exp3$protid, featuretable$protid), ]

# Attach labels + keys + gene info
feat_all$Class <- factor(var_exp3$prd_mech_revised, levels = c("gof","lof"))  # GOF positive (first level)
feat_all$variant_key <- var_exp3$variant_key
feat_all$gene <- var_exp3$gene
feat_all$gene_group <- var_exp3$gene_group


# OUR feature mapping

# 1) S4 = IIIS4 only
if (!("S4M3" %in% colnames(feat_all))) stop("Missing S4M3 in Heyne featuretable (needed for IIIS4).")
feat_all$S4 <- feat_all$S4M3

# 2) bury/medium/exposed continuous -> binary one-hot (argmax)
need_bme <- c("bury","medium","exposed")
if (!all(need_bme %in% colnames(feat_all))) stop("Missing bury/medium/exposed in Heyne featuretable.")
bme <- as.matrix(feat_all[, need_bme, drop=FALSE])
max_idx <- max.col(bme, ties.method = "first")
feat_all$bury    <- as.integer(max_idx == 1)
feat_all$medium  <- as.integer(max_idx == 2)
feat_all$exposed <- as.integer(max_idx == 3)

# 3) gofsotherKCNQgenes analogue
need_gofs <- c("scn","gofsotherSCNgenes","gofsotherCACgenes")
if (!all(need_gofs %in% colnames(feat_all))) stop("Missing scn/gofsotherSCNgenes/gofsotherCACgenes in Heyne featuretable.")
feat_all$gofsotherKCNQgenes <- ifelse(feat_all$scn == 1, feat_all$gofsotherSCNgenes, feat_all$gofsotherCACgenes)

# 4) kcnqcon analogue
need_con <- c("scn","scncon","caccon")
if (!all(need_con %in% colnames(feat_all))) stop("Missing scn/scncon/caccon in Heyne featuretable.")
feat_all$kcnqcon <- ifelse(feat_all$scn == 1, feat_all$scncon, feat_all$caccon)

# 5) 50/95 mapping into your names (50.kcnq / 95.kcnq)
has_50scn <- "50.scn" %in% colnames(feat_all)
has_95scn <- "95.scn" %in% colnames(feat_all)
has_50scncac <- "50.scncac" %in% colnames(feat_all)
has_95scncac <- "95.scncac" %in% colnames(feat_all)

if (has_50scn && has_50scncac) {
  feat_all[["50.kcnq"]] <- ifelse(feat_all$scn == 1, feat_all[["50.scn"]], feat_all[["50.scncac"]])
} else if (has_50scn) {
  feat_all[["50.kcnq"]] <- feat_all[["50.scn"]]
} else if (has_50scncac) {
  feat_all[["50.kcnq"]] <- feat_all[["50.scncac"]]
} else {
  stop("Missing both 50.scn and 50.scncac in Heyne featuretable.")
}

if (has_95scn && has_95scncac) {
  feat_all[["95.kcnq"]] <- ifelse(feat_all$scn == 1, feat_all[["95.scn"]], feat_all[["95.scncac"]])
} else if (has_95scn) {
  feat_all[["95.kcnq"]] <- feat_all[["95.scn"]]
} else if (has_95scncac) {
  feat_all[["95.kcnq"]] <- feat_all[["95.scncac"]]
} else {
  stop("Missing both 95.scn and 95.scncac in Heyne featuretable.")
}


# Define OUR feature list

our_features <- c(
  "S4",
  "H","G","I","E","B","T","S","L",
  "bury","medium","exposed",
  "gofsotherKCNQgenes",
  "grantham",
  "ref_hydrophobic-side-chain:aliphatic",
  "ref_negatively-charged",
  "hydrophobic-side-chain:aliphatic",
  "negatively-charged",
  "hydrophobicity_altAA","hydrophobicity_refAA",
  "hydrophobicity3_altAA","hydrophobicity3_refAA",
  "hydrophobicity4_altAA","hydrophobicity4_refAA",
  "densgof","densgof3aa",
  "kcnqcon","50.kcnq","95.kcnq"
)

missing <- setdiff(our_features, colnames(feat_all))
if (length(missing) > 0) stop(paste("Missing required features:", paste(missing, collapse=", ")))


# Build matrices

train_idx <- which(feat_all$gene_group == "mixed")

X_all <- feat_all[, our_features, drop=FALSE]
y_all <- feat_all$Class

X_train <- X_all[train_idx, , drop=FALSE]
y_train <- y_all[train_idx]

nzv <- caret::nearZeroVar(X_train)
if (length(nzv) > 0) {
  X_train <- X_train[, -nzv, drop=FALSE]
  X_all   <- X_all[, -nzv, drop=FALSE]
}

cat("\nDesign matrices:\n")
cat("TRAIN (MIXED): n_rows:", nrow(X_train), " | n_features:", ncol(X_train), "\n")
print(table(y_train))
cat("ALL kept groups: n_rows:", nrow(X_all), " | n_features:", ncol(X_all), "\n")
print(table(y_all))


# Train the model on "mixed" genes only

ctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

n_cores <- max(1, parallel::detectCores() - 1)
cl <- makePSOCKcluster(n_cores)
doParallel::registerDoParallel(cl)

gbm_grid <- expand.grid(
  n.trees = c(1000),
  interaction.depth = c(1,2,3),
  shrinkage = c(0.01),
  n.minobsinnode = c(5,10)
)

model <- caret::train(
  x = X_train, y = y_train,
  method = "gbm",
  trControl = ctrl,
  tuneGrid = gbm_grid,
  metric = "ROC",
  verbose = FALSE
)

stopCluster(cl)

cat("\nBest tune:\n")
print(model$bestTune)
cat("\nCV ROC (caret):\n")
print(max(model$results$ROC, na.rm = TRUE))

# OOF predictions (MIXED only)

pred <- model$pred
best <- model$bestTune

pred_best <- pred[
  pred$n.trees == best$n.trees &
    pred$interaction.depth == best$interaction.depth &
    pred$shrinkage == best$shrinkage &
    pred$n.minobsinnode == best$n.minobsinnode, ]

feat_train <- feat_all[train_idx, , drop=FALSE]

pred_best$variant_key <- feat_train$variant_key[pred_best$rowIndex]
pred_best$y_true <- feat_train$Class[pred_best$rowIndex]

if (!("gof" %in% colnames(pred_best))) stop("Missing 'gof' probability column in caret predictions.")
pred_best$p_gof <- pred_best$gof

pred_best <- data.table::as.data.table(pred_best)

agg <- pred_best[, .(
  y_true = as.character(stats::na.omit(y_true)[1]),
  p_gof = mean(p_gof, na.rm = TRUE)
), by = variant_key]

agg$y_true <- factor(agg$y_true, levels=c("gof","lof"))

cat("\n[MIXED] Unique variants after OOF aggregation:", nrow(agg), "\n")
print(table(agg$y_true))


# Metrics

mcc <- function(tp, tn, fp, fn) {
  tp <- as.numeric(tp); tn <- as.numeric(tn); fp <- as.numeric(fp); fn <- as.numeric(fn)
  num <- (tp * tn) - (fp * fn)
  den_sq <- (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
  if (is.na(den_sq) || den_sq <= 0) return(NA_real_)
  den <- sqrt(den_sq)
  if (den == 0) return(NA_real_)
  num / den
}

metrics_at_threshold <- function(y_true, p_gof, thr) {
  pred_lab <- ifelse(p_gof >= thr, "gof", "lof")
  y_true_chr <- as.character(y_true)
  
  tp <- sum(pred_lab == "gof" & y_true_chr == "gof")
  tn <- sum(pred_lab == "lof" & y_true_chr == "lof")
  fp <- sum(pred_lab == "gof" & y_true_chr == "lof")
  fn <- sum(pred_lab == "lof" & y_true_chr == "gof")
  
  tp <- as.numeric(tp); tn <- as.numeric(tn); fp <- as.numeric(fp); fn <- as.numeric(fn)
  
  sens <- ifelse((tp + fn) == 0, NA_real_, tp / (tp + fn)) # Recall GOF
  spec <- ifelse((tn + fp) == 0, NA_real_, tn / (tn + fp)) # Spec LOF
  ba <- mean(c(sens, spec), na.rm = TRUE)
  
  cm <- caret::confusionMatrix(
    factor(pred_lab, levels=c("gof","lof")),
    factor(y_true_chr, levels=c("gof","lof"))
  )
  kappa <- unname(cm$overall["Kappa"])
  
  data.table::data.table(
    thr = thr,
    TP = tp, FN = fn, FP = fp, TN = tn,
    BA = ba,
    Kappa = kappa,
    MCC = mcc(tp, tn, fp, fn),
    Recall_GOF = sens,
    Spec_LOF = spec
  )
}

compute_roc_pr <- function(df, name) {
  roc_obj <- pROC::roc(response = df$y_true, predictor = df$p_gof, levels = c("lof","gof"), direction = "<")
  roc_auc <- as.numeric(pROC::auc(roc_obj))
  
  pos_scores <- df$p_gof[df$y_true == "gof"]
  neg_scores <- df$p_gof[df$y_true == "lof"]
  pr <- PRROC::pr.curve(scores.class0 = pos_scores, scores.class1 = neg_scores, curve = FALSE)
  prauc <- pr$auc.integral
  
  cat("\n---", name, "---\n")
  cat("ROC AUC:", round(roc_auc, 4), "\n")
  cat("prAUC :", round(prauc, 4), "\n")
  
  list(roc_auc=roc_auc, prauc=prauc)
}


# Build unique table for all kept variants (Hybrid)

truth_all <- unique(feat_all[, c("variant_key","Class","gene_group")])
truth_all$y_true <- truth_all$Class
truth_all$Class <- NULL

truth_all <- data.table::as.data.table(truth_all)
agg_dt <- data.table::as.data.table(agg)

# Merge OOF p_gof for MIXED variants
hyb <- merge(truth_all, agg_dt[, .(variant_key, p_gof)], by="variant_key", all.x=TRUE)

# Deterministic p_gof for SINGLE genes
hyb$p_gof <- ifelse(hyb$gene_group == "single_lof", 0,
                    ifelse(hyb$gene_group == "single_gof", 1, hyb$p_gof))

cat("\n[HYBRID] Unique variants table:\n")
cat("n_unique:", nrow(hyb), "\n")
print(table(hyb$gene_group))
cat("Missing p_gof after hybrid (should be 0):", sum(is.na(hyb$p_gof)), "\n")

if (sum(is.na(hyb$p_gof)) > 0) {
  cat("\nVariants with missing p_gof (unexpected):\n")
  print(hyb[is.na(p_gof)])
  stop("Some hybrid p_gof are NA. This should not happen if all MIXED variants got OOF predictions.")
}


mix_only <- agg_dt
mix_only$y_true <- factor(mix_only$y_true, levels=c("gof","lof"))
hyb$y_true <- factor(hyb$y_true, levels=c("gof","lof"))

mix_auc <- compute_roc_pr(mix_only, "MIXED only (ML)  [EXP2-like]")
hyb_auc <- compute_roc_pr(hyb,      "ALL genes (Hybrid ML+Rule)  [EXP3]")

ths <- c(0.3, 0.5, 0.7)

out_mix <- data.table::rbindlist(lapply(ths, function(t) metrics_at_threshold(mix_only$y_true, mix_only$p_gof, t)))
out_hyb <- data.table::rbindlist(lapply(ths, function(t) metrics_at_threshold(hyb$y_true, hyb$p_gof, t)))

cat("\nMetrics @ thresholds (MIXED only):\n")
print(out_mix)

cat("\nMetrics @ thresholds (ALL Hybrid):\n")
print(out_hyb)