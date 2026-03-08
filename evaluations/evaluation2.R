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


# EVALUATION2 CHANGE: Filter to selected genes ONLY

keep_genes <- c("SCN5A","SCN1A","SCN2A","CACNA1A","SCN4A")
varall <- varall[varall$gene %in% keep_genes, ]

# SANITY CHECK (short): ensure ONLY these genes remain
genes_present <- sort(unique(varall$gene))
if (!all(genes_present %in% keep_genes)) {
  stop(paste0("Sanity check failed: found unexpected genes: ",
              paste(setdiff(genes_present, keep_genes), collapse=", ")))
}
cat("\n[Sanity check] Genes included (should be exactly the 5 selected):\n")
print(table(varall$gene))
cat("[Sanity check] n_genes:", length(genes_present), "\n")
cat("[Sanity check] genes_present:", paste(genes_present, collapse=", "), "\n")

# Deduplicate: unique by (gene, pos, altAA)
varall <- varall[!duplicated(varall[, c("gene","altAA","pos")]), ]

# key for aggregation
varall$variant_key <- paste(varall$gene, varall$pos, varall$altAA, sep=":")

cat("\nAfter Heyne filtering (+ gene subset):\n")
cat("n_rows:", nrow(varall), "\n")
print(table(varall$prd_mech_revised))


# Read features
featuretable <- data.table::fread(
  featuretable_path,
  sep = "\t",
  data.table = FALSE,
  showProgress = TRUE
)

if (!("protid" %in% colnames(featuretable))) stop("featuretable must contain column 'protid'")

# Join by protid
feat <- featuretable[match(varall$protid, featuretable$protid), ]

# Attach labels
feat$Class <- factor(varall$prd_mech_revised, levels = c("gof","lof"))  # GOF positive
feat$variant_key <- varall$variant_key


# OUR feature mapping

# 1) S4 = IIIS4 only
if (!("S4M3" %in% colnames(feat))) stop("Missing S4M3 in Heyne featuretable (needed for IIIS4).")
feat$S4 <- feat$S4M3

# 2) bury/medium/exposed continuous -> binary one-hot (argmax)
need_bme <- c("bury","medium","exposed")
if (!all(need_bme %in% colnames(feat))) stop("Missing bury/medium/exposed in Heyne featuretable.")
bme <- as.matrix(feat[, need_bme, drop=FALSE])
max_idx <- max.col(bme, ties.method = "first")  # deterministic
feat$bury    <- as.integer(max_idx == 1)
feat$medium  <- as.integer(max_idx == 2)
feat$exposed <- as.integer(max_idx == 3)

# 3) gofsotherKCNQgenes analogue
need_gofs <- c("scn","gofsotherSCNgenes","gofsotherCACgenes")
if (!all(need_gofs %in% colnames(feat))) stop("Missing scn/gofsotherSCNgenes/gofsotherCACgenes in Heyne featuretable.")
feat$gofsotherKCNQgenes <- ifelse(feat$scn == 1, feat$gofsotherSCNgenes, feat$gofsotherCACgenes)

# 4) kcnqcon analogue
need_con <- c("scn","scncon","caccon")
if (!all(need_con %in% colnames(feat))) stop("Missing scn/scncon/caccon in Heyne featuretable.")
feat$kcnqcon <- ifelse(feat$scn == 1, feat$scncon, feat$caccon)

# 5) 50/95 mapping into your names (50.kcnq / 95.kcnq)
has_50scn <- "50.scn" %in% colnames(feat)
has_95scn <- "95.scn" %in% colnames(feat)
has_50scncac <- "50.scncac" %in% colnames(feat)
has_95scncac <- "95.scncac" %in% colnames(feat)

if (has_50scn && has_50scncac) {
  feat[["50.kcnq"]] <- ifelse(feat$scn == 1, feat[["50.scn"]], feat[["50.scncac"]])
} else if (has_50scn) {
  feat[["50.kcnq"]] <- feat[["50.scn"]]
} else if (has_50scncac) {
  feat[["50.kcnq"]] <- feat[["50.scncac"]]
} else {
  stop("Missing both 50.scn and 50.scncac in Heyne featuretable.")
}

if (has_95scn && has_95scncac) {
  feat[["95.kcnq"]] <- ifelse(feat$scn == 1, feat[["95.scn"]], feat[["95.scncac"]])
} else if (has_95scn) {
  feat[["95.kcnq"]] <- feat[["95.scn"]]
} else if (has_95scncac) {
  feat[["95.kcnq"]] <- feat[["95.scncac"]]
} else {
  stop("Missing both 95.scn and 95.scncac in Heyne featuretable.")
}

# Define our feature list
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

missing <- setdiff(our_features, colnames(feat))
if (length(missing) > 0) stop(paste("Missing required features:", paste(missing, collapse=", ")))

X <- feat[, our_features, drop=FALSE]
y <- feat$Class

# Drop near-zero variance columns (optional)
nzv <- caret::nearZeroVar(X)
if (length(nzv) > 0) X <- X[, -nzv, drop=FALSE]

cat("\nDesign matrix:\n")
cat("n_rows:", nrow(X), " | n_features:", ncol(X), "\n")
print(table(y))


# Run the model
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
  x = X, y = y,
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


# OOF predictions

pred <- model$pred
best <- model$bestTune

pred_best <- pred[
  pred$n.trees == best$n.trees &
    pred$interaction.depth == best$interaction.depth &
    pred$shrinkage == best$shrinkage &
    pred$n.minobsinnode == best$n.minobsinnode, ]

pred_best$variant_key <- feat$variant_key[pred_best$rowIndex]
pred_best$y_true <- feat$Class[pred_best$rowIndex]

# GOF prob
if (!("gof" %in% colnames(pred_best))) stop("Missing 'gof' probability column in caret predictions.")
pred_best$p_gof <- pred_best$gof

# IMPORTANT FIX: convert to data.table so "by=" works
pred_best <- data.table::as.data.table(pred_best)

agg <- pred_best[, .(
  y_true = as.character(stats::na.omit(y_true)[1]),
  p_gof = mean(p_gof, na.rm = TRUE)
), by = variant_key]

agg$y_true <- factor(agg$y_true, levels=c("gof","lof"))

cat("\nUnique variants after aggregation:", nrow(agg), "\n")
print(table(agg$y_true))


# Evaluation Metrics

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

roc_obj <- pROC::roc(response = agg$y_true, predictor = agg$p_gof, levels = c("lof","gof"), direction = "<")
roc_auc <- as.numeric(pROC::auc(roc_obj))

pos_scores <- agg$p_gof[agg$y_true == "gof"]
neg_scores <- agg$p_gof[agg$y_true == "lof"]
pr <- PRROC::pr.curve(scores.class0 = pos_scores, scores.class1 = neg_scores, curve = FALSE)
prauc <- pr$auc.integral

cat("\nThreshold-independent:\n")
cat("ROC AUC:", round(roc_auc, 4), "\n")
cat("prAUC :", round(prauc, 4), "\n")

ths <- c(0.3, 0.5, 0.7)
out <- data.table::rbindlist(lapply(ths, function(t) metrics_at_threshold(agg$y_true, agg$p_gof, t)))

cat("\nMetrics @ thresholds:\n")
print(out)