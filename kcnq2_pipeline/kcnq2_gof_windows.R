# GOF density features (10aa + 3aa) calculations on KCNQ2

# INPUT
gof_pos <- c(106,141,144,144,144,153,175,192,194,198,201,201,317,318,335)
lof_pos <- c(
  1,107,113,114,126,127,130,178,185,196,203,205,207,213,213,214,239,247,253,254,
  254,256,258,265,265,265,265,267,268,268,268,268,269,272,274,274,276,277,277,
  281,284,285,285,287,290,290,294,294,295,296,301,304,306,309,313,315,325,333,
  337,339,341,344,352,353,355,358,359,362,386,489,541,541,546,546,547,553,553,
  554,556,556,558,560,563,567,581,581,588,592,637,764
)

L <- 872  # KCNQ2 sequence length

# Build full variant table (one row per variant occurrence)
varallmod <- rbind(
  data.frame(gene="KCNQ2", pos=gof_pos, Class="gof", stringsAsFactors=FALSE),
  data.frame(gene="KCNQ2", pos=lof_pos, Class="lof", stringsAsFactors=FALSE)
)

# Add a row id to keep duplicates distinguishable
varallmod$row_id <- seq_len(nrow(varallmod))

# Data Splits
suppressWarnings(RNGversion("3.5.3"))
set.seed(999)

inTraining  <- caret::createDataPartition(as.factor(varallmod$Class), p = 0.9, list = FALSE)
trainingall <- varallmod[inTraining, , drop=FALSE]
testing     <- varallmod[-inTraining, , drop=FALSE]

set.seed(989)
inTraining1 <- caret::createDataPartition(as.factor(trainingall$Class), p = 0.5, list = FALSE)
training1   <- trainingall[inTraining1, , drop=FALSE]
training2   <- trainingall[-inTraining1, , drop=FALSE]

# Density machinery (moving average with circular wrap)
ma <- function(x, w) {
  as.numeric(stats::filter(x, rep(1/w, w), circular = TRUE))
}

dens_at_pos <- function(pos, y, w, subtract_self = 0L) {
  y2 <- y
  y2[pos] <- y2[pos] - subtract_self
  dens <- ma(y2, w)
  dens[pos]
}

# Build GOF counts vector from TRAINING1 GOF only
training1_gof <- training1[training1$Class == "gof", , drop=FALSE]

y_train1_gof <- integer(L)
tab_gof <- table(training1_gof$pos)
y_train1_gof[as.integer(names(tab_gof))] <- as.integer(tab_gof)

# Add GOF density features to any dataset, mapped from training1 GOF
add_gof_density_features <- function(df, y_train1_gof, do_wovar = FALSE) {
  subtract_self_vec <- if (do_wovar) rep(1L, nrow(df)) else rep(0L, nrow(df))
  
  raw10 <- mapply(dens_at_pos, pos = df$pos,
                  MoreArgs = list(y = y_train1_gof, w = 10),
                  subtract_self = subtract_self_vec)
  
  raw3  <- mapply(dens_at_pos, pos = df$pos,
                  MoreArgs = list(y = y_train1_gof, w = 3),
                  subtract_self = subtract_self_vec)
  
  df$densgofwovar_raw10 <- as.numeric(raw10)
  df$densgofwovar_raw3  <- as.numeric(raw3)
  df
}


training1_gof <- add_gof_density_features(training1_gof, y_train1_gof, do_wovar = TRUE)

training1_lof <- training1[training1$Class == "lof", , drop=FALSE]
if (nrow(training1_lof) > 0) {
  training1_lof <- add_gof_density_features(training1_lof, y_train1_gof, do_wovar = FALSE)
}

training2 <- add_gof_density_features(training2, y_train1_gof, do_wovar = FALSE)
testing   <- add_gof_density_features(testing,   y_train1_gof, do_wovar = FALSE)

training1 <- rbind(training1_gof, training1_lof)
training1 <- training1[order(training1$row_id), , drop=FALSE]

# Z-score calculations
zfit <- function(x, mu, sd) (x - mu) / sd

mu10 <- mean(training1_gof$densgofwovar_raw10)
sd10 <- sd(training1_gof$densgofwovar_raw10)

mu3  <- mean(training1_gof$densgofwovar_raw3)
sd3  <- sd(training1_gof$densgofwovar_raw3)

if (is.na(sd10) || sd10 == 0) stop("sd10 is 0/NA: training1_gof has too little variation for 10aa density.")
if (is.na(sd3)  || sd3  == 0) stop("sd3 is 0/NA: training1_gof has too little variation for 3aa density.")

training1$densgofwovar    <- round(zfit(training1$densgofwovar_raw10, mu10, sd10), 2)
training1$densgofwovar3aa <- round(zfit(training1$densgofwovar_raw3,  mu3,  sd3),  2)

training2$densgofwovar    <- round(zfit(training2$densgofwovar_raw10, mu10, sd10), 2)
training2$densgofwovar3aa <- round(zfit(training2$densgofwovar_raw3,  mu3,  sd3),  2)

testing$densgofwovar      <- round(zfit(testing$densgofwovar_raw10,   mu10, sd10), 2)
testing$densgofwovar3aa   <- round(zfit(testing$densgofwovar_raw3,    mu3,  sd3),  2)


# FINAL MERGE + SANITY CHECKS + SORTED EXPORT
final_all <- rbind(training1, training2, testing)

expected_n <- length(gof_pos) + length(lof_pos)

cat("Expected rows:", expected_n, "\n")
cat("Actual rows   :", nrow(final_all), "\n")

# Check missing row_ids (should be empty)
missing_ids <- setdiff(1:expected_n, final_all$row_id)
cat("Missing row_ids:", ifelse(length(missing_ids)==0, "NONE", paste(missing_ids, collapse=", ")), "\n")

# Check duplicated row_ids (should be empty)
dup_ids <- final_all$row_id[duplicated(final_all$row_id)]
cat("Duplicated row_ids:", ifelse(length(dup_ids)==0, "NONE", paste(dup_ids, collapse=", ")), "\n")

# If something is missing, print the original rows that correspond to the missing ids:
if (length(missing_ids) > 0) {
  cat("\nMissing rows from original varallmod:\n")
  print(varallmod[varallmod$row_id %in% missing_ids, ])
}

# Sort by position
final_all <- final_all[order(final_all$pos, final_all$Class, final_all$row_id), ]

final_features <- final_all[, c("row_id","gene","pos","Class","densgofwovar","densgofwovar3aa")]

write.table(final_features,
            file = "KCNQ2_GOF_density_features_3aa_10aa.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Return final table in the console
final_features