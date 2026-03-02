# =========================
# 03_modeling.R
# Baseline vs Delta models
# =========================

install.packages("pROC")
library(pROC)

run_logistic_model <- function(X, y, out_prefix, seed = 123, topK = 100) {
  set.seed(seed)
  stopifnot(nrow(X) == length(y))
  n <- nrow(X)
  train_idx <- sample(seq_len(n), size = floor(0.7 * n))
  X_train <- X[train_idx, , drop = FALSE]
  X_test  <- X[-train_idx, , drop = FALSE]
  y_train <- y[train_idx]
  y_test  <- y[-train_idx]
  if (length(unique(y_train)) < 2 || length(unique(y_test)) < 2) {
    stop("Train/test split resulted in only one class. Re-run with a different seed.")
  }
  v <- apply(X_train, 2, var)
  top_features <- order(v, decreasing = TRUE)[1:min(topK, length(v))]
  X_train_small <- X_train[, top_features, drop = FALSE]
  X_test_small  <- X_test[, top_features, drop = FALSE]
  df_train <- data.frame(y = y_train, X_train_small)
  m <- glm(y ~ ., data = df_train, family = binomial)
  prob <- predict(m, newdata = data.frame(X_test_small), type = "response")
  pred <- ifelse(prob >= 0.5, 1, 0)
  cm <- table(Predicted = pred, Actual = y_test)
  roc_obj <- roc(y_test, prob, quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  png(paste0("results/ROC_", out_prefix, ".png"), width = 800, height = 800)
  plot(roc_obj, main = paste0("ROC — ", out_prefix, " model (AUC=", round(auc_val, 3), ")"))
  dev.off()
  metrics_file <- paste0("results/model_", out_prefix, "_metrics.txt")
  con2 <- file(metrics_file, open = "wt")
  writeLines(paste0("Model: ", out_prefix), con2)
  writeLines(paste0("Seed: ", seed), con2)
  writeLines(paste0("Train size: ", length(train_idx)), con2)
  writeLines(paste0("Test size: ", n - length(train_idx)), con2)
  writeLines(paste0("TopK features (variance, train-only): ", topK), con2)
  writeLines(paste0("AUC: ", round(auc_val, 4)), con2)
  writeLines("\nConfusion Matrix (threshold=0.5):", con2)
  writeLines(capture.output(cm), con2)
  close(con2)
  return(list(auc = auc_val, cm = cm, model = m))
}

dir.create("results", showWarnings = FALSE)

obj_base <- readRDS("data/DLX_baseline_modeling_ready.rds")
X_base <- obj_base$X
y_base <- obj_base$y
if (is.factor(y_base)) y_base <- as.integer(as.character(y_base))
if (all(sort(unique(y_base)) == c(1,2))) y_base <- ifelse(y_base == 2, 1, 0)
if (ncol(X_base) == length(y_base) && nrow(X_base) != length(y_base)) X_base <- t(X_base)
res_base <- run_logistic_model(X_base, y_base, out_prefix = "baseline", seed = 123, topK = 100)

obj_delta <- readRDS("data/processed/DLX_delta_data.rds")
delta_top <- obj_delta$delta_top
y_delta <- obj_delta$y_delta
if (is.factor(y_delta)) y_delta <- as.integer(as.character(y_delta))
if (all(sort(unique(y_delta)) == c(1,2))) y_delta <- ifelse(y_delta == 2, 1, 0)
X_delta <- t(delta_top)
res_delta <- run_logistic_model(X_delta, y_delta, out_prefix = "delta", seed = 123, topK = 100)

cat("Baseline AUC:", round(res_base$auc, 3), "\n")
cat("Delta AUC   :", round(res_delta$auc, 3), "\n")
