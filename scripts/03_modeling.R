# =========================
# 03_modeling.R (Upgraded)
# Baseline vs Delta models
#   - v1: 70/30 split
#   - v2: 5-fold CV (out-of-fold ROC/AUC)
# =========================

# Avoid X11/XQuartz dependency on macOS when saving PNG

dir.create("results", showWarnings = FALSE)

# ---------- Packages ----------
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
library(pROC)
if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")

# ---------- Helpers ----------
normalize_y01 <- function(y) {
  if (is.factor(y)) y <- as.integer(as.character(y))
  if (!is.numeric(y)) y <- as.integer(y)
  u <- sort(unique(y))
  if (identical(u, c(1, 2))) return(ifelse(y == 2, 1, 0))
  if (identical(u, c(0, 1))) return(y)
  if (length(u) == 2) return(ifelse(y == max(u), 1, 0))
  stop("y has unexpected classes: ", paste(u, collapse = ", "))
}

ensure_samples_by_features <- function(X, y) {
  if (nrow(X) == length(y)) return(X)
  if (ncol(X) == length(y)) return(t(X))
  stop("X shape doesn't match y. dim(X)=",
       paste(dim(X), collapse="x"), " length(y)=", length(y))
}

make_stratified_folds <- function(y01, k = 5, seed = 123) {
  set.seed(seed)
  idx0 <- sample(which(y01 == 0))
  idx1 <- sample(which(y01 == 1))

  folds <- vector("list", k)
  for (i in seq_len(k)) folds[[i]] <- integer(0)

  for (i in seq_along(idx0)) folds[[ (i - 1) %% k + 1 ]] <- c(folds[[ (i - 1) %% k + 1 ]], idx0[i])
  for (i in seq_along(idx1)) folds[[ (i - 1) %% k + 1 ]] <- c(folds[[ (i - 1) %% k + 1 ]], idx1[i])

  folds
}

# ---------- v1: 70/30 split logistic ----------
run_logistic_split <- function(X, y01, out_prefix, seed = 123, topK = 100) {
  set.seed(seed)
  stopifnot(nrow(X) == length(y01))

  n <- nrow(X)
  train_idx <- sample(seq_len(n), size = floor(0.7 * n))
  test_idx  <- setdiff(seq_len(n), train_idx)

  X_train <- X[train_idx, , drop = FALSE]
  X_test  <- X[test_idx, , drop = FALSE]
  y_train <- y01[train_idx]
  y_test  <- y01[test_idx]

  if (length(unique(y_train)) < 2 || length(unique(y_test)) < 2) {
    stop("Train/test split resulted in only one class. Try a different seed.")
  }

  # Feature selection on TRAIN only (avoid leakage)
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

ragg::agg_png(paste0("results/ROC_", out_prefix, ".png"), width = 900, height = 900, units = "px")
  plot(roc_obj, main = paste0("ROC — ", out_prefix, " (Split)  AUC=", round(auc_val, 3)))
  dev.off()

  metrics_file <- paste0("results/model_", out_prefix, "_metrics.txt")
  con <- file(metrics_file, open = "wt")
  writeLines(paste0("Model: ", out_prefix, " (70/30 split)"), con)
  writeLines(paste0("Seed: ", seed), con)
  writeLines(paste0("Train size: ", length(train_idx)), con)
  writeLines(paste0("Test size: ", length(test_idx)), con)
  writeLines(paste0("TopK (variance on TRAIN only): ", topK), con)
  writeLines(paste0("AUC: ", round(auc_val, 4)), con)
  writeLines("\nConfusion Matrix (threshold=0.5):", con)
  writeLines(capture.output(cm), con)
  close(con)

  list(auc = auc_val, cm = cm, model = m)
}

# ---------- v2: 5-fold CV logistic (OOF ROC/AUC) ----------
run_logistic_cv <- function(X, y01, out_prefix, seed = 123, k = 5, topK = 100) {
  stopifnot(nrow(X) == length(y01))
  folds <- make_stratified_folds(y01, k = k, seed = seed)

  oof_prob <- rep(NA_real_, length(y01))
  fold_auc <- numeric(k)

  for (i in seq_len(k)) {
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_len(nrow(X)), test_idx)

    X_train <- X[train_idx, , drop = FALSE]
    X_test  <- X[test_idx, , drop = FALSE]
    y_train <- y01[train_idx]
    y_test  <- y01[test_idx]

    v <- apply(X_train, 2, var)
    top_features <- order(v, decreasing = TRUE)[1:min(topK, length(v))]
    X_train_small <- X_train[, top_features, drop = FALSE]
    X_test_small  <- X_test[, top_features, drop = FALSE]

    df_train <- data.frame(y = y_train, X_train_small)
    m <- glm(y ~ ., data = df_train, family = binomial)

    prob <- predict(m, newdata = data.frame(X_test_small), type = "response")
    oof_prob[test_idx] <- prob

    roc_obj <- roc(y_test, prob, quiet = TRUE)
    fold_auc[i] <- as.numeric(auc(roc_obj))
  }

  roc_oof <- roc(y01, oof_prob, quiet = TRUE)
  auc_oof <- as.numeric(auc(roc_oof))

ragg::agg_png(paste0("results/ROC_", out_prefix, "_cv.png"), width = 900, height = 900, units = "px")
  plot(roc_oof, main = paste0("ROC — ", out_prefix, " (5-fold CV OOF)  AUC=", round(auc_oof, 3)))
  dev.off()

  metrics_file <- paste0("results/model_", out_prefix, "_cv_metrics.txt")
  con <- file(metrics_file, open = "wt")
  writeLines(paste0("Model: ", out_prefix, " (", k, "-fold CV, out-of-fold)"), con)
  writeLines(paste0("Seed: ", seed), con)
  writeLines(paste0("TopK (variance on TRAIN only): ", topK), con)
  writeLines(paste0("OOF AUC: ", round(auc_oof, 4)), con)
  writeLines(paste0("Fold AUC mean: ", round(mean(fold_auc), 4)), con)
  writeLines(paste0("Fold AUC sd: ", round(sd(fold_auc), 4)), con)
  writeLines("\nFold AUCs:", con)
  writeLines(paste0("  ", paste(round(fold_auc, 4), collapse = ", ")), con)
  close(con)

  list(auc_oof = auc_oof, fold_auc = fold_auc, oof_prob = oof_prob)
}

# =========================
# Load Baseline data
# =========================
obj_base <- readRDS("data/DLX_baseline_modeling_ready.rds")
X_base <- obj_base$X
y_base <- obj_base$y
y_base01 <- normalize_y01(y_base)
X_base <- ensure_samples_by_features(X_base, y_base01)

# =========================
# Load Delta data
# =========================
obj_delta <- readRDS("data/processed/DLX_delta_data.rds")
delta_top <- obj_delta$delta_top
y_delta <- obj_delta$y_delta
y_delta01 <- normalize_y01(y_delta)
X_delta <- ensure_samples_by_features(t(delta_top), y_delta01)

# =========================
# Run split + CV
# =========================
seed <- 123
topK <- 100
kfold <- 5

res_base_split <- run_logistic_split(X_base, y_base01, out_prefix = "baseline", seed = seed, topK = topK)
res_delta_split <- run_logistic_split(X_delta, y_delta01, out_prefix = "delta", seed = seed, topK = topK)

res_base_cv <- run_logistic_cv(X_base, y_base01, out_prefix = "baseline", seed = seed, k = kfold, topK = topK)
res_delta_cv <- run_logistic_cv(X_delta, y_delta01, out_prefix = "delta", seed = seed, k = kfold, topK = topK)

cat("=== Split (70/30) AUC ===\n")
cat("Baseline AUC:", round(res_base_split$auc, 3), "\n")
cat("Delta AUC   :", round(res_delta_split$auc, 3), "\n\n")

cat("=== 5-fold CV (OOF) AUC ===\n")
cat("Baseline OOF AUC:", round(res_base_cv$auc_oof, 3),
    " | mean±sd:", round(mean(res_base_cv$fold_auc), 3), "±", round(sd(res_base_cv$fold_auc), 3), "\n")
cat("Delta OOF AUC   :", round(res_delta_cv$auc_oof, 3),
    " | mean±sd:", round(mean(res_delta_cv$fold_auc), 3), "±", round(sd(res_delta_cv$fold_auc), 3), "\n")
