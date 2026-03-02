#!/usr/bin/env Rscript
# ==================================================
# 01_cleaning.R
# Data download/load + clinical cleaning +
# baseline modeling-ready dataset + delta dataset
# ==================================================

# ---------- Setup ----------
dir.create("data", showWarnings = FALSE)
dir.create("data/raw", showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE)

# ---------- Packages ----------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
library(GEOquery)

# ---------- Load or Download GEO ----------
gse_id <- "GSE146446"
raw_rds <- file.path("data", "raw", paste0(gse_id, "_exprSet.rds"))

if (file.exists(raw_rds)) {
  message("Loading cached ExpressionSet: ", raw_rds)
  exprSet <- readRDS(raw_rds)
} else {
  message("Downloading GEO dataset: ", gse_id)
  exprSet <- getGEO(gse_id, GSEMatrix = TRUE)[[1]]
  saveRDS(exprSet, raw_rds)
  message("Saved ExpressionSet to: ", raw_rds)
}

# ---------- Extract core objects ----------
expr_matrix <- exprs(exprSet)
pheno_data  <- pData(exprSet)

message("expr_matrix dim: ", paste(dim(expr_matrix), collapse=" x "))
message("pheno_data  dim: ", paste(dim(pheno_data), collapse=" x "))

# ---------- Step 1: Align sample order ----------
pheno_data <- pheno_data[colnames(expr_matrix), ]
stopifnot(identical(rownames(pheno_data), colnames(expr_matrix)))

# ---------- Step 2: Standardize clinical fields ----------
# Required GEO columns:
# visit:ch1, response:ch1, treatment:ch1, subject_id:ch1
req_cols <- c("visit:ch1","response:ch1","treatment:ch1","subject_id:ch1")
missing <- setdiff(req_cols, colnames(pheno_data))
if (length(missing) > 0) stop("Missing required columns in pheno_data: ", paste(missing, collapse=", "))

pheno_data$visit      <- pheno_data$`visit:ch1`
pheno_data$response   <- pheno_data$`response:ch1`
pheno_data$treatment  <- pheno_data$`treatment:ch1`
pheno_data$subject_id <- pheno_data$`subject_id:ch1`

# Factors (optional but recommended)
pheno_data$visit     <- factor(pheno_data$visit, levels = c("T0","T8"))
pheno_data$treatment <- factor(pheno_data$treatment)
pheno_data$subject_id <- factor(pheno_data$subject_id)

# Response may be stored as character; convert to integer if possible
# (We keep original coding; modeling script will normalize 1/2 -> 0/1)
if (!is.numeric(pheno_data$response)) {
  suppressWarnings({
    pheno_data$response <- as.integer(as.character(pheno_data$response))
  })
}

message("Visit counts:"); print(table(pheno_data$visit))
message("Treatment counts:"); print(table(pheno_data$treatment))
message("Response counts:"); print(table(pheno_data$response, useNA="ifany"))

# ==================================================
# BASELINE DATASET (DLX, T0) -> data/DLX_baseline_modeling_ready.rds
# ==================================================
message("---- Building baseline modeling-ready dataset (DLX, T0) ----")
pd0   <- subset(pheno_data, visit == "T0")
expr0 <- expr_matrix[, rownames(pd0)]

# Keep only with response not NA
keep <- !is.na(pd0$response)
pd0 <- pd0[keep, ]
expr0 <- expr0[, rownames(pd0)]

# Remove duplicate subject_id (keep first occurrence)
dup <- duplicated(pd0$subject_id)
pd0 <- pd0[!dup, ]
expr0 <- expr0[, rownames(pd0)]

# DLX baseline only
pd_dlx <- subset(pd0, treatment == "DLX")
expr_dlx <- expr0[, rownames(pd_dlx)]

message("DLX baseline samples: ", nrow(pd_dlx))
message("Table response x treatment (baseline T0):")
print(table(pd0$response, pd0$treatment))

# y baseline (keep as integer; may be 0/1 or 1/2)
y_base <- as.integer(pd_dlx$response)

# Variance filter to top 5000 genes/probes
v <- apply(expr_dlx, 1, var)
topN <- 5000
keep_genes <- order(v, decreasing = TRUE)[1:min(topN, length(v))]
X_gs <- expr_dlx[keep_genes, , drop = FALSE]  # genes x samples
X_base <- t(X_gs)                              # samples x genes

stopifnot(nrow(X_base) == length(y_base))

base_out <- "data/DLX_baseline_modeling_ready.rds"
saveRDS(list(X = X_base, y = y_base, clinical = pd_dlx), base_out)
message("Saved baseline modeling-ready dataset: ", base_out)
message("X_base dim (samples x genes): ", paste(dim(X_base), collapse=" x "))

# ==================================================
# DELTA DATASET (DLX, paired T8-T0) -> data/processed/DLX_delta_data.rds
# ==================================================
message("---- Building delta dataset (DLX, paired T8-T0) ----")
pd_dlx_all <- subset(pheno_data, treatment == "DLX")
expr_dlx_all <- expr_matrix[, rownames(pd_dlx_all)]

message("DLX all samples: ", nrow(pd_dlx_all))
print(table(pd_dlx_all$visit))

# keep only subjects with exactly 2 timepoints
tab <- table(pd_dlx_all$subject_id)
subjects_two <- names(which(tab == 2))

pd_pair <- subset(pd_dlx_all, subject_id %in% subjects_two)
expr_pair <- expr_dlx_all[, rownames(pd_pair)]

# order by subject then visit (T0 then T8 due to factor levels)
pd_pair <- pd_pair[order(pd_pair$subject_id, pd_pair$visit), ]
expr_pair <- expr_pair[, rownames(pd_pair)]

# compute delta per subject: T8 - T0
subjects <- unique(pd_pair$subject_id)
delta_matrix <- matrix(NA_real_,
                       nrow = nrow(expr_pair),
                       ncol = length(subjects),
                       dimnames = list(rownames(expr_pair), as.character(subjects)))

for (i in seq_along(subjects)) {
  sid <- subjects[i]
  idx <- which(pd_pair$subject_id == sid)
  idx_T0 <- idx[pd_pair$visit[idx] == "T0"]
  idx_T8 <- idx[pd_pair$visit[idx] == "T8"]
  delta_matrix[, i] <- expr_pair[, idx_T8] - expr_pair[, idx_T0]
}

message("delta_matrix dim (genes x subjects): ", paste(dim(delta_matrix), collapse=" x "))

# response label from T0 rows (one per subject)
pd_T0 <- subset(pd_pair, visit == "T0")
pd_T0 <- pd_T0[match(colnames(delta_matrix), pd_T0$subject_id), ]
stopifnot(identical(as.character(pd_T0$subject_id), colnames(delta_matrix)))

y_delta <- as.integer(pd_T0$response)
message("y_delta counts:"); print(table(y_delta, useNA="ifany"))

# variance filter for delta to top 5000
v_delta <- apply(delta_matrix, 1, var)
topN <- 5000
top_genes_delta <- order(v_delta, decreasing = TRUE)[1:min(topN, length(v_delta))]
delta_top <- delta_matrix[top_genes_delta, , drop = FALSE]  # genes x subjects
message("delta_top dim (genes x subjects): ", paste(dim(delta_top), collapse=" x "))

delta_out <- "data/processed/DLX_delta_data.rds"
saveRDS(list(delta_matrix = delta_matrix, delta_top = delta_top, y_delta = y_delta, clinical_T0 = pd_T0), delta_out)
message("Saved delta dataset: ", delta_out)

message("DONE: 01_cleaning.R")

