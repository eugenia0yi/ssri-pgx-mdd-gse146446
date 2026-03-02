#!/usr/bin/env Rscript
# ==================================================
# 02_DE_analysis.R
# Baseline DE (DLX, T0) + Delta DE (DLX, T8-T0)
# Output CSV + Volcano + Heatmap
# ==================================================

dir.create("results", showWarnings = FALSE)

# ---------- Packages ----------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
library(limma)
library(pheatmap)

# ==================================================
# Part A: Baseline DE (DLX, T0) - exploratory
# ==================================================
message("---- Baseline DE (DLX, T0) ----")
obj_base <- readRDS("data/DLX_baseline_modeling_ready.rds")
X_base <- obj_base$X   # samples x genes
y_base <- obj_base$y

# normalize y to 0/1 if it is 1/2
if (is.factor(y_base)) y_base <- as.integer(as.character(y_base))
if (all(sort(unique(y_base)) == c(1,2))) y_base <- ifelse(y_base == 2, 1, 0)

# limma expects genes x samples
X_gs <- t(X_base)
group <- factor(ifelse(y_base == 1, "Responder", "Nonresponder"),
                levels = c("Nonresponder","Responder"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- lmFit(X_gs, design)
contr <- makeContrasts(RvsNR = Responder - Nonresponder, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contr))
res_base <- topTable(fit2, coef="RvsNR", number=Inf, adjust.method="BH")

write.csv(res_base, "results/DE_baseline_all.csv")
write.csv(res_base[1:100, ], "results/DE_signature_top100.csv")

message("Baseline min adj.P.Val: ", signif(min(res_base$adj.P.Val), 4))

# ==================================================
# Part B: Delta DE (DLX, paired T8-T0)
# ==================================================
message("---- Delta DE (DLX, T8-T0) ----")
obj_delta <- readRDS("data/processed/DLX_delta_data.rds")
delta_top <- obj_delta$delta_top  # genes x subjects
y_delta <- obj_delta$y_delta

# normalize y_delta to 0/1 if it is 1/2
if (is.factor(y_delta)) y_delta <- as.integer(as.character(y_delta))
if (all(sort(unique(y_delta)) == c(1,2))) y_delta <- ifelse(y_delta == 2, 1, 0)

group_d <- factor(ifelse(y_delta == 1, "Responder", "Nonresponder"),
                  levels = c("Nonresponder","Responder"))
design_d <- model.matrix(~ 0 + group_d)
colnames(design_d) <- levels(group_d)

fit <- lmFit(delta_top, design_d)
contr <- makeContrasts(RvsNR = Responder - Nonresponder, levels = design_d)
fit2 <- eBayes(contrasts.fit(fit, contr))
res_delta <- topTable(fit2, coef="RvsNR", number=Inf, adjust.method="BH")

write.csv(res_delta, "results/DE_delta_all.csv")

# signature file: if none pass FDR, still save top100 for downstream use
sig_delta <- subset(res_delta, adj.P.Val < 0.05 & abs(logFC) > 0.3)
write.csv(sig_delta, "results/DE_delta_signature.csv")
write.csv(res_delta[1:100, ], "results/DE_delta_top100.csv")

message("Delta min adj.P.Val: ", signif(min(res_delta$adj.P.Val), 4))
message("Delta significant genes (FDR<0.05 & |logFC|>0.3): ", nrow(sig_delta))

# ---------- Volcano (delta) ----------
png("results/volcano_delta.png", width = 1000, height = 800)
plot(res_delta$logFC, -log10(res_delta$P.Value),
     pch = 16,
     xlab = "logFC (ΔResponder vs ΔNonresponder)",
     ylab = "-log10(P)",
     main = "Volcano: Delta Expression (DLX, T8-T0)")
abline(h = -log10(0.05), lty = 2)
abline(v = c(-0.3, 0.3), lty = 2)
dev.off()

# ---------- Heatmap top30 (delta) ----------
top30 <- rownames(res_delta)[1:30]
mat <- delta_top[top30, , drop = FALSE]
mat_z <- t(scale(t(mat)))

ann <- data.frame(Response = group_d)
rownames(ann) <- colnames(mat_z)

png("results/heatmap_delta.png", width = 1000, height = 900)
pheatmap(mat_z, annotation_col = ann, show_colnames = FALSE,
         main = "Top 30 Delta Genes (DLX, T8-T0)")
dev.off()

message("DONE: 02_DE_analysis.R")

