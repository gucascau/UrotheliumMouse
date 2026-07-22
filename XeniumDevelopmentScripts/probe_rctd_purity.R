suppressPackageStartupMessages({library(Seurat); library(dplyr)})
OUT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/XeniumDevelopmentScripts/output"
object <- readRDS(file.path(OUT_DIR, "XeniumDev_RCTD_deconvolved.rds"))

wcols <- grep("^rctd_", colnames(object@meta.data), value = TRUE)
wcols <- setdiff(wcols, "rctd_dominant_celltype")
message("Weight columns (", length(wcols), "): ", paste(head(wcols, 10), collapse=", "), " ...")

md <- object@meta.data
uro_idx <- which(md$rctd_dominant_celltype == "Urothelium")
message("Urothelium-dominant cells: ", length(uro_idx))

wmat <- as.matrix(md[uro_idx, wcols])
top1 <- apply(wmat, 1, max)
top2 <- apply(wmat, 1, function(x) sort(x, decreasing = TRUE)[2])
second_type <- apply(wmat, 1, function(x) wcols[order(x, decreasing = TRUE)[2]])

message("Urothelium weight (top1) summary:")
print(summary(top1))
message("Gap to 2nd-best type summary (top1 - top2):")
print(summary(top1 - top2))

message("Distribution of 2nd-highest-weight type among Urothelium-dominant cells:")
print(sort(table(sub("^rctd_", "", second_type)), decreasing = TRUE))

message("\nFraction of Urothelium cells with top weight > various thresholds:")
for (thr in c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)) {
  message(sprintf("  weight > %.1f : %d / %d (%.1f%%)", thr, sum(top1 > thr), length(top1), 100*mean(top1 > thr)))
}
