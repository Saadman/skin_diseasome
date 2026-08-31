# GO biological-process enrichment on the harmonized (equal-N) gene sets.
# Tests the published claim that neutrophil activation is enriched in every
# disease except SLS, and that SLS's profile is keratinisation-dominated.
suppressMessages({library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi)})

N <- 582
sets <- readRDS(sprintf("sensitivity/harmonized_symbols_N%d.rds", N))
universe <- readRDS("sensitivity/universe_symbols.rds")
uni_eg <- unique(na.omit(as.character(suppressMessages(
  mapIds(org.Hs.eg.db, keys = universe, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")))))

res <- list()
for (a in names(sets)) {
  eg <- unique(na.omit(as.character(suppressMessages(
    mapIds(org.Hs.eg.db, keys = sets[[a]], column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")))))
  e <- enrichGO(gene = eg, universe = uni_eg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                readable = TRUE)
  df <- as.data.frame(e)
  res[[a]] <- df
  cat(sprintf("\n=== %s : %d significant BP terms ===\n", a, nrow(df)))
  if (nrow(df)) {
    for (i in 1:min(6, nrow(df)))
      cat(sprintf("   %-58s p.adj=%.2e\n", substr(df$Description[i], 1, 57), df$p.adjust[i]))
    hit <- grep("neutrophil", df$Description, ignore.case = TRUE)
    cat("   neutrophil terms: ", if (length(hit)) paste(sprintf("%s (rank %d, p.adj=%.1e)",
        df$Description[hit][1], hit[1], df$p.adjust[hit][1])) else "NONE", "\n")
    ker <- grep("keratin|cornif|epiderm", df$Description, ignore.case = TRUE)
    cat("   keratinisation/epidermis terms: ", if (length(ker))
        sprintf("%d terms, top = %s (rank %d)", length(ker), df$Description[ker][1], ker[1]) else "NONE", "\n")
  }
}
saveRDS(res, sprintf("sensitivity/harmonized_GO_N%d.rds", N))
out <- do.call(rbind, lapply(names(res), function(a)
  if (nrow(res[[a]])) cbind(disease = a, res[[a]][, c("ID","Description","GeneRatio","p.adjust","Count")]) else NULL))
write.csv(out, sprintf("sensitivity/harmonized_GO_N%d.csv", N), row.names = FALSE)
cat("\nwrote sensitivity/harmonized_GO_N582.csv\n")
