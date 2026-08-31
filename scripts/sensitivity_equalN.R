# Sensitivity analysis for the mixed-threshold critique.
# Removes both the threshold and the set-size confound by ranking each disease's
# genes on p-value and taking equal-sized top-N sets across all six diseases.
# Also reports size-normalised overlap (Jaccard) for the published DEG sets.
suppressMessages({library(readxl); library(org.Hs.eg.db); library(AnnotationDbi)})

RAW <- list(
  PS  = list(f = "DEG/GEX_datasets/3-GSE30999_Psoriasis.xlsx",        s = "GSE30999_Psoriasis",   id = "Gene.symbol"),
  AD  = list(f = "DEG/GEX_datasets/2-GSE32924-AD-LvsNL.xlsx",         s = "GSE32924.top.table",   id = "Gene.symbol"),
  HS  = list(f = "DEG/GEX_datasets/5-GSE148027_hidradenitis.xlsx",    s = "Sheet1",               id = "Gene.symbol"),
  VL  = list(f = "DEG/GEX_datasets/6-GSE75819_vitilago.xlsx",         s = "Sheet1",               id = "Gene.symbol"),
  CSU = list(f = "DEG/GEX_datasets/4-GSE57178_urticaria.xlsx",        s = "Sheet1",               id = "Gene.symbol"),
  SLS = list(f = "DEG/GEX_datasets/1-sls_vs_control_GSE168735.xlsx",  s = "GSE168735.top.table",  id = "GENEID"))
ABB <- names(RAW)

# ranked gene lists: one row per gene symbol, best (smallest) p-value, ascending
ranked <- lapply(ABB, function(a) {
  x <- RAW[[a]]
  d <- suppressMessages(read_excel(x$f, sheet = x$s))
  sym <- as.character(d[[x$id]])
  if (x$id == "GENEID") {                      # SLS is Entrez-only: map to symbol
    sym <- suppressMessages(AnnotationDbi::mapIds(
      org.Hs.eg.db, keys = sym, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"))
  }
  sym <- trimws(sym)
  sym <- sub("///.*$", "", sym)                # multi-mapping probes: keep first symbol
  p <- suppressWarnings(as.numeric(d$P.Value))
  ok <- !is.na(sym) & sym != "" & !is.na(p)
  df <- data.frame(sym = sym[ok], p = p[ok], stringsAsFactors = FALSE)
  df <- df[order(df$p), ]
  df[!duplicated(df$sym), ]
})
names(ranked) <- ABB
cat("Ranked genes available per disease:\n"); print(sapply(ranked, nrow)); cat("\n")

pairs <- t(combn(ABB, 2))
jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))

# ---- published DEG sets: size-normalised overlap -------------------------
pub <- read_excel("DEG/GEX_datasets/Venn_allDisease_final.xlsx")
psets <- lapply(colnames(pub), function(c) unique(na.omit(trimws(as.character(pub[[c]])))))
names(psets) <- colnames(pub)
pj <- data.frame(pair = paste(pairs[, 1], pairs[, 2], sep = "-"),
                 n_overlap = NA, jaccard = NA, stringsAsFactors = FALSE)
for (i in seq_len(nrow(pairs))) {
  A <- psets[[pairs[i, 1]]]; B <- psets[[pairs[i, 2]]]
  pj$n_overlap[i] <- length(intersect(A, B)); pj$jaccard[i] <- round(jac(A, B), 4)
}
cat("PUBLISHED SETS - overlap ranked by raw count vs by Jaccard:\n")
print(head(pj[order(-pj$n_overlap), ], 6), row.names = FALSE)
cat("\n"); print(head(pj[order(-pj$jaccard), ], 6), row.names = FALSE); cat("\n")

# ---- equal-N top ranked sets -------------------------------------------
for (N in c(250, 500, 1000, 2000)) {
  tops <- lapply(ranked, function(d) head(d$sym, N))
  res <- data.frame(pair = paste(pairs[, 1], pairs[, 2], sep = "-"), n = NA, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(pairs))) res$n[i] <- length(intersect(tops[[pairs[i, 1]]], tops[[pairs[i, 2]]]))
  res <- res[order(-res$n), ]
  sls <- res[grepl("SLS", res$pair), ]
  cat(sprintf("--- top-%d ---\n", N))
  cat("  rank1:", res$pair[1], "=", res$n[1],
      "| rank2:", res$pair[2], "=", res$n[2],
      "| rank3:", res$pair[3], "=", res$n[3], "\n")
  cat("  PS-HS rank:", which(res$pair == "PS-HS"), "of 15 (n =", res$n[res$pair == "PS-HS"], ")\n")
  cat("  SLS pairs (n):", paste(sprintf("%s=%d", sls$pair, sls$n), collapse = " "), "\n")
  cat("  SLS mean =", round(mean(sls$n), 1), "| non-SLS mean =", round(mean(res$n[!grepl("SLS", res$pair)]), 1), "\n\n")
}
