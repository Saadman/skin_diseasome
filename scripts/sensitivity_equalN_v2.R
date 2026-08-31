# Extended sensitivity analysis: does the choice of RANKING STATISTIC change the
# conclusion? Ranking on raw p-value conflates effect size with sample size
# (PS has 170 samples, CSU 18), so it is tested against effect-size and
# t-statistic rankings, on a gene universe common to all six platforms.
suppressMessages({library(readxl); library(org.Hs.eg.db); library(AnnotationDbi)})

RAW <- list(
  PS  = list(f = "DEG/GEX_datasets/3-GSE30999_Psoriasis.xlsx",       s = "GSE30999_Psoriasis",  id = "Gene.symbol"),
  AD  = list(f = "DEG/GEX_datasets/2-GSE32924-AD-LvsNL.xlsx",        s = "GSE32924.top.table",  id = "Gene.symbol"),
  HS  = list(f = "DEG/GEX_datasets/5-GSE148027_hidradenitis.xlsx",   s = "Sheet1",              id = "Gene.symbol"),
  VL  = list(f = "DEG/GEX_datasets/6-GSE75819_vitilago.xlsx",        s = "Sheet1",              id = "Gene.symbol"),
  CSU = list(f = "DEG/GEX_datasets/4-GSE57178_urticaria.xlsx",       s = "Sheet1",              id = "Gene.symbol"),
  SLS = list(f = "DEG/GEX_datasets/1-sls_vs_control_GSE168735.xlsx", s = "GSE168735.top.table", id = "GENEID"))
ABB <- names(RAW)

tabs <- lapply(ABB, function(a) {
  x <- RAW[[a]]
  d <- suppressMessages(read_excel(x$f, sheet = x$s))
  sym <- as.character(d[[x$id]])
  if (x$id == "GENEID")
    sym <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db, keys = sym,
             column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"))
  sym <- sub("///.*$", "", trimws(sym))
  df <- data.frame(sym = sym,
                   p  = suppressWarnings(as.numeric(d$P.Value)),
                   fc = suppressWarnings(as.numeric(d$logFC)),
                   t  = suppressWarnings(as.numeric(d$t)), stringsAsFactors = FALSE)
  df <- df[!is.na(df$sym) & df$sym != "" & !is.na(df$p) & !is.na(df$fc), ]
  df <- df[order(df$p), ]
  df[!duplicated(df$sym), ]
})
names(tabs) <- ABB

universe <- Reduce(intersect, lapply(tabs, function(d) d$sym))
cat("Common gene universe across all six platforms:", length(universe), "genes\n\n")
tabs <- lapply(tabs, function(d) d[d$sym %in% universe, ])

pairs <- t(combn(ABB, 2))
report <- function(tops, label) {
  res <- data.frame(pair = paste(pairs[,1], pairs[,2], sep = "-"), n = NA, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(pairs))) res$n[i] <- length(intersect(tops[[pairs[i,1]]], tops[[pairs[i,2]]]))
  res <- res[order(-res$n), ]
  sls <- res[grepl("SLS", res$pair), ]; nonsls <- res[!grepl("SLS", res$pair), ]
  cat(sprintf("%-34s top1=%-7s PS-HS rank %2d/15 | SLS mean %6.1f vs other %6.1f | SLS worst-rank %d\n",
      label, paste0(res$pair[1], "(", res$n[1], ")"), which(res$pair == "PS-HS"),
      mean(sls$n), mean(nonsls$n), max(which(res$pair %in% sls$pair))))
  invisible(res)
}

for (N in c(500, 1000, 2000)) {
  cat("=== N =", N, "===\n")
  report(lapply(tabs, function(d) head(d$sym[order(d$p)], N)),                       "rank by p-value")
  report(lapply(tabs, function(d) head(d$sym[order(-abs(d$fc))], N)),                "rank by |log2FC|")
  report(lapply(tabs, function(d) head(d$sym[order(-abs(d$t))], N)),                 "rank by |t|")
  report(lapply(tabs, function(d) { e <- d[d$p < 0.05, ]; head(e$sym[order(-abs(e$fc))], N) }),
                                                                                     "p<0.05 then |log2FC|")
  cat("\n")
}
