# Build harmonized (equal-N) DEG sets for the mixed-threshold sensitivity analysis.
# Rule applied identically to all six diseases:
#   (1) restrict to genes measurable on all six platforms
#   (2) nominal p < 0.05
#   (3) rank by |log2FC| descending, take top N
# N is capped at 582 (the most SLS can supply under step 2).
# Writes Entrez-ID gene lists for separation.py and symbol lists for clusterProfiler.
suppressMessages({library(readxl); library(org.Hs.eg.db); library(AnnotationDbi)})

OUT <- "sensitivity"; dir.create(OUT, showWarnings = FALSE)
RAW <- list(
  PS  = c("DEG/GEX_datasets/3-GSE30999_Psoriasis.xlsx",       "GSE30999_Psoriasis",  "Gene.symbol"),
  AD  = c("DEG/GEX_datasets/2-GSE32924-AD-LvsNL.xlsx",        "GSE32924.top.table",  "Gene.symbol"),
  HS  = c("DEG/GEX_datasets/5-GSE148027_hidradenitis.xlsx",   "Sheet1",              "Gene.symbol"),
  VL  = c("DEG/GEX_datasets/6-GSE75819_vitilago.xlsx",        "Sheet1",              "Gene.symbol"),
  CSU = c("DEG/GEX_datasets/4-GSE57178_urticaria.xlsx",       "Sheet1",              "Gene.symbol"),
  SLS = c("DEG/GEX_datasets/1-sls_vs_control_GSE168735.xlsx", "GSE168735.top.table", "GENEID"))

tabs <- lapply(names(RAW), function(a) {
  x <- RAW[[a]]
  d <- suppressMessages(read_excel(x[1], sheet = x[2]))
  s <- as.character(d[[x[3]]])
  if (x[3] == "GENEID")
    s <- suppressMessages(mapIds(org.Hs.eg.db, keys = s, column = "SYMBOL",
                                 keytype = "ENTREZID", multiVals = "first"))
  s <- sub("///.*$", "", trimws(s))
  df <- data.frame(sym = s,
                   p  = suppressWarnings(as.numeric(d$P.Value)),
                   fc = suppressWarnings(as.numeric(d$logFC)), stringsAsFactors = FALSE)
  df <- df[!is.na(df$sym) & df$sym != "" & !is.na(df$p) & !is.na(df$fc), ]
  df <- df[order(df$p), ]
  df[!duplicated(df$sym), ]
})
names(tabs) <- names(RAW)

universe <- Reduce(intersect, lapply(tabs, function(d) d$sym))
cat("common gene universe:", length(universe), "\n")
saveRDS(universe, file.path(OUT, "universe_symbols.rds"))

for (N in c(250, 400, 582)) {
  sets <- lapply(tabs, function(d) {
    e <- d[d$sym %in% universe & d$p < 0.05, ]
    e <- e[order(-abs(e$fc)), ]
    head(e$sym, N)
  })
  stopifnot(all(sapply(sets, length) == N))          # genuine equal-N or fail loudly
  saveRDS(sets, file.path(OUT, sprintf("harmonized_symbols_N%d.rds", N)))
  for (a in names(sets)) {
    eg <- suppressMessages(mapIds(org.Hs.eg.db, keys = sets[[a]], column = "ENTREZID",
                                  keytype = "SYMBOL", multiVals = "first"))
    eg <- unique(na.omit(as.character(eg)))
    writeLines(eg, file.path(OUT, sprintf("%s_N%d.txt", a, N)))
  }
  cat(sprintf("N=%d: sets built (Entrez mapped: %s)\n", N,
      paste(sapply(names(sets), function(a)
        sprintf("%s=%d", a, length(readLines(file.path(OUT, sprintf("%s_N%d.txt", a, N)))))),
        collapse = " ")))
}
cat("done\n")
