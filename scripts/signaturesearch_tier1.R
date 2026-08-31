# Phase 3 Tier 1: per-disease LINCS drug-repurposing reconstruction (signatureSearch).
# Builds an up/down DEG signature per disease from its raw GSE file, using each
# disease's established significance criterion (see Supplementary Table S2), maps to
# Entrez IDs, and runs gess_lincs() against the cached `lincs`
# reference (Subramanian et al. 2017 LINCS L1000 Level5 z-scores, EH3226).
#
# Run with working directory = repo root.

suppressMessages(library(readxl))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(AnnotationDbi))
suppressMessages(library(signatureSearch))
suppressMessages(library(signatureSearchData))

OUT_DIR <- "signaturesearch"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

sym2entrez <- function(symbols) {
  m <- AnnotationDbi::select(org.Hs.eg.db, keys = symbols, keytype = "SYMBOL", columns = "ENTREZID")
  unique(na.omit(m$ENTREZID))
}

build_signature_symbol <- function(rawfile, pcol) {
  d <- read_excel(rawfile)
  d <- d[!is.na(d[[pcol]]) & !is.na(d$logFC) & !is.na(d$Gene.symbol), ]
  d <- d[d[[pcol]] < 0.05 & abs(d$logFC) > 0.585, ]
  list(
    up = sym2entrez(unique(d$Gene.symbol[d$logFC > 0.585])),
    down = sym2entrez(unique(d$Gene.symbol[d$logFC < -0.585]))
  )
}

build_signature_entrez <- function(rawfile, pcol) {
  d <- read_excel(rawfile)
  d <- d[!is.na(d[[pcol]]) & !is.na(d$logFC) & !is.na(d$GENEID), ]
  d <- d[d[[pcol]] < 0.05 & abs(d$logFC) > 0.585, ]
  list(
    up = unique(as.character(d$GENEID[d$logFC > 0.585])),
    down = unique(as.character(d$GENEID[d$logFC < -0.585]))
  )
}

diseases <- list(
  PS  = list(file = "DEG/GEX_datasets/3-GSE30999_Psoriasis.xlsx", pcol = "P.Value", type = "symbol"),
  AD  = list(file = "DEG/GEX_datasets/2-GSE32924-AD-LvsNL.xlsx", pcol = "P.Value", type = "symbol"),
  HS  = list(file = "DEG/GEX_datasets/5-GSE148027_hidradenitis.xlsx", pcol = "adj.P.Val", type = "symbol"),
  VL  = list(file = "DEG/GEX_datasets/6-GSE75819_vitilago.xlsx", pcol = "adj.P.Val", type = "symbol"),
  CSU = list(file = "DEG/GEX_datasets/4-GSE57178_urticaria.xlsx", pcol = "P.Value", type = "symbol"),
  SLS = list(file = "DEG/GEX_datasets/1-sls_vs_control_GSE168735.xlsx", pcol = "P.Value", type = "entrez")
)

for (abb in names(diseases)) {
  spec <- diseases[[abb]]
  message("=== ", abb, " ===")
  sig <- if (spec$type == "symbol") {
    build_signature_symbol(spec$file, spec$pcol)
  } else {
    build_signature_entrez(spec$file, spec$pcol)
  }
  message(abb, " up=", length(sig$up), " down=", length(sig$down))
  if (length(sig$up) < 5 || length(sig$down) < 5) {
    message(abb, ": signature too small (<5 genes in a direction), skipping LINCS search")
    next
  }
  qsig <- qSig(query = list(upset = sig$up, downset = sig$down), gess_method = "LINCS", refdb = "lincs")
  res <- gess_lincs(qsig, sortby = "NCS", tau = FALSE, workers = 2)
  df <- result(res)
  saveRDS(df, file.path(OUT_DIR, paste0(abb, "_lincs_result.rds")))
  write.csv(df, file.path(OUT_DIR, paste0(abb, "_lincs_result.csv")), row.names = FALSE)
  message(abb, " done: ", nrow(df), " rows written")
}
message("Tier 1 complete.")
