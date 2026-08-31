# Tier 2: cross-disease comparison of LINCS reversal candidates.
#
# Reads the per-disease gess_lincs() output written by scripts/signaturesearch_tier1.R,
# reduces each to its top-N distinct reversal compounds, counts the compounds shared by
# each disease pair, and joins the counts to the interactome separation values so that
# drug-signature overlap can be compared with network proximity.
#
# Selection rule (the one used for the values reported in the manuscript):
#   1. keep only signatures with trend == "down" (i.e. the compound moves the
#      disease-associated genes opposite to the disease itself: candidate reversal);
#   2. collapse to one row per compound, keeping its best-scoring cell line
#      (gess_lincs sorts by |NCS| descending, so the first occurrence is the best);
#   3. take the top TOP_N compounds.
#
# Run with working directory = repo root.

suppressMessages(library(readxl))
suppressMessages(library(dplyr))

IN_DIR   <- "signaturesearch"
SAB_FILE <- "networksim/sab_values.csv"
TOP_N    <- 50
DISEASES <- c("PS", "AD", "HS", "VL", "CSU", "SLS")

top_reversal <- function(abb, n = TOP_N) {
  d <- read.csv(file.path(IN_DIR, paste0(abb, "_lincs_result.csv")), stringsAsFactors = FALSE)
  d <- d[d$trend == "down", ]
  head(d$pert[!duplicated(d$pert)], n)
}

tops <- lapply(DISEASES, top_reversal)
names(tops) <- DISEASES

pairs <- t(combn(DISEASES, 2))
overlap <- data.frame(
  A = pairs[, 1],
  B = pairs[, 2],
  overlap_n = apply(pairs, 1, function(p) length(intersect(tops[[p[1]]], tops[[p[2]]]))),
  overlap_drugs = apply(pairs, 1, function(p) paste(intersect(tops[[p[1]]], tops[[p[2]]]), collapse = ";")),
  stringsAsFactors = FALSE
)
overlap <- overlap[order(-overlap$overlap_n), ]
write.csv(overlap, file.path(IN_DIR, "tier2_pairwise_overlap.csv"), row.names = FALSE)

# join to interactome separation (S_AB), which is stored one row per unordered pair
sab <- read.csv(SAB_FILE, stringsAsFactors = FALSE)
names(sab)[names(sab) == "Network_separation_s_AB"] <- "Network_separation_s_AB"
key <- function(a, b) apply(cbind(a, b), 1, function(v) paste(sort(v), collapse = "-"))
sab$key     <- key(sab$from, sab$to)
overlap$key <- key(overlap$A, overlap$B)

merged <- merge(sab[, c("key", "from", "to", "Network_separation_s_AB")],
                overlap[, c("key", "overlap_n")], by = "key")
merged <- merged[order(merged$Network_separation_s_AB), ]
write.csv(merged, file.path(IN_DIR, "tier2_sab_vs_overlap.csv"), row.names = FALSE)

ct <- suppressWarnings(cor.test(merged$Network_separation_s_AB, merged$overlap_n, method = "spearman"))
message(sprintf("Tier 2 complete: %d pairs; Spearman rho = %.3f, P = %.3f",
                nrow(merged), ct$estimate, ct$p.value))
