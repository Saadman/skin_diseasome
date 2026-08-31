# Tier 3: for each DisGeNET-predicted disease with an observed TriNetX co-occurrence
# (see trinetx/genetic_rwe_fdr_relation.xlsx), pull its DisGeNET
# gene list directly (via api.disgenet.com REST) and run gess_fisher (unsigned
# gene-set enrichment) against the same `lincs` reference used for Tier 1.
# Requires DISGENET_API_KEY set in the environment.

suppressMessages(library(httr))
suppressMessages(library(jsonlite))
suppressMessages(library(signatureSearch))
suppressMessages(library(signatureSearchData))

KEY <- Sys.getenv("DISGENET_API_KEY")
if (KEY == "") stop("Set DISGENET_API_KEY in your environment before running this script.")

OUT_DIR <- "signaturesearch/tier3"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

lookup <- read.csv("signaturesearch/tier3_disease_list.csv")
lookup <- lookup[!is.na(lookup$CUI), ]
cat("Diseases to process:", nrow(lookup), "\n")

get_disease_genes <- function(cui) {
  all_rows <- list(); page <- 0
  repeat {
    resp <- NULL
    for (attempt in 1:3) {
      resp <- tryCatch(
        httr::GET(paste0("https://api.disgenet.com/api/v1/gda/summary?disease=UMLS_", cui,
                   "&source=CURATED&page_number=", page),
            httr::add_headers(Authorization = KEY), httr::timeout(60)),
        error = function(e) NULL)
      if (!is.null(resp)) break
      Sys.sleep(5 * attempt)
    }
    if (is.null(resp)) return(character(0))
    parsed <- tryCatch(fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"), flatten = TRUE),
                        error = function(e) NULL)
    if (is.null(parsed)) break
    payload <- parsed$payload
    if (is.null(payload)) break
    if (!is.data.frame(payload)) payload <- tryCatch(as.data.frame(payload, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(payload) || nrow(payload) == 0) break
    if (!all(c("symbolOfGene", "geneNcbiID") %in% colnames(payload))) break
    all_rows[[length(all_rows) + 1]] <- payload[, c("symbolOfGene", "geneNcbiID")]
    total <- parsed$paging$totalElements
    got <- sum(vapply(all_rows, nrow, integer(1)))
    if (is.null(total) || length(total) == 0 || got >= total) break
    page <- page + 1
  }
  if (length(all_rows) == 0) return(character(0))
  df <- do.call(rbind, all_rows)
  unique(as.character(df$geneNcbiID))
}

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

for (i in seq_len(nrow(lookup))) {
  disease <- lookup$disease[i]
  cui <- lookup$CUI[i]
  fname <- file.path(OUT_DIR, paste0(safe_name(disease), "_fisher_result.rds"))
  if (file.exists(fname)) { message("skip (exists): ", disease); next }
  process_one <- function() {
    message("=== ", disease, " (", cui, ") ===")
    genes <- get_disease_genes(cui)
    message("  ", length(genes), " genes")
    if (length(genes) < 5) { message("  too few genes, skipping"); return(invisible(NULL)) }
    qsig <- qSig(query = list(upset = genes, downset = character(0)), gess_method = "Fisher", refdb = "lincs")
    res <- gess_fisher(qsig, higher = 2, lower = -2, workers = 2)
    df <- result(res)
    saveRDS(list(disease = disease, cui = cui, n_genes = length(genes), result = df), fname)
    message("  done: ", nrow(df), " rows")
  }
  tryCatch(process_one(), error = function(e) message("  SKIPPED (error): ", conditionMessage(e)))
}
message("Tier 3 complete.")
