# Figure 5: drug-repurposing tiers (not in the original figure-generation skill
# spec, added after the Tier 1/2/3 signatureSearch drug-repurposing analysis).
# Panel A: cross-skin-disease drug-overlap heatmap (Tier 2), same visual
#          convention as Figure 1's S_AB heatmap, for direct comparison.
# Panel B: top connected-disease <-> skin-disease drug overlaps (Tier 3),
#          i.e. drugs found for a skin disease that are independently also
#          found for one of its DisGeNET+TriNetX-corroborated connected diseases.

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(patchwork))

PALETTE <- c(PS = "#E64B35", AD = "#4DBBD5", HS = "#00A087",
             VL = "#3C5488", CSU = "#F39B7F", SLS = "#8491B4")
DISEASE_ORDER <- c("HS", "PS", "VL", "AD", "CSU", "SLS")

## Panel A: 6x6 drug-overlap heatmap (Tier 2)
ov <- read.csv("signaturesearch/tier2_pairwise_overlap.csv")
pairs <- expand.grid(from = DISEASE_ORDER, to = DISEASE_ORDER, stringsAsFactors = FALSE)
lookup_n <- function(a, b) {
  hit <- ov$overlap_n[(ov$A == a & ov$B == b) | (ov$A == b & ov$B == a)]
  if (length(hit) == 0) return(NA_integer_)
  hit[1]
}
pairs$n <- mapply(lookup_n, pairs$from, pairs$to)
pairs$n[pairs$from == pairs$to] <- NA
pairs$from <- factor(pairs$from, levels = DISEASE_ORDER)
pairs$to <- factor(pairs$to, levels = rev(DISEASE_ORDER))

panelA <- ggplot(pairs, aes(x = from, y = to, fill = n)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = ifelse(is.na(n), "", n)), size = 3.2, color = "black") +
  scale_fill_gradient(low = "white", high = "#2166AC", na.value = "grey90", name = "Shared\ntop-50 drugs") +
  coord_fixed() +
  labs(title = "A. Cross-disease drug-repurposing overlap",
       subtitle = "Number of shared compounds among each pair's top 50 LINCS reversal hits",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(color = PALETTE[levels(pairs$from)], face = "bold", size = 12),
    axis.text.y = element_text(color = PALETTE[rev(levels(pairs$to))], face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 8.5, color = "grey30")
  )

## Panel B: top Tier-3 connected-disease <-> skin-disease overlaps
t3 <- read.csv("signaturesearch/tier3_overlap_results.csv")
topB <- t3 %>% arrange(desc(overlap_n)) %>% slice_head(n = 10)
topB$label <- paste0(topB$connected_disease, " x ", topB$skin_disease)
topB$label <- factor(topB$label, levels = rev(topB$label))
# show the top 3 drug names per bar as a compact annotation
topB$drug_preview <- sapply(strsplit(topB$drugs, ";"), function(x) paste(head(x, 3), collapse = ", "))
topB$bar_label <- paste0(topB$overlap_n, "  (", topB$drug_preview, ", ...)")

panelB <- ggplot(topB, aes(x = overlap_n, y = label, fill = skin_disease)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = bar_label), hjust = -0.03, size = 2.9) +
  scale_fill_manual(values = PALETTE, name = "Skin disease") +
  scale_x_continuous(expand = expansion(mult = c(0, 1.1))) +
  labs(title = "B. Top connected-disease vs. skin-disease drug overlaps",
       subtitle = "Diseases both DisGeNET-predicted and TriNetX-corroborated as related to a skin disease,\nindependently sharing top-100 drug hits with that skin disease's LINCS signature (top 3 shared drugs shown)",
       x = "Shared drugs (of top 100 each)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 8, color = "grey30"),
    legend.position = "bottom"
  )

combined <- panelA + panelB + plot_layout(widths = c(0.8, 1.6))

dir.create("images/high_res", showWarnings = FALSE, recursive = TRUE)
ggsave("images/high_res/drug_repurposing_tiers.png", combined, width = 17, height = 6.5, dpi = 300, bg = "white")
ggsave("images/high_res/drug_repurposing_tiers.tiff", combined, width = 17, height = 6.5, dpi = 300, compression = "lzw", bg = "white")
message("Saved images/high_res/drug_repurposing_tiers.png and .tiff")
