# Figure 3: GO biological-process dot plots (figure-generation skill).
# Source: go_enrichment/disease_GO_enrichment_BiologicalProcess.csv
# (the correct 6-disease GO table for this manuscript; NOT the old EAACI/ 4-disease file).

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

PALETTE <- c(PS = "#E64B35", AD = "#4DBBD5", HS = "#00A087",
             VL = "#3C5488", CSU = "#F39B7F", SLS = "#8491B4")
CLUSTER_TO_ABB <- c(
  Psoriasis = "PS", Atopic_Dermatitis = "AD", Hidradenitis = "HS",
  Vitiligo = "VL", Chronic_Urticaria = "CSU", SLS = "SLS"
)
DISEASE_LABEL <- c(PS = "Psoriasis", AD = "Atopic Dermatitis", HS = "Hidradenitis Suppurativa",
                    VL = "Vitiligo", CSU = "Chronic Urticaria", SLS = "SLS")

go <- read.csv("go_enrichment/disease_GO_enrichment_BiologicalProcess.csv")
go$abb <- CLUSTER_TO_ABB[go$Cluster]
go$abb <- factor(go$abb, levels = c("PS", "AD", "HS", "VL", "CSU", "SLS"))

# parse GeneRatio "211/5313" -> 211/5313
parse_ratio <- function(x) {
  parts <- strsplit(x, "/")
  sapply(parts, function(p) as.numeric(p[1]) / as.numeric(p[2]))
}
go$gene_ratio <- parse_ratio(go$GeneRatio)

TOP_N <- 10
top <- go %>%
  group_by(abb) %>%
  arrange(p.adjust) %>%
  slice_head(n = TOP_N) %>%
  ungroup()

# order terms within each disease facet by gene_ratio for a clean look
top$Description <- factor(top$Description, levels = rev(unique(top$Description[order(top$abb, top$gene_ratio)])))

p <- ggplot(top, aes(x = gene_ratio, y = Description, size = Count, color = p.adjust)) +
  geom_point() +
  facet_wrap(~ abb, scales = "free_y", ncol = 2,
             labeller = labeller(abb = DISEASE_LABEL)) +
  scale_color_gradient(low = "#B2182B", high = "#2166AC", trans = "log10", name = "adj. p-value") +
  scale_size_continuous(name = "Gene count", range = c(2, 7)) +
  labs(title = "Top enriched GO biological processes per skin disease",
       x = "Gene ratio", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold", size = 13),
    panel.grid.minor = element_blank()
  )

dir.create("images/high_res", showWarnings = FALSE, recursive = TRUE)
ggsave("images/high_res/go_dotplot.png", p, width = 11, height = 10, dpi = 300, bg = "white")
ggsave("images/high_res/go_dotplot.tiff", p, width = 11, height = 10, dpi = 300, compression = "lzw", bg = "white")
message("Saved images/high_res/go_dotplot.png and .tiff")
