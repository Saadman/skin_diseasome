# Figure 3: S_AB network-separation heatmap.
# Source: networksim/sab_values.csv -- the August 2022 separation.py run, computed from
# the differential-expression gene sets (see networksim/README.md).
# Sign convention (Menche et al. 2015): S_AB < 0 = disease modules closer than
# random expectation (i.e. more related); S_AB > 0 = more separated than random.

suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))

PALETTE <- c(PS = "#E64B35", AD = "#4DBBD5", HS = "#00A087",
             VL = "#3C5488", CSU = "#F39B7F", SLS = "#8491B4")
DISEASE_ORDER <- c("HS", "PS", "VL", "AD", "CSU", "SLS")

sab <- read.csv("networksim/sab_values.csv", stringsAsFactors = FALSE)
sab <- sab[, c("from", "to", "Network_separation_s_AB")]

# build symmetric matrix
pairs <- expand.grid(from = DISEASE_ORDER, to = DISEASE_ORDER, stringsAsFactors = FALSE)
lookup <- function(a, b) {
  hit <- sab$Network_separation_s_AB[(sab$from == a & sab$to == b) | (sab$from == b & sab$to == a)]
  if (length(hit) == 0) return(NA_real_)
  hit[1]
}
pairs$s_ab <- mapply(lookup, pairs$from, pairs$to)
pairs$s_ab[pairs$from == pairs$to] <- NA

pairs$from <- factor(pairs$from, levels = DISEASE_ORDER)
pairs$to <- factor(pairs$to, levels = rev(DISEASE_ORDER))

p <- ggplot(pairs, aes(x = from, y = to, fill = s_ab)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = ifelse(is.na(s_ab), "", sprintf("%.3f", s_ab))),
            size = 3.2, color = "black") +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0,
                        na.value = "grey90", name = expression(S[AB])) +
  coord_fixed() +
  labs(title = expression("Network separation (" * S[AB] * ") among the six skin diseases"),
       subtitle = expression(S[AB] * " < 0: modules closer than random (related)  |  " * S[AB] * " > 0: more separated than random"),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(color = PALETTE[levels(pairs$from)], face = "bold", size = 12),
    axis.text.y = element_text(color = PALETTE[levels(pairs$to)], face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "grey30")
  )

dir.create("images/high_res", showWarnings = FALSE, recursive = TRUE)
ggsave("images/high_res/sab_heatmap.png", p, width = 6.5, height = 5.5, dpi = 300, bg = "white")
ggsave("images/high_res/sab_heatmap.tiff", p, width = 6.5, height = 5.5, dpi = 300, compression = "lzw", bg = "white")
message("Saved images/high_res/sab_heatmap.png and .tiff")
