# Figure 4 (static, print-ready): the diseasome network for the manuscript,
# styled after Goh et al. 2007's classic "human disease network" — comorbidity
# nodes colored by disorder class (not a flat grey), skin diseases as labeled hubs.
#
# Two edge types, drawn distinctly:
#  (a) disease-disease "backbone" edges among the 6 skin diseases themselves,
#      weighted by S_AB network separation (networksim/sab_values.csv, the same
#      data as the Figure 3 heatmap) — drawn only where S_AB < 0 (modules closer
#      than random). SLS-ICD therefore has no backbone edge, and CSU connects
#      only to AD and VL; both are real results, not missing data.
#  (b) predicted-comorbidity edges from each skin disease to its top-10 DisGeNET
#      associations (disgenet/diseasome_final_network_connections.xlsx, 2022 dataset),
#      colored by the comorbid disease's disorder class.

suppressMessages(library(readxl))
suppressMessages(library(dplyr))
suppressMessages(library(igraph))
suppressMessages(library(ggraph))
suppressMessages(library(ggplot2))

PALETTE <- c(PS = "#E64B35", AD = "#4DBBD5", HS = "#00A087",
             VL = "#3C5488", CSU = "#F39B7F", SLS = "#8491B4")
SRC_TO_ABB <- c("Psoriasis"="PS","Atopic Dermatitis"="AD","Hidradenitis Suppurativa"="HS",
                "Vitiligo"="VL","Chronic Urticaria"="CSU","SLS"="SLS")
FULL_NAME <- c(PS="Psoriasis", AD="Atopic Dermatitis", HS="Hidradenitis Suppurativa",
               VL="Vitiligo", CSU="Chronic Urticaria", SLS="SLS")

CLASS_PALETTE <- c(
  "Infections" = "#D95F02",
  "Neoplasms" = "#7570B3",
  "Nervous System Diseases" = "#1F78B4",
  "Musculoskeletal Diseases" = "#66A61E",
  "Cardiovascular Diseases" = "#E7298A",
  "Female Urogenital Diseases and Pregnancy Complications" = "#E6AB02",
  "Skin and Connective Tissue Diseases" = "#A6761D",
  "Digestive System Diseases" = "#B15928",
  "Hemic and Lymphatic Diseases" = "#CAB2D6",
  "Nutritional and Metabolic Diseases" = "#1B9E77",
  "Pathological Conditions, Signs and Symptoms" = "#B2DF8A",
  "Other / unclassified" = "#999999"
)

TOP_N <- 10  # per disease; kept small so every target label stays legible

## (a) disease-disease backbone from S_AB
sab <- read.csv("networksim/sab_values.csv", stringsAsFactors = FALSE)
backbone <- sab %>%
  transmute(source = FULL_NAME[from], target = FULL_NAME[to], weight = -Network_separation_s_AB,
            etype = "backbone")
backbone <- backbone[backbone$weight > 0, ]  # only genuinely related pairs (S_AB<0)

## (b) predicted comorbidities from DisGeNET, with primary disorder class
d <- read_excel("disgenet/diseasome_final_network_connections.xlsx")
d$abb <- SRC_TO_ABB[d$source]
d$primary_class <- sapply(strsplit(d$target_class, ";"), function(x) trimws(x[1]))
d$primary_class[is.na(d$primary_class) | d$primary_class == "NA"] <- "Other / unclassified"
top <- d %>% group_by(source) %>% arrange(FDR) %>% slice_head(n = TOP_N) %>% ungroup()
comorbid <- top %>% transmute(source = source, target = target, weight = -log10(pmax(FDR, 1e-30)), etype = "comorbidity")

edges <- bind_rows(backbone, comorbid)
edges$edge_col <- ifelse(edges$etype == "backbone", "#08306B", "#B0ABA0")
# backbone width scaled linearly across its own actual weight range (0.018-0.594),
# NOT a fixed floor -- otherwise all but the strongest edge collapse to one width.
bb_range <- range(backbone$weight)
edges$edge_w <- ifelse(
  edges$etype == "backbone",
  0.6 + (edges$weight - bb_range[1]) / diff(bb_range) * (3.6 - 0.6),
  0.35
)
edges$edge_a <- ifelse(edges$etype == "backbone", 0.75, 0.55)
edges$etype <- factor(edges$etype, levels = c("backbone", "comorbidity"),
                       labels = c("Disease-disease relatedness (S_AB)", "Predicted comorbidity (DisGeNET)"))

nodes <- data.frame(name = unique(c(edges$source, edges$target)), stringsAsFactors = FALSE)
nodes$is_source <- nodes$name %in% FULL_NAME
nodes$abb <- SRC_TO_ABB[nodes$name]
nodes$node_col <- ifelse(nodes$is_source, PALETTE[nodes$abb], NA)  # identity color for sources only
class_lookup <- top$primary_class[match(nodes$name, top$target)]
nodes$primary_class <- ifelse(nodes$is_source, NA, class_lookup)
degree_tab <- table(comorbid$target)
nodes$degree <- ifelse(nodes$is_source, NA, degree_tab[nodes$name])
nodes$size <- ifelse(nodes$is_source, 9, pmin(2.6 + ifelse(is.na(nodes$degree), 1, nodes$degree) * 0.8, 5.5))

g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

# The force-directed layout places disconnected components (SLS has zero edges
# to the rest) at an arbitrary, meaningless distance apart. Compute the layout,
# then rigidly translate each disconnected component closer to the main one's
# centroid -- internal relative positions within a component are untouched, only
# the whole component is moved, purely for legible display.
set.seed(42)
lay <- create_layout(g, layout = "fr", niter = 3000)
comp <- components(g)$membership
main_comp <- which.max(table(comp))
main_centroid <- colMeans(lay[comp[lay$name] == main_comp, c("x", "y")])
for (cid in unique(comp)) {
  if (cid == main_comp) next
  idx <- comp[lay$name] == cid
  this_centroid <- colMeans(lay[idx, c("x", "y")])
  # pull the component to within a fixed radius of the main centroid
  target_radius <- 4.2
  dir_vec <- this_centroid - main_centroid
  dir_len <- sqrt(sum(dir_vec^2))
  if (dir_len == 0) dir_vec <- c(1, 0) else dir_vec <- dir_vec / dir_len
  new_centroid <- main_centroid + dir_vec * target_radius
  shift <- new_centroid - this_centroid
  lay$x[idx] <- lay$x[idx] + shift[1]
  lay$y[idx] <- lay$y[idx] + shift[2]
}

p <- ggraph(lay) +
  geom_edge_link(aes(edge_width = edge_w, edge_alpha = edge_a, edge_color = edge_col,
                      edge_linetype = etype),
                 show.legend = TRUE) +
  scale_edge_width_identity() +
  scale_edge_alpha_identity() +
  scale_edge_color_identity() +
  scale_edge_linetype_manual(values = c("Disease-disease relatedness (S_AB)" = "solid",
                                         "Predicted comorbidity (DisGeNET)" = "22"),
                              name = "Connection type") +
  guides(edge_linetype = guide_legend(override.aes = list(edge_width = c(2.1, 0.35),
                                                           edge_color = c("#08306B", "#B0ABA0")))) +
  # target nodes: colored by disorder class, with a legend
  geom_node_point(data = function(x) dplyr::filter(x, !is_source),
                   aes(size = size, color = primary_class)) +
  scale_color_manual(values = CLASS_PALETTE, name = "Disorder class\n(comorbidity nodes)", na.translate = FALSE) +
  ggnewscale::new_scale_color() +
  # source nodes: fixed disease palette, drawn on top, identity color
  geom_node_point(data = function(x) dplyr::filter(x, is_source),
                   aes(size = size, color = node_col)) +
  scale_color_identity() +
  geom_node_text(aes(label = ifelse(is_source, abb, NA)), fontface = "bold", size = 3.6, color = "white") +
  geom_node_text(aes(label = ifelse(!is_source, name, NA)), repel = TRUE, size = 2.3, color = "grey20",
                 max.overlaps = 100, segment.size = 0.15, segment.color = "grey70",
                 bg.color = "white", bg.r = 0.15) +
  scale_size_identity() +
  labs(title = "Skin diseasome — network of disease relationships and predicted comorbidities",
       subtitle = paste0("Solid lines: disease-disease relatedness among the 6 skin diseases (network separation S_AB < 0, Figure 3).\n",
                          "Dashed lines: top ", TOP_N, " DisGeNET-predicted comorbidities per disease, colored by disorder class.")) +
  guides(size = "none") +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 8, color = "grey30", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    plot.margin = margin(10, 10, 10, 10)
  )

dir.create("images/high_res", showWarnings = FALSE, recursive = TRUE)
ggsave("images/high_res/diseasome_network_static.png", p, width = 13, height = 9.5, dpi = 300, bg = "white")
ggsave("images/high_res/diseasome_network_static.tiff", p, width = 13, height = 9.5, dpi = 300, compression = "lzw", bg = "white")
message("Saved images/high_res/diseasome_network_static.png and .tiff")
