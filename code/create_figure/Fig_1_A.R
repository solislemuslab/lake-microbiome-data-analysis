rm(list=ls())
# === Required Libraries ===
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(patchwork)  # For combining plots
library(gridExtra)
library(grid)
# === Load Data ===
presence_matrix <- read_csv("presence_matrix.csv")
taxonomy <- read_delim("MAG_taxonomy.tsv", delim = "\t", show_col_types = FALSE)

# === Prepare presence & percentage info ===
presence_df <- presence_matrix %>%
  mutate(Presence = rowSums(across(-MAG))) %>%
  mutate(Percent = Presence / 11 * 100) %>%
  left_join(taxonomy %>% select(Genome, Phylum), by = c("MAG" = "Genome")) %>%
  mutate(Phylum = str_remove(Phylum, "p__")) %>%
  arrange(Phylum, desc(Percent))

# === Reorder matrix ===
presence_mat <- presence_matrix %>%
  column_to_rownames("MAG") %>%
  .[presence_df$MAG, ]

# === Annotations ===
annotation_row <- presence_df %>%
  select(MAG, Phylum) %>%
  column_to_rownames("MAG")

# === Define colors ===
phylum_colors <- c(
  "Proteobacteria"      = "#D8BFD8",
  "Desulfobacterota"    = "#F4C28C",
  "Bacteroidota"        = "#9491D9",
  "Myxococcota"         = "#82C341",
  "Verrucomicrobiota"   = "#A2C4EC",
  "Actinobacteriota"    = "#B6D7A8",
  "Planctomycetota"     = "#A94442",
  "Cyanobacteria"       = "#F4CCCC",
  "Krumholzibacteriota" = "#D1EEEE",
  "Chloroflexota"       = "#DDCC77",
  "Firmicutes_A"        = "#44AA99"
)

annotation_colors <- list(Phylum = phylum_colors)

# === Bar chart for % presence ===
bar_df <- presence_df %>%
  select(MAG, Percent, Phylum)

bar_plot <- ggplot(bar_df, aes(x = Percent, y = fct_rev(factor(MAG, levels = presence_df$MAG)), fill = Phylum)) +
  geom_col(width = 0.4) +  # 🔻 Thinner bars
  scale_fill_manual(values = phylum_colors) +
  theme_minimal() +
  labs(x = "Presence in the top 10 ranked taxa across environmental conditions", y = NULL) +
  theme(
    legend.position = "right",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(17, 11, 30, 7)  # 🔻 Tighter spacing
  )



# === Heatmap ===
pheatmap_out <- pheatmap(as.matrix(presence_mat),
                         cluster_rows = FALSE,
                         cluster_cols = FALSE,
                         color = c("white", "gray"),
                         annotation_row = annotation_row,
                         annotation_colors = annotation_colors,
                         labels_row = rownames(presence_mat),
                         show_rownames = TRUE,
                         fontsize_row = 8,
                         fontsize_col = 10,
                         angle_col = 90,
                         legend = FALSE,
                         annotation_legend = FALSE,
                         border_color = "grey80",
                         gaps_col = c(4, 8),
                         main = "")

# === Combine plots ===
heatmap_grob <- grid::grid.grabExpr(pheatmap_out$gtable)

# === Combine barplot and heatmap (side-by-side) ===
pdf("Fig_1_A.pdf", , width = 14, height = 10)
grid.arrange(
  pheatmap_out$gtable,
  ggplotGrob(bar_plot),
  ncol = 2,
  widths = c(4.5, 2)  # Adjust for spacing between bar and heatmap
)
dev.off()
