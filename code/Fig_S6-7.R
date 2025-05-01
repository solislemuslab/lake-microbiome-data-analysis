rm(list=ls())
set.seed(1234)
# Load required libraries
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# === Load data ===
abundance <- read_csv("coverm_431_MAGS_metagenomes_reads_count.csv")
samples <- read_csv("Samples-mendota.csv")

# === Transpose abundance: rows = samples ===
abundance_t <- abundance %>%
  pivot_longer(-Genome, names_to = "Sample.Name", values_to = "Count") %>%
  pivot_wider(names_from = Genome, values_from = Count)

# === Merge with metadata ===
merged_df <- left_join(samples, abundance_t, by = "Sample.Name")

# === Extract sample metadata and labels ===
merged_df <- merged_df %>%
  mutate(
    Month_Num = str_sub(Sample.Name, 6, 7),
    Day = str_sub(Sample.Name, 9, 10),
    Month = case_when(
      Month_Num == "07" ~ "July",
      Month_Num == "08" ~ "August",
      Month_Num == "09" ~ "September",
      Month_Num == "10" ~ "October",
      TRUE ~ Month_Num
    ),
    Month = factor(Month, levels = c("July", "August", "September", "October")),
    Depth = as.character(Depth),
    Sample_Label = paste(Month, Day, paste0(Depth, "m"), Oxygen, sep = "-")
  )

# === Get correct sample order ===
sample_order <- merged_df %>%
  mutate(Depth_num = as.numeric(Depth)) %>%
  arrange(Month, Depth_num) %>%
  pull(Sample_Label)

# === Prepare abundance matrix for correlation ===
mag_cols <- setdiff(colnames(merged_df), colnames(samples))

# abundance_matrix <- merged_df %>%
#   select(Sample_Label, all_of(mag_cols)) %>%
#   column_to_rownames("Sample_Label") %>%
#   mutate(across(everything(), as.numeric)) %>%  # <- this ensures numeric columns
#   as.matrix()
# Normalize with Total Sum Scaling (TSS)
abundance_matrix <- merged_df %>%
  select(Sample_Label, all_of(mag_cols)) %>%
  column_to_rownames("Sample_Label") %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

# Apply Total Sum Scaling row-wise
abundance_matrix <- abundance_matrix / rowSums(abundance_matrix)
rowSums(abundance_matrix )



# === Correlation matrices (samples in rows) ===
cor_pearson <- cor(t(abundance_matrix), method = "pearson")
cor_spearman <- cor(t(abundance_matrix), method = "spearman")

# === Reorder matrices by sample_order ===
cor_pearson <- cor_pearson[sample_order, sample_order]
cor_spearman <- cor_spearman[sample_order, sample_order]

# Convert correlation matrix to long format and filter
high_corr_pairs <- which(cor_pearson > 0.8 & cor_pearson < 1, arr.ind = TRUE)
# Create a data frame of results
high_corr_df <- data.frame(
  Sample1 = rownames(cor_pearson)[high_corr_pairs[, 1]],
  Sample2 = colnames(cor_pearson)[high_corr_pairs[, 2]],
  Correlation = cor_pearson[high_corr_pairs]
)
# Remove duplicate pairs (e.g., A-B and B-A)
high_corr_df <- high_corr_df[as.character(high_corr_df$Sample1) < as.character(high_corr_df$Sample2), ]
print(high_corr_df)

# Convert correlation matrix to long format and filter
high_corr_spearman_pairs <- which(cor_spearman > 0.8 & cor_spearman < 1, arr.ind = TRUE)
# Create a data frame of results
high_corr_spearman_df <- data.frame(
  Sample1 = rownames(cor_spearman)[high_corr_spearman_pairs[, 1]],
  Sample2 = colnames(cor_spearman)[high_corr_spearman_pairs[, 2]],
  Correlation_spearman = cor_spearman[high_corr_spearman_pairs]
)
# Remove duplicate pairs (e.g., A-B and B-A)
high_corr_spearman_df <- high_corr_spearman_df[as.character(high_corr_spearman_df$Sample1) < as.character(high_corr_spearman_df$Sample2), ]
print(high_corr_spearman_df)




# Create identifier columns for matching (ignoring correlation values)
high_corr_df$pair_id <- paste(high_corr_df$Sample1, high_corr_df$Sample2, sep = "___")
high_corr_spearman_df$pair_id <- paste(high_corr_spearman_df$Sample1, high_corr_spearman_df$Sample2, sep = "___")

# Find intersecting pair IDs
common_ids <- intersect(high_corr_df$pair_id, high_corr_spearman_df$pair_id)

# Subset original data frames to keep only the intersecting pairs
intersect_df <- high_corr_df[high_corr_df$pair_id %in% common_ids, ]
intersect_df <- merge(intersect_df,
                      high_corr_spearman_df[, c("Sample1", "Sample2", "Correlation_spearman", "pair_id")],
                      by = "pair_id")

# Clean up: remove the 'pair_id' column if you don’t want it in the final result
intersect_df <- intersect_df[, c("Sample1.x", "Sample2.x", "Correlation", "Correlation_spearman")]
colnames(intersect_df) <- c("Sample1", "Sample2", "Pearson", "Spearman")

# View result
print(intersect_df)


# === Plot and save heatmaps ===
pdf("correlation_Pearson.pdf", width = 11, height = 9)

# Pearson
Heatmap(cor_pearson,
        name = "Pearson",
        row_order = sample_order,
        column_order = sample_order,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_names_rot = 90,
        row_names_rot = 0,
        column_names_gp = gpar(fontsize = 9),
        row_names_gp = gpar(fontsize = 9),
        heatmap_legend_param = list(title = "Pearson"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", cor_pearson[i, j]), x, y, gp = gpar(fontsize = 6))
        },
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
)
#dev.off()

pdf("correlation_Spearman.pdf", width = 11, height = 9)
Heatmap(cor_spearman,
        name = "Spearman",
        row_order = sample_order,
        column_order = sample_order,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_names_rot = 90,
        row_names_rot = 0,
        column_names_gp = gpar(fontsize = 9),
        row_names_gp = gpar(fontsize = 9),
        heatmap_legend_param = list(title = "Spearman"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", cor_spearman[i, j]), x, y, gp = gpar(fontsize = 6))
        },
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
)

dev.off()














# === Create binary matrices for Pearson, Spearman, and intersection ===

# Helper: zero diagonal
zero_diag <- function(mat) {
  mat[upper.tri(mat, diag = TRUE)] <- 0
  return(mat)
}

# Binary Pearson: 1 if r > 0.8 and < 1
bin_pearson <- (cor_pearson > 0.8 & cor_pearson < 1) * 1
bin_pearson <- zero_diag(bin_pearson)

# Binary Spearman
bin_spearman <- (cor_spearman > 0.8 & cor_spearman < 1) * 1
bin_spearman <- zero_diag(bin_spearman)

# Binary Intersection: both Pearson & Spearman > 0.8
bin_intersection <- ((bin_pearson == 1) & (bin_spearman == 1)) * 1
bin_intersection <- zero_diag(bin_intersection)

# === Plot binary heatmaps ===
pdf("FigS6-7.pdf", width = 10, height = 10)

# Pearson binary
Heatmap(bin_pearson,
        name = "Pearson > 0.8",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_order = sample_order,
        column_order = sample_order,
        column_names_rot = 90,
        row_names_rot = 0,
        column_names_gp = gpar(fontsize = 9),
        row_names_gp = gpar(fontsize = 9),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", cor_spearman[i, j]), x, y, gp = gpar(fontsize = 6))
        },      col = c("0" = "white", "1" = "black")
        
)

# Spearman binary
Heatmap(bin_spearman,
        name = "Spearman > 0.8",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_order = sample_order,
        column_order = sample_order,
        column_names_rot = 90,
        row_names_rot = 0,
        column_names_gp = gpar(fontsize = 9),
        row_names_gp = gpar(fontsize = 9),
        col = c("0" = "white", "1" = "black")
)

# Intersection
Heatmap(bin_intersection,
        name = "Pearson & Spearman > 0.8",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_order = sample_order,
        column_order = sample_order,
        column_names_rot = 90,
        row_names_rot = 0,
        column_names_gp = gpar(fontsize = 9),
        row_names_gp = gpar(fontsize = 9),
        col = c("0" = "white", "1" = "black")
)

dev.off()

