# =================== Load Required Packages =====================
library(dplyr)
library(readr)
library(pheatmap)
library(stringr)
library(RColorBrewer)

# =================== Load Your Data =====================
abundance <- read_csv("coverm_431_MAGS_metagenomes_reads_count.csv")
taxonomy <- read_delim("MAG_taxonomy.tsv", delim = "\t", show_col_types = FALSE)
sample_meta <- read_csv("Samples-mendota.csv")

# =================== Normalize Abundance =====================
abundance_norm <- abundance
abundance_norm[ , -1] <- sweep(abundance[ , -1], 2, colSums(abundance[ , -1], na.rm = TRUE), FUN = "/")

# =================== Grouping Samples =====================
group_samples <- function(var) {
  split(sample_meta$Sample.Name, sample_meta[[var]])
}

depth_groups <- group_samples("Depth")
month_groups <- group_samples("Month")
oxygen_groups <- group_samples("Oxygen")

# =================== Top 10 MAGs per Group =====================
get_top_10_MAGs <- function(group_list, abundance_data) {
  result <- list()
  for (grp in names(group_list)) {
    samples <- group_list[[grp]]
    if (length(samples) == 0) next
    top_mags <- abundance_data %>%
      select(Genome, all_of(samples)) %>%
      mutate(avg = rowMeans(across(-Genome), na.rm = TRUE)) %>%
      arrange(desc(avg)) %>%
      slice_head(n = 10) %>%
      pull(Genome)
    result[[grp]] <- top_mags
  }
  return(result)
}

top_depth <- get_top_10_MAGs(depth_groups, abundance_norm)
top_month <- get_top_10_MAGs(month_groups, abundance_norm)
top_oxygen <- get_top_10_MAGs(oxygen_groups, abundance_norm)

# =================== Build Presence Matrix =====================
all_conditions <- c(names(top_depth), names(top_month), names(top_oxygen))
all_mags <- unique(unlist(c(top_depth, top_month, top_oxygen)))

presence_matrix <- matrix(0, nrow = length(all_mags), ncol = length(all_conditions),
                          dimnames = list(all_mags, all_conditions))

# Fill matrix
for (cond in names(top_depth)) {
  presence_matrix[top_depth[[cond]], cond] <- 1
}
for (cond in names(top_month)) {
  presence_matrix[top_month[[cond]], cond] <- 1
}
for (cond in names(top_oxygen)) {
  presence_matrix[top_oxygen[[cond]], cond] <- 1
}

# =================== Map Phylum to MAGs =====================
phylum_map <- taxonomy %>%
  filter(Genome %in% rownames(presence_matrix)) %>%
  mutate(Phylum = str_remove(Phylum, "p__")) %>%
  select(Genome, Phylum)

# Create named vector for row annotation
row_annotation <- data.frame(Phylum = setNames(phylum_map$Phylum, phylum_map$Genome)[rownames(presence_matrix)])
rownames(row_annotation) <- rownames(presence_matrix)

# Create color palette for phylum
phylum_levels <- unique(row_annotation$Phylum)
phylum_colors <- setNames(colorRampPalette(brewer.pal(8, "Set2"))(length(phylum_levels)), phylum_levels)
ann_colors <- list(Phylum = phylum_colors)

# =================== Plot Heatmap =====================
pheatmap(presence_matrix,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         annotation_row = row_annotation,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("white", "blue"))(100),
         fontsize_row = 8,
         fontsize_col = 10,
         main = "Presence of Top 10 MAGs Across Environmental Conditions")

write.csv(presence_matrix ,"presence_matrix.csv")
