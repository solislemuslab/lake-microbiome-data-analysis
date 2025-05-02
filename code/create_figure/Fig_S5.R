rm(list=ls())
# Load required libraries
library(vegan)
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
# Load abundance data
df <- read.csv("coverm_431_MAGS_metagenomes_reads_count.csv", row.names = 1, check.names = FALSE)
df_t <- as.data.frame(t(df))  # Transpose: samples as rows

# Create Sample column from row names
df_t$Sample <- rownames(df_t)

# Clean up sample names
df_t$Sample <- gsub("^X", "", df_t$Sample)

# Load sample metadata
metadata <- read_csv("Samples-mendota.csv")

# Merge metadata into the abundance data
df_merged <- merge(df_t, metadata[, c("Sample.Name", "Month", "Depth", "Oxygen")],
                   by.x = "Sample", by.y = "Sample.Name")

# Reorder months and convert to factors
df_merged$Month <- factor(df_merged$Month, levels = c("July", "August", "September", "October"))
df_merged$Depth <- as.factor(df_merged$Depth)
df_merged$Oxygen <- as.factor(df_merged$Oxygen)

# Separate abundance data
mag_columns <- setdiff(colnames(df_merged), c("Sample", "Month", "Depth", "Oxygen"))
abundance_data <- df_merged[, mag_columns]

# Convert to numeric
abundance_data <- as.data.frame(lapply(abundance_data, as.numeric))
abundance_data[is.na(abundance_data)] <- 0  # Replace NAs with 0

####IMP change the resulttttttt in plot 
#abundance_data <- abundance_data[order(rownames(abundance_data)), ]

# Hellinger transformation
hellinger_transform <- decostand(abundance_data, method = "hellinger")

# Bray-Curtis dissimilarity
bray_dist <- vegdist(hellinger_transform, method = "bray")

# PCoA
pcoa_result <- cmdscale(bray_dist, eig = TRUE, k = 2)
df_merged$PCoA1 <- pcoa_result$points[,1]
df_merged$PCoA2 <- pcoa_result$points[,2]

# Variance explained
eigenvalues <- pcoa_result$eig
variance_explained <- eigenvalues / sum(eigenvalues) * 100

# Format axis labels
pcoa_x_label <- paste0("PCoA Axis 1 (", round(variance_explained[1], 2), "%)")
pcoa_y_label <- paste0("PCoA Axis 2 (", round(variance_explained[2], 2), "%)")

# Define custom colors for depth
#depth_colors <- c("5" = "#9491D9", "10" = "#A2C4EC", "15" = "#B6D7A8", "23.5" = "#DDCC77")

depth_colors <- c("5" = "#9491D9", "10" = "#A2C4EC", "15" = "#44AA99", "23.5" = "#DDCC77")

pcoa_plot <- ggplot(df_merged, aes(x = PCoA1, y = PCoA2, color = Depth, shape = Month)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text(aes(label = Oxygen), vjust = -1.1, size = 3.2, show.legend = FALSE) +  # 👈 suppress unwanted legend
  scale_color_manual(values = depth_colors, name = "Depth (m)") +
  theme_minimal() +
  labs(title = "PCoA Ordination (Bray-Curtis Dissimilarity)",
       x = pcoa_x_label, y = pcoa_y_label,
       shape = "Month") +
  theme(legend.position = "right")



# Save the plot
ggsave("FigS5.png", plot = pcoa_plot, width = 9, height = 6.5, dpi = 300)
