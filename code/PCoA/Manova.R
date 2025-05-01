# Clear workspace and load packages
rm(list=ls())
library(vegan)
library(pairwiseAdonis)
library(tidyverse)
library(readr)
set.seed(1234)
# Load abundance and metadata
abundance <- read_csv("coverm_431_MAGS_metagenomes_reads_count.csv")
metadata <- read_csv("Samples-mendota.csv")

# Transpose abundance data
abundance_t <- abundance %>%
  pivot_longer(-Genome, names_to = "Sample", values_to = "Count") %>%
  pivot_wider(names_from = Genome, values_from = Count)

# Merge with metadata
df <- left_join(metadata, abundance_t, by = c("Sample.Name" = "Sample"))

# Extract abundance matrix and ensure numeric
# Use only MAG columns from abundance matrix (avoid metadata columns)
mag_cols <- setdiff(colnames(abundance_t), "Sample")
abund_mat <- df[, mag_cols] %>%
  mutate_all(as.numeric)
abund_mat[is.na(abund_mat)] <- 0


# Hellinger transformation
hellinger <- decostand(abund_mat, method = "hellinger")

# Bray-Curtis dissimilarity
bray <- vegdist(hellinger, method = "bray")

# Make sure metadata variables are factors
df$Month <- factor(df$Month, levels = c("July", "August", "September", "October"))
df$Depth <- factor(df$Depth, levels = c(5, 10, 15, 23.5))
df$Oxygen <- factor(df$Oxygen)

# ----------- PERMANOVA main effects -----------
cat("\n--- PERMANOVA: Depth ---\n")
print(adonis2(bray ~ Depth, data = df, permutations = 999))

cat("\n--- PERMANOVA: Month ---\n")
print(adonis2(bray ~ Month, data = df, permutations = 999))

cat("\n--- PERMANOVA: Oxygen ---\n")
print(adonis2(bray ~ Oxygen, data = df, permutations = 999))

# ----------- PERMANOVA interaction and combined models -----------
cat("\n--- PERMANOVA: Depth + Month + Oxygen ---\n")
print(adonis2(bray ~ Depth + Month + Oxygen, data = df, permutations = 999))

# ----------- PERMANOVA interaction and combined models -----------
cat("\n--- PERMANOVA: Depth * Month * Oxygen ---\n")
print(adonis2(bray ~ Depth * Month * Oxygen, data = df, permutations = 999))


cat("\n--- PERMANOVA: Depth * Month ---\n")
print(adonis2(bray ~ Depth * Month, data = df, permutations = 999))

cat("\n--- PERMANOVA: Month * Oxygen ---\n")
print(adonis2(bray ~ Month * Oxygen, data = df, permutations = 999))

cat("\n--- PERMANOVA: Depth * Oxygen ---\n")
print(adonis2(bray ~ Depth * Oxygen, data = df, permutations = 999))


# ----------- Pairwise PERMANOVA -----------
cat("\n--- Pairwise PERMANOVA: Depth ---\n")
print(pairwise.adonis2(bray ~ Depth, data = df, permutations = 999))

cat("\n--- Pairwise PERMANOVA: Month ---\n")
print(pairwise.adonis2(bray ~ Month, data = df, permutations = 999))

cat("\n--- Pairwise PERMANOVA: Oxygen ---\n")
print(pairwise.adonis2(bray ~ Oxygen, data = df, permutations = 999))
