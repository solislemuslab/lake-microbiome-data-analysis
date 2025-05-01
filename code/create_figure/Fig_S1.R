# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)

# Load data
abundance <- read_csv("coverm_431_MAGS_metagenomes_reads_count.csv")
samples <- read_csv("Samples-mendota.csv")

# Transpose abundance: rows = samples
abundance_t <- abundance %>%
  pivot_longer(-Genome, names_to = "Sample.Name", values_to = "Count") %>%
  pivot_wider(names_from = Genome, values_from = Count)

# Merge with sample metadata
merged_df <- left_join(samples, abundance_t, by = "Sample.Name")

# Identify MAG columns (exclude metadata)
mag_cols <- setdiff(colnames(merged_df), colnames(samples))

# ✅ Extract Month and Day from Sample.Name (format: 2020-08-05)
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
    Month = factor(Month, levels = c("July", "August", "September", "October"))
  )

# ✅ Create label: Month-Day–Depth–Oxygen
merged_df <- merged_df %>%
  mutate(
    Depth = as.character(Depth),
    Sample_Label = paste(Month, Day, paste0(Depth, "m"), Oxygen, sep = "-")
  )

# ✅ Sort x-axis labels by Month then numeric Depth
sample_order <- merged_df %>%
  mutate(Depth_num = as.numeric(Depth)) %>%
  arrange(Month, Depth_num) %>%
  pull(Sample_Label)

# Reshape to long format (raw counts)
long_df_samples <- merged_df %>%
  select(Sample_Label, all_of(mag_cols)) %>%
  pivot_longer(cols = all_of(mag_cols), names_to = "MAG", values_to = "Abundance") %>%
  mutate(Sample_Label = factor(Sample_Label, levels = unique(sample_order)))

# Plot
plot_samples <- ggplot(long_df_samples, aes(x = Sample_Label, y = Abundance)) +
  geom_boxplot(fill = "#D9CCE3", outlier.size = 0.5) +
  theme_minimal(base_size = 13) +
  labs(x = "Sample (Month-Day-Depth-Oxygen)", y = "Abundance of MAGs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

# Save plot
ggsave("FigS1.pdf", plot = plot_samples, width = 12, height = 5, dpi = 400)
