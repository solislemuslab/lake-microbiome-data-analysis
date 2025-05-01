# Load necessary libraries
library(tidyverse)
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

# Normalize (TSS)
merged_df[mag_cols] <- merged_df[mag_cols] / rowSums(merged_df[mag_cols])

# Reshape to long format
long_df <- merged_df %>%
  pivot_longer(cols = all_of(mag_cols), names_to = "MAG", values_to = "Relative_Abundance")

# ----------- Add sample count and custom sorted labels -----------

# DEPTH: order by 5, 10, 15, 23.5
depth_counts <- merged_df %>%
  mutate(Depth = as.character(Depth)) %>%
  count(Depth) %>%
  mutate(label = factor(paste0(Depth, "m (n=", n, ")"),
                        levels = paste0(c(5, 10, 15, 23.5), "m (n=", 
                                        c(
                                          n[Depth == "5"],
                                          n[Depth == "10"],
                                          n[Depth == "15"],
                                          n[Depth == "23.5"]
                                        ), ")")))

# MONTH: order by July → August → September → October
month_counts <- merged_df %>%
  mutate(Month = as.character(Month)) %>%
  count(Month) %>%
  mutate(label = factor(paste0(Month, " (n=", n, ")"),
                        levels = paste0(c("July", "August", "September", "October"), " (n=",
                                        c(
                                          n[Month == "July"],
                                          n[Month == "August"],
                                          n[Month == "September"],
                                          n[Month == "October"]
                                        ), ")")))

# OXYGEN: default order or customize if needed
oxygen_counts <- merged_df %>%
  mutate(Oxygen = as.character(Oxygen)) %>%
  count(Oxygen) %>%
  mutate(label = paste0(Oxygen, " (n=", n, ")"))

# ----------- Join labels back -----------

long_df_depth <- long_df %>%
  mutate(Depth = as.character(Depth)) %>%
  left_join(depth_counts, by = "Depth")

long_df_month <- long_df %>%
  mutate(Month = as.character(Month)) %>%
  left_join(month_counts, by = "Month")

long_df_oxygen <- long_df %>%
  mutate(Oxygen = as.character(Oxygen)) %>%
  left_join(oxygen_counts, by = "Oxygen")

# ----------- Plots -----------

# Plot: by Depth
plot_depth <- ggplot(long_df_depth, aes(x = label, y = Relative_Abundance)) +
  geom_boxplot(fill = "#A2C4EC", outlier.size = 0.5) +
  theme_minimal(base_size = 13) +
  labs(x = "(a) Depth", y = "Relative Abundance of MAGs") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# Plot: by Month
plot_month <- ggplot(long_df_month, aes(x = label, y = Relative_Abundance)) +
  geom_boxplot(fill = "#B6D7A8", outlier.size = 0.5) +
  theme_minimal(base_size = 13) +
  labs(x = "(b) Month", y = "") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# Plot: by Oxygen
plot_oxygen <- ggplot(long_df_oxygen, aes(x = label, y = Relative_Abundance)) +
  geom_boxplot(fill = "#DDCC77", outlier.size = 0.5) +
  theme_minimal(base_size = 13) +
  labs(x = "(c) Oxygen Level", y = "") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# Combine
final_plot <- ggarrange(plot_depth, plot_month, plot_oxygen, ncol = 3)

# Add title
final_plot <- annotate_figure(
  final_plot,
  top = text_grob("",
                  face = "bold", size = 14)
)

# Save
ggsave("Fig_S4.pdf", plot = final_plot, width = 13, height = 6, dpi = 400)
