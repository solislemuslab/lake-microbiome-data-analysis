# Load necessary libraries
library(tidyverse)
library(ggpubr)
library(forcats)

# Load data
abundance <- read_csv("coverm_431_MAGS_metagenomes_reads_count.csv")
samples <- read_csv("Samples-mendota.csv")

# Transpose abundance data: samples as rows
abundance_t <- abundance %>%
  pivot_longer(-Genome, names_to = "Sample.Name", values_to = "Count") %>%
  pivot_wider(names_from = Genome, values_from = Count)

# Merge with sample metadata
merged_df <- left_join(samples, abundance_t, by = "Sample.Name")

# Identify MAG columns (exclude metadata)
mag_cols <- setdiff(colnames(merged_df), colnames(samples))

# Normalize abundance (TSS)
merged_df[mag_cols] <- merged_df[mag_cols] / rowSums(merged_df[mag_cols])
# Pivot to long format
long_df <- merged_df %>%
  pivot_longer(cols = all_of(mag_cols), names_to = "MAG", values_to = "Relative_Abundance")

# ---------- Add Sample Counts and Ordered Labels ----------

# Set Depth as ordered factor and make counts
depth_counts <- merged_df %>%
  mutate(Depth = factor(Depth, levels = c(5, 10, 15, 23.5))) %>%
  count(Depth) %>%
  mutate(Depth_Label = paste0(Depth, "m (n=", n, ")"))

# Month in chronological order
month_counts <- merged_df %>%
  mutate(Month = factor(Month, levels = c("July", "August", "September", "October"))) %>%
  count(Month) %>%
  mutate(Month_Label = paste0(Month, " (n=", n, ")"))


# Fix ordering in label columns

# Sort Depth_Label manually
depth_counts <- depth_counts %>%
  mutate(Depth_Label = factor(Depth_Label, 
                              levels = paste0(c(5, 10, 15, 23.5), "m (n=", c(
                                depth_counts$n[depth_counts$Depth == 5],
                                depth_counts$n[depth_counts$Depth == 10],
                                depth_counts$n[depth_counts$Depth == 15],
                                depth_counts$n[depth_counts$Depth == 23.5]
                              ), ")")))

# Sort Month_Label manually
month_counts <- month_counts %>%
  mutate(Month_Label = factor(Month_Label, 
                              levels = paste0(c("July", "August", "September", "October"), 
                                              " (n=", c(
                                                month_counts$n[month_counts$Month == "July"],
                                                month_counts$n[month_counts$Month == "August"],
                                                month_counts$n[month_counts$Month == "September"],
                                                month_counts$n[month_counts$Month == "October"]
                                              ), ")")))



# Oxygen (default order or customize if needed)
oxygen_counts <- merged_df %>%
  count(Oxygen) %>%
  mutate(Oxygen_Label = paste0(Oxygen, " (n=", n, ")"))

# ---------- Mean MAG abundance per group (431 values per group) ----------

# By Depth (convert to character to avoid join errors)
depth_mean <- long_df %>%
  mutate(Depth = as.character(Depth)) %>%
  group_by(Depth, MAG) %>%
  summarise(Mean_Relative_Abundance = mean(Relative_Abundance), .groups = "drop") %>%
  left_join(depth_counts %>% mutate(Depth = as.character(Depth)), by = "Depth")

# By Month
month_mean <- long_df %>%
  mutate(Month = as.character(Month)) %>%
  group_by(Month, MAG) %>%
  summarise(Mean_Relative_Abundance = mean(Relative_Abundance), .groups = "drop") %>%
  left_join(month_counts %>% mutate(Month = as.character(Month)), by = "Month")

# By Oxygen
oxygen_mean <- long_df %>%
  mutate(Oxygen = as.character(Oxygen)) %>%
  group_by(Oxygen, MAG) %>%
  summarise(Mean_Relative_Abundance = mean(Relative_Abundance), .groups = "drop") %>%
  left_join(oxygen_counts %>% mutate(Oxygen = as.character(Oxygen)), by = "Oxygen")

# ---------- Plotting ----------

# Depth
plot_depth <- ggplot(depth_mean, aes(x = Depth_Label, y = Mean_Relative_Abundance)) +
  geom_boxplot(fill = "#A2C4EC", outlier.size = 0.5) +
  theme_minimal(base_size = 13) +
  labs(x = "(a) Depth", y = "Mean Relative Abundance of MAGs") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# Month
plot_month <- ggplot(month_mean, aes(x = Month_Label, y = Mean_Relative_Abundance)) +
  geom_boxplot(fill = "#B6D7A8", outlier.size = 0.5) +
  theme_minimal(base_size = 13) +
  labs(x = "(b) Month", y = "") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# Oxygen
plot_oxygen <- ggplot(oxygen_mean, aes(x = fct_inorder(Oxygen_Label), y = Mean_Relative_Abundance)) +
  geom_boxplot(fill = "#DDCC77", outlier.size = 0.5) +
  theme_minimal(base_size = 13) +
  labs(x = "(c) Oxygen Level", y = "") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# ---------- Combine and Save ----------

final_plot <- ggarrange(plot_depth, plot_month, plot_oxygen, ncol = 3)

final_plot <- annotate_figure(
  final_plot,
  top = text_grob("",
                  face = "bold", size = 14)
)

ggsave("FigS3.pdf",
       plot = final_plot, width = 13, height = 6, dpi = 400)
