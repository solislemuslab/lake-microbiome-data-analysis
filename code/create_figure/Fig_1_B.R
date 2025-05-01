# === Load required libraries ===
library(tidyverse)
library(ggvenn)
library(readr)

# === Load presence matrix (MAG x 11 conditions) ===
presence_matrix <- read_csv("presence_matrix.csv")

# === Define groups (columns) ===
depth_cols <- c("5", "10", "15","23.5")
month_cols <- c("July", "August", "September","October")
oxygen_cols <- c("Oxic", "Anoxic", "Oxycline")

# === Convert to presence lists for each group ===
get_mag_lists <- function(cols) {
  lapply(cols, function(col) {
    presence_matrix %>% filter(!!sym(col) == 1) %>% pull(MAG)
  }) %>% setNames(cols)
}

depth_sets <- get_mag_lists(depth_cols)
month_sets <- get_mag_lists(month_cols)
oxygen_sets <- get_mag_lists(oxygen_cols)

# Define colors
venn_colors <- c("#DDCC77", "#11A0D9", "#44AA99", "#9491D9")


# === Plot Venn: Depth ===
png("venn_depths.png", width = 2500, height = 2200, res = 300)
ggvenn(depth_sets,
       fill_color = venn_colors,
       show_percentage = TRUE,
       stroke_size = 1,
       set_name_size = 6,
       text_size = 5)
dev.off()

# === Plot Venn: Month ===
png("venn_months.png", width = 2500, height = 2200, res = 300)
ggvenn(month_sets,
       fill_color = venn_colors,
       show_percentage = TRUE,
       stroke_size = 1,
       set_name_size = 6,
       text_size = 5)
dev.off()

# === Plot Venn: Oxygen ===
png("Fig_1_B.png", width = 2500, height = 2200, res = 300)
ggvenn(oxygen_sets,
       fill_color = venn_colors,
       show_percentage = TRUE,
       stroke_size = 1,
       set_name_size = 6,
       text_size = 5)
dev.off()
