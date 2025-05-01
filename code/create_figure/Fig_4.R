rm = list(ls())
# Load required libraries
library(dplyr)
library(readr)
library(ggplot2)
library(ggalluvial)
library(stringr)

# === Load data ===
abundance_data <- read_csv("coverm_431_MAGS_metagenomes_reads_count.csv", show_col_types = FALSE)
# Convert to relative abundance (row-wise TSS)
abundance_data_rel <- abundance_data %>%
  mutate(across(-Genome, ~ . / sum(., na.rm = TRUE)))

relative_abundance = abundance_data_rel

taxonomy_data <- read_delim("MAG_taxonomy.tsv", delim = "\t", show_col_types = FALSE)

# === Define MAGs in each table ===
table1_taxa <- c("Ga0485158_metabat2_ours.098", "Ga0485159_metabat2_ours.079", "Ga0485171_metabat2_ours.127_sub",
                 "Ga0485171_metabat2_ours.004", "Ga0485170_maxbin.090", "Ga0485157_metabat2_ours.019", 
                 "Ga0485166_metabat2_ours.038", "Ga0485167_metabat2_ours.023", "Ga0485170_maxbin.059_sub",
                 "Ga0485171_maxbin.130_sub", "Ga0485160_metabat2_ours.158_sub", "Ga0485161_maxbin.110", 
                 "Ga0485161_metabat1.096", "Ga0485161_metabat2_ours.167_sub", "Ga0485171_metabat1.063")

table2_taxa <- c("Ga0485159_maxbin.004_sub", "Ga0485169_maxbin.153", "Ga0485163_maxbin.002_sub", 
                 "Ga0485159_metabat2_ours.130_sub", "Ga0485164_metabat1.042_sub", "Ga0485162_metabat2_ours.088", 
                 "Ga0485161_maxbin.064", "Ga0485172_metabat2_ours.152", "Ga0485161_metabat2_jgi.003_sub", 
                 "Ga0485168_metabat2_ours.010", "Ga0485157_metabat2_ours.068_sub", "Ga0485169_maxbin.226_sub", 
                 "Ga0485169_metabat2_ours.154", "Ga0485160_maxbin.164_sub", "Ga0485171_metabat1.063")

# === Combine 19 MAGs ===
all_taxa <- union(table1_taxa, table2_taxa)  # ensures shared MAG appears once

# === Compute average abundance only for those 29 MAGs ===
avg_abundance <- relative_abundance %>%
  filter(Genome %in% all_taxa) %>%
  mutate(Avg_Abundance = rowMeans(select(., -Genome), na.rm = TRUE)) %>%
  select(Genome, Avg_Abundance)

# === Assign groups — this will now include "Both" ===
group_df <- tibble(
  Genome = all_taxa,
  Group = case_when(
    Genome %in% table1_taxa & Genome %in% table2_taxa ~ "Both",
    Genome %in% table1_taxa ~ "Table 1",
    Genome %in% table2_taxa ~ "Table 2"
  )
)

# === Clean taxonomy ===
taxonomy_clean <- taxonomy_data %>%
  filter(Genome %in% all_taxa) %>%
  mutate(
    Phylum = str_remove(Phylum, "p__"),
    Class = str_remove(Class, "c__"),
    Family = str_remove(Family, "f__")
  ) %>%
  select(Genome, Phylum, Class, Family)

# === Merge everything ===
plot_data <- avg_abundance %>%
  left_join(group_df, by = "Genome") %>%
  left_join(taxonomy_clean, by = "Genome") %>%
  filter(!is.na(Phylum) & !is.na(Class) & !is.na(Family))

# === Prepare for ggalluvial ===
alluvial_data <- plot_data %>%
  select(Group, Genome, Phylum, Class, Family, Avg_Abundance)

# Set group order
alluvial_data$Group <- factor(alluvial_data$Group,
                              levels = c("Both", "Table 1", "Table 2"),
                              labels = c("Shared MAG", "Table 1", "Table 2"))


# === Custom Phylum colors ===
phylum_colors <- c(
  "Proteobacteria"      = "#D8BFD8",
  "Desulfobacterota"    = "#F4C28C",
  "Bacteroidota"        = "#9491D9",
  "Myxococcota"         = "#82C341",
  "Verrucomicrobiota"   = "#A2C4EC",
  "Actinobacteriota"    = "#B6D7A8",
  "Planctomycetota"     = "#A94442",
  "Cyanobacteria"       = "#F4CCCC",
  "Krumholzibacteriota" = "#D1EEEE",#
  "Chloroflexota"       = "#DDCC77",
  "Firmicutes_A"        = "#44AA99"
)



ggplot(alluvial_data,
       aes(axis1 = Group, axis2 = Phylum, axis3 = Class, axis4 = Family, axis5 = Genome, y = Avg_Abundance)) +
  
  # BORDER layer (slightly larger, darker version behind)
  geom_alluvium(aes(fill = Phylum), width = 0.06, alpha = 1, color = "gray30", size = 0.2) +
  
  # MAIN layer (slightly narrower and on top)
  geom_alluvium(aes(fill = Phylum), width = 0.045, size = 0.4,alpha = 1) +
  
  geom_stratum(fill = "white", color = "black",linetype = "solid", width = 1/14) +
  #geom_stratum(fill = "white", color = "black", linetype = "solid", size = 0.4, width = 1/14)

  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4) +
  scale_x_discrete(
    limits = c("Group", "Phylum", "Class", "Family", "MAGs"),
    expand = c(0.1, 0.05)
  ) +
  scale_fill_manual(
    values = c(
      phylum_colors,
      "Table 1" = "gray90",
      "Table 2" = "gray85",
      "Shared MAG" = "white"
    )
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 13, face = "bold")
  )
# === Save ===
ggsave("Figure4.pdf", width = 15, height = 20)
