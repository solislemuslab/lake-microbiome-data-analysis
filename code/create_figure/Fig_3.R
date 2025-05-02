library(dplyr)
library(readr)
library(CARlasso)

# Load data
data <- read_csv("Fw_ lake mendota data/coverm_431_MAGS_metagenomes_reads_count.csv")
data_t <- t(data)
data_df <- as.data.frame(data_t)
colnames(data_df) <- data_df[1, ]
data_df <- data_df[-1, ]

# Exclude certain genomes
excluded_columns <- c(
  "Ga0485158_metabat2_ours.098", "Ga0485159_metabat2_ours.079", "Ga0485171_metabat2_ours.127_sub",
  "Ga0485171_metabat2_ours.004", "Ga0485170_maxbin.090", "Ga0485171_metabat1.063",
  "Ga0485157_metabat2_ours.019", "Ga0485166_metabat2_ours.038", "Ga0485167_metabat2_ours.023",
  "Ga0485170_maxbin.059_sub", "Ga0485171_maxbin.130_sub", "Ga0485160_metabat2_ours.158_sub",
  "Ga0485161_maxbin.110", "Ga0485161_metabat1.096", "Ga0485161_metabat2_ours.167_sub"
)
available_columns <- setdiff(colnames(data_df), excluded_columns)
random_column_names <- sample(available_columns, 14)

# Read taxonomy and t_count data
taxonomy_df <- read_tsv("Fw_ lake mendota data/MAG_taxonomy.tsv")
t_count_df <- read_csv("t_count_columnC.csv")

# Map genomes to phylum and rename columns
genome_to_phylum_map <- setNames(taxonomy_df$Phylum, taxonomy_df$Genome)
genomes_list <- c(
  "Ga0485159_maxbin.004_sub","Ga0485169_maxbin.153","Ga0485163_maxbin.002_sub","Ga0485159_metabat2_ours.130_sub",
  "Ga0485164_metabat1.042_sub","Ga0485171_metabat1.063","Ga0485162_metabat2_ours.088",
  "Ga0485161_maxbin.064","Ga0485172_metabat2_ours.152","Ga0485161_metabat2_jgi.003_sub",
  "Ga0485168_metabat2_ours.010","Ga0485157_metabat2_ours.068_sub","Ga0485169_maxbin.226_sub",
  "Ga0485169_metabat2_ours.154","Ga0485160_maxbin.164_sub"
)

new_column_names_list <- vector("list", length(genomes_list))
names(new_column_names_list) <- genomes_list
seq_num <- 1

for (genome in genomes_list) {
  phylum <- sub("^p__", "", genome_to_phylum_map[genome])
  if (!is.null(phylum) && genome %in% colnames(t_count_df)) {
    new_column_name <- paste(phylum, seq_num, sep = "_")
    new_column_names_list[genome] <- new_column_name
    colnames(t_count_df)[colnames(t_count_df) == genome] <- new_column_name
    seq_num <- seq_num + 1
  }
}

# Normalize
set.seed(42)
t_count_df$DO_mg_L <- ifelse(t_count_df$DO_mg_L < 0, 0, t_count_df$DO_mg_L)
last_10_cols <- tail(names(t_count_df), 10)
t_count_df <- t_count_df %>%
  mutate(across(all_of(last_10_cols), ~scale(., center = min(.), scale = max(.) - min(.))[, 1]))

# Build formula and fit CARlasso
predictor_formula_part <- paste(unlist(new_column_names_list), collapse = " + ")
full_formula <- as.formula(paste(predictor_formula_part, "~ depth + wtemp_in_celsius + specific_conductivity + chlorophyll_RFU + phycocyanin_RFU + Dissolved_organic_matter_RFU + Turbidity_RFU + DO_mg_L + PH"))
otu_res <- CARlasso(full_formula, data = t_count_df, adaptive = TRUE, link = "log", n_iter = 5000, n_burn_in = 1000, thin_by = 10)
otu_res <- horseshoe(otu_res)

# Plot
png(filename = "Permutation_Test.png", width = 10, height = 8, units = "in", res = 300)
plot(otu_res)
dev.off()
