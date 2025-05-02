library(dplyr)
library(CARlasso)

# Load taxonomy and create mapping
taxonomy_file_path <- 'Fw_ lake mendota data/MAG_taxonomy.tsv'
taxonomy_df <- read.csv(taxonomy_file_path, sep = '\t')
genome_to_phylum_map <- setNames(taxonomy_df$Phylum, taxonomy_df$Genome)

# Load and process input data
set.seed(42)
changed_col <- read.csv("Network files/tax_Degree.csv")
colnames(changed_col) <- sub("^p__", "", colnames(changed_col))
changed_col$DO_mg_L <- ifelse(changed_col$DO_mg_L < 0, 0, changed_col$DO_mg_L)
last_10_cols <- tail(names(changed_col), 10)
changed_col <- changed_col %>%
  mutate(across(all_of(last_10_cols), ~scale(., center = min(.), scale = max(.) - min(.))[,1]))

# Fit CARlasso model
otu_res <- CARlasso(
  Krumholzibacteriota_1 + Bacteroidota_2 + Bacteroidota_3 + Proteobacteria_4 +
    Planctomycetota_5 + Bacteroidota_6 + Chloroflexota_7 + Bacteroidota_8 +
    Verrucomicrobiota_9 + Actinobacteriota_10 + Planctomycetota_11 +
    Chloroflexota_12 + Proteobacteria_13 + Firmicutes_A_14 + Verrucomicrobiota_15 ~
    depth + wtemp_in_celsius + specific_conductivity + chlorophyll_RFU +
    phycocyanin_RFU + Dissolved_organic_matter_RFU + Turbidity_RFU + DO_mg_L + PH,
  data = changed_col, adaptive = TRUE, link = "log", n_iter = 5000, n_burn_in = 1000, thin_by = 10
)
otu_res <- horseshoe(otu_res)

# Plot
png(filename = "Top_15_Nodes_Degree.png", width = 10, height = 8, units = "in", res = 300)
plot(otu_res)
dev.off()
