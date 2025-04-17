# Set CRAN repository
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# Install packages
install.packages("ggplot2")
install.packages("dplyr")
install.packages("readr")
install.packages("CARlasso")

# Proceed with the rest of the script
library(dplyr)
library(CARlasso)

#Figure 10
taxonomy_file_path = 'Fw_ lake mendota data/MAG_taxonomy.tsv'
taxonomy_df = read.csv(taxonomy_file_path, sep='\t')
genome_to_phylum_map <- setNames(taxonomy_df$Phylum, taxonomy_df$Genome)

set.seed(42)
changed_col<-read.csv("Network files/tax_Degree.csv")
colnames(changed_col) <- sub("^p__", "", colnames(changed_col))
changed_col$DO_mg_L <- ifelse(changed_col$DO_mg_L < 0, 0, changed_col$DO_mg_L)
last_10_cols <- tail(names(changed_col), 10)
changed_col <- changed_col%>%
  mutate(across(all_of(last_10_cols), ~scale(., center = min(.), scale = max(.) - min(.))[,1]))

otu_res <- CARlasso(Krumholzibacteriota_1
                    +Bacteroidota_2
                    +Bacteroidota_3
                    +Proteobacteria_4
                    +Planctomycetota_5
                    +Bacteroidota_6
                    +Chloroflexota_7
                    +Bacteroidota_8
                    +Verrucomicrobiota_9
                    +Actinobacteriota_10
                    +Planctomycetota_11
                    +Chloroflexota_12
                    +Proteobacteria_13
                    +Firmicutes_A_14
                    +Verrucomicrobiota_15~depth+wtemp_in_celsius+specific_conductivity+chlorophyll_RFU+phycocyanin_RFU+Dissolved_organic_matter_RFU+Turbidity_RFU+DO_mg_L+PH, data = changed_col, adaptive = TRUE, link="log", n_iter = 5000, n_burn_in = 1000, thin_by = 10)
otu_res <- horseshoe(otu_res)


genomes_list <- c(
  "Ga0485158_metabat2_ours.098",
  "Ga0485159_metabat2_ours.079",
  "Ga0485171_metabat2_ours.127_sub",
  "Ga0485171_metabat2_ours.004",
  "Ga0485170_maxbin.090",
  "Ga0485171_metabat1.063",
  "Ga0485157_metabat2_ours.019",
  "Ga0485166_metabat2_ours.038",
  "Ga0485167_metabat2_ours.023",
  "Ga0485170_maxbin.059_sub",
  "Ga0485171_maxbin.130_sub",
  "Ga0485160_metabat2_ours.158_sub",
  "Ga0485161_maxbin.110",
  "Ga0485161_metabat1.096",
  "Ga0485161_metabat2_ours.167_sub"
  
)


phylum_names <- sapply(genomes_list, function(genome) sub("^p__", "", genome_to_phylum_map[[genome]]))
new_column_names <- paste(phylum_names, seq_along(genomes_list), sep="_")

reference_table <- data.frame(
  Original_Genome = genomes_list,
  Phylum_Level = new_column_names
)

# Set up parameters for high-resolution output
output_file_path <- "OutputFigures/Top 15 Nodes for Degree with All Environmental Features.png"
width_in_inches <- 10
height_in_inches <- 8
res_dpi <- 300  # Resolution in DPI

# Create output directory if it doesn't exist
if (!dir.exists("OutputFigures")) {
  dir.create("OutputFigures")
}

# Open a PNG device with specified resolution
png(filename = output_file_path, width = width_in_inches, height = height_in_inches, units = "in", res = res_dpi)

# Plot the horseshoe result
plot(otu_res)

# Close the device to finish writing to the file
dev.off()


#Figure 11
library(readr)
data <- read_csv("Fw_ lake mendota data/coverm_431_MAGS_metagenomes_reads_count.csv")
data_t <- t(data)
data_df <- as.data.frame(data_t)
colnames(data_df) <- data_df[1, ]
data_df <- data_df[-1, ]

excluded_columns <- c("Ga0485158_metabat2_ours.098",
                      "Ga0485159_metabat2_ours.079",
                      "Ga0485171_metabat2_ours.127_sub",
                      "Ga0485171_metabat2_ours.004",
                      "Ga0485170_maxbin.090",
                      "Ga0485171_metabat1.063",
                      "Ga0485157_metabat2_ours.019",
                      "Ga0485166_metabat2_ours.038",
                      "Ga0485167_metabat2_ours.023",
                      "Ga0485170_maxbin.059_sub",
                      "Ga0485171_maxbin.130_sub",
                      "Ga0485160_metabat2_ours.158_sub",
                      "Ga0485161_maxbin.110",
                      "Ga0485161_metabat1.096",
                      "Ga0485161_metabat2_ours.167_sub")

set.seed(42)
available_columns <- setdiff(colnames(data_df), excluded_columns)
random_column_names <- sample(available_columns, 14)
#print(random_column_names)
taxonomy_file_path <- 'Fw_ lake mendota data/MAG_taxonomy.tsv'
taxonomy_df <- read_tsv(taxonomy_file_path)
t_count_file_path <- 't_count_columnC.csv'
t_count_df <- read_csv(t_count_file_path)

# Create a mapping of genomes to phylum
genome_to_phylum_map <- setNames(taxonomy_df$Phylum, taxonomy_df$Genome)

genomes_list <- c(
  "Ga0485159_maxbin.004_sub","Ga0485169_maxbin.153","Ga0485163_maxbin.002_sub","Ga0485159_metabat2_ours.130_sub","Ga0485164_metabat1.042_sub","Ga0485171_metabat1.063","Ga0485162_metabat2_ours.088" ,"Ga0485161_maxbin.064","Ga0485172_metabat2_ours.152" ,"Ga0485161_metabat2_jgi.003_sub","Ga0485168_metabat2_ours.010","Ga0485157_metabat2_ours.068_sub","Ga0485169_maxbin.226_sub","Ga0485169_metabat2_ours.154","Ga0485160_maxbin.164_sub"   
)

new_column_names_list <- vector("list", length(genomes_list))
names(new_column_names_list) <- genomes_list

seq_num <- 1

for (genome in genomes_list) {
  phylum <- sub("^p__", "", genome_to_phylum_map[genome])
  
  if (!is.null(phylum) && genome %in% colnames(t_count_df)) {
    # Construct the new column name using the phylum and the current sequence number
    new_column_name <- paste(phylum, seq_num, sep = "_")
    
    # Update the list with the new column name
    new_column_names_list[genome] <- new_column_name
    
    # Rename the column in t_count_df
    colnames(t_count_df)[colnames(t_count_df) == genome] <- new_column_name
    
    # Increment the sequence number for the next iteration
    seq_num <- seq_num + 1
  }
}
set.seed(42)
t_count_df$DO_mg_L <- ifelse(t_count_df$DO_mg_L < 0, 0, t_count_df$DO_mg_L)
last_10_cols <- tail(names(t_count_df), 10)
t_count_df <- t_count_df%>%
  mutate(across(all_of(last_10_cols), ~scale(., center = min(.), scale = max(.) - min(.))[,1]))

predictor_formula_part <- paste(new_column_names_list, collapse = " + ")
# Construct the full formula by appending the response variables and predictors
full_formula <- as.formula(paste(predictor_formula_part,"~", "depth + wtemp_in_celsius + specific_conductivity + chlorophyll_RFU + phycocyanin_RFU + Dissolved_organic_matter_RFU + Turbidity_RFU + DO_mg_L + PH "))
otu_res <- CARlasso(full_formula, data = t_count_df, adaptive = TRUE, link="log", n_iter = 5000, n_burn_in = 1000, thin_by = 10)


otu_res <- horseshoe(otu_res)

# Set up parameters for high-resolution output
output_file_path <- "OutputFigures/permutation.png"
width_in_inches <- 10
height_in_inches <- 8
res_dpi <- 300  # Resolution in DPI

# Create output directory if it doesn't exist
if (!dir.exists("OutputFigures")) {
  dir.create("OutputFigures")
}

# Open a PNG device with specified resolution
png(filename = output_file_path, width = width_in_inches, height = height_in_inches, units = "in", res = res_dpi)

# Plot the horseshoe result
plot(otu_res)

# Close the device to finish writing to the file
dev.off()

## figure for Metatranscriptome in the supplementary file


set.seed(42)

# Load data
changed_col <- read.csv("metatransciptome_Networkfile/tax_Degree.csv")
taxonomy_df <- read.csv("Fw_ lake mendota data/MAG_taxonomy.tsv", sep = "\t")

colnames(changed_col) <- sub("^p__", "", colnames(changed_col))

# Clean and normalize DO column
changed_col$DO_mg_L <- ifelse(changed_col$DO_mg_L < 0, 0, changed_col$DO_mg_L)

# Normalize the last 10 columns
last_10_cols <- tail(names(changed_col), 10)
changed_col <- changed_col %>%
  mutate(across(all_of(last_10_cols), ~scale(., center = min(.), scale = max(.) - min(.))[, 1]))

# Run CARlasso model
otu_res <- CARlasso(
  Krumholzibacteriota_1 + Bacteroidota_2 + Bacteroidota_3 + Proteobacteria_4 +
    Planctomycetota_5 + Bacteroidota_6 + Chloroflexota_7 + Bacteroidota_8 +
    Verrucomicrobiota_9 + Actinobacteriota_10 + Planctomycetota_11 +
    Chloroflexota_12 + Proteobacteria_13 + Firmicutes_A_14 + Verrucomicrobiota_15 ~
    depth + wtemp_in_celsius + specific_conductivity + chlorophyll_RFU +
    phycocyanin_RFU + Dissolved_organic_matter_RFU + Turbidity_RFU + DO_mg_L + PH,
  data = changed_col,
  adaptive = TRUE, link = "log", n_iter = 5000, n_burn_in = 1000, thin_by = 10
)

# Apply horseshoe shrinkage
otu_res <- horseshoe(otu_res)

# Genomes of interest
genomes_list <- c(
  "Ga0485158_metabat2_ours.098", "Ga0485159_metabat2_ours.079", "Ga0485171_metabat2_ours.127_sub",
  "Ga0485171_metabat2_ours.004", "Ga0485170_maxbin.090", "Ga0485171_metabat1.063",
  "Ga0485157_metabat2_ours.019", "Ga0485166_metabat2_ours.038", "Ga0485167_metabat2_ours.023",
  "Ga0485170_maxbin.059_sub", "Ga0485171_maxbin.130_sub", "Ga0485160_metabat2_ours.158_sub",
  "Ga0485161_maxbin.110", "Ga0485161_metabat1.096", "Ga0485161_metabat2_ours.167_sub"
)

# Create taxonomy mappings
genome_to_domain_map <- setNames(taxonomy_df$Domain, taxonomy_df$Genome)
genome_to_phylum_map <- setNames(taxonomy_df$Phylum, taxonomy_df$Genome)
genome_to_class_map <- setNames(taxonomy_df$Class, taxonomy_df$Genome)
genome_to_order_map <- setNames(taxonomy_df$Order, taxonomy_df$Genome)
genome_to_family_map <- setNames(taxonomy_df$Family, taxonomy_df$Genome)
genome_to_genus_map <- setNames(taxonomy_df$Genus, taxonomy_df$Genome)
genome_to_species_map <- setNames(taxonomy_df$Species, taxonomy_df$Genome)

# Extract taxonomy names
domain_names <- sapply(genomes_list, function(genome) genome_to_domain_map[[genome]])
phylum_names <- sapply(genomes_list, function(genome) genome_to_phylum_map[[genome]])
class_names <- sapply(genomes_list, function(genome) genome_to_class_map[[genome]])
order_names <- sapply(genomes_list, function(genome) genome_to_order_map[[genome]])
family_names <- sapply(genomes_list, function(genome) genome_to_family_map[[genome]])
genus_names <- sapply(genomes_list, function(genome) genome_to_genus_map[[genome]])
species_names <- sapply(genomes_list, function(genome) genome_to_species_map[[genome]])

# Clean up phylum names (remove "p__" if needed)
phylum_names <- sub("^p__", "", phylum_names)

# Create readable node names
new_column_names <- paste(phylum_names, seq_along(genomes_list), sep = "_")

# Create output directory if it doesn't exist
if (!dir.exists("OutputFigures")) {
  dir.create("OutputFigures")
}

# Save plot as PNG
png(filename = "OutputFigures/Metatranscriptome_supp.png",
    width = 10, height = 8, units = "in", res = 300)

plot(otu_res)

dev.off()



