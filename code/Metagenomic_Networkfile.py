import os
import pandas as pd

# Directory for saving output files
folder_path = 'Metagenomic_Networkfile'
os.makedirs(folder_path, exist_ok=True) 

taxonomy_file_path = 'MAG_taxonomy.tsv'
taxonomy_df = pd.read_csv(taxonomy_file_path, sep='\t')

t_count_file_path = 'Metagenomic_columnC.csv'
t_count_df = pd.read_csv(t_count_file_path)


#top10 OTUs
genome_to_phylum_map = dict(zip(taxonomy_df['Genome'], taxonomy_df['Phylum']))

# List of genome names to map to Phylum
genomes_list = [
    'Ga0485157_metabat1.059',
    'Ga0485167_maxbin.109',
    'Ga0485162_maxbin.089',
    'Ga0485161_maxbin.110',
    'Ga0485161_maxbin.075',
    'Ga0485157_metabat1.036',
    'Ga0485165_metabat2_ours.012_sub',
    'Ga0485169_maxbin.201_sub',
    'Ga0485172_maxbin.081_sub',
    'Ga0485168_maxbin.153'
    
]
# Resetting the DataFrame to its original state before reapplying changes with the new requirement
t_count_df = pd.read_csv(t_count_file_path)

# Applying the updated mapping with sequence numbers
for i, genome in enumerate(genomes_list, start=1):
    phylum = genome_to_phylum_map.get(genome)
    if phylum and genome in t_count_df.columns:
        new_column_name = f"{phylum}_{i}"  # Append underscore and sequence number
        t_count_df = t_count_df.rename(columns={genome: new_column_name})

t_count_df.to_csv(os.path.join(folder_path, 'tax_10.csv'), index=True)

#top15 OTUs
# List of genome names to map to Phylum
genomes_list = [
    'Ga0485157_metabat1.059',
    'Ga0485167_maxbin.109',
    'Ga0485162_maxbin.089',
    'Ga0485161_maxbin.110',
    'Ga0485161_maxbin.075',
    'Ga0485157_metabat1.036',
    'Ga0485165_metabat2_ours.012_sub',
    'Ga0485169_maxbin.201_sub',
    'Ga0485172_maxbin.081_sub',
    'Ga0485168_maxbin.153',
    'Ga0485172_metabat2_ours.083',
    'Ga0485162_maxbin.023',
    'Ga0485160_maxbin.092',
    'Ga0485158_metabat2_jgi.024',
    'Ga0485162_metabat1.001'
]


# Applying the updated mapping with sequence numbers
for i, genome in enumerate(genomes_list, start=1):
    phylum = genome_to_phylum_map.get(genome)
    if phylum and genome in t_count_df.columns:
        new_column_name = f"{phylum}_{i}"  # Append underscore and sequence number
        t_count_df = t_count_df.rename(columns={genome: new_column_name})
        
t_count_df.to_csv(os.path.join(folder_path, 'tax_15.csv'), index=True)


# betweenness

# List of genome names to map to Phylum
genomes_list = [
    "Ga0485158_maxbin.120_sub",
"Ga0485170_maxbin.090",
"Ga0485159_metabat1.040",
"Ga0485167_metabat2_ours.163",
"Ga0485166_metabat2_ours.038",
"Ga0485162_metabat2_ours.116",
"Ga0485162_metabat2_ours.050",
"Ga0485163_metabat1.181",
"Ga0485158_metabat2_ours.098",
"Ga0485169_metabat2_ours.035",
"Ga0485159_metabat2_ours.130_sub",
"Ga0485163_metabat1.115_sub",
"Ga0485159_maxbin.025",
"Ga0485161_metabat1.096",
"Ga0485167_metabat1.127_sub"

]
# Load the taxonomy file again (repeating due to re-upload)
taxonomy_file_path = 'MAG_taxonomy.tsv'
taxonomy_df = pd.read_csv(taxonomy_file_path, sep='\t')

# Load the T-count file
t_count_file_path = 't_count_columnC.csv'
t_count_df = pd.read_csv(t_count_file_path)
genome_to_phylum_map = dict(zip(taxonomy_df['Genome'], taxonomy_df['Phylum']))
t_count_df = pd.read_csv(t_count_file_path)


# Applying the updated mapping with sequence numbers
for i, genome in enumerate(genomes_list, start=1):
    phylum = genome_to_phylum_map.get(genome)
    if phylum and genome in t_count_df.columns:
        new_column_name = f"{phylum}_{i}"  # Append underscore and sequence number
        #print(new_column_name)
        t_count_df = t_count_df.rename(columns={genome: new_column_name})

t_count_df.to_csv(os.path.join(folder_path, 'tax_Betweenness.csv'), index=True)



#Closeness
# List of genome names to map to Phylum
genomes_list = [
    "Ga0485157_metabat1.069",
"Ga0485157_metabat2_ours.025",
"Ga0485157_metabat2_ours.091",
"Ga0485158_metabat2_ours.184",
"Ga0485159_maxbin.015",
"Ga0485160_maxbin.155_sub",
"Ga0485160_metabat2_ours.051",
"Ga0485160_metabat2_ours.071",
"Ga0485161_metabat1.109_sub",
"Ga0485162_metabat2_ours.036",
"Ga0485163_maxbin.002_sub",
"Ga0485163_maxbin.104_sub",
"Ga0485163_maxbin.129",
"Ga0485163_maxbin.130_sub",
"Ga0485163_maxbin.183_sub",

]
# Load the taxonomy file again (repeating due to re-upload)
taxonomy_file_path = 'MAG_taxonomy.tsv'
taxonomy_df = pd.read_csv(taxonomy_file_path, sep='\t')

# Load the T-count file
t_count_file_path = 't_count_columnC.csv'
t_count_df = pd.read_csv(t_count_file_path)
genome_to_phylum_map = dict(zip(taxonomy_df['Genome'], taxonomy_df['Phylum']))
t_count_df = pd.read_csv(t_count_file_path)


# Applying the updated mapping with sequence numbers
for i, genome in enumerate(genomes_list, start=1):
    phylum = genome_to_phylum_map.get(genome)
    if phylum and genome in t_count_df.columns:
        new_column_name = f"{phylum}_{i}"  # Append underscore and sequence number
        #print(new_column_name)
        t_count_df = t_count_df.rename(columns={genome: new_column_name})

t_count_df.to_csv(os.path.join(folder_path, 'tax_Closeness.csv'), index=True)


#Page Rank
# List of genome names to map to Phylum
genomes_list = [
    "Ga0485158_metabat2_ours.098",
"Ga0485166_metabat2_ours.038",
"Ga0485167_metabat2_ours.023",
"Ga0485170_maxbin.090",
"Ga0485161_metabat2_ours.167_sub",
"Ga0485171_maxbin.130_sub",
"Ga0485161_maxbin.110",
"Ga0485164_metabat2_ours.069_sub",
"Ga0485171_metabat1.063",
"Ga0485159_metabat2_ours.079",
"Ga0485171_metabat2_ours.127_sub",
"Ga0485159_metabat2_ours.155_sub",
"Ga0485158_metabat1.076",
"Ga0485171_metabat1.030",
"Ga0485168_metabat2_ours.135_sub"
]
# Load the taxonomy file again (repeating due to re-upload)
taxonomy_file_path = 'MAG_taxonomy.tsv'
taxonomy_df = pd.read_csv(taxonomy_file_path, sep='\t')

# Load the T-count file
t_count_file_path = 't_count_columnC.csv'
t_count_df = pd.read_csv(t_count_file_path)
genome_to_phylum_map = dict(zip(taxonomy_df['Genome'], taxonomy_df['Phylum']))
t_count_df = pd.read_csv(t_count_file_path)


# Applying the updated mapping with sequence numbers
for i, genome in enumerate(genomes_list, start=1):
    phylum = genome_to_phylum_map.get(genome)
    if phylum and genome in t_count_df.columns:
        new_column_name = f"{phylum}_{i}"  # Append underscore and sequence number
        #print(new_column_name)
        t_count_df = t_count_df.rename(columns={genome: new_column_name})

t_count_df.to_csv(os.path.join(folder_path, 'tax_Page_Rank.csv'), index=True)



#Degree
# List of genome names to map to Phylum
genomes_list = [
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

]
# Load the taxonomy file again (repeating due to re-upload)
taxonomy_file_path = 'MAG_taxonomy.tsv'
taxonomy_df = pd.read_csv(taxonomy_file_path, sep='\t')

# Load the T-count file
t_count_file_path = 't_count_columnC.csv'
t_count_df = pd.read_csv(t_count_file_path)
genome_to_phylum_map = dict(zip(taxonomy_df['Genome'], taxonomy_df['Phylum']))
t_count_df = pd.read_csv(t_count_file_path)


# Applying the updated mapping with sequence numbers
for i, genome in enumerate(genomes_list, start=1):
    phylum = genome_to_phylum_map.get(genome)
    if phylum and genome in t_count_df.columns:
        new_column_name = f"{phylum}_{i}"  # Append underscore and sequence number
        #print(new_column_name)
        t_count_df = t_count_df.rename(columns={genome: new_column_name})

t_count_df.to_csv(os.path.join(folder_path, 'tax_Degree.csv'), index=True)




#intersection
# List of genome names to map to Phylum
genomes_list = [
    "Ga0485158_metabat2_ours.098",
"Ga0485170_maxbin.090",
"Ga0485166_metabat2_ours.038",
"Ga0485161_metabat1.096",
"Ga0485161_metabat2_ours.167_sub",
"Ga0485172_maxbin.081_sub",
"Ga0485163_maxbin.156_sub",
"Ga0485171_metabat1.030",
"Ga0485161_metabat2_ours.042",
"Ga0485168_metabat2_ours.045_sub",
"Ga0485168_metabat2_ours.050_sub",
"Ga0485170_metabat2_ours.028",
"Ga0485170_metabat2_ours.020_sub",
"Ga0485171_metabat2_ours.185_sub",
"Ga0485165_metabat2_ours.031",
"Ga0485167_metabat2_ours.163"


]
# Load the taxonomy file again (repeating due to re-upload)
taxonomy_file_path = 'MAG_taxonomy.tsv'
taxonomy_df = pd.read_csv(taxonomy_file_path, sep='\t')

# Load the T-count file
t_count_file_path = 't_count_columnC.csv'
t_count_df = pd.read_csv(t_count_file_path)
genome_to_phylum_map = dict(zip(taxonomy_df['Genome'], taxonomy_df['Phylum']))
t_count_df = pd.read_csv(t_count_file_path)


# Applying the updated mapping with sequence numbers
for i, genome in enumerate(genomes_list, start=1):
    phylum = genome_to_phylum_map.get(genome)
    if phylum and genome in t_count_df.columns:
        new_column_name = f"{phylum}_{i}"  # Append underscore and sequence number
        #print(new_column_name)
        t_count_df = t_count_df.rename(columns={genome: new_column_name})

t_count_df.to_csv(os.path.join(folder_path, 'tax_D-B-C-E-P.csv'), index=True)



