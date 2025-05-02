# code

This document provides explanations of the Python and R scripts in this repository, along with a tutorial on how to reproduce the experiment results. The primary focus of these scripts is on processing lake microbiome data and visualizing the results.

Our data is available on Zenodo: [10.5281/zenodo.14046929](https://doi.org/10.5281/zenodo.14046929)

## Scripts Overview

| Script File                          | Description                                                                                                           |
|--------------------------------------|-----------------------------------------------------------------------------------------------------------------------|
| **top10_15.py**       | Generates a directory containing files using the fraction feature selection methods for the MAGs. |
| **SE_features.R**    | Generates files using SpiecEasi Algorithm for the feature selection of the MAGs. |
| **Metagenomic_Networkfile.py**       | Generates a directory containing files using the SPIEC-EASI algorithm as well as fraction method for metagenomic data analysis. Serves as input for `Metagenomic_Network.Rmd`. |
| **Metatransciptome_Networkfile.py**  | Similar function for metatranscriptomic data, generating input files for `Metatransciptomic_Network.Rmd`.             |
| **Metagenomic_Network.Rmd**          | R Markdown file that explores connections between MAGs and environmental parameters that includes seven different feature selection methods, four permutation tests, and three percentage analyses. Outputs CarLasso models.         |
| **Metatransciptomic_Network.Rmd**    | Complements the Metagenomic analysis by focusing on metatranscriptomic data to produce comparative CarLasso models.  |


## Prerequisites

Ensure you have the following software installed before running the scripts:

- **Python**: Version 3.8 or higher recommended.
- **R**: Version 4.0 or higher recommended.
- **CARlasso**: Conditional Auto-Regressive LASSO in R: [CARlasso](https://yunyishen.github.io/CAR-LASSO/) (please follow the installation process)
- **SpiecEasi**: An R package for microbial network inference: [SpiecEasi](https://github.com/zdk123/SpiecEasi.git) (please follow the installation process)

## Running the Scripts


### Execute the Code

Run the Python and R scripts by executing these commands in your terminal:

```bash
# Generate selected features(MAGs) by using different methodologies
python3 top10_15.py
Rscript SE_features.R

# Generate network files
python3 Metagenomic_Networkfile.py
python3 Metatransciptome_Networkfile.py

# Run R Markdown analyses
Rscript -e "rmarkdown::render('Metagenomic_Network.Rmd')"
Rscript -e "rmarkdown::render('Metatransciptomic_Network.Rmd')"

```
For running these code, automatically create the following folder as input of algorithms.

| Folder Name                     | File Name                   | Description                                                                                      |
|---------------------------------|-----------------------------|--------------------------------------------------------------------------------------------------|
| **feature_selection**     | `top_10_genomes.csv`                | Using the fraction method to select top 10 OTUs based on consistent presence over time and depth from `Metagenomic_columnC.csv`.              |
|      | `top_15_genomes.csv`                | Using the fraction method to select top 15 OTUs based on consistent presence over time and depth from `Metagenomic_columnC.csv`.                         |
|      | `B.csv`       | Using SpiecEasi Algorithm to select the top 15 nodes based on betweenness centrality metrics from `Metagenomic_columnC.csv`.|
|      | `C.csv`         | Using SpiecEasi Algorithm to select the top 15 nodes based on closeness centrality metrics from `Metagenomic_columnC.csv`.    |
|     | `D.csv`            | Using SpiecEasi Algorithm to select the top 15 nodes based on degree centrality metrics from `Metagenomic_columnC.csv`.    |
|      | `P.csv`         | Using SpiecEasi Algorithm to select the top 15 nodes based on Page Rank centrality metrics from `Metagenomic_columnC.csv`. |
| **Metagenomic_Networkfile**     | `tax_10.csv`                | Features the top 10 OTUs selected using the fraction method based on consistent presence over time and depth from `Metagenomic_columnC.csv`, with names updated to their phylum level.              |
|      | `tax_15.csv`                | Features the top 15 OTUs selected using the fraction method based on consistent presence over time and depth from `Metagenomic_columnC.csv`, with names updated to their phylum level.                         |
|      | `tax_Betweenness.csv`       | Reflects betweenness centrality metrics from `Metagenomic_columnC.csv`, with selected nodes renamed to their phylum level.|
|      | `tax_Closeness.csv`         | Contains closeness centrality metrics from `Metagenomic_columnC.csv`, with selected nodes renamed to their phylum level.     |
|     | `tax_Degree.csv`            | Offers degree centrality metrics from `Metagenomic_columnC.csv` with selected nodes renamed to their phylum level.     |
|      | `tax_Page_Rank.csv`         | Offers Page Rank centrality metrics from `Metagenomic_columnC.csv` with selected nodes renamed to their phylum level. |
| **Metatransciptome_Networkfile**| `tax_10.csv`                | Features the top 10 OTUs selected using the fraction method based on consistent presence over time and depth from `metatranscriptome_columnC.csv`, with names updated to their phylum level.       |
| | `tax_15.csv`                | Features the top 15 OTUs selected using the fraction method based on consistent presence over time and depth from `metatranscriptome_columnC.csv`, with names updated to their phylum level.          |
| | `tax_Betweenness.csv`       | Reflects betweenness centrality metrics from `metatranscriptome_columnC.csv`, with selected nodes renamed to their phylum level.|
| | `tax_Closeness.csv`         | Features closeness centrality metrics from `metatranscriptome_columnC.csv`, with selected node names at the phylum level. |
|| `tax_Degree.csv`            | Offers degree centrality metrics from `metatranscriptome_columnC.csv` with selected nodes renamed to their phylum level. |
| | `tax_Page_Rank.csv`         | Offers Page Rank centrality metrics from `metatranscriptome_columnC.csv` with selected nodes renamed to their phylum level. |



For running **Metagenomic_Network.Rmd**, **Metatransciptomic_Network.Rmd**, the following folders are required: `Metagenomic_Networkfile`, `Metagenomic_Networkfile`.





