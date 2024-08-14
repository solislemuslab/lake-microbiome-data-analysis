# code

This document provides explanations of the Python and R scripts in this repository, along with a tutorial on how to reproduce the experiment results. The primary focus of these scripts is on processing lake microbiome data and visualizing the results.

## Scripts Overview

| Script File                          | Description                                                                                                           |
|--------------------------------------|-----------------------------------------------------------------------------------------------------------------------|
| **Metagenomic_Networkfile.py**       | Generates a directory containing files using the SPIEC-EASI algorithm as well as fraction method for metagenomic data analysis. Serves as input for `Metagenomic_Network.Rmd`. |
| **Metatransciptome_Networkfile.py**  | Similar function for metatranscriptomic data, generating input files for `Metatransciptomic_Network.Rmd`.             |
| **lake-microbiome-figures.py**       | Analyzes microbial abundance in varied environmental conditions and generates visualizations.                         |
| **Metagenomic_Network.Rmd**          | R Markdown file that explores connections between MAGs and environmental parameters that includes seven different feature selection methods, four permutation tests, and three percentage analyses. Outputs CarLasso models.         |
| **Metatransciptomic_Network.Rmd**    | Complements the Metagenomic analysis by focusing on metatranscriptomic data to produce comparative CarLasso models.  |
| **Car-figures.R**                    | Uses the CARlasso model to identify and visualize connections between MAGs and environmental parameters.               |


## Prerequisites

Ensure you have the following software installed before running the scripts:

- **Python**: Version 3.8 or higher recommended.
- **R**: Version 4.0 or higher recommended.
- CARlasso: Conditional Auto-Regressive LASSO in R: [CARlasso](https://yunyishen.github.io/CAR-LASSO/) (please follow the installation process)
- SpiecEasi: An R package for microbial network inference: [SpiecEasi](https://github.com/zdk123/SpiecEasi.git) (please follow the installation process)

## Running the Scripts


### Execute the Code

Run the Python and R scripts by executing these commands in your terminal:

```bash
# Generate network files
python3 Metagenomic_Networkfile.py
python3 Metatransciptome_Networkfile.py

# Run R Markdown analyses
Rscript -e "rmarkdown::render('Metagenomic_Network.Rmd')"
Rscript -e "rmarkdown::render('Metatransciptomic_Network.Rmd')"

# Generate figures
python3 lake-microbiome-figures.py
Rscript Car-figures.R

```
For running these code, automatically create the following folder as input of algorithms.

| Folder Name                     | File Name                   | Description                                                                                      |
|---------------------------------|-----------------------------|--------------------------------------------------------------------------------------------------|
| **Metagenomic_Networkfile**     | `tax_10.csv`                | Features the top 10 OTUs selected using the fraction method based on consistent presence over time and depth from `Metagenomic_columnC.csv`, with names updated to their phylum level.              |
|      | `tax_15.csv`                | Features the top 15 OTUs selected using the fraction method based on consistent presence over time and depth from `Metagenomic_columnC.csv`, with names updated to their phylum level.                         |
|      | `tax_Betweenness.csv`       | Reflects betweenness centrality metrics from `Metagenomic_columnC.csv`, with selected nodes renamed to their phylum level.|
|      | `tax_Closeness.csv`         | Contains closeness centrality metrics from `Metagenomic_columnC.csv`, with selected nodes renamed to their phylum level.     |
|      | `tax_D-B-C-E-P.csv`         | Contains an intersection of "degree", "betweenness", "closeness", "eigenvector centrality", and "PageRank" metrics from `Metagenomic_columnC.csv`, but with node names adjusted to their phylum level. |
|     | `tax_Degree.csv`            | Offers degree centrality metrics from `Metagenomic_columnC.csv` with selected nodes renamed to their phylum level.     |
|      | `tax_Page_Rank.csv`         | Offers Page Rank centrality metrics from `Metagenomic_columnC.csv` with selected nodes renamed to their phylum level. |
| **Metatransciptome_Networkfile**| `tax_10.csv`                | Features the top 10 OTUs selected using the fraction method based on consistent presence over time and depth from `metatranscriptome_columnC.csv`, with names updated to their phylum level.       |
| | `tax_15.csv`                | Features the top 15 OTUs selected using the fraction method based on consistent presence over time and depth from `metatranscriptome_columnC.csv`, with names updated to their phylum level.          |
| | `tax_Betweenness.csv`       | Reflects betweenness centrality metrics from `metatranscriptome_columnC.csv`, with selected nodes renamed to their phylum level.|
| | `tax_Closeness.csv`         | Features closeness centrality metrics from `metatranscriptome_columnC.csv`, with selected node names at the phylum level. |
|| `tax_Degree.csv`            | Offers degree centrality metrics from `metatranscriptome_columnC.csv` with selected nodes renamed to their phylum level. |
| | `tax_Page_Rank.csv`         | Offers Page Rank centrality metrics from `metatranscriptome_columnC.csv` with selected nodes renamed to their phylum level. |
| **OutputFigures**               | -                           | Contains all figures generated by the scripts, showcasing visual representations of data analyses conducted in the research paper. For more detailed files generated inside the OutputFigures folder, please refer to the `output` folder.|




For running **Metagenomic_Network.Rmd**, **Metatransciptomic_Network.Rmd**, **Car-figures.R**, the following folders are required: `Metagenomic_Networkfile`, `Metagenomic_Networkfile`.



