###################################### 
# Load Required Libraries 
######################################
library(SpiecEasi)
library(Matrix)
library(dplyr)
library(igraph)

######################################
# Set Seed for Reproducibility
######################################
set.seed(100)
######################################
# Load and Preprocess Data
######################################
# Read the ASV table from a CSV file
asvtable <- read.csv("Metagenomic_columnC.csv", row.names = 1)
# Remove the last 10 columns from the ASV table
# The last 10 columns are assumed to be metadata or unnecessary for the analysis
asvtable <- asvtable[, 1:(ncol(asvtable) - 10)]
######################################
# Apply SpiecEasi to Construct a Network
######################################
# Run SpiecEasi using the glasso method
se <- spiec.easi(as.matrix(asvtable), method = 'glasso', lambda.min.ratio = 1e-2, nlambda = 100)
# Extract the adjacency matrix of the network
se_adj <- se$refit$stars
# Convert the adjacency matrix to a format suitable for graph analysis
SE_glasso <- as.matrix(se_adj)
######################################
# Network Analysis
######################################
# Create an igraph object from the adjacency matrix
graph <- graph_from_adjacency_matrix(SE_glasso, mode = "undirected", diag = FALSE)
# Calculate degree centrality for each node in the graph
degree <- degree(graph)
# Calculate betweenness centrality for each node in the graph
betweenness <- betweenness(graph)
# Calculate closeness centrality for each node in the graph
closeness <- closeness(graph)
# Calculate eigenvalue centrality for each node in the graph
evcent <- evcent(graph)$vector
# Calculate page.rank centrality for each node in the graph
page_rank <- page.rank(graph)$vector

# Sort nodes based on  centrality in decreasing order
sorted_degree <- order(degree, decreasing = TRUE)
sorted_betweenness <- order(betweenness, decreasing = TRUE)
sorted_closeness <- order(closeness, decreasing = TRUE)
sorted_evcent <- order(evcent, decreasing = TRUE)
sorted_page_rank <- order(page_rank, decreasing = TRUE)
# Define the proportion of top nodes to be selected

p <- 0.2

# Calculate the top p percent of nodes based on centrality measures
top_degree <- sorted_degree[1:round(p * dim(asvtable)[2])]
top_betweenness <- sorted_betweenness[1:round(p*dim(asvtable)[2])]
top_closeness <- sorted_closeness[1:round(p*dim(asvtable)[2])]
top_evcent <- sorted_evcent[1:round(p*dim(asvtable)[2])]
top_page_rank <- sorted_page_rank[1:round(p*dim(asvtable)[2])]

######################################
# Output Results
######################################
# Save the names of the top 15 nodes (based on degree centrality) to a CSV file
write.csv(colnames(asvtable)[top_degree[1:15]], file = "D.csv")
write.csv(colnames(asvtable)[top_betweenness[1:15]],file="B.csv")
write.csv(colnames(asvtable)[top_closeness[1:15]],file="C.csv")
write.csv(colnames(asvtable)[top_page_rank[1:15]],file="P.csv")
