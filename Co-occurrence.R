#!/usr/bin/env Rscript

# Load required libraries
library(optparse)
library(igraph)
library(ggplot2)
library(viridis)
library(ggraph)
library(tidygraph)

# Function to compute Jaccard distances
compute_jaccard_distances <- function(binary_matrix) {
  gene_list <- rownames(binary_matrix)
  results <- data.frame(Gene1 = character(), Gene2 = character(), Jaccard_Distance = numeric())
  
  for (i in 1:(length(gene_list) - 1)) {
    for (j in (i + 1):length(gene_list)) {
      gene1 <- gene_list[i]
      gene2 <- gene_list[j]
      
      # Compute Jaccard Distance
      A <- binary_matrix[gene1, ] == 1
      B <- binary_matrix[gene2, ] == 1
      jaccard_dist <- 1 - (sum(A & B) / sum(A | B))
      
      # Store results
      results <- rbind(results, data.frame(Gene1 = gene1, Gene2 = gene2, Jaccard_Distance = jaccard_dist))
    }
  }
  
  return(results[order(results$Jaccard_Distance), ])
}

# Function to generate network graph
create_network_graph <- function(binary_matrix, jaccard_threshold, use_prevalence) {
  jaccard_results <- compute_jaccard_distances(binary_matrix)
  filtered_jaccard <- subset(jaccard_results, Jaccard_Distance < jaccard_threshold)
  filtered_jaccard$Jaccard_Similarity <- 1 - filtered_jaccard$Jaccard_Distance
  
  # Create network
  g <- graph_from_data_frame(filtered_jaccard[, c("Gene1", "Gene2")], directed = FALSE)
  E(g)$weight <- filtered_jaccard$Jaccard_Similarity
  
  # Compute gene prevalence if required
  if (use_prevalence) {
    gene_prevalence <- rowSums(binary_matrix)
    V(g)$size <- log10(gene_prevalence[V(g)$name] + 1)
  } else {
    V(g)$size <- 5  # Default size
  }
  
  # Convert to tbl_graph and plot
  g_tidy <- as_tbl_graph(g)
  ggraph(g_tidy, layout = "fr") +
    geom_edge_link(aes(width = weight), color = "gray70", alpha = 0.8) +
    geom_node_point(aes(size = size), color = "skyblue", alpha = 0.9) +
    geom_node_text(aes(label = name), size = 4, vjust = 1.5, hjust = 0.5, fontface = "italic") +
    scale_size_continuous(range = c(3, 12)) +
    scale_edge_width_continuous(range = c(0.5, 5)) +
    theme_void()
}

# Define command-line options
option_list <- list(
  make_option(c("-f", "--file"), type = "character", help = "Input CSV file (gene presence/absence matrix)", metavar = "FILE"),
  make_option(c("-t", "--threshold"), type = "double", default = 0.5, help = "Jaccard distance threshold (default: 0.5)", metavar = "VALUE"),
  make_option(c("-p", "--prevalence"), action = "store_true", default = FALSE, help = "Scale node size by gene prevalence")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Display help if no file is provided
if (is.null(opt$file)) {
  print_help(opt_parser)
  quit(status = 1)
}

# Load data
Data <- read.csv(opt$file, row.names = 1)
binary_matrix <- as.matrix(Data)
binary_matrix[binary_matrix > 1] <- 1
binary_matrix[binary_matrix < 0] <- 0

# Generate network graph
create_network_graph(binary_matrix, opt$threshold, opt$prevalence)
