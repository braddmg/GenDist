#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(ggraph))
suppressPackageStartupMessages(library(tidygraph))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(viridis))

# Define command-line options
option_list <- list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="Input TSV file containing gene presence data", metavar="file"),
  make_option(c("-t", "--threshold"), type="numeric", default=0.5,
              help="Jaccard distance threshold for network edges (default: 0.5)", metavar="threshold"),
  make_option(c("-s", "--scale_node"), action="store_true", default=FALSE,
              help="Scale node size by gene prevalence"),
  make_option(c("-o", "--output"), type="character", default="network_plot.pdf",
              help="Output filename for the network plot (default: network_plot.pdf)", metavar="output"),
  make_option(c("-p", "--prevalence"), type="numeric", default=0,
              help="Prevalence threshold (default: 0, no filtering). Genes must be present in at least X fraction of sequences")
)

# Parse command-line arguments
opt <- parse_args(OptionParser(option_list=option_list))

# Check if input file is provided
if (is.null(opt$file)) {
  stop("Error: You must provide an input file using -f or --file")
}

# Load data
data <- read_tsv(opt$file, col_types = cols())

# Create presence-absence matrix
gene_matrix <- data %>%
  group_by(GENE, SEQUENCE) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = SEQUENCE, values_from = count, values_fill = list(count = 0))

# Convert to binary matrix
binary_matrix <- as.matrix(gene_matrix[, -1])
rownames(binary_matrix) <- gene_matrix$GENE
binary_matrix[binary_matrix > 1] <- 1  # Convert counts to binary presence/absence

# Filter by prevalence threshold
if (opt$prevalence > 0) {
  gene_prevalence <- rowSums(binary_matrix) / ncol(binary_matrix)
  binary_matrix <- binary_matrix[gene_prevalence >= opt$prevalence, , drop = FALSE]
}

# Compute Jaccard distances
compute_jaccard_distances <- function(matrix) {
  gene_list <- rownames(matrix)
  results <- expand.grid(Gene1 = gene_list, Gene2 = gene_list, stringsAsFactors = FALSE)
  results <- results[results$Gene1 < results$Gene2, ]
  results$Jaccard_Distance <- apply(results, 1, function(row) {
    A <- matrix[row["Gene1"], ] == 1
    B <- matrix[row["Gene2"], ] == 1
    1 - (sum(A & B) / sum(A | B))
  })
  return(results[order(results$Jaccard_Distance), ])
}

# Compute Jaccard distances
jaccard_results <- compute_jaccard_distances(binary_matrix)
write_tsv(jaccard_results, "jaccard_distance_results.tsv")

# Filter based on threshold
filtered_jaccard <- subset(jaccard_results, Jaccard_Distance < opt$threshold)
filtered_jaccard$Jaccard_Similarity <- 1 - filtered_jaccard$Jaccard_Distance

# Compute gene prevalence (ensure it's a named vector)
gene_prevalence <- rowSums(binary_matrix)
names(gene_prevalence) <- rownames(binary_matrix)

# Create network graph
g <- graph_from_data_frame(filtered_jaccard[, c("Gene1", "Gene2")], directed = FALSE)
E(g)$weight <- filtered_jaccard$Jaccard_Similarity

# Fix the node size issue
V(g)$size <- if (opt$scale_node) {
  ifelse(V(g)$name %in% names(gene_prevalence), gene_prevalence[V(g)$name], 5)
} else {
  5
}

# Convert to tidygraph
g_tidy <- as_tbl_graph(g)

# Generate network plot
p <- ggraph(g_tidy, layout = "fr") +
  geom_edge_link(aes(width = weight), color = "gray70", alpha = 0.8) +
  geom_node_point(aes(size = size), color = "skyblue", alpha = 0.9) +
  geom_node_text(aes(label = name), size = 4, vjust = 1.5, hjust = 0.5) +
  scale_size_continuous(range = c(6, 20)) +
  scale_edge_width_continuous(range = c(1, 6)) +
  theme_void() +
  theme(plot.margin = margin(10, 10, 10, 10))

# Save the plot as a PDF
ggsave(opt$output, plot = p, width = 12, height = 10, dpi = 300, device = "pdf")

message("Gene presence-absence matrix and network plot saved!")
