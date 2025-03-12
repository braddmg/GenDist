library(optparse)
library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)
library(viridis)

# Define command-line options
option_list <- list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Input CSV file containing the gene presence-absence matrix", metavar="file"),
  make_option(c("-t", "--threshold"), type="numeric", default=0.5, 
              help="Jaccard distance threshold for network edges (default: 0.5)", metavar="threshold"),
  make_option(c("-p", "--prevalence"), action="store_true", default=FALSE, 
              help="Scale node size by gene prevalence"),
  make_option(c("-o", "--output"), type="character", default="network_plot.pdf", 
              help="Output filename for the plot (default: network_plot.pdf)", metavar="output")
)

# Parse command-line arguments
opt <- parse_args(OptionParser(option_list=option_list))

# Load input data
data <- read.csv(opt$file, row.names = 1)
binary_matrix <- as.matrix(data)

# Ensure binary values
binary_matrix[binary_matrix > 1] <- 1
binary_matrix[binary_matrix < 0] <- 0

# Function to compute Jaccard distances
compute_jaccard_distances <- function(matrix) {
  gene_list <- rownames(matrix)
  results <- data.frame(Gene1 = character(), Gene2 = character(), Jaccard_Distance = numeric())
  
  for (i in 1:(length(gene_list) - 1)) {
    for (j in (i + 1):length(gene_list)) {
      gene1 <- gene_list[i]
      gene2 <- gene_list[j]
      
      A <- matrix[gene1, ] == 1
      B <- matrix[gene2, ] == 1
      jaccard_dist <- 1 - (sum(A & B) / sum(A | B))
      
      results <- rbind(results, data.frame(Gene1 = gene1, Gene2 = gene2, Jaccard_Distance = jaccard_dist))
    }
  }
  
  return(results[order(results$Jaccard_Distance), ])
}

# Compute Jaccard distances
jaccard_results <- compute_jaccard_distances(binary_matrix)
write.csv(jaccard_results, "jaccard_distance_results.csv", row.names = FALSE)

# Filter based on threshold
filtered_jaccard <- subset(jaccard_results, Jaccard_Distance < opt$threshold)
filtered_jaccard$Jaccard_Similarity <- 1 - filtered_jaccard$Jaccard_Distance

# Compute gene prevalence
gene_prevalence <- rowSums(binary_matrix)

# Create network graph
g <- graph_from_data_frame(filtered_jaccard[, c("Gene1", "Gene2")], directed = FALSE)
E(g)$weight <- filtered_jaccard$Jaccard_Similarity
V(g)$size <- if (opt$prevalence) gene_prevalence[V(g)$name] else 5  # Default size if not scaling

# Convert to tidygraph
g_tidy <- as_tbl_graph(g)

# Generate wider network plot
p <- ggraph(g_tidy, layout = "fr") +
  geom_edge_link(aes(width = weight), color = "gray70", alpha = 0.8) +
  geom_node_point(aes(size = size), color = "skyblue", alpha = 0.9) +
  geom_node_text(aes(label = name), size = 4, vjust = 1.5, hjust = 0.5) +
  scale_size_continuous(range = c(6, 20)) +
  scale_edge_width_continuous(range = c(1, 6)) +
  theme_void() +
  theme(plot.margin = margin(10, 10, 10, 10))

# Save the plot as a PDF with a wider aspect ratio
ggsave(opt$output, plot = p, width = 12, height = 10, dpi = 300, device = "pdf")

message("Plot saved as: ", opt$output)