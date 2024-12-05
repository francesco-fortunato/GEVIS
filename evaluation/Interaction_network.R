# Load necessary libraries
library(multiMiR)
library(igraph)
library(RColorBrewer)

# Define the miRNA_target function
miRNA_target <- function(miRNAs) {
  # Retrieve the target genes using multiMiR
  targets_multiMiR <- get_multimir(mirna = miRNAs, table = "predicted")

  # Access the data slot
  data_slot <- targets_multiMiR@data

  # Convert to data frame
  targets_df <- as.data.frame(data_slot)

  # Check if targets_df is empty
  if (nrow(targets_df) == 0) {
    warning("No targets found for the specified miRNAs.")
    return(NULL)
  }

  # Create a data frame for edges (miRNA-gene_symbol pairs)
  edges_df <- data.frame(
    miRNA = targets_df$mature_mirna_id,  # Assuming 'mature_mirna_id' holds the miRNA IDs
    gene_symbol = targets_df$target_symbol,  # Use 'target_symbol' for gene symbols
    stringsAsFactors = FALSE
  )

  edges_df <- unique(edges_df)

  # Create the graph using igraph
  g <- graph_from_data_frame(edges_df, directed = TRUE)

  # Get the list of unique miRNAs
  miRNAs_unique <- unique(edges_df$miRNA)

  # Create vertex attributes for color and size
  V(g)$color <- ifelse(V(g)$name %in% miRNAs_unique, "blue", "grey")
  V(g)$size <- ifelse(V(g)$name %in% miRNAs_unique, 30, 10)  # miRNAs are larger

  # Plotting the network
  pdf("miRNA_target_network.pdf", width = 10, height = 10)  # Save the plot as a PDF
  plot(g,
       vertex.label = NA,  # Remove labels
       vertex.color = V(g)$color,
       vertex.size = V(g)$size,
       edge.arrow.size = 0.5,
       main = "miRNA-Target Gene Network"
  )
  dev.off()  # Close the PDF device

  # View the first few rows of the edges data frame (miRNA-gene pairs)
  return(edges_df)
}

# Example miRNAs
miRNAs <- c("hsa-miR-21-5p", "hsa-miR-155-5p")

# Call the function to generate the network and save as a PDF
miRNA_target(miRNAs)
