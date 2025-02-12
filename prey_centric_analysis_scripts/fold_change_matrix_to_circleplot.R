# Load necessary libraries
library(tidygraph)    # For graph manipulation
library(tidyverse)    # For data wrangling and visualization
library(ggraph)       # For network visualization using ggplot
library(dendextend)   # For dendrogram manipulation
library(igraph)       # For network analysis

# Function to cut a dendrogram at a specific height and assign unique cluster IDs
cut_df <- function(dendrogram, height, c){
  # Cut the dendrogram at the specified height
  cd <- cutree(dendrogram, h = height) %>% as.data.frame()
  cd$ID <- row.names(cd)
  cd <- cd %>% as_tibble()
  
  # Generate a column name based on the cut value
  colname <- paste("cut", as.character(round(c, 2)), sep = "_")
  names(cd) <- c(colname, "ID")  # Rename columns for clarity
  return(cd)
}

# Function to cut the dendrogram at multiple heights
cut_dend <- function(dendrogram, cuts){
  ht <- max(get_nodes_attr(dendrogram, "height"))  # Get maximum height of the dendrogram
  cut_clusters <- data.frame(ID = as.character())  # Initialize an empty dataframe for clusters
  
  # Loop through each cut point and merge the results
  for (c in cuts){
    cut_clusters <- merge(cut_clusters, cut_df(dendrogram, c*ht, c), all=TRUE)
  }
  return(cut_clusters)  # Return the clusters for each cut
}

# Specify the path to the network matrix file (replace with your actual file path)
network_matrix_file = 'path_to_fold_change_matrix_file.csv'  # Placeholder path
network_matrix <- read.csv(network_matrix_file)  # Read the network matrix into a data frame

# Create a dendrogram from the hierarchical clustering of the network matrix
d_iris <- as.dendrogram(hclust(dist(network_matrix))) 
plot(d_iris)  # Plot the dendrogram

# Specify the path to the annotation file (replace with your actual file path)
annotation_file <- 'path_to_annotation_file.csv'  # Placeholder path
annotation <- read.csv(annotation_file)  # Read the annotation data

# Select relevant columns from the annotation file
annotation <- annotation[, c("Protein.ID", "characterization", "annotate_0.4")]

# Merge annotation data with the network matrix based on "Protein.ID"
network_matrix <- left_join(network_matrix, annotation, by = c("Protein.ID" = "Protein.ID"))

# Define the heights at which to cut the dendrogram
cut_clusters <- cut_dend(d_iris, c(0.6, 0.5, 0.4))  # Modify the cuts as needed

# Create a unique cluster ID for each cluster at each height
clusters_uniqued <- cut_clusters  %>%
  gather(clusterset, clusternum, -ID) %>%
  mutate(clusterid = paste0(clusternum, clusterset)) %>%
  select(-clusternum) %>%
  spread(clusterset, clusterid)

# Assign cluster information to the original data frame using left_join
network_matrix_with_clusters <- network_matrix %>% 
  rownames_to_column("Protein ID") %>% 
  left_join(cut_clusters %>% gather(clusterset, clusternum, -ID) %>% 
              mutate(clusterid = paste0(clusternum, clusterset)) %>% 
              select(-clusternum) %>% 
              spread(clusterset, clusterid), 
            by = c("Protein ID" = "ID"))

# The following block of code generates edges for each cluster level
# These edges represent the hierarchical relationships between clusters at different levels

# Create a link column for all nodes (using "origin" as the placeholder)
clusters_uniqued$link <- "origin" 

# Create edges between clusters at different levels (0.6, 0.5, and 0.4 cuts)
edges_level0_1 = clusters_uniqued %>% select(link, cut_0.6) %>% unique %>% rename(from = link, to = cut_0.6)
edges_level2_3 = clusters_uniqued %>% select( cut_0.6, cut_0.5) %>% unique %>% rename(from = cut_0.6, to = cut_0.5)
edges_level3_4 = clusters_uniqued %>% select( cut_0.5, cut_0.4) %>% unique %>% rename(from = cut_0.5, to = cut_0.4)
edges_level5_final = clusters_uniqued %>% select( cut_0.4, ID) %>% unique %>% rename(from = cut_0.4, to = ID)

# Combine all edges into a single list
edge_list = rbind(edges_level0_1, edges_level2_3, edges_level3_4, edges_level5_final)

# Create a graph object from the edge list
mygraph <- as_tbl_graph(edge_list)

# Add the row names (Protein ID) from the network matrix as a new column
csm_plain <- network_matrix_with_clusters %>% mutate(name = rownames(network_matrix_with_clusters))

# Add species labels to the nodes in the graph
mygraph <- mygraph %>% activate("nodes") %>% 
  left_join(csm_plain, by = "name")

# Create a circular plot of the network with specified colors
csm_circleplot <- ggraph(mygraph, layout = 'circlepack') + 
  geom_node_circle(aes(fill = characterization), linewidth = 1) +  # Visualize nodes with color based on characterization
  theme_void() +  # Remove axis and gridlines for a clean plot
  scale_fill_brewer(palette = 'Paired') +  # Choose color palette
  guides(fill = guide_legend(title = "Characterized state")) +  # Add legend for color scale
  theme(text = element_text(size = 25),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())  # Customize the theme for the plot

# Display the plot
csm_circleplot


# Save the circular plot as a PDF file (replace with your actual file path)
# ggsave(csm_circleplot,
#        filename = 'path_to_output_file.pdf',  # Placeholder path
#        width = 12, height = 7, units = 'in', device = cairo_pdf)  # Save plot to file
