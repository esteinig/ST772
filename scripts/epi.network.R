# This script computes the connected components in a thresholded graph, where
# the threshold corresponds to a maximum SNP distance between isolates, e.g. based on
# nucleotide substitution rate (4 SNPs per year) and use for:
# Supplementary File 2 and Supplementary Figure 6

library(ggplot2)
library(network)
library(sna)
library(ggnet)
library(igraph)
library(intergraph)

# Check names is FALSE since row isolate IDs start with number:

distance = as.matrix(read.table("snippy.hamming.counts.tab", sep="\t", header=T, row.names=1, check.names=F))

# Subset distance matrix by reducing the veterinary cluster to one rperesentative isolate:

meta = read.csv("meta.csv", header=T, stringsAsFactors=F)

meta = meta[meta$Include == "Yes",]

# Distance matrix is Hamming distance from Snippy core alignment, replace 'DAR_4145' with 'Reference':

meta[meta$ID == "DAR_4145", "ID"] <- "Reference"

distance = distance[rownames(distance) %in% meta$ID, colnames(distance) %in% meta$ID]

# Dimensions should be 302, 302

dim(distance)

# Define subsetting outbreak function to generate network plot:

get.outbreaks <- function(snp.distances, snp.threshold=4){
  
  
  # Do not exclude clones i.e. where distance == 0 - replace with dummy value and so that
  # weighted edges can be read by adjacency matrix function:
  
  # Includes diagnonal so reset diagonal after:
  snp.distances[snp.distances == 0] <- 0.1
  diag(snp.distances) <- 0
  
  # Subset distance matrix where pairs <= snp.threshold,
  # based on substitution rate of 1.613707e-06 (Gubbins, LSD)
  
  snp.distances[snp.distances > snp.threshold] <- 0
  
  net <- graph_from_adjacency_matrix(snp.distances, mode="undirected", weighted=TRUE, add.rownames=T)
  
  # Reset clone edge weights so they are displayed correctly:
  E(net)$weight[E(net)$weight == 0.1] <- 0
  
  g <- delete.vertices(simplify(net), degree(net)==0)
  
  return(g)
  
}

g <- get.outbreaks(distance)

# Get epidemiological Link status for nodes in network:

node.links <- meta[match(V(g)$name, meta$ID), "Link"]
names(node.links) <- V(g)$name

# Translate to color mapping:

color.link <- c("Yes" = "#762a83", "No" = "#5aae61", "Unknown" = "#d3d3d3")

V(g)$color <- as.character(plyr::revalue(as.factor(node.links), color.link))

pdf(file="epi.network.pdf", width=41, height=41)

ggnet2(g, label = TRUE, edge.label = "weight", label.size=8, node.size=15, node.color="color", edge.size = 1.5, edge.label.size=3)

dev.off()