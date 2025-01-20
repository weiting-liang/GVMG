# Set working directory
workdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
rm(list = ls())

# Install necessary packages
#install.packages("ape")
#install.packages("phangorn")
library(ape)
library(phangorn)

# Read fastANI result file
ani_data <- read.table("ANI.name.txt", header=FALSE, stringsAsFactors=FALSE)

# Extract genome names and ANI values
query_genomes <- ani_data$V1
reference_genomes <- ani_data$V2
ani_values <- ani_data$V3

# Get all unique genome names
genomes <- unique(c(query_genomes, reference_genomes))

# Initialize distance matrix
distance_matrix <- matrix(NA, nrow=length(genomes), ncol=length(genomes))
rownames(distance_matrix) <- genomes
colnames(distance_matrix) <- genomes

# Fill the distance matrix
for (i in 1:nrow(ani_data)) {
  query <- ani_data[i, 1]
  reference <- ani_data[i, 2]
  ani <- as.numeric(ani_data[i, 3])
  distance <- 1 - ani / 100  # Convert ANI value to distance
  distance_matrix[query, reference] <- distance
  distance_matrix[reference, query] <- distance
}

# Set diagonal to 0 (distance to itself is 0)
diag(distance_matrix) <- 0

# Replace NA values with maximum distance (1) to ensure matrix completeness
distance_matrix[is.na(distance_matrix)] <- 1

# Construct phylogenetic tree using maximum likelihood method
dm <- as.dist(distance_matrix)
tree <- upgma(dm)

# Ladderize the tree
tree <- ladderize(tree)

## Construct phylogenetic tree using neighbor-joining method
#tree <- nj(as.dist(distance_matrix))

# Ladderize the tree
#tree <- ladderize(nj_tree)

# Save as Newick format tree file
write.tree(tree, file="phylogenetic_tree.nwk")

# Read Newick format tree file
tree <- read.tree("phylogenetic_tree.nwk")

# Create PDF device and plot phylogenetic tree
pdf("phylogenetic_tree.pdf", width=6, height=10)
#plot(tree, type="cladogram", cex=0.6, no.margin=TRUE, align.tip.label=TRUE)
plot(tree, type="phylogram", cex=0.6, no.margin=TRUE, align.tip.label=TRUE)
#ggtree(tree) + geom_tiplab(size=3, align=TRUE) + theme_tree2()

# Add scale bar
add.scale.bar(length=0.1, lwd=2, col="black", cex=0.8)

dev.off()  # Close PDF device
