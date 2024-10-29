# Create a directory
dir.create("multi20_list")

# Read the Mash distance file
mash_data <- read.table("mash_distances.tsv", header=FALSE)

# Assume the columns follow this format:
# ReferenceID, QueryID, Distance, PValue, SharedHashes
colnames(mash_data) <- c("ReferenceID", "QueryID", "Distance", "PValue", "SharedHashes")

# Trim ReferenceID and QueryID to keep only the characters before the first dot after the last slash
mash_data$ReferenceID <- sub(".*/", "", mash_data$ReferenceID)
mash_data$QueryID <- sub(".*/", "", mash_data$QueryID)
mash_data$ReferenceID <- sub("\\..*$", "", mash_data$ReferenceID)
mash_data$QueryID <- sub("\\..*$", "", mash_data$QueryID)

# Get a list of all unique samples
all_samples <- unique(c(mash_data$ReferenceID, mash_data$QueryID))

# Initialize the total table data frame
total_table <- data.frame(SampleID=character(), Nearest20=character(), stringsAsFactors=FALSE)

# Operate on each sample
for(sample in all_samples) {
  # Filter the comparison results related to the current sample
  related_data <- subset(mash_data, ReferenceID == sample)

  # Sort by Mash distance and select the nearest 20 samples
  top20 <- head(related_data[order(related_data$Distance), ], 20)

  # Extract the IDs of the nearest 20 samples
  nearest_sample_ids <- top20$QueryID

  # Output to file
  file_name <- paste0("./multi20_list/", sample, ".multi20.bwa.txt")
  writeLines(nearest_sample_ids, file_name)

  # Add the current sample and the corresponding nearest 20 sample IDs to the total table
  total_table <- rbind(total_table, data.frame(SampleID=sample, Nearest20=paste(nearest_sample_ids, collapse=",")))
}

# Output the total table to a TSV file
write.table(total_table, "total_nearest_samples.tsv", sep="\t", row.names=FALSE, quote=FALSE)
