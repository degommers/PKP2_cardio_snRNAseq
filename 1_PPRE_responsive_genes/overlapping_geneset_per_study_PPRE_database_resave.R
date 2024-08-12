install.packages("VennDiagram")
library(VennDiagram)
library(tidyverse)
library(readr)
library(readxl)
library(RColorBrewer)
library(pheatmap)

# External variables
myCol <- brewer.pal(3, "Pastel2")


#### Load PRRE database ####
database_PPRE <- read_delim("Database/PPRE_DATABASE_17042024.csv", delim = ",")

#### Load PPRE database meta data ####
PPRE_metadata <- read_excel("Database/PPRE_database_metadata_annotation.xlsx")

#### Load FAO genes ####
FAO_genes <- read_delim("Database/FAO_genes_ensembl_ID.txt", delim = "\t", col_names = "FAO")
gene_list <- c("ABCB11", "ACACB", "SLC25A17", "DECR1",  "CPT1A",  "EHHADH",  "ACADS",  "SESN2",
               "HADH", "HADHB", "HADHA", "SLC27A2", "FABP1", "ETFDH", 
               "MTLN", "ADIPOQ", "PPARA", "ACAD11", "DECR2", "FABP3", "PPARG", 
               "PPARD", "LEP", "RXRA", "RXRB", "RXRG")


# Add FAO genes to database
database_PPRE_FAO <- database_PPRE %>%
  mutate(
    FAO = ifelse(ensembl_gene_id %in% FAO_genes$FAO, 1, 0),
    gene_list = ifelse(gene_name %in% gene_list, 1, 0)
    
  )

# Function to filter a given column
filter_ensembl_ids <- function(database, column_name) {
  filtered_ids <- database$ensembl_gene_id[database[, column_name] == 1]
  return(filtered_ids)
}

#### Extract the column names of the database ####
# Picks from column 8 to the end only the values of the presence matrix of the database
cols <- colnames(database_PPRE_FAO)[8:41]

# Initialize an empty list to store filtered IDs for each column
filtered_ids_list <- list()

#### Making a list of lists including the gene IDs that are 'present' ####
# present ==1
for (col in cols) {
  # Filter the IDs for the current column
  filtered_ids <- filter_ensembl_ids(database_PPRE_FAO, col)
  
  # Store the filtered IDs in the list with the column name as key
  filtered_ids_list[[col]] <- filtered_ids
  
}

# Initialize an empty matrix to store the presence/absence of Ensembl gene IDs
ids_matrix <- matrix(0, nrow = length(database_PPRE_FAO$ensembl_gene_id), ncol = length(filtered_ids_list))

# Populate the matrix with 1 where the Ensembl gene ID is present in the corresponding category
for (i in seq_along(filtered_ids_list)) {
  ids_matrix[, i] <- as.integer(database_PPRE_FAO$ensembl_gene_id %in% filtered_ids_list[[i]])
}

#### Set row names and column names to Ensembl gene IDs and the column names of the matrix ####
rownames(ids_matrix) <- database_PPRE_FAO$ensembl_gene_id

# Set column names to category names
colnames(ids_matrix) <- names(filtered_ids_list)

#### Combine the ATACseq replicates per experiment/modulation ####
# Create a list of matrices
matrices <- list(
  atac_matrix_control = ids_matrix[, 13:15],
  atac_matrix_GSK0660 = ids_matrix[, 16:18],
  atac_matrix_GW0742 = ids_matrix[, 19:21],
  atac_matrix_LCFA = ids_matrix[, 22:24],
  atac_matrix_LCFA_GSK0660 = ids_matrix[, 25:27],
  atac_matrix_LCFA_GW0742 = ids_matrix[, 28:30]
)

# Loop through the list
for (i in seq_along(matrices)) {
  # Extract the matrix from the list
  matrix <- matrices[[i]]
  
  matrix <- as.data.frame(matrix)
  
  # Create a new combined column initialized with zeros
  matrix$combined_replicates <- 0
  
  # Check if there are any 1s in each row
  row_sums <- rowSums(matrix[, 1:3]) 
  
  # Assign 1 to the combined column if any row sum is greater than 0
  matrix$combined_replicates <- ifelse(row_sums > 2, 1, 0)
  
  # Assign the modified matrix back to the list
  matrices[[i]] <- matrix
}

# Create an empty matrix to store the combined data
combined_matrix <- matrix(0, nrow = nrow(ids_matrix), ncol = 0)

# Iterate over each matrix
for (matrix_name in names(matrices)) {
  # Extract the "combined_replicates" column from the current matrix
  combined_column <- matrices[[matrix_name]]$combined_replicates
  
  # Add the combined column to the combined matrix
  combined_matrix <- cbind(combined_matrix, combined_column)
}

# Set row names of the combined matrix
rownames(combined_matrix) <- rownames(ids_matrix)
colnames(combined_matrix) <- names(matrices)


# Replace combined matrix columns in the ids_matrix
ids_matrix <- ids_matrix[,c(1:12,31:34)]
ids_matrix_final <- cbind(ids_matrix, combined_matrix)


#### Filter matrix on given gene list ####
# Filter row names based on the list
FAO_ids <- rownames(ids_matrix_final)[rownames(ids_matrix_final) %in% filtered_ids_list[["FAO"]]]
genelist_ids <- rownames(ids_matrix_final)[rownames(ids_matrix_final) %in% filtered_ids_list[["gene_list"]]]
# Filter the matrix based on the filtered row names
filtered_matrix_fao_genes <- ids_matrix_final[FAO_ids, , drop = FALSE]
filtered_matrix_genelist_genes <- ids_matrix_final[genelist_ids, , drop = FALSE]


#### Check presence of columns in matrix with metadata for annotation ####
studies <- sort(colnames(filtered_matrix_genelist_genes))
metadata_studies <- sort(PPRE_metadata[8:28,]$Col_name)
studies %in% metadata_studies

modify_matrix_values <- function(matrix, metadata, column_name, study_type_column) {
  # Loop over the columns of the matrix
  for (col in colnames(matrix)) {
    # Skip the "gene_list" column
    if (col == "gene_list") {
      next
    }
    
    print(paste("Processing column:", col))
    
    # Get the study type for the current column
    study_type <- metadata[metadata[[column_name]] == col, study_type_column]
    
    # Check if study_type is not found or if it's the last column
    if (length(study_type) == 0 || col == tail(colnames(matrix), 1)) {
      if (length(study_type) == 0) {
        warning(paste("Study type for column", col, "not found in metadata"))
      }
      if (col == tail(colnames(matrix), 1)) {
        print("No more columns to process, exiting loop")
      } else {
        print("Skipping to next column")
      }
      next
    }
    
    # Print the study type for debugging
    print(paste("Study type for column", col, "is:", study_type))
    
    # Apply modification based on study type, only if the value is not 0
    if (length(study_type) > 0 && any(matrix[, col] != 0)) { # Check if study type is not empty and any non-zero values exist
      if (study_type == "ChIP-seq") {
        matrix[matrix[, col] != 0, col] <- matrix[matrix[, col] != 0, col] + 0.5 # these values can be altered
      } 
      else if (study_type == "Motif") {
        matrix[matrix[, col] != 0, col] <- matrix[matrix[, col] != 0, col] + 1 # these values can be altered
      }
    }
  }
  return(matrix)
}

modified_heatmap_data <- modify_matrix_values(filtered_matrix_genelist_genes, PPRE_metadata, "Col_name", "Type_of_data")

### Cummulative score per gene
cum_matrix <- rowSums(modified_heatmap_data)
row_sums <- data.frame(`Cumulative Score` = cum_matrix)
rownames(row_sums) <- rownames(modified_heatmap_data)

## setting gene names as rownames
genes_ensembl_id <- rownames(modified_heatmap_data)

genenames_for_matrix <- database_PPRE_FAO %>% 
  select(ensembl_gene_id, gene_name) %>%
  filter(ensembl_gene_id %in% genes_ensembl_id) %>%
  select(gene_name)

genenames_for_matrix <- as.list(genenames_for_matrix)

rownames(modified_heatmap_data) <- genenames_for_matrix[["gene_name"]]


#### Visualize the overlap between the given gene list and other studies ####
#my_color_ramp <- colorRampPalette(c("white", "lightblue", "darkblue"))(n=10)

# annotation
study_type <- PPRE_metadata[8:29,] %>% select(Col_name, Type_of_data)

target <- PPRE_metadata[8:29,] %>% select(Col_name, Target)

# Create the annotation DataFrame for columns
annotation_columns <- data.frame(
  Study_type = study_type,
  Target = target
)

# Set first column to rownames
rownames(annotation_columns) <- annotation_columns$Study_type.Col_name

# Select and rename the columns required for annotation (data type and target)
annotation_columns <- annotation_columns %>% 
  select(Study_type.Type_of_data, Target.Target) %>% 
  rename(Type_of_data ="Study_type.Type_of_data", Target = "Target.Target")

ann_colors <- list(
  Type_of_data = c(
    "ChIP-seq" = "#377EB8",     # Blue
    "Motif" = "#E41A1C",        # Red
    "RNA-seq" = "#4DAF4A",      # Green
    "ATAC-seq" = "#984EA3",     # Purple
    "GO" = "#FF7F00",           # Orange
    "All" = "#FFFF33"           # Yellow
  ), 
  Target = c(
    "PPARa" = "#3CAEA3",        # Teal
    "RXR" = "#FF5733",          # Coral
    "PPARa-RXRa" = "#D7E8BA",   # Pale Green
    "PPARy-RXRa" = "#BDD7EE",   # Sky Blue
    "PPARg" = "#FFC700",        # Gold
    "RNA" = "#1B9E77",          # Dark Green
    "Open chromatin" = "#666666",  # Gray
    "H3K27ac" = "#FFA500",      # Orange
    "Proteins" = "#A6CEE3",     # Light Blue
    "NA" = "grey"               # Gray
  )
)



# Create the heatmap
pheatmap(
  modified_heatmap_data,
  fontsize_row = 10,
  fontsize_col = 12,
  legend = TRUE,
  # annotation_col = annotation_columns,
  # annotation_colors = ann_colors,
  # annotation_row = row_sums,
  cluster_cols = TRUE,  # Don't cluster rows
  cluster_rows = TRUE,  # Cluster columns
  show_rownames = TRUE,  # Show row names (Ensembl gene IDs)
  show_colnames = TRUE,  # Show column names (category names)
  color = my_color_ramp  # Color scheme (white for absence, blue for presence)
)
