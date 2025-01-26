library(vcfR)

vcf <- read.vcfR("/Users/swethayadavalli/Downloads/overlaped_file.vcf")#reading in the overlapped vcf file only with the block regions from 1000 genome blocks

blocks <- read.table('/Users/swethayadavalli/Downloads/LD_blocks_ADNI_19_2000kb.blocks.det', header = TRUE) # reading in the blocks file

# Creating a new directory for the .rds files
dir.create("rds_files")
dir.create("counts")
# Looping over all the rows in the blocks dataframe
for (i in 1:nrow(blocks)) {
  # Extracting the i-th block
  block <- blocks[i, ]
  
  # Splitting the SNPs string into a vector of SNP IDs
  snps <- strsplit(block$SNPS, "\\|")[[1]]
  
  # Extracting the SNPs from the VCF file
  block_snps <- vcf[vcf@fix[, "ID"] %in% snps, ]
  
  samples <- data.frame(block_snps@gt)
  rsIDs <- data.frame(block_snps@fix)
  
  # Selecting the "ID", "REF", and "ALT" columns from the rsids data frame
  selected_columns <- rsIDs[, c("ID", "REF", "ALT")]
  
  combined_df <- cbind(samples, selected_columns)
  
  # Getting the number of columns in the data frame
  num_columns <- ncol(combined_df)
  
  # Creating a vector of column indices in the desired order
  new_order <- c((num_columns-2):num_columns, 1:(num_columns-3))
  
  # Rearranging the columns
  rearranged_df <- combined_df[, new_order]
  
  # Removing the fourth column
  Block_df <- rearranged_df[,-4]
  
  # Saving the dataframe as an .rds file in the new directory
  saveRDS(Block_df, file = paste0("rds_files/Block", i, ".rds"))
}


# code for unique sequence (loop)
library(dplyr)
library(tidyr)
library(reshape2)

process_block <- function(Block_df) {
  # Subsetting the first three columns of Block_df
  ref_alt_rsid <- Block_df %>% select(1:3)
  
  # Removing the first three columns of Block_df
  all_blocks <- Block_df %>% select(-1:-3)
  
  # Getting the column names
  col_names <- colnames(all_blocks)
  
  # Initializing the dataframe
  separated_all_blocks <- all_blocks
  
  # Looping over each column of all_blocks
  for (i in seq_along(col_names)) {
    # Separating the "0|0" into two separate columns titled "REF_sample" and "ALT_sample"
    separated_all_blocks <- separated_all_blocks %>%
      separate(col_names[i], into = c(paste0("MAT_", col_names[i]), paste0("PAT_", col_names[i])), sep = "\\|")
  }
  
  # Combining the dataframes
  separated_block_df <- cbind(ref_alt_rsid, separated_all_blocks)
  
  # Getting the column names
  col_names <- colnames(separated_block_df)[4:ncol(separated_block_df)]  # Adjust this to match the columns you want to change
  # Replacing 0s and 1s with REF and ALT 
  separated_block_df2 <- separated_block_df %>%
    mutate_at(vars(col_names), ~ifelse(. == 0, separated_block_df$REF, ifelse(. == 1, separated_block_df$ALT, .)))
  
  
  # Creating a dataframe with only the 'MAT' Maternal columns
  df_MAT <- separated_block_df2 %>% select(ID, REF, ALT, starts_with("MAT_"))
  
  # Creating a dataframe with only the 'PAT' Paternal columns
  df_PAT <- separated_block_df2 %>% select(ID, REF, ALT, starts_with("PAT_"))
  
  # Removing the 'REF' and 'ALT' columns
  df_MAT2 <- df_MAT %>% select(-REF, -ALT)
  
  # Transposing the dataframe
  df_MAT_transposed <- as.data.frame(t(df_MAT2))
  
  # Setting the first row as column names
  colnames(df_MAT_transposed) <- df_MAT_transposed[1,]
  
  # Removing the first row
  df_MAT_transposed <- df_MAT_transposed[-1,]
  
  # Combining all the values from different columns into a single string for each row
  df_MAT_transposed$Combined <- apply(df_MAT_transposed, 1, paste, collapse = "")
  
  # Creating a new dataframe, preserving the row names
  df_MAT_sequences <- data.frame(RowNames = row.names(df_MAT_transposed), Sequence1_block1 = df_MAT_transposed$Combined)
  
  # Setting the row names of the dataframe to the 'RowNames' column
  row.names(df_MAT_sequences) <- df_MAT_sequences$RowNames
  
  # Removing the 'RowNames' column
  df_MAT_sequences$RowNames <- NULL
  
  # Removing the 'REF' and 'ALT' columns
  df_PAT2 <- df_PAT %>% select(-REF, -ALT)
  
  # Transposing the dataframe
  df_PAT_transposed <- as.data.frame(t(df_PAT2))
  
  # Setting the first row as column names
  colnames(df_PAT_transposed) <- df_PAT_transposed[1,]
  
  # Removing the first row
  df_PAT_transposed <- df_PAT_transposed[-1,]
  
  # Combining all the values from different columns into a single string for each row
  df_PAT_transposed$Combined <- apply(df_PAT_transposed, 1, paste, collapse = "")
  
  # Creating a new dataframe, preserving the row names
  df_PAT_sequences <- data.frame(RowNames = row.names(df_PAT_transposed), Sequence1_block1 = df_PAT_transposed$Combined)
  
  # Setting the row names of the dataframe to the 'RowNames' column
  row.names(df_PAT_sequences) <- df_PAT_sequences$RowNames
  
  # Removing the 'RowNames' column
  df_PAT_sequences$RowNames <- NULL
  
  combined_MAT_PAT_df <- cbind(df_MAT_sequences, df_PAT_sequences)
  colnames(combined_MAT_PAT_df)[1] <- "MAT"
  colnames(combined_MAT_PAT_df)[2] <- "PAT"
  
  # Converting row names to a column
  combined_MAT_PAT_df$ID <- rownames(combined_MAT_PAT_df)
  
  # Reshaping the dataframe from wide to long format
  df_long <- melt(combined_MAT_PAT_df, id.vars = 'ID')
  
  # Creating a function to append a suffix to duplicate row names
  make_unique <- function(x){
    return(paste0(x, ".", ave(seq_along(x), x, FUN = seq_along)))
  }
  
  # Applying the function to the 'ID' column
  df_long$ID <- make_unique(df_long$ID)
  
  # Setting the 'ID' column as the row names of the dataframe
  rownames(df_long) <- df_long$ID
  
  # Removing the 'ID' column
  df_long$ID <- NULL
  df_long$variable <- NULL
  
  # Now proceeding with creating the counts matrix
  unique_sequences <- unique(df_long$value)
  
  # Creating a counts matrix
  counts_matrix <- table(rownames(df_long), df_long$value)
  
  # Converting the counts matrix to a dataframe
  counts_df <- as.data.frame.matrix(counts_matrix)
  
  # Checking the sum of all values in each row before removing the rows with all zero counts
  row_sums_before <- rowSums(counts_df)
  
  # Removeing rows with all zero counts
  counts_df <- counts_df[apply(counts_df, 1, function(x) any(x != 0)), ]
  
  # Converting the counts matrix to a dataframe and reset the row names
  counts_df$ID <- rownames(counts_df)
  rownames(counts_df) <- NULL
  
  # Splitting the ID and keep only the part before the suffix
  counts_df$ID <- sub("\\..*", "", counts_df$ID)
  
  
  # Grouping by the ID and sum the counts for each group
  counts_df <- counts_df %>%
    group_by(ID) %>%
    summarise_all(sum)
  
  
  # Converting the tibble to a data frame
  counts_df <- as.data.frame(counts_df)
  
  # Setting the ID column as row names
  rownames(counts_df) <- counts_df$ID
  
  # Removing the ID column
  counts_df$ID <- NULL
  
  # Converting back to a matrix
  counts_matrix <- as.matrix(counts_df)
}
# Defining the directory where the .rds files are stored
rds_directory <- "/Users/swethayadavalli/rds_files"

# Defining the directory where the output .rds files should be saved
counts_directory <- "/Users/swethayadavalli/counts"

# Defining the sequence of blocks
blocks <- 1:4066  # replace with your actual sequence

# Initializing an empty list to store the results
results <- list()

# Looping over each block
for (i in blocks) {
  # Loading the block from an .rds file
  Block_df <- readRDS(file.path(rds_directory, paste0("Block", i, ".rds")))
  
  # Processing the block and store the result
  results[[i]] <- process_block(Block_df)
  
  # Saving the result to an .rds file in the counts directory
  saveRDS(results[[i]], file.path(counts_directory, paste0("counts_df_Block", i, ".rds")))
}

######### Code to append the Block number to haplotypes and generating a combined matrix ##########

# Setting the directory where your .rds files are stored
directory <- "/Users/swethayadavalli/counts"

# Creating a directory to store the modified result
modified_directory <- "/Users/swethayadavalli/combined_matrices"
dir.create(modified_directory, showWarnings = FALSE)

# Getting the list of .rds files in the directory
files <- Sys.glob(file.path(directory, "counts_df_Block*.rds"))

# Extracting the block numbers from the filenames, convert them to integers, and sort the filenames based on these integers
files <- files[order(as.integer(gsub(".*counts_df_Block(\\d+)\\.rds$", "\\1", basename(files))))]

# Initializing a list to store modified matrices
modified_matrices <- list()

# Looping through each file
for (i in 1:length(files)) {
  # Reading the count matrix from the .rds file
  count_matrix <- readRDS(files[i])
  
  # Extracting the block number from the filename
  block_num <- stringr::str_extract(basename(files[i]), "\\d+")
  
  # Generating the block label
  block_label <- paste0("B", block_num)
  
  # Appending the block label to each column name
  colnames(count_matrix) <- paste0(block_label, "_", colnames(count_matrix))
  
  # Storing the modified matrix
  modified_matrices[[i]] <- count_matrix
  
  # Saving the modified matrix as .rds file
  saveRDS(count_matrix, file.path(modified_directory, paste0("modified_count_matrix_", block_num, ".rds")))
}

# Combining all the modified matrices into a single matrix
combined_matrix <- do.call(cbind, modified_matrices)

# Saving the combined matrix as .rds file
saveRDS(combined_matrix, file.path(modified_directory, "combined_matrix.rds"))


#################     Haplotype Association analysis ###############

# Loading the .fam file
fam_data <- read.table("/Users/swethayadavalli/Downloads/ADNI_HRC_1GO23_V3.fam")

# Set the column names
colnames(fam_data) <- c("FamilyID", "IndividualID", "PaternalID", "MaternalID", "Sex", "Phenotype")

# Converting the phenotype data to a factor
phenotype_factor <- factor(fam_data$Phenotype, levels = c(1, 2), labels = c("control", "AD"))
names(phenotype_factor) <- fam_data$IndividualID

# Continuing with the Chi-square test as before...
# Replace 'MAT_X2_' in the row names of the combined matrix
rownames(combined_matrix) <- sub("MAT_X2_", "", rownames(combined_matrix))

# Making sure the names of the phenotype factor match the row names of the combined matrix
phenotype_factor <- phenotype_factor[rownames(combined_matrix)]

# Initializing a list to store the results
chi_square_results <- list()

# Opening the file for writing
sink("chi_square_results.txt")

# Perform the Chi-square test for each column of combined_matrix
for (i in seq_len(ncol(combined_matrix))) {
  contingency_table <- table(phenotype_factor, combined_matrix[, i])
  chi_square_result <- chisq.test(contingency_table)
  
  # Printing the results
  cat("Column:", i, "\n")
  print(paste("Chi-square statistic:", chi_square_result$statistic))
  print(paste("Degrees of freedom:", chi_square_result$parameter))
  print(paste("p-value:", chi_square_result$p.value))
  print(paste("Method used:", chi_square_result$method))
  cat("\n")
}

# Closing the file
sink()

####### Converting the results if the chi-square to cvs ###########
# Initializing a data frame to store the results
chi_square_results <- data.frame(Column = integer(), Statistic = numeric(), DF = integer(), PValue = numeric(), Method = character())

# Performing the Chi-square test for each column of combined_matrix
for (i in seq_len(ncol(combined_matrix))) {
  contingency_table <- table(phenotype_factor, combined_matrix[, i])
  chi_square_result <- chisq.test(contingency_table)
  
  # Storing the results in the data frame
  chi_square_results <- rbind(chi_square_results, data.frame(Column = i, Statistic = chi_square_result$statistic, DF = chi_square_result$parameter, PValue = chi_square_result$p.value, Method = chi_square_result$method))
}

# Printing the results
print(chi_square_results)


# Saving the data frame to a CSV file
write.csv(chi_square_results, file = "/Users/swethayadavalli/Downloads/chi_square_results.csv", row.names = FALSE)

#reading in th efile that we got above here
df <- read.csv("/Users/swethayadavalli/Downloads/chi_square_results.csv", header = TRUE)

# Loading the necessary libraries
library(stringr)

# Reading the .csv file
df <- read.csv("/Users/swethayadavalli/Downloads/chi_square_results.csv", header = TRUE)

# Extracting the column numbers and p-values
columns <- df$Column
p_values <- df$PValue

# Creating a data frame with the column numbers and p-values
df <- data.frame(Column = columns, PValue = p_values)

# Removing rows with NA values
df <- df[complete.cases(df), ]

# Filtering the rows where the p-value is less than 0.05
filtered_df <- df[df$PValue < 0.05, ]

# Printing the filtered data frame
print(filtered_df)

# Extracting the column names using the column numbers in filtered_df
filtered_df$ColumnName <- colnames(combined_matrix)[filtered_df$Column]


# Adding a new column with the -log10 of the p-values
filtered_df$NegLogPValue <- -log10(filtered_df$PValue)

# Creating a subset of the data frame with the most significant p-values
top_df <- filtered_df[filtered_df$NegLogPValue > -log10(0.01), ]  # Adjust the threshold as needed

library(ggplot2)
# Creating the Manhattan plot
ggplot(filtered_df, aes(x = Column, y = NegLogPValue)) +
  geom_point() +
  geom_text(data = top_df, aes(label = Column), vjust = -0.5, size = 3) +
  theme_minimal() +
  labs(x = "Column Number", y = "-log10(P Value)") +
  ggtitle("Manhattan Plot of P Values by Column Number")






