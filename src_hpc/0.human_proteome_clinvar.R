
# Set up parallel processing
plan(multisession, workers = availableCores())  # Or set a specific number

csv_dir <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/meta"
csv_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)
csv_files

total_files <- length(csv_files) #20144
fread(csv_files[1])

# Define the clinical significance values to filter
clinical_significance_filter <- c("benign", "likely_benign", "pathogenic", "likely_pathogenic")

# Function to read a single file, filter rows, and keep all columns
read_file <- function(file, index) {
  dt <- fread(file)  # Read the entire file
  filtered_dt <- dt[clinvar_clinical_significance %in% clinical_significance_filter]  # Filter rows but keep all columns
  
  # Manual progress bar
  cat(sprintf("\rProcessing file %d of %d (%.2f%%)", index, total_files, (index / total_files) * 100))
  flush.console()  # Ensure the progress updates in the console
  
  return(filtered_dt)
}

# Read and combine all CSV files in parallel
# this takes a while, be patient 
combined_df_new <- future_lapply(seq_along(csv_files), function(i) read_file(csv_files[i], i), future.seed = TRUE) |> rbindlist(fill = TRUE)

# Print completion message
cat("\nAll files processed successfully!\n")

head(combined_df_new)
nrow(combined_df_new) #76853

write.csv(combined_df, "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/clinvar_all_protein_filtered_new.csv", row.names = FALSE)

length(unique(combined_df_new$uniprot))

# Paths
esm1v_dir <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_inputs/unzipped/all_ESM1v"

# List all ESM1v CSV files and create a mapping of UniProt IDs
esm1v_files <- list.files(esm1v_dir, pattern = "*.csv", full.names = TRUE)
esm1v_map <- setNames(esm1v_files, str_replace(basename(esm1v_files), ".csv$", ""))

# Function to merge ESM1v data into combined_df
merge_esm1v_to_combined_df <- function(combined_df, esm1v_map) {
  
  # Split combined_df by uniprot for individual processing
  combined_split <- split(combined_df, combined_df$uniprot)
  
  # Initialize list to store merged data
  merged_data_list <- list()
  
  # Track all columns that appear during merging
  all_columns <- colnames(combined_df)
  
  # Process each protein's data
  for (uniprot_id in names(combined_split)) {
    
    protein_data <- combined_split[[uniprot_id]]
    
    # Check if the uniprot_id has a corresponding ESM1v file
    if (uniprot_id %in% names(esm1v_map)) {
      esm1v_file <- esm1v_map[[uniprot_id]]
      df_esm1v <- read.csv(esm1v_file, stringsAsFactors = FALSE)
      
      # Check if the ESM1v file contains the "variant" column
      if ("variant" %in% colnames(df_esm1v)) {
        
        # Merge ESM1v scores by "variant"
        protein_data <- protein_data %>%
          left_join(df_esm1v, by = "variant")  # Merging by "variant"
        
        # Update the list of all columns encountered
        all_columns <- unique(c(all_columns, colnames(protein_data)))
        
        message("Merged ESM1v data for: ", uniprot_id)
        
      } else {
        message("No 'variant' column found in ESM1v file for: ", uniprot_id)
      }
    } else {
      message("No ESM1v file found for: ", uniprot_id)
    }
    
    # Save the processed data for each uniprot_id
    merged_data_list[[uniprot_id]] <- protein_data
  }
  
  # Ensure all data frames have the same columns by filling missing columns with NA
  for (name in names(merged_data_list)) {
    missing_cols <- setdiff(all_columns, colnames(merged_data_list[[name]]))
    if (length(missing_cols) > 0) {
      merged_data_list[[name]][, missing_cols] <- NA
    }
  }
  
  # Combine all processed data back into a single DataFrame
  combined_merged_df <- do.call(rbind, merged_data_list)
  
  return(combined_merged_df)
}


# Merge ESM1v scores into combined_df
combined_merged_df_new <- merge_esm1v_to_combined_df(combined_df_new, esm1v_map)
nrow(combined_merged_df_new) #76853
length(unique(combined_merged_df_new$uniprot)) #10899
head(combined_merged_df_new)

# Save the final DataFrame to file
write.csv(combined_merged_df_new, file.path("/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/clinvar_all_protein_filtered_with_esm1v_new.csv"), row.names = FALSE)

























