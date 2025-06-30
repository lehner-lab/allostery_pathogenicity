library(data.table)
library(future.apply)
library(performance)
library(readr)
library(utils)
library(furrr)
library(arrow)

# Set up parallelization
future::plan(multisession)

process_proteins_parallel <- function(ddg_dir, var_dir, output_dir) {
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Get CSV files
  ddg_files <- list.files(ddg_dir, full.names = TRUE, pattern = "\\.csv$")
  var_files <- list.files(var_dir, full.names = TRUE, pattern = "\\.csv$")
  
  # Extract UniProt IDs
  ddg_ids <- sub(".*ThermoMPNN_inference_AF-([A-Z0-9]+)-.*", "\\1", basename(ddg_files))
  var_ids <- sub("^([A-Z0-9]+)-.*", "\\1", basename(var_files))
  
  # Match only shared UniProt IDs
  matched_ids <- intersect(ddg_ids, var_ids)
  ddg_files <- ddg_files[ddg_ids %in% matched_ids]
  var_files <- var_files[var_ids %in% matched_ids]
  
  total_tasks <- length(ddg_files)
  progress <- 0
  
  # Single pair processor
  process_single_pair <- function(ddg_file, var_file) {
    ddg_df <- fread(ddg_file)
    var_df <- fread(var_file)
    
    # Construct variant column in ddg_df
    setDT(ddg_df)
    ddg_df[, new_position := position + 1]
    ddg_df[, variant := paste0(wildtype, new_position, mutation)]
    
    # Filter relevant columns from var_df (use only existing ones safely)
    keep_cols <- c("uniprot", "variant", "exposure_ss", "exposure_rasa", "spot_disorder",
                   "clinvar_clinical_significance", "aa_sequence")
    
    keep_cols <- keep_cols[keep_cols %in% names(var_df)]
    var_df_fil <- var_df[, ..keep_cols]
    
    # Merge on variant
    merged_df <- merge(ddg_df, var_df_fil, by = "variant", all = FALSE)
    return(merged_df)
  }
  
  update_progress <- function(current, total) {
    cat(sprintf("\rProgress: %d/%d (%.2f%%)", current, total, (current / total) * 100))
    flush.console()
  }
  
  # Parallel execution
  furrr::future_map2(ddg_files, var_files, ~{
    protein_id <- sub(".*ThermoMPNN_inference_AF-([A-Z0-9]+)-.*", "\\1", basename(.x))
    result <- process_single_pair(.x, .y)
    
    output_file <- file.path(output_dir, paste0(protein_id, "_clinvar.csv"))
    fwrite(result, output_file)
    
    progress <<- progress + 1
    update_progress(progress, total_tasks)
  })
  
  cat("\nProcessing complete! All results saved.\n")
}


# Define directories
ddg_dir <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/thermompnn_outs"
var_dir <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_inputs/cagiada_csv"
output_dir <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/meta"

# Run the processing function
process_proteins_parallel(ddg_dir, var_dir, output_dir)

###############################################################################################################

bi_df <- read.table("/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_inputs/pioneer_interactome/human_high_very_high.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)  
nrow(bi_df) #146138
bi_df <- bi_df %>% dplyr::select (uniprot1, uniprot2, interface_residues1, interface_residues2)
head(bi_df)

# Step 1: Separate the DataFrame into two DataFrames
df1 <- bi_df %>%
  dplyr::select(uniprot = uniprot1, interface_residues = interface_residues1)

df2 <- bi_df %>%
  dplyr::select(uniprot = uniprot2, interface_residues = interface_residues2)

# Step 2: Combine the two DataFrames
bi_df_merged <- bind_rows(df1, df2)

# Step 3: Group by UniProt ID and aggregate interface residues
flattened_df <- bi_df_merged %>%
  group_by(uniprot) %>%
  summarise(interface_residues = list(unique(unlist(interface_residues))))

flattened_df_cleaned <- flattened_df %>%
  mutate(
    interface_residues = str_replace_all(interface_residues, "c\\(|\\)|\\[|\\]|\"", ""), # Remove c(), [], and quotes
    interface_residues = str_split(interface_residues, ",\\s*") # Split into individual residues
  ) %>%
  unnest(interface_residues) %>% # Unnest the list of residues
  filter(interface_residues != "") %>% # Remove empty strings
  mutate(interface_residues = as.numeric(interface_residues)) # Convert to numeric

final_df <- flattened_df_cleaned %>%
  group_by(uniprot) %>%
  summarise(interface_residues = list(unique(interface_residues)))

head(final_df)
nrow(final_df) #12569

final_df <- final_df %>%
  mutate(interface_residues = sapply(interface_residues, function(x) paste(x, collapse = ",")))

write.csv(final_df, "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_inputs/pioneer_interactome/all_bi_residues.csv", row.names = FALSE)

###############################################################################################################
# Define directories
input_dir <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/meta"
output_dir <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/meta_core"

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load interface residue data
bi_df_path <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_inputs/pioneer_interactome/all_bi_residues.csv"
bi_df <- read.csv(bi_df_path, stringsAsFactors = FALSE)

# Convert interface_residues to a list of numeric values
bi_df$interface_residues <- lapply(strsplit(bi_df$interface_residues, ","), as.numeric)
head(bi_df)
nrow(bi_df) #12569
sum(is.na(bi_df))

# Lisbi_df# List all CSV files in the directory
file_list <- list.files(input_dir, pattern = "*.csv", full.names = TRUE)

# Function to filter a dataframe based on interface residues
filter_interface_residues <- function(df, uniprot_id) {
  match_row <- bi_df %>% filter(uniprot == uniprot_id)
  if (nrow(match_row) > 0) {
    interface_positions <- unlist(match_row$interface_residues)
    df <- df %>% filter(!(new_position %in% interface_positions))  # Remove matching positions
  }
  return(df)
}

# Process each CSV file
for (file in file_list) {
  # Extract Uniprot ID from the filename
  uniprot_id <- str_replace(basename(file), "_clinvar.csv$", "")
  
  # Read CSV file
  df <- read.csv(file, stringsAsFactors = FALSE)
  
  # Drop 'sequence' column if it exists
  if ("sequence" %in% colnames(df)) {
    df <- df %>% select(-sequence)
  }
  
  # Filter rows based on interface residues
  df_filtered <- filter_interface_residues(df, uniprot_id)
  df_filtered <- df_filtered %>% filter (spot_disorder == 0) %>% filter (exposure_rasa <= 0.2)
  
  # Save filtered file
  output_file <- file.path(output_dir, paste0(uniprot_id, "_filtered.csv"))
  write.csv(df_filtered, output_file, row.names = FALSE)
  
  print(paste("Processed:", uniprot_id, "-> Saved to:", output_file))
}

print("All files processed successfully.")


########################################################################################
# find empty files and fix 
# Set the directory path
directory_path <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/meta_core"

# List all CSV files in the directory
csv_files <- list.files(path = directory_path, pattern = "*.csv", full.names = TRUE)
length(csv_files)

# Fast function to check if a CSV has data rows
is_empty_file_fast <- function(file) {
  tryCatch({
    data <- fread(file, nrows = 2)  # Only read 2 rows to check if data exists
    nrow(data) == 0
  }, error = function(e) {
    TRUE  # If any error, assume empty
  })
}

# Apply in parallel if needed
empty_files <- csv_files[vapply(csv_files, is_empty_file_fast, logical(1L))]
length(empty_files) #2620

# Extract UniProt IDs assuming format: <uniprot_id>_filtered.csv
empty_uniprot_ids <- gsub(".*\\/|_filtered\\.csv", "", empty_files)

# Filter bi_df
matching_uniprots <- bi_df %>% filter(uniprot %in% empty_uniprot_ids)
nrow(matching_uniprots) #1475

###############################
# Function to filter a dataframe based on interface residues
filter_interface_residues <- function(df, uniprot_id) {
  match_row <- bi_df %>% filter(uniprot == uniprot_id)
  
  if (nrow(match_row) > 0) {
    # If matching uniprot found, filter based on interface residues
    interface_positions <- unlist(match_row$interface_residues)
    df <- df %>% filter(!(new_position %in% interface_positions))  # Remove matching positions
  } else {
    # If no match found, filter based on spot_disorder and exposure_rASA
    df <- df %>% filter(spot_disorder == 0, exposure_rasa <= 0.2)
  }
  
  return(df)
}

# Process each CSV file
for (file in empty_files) {
  # Extract Uniprot ID from the filename
  uniprot_id <- str_replace(basename(file), "_filtered.csv", "")
  
  # Read CSV file
  df <- read.csv(file, stringsAsFactors = FALSE)
  
  # Drop 'sequence' column if it exists
  if ("sequence" %in% colnames(df)) {
    df <- df %>% select(-sequence)
  }
  
  # Filter rows using the updated logic
  df_filtered <- filter_interface_residues(df, uniprot_id)
  
  # Save filtered file
  output_file <- file.path(directory_path, paste0(uniprot_id, "_filtered.csv"))
  write.csv(df_filtered, output_file, row.names = FALSE)
  
  print(paste("Processed:", uniprot_id, "-> Saved to:", output_file))
}

print("All files processed successfully.")

tail(empty_files)
test <- read.csv("/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/meta/X6R8D5_clinvar.csv")
test %>% filter(spot_disorder == 0, exposure_rasa <= 0.2)

# empty files don't have residues meeting the condition

# Remove empty files
file.remove(empty_files)

########################################################################################
# Define directories
filtered_dir <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/meta_core/"
esm1v_dir <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_inputs/unzipped/all_ESM1v"
output_dir <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/meta_core_with_esm1v/"

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# List all filtered ClinVar CSV files
filtered_files <- list.files(filtered_dir, pattern = "*.csv", full.names = TRUE)
length(filtered_files) #17524

# List all ESM1v CSV files and create a mapping of UniProt IDs
esm1v_files <- list.files(esm1v_dir, pattern = "*.csv", full.names = TRUE)
esm1v_map <- setNames(esm1v_files, str_replace(basename(esm1v_files), ".csv$", ""))

# Process each filtered ClinVar CSV file
for (file in filtered_files) {
  # Extract UniProt ID from the filename
  uniprot_id <- str_replace(basename(file), "_filtered.csv$", "")
  
  print(paste("Processing:", uniprot_id))
  
  # Read filtered ClinVar CSV file
  df_clinvar <- read.csv(file, stringsAsFactors = FALSE)
  
  # Find the corresponding ESM1v file
  if (uniprot_id %in% names(esm1v_map)) {
    esm1v_file <- esm1v_map[[uniprot_id]]
    df_esm1v <- read.csv(esm1v_file, stringsAsFactors = FALSE)
    
    # Merge ESM1v scores with the filtered ClinVar data
    df_merged <- df_clinvar %>%
      left_join(df_esm1v, by = "variant")  # Ensure correct merging by variant
    
    print(paste("Merged ESM1v data for:", uniprot_id))
  } else {
    df_merged <- df_clinvar
    print(paste("No ESM1v file found for:", uniprot_id))
  }
  
  # Save the merged file
  output_file <- file.path(output_dir, paste0(uniprot_id, "_filtered_with_esm1v.csv"))
  write.csv(df_merged, output_file, row.names = FALSE)
  
  print(paste("Saved:", output_file))
}

print("All files processed successfully.")


########################################################################################
# Set up parallelization
future::plan(multisession)

# Define paths
input_dir <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/meta_core_with_esm1v/"
output_dir <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/meta_core_with_esm1v_out/"
dir.create(output_dir, showWarnings = FALSE)

# Function to process each protein file
process_protein <- function(file_path, output_dir) {
  # Extract protein ID from filename
  protein_id <- gsub("\\_core_with_esm1v.csv$", "", basename(file_path))
  
  # Read data
  protein_data <- tryCatch(fread(file_path), error = function(e) {
    message("Error reading file: ", file_path)
    return(NULL)
  })
  
  if (is.null(protein_data) || nrow(protein_data) == 0) {
    message("Skipping empty or unreadable file: ", file_path)
    return(NULL)
  }
  
  # Ensure required columns exist
  required_columns <- c("variant", "ddG_pred", "ESM.1v")
  if (!all(required_columns %in% names(protein_data))) {
    message("Skipping file due to missing columns: ", file_path)
    return(NULL)
  }
  
  ### Linear Model (LM)
  lm_fit <- tryCatch({
    lm(ESM.1v ~ ddG_pred, data = protein_data)
  }, error = function(e) {
    message("LM failed for protein: ", protein_id, " - ", e$message)
    return(NULL)
  })
  
  if (!is.null(lm_fit)) {
    # Add predictions to the data frame
    protein_data[, predicted_ESM1v_lm := predict(lm_fit, newdata = protein_data)]
    
    # Save LM results
    lm_results <- list(
      fit = lm_fit,
      r_squared = summary(lm_fit)$r.squared,
      adjusted_r_squared = summary(lm_fit)$adj.r.squared,
      predicted_data = protein_data[, .(variant, position, ddG_pred, ESM.1v, predicted_ESM1v_lm)]
    )
    saveRDS(lm_results, file = file.path(output_dir, paste0(protein_id, "_lm.rds")))
  }
  
  # Return a summary of the results
  return(list(protein_id = protein_id, lm_success = !is.null(lm_fit)))
}

# Get the list of CSV files
file_list <- list.files(path = input_dir, pattern = "\\.csv$", full.names = TRUE)
length(file_list) #17524
# Process 
results <- future_lapply(file_list, process_protein, output_dir = output_dir)

# Summarize results
results_dt <- rbindlist(results, fill = TRUE)
results_dt
table(results_dt$lm_success)
# FALSE  TRUE 
# 8 17420
########################################################################################
# Set upggplot2# Set up parallelization
future::plan(multisession)

# Directory containing .rds files
rds_dir <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/meta_core_with_esm1v_out/"

# List all _lmm.rds files

rds_files <- list.files(path = rds_dir, pattern = "_lm.rds$", full.names = TRUE)
length(rds_files) #17420
#test_file <- rds_files[1]
#readRDS(test_file)

extract_r_squared <- function(file) {
  protein_id <- sub("_lm.rds$", "", basename(file))
  
  adj_r2 <- tryCatch({
    lm_rds <- readRDS(file)
    lm_rds$adjusted_r_squared
  }, error = function(e) {
    warning(paste("Skipping", file, ":", e$message))
    NA
  })
  
  data.table(protein_id = protein_id, adj_r2 = adj_r2)
}

all_r2 <- rbindlist(lapply(rds_files, extract_r_squared), fill = TRUE)
nrow(all_r2) #17420
head(all_r2)
tail(all_r2)

sum(!is.na(all_r2$adj_r2))  #17420


library(furrr)
library(purrr)

# Use multiple sessions (parallel workers)
future::plan(multisession)

# Set the directory path where your RDS files are stored
dir_path <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/meta_core_with_esm1v_out/"

# List all RDS files in the directory
rds_files <- list.files(path = dir_path, pattern = "\\.rds$", full.names = TRUE)
length(rds_files) #17420

# Parallelize the row counting
row_counts <- future_map_int(rds_files, function(file) {
  lm_rds <- readRDS(file)
  if (!is.null(lm_rds$predicted_data)) {
    nrow(lm_rds$predicted_data)
  } else {
    0
  }
})

# Sum up the total rows
total_rows <- sum(row_counts)

# Print the result
cat("Total number of rows in predicted_data across all RDS files:", total_rows, "\n")
#Total number of rows in predicted_data across all RDS files: 47841329 

# Calculate average marginal R-squared
median_r_squared <- median(all_r2$adj_r2, na.rm = TRUE)
median_r_squared #0.1581737

range(all_r2$adj_r2) #-0.05881657  0.86482045

write_csv(r_squared_dt_lm, "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/core_lm_r2.csv")

####################################################
####################################################################################
# Directory containing .rds files
rds_dir <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/lm_rds"

# List all _lmm.rds files
# Parallel setup
plan(multisession, workers = availableCores())  # Or set a specific number

# List RDS files
rds_files <- list.files(path = rds_dir, pattern = "_lm.rds$", full.names = TRUE)
length(rds_files) #19852

# Parallel function
extract_r_squared <- function(file) {
  protein_id <- sub("_lm.rds$", "", basename(file))
  adj_r2 <- tryCatch({
    lm_rds <- readRDS(file)
    lm_rds$adjusted_r_squared
  }, error = function(e) {
    NA_real_  # Return numeric NA
  })
  data.table(protein_id = protein_id, adj_r2 = adj_r2)
}

# Run in parallel
all_var_r2 <- rbindlist(future_lapply(rds_files, extract_r_squared), fill = TRUE)

# Result
nrow(all_var_r2)  # 19852
head(all_var_r2)
tail(all_r2)

sum(!is.na(all_r2$adj_r2))  #17420

# Save the results to a CSV 
fwrite(all_var_r2, file.path("/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/all_var_lm_r2.csv"))














