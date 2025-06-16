### module load HGI/softpack/groups/alpha-allostery-global/xl7_proteome_lmm-2/1
### rstudio start -q hugemem -M 1000000 -c 6

# Load necessary libraries
library(lme4)
library(data.table)
library(future.apply)
library(performance)

# Set up parallelization
future::plan(multisession)

# Paths to data files
ddg_path <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/thermompnn_outs"
esm_path <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_inputs/unzipped/all_ESM1v"

# Output directory for results
output_dir <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/lm_rds"
#dir.create(output_dir, showWarnings = FALSE)

# Function to process each protein
process_protein <- function(protein_id, ddg_path, esm_path, output_dir) {
  # File paths
  ddg_file <- file.path(ddg_path, paste0("ThermoMPNN_inference_AF-", protein_id, "-F1-model_v4.csv"))
  esm_file <- file.path(esm_path, paste0(protein_id, ".csv"))
  
  # Skip if files are missing
  if (!file.exists(ddg_file) || !file.exists(esm_file)) {
    return(NULL)
  }
  
  # Read data
  test_ddg <- fread(ddg_file)
  test_esm <- fread(esm_file)
  colnames(test_esm)[2] <- "ESM1v"  # Rename second column to ESM1v
  
  # Transform and merge
  test_ddg[, new_position := position + 1]
  test_ddg[, variant := paste0(wildtype, new_position, mutation)]
  test_df <- merge(test_ddg, test_esm, by = "variant")
  
  # Skip if no data after merging
  if (nrow(test_df) == 0) {
    return(NULL)
  }
  
  ### Linear Mixed Model (LMM)
  lmm_fit <- tryCatch({
    lmer(ESM1v ~ ddG_pred + (1 | new_position), data = test_df)
  }, error = function(e) {
    message("LMM failed for protein: ", protein_id, " - ", e$message)
    NULL
  })
  
  if (!is.null(lmm_fit)) {
    # Add predictions to the data frame
    test_df[, predicted_ESM1v_lmm := predict(lmm_fit, newdata = test_df, re.form = NULL)]
    
    # Save LMM results
    lmm_results <- list(
      fit = lmm_fit,
      marginal_r_squared = r2(lmm_fit)$R2_marginal,
      conditional_r_squared = r2(lmm_fit)$R2_conditional,
      predicted_data = test_df[, .(variant, new_position, ddG_pred, ESM1v, predicted_ESM1v_lmm, pdb)]
    )
    saveRDS(lmm_results, file = file.path(output_dir, paste0(protein_id, "_lmm.rds")))
  }
  
  ### Linear Model (LM)
  lm_fit <- tryCatch({
    lm(ESM1v ~ ddG_pred, data = test_df)
  }, error = function(e) {
    message("LM failed for protein: ", protein_id, " - ", e$message)
    NULL
  })
  
  if (!is.null(lm_fit)) {
    # Add predictions to the data frame
    test_df[, predicted_ESM1v_lm := predict(lm_fit, newdata = test_df)]
    
    # Save LM results
    lm_results <- list(
      fit = lm_fit,
      r_squared = summary(lm_fit)$r.squared,
      adjusted_r_squared = summary(lm_fit)$adj.r.squared,
      predicted_data = test_df[, .(variant, new_position,  ddG_pred, ESM1v, predicted_ESM1v_lm, pdb)]
    )
    saveRDS(lm_results, file = file.path(output_dir, paste0(protein_id, "_lm.rds")))
  }
  
  ### Linear Model 2 (LM2): ddG_pred + ddG_pred^2
  test_df[, ddG_pred_sq := ddG_pred^2]
  
  lm2_fit <- tryCatch({
    lm(ESM1v ~ ddG_pred + ddG_pred_sq, data = test_df)
  }, error = function(e) {
    message("LM2 failed for protein: ", protein_id, " - ", e$message)
    NULL
  })
  
  if (!is.null(lm2_fit)) {
    test_df[, predicted_ESM1v_lm2 := predict(lm2_fit, newdata = test_df)]
    
    lm2_results <- list(
      fit = lm2_fit,
      r_squared = summary(lm2_fit)$r.squared,
      adjusted_r_squared = summary(lm2_fit)$adj.r.squared,
      predicted_data = test_df[, .(variant, new_position, ddG_pred, ddG_pred_sq, ESM1v, predicted_ESM1v_lm2, pdb)]
    )
    saveRDS(lm2_results, file = file.path(output_dir, paste0(protein_id, "_lm2.rds")))
  }
  
  # Return a summary of the results
  return(list(
    protein_id = protein_id, 
    lmm_success = !is.null(lmm_fit), 
    lm_success = !is.null(lm_fit),
    lm2_success = !is.null(lm2_fit)
  ))
}

# # List of protein IDs
# protein_ids <- c("Q96HR3") 
# # Example protein IDs; replace with the actual list of 19,000 IDs
# 
# # Run processing in parallel
# results <- future_lapply(protein_ids, process_protein, ddg_path = ddg_path, esm_path = esm_path, output_dir = output_dir)
# results
# 
# # Load the .rds file
# lmm_rds <- readRDS("/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/lm_rds/Q96HR3_lmm.rds")  # Replace with the actual file path
# 
# # Print the object to inspect contents
# print(lmm_rds)
# 
# 
# # Load the .rds file
# lm_rds <- readRDS("/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/lm_rds/Q96HR3_lm.rds")  # Replace with the actual file path
# 
# # Print the object to inspect contents
# print(lm_rds)
# 
# # Load the .rds file
# lm2_rds <- readRDS("/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/lm_rds/Q96HR3_lm2.rds")  # Replace with the actual file path
# 
# # Print the object to inspect contents
# print(lm2_rds)

###########################
# List all CSV files in the directory
file_list <- list.files(path = ddg_path, pattern = "\\.csv$", full.names = FALSE)

# Extract the PDB IDs using regex
pdb_ids <- sub(".*AF-([A-Z0-9]+)-F1-model_v4\\.csv", "\\1", file_list)

# Save the PDB IDs to a text file
output_file <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/pdb_ids_batch1.txt"
writeLines(pdb_ids, output_file)

# Run processing in parallel
results_all <- future_lapply(pdb_ids, process_protein, ddg_path = ddg_path, esm_path = esm_path, output_dir = output_dir)

# Summarize results
results_dt <- rbindlist(results_all, fill = TRUE)
print(results_dt)

################################################################################
# List all CSV files in the directory
file_list <- list.files(path = ddg_path, pattern = "\\.csv$", full.names = FALSE)
length(file_list)

# Extract the PDB IDs using regex
pdb_ids <- sub(".*AF-([A-Z0-9]+)-F1-model_v4\\.csv", "\\1", file_list)
length(pdb_ids) #20144

# Define the directory containing the files
out_dir_path <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/lm_rds"

# List all files ending with '_lmm.rds' in the directory
out_file_list <- list.files(path = out_dir_path, pattern = "_lm2.rds$", full.names = FALSE)

# Extract the PDB IDs using regex
out_ids <- sub("_lm2.rds$", "", out_file_list)

length(out_ids) #17259

remaining_id <- setdiff(pdb_ids, out_ids)
length(remaining_id) #2575

process_protein <- function(protein_id, ddg_path, esm_path, output_dir) {
  # File paths
  ddg_file <- file.path(ddg_path, paste0("ThermoMPNN_inference_AF-", protein_id, "-F1-model_v4.csv"))
  esm_file <- file.path(esm_path, paste0(protein_id, ".csv"))
  
  # Skip if files are missing
  if (!file.exists(ddg_file) || !file.exists(esm_file)) {
    return(NULL)
  }
  
  # Read data
  test_ddg <- fread(ddg_file)
  test_esm <- fread(esm_file)
  colnames(test_esm)[2] <- "ESM1v"  # Rename second column to ESM1v
  
  # Transform and merge
  test_ddg[, new_position := position + 1]
  test_ddg[, variant := paste0(wildtype, new_position, mutation)]
  test_df <- merge(test_ddg, test_esm, by = "variant")
  
  # Skip if no data after merging
  if (nrow(test_df) == 0) {
    return(NULL)
  }
  
  # Ensure new_position is a factor and remove rows with NA
  test_df <- test_df[!is.na(new_position)]
  test_df[, new_position := as.factor(new_position)]
  
  ### Linear Mixed Model (LMM)
  lmm_fit <- tryCatch({
    lmer(ESM1v ~ ddG_pred + (1 | new_position), data = test_df)
  }, error = function(e) {
    message("LMM failed for protein: ", protein_id, " - ", e$message)
    NULL
  })
  
  if (!is.null(lmm_fit)) {
    # Align levels for predictions
    test_df[, new_position := factor(new_position, levels = levels(lmm_fit@frame$new_position))]
    
    # Add predictions to the data frame
    test_df[, predicted_ESM1v_lmm := predict(lmm_fit, newdata = test_df, allow.new.levels = TRUE)]
    
    # Save LMM results
    lmm_results <- list(
      fit = lmm_fit,
      marginal_r_squared = r2(lmm_fit)$R2_marginal,
      conditional_r_squared = r2(lmm_fit)$R2_conditional,
      predicted_data = test_df[, .(variant, new_position, ddG_pred, ESM1v, predicted_ESM1v_lmm, pdb)]
    )
    saveRDS(lmm_results, file = file.path(output_dir, paste0(protein_id, "_lmm.rds")))
  }
  
  ### Linear Model (LM)
  lm_fit <- tryCatch({
    lm(ESM1v ~ ddG_pred, data = test_df)
  }, error = function(e) {
    message("LM failed for protein: ", protein_id, " - ", e$message)
    NULL
  })
  
  if (!is.null(lm_fit)) {
    # Add predictions to the data frame
    test_df[, predicted_ESM1v_lm := predict(lm_fit, newdata = test_df)]
    
    # Save LM results
    lm_results <- list(
      fit = lm_fit,
      r_squared = summary(lm_fit)$r.squared,
      adjusted_r_squared = summary(lm_fit)$adj.r.squared,
      predicted_data = test_df[, .(variant, new_position, ddG_pred, ESM1v, predicted_ESM1v_lm, pdb)]
    )
    saveRDS(lm_results, file = file.path(output_dir, paste0(protein_id, "_lm.rds")))
  }
  
  ### Linear Model 2 (LM2): ddG_pred + ddG_pred^2
  test_df[, ddG_pred_sq := ddG_pred^2]
  
  lm2_fit <- tryCatch({
    lm(ESM1v ~ ddG_pred + ddG_pred_sq, data = test_df)
  }, error = function(e) {
    message("LM2 failed for protein: ", protein_id, " - ", e$message)
    NULL
  })
  
  if (!is.null(lm2_fit)) {
    test_df[, predicted_ESM1v_lm2 := predict(lm2_fit, newdata = test_df)]
    
    lm2_results <- list(
      fit = lm2_fit,
      r_squared = summary(lm2_fit)$r.squared,
      adjusted_r_squared = summary(lm2_fit)$adj.r.squared,
      predicted_data = test_df[, .(variant, new_position, ddG_pred, ddG_pred_sq, ESM1v, predicted_ESM1v_lm2, pdb)]
    )
    saveRDS(lm2_results, file = file.path(output_dir, paste0(protein_id, "_lm2.rds")))
  }
  
  # Return a summary of the results
  return(list(
    protein_id = protein_id, 
    lmm_success = !is.null(lmm_fit), 
    lm_success = !is.null(lm_fit),
    lm2_success = !is.null(lm2_fit)
  ))
}

# Run processing in parallel
results_remaining <- future_lapply(remaining_id, process_protein, ddg_path = ddg_path, esm_path = esm_path, output_dir = output_dir)
















