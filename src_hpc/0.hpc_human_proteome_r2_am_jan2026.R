# Load necessary libraries
library(lme4)
library(data.table)
library(future.apply)
library(performance)

# Set up parallelization
future::plan(multisession)

# Input directory containing *_merged.csv
merged_dir <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/00.human_proteome_unzipped/results_unzipped/merged_am_ddgf_jan2026"

# Output directory for LM rds
output_dir <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/00.human_proteome_results/am_lm_rds_jan2026"

process_protein_merged <- function(protein_id, merged_dir, output_dir, epsilon = 1e-10) {
  message("Processing ", protein_id)
  
  merged_file <- file.path(merged_dir, paste0(protein_id, "_merged.csv"))
  if (!file.exists(merged_file)) return(NULL)
  
  # Read merged data
  dt <- tryCatch(fread(merged_file, showProgress = FALSE),
                 error = function(e) {
                   message("  fread failed for ", protein_id, ": ", e$message)
                   return(NULL)
                 })
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  
  # ---- required columns ----
  # UPDATED: Check for AM instead of ESM1v
  req <- c("variant", "ddG_pred", "am_pathogenicity")
  miss <- setdiff(req, names(dt))
  if (length(miss) > 0) {
    message("  Missing columns for ", protein_id, ": ", paste(miss, collapse = ", "))
    return(NULL)
  }
  
  # Ensure numeric
  dt[, ddG_pred := as.numeric(ddG_pred)]
  dt[, am_pathogenicity := as.numeric(am_pathogenicity)]
  
  # Drop rows where key vars are NA
  dt <- dt[!is.na(ddG_pred) & !is.na(am_pathogenicity)]
  if (nrow(dt) == 0) return(NULL)
  
  # ---- Transformation (Specific to AM) ----
  # Clamp AM scores and compute logits (from your original snippet)
  dt[, AM_clamped := pmin(pmax(am_pathogenicity, epsilon), 1 - epsilon)]
  dt[, AM_logit := log(AM_clamped / (1 - AM_clamped))]
  
  # If you want wt_aa and position parsing from "A126C"
  if (!("wt_aa" %in% names(dt))) dt[, wt_aa := substr(variant, 1, 1)]
  if (!("new_position" %in% names(dt))) {
    dt[, new_position := as.integer(sub("^[A-Z](\\d+)[A-Z]$", "\\1", variant))]
  }
  
  # ---- Linear Model (Predicting Logits) ----
  lm_fit <- tryCatch(
    lm(AM_logit ~ ddG_pred, data = dt),
    error = function(e) {
      message("  LM failed for ", protein_id, " - ", e$message)
      NULL
    }
  )
  
  if (is.null(lm_fit)) return(list(protein_id = protein_id, lm_success = FALSE))
  
  # Add predictions (on the logit scale)
  dt[, predicted_AM_lm := predict(lm_fit, newdata = dt)]
  
  # choose columns to save if present
  keep <- c("variant", "new_position", "ddG_pred", 
            "Model", "Dataset", "uniprot", "chain", 
            "am_pathogenicity", "AM_logit", "predicted_AM_lm") # Updated kept columns
  keep <- intersect(keep, names(dt))
  
  lm_results <- list(
    fit = lm_fit,
    r_squared = summary(lm_fit)$r.squared,
    adjusted_r_squared = summary(lm_fit)$adj.r.squared,
    n = nrow(dt),
    predicted_data = dt[, ..keep]
  )
  
  # Save with _am suffix
  saveRDS(lm_results, file = file.path(output_dir, paste0(protein_id, "_lm_am.rds")))
  
  list(
    protein_id = protein_id,
    lm_success = TRUE,
    n = nrow(dt),
    r2 = lm_results$r_squared
  )
}

files <- list.files(merged_dir, pattern = "_merged\\.csv$", full.names = FALSE)
protein_ids <- sub("_merged\\.csv$", "", files)
length(protein_ids)

results <- rbindlist(lapply(protein_ids, process_protein_merged,
                            merged_dir = merged_dir, output_dir = output_dir),
                     fill = TRUE)


plan(multisession, workers = availableCores())  # Or set a specific number

# List RDS files
rds_files <- list.files(path = output_dir, pattern = "_lm_am.rds$", full.names = TRUE)

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
nrow(all_var_r2) # 20049
sum(!is.na(all_var_r2$adj_r2))  #20049

# Save the results to a CSV 
fwrite(all_var_r2, file.path("/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/00.human_proteome_results/all_var_am_lm_r2_jan2026.csv"))


#############################################







