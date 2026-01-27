library(data.table)
library(parallel)
library(purrr)
library(furrr)

# Directory with CSVs
dir_path <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/00.human_proteome_unzipped/results_unzipped/merged_out_dec2025"

# Get list of non-empty CSV files
files <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE)
non_empty_files <- files[file.info(files)$size > 0]
length(non_empty_files)

# Function to classify file issues
check_file <- function(f) {
  dat <- tryCatch(fread(f, showProgress = FALSE), error = function(e) return("error"))
  if (is.character(dat) && dat == "error") return("read_error")
  if (!"am_pathogenicity" %in% colnames(dat)) return("missing_column")
  if (any(is.na(dat$am_pathogenicity))) return("has_na")
  return("ok")
}

# Run in parallel
num_cores <- max(1, detectCores() - 1)
statuses <- mclapply(non_empty_files, check_file, mc.cores = num_cores)

# Map results
status_table <- data.table(file = non_empty_files, status = unlist(statuses))

# Summarize
summary_counts <- status_table[, .N, by = status]
print(summary_counts)
#status     N
#<char> <int>
#1:     ok 20422
#2: has_na   238

status_table %>% filter (status == "has_na")

# destination directory
out_dir <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/00.human_proteome_unzipped/results_unzipped/merged_no_am_jan2026"

# files to copy
files_ok <- status_table[status != "has_na", file]
length(files_ok)

# copy
file.copy(
  from = files_ok,
  to   = file.path(out_dir, basename(files_ok)),
  overwrite = TRUE
)

##################################################################################
# check for NA ddG_pred 
# Directory with CSVs
dir_path <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/00.human_proteome_unzipped/results_unzipped/merged_no_am_jan2026"

# Get list of non-empty CSV files
files <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE)
non_empty_files <- files[file.info(files)$size > 0]

# Function to classify file issues
check_file <- function(f) {
  dat <- tryCatch(fread(f, showProgress = FALSE), error = function(e) return("error"))
  if (is.character(dat) && dat == "error") return("read_error")
  if (!"ddG_pred" %in% colnames(dat)) return("missing_column")
  if (any(is.na(dat$ddG_pred))) return("has_na")
  return("ok")
}

# Run in parallel
num_cores <- max(1, detectCores() - 1)
statuses <- mclapply(non_empty_files, check_file, mc.cores = num_cores)

# Map results
status_table <- data.table(file = non_empty_files, status = unlist(statuses))

# Summarize
summary_counts <- status_table[, .N, by = status]
print(summary_counts)
#1:             ok 19532
#2:         has_na   517
#3: missing_column   373

status_table %>% filter (status == "has_na")

# destination directory
out_dir <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/00.human_proteome_unzipped/results_unzipped/merged_am_ddgf_jan2026"

# files to copy
files_ok <- status_table[status != "missing_column", file]
length(files_ok) #20049

# copy
file.copy(
  from = files_ok,
  to   = file.path(out_dir, basename(files_ok)),
  overwrite = TRUE
)

################# check AM all
# Set the directory path where your RDS files are stored
dir_path <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/00.human_proteome_results/am_lm_rds_jan2026"

# List all RDS files in the directory
rds_files <- list.files(path = dir_path, pattern = "\\_lm_am.rds$", full.names = TRUE)

safe_read_and_count <- safely(function(file) {
  lm_rds <- readRDS(file)
  
  # Check if 'n' exists in the list and return it directly
  if (!is.null(lm_rds$n)) {
    lm_rds$n
  } else {
    0
  }
})


# 1. Run on ALL files (removed the [1] index)
results <- future_map(rds_files, safe_read_and_count)

# 2. Extract row counts
# We use map_dbl instead of map_int to avoid errors if 'n' is stored as a numeric type
row_counts <- map_dbl(results, ~ if (is.null(.x$error)) .x$result else 0)

# 3. Sum the total
total_count <- sum(row_counts)
total_count
#194791762

######## check AM core
# Set the directory path where your RDS files are stored
dir_path <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/00.human_proteome_unzipped/results_unzipped/AM_results/core_lm_rds"

# List all RDS files in the directory
rds_files <- list.files(path = dir_path, pattern = "\\_lm.rds$", full.names = TRUE)

get_row_count <- function(file) {
  # Read the file
  lm_rds <- readRDS(file)
  
  # Check if predicted_data exists and return count
  if (!is.null(lm_rds$predicted_data)) {
    return(nrow(lm_rds$predicted_data))
  } else {
    return(0)
  }
}

results <- future_map_int(rds_files, get_row_count)
readRDS(rds_files[1])
results[1]

sum(results = 0)
sum(results > 0) #17519

# 3. Sum the total
total_count <- sum(results)
total_count
#47982923










