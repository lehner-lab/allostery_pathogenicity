library(readr)
library(writexl)
library(purrr)
library(stringr)

# 1. Define the directory path (using the path you provided)
data_dir <- "~/Documents/0.Projects/01.protein-seq-evo-v1/data/source_data"

# 2. List all .csv files in that directory
file_list <- list.files(path = data_dir, pattern = "\\.csv$", full.names = TRUE)

# 3. Read the files into a named list
# The names of the list elements will become the Excel tab names
data_list <- file_list %>%
  # Create names by taking the filename and removing the path and ".csv" extension
  set_names(nm = basename(.) %>% str_replace("\\.csv$", "")) %>%
  # Read each CSV into a dataframe
  map(~read_csv(.x, show_col_types = FALSE))

# 4. Write the list to a single .xlsx file
output_file <- "~/Documents/0.Projects/01.protein-seq-evo-v1/data/source_data/consolidated.xlsx"
write_xlsx(data_list, path = output_file)

cat("Success! All CSVs have been merged into:", output_file, "\n")

R.version.string

# List of packages used
pkgs <- c("ggplot2")

# Check versions for all of them
sapply(pkgs, function(x) as.character(packageVersion(x)))
