
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



# Set up parallel processing
plan(multisession, workers = availableCores())  # Or set a specific number

process_and_merge_clinvar_with_esm1v <- function(
    csv_dir = "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/meta",
    esm1v_dir = "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_inputs/unzipped/all_ESM1v",
    output_file = "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/clinvar_all_protein_filtered_with_esm1v_new.csv"
) {
  # Set up parallel processing
  plan(multisession, workers = availableCores())
  
  # List CSV files to process
  csv_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)
  total_files <- length(csv_files)
  
  clinical_significance_filter <- c("benign", "likely_benign", "pathogenic", "likely_pathogenic")
  
  # Function to read and filter each file
  read_file <- function(file, index) {
    dt <- fread(file)
    filtered_dt <- dt[clinvar_clinical_significance %in% clinical_significance_filter]
    cat(sprintf("\rProcessing file %d of %d (%.2f%%)", index, total_files, (index / total_files) * 100))
    flush.console()
    return(filtered_dt)
  }
  
  # Read and combine all files in parallel
  combined_df <- future_lapply(seq_along(csv_files), function(i) read_file(csv_files[i], i), future.seed = TRUE) |> rbindlist(fill = TRUE)
  cat("\nAll files processed successfully!\n")
  
  # Map ESM1v files
  esm1v_files <- list.files(esm1v_dir, pattern = "*.csv", full.names = TRUE)
  esm1v_map <- setNames(esm1v_files, str_replace(basename(esm1v_files), ".csv$", ""))
  
  # Merge function
  merge_esm1v_to_combined_df <- function(combined_df, esm1v_map) {
    combined_split <- split(combined_df, combined_df$uniprot)
    merged_data_list <- list()
    all_columns <- colnames(combined_df)
    
    for (uniprot_id in names(combined_split)) {
      protein_data <- combined_split[[uniprot_id]]
      if (uniprot_id %in% names(esm1v_map)) {
        df_esm1v <- read.csv(esm1v_map[[uniprot_id]], stringsAsFactors = FALSE)
        if ("variant" %in% colnames(df_esm1v)) {
          protein_data <- protein_data %>%
            left_join(df_esm1v, by = "variant")
          all_columns <- unique(c(all_columns, colnames(protein_data)))
          message("Merged ESM1v data for: ", uniprot_id)
        } else {
          message("No 'variant' column found in ESM1v file for: ", uniprot_id)
        }
      } else {
        message("No ESM1v file found for: ", uniprot_id)
      }
      merged_data_list[[uniprot_id]] <- protein_data
    }
    
    for (name in names(merged_data_list)) {
      missing_cols <- setdiff(all_columns, colnames(merged_data_list[[name]]))
      if (length(missing_cols) > 0) {
        merged_data_list[[name]][, missing_cols] <- NA
      }
    }
    
    combined_merged_df <- do.call(rbind, merged_data_list)
    return(combined_merged_df)
  }
  
  # Merge and save
  combined_merged_df <- merge_esm1v_to_combined_df(combined_df, esm1v_map)
  write.csv(combined_merged_df, output_file, row.names = FALSE)
  cat("\nMerged data saved to: ", output_file, "\n")
  
  return(combined_merged_df)
}

final_df <- process_and_merge_clinvar_with_esm1v()

final_df <- read.csv("/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/clinvar_all_protein_filtered_with_esm1v_new.csv")
nrow(final_df) #114903
length(unique(final_df$uniprot))
head(final_df)

# check <- fread("/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_inputs/cagiada_2024_supp/20250128_functional_proteome_clinvar.csv")
# nrow(check)
# head(check)
# check
# 
# nrow(check %>% filter (clinvar_clinical_significance %in% c("pathogenic","likely_pathogenic","benign","likely_benign"))) #114903
# 
total_patho <- final_df %>%
  dplyr::filter(clinvar_clinical_significance %in% c("pathogenic", "likely_pathogenic"))
nrow(total_patho) #43014

total_benign <- final_df %>%
  dplyr::filter(clinvar_clinical_significance %in% c("benign", "likely_benign"))
nrow(total_benign) #71889

final_df <- final_df %>%
  mutate(stability = case_when(
    ddG_pred > 1 ~ "destabilizing",
    ddG_pred < -0.5 ~ "stabilizing",
    TRUE ~ "WT-like"  # Default case
  ))


# Summarize data by categories and pathogenicity
summary_df <- final_df %>%
  mutate(clinvar_clinical_significance = case_when(
    clinvar_clinical_significance %in% c("likely_benign", "benign") ~ "benign",
    clinvar_clinical_significance %in% c("likely_pathogenic", "pathogenic") ~ "pathogenic",
    TRUE ~ clinvar_clinical_significance  # Keep other categories unchanged
  )) %>%
  # Group by the new clinical significance categories and stability
  group_by(clinvar_clinical_significance, stability) %>%
  summarise(count = n(), .groups = "drop") %>%
  # Recalculate percentages after merging
  group_by(clinvar_clinical_significance) %>%
  mutate(percent = (count / sum(count)) * 100) %>%
  ungroup()

summary_df


# Creating a contingency table
contingency_table <- matrix(c(
  summary_df$count[summary_df$clinvar_clinical_significance == "pathogenic" & summary_df$stability == "destabilizing"], 
  summary_df$count[summary_df$clinvar_clinical_significance == "benign" & summary_df$stability == "destabilizing"],
  sum(summary_df$count[summary_df$clinvar_clinical_significance == "pathogenic" & summary_df$stability != "destabilizing"]),
  sum(summary_df$count[summary_df$clinvar_clinical_significance == "benign" & summary_df$stability != "destabilizing"])
), nrow = 2, byrow = TRUE)

# Adding row and column names for clarity
colnames(contingency_table) <- c("Pathogenic", "Benign")
rownames(contingency_table) <- c("Destabilizing", "Non-Destabilizing")

# Displaying the table
print(contingency_table)
# Pathogenic Benign
# Destabilizing          20815   6319
# Non-Destabilizing      22199  65570

# Calculate the Odds Ratio using Fisher's Test
result <- fisher.test(contingency_table)
print(result$estimate)  # This is the Odds Ratio
# odds ratio 
# 9.727316 


table(combined_df$stability)

summary_df$stability <- factor(summary_df$stability, levels = c("WT-like", "destabilizing", "stabilizing"))
summary_df$clinvar_clinical_significance <- factor(summary_df$clinvar_clinical_significance, levels = c("benign", "pathogenic"))

nrow(final_df) #114903
length(unique(final_df$uniprot)) #12372

ggplot(summary_df, aes(x = clinvar_clinical_significance, y = percent, fill = stability)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  # Stacked bars
  geom_text(aes(label = paste0(round(percent, 1), "% (", count, ")")), 
            position = position_stack(vjust = 0.5), size = 6.5) +  # Center labels in stacks
  scale_fill_brewer(palette = "Set2") +
  labs(title = "114,903 ClinVar Variants", subtitle = "12,372 proteins",
       x = NULL, y = "Variant percentage (%)", fill = NULL) +
  theme_classic() + ylim(0, 100) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 24, hjust = 0, face = "bold"),
        plot.subtitle = element_text(size = 22, hjust = 0),
        axis.text.y = element_text(size = 20) ,
        axis.text.x = element_text(size = 20) ,
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 20),       # Increase legend text size
        legend.title = element_text(size = 20),      # Increase legend title size
        legend.key.size = unit(1.2, "cm"))            # Increase legend box size


bi_df <- read.csv("/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_inputs/pioneer_interactome/all_bi_residues.csv")
head(bi_df)

nrow(final_df) #114903
final_df2 <- merge(final_df, bi_df, by = "uniprot", all.x = TRUE)
nrow(final_df2) #114903

clinvar_core <- final_df2 %>% filter (exposure_rasa <= 0.2) %>% filter (spot_disorder == 0)
nrow(clinvar_core) #32727
head(clinvar_core)

# Step 1: Convert interface_residues to a list of numeric values
clinvar_core <- clinvar_core %>%
  mutate(
    interface_residues = str_split(interface_residues, ",\\s*"), # Split into individual residues
    interface_residues = lapply(interface_residues, as.numeric) # Convert to numeric
  )

# Step 2: Filter rows where new_position is not in interface_residues
clinvar_core_fil <- clinvar_core %>%
  rowwise() %>%
  filter(!(new_position %in% interface_residues)) %>%
  ungroup()

clinvar_core_fil <- clinvar_core_fil %>%
  mutate(interface_residues = sapply(interface_residues, function(x) paste(x, collapse = ",")))

nrow(clinvar_core_fil) #30744
head(clinvar_core_fil)

write.csv(clinvar_core_fil, "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/clinvar_all_protein_filtered_with_esm1v_core_new.csv", row.names = FALSE)


nrow(clinvar_core_fil) #19491

summary_df <- clinvar_core_fil %>%
  mutate(clinvar_clinical_significance = case_when(
    clinvar_clinical_significance %in% c("likely_benign", "benign") ~ "benign",
    clinvar_clinical_significance %in% c("likely_pathogenic", "pathogenic") ~ "pathogenic",
    TRUE ~ clinvar_clinical_significance  # Keep other categories unchanged
  )) %>%
  # Group by the new clinical significance categories and stability
  group_by(clinvar_clinical_significance, stability) %>%
  summarise(count = n(), .groups = "drop") %>%
  # Recalculate percentages after merging
  group_by(clinvar_clinical_significance) %>%
  mutate(percent = (count / sum(count)) * 100) %>%
  ungroup()

summary_df

# Set factor levels for consistent ordering
summary_df$stability <- factor(summary_df$stability, levels = c("WT-like", "destabilizing", "stabilizing"))
summary_df$clinvar_clinical_significance <- factor(summary_df$clinvar_clinical_significance, levels = c("benign", "pathogenic"))

nrow(clinvar_core_fil) #30744
length(unique(clinvar_core_fil$uniprot)) #5300

ggplot(summary_df, aes(x = clinvar_clinical_significance, y = percent, fill = stability)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  # Stacked bars
  geom_text(aes(label = paste0(round(percent, 1), "% (", count, ")")), 
            position = position_stack(vjust = 0.5), size = 6.5) +  # Center labels in stacks
  scale_fill_brewer(palette = "Set2") +
  labs(title = "30,744 ClinVar Variants", subtitle = "5,300 proteins",
       x = NULL, y = "Variant percentage (%)", fill = NULL) +
  theme_classic() + ylim(0, 100) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 24, hjust = 0, face = "bold"),
        plot.subtitle = element_text(size = 22, hjust = 0),
        axis.text.y = element_text(size = 20) ,
        axis.text.x = element_text(size = 20) ,
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 20),       # Increase legend text size
        legend.title = element_text(size = 20),      # Increase legend title size
        legend.key.size = unit(1.2, "cm"))            # Increase legend box size


################################################################################
protein_atlas <- read_csv("/lustre/scratch126/gengen/projects/alpha-allostery-global/00.human_proteome_protein_atlas/combined_data.csv")

# Read the UniProt ID mapping file
id_mapping <- read.delim("/lustre/scratch126/gengen/projects/alpha-allostery-global/00.human_proteome_protein_atlas/uniprot_id_mapping.dat", header = TRUE, sep = "\t", stringsAsFactors = FALSE)






























