# Load required libraries
library(dplyr)
library(tidyr)
library(data.table)
library(future.apply)
library(ggplot2)
library(cowplot)
library(parallel)
library(furrr)
library(jsonlite)
library(viridis)

# Set up parallel processing
future::plan(multisession)

combined_df <- read.csv("/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/clinvar_all_protein_filtered_with_esm1v_new.csv")
head(combined_df)
nrow(combined_df)
length(unique(combined_df$uniprot))

# List all JSON files in the directory (adjust the path as needed)
json_files <- list.files(path = "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_inputs/protein_atlas", pattern = "*.json", full.names = TRUE)

# Initialize an empty list to store data frames
df_list <- list()

# Loop through each JSON file
for (file in json_files) {
  # Read the JSON file
  json_data <- fromJSON(file)
  
  # Extract the 'Gene' and 'Gene description' fields
  gene_data <- json_data %>% dplyr::select(Gene, `Gene description`)
  
  # Add a column to specify the source JSON file
  gene_data <- gene_data %>% mutate(source_file = basename(file))
  
  # Append the data frame to the list
  df_list[[file]] <- gene_data
}

# Combine all data frames into one
protein_atlas <- bind_rows(df_list)
protein_atlas$source_file <- gsub("\\.json$", "", protein_atlas$source_file)

nrow(protein_atlas) #40031
head(protein_atlas)

# Save the combined data frame as a CSV file (optional)
write.csv(protein_atlas, "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_inputs/human_protein_atlas_anno.csv", row.names = FALSE)



id_mapping <- fread(
  cmd = "grep -P '\tGene_Name\t' /lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_inputs/uniprot_id_mapping/idmapping.dat",
  sep = "\t",
  header = FALSE,
  col.names = c("uniprot", "Type", "ID")
) %>%
  rename(
    gene_name = ID
  )

# Check the first few rows
head(id_mapping)
nrow(id_mapping) #58690748
length(unique(id_mapping$uniprot)) #4614909
tail(id_mapping)

fwrite(
  id_mapping,
  file = "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_inputs/uniprot_id_mapping/idmapping_gene_name.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  compress = "gzip"  # Optional: adds .gz and compresses on the fly
)

head(id_mapping)
merged_df <- merge(protein_atlas, id_mapping, by.x="Gene", by.y = "gene_name", all.x = TRUE)

head(merged_df)
nrow(merged_df) #12607374

write.csv(merged_df, "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_inputs/protein_atlas_meta/idmapping_gene_name.tsv", row.names = FALSE)


clinvar_core_fil <- read.csv("/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/00.human_proteome_results/clinvar_all_protein_filtered_with_esm1v_core_new.csv")
head(clinvar_core_fil)


core_r2_df <- merge(r_squared_dt_lm,merged_df, by.x="protein_id", by.y="uniprot", all.x = TRUE )
nrow(core_r2_df) #34174


