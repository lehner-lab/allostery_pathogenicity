library(stringr)
library(biomaRt)
library(data.table)

# Path containing the merged CSVs
data_dir <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/00.human_proteome_unzipped/results_unzipped/merged_out_dec2025"

# List all *_merged.csv files
files <- list.files(
  path = data_dir,
  pattern = "_merged\\.csv$",
  full.names = TRUE
)

length(files)  # sanity check (~20,000)

# Read and combine
master_dt <- rbindlist(
  lapply(files, function(f) {
    dt <- fread(f)
    dt[, source_file := basename(f)]  # optional but VERY useful
    dt
  }),
  use.names = TRUE,
  fill = TRUE
)

# Write combined master CSV
fwrite(
  master_dt,
  file = file.path("/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/00.human_proteome_results/clinvar_all_protein_filtered_with_esm1v_2026.csv")
)

combined_df <- read.csv("/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/00.human_proteome_results/clinvar_all_protein_filtered_with_esm1v_2026.csv")
head(combined_df)

nrow(combined_df) #217813060

unique(combined_df$clinvar_review_status)

combined_df_clean <- combined_df %>% filter (!clinvar_review_status %in% c("0", "0,0", "0,0,0"))
nrow(combined_df_clean) #217792789

combined_df_by_gene <- combined_df_clean %>%
  filter(clinvar_clinical_significance %in% c("pathogenic", "likely_pathogenic")) %>% # Filter pathogenic variants
  group_by(Dataset) %>% # Group by gene
  summarise(
    pathogenic_count = n(), # Count number of pathogenic mutations per gene
    variants = paste(unique(variant), collapse = ", ")
  ) %>%
  arrange(desc(pathogenic_count)) # Rank genes by number of pathogenic mutations


combined_df_by_gene

length(unique(combined_df_by_gene$Dataset)) #2749

# Extract UniProt ID from Dataset
combined_df_by_gene<- combined_df_by_gene %>%
  mutate(UniProt_ID = str_extract(Dataset, "(?<=AF-)[^-]+")) %>%
  dplyr::rename(AFDB_id = Dataset)

# ensembl <- useEnsembl(
#   biomart = "ensembl",
#   dataset = "hsapiens_gene_ensembl",
#   mirror = "useast"  # Use the US East mirror site
# )

# hgnc_mapping <- getBM(
#   attributes = c("uniprot_gn_id", "hgnc_symbol"),
#   filters = "uniprot_gn_id",
#   values = combined_df_by_gene$UniProt_ID,
#   mart = ensembl
# )

# Merge HGNC symbols into combined_df
# combined_df_by_gene <- combined_df_by_gene %>%
#   left_join(hgnc_mapping, by = c("UniProt_ID" = "uniprot_gn_id"))

nrow(combined_df_by_gene) #3773

head(combined_df_by_gene)

write.csv(combined_df_by_gene, "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/git_allostery_pathogenicity/data/gene_ranked_by_likely_pathogenic_variants_2026.csv", row.names = FALSE)


