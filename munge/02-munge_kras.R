
########## DATA SANITY-CHECKS ##########

########################################
######### KRAS DDG check alignment #####
########################################
# kras_ddg <- read_excel('./data/paper_supplements/kras_chenchun/supplement_tables/media-5.xlsx', sheet = 2)
# kras_ddg
# min(kras_ddg$`mean_kcal/mol`)
# 
# 
# # 
# # # Function to process ddG data for each assay type
# process_ddg <- function(df, assay_type, new_ddg_name, new_ddg_std_name) {
#   df %>%
#     filter(assay == assay_type) %>%
#     dplyr::select(id, Pos_real, `mean_kcal/mol`, `std_kcal/mol`, wt_codon)
# }
# 
# # Ensure you drop NA in 'id' column before merging
# df_ddg_folding <- process_ddg(kras_ddg, "folding", "abundance_ddg", "abundance_ddg_std")
# df_ddg_folding <- df_ddg_folding[df_ddg_folding$id != "NA"& df_ddg_folding$id != "WT", ]
# min(as.numeric(df_ddg_folding$`mean_kcal/mol`))
# max(as.numeric(df_ddg_folding$`mean_kcal/mol`))
# 
# ggplot(df_ddg_folding, aes(x = `mean_kcal/mol`)) +
#   geom_density() +
#   labs(title = "Density Plot", x = "Value", y = "Density")

# # Process each dataset and create sequences
# kras_ddg_apca_unique <- ali_process_sequence(df_ddg_folding, "Pos_real", "wt_codon")
# kras_am_unique <- ali_process_sequence(kras_data$am, "position", "wt_aa")
# kras_esm1v_unique <- ali_process_sequence(kras_data$esm1v, "position", "wt_aa", "id")
# kras_esm1b_unique <- ali_process_sequence(kras_data$esm1b, "position", "wt_aa", "id")
# kras_popeve_unique <- ali_process_sequence(kras_data$popeve, "position", "wt_aa", "id")
# kras_eve_unique <- ali_process_sequence(kras_data$eve, "position", "wt_aa", "id")
# 
# # Find overlaps for each sequence
# ali_find_sequence_overlap(substr(kras_am_unique, 1, 150), substr(kras_ddg_apca_unique,1,148), "AM sequence")
# #"AM sequence starts at index 2 in the reference sequence."
# ali_find_sequence_overlap(kras_esm1v_unique, kras_ddg_apca_unique, "ESM1v sequence")
# #"ESM1v sequence starts at index 2 in the reference sequence."
# ali_find_sequence_overlap(substr(kras_esm1b_unique, 1, 150), substr(kras_ddg_apca_unique,1,148), "ESM1b sequence")
# #"ESM1b sequence starts at index 2 in the reference sequence."
# ali_find_sequence_overlap(substr(kras_popeve_unique, 1, 150), substr(kras_ddg_apca_unique,1,148), "popeve sequence")
# #"popeve sequence starts at index 2 in the reference sequence."
# ali_find_sequence_overlap(substr(kras_eve_unique, 1, 150), substr(kras_ddg_apca_unique,1,148), "eve sequence")
# #"eve sequence starts at index 2 in the reference sequence."
# 
# ########################################
# ######### KRAS DDG check alignment #####
# ########################################
# 
#kras_fit <- read_csv('./data/paper_supplements/kras_chenchun/supplement_tables/kras_fit_cleaned.csv')
# 
# process_fit <- function(df, assay_type, new_fit_name, new_fit_std_name) {
#   df %>%
#     filter(assay == assay_type) %>%
#     dplyr::select(id, nor_fitness, nor_fitness_sigma) %>%
#     dplyr::rename(!!new_fit_name := nor_fitness, !!new_fit_std_name := nor_fitness_sigma) %>%
#     mutate (Pos_real = as.integer(str_extract(id, '\\d+')),
#             wt_aa = substr(id, 1, 1),)
# }
# 
# # Define the assays and their new column names
# df_fit_folding <- process_fit(kras_fit, "AbundancePCA", "abundance_nor_fit", "abundance_nor_fit_std")
# 
# # Process each dataset and create sequences
# kras_fit_apca_unique <- ali_process_sequence(df_fit_folding, "Pos_real", "wt_aa")
# 
# # Find overlaps for each sequence
# ali_find_sequence_overlap(substr(kras_am_unique, 1, 70), substr(kras_fit_apca_unique,1,68), "AM sequence")
# #"AM sequence starts at index 2 in the reference sequence."
# # check alignment fig, there is a N missing in the fitness data but index is the same

########## DATA PREPROCESSING ##########

########################################
######### KRAS ddg Data ################
########################################
# kras_ddg <- read_excel('./data/paper_supplements/kras_chenchun/supplement_tables/media-5.xlsx', sheet = 2)
# 
# #Function to process ddG data for each assay type
# process_ddg <- function(df, assay_type, new_ddg_name, new_ddg_std_name) {
#   df %>%
#     filter(assay == assay_type) %>%
#     dplyr::select(id, Pos_real, `mean_kcal/mol`, `std_kcal/mol`, wt_codon)
# }
# 
# # Ensure you drop NA in 'id' column before merging
# df_ddg_folding <- process_ddg(kras_ddg, "folding", "abundance_ddg", "abundance_ddg_std")
# df_ddg_pik <- process_ddg(kras_ddg, "PIK3CG", "pik3cg_ddg", "pik3cg_ddg_std")
# df_ddg_raf1 <- process_ddg(kras_ddg, "RAF1", "raf1_ddg", "raf1_ddg_std")
# df_ddg_ral <- process_ddg(kras_ddg, "RALGDS", "ralgds_ddg", "ralgds_ddg_std")
# df_ddg_sos <- process_ddg(kras_ddg, "SOS1", "sos1_ddg", "sos1_ddg_std")
# 
# df_ddg_folding <- df_ddg_folding[df_ddg_folding$id != "NA"& df_ddg_folding$id != "WT", ]
# df_ddg_pik <- df_ddg_pik[df_ddg_pik$id != "NA" & df_ddg_pik$id != "WT",]
# df_ddg_raf1 <- df_ddg_raf1[df_ddg_raf1$id != "NA" & df_ddg_raf1$id != "WT",]
# df_ddg_ral <- df_ddg_ral[df_ddg_ral$id != "NA" & df_ddg_ral$id != "WT",]
# df_ddg_sos <- df_ddg_sos[df_ddg_sos$id != "NA" & df_ddg_sos$id != "WT",]
# 
# # Example call for one protein
kras_data <- process_protein_data(
  consurf_file = './data/scores/3_consurf/kras_consurf/new_cleaned_file.txt',
  am_file = './data/scores/0_am/alphamissense_kras.tsv',
  esm1b_file = './data/scores/1_esm1b/KRAS_P01116_esm1b.csv',
  popeve_file = './data/scores/2_popeve/KRAS_NP_004976.2.csv',
  eve_file = './data/scores/4_eve/EVE_RASK_HUMAN/RASK_HUMAN.csv'
)
# 
process_and_merge_protein_data_esm1v <- function(protein_data, merged_df) {

  # Merge with ESM1v data
  merged_df_ddg <- merged_df %>%
    merge(protein_data$esm1v, by = "id", all.x = FALSE)

  return(merged_df_ddg)
}
# 
# kras_ddg_apca_esm1v <- process_and_merge_protein_data_esm1v(kras_data, df_ddg_folding)
# write.csv(kras_ddg_apca_esm1v, './data/cleaned/kras_ddg_apca_esm1v_full.csv', row.names = FALSE)
# 
# # Merge the processed DataFrames into one final DataFrame
# kras_merged_df <- df_ddg_folding %>%
#  inner_join(df_ddg_pik, by = "id") %>%
#  inner_join(df_ddg_raf1, by = "id") %>%
#  inner_join(df_ddg_ral, by = "id") %>%
#  inner_join(df_ddg_sos, by = "id") %>%
#  mutate(Pos_real = coalesce(Pos_real.x, Pos_real.y, Pos_real.x.x, Pos_real.y.y)) %>%
# Select desired columns and remove extra Pos_real columns
#  select(id, Pos_real, abundance_ddg, abundance_ddg_std, pik3cg_ddg, pik3cg_ddg_std, raf1_ddg, raf1_ddg_std, ralgds_ddg, ralgds_ddg_std, sos1_ddg, sos1_ddg_std)
# 
# # Inspect the final merged dataframe
# #print(kras_merged_df)
# 
# 
# # Example call 
# kras_final_df_ddg <- process_and_merge_protein_data(kras_data, kras_merged_df)
# dim(kras_final_df)
# 
# # Inspect the final merged dataframe
# # kras_final_df
# 
# # Check for NA values in each column
# #colSums(is.na(kras_final_df))
# 
# # Load anno
# kras_anno = read_csv('./data/paper_supplements/kras_chenchun/DATA/anno_final3.csv')
# kras_anno_subset <- kras_anno[c("Pos","bp_interface","binding_RAF","binding_RAL","binding_PI3","binding_SOS")]
# kras_final_df_ddg <- merge(kras_final_df_ddg, kras_anno_subset, by.x="Pos_real", by.y="Pos")
# 
# write.csv(kras_final_df, './data/cleaned/kras_ddg_cleaned.csv', row.names = FALSE)

kras_final_df_ddg <- read_csv('./data/cleaned_ddg/kras_ddg_cleaned.csv')

# Convert relevant columns to numeric, coercing non-numeric values to NA
kras_final_df_ddg <- kras_final_df_ddg %>%
  mutate(
    abundance_ddg = as.numeric(abundance_ddg),
    pik3cg_ddg = as.numeric(pik3cg_ddg),
    raf1_ddg = as.numeric(raf1_ddg),
    ralgds_ddg = as.numeric(ralgds_ddg),
    sos1_ddg = as.numeric(sos1_ddg),
    abundance_ddg_std = as.numeric(abundance_ddg_std),
    pik3cg_ddg_std = as.numeric(pik3cg_ddg_std),
    raf1_ddg_std = as.numeric(raf1_ddg_std),
    ralgds_ddg_std = as.numeric(ralgds_ddg_std),
    sos1_ddg_std = as.numeric(sos1_ddg_std)
  )

kras_ddg_apca_esm1v_full <- read_csv('./data/cleaned_ddg/kras_ddg_apca_esm1v_full.csv')
#dim(kras_ddg_apca_esm1v_full) # 3453 6 
# only KRAS needs this as the other proteins align fully 

########################################
######### KRAS fitness Data ############
########################################

# # Read Excel file 
# kras_df_fitness <- read_excel('./data/paper_supplements/kras_chenchun/supplement_tables/media-4.xlsx', sheet = 2)
# 
# # Extract the WT sequence (equivalent to df_fitness.loc)
# wt_seq <- kras_df_fitness %>%
#   filter(WT == TRUE) %>%
#   pull(aa_seq) %>%
#   .[1]  # Select the first element
# 
# # Filter rows where 'Nham_aa' equals 1
# kras_df_fitness <- kras_df_fitness %>%
#   filter(Nham_aa == 1)
# 
# # Test the function
# # compare_sequences("AA", "A*")
# # compare_sequences(wt_seq, kras_df_fitness[grepl("\\*", kras_df_fitness$aa_seq), ][1,]$aa_seq)
# 
# # Apply the comparison function to each row in the data frame
# kras_df_fitness <- kras_df_fitness %>%
#  mutate(id = sapply(aa_seq, function(x) compare_sequences_kras(wt_seq, x)))
# # 
# # kras_df_fitness[grepl("\\*", kras_df_fitness$aa_seq), ]
# 
# #dim(kras_df_fitness) #21811 44
# #colSums(is.na(kras_df_fitness)) # no NA!
# 
# #min(kras_df_fitness$nor_fitness) #-1.581651
# #max(kras_df_fitness$nor_fitness) #0.8920027
# 
# write.csv(kras_df_fitness, './data/paper_supplements/kras_chenchun/supplement_tables/kras_fit_cleaned.csv', row.names = FALSE)

# kras_fit <- read_csv('./data/paper_supplements/kras_chenchun/supplement_tables/kras_fit_cleaned.csv')
# # kras_fit[grepl("\\*", kras_fit$aa_seq), ]
# 
# # Function to process fitness data for each assay type
# process_fit <- function(df, assay_type, new_fit_name, new_fit_std_name) {
#  df %>%
#    filter(assay == assay_type) %>%
#    dplyr::select(id, nor_fitness, nor_fitness_sigma) %>%
#    dplyr::rename(!!new_fit_name := nor_fitness, !!new_fit_std_name := nor_fitness_sigma) %>%
#    distinct()
# }
# 
# # Define the assays and their new column names
# df_fit_folding <- process_fit(kras_fit, "AbundancePCA", "abundance_nor_fit", "abundance_nor_fit_std")
# df_fit_pik <- process_fit(kras_fit, "BindingPCA PIK3CG", "pik3cg_nor_fit", "pik3cg_nor_fit_std")
# df_fit_raf1 <- process_fit(kras_fit, "BindingPCA RAF1", "raf1_nor_fit", "raf1_nor_fit_std")
# df_fit_ral <- process_fit(kras_fit, "BindingPCA RALGDS", "ralgds_nor_fit", "ralgds_nor_fit_std")
# df_fit_sos <- process_fit(kras_fit, "BindingPCA SOS1", "sos1_nor_fit", "sos1_nor_fit_std")
# 
# 
# kras_fit_apca_esm1v <- process_and_merge_protein_data_esm1v(kras_data, df_fit_folding)
# kras_fit_apca_esm1v <- kras_fit_apca_esm1v  %>%
#     mutate (Pos_real = as.integer(str_extract(id, '\\d+')),
#       wt_aa = substr(id, 1, 1),) %>%
#     arrange(as.numeric(Pos_real))
# 
# dim(kras_fit_apca_esm1v) # 3135    4
# write.csv(kras_fit_apca_esm1v, './data/cleaned/kras_fit_apca_esm1v_full.csv', row.names = FALSE)
# 
# 
# # List of datasets
# fit_datasets <- list(
#  df_fit_folding = df_fit_folding,
#  df_fit_pik = df_fit_pik,
#  df_fit_raf1 = df_fit_raf1,
#  df_fit_ral = df_fit_ral,
#  df_fit_sos = df_fit_sos
# )
# 
# # Corresponding column names for fitness and standard deviation
# fit_columns <- list(
#  c("abundance_nor_fit", "abundance_nor_fit_std"),
#  c("pik3cg_nor_fit", "pik3cg_nor_fit_std"),
#  c("raf1_nor_fit", "raf1_nor_fit_std"),
#  c("ralgds_nor_fit", "ralgds_nor_fit_std"),
#  c("sos1_nor_fit", "sos1_nor_fit_std")
# )
# 
# # Function to summarize fitness and standard deviation
# summarize_fitness <- function(df, fit_col, std_col) {
#  df %>%
#    group_by(id) %>%
#    summarise(
#      !!fit_col := mean(!!sym(fit_col)), # need to take average as certain positions were covered by both blocks
#      !!std_col := mean(!!sym(std_col))
#    )
# }
# 
# # Automate the summarization for all datasets
# summarized_fit <- map2(
#  fit_datasets,
#  fit_columns,
#  ~ summarize_fitness(.x, .y[1], .y[2])
# )
# 
# kras_merged_df_fit <- summarized_fit$df_fit_folding %>%
#  inner_join(summarized_fit$df_fit_pik, by = "id") %>%
#  inner_join(summarized_fit$df_fit_raf1, by = "id") %>%
#  inner_join(summarized_fit$df_fit_ral, by = "id") %>%
#  inner_join(summarized_fit$df_fit_sos, by = "id") %>%
#  dplyr::select(id, abundance_nor_fit, abundance_nor_fit_std,
#         pik3cg_nor_fit, pik3cg_nor_fit_std,
#         raf1_nor_fit, raf1_nor_fit_std,
#         ralgds_nor_fit, ralgds_nor_fit_std,
#         sos1_nor_fit, sos1_nor_fit_std) %>%
#   mutate (Pos_real = as.integer(str_extract(id, '\\d+')),
#           wt_aa = substr(id, 1, 1),) %>%
#   arrange(as.numeric(Pos_real))
# 
# dim(kras_merged_df_fit) #2563   13
# 
# # kras_merged_df_fit[grepl("\\*", kras_merged_df_fit$id), ]
# 
# #fit_unique <- kras_merged_df_fit %>%  distinct(Pos_real, .keep_all = TRUE) %>%  filter(!is.na(Pos_real)) %>%  arrange(Pos_real)
# #sequence <- paste(fit_unique$wt_aa, collapse = "")
# #print(sequence)
# 
# kras_final_df_fit <- process_and_merge_protein_data(kras_data, kras_merged_df_fit)
# dim(kras_final_df_fit) #2339   29
# # # the merge will get rid of all stop codon muts
# kras_final_df_fit <- merge(kras_final_df_fit, kras_anno_subset, by.x="Pos_real", by.y="Pos")
# 
# write.csv(kras_final_df_fit, './data/cleaned/kras_fit_cleaned.csv', row.names = FALSE)

# kras_final_df_fit <- read_csv('./data/cleaned_ddg/kras_fit_cleaned.csv')
# 
# # Convert relevant columns to numeric, coercing non-numeric values to NA
# kras_final_df_fit <- kras_final_df_fit %>%
#   mutate(
#     abundance_nor_fit = as.numeric(abundance_nor_fit),
#     pik3cg_nor_fit = as.numeric(pik3cg_nor_fit),
#     raf1_nor_fit = as.numeric(raf1_nor_fit),
#     ralgds_nor_fit = as.numeric(ralgds_nor_fit),
#     sos1_nor_fit = as.numeric(sos1_nor_fit),
#     abundance_nor_fit_std = as.numeric(abundance_nor_fit_std),
#     pik3cg_nor_fit_std = as.numeric(pik3cg_nor_fit_std),
#     raf1_nor_fit_std = as.numeric(raf1_nor_fit_std),
#     ralgds_nor_fit_std = as.numeric(ralgds_nor_fit_std),
#     sos1_nor_fit_std = as.numeric(sos1_nor_fit_std)
#   )
# 
# kras_fit_apca_esm1v_full <- read_csv('./data/cleaned_ddg/kras_fit_apca_esm1v_full.csv')
# only KRAS needs this as the other proteins align fully 

########################################
######### KRAS RSA Data ################
########################################
# Use the function to read your file
kras_rsa_file <- "./data/sasa/6vjj_rsa.txt"
kras_rsa_df <- read_sasa_file(kras_rsa_file)
kras_rsa_df <- kras_rsa_df[kras_rsa_df$Chain == "A",]






