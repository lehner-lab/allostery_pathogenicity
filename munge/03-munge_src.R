########## DATA PREPROCESSING ##########

# ########################################
# ######### SRC DDG check alignment #####
# ########################################
# src_ddg <- read_delim('./data/paper_supplements/src_toni/Supplementary_table_2_activity_folding_ddGs.txt', delim = ' ')
# 
# src_ddg_prep <- src_ddg %>%  
#   mutate (Pos_real = as.integer(str_extract(id, '\\d+')),
#           wt_aa = substr(id, 1, 1),) 
# 
# # Process each dataset and create sequences
# src_ddg_unique <- ali_process_sequence(src_ddg_prep, "Pos_real", "wt_aa")
# src_am_unique <- ali_process_sequence(src_data$am, "position", "wt_aa")
# src_esm1v_unique <- ali_process_sequence(src_data$esm1v, "position", "wt_aa", "id")
# src_esm1b_unique <- ali_process_sequence(src_data$esm1b, "position", "wt_aa", "id")
# src_popeve_unique <- ali_process_sequence(src_data$popeve, "position", "wt_aa", "id")
# src_eve_unique <- ali_process_sequence(src_data$eve, "position", "wt_aa", "id")
# 
# 
# # Find overlaps for each sequence
# ali_find_sequence_overlap(src_am_unique, src_ddg_unique,"AM sequence")
# #"AM sequence starts at index 268 in the reference sequence."
# ali_find_sequence_overlap(src_esm1v_unique, src_ddg_unique,"ESM1v sequence")
# #"ESM1v sequence starts at index 268 in the reference sequence."
# ali_find_sequence_overlap(src_esm1b_unique, src_ddg_unique,"ESM1b sequence")
# #"ESM1b sequence starts at index 268 in the reference sequence."
# ali_find_sequence_overlap(src_popeve_unique, src_ddg_unique,"popEVE sequence")
# #"popEVE sequence starts at index 268 in the reference sequence."
# ali_find_sequence_overlap(src_eve_unique, src_ddg_unique,"EVE sequence")
# #"EVE sequence starts at index 268 in the reference sequence."
# 
# 
# ########################################
# ###### SRC fitness check alignment #####
# ########################################
# src_df_fitness_cleaned <- read_csv('./data/cleaned/src_fit_cleaned')
# src_fit_prep <- src_df_fitness_cleaned %>%  
#   mutate (Pos_real = as.integer(str_extract(id, '\\d+')),
#           wt_aa = substr(id, 1, 1),) 
# 
# # Process each dataset and create sequences
# src_fit_unique <- ali_process_sequence(src_fit_prep, "Pos_real", "wt_aa")
# 
# # Find overlaps for each sequence
# ali_find_sequence_overlap(src_am_unique, src_fit_unique,"AM sequence")
# #"AM sequence starts at index 268 in the reference sequence."
# ali_find_sequence_overlap(src_esm1v_unique, src_fit_unique,"ESM1v sequence")
# #"ESM1v sequence starts at index 268 in the reference sequence."
# ali_find_sequence_overlap(src_esm1b_unique, src_fit_unique,"ESM1b sequence")
# #"ESM1b sequence starts at index 268 in the reference sequence."
# ali_find_sequence_overlap(src_popeve_unique, src_fit_unique,"popEVE sequence")
# #"popEVE sequence starts at index 268 in the reference sequence."
# ali_find_sequence_overlap(src_eve_unique, src_fit_unique,"EVE sequence")
# #"EVE sequence starts at index 268 in the reference sequence."


########################################
######### SRC DDG Data #################
########################################
src_ddg <- read_delim('./data/paper_supplements/src_toni/Supplementary_table_2_activity_folding_ddGs.txt', delim = ' ')

src_merged_df <- src_ddg %>%
  mutate(wt_aa = str_extract(id, '^[A-Za-z]+'),
         position_old = str_extract(id, '\\d+'),
         mt_aa = str_extract(id, '[A-Za-z]+$')) %>%
  mutate(position_old = as.numeric(position_old)) %>% 
  rename(old_id = id) %>% 
  mutate(Pos_real = position_old + 267) %>% 
  mutate(id = paste0(wt_aa, Pos_real, mt_aa))

src_data <- process_protein_data(
  consurf_file = './data/scores/3_consurf/src_consurf/new_cleaned_file.txt',
  am_file = './data/scores/0_am/alphamissense_src_p12931.tsv',
  esm1b_file = './data/scores/1_esm1b/SRC_P12931_esm1b.csv',
  popeve_file = './data/scores/2_popeve/SRC_NP_938033.1.csv'
)

src_final_df_ddg <- process_and_merge_protein_data(src_data,src_merged_df)

# Load anno
src_anno <- read_csv('./data/paper_supplements/src_toni/c-SRC_anno.csv')

src_final_df_ddg_anno <- merge(src_final_df_ddg, src_anno, by.x="Pos_real", by.y = "position")

# Define the active site positions
active_site_pos <- c(277, 279, 281, 282, 284, 296, 298, 389, 391, 394, 396, 407, 428)

# Initialize the 'active site' column with NA
src_final_df_ddg_anno$active_site <- NA

# Update the 'active site' column for the specified positions to TRUE
src_final_df_ddg_anno$active_site[src_final_df_ddg_anno$Pos_real %in% active_site_pos] <- TRUE

# Fill the remaining with FALSE
src_final_df_ddg_anno$active_site[is.na(src_final_df_ddg_anno$active_site)] <- FALSE

########################################
######### SRC fitness Data #############
########################################
# src_df_fitness <- read_delim('./data/paper_supplements/src_toni/media-3.txt', delim = '\t')
# #dim(src_df_fitness) # 56855 16
# src_df_fitness <- src_df_fitness %>% filter(Nham_aa == 1)
# #dim(src_df_fitness) # 5310 16
# 
# src_df_fitness_cleaned <- src_df_fitness %>%
#   dplyr::select(variant_name, FL_kinase_fitness_scaled, FL_kinase_sigma_scaled,
#          FL_abundance_fitness_scaled, FL_abundance_sigma_scaled, KD_kinase_fitness_scaled,
#          KD_kinase_sigma_scaled, KD_abundance_fitness_scaled, KD_abundance_sigma_scaled) %>%
#   group_by(variant_name) %>%
#   summarise(across(everything(), mean, na.rm = TRUE)) %>%
#   dplyr::rename(id = variant_name) %>%
#   mutate (Pos_real = as.integer(str_extract(id, '\\d+')),
#             wt_aa = substr(id, 1, 1),)

# src_df_fitness_cleaned_unique <- src_df_fitness_cleaned %>%  
#   distinct(Pos_real, .keep_all = TRUE) %>%  
#   filter(!is.na(Pos_real)) %>%  arrange(Pos_real)
# sequence <- paste(src_df_fitness_cleaned_unique$wt_aa, collapse = "")
# print(sequence)

#write.csv(src_df_fitness_cleaned, './data/cleaned/src_fit_cleaned.csv', row.names = FALSE)

src_df_fitness_cleaned <- read_csv('./data/cleaned_ddg/src_fit_cleaned.csv')
src_final_df_fit <- process_and_merge_protein_data(src_data, src_df_fitness_cleaned)
src_final_df_fit_anno <- merge(src_final_df_fit, src_anno, by.x="Pos_real", by.y = "position")

# Initialize the 'active site' column with NA
src_final_df_fit_anno$active_site <- NA

# Update the 'active site' column for the specified positions to TRUE
src_final_df_fit_anno$active_site[src_final_df_fit_anno$Pos_real %in% active_site_pos] <- TRUE

# Fill the remaining with FALSE
src_final_df_fit_anno$active_site[is.na(src_final_df_fit_anno$active_site)] <- FALSE


########################################
######### SRC RSA Data #################
########################################
# Use the function to read your file
src_rsa_file <- "./data/sasa/2src_rsa.txt"
src_rsa_df <- read_sasa_file(src_rsa_file)
src_rsa_df <- src_rsa_df[src_rsa_df$Chain == "A",]











