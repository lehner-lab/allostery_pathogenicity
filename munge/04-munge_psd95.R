########## DATA PREPROCESSING ##########

# ddg_df <- read.csv("./data/paper_supplements/psd95_mochi_refit.csv")
# ddg_df <- ddg_df %>% filter(library %in% c("761_abundance", "761_808"))
# dim(ddg_df)
# 
# pdz3_data <- process_protein_data(
#   consurf_file = './data/scores/3_consurf/psd95_consurf/new_cleaned_file.txt',
#   am_file = './data/scores/0_am/alphamissense_psd95.tsv',
#   esm1b_file = './data/scores/1_esm1b/DLG4_P78352_esm1b.csv',
#   popeve_file = './data/scores/2_popeve/DLG4_NP_001356.1.csv'
# )
# 
# #pdz_ddg_unique <- ali_process_sequence(ddg_df, "Pos", "WT_aa.1")
# #pdz_esm1v_unique <- ali_process_sequence(pdz3_data$esm1v, "position", "wt_aa", "id")
# #ali_find_sequence_overlap(pdz_esm1v_unique, pdz_ddg_unique ,"ESM1v sequence")
# #ESM1v sequence starts at index 354 in the reference sequence.
# 
# pdz3_ddg <- ddg_df %>%
#   mutate(wt_aa = str_extract(id, '^[A-Za-z]+'),
#          mt_aa = str_extract(id, '[A-Za-z]+$')) %>%
#   mutate(pos_eve = as.numeric(Pos) + 353) %>%
#   rename(pos_am = Pos)  %>%
#   rename(id_old = id)  %>%
#   mutate(id_eve = paste0(wt_aa, pos_eve, mt_aa))
# head(pdz3_ddg)
# 
# process_and_merge_protein_data_pdz <- function(protein_data, merged_df) {
#   
#   # Merge with popEVE and ESM1v data
#   merged_df_ddg <- merged_df %>%
#     merge(protein_data$popeve, by.x = "id_eve", by.y = "id", all.x = FALSE) %>%
#     merge(protein_data$esm1v, by.x = "id_eve", by.y = "id", all.x = FALSE)
# 
#   return(merged_df_ddg)
# }
# 
# pdz3_final_df_ddg <-process_and_merge_protein_data_pdz(pdz3_data, pdz3_ddg)
# dim(pdz3_final_df_ddg) # 3154   41
# 
# write.csv(pdz3_final_df_ddg, './data/cleaned_ddg/pdz3_ddg_cleaned_mochi_refit.csv', row.names = FALSE)

pdz3_final_df_ddg <- read_csv('./data/cleaned_ddg/pdz3_ddg_cleaned_mochi_refit.csv')
#head(pdz3_final_df_ddg)
#colnames(pdz3_final_df_ddg)

########################################
######### PDZ3 RSA Data #################
########################################
# Use the function to read your file
# pdz3_rsa_file <- "./data/sasa/1be9_rsa.txt"
# pdz3_rsa_df <- read_sasa_file(pdz3_rsa_file)
# pdz3_rsa_df <- pdz3_rsa_df[pdz3_rsa_df$Chain == "A",]

##### ARCHIVED
########################################
### PSD95-PDZ3 DDG check alignment #####
########################################
# domain_ddg <- read_excel("./data/paper_supplements/domains_andre/41586_2022_4586_MOESM10_ESM.xlsx", sheet = 2)
# 
# pdz3_ddg <- domain_ddg %>% filter(protein == "PSD95-PDZ3") %>% filter(mut_order==1) 
# 
# pdz3_ddg_prep <- pdz3_ddg %>%
#   mutate (wt_aa = substr(id, 1, 1),)
# 
# # Process each dataset and create sequences
# pdz_ddg_unique <- ali_process_sequence(pdz3_ddg_prep, "Pos_ref", "wt_aa")
# pdz_am_unique <- ali_process_sequence(pdz3_data$am, "position", "wt_aa")
# pdz_esm1v_unique <- ali_process_sequence(pdz3_data$esm1v, "position", "wt_aa", "id")
# pdz_esm1b_unique <- ali_process_sequence(pdz3_data$esm1b, "position", "wt_aa", "id")
# pdz_popeve_unique <- ali_process_sequence(pdz3_data$popeve, "position", "wt_aa", "id")
# pdz_eve_unique <- ali_process_sequence(pdz3_data$eve, "position", "wt_aa", "id")
# 
# 
# # Find overlaps for each sequence
# ali_find_sequence_overlap(pdz_am_unique, pdz_ddg_unique ,"AM sequence")
# #"AM sequence starts at index 311 in the reference sequence."
# ali_find_sequence_overlap(pdz_esm1v_unique, pdz_ddg_unique,"ESM1v sequence")
# #"ESM1v sequence starts at index 354 in the reference sequence."
# ali_find_sequence_overlap(pdz_esm1b_unique, pdz_ddg_unique,"ESM1b sequence")
# #"ESM1b sequence starts at index 311 in the reference sequence."
# ali_find_sequence_overlap(pdz_popeve_unique, pdz_ddg_unique,"popEVE sequence")
# #"popEVE sequence starts at index 354 in the reference sequence."
# ali_find_sequence_overlap(pdz_eve_unique, pdz_ddg_unique,"EVE sequence")
# #"EVE sequence starts at index 354 in the reference sequence."
# #

# ########################################
# # PSD95-PDZ3 fitness check alignment ###
# ########################################
# domain_fit <- read_excel('./data/paper_supplements/domains_andre/41586_2022_4586_MOESM9_ESM.xlsx', sheet = "TableS6")
# pdz3_fit <- domain_fit %>% filter(protein == "PSD95-PDZ3")
# pdz3_wt <- pdz3_fit %>% filter(WT == TRUE) %>% pull(aa_seq) %>% .[1]
# #dim(pdz3_fit) #15229     9
# 
# pdz3_fit <- pdz3_fit %>% filter(Nham_aa == 1)
# pdz3_fit <- pdz3_fit %>%
#   mutate(id = sapply(aa_seq, function(x) compare_sequences(pdz3_wt, x))) 
# 
# pdz3_fit_prep <- pdz3_fit %>%
#   mutate (Pos_real = as.integer(str_extract(id, '\\d+')),
#           wt_aa = substr(id, 1, 1),) 
# 
# pdz_fit_unique <- ali_process_sequence(pdz3_fit_prep, "Pos_real", "wt_aa")
# pdz_fit_unique
# 
# # Find overlaps for each sequence
# ali_find_sequence_overlap(pdz_am_unique, pdz_fit_unique ,"AM sequence")
# #"AM sequence starts at index 311 in the reference sequence."
# ali_find_sequence_overlap(pdz_esm1v_unique, pdz_fit_unique,"ESM1v sequence")
# #"ESM1v sequence starts at index 354 in the reference sequence."
# ali_find_sequence_overlap(pdz_esm1b_unique, pdz_fit_unique,"ESM1b sequence")
# #"ESM1b sequence starts at index 311 in the reference sequence."
# ali_find_sequence_overlap(pdz_popeve_unique, pdz_fit_unique,"popEVE sequence")
# #"popEVE sequence starts at index 354 in the reference sequence."
# ali_find_sequence_overlap(pdz_eve_unique, pdz_fit_unique,"EVE sequence")
# #"EVE sequence starts at index 354 in the reference sequence."

########## DATA PREPROCESSING ##########

########################################
######### PSD95-PDZ3 DDG Data ##########
########################################

# Read in ddg
#domain_ddg <- read_excel("./data/paper_supplements/domains_andre/41586_2022_4586_MOESM10_ESM.xlsx", sheet = 2)
# 
#pdz3_ddg <- domain_ddg %>% filter(protein == "PSD95-PDZ3") %>% filter(mut_order==1)
#unique(pdz3_ddg[pdz3_ddg$orthosteric==TRUE,]$Pos_ref)
# pdz3_ddg <- pdz3_ddg %>%
#   mutate(wt_aa = str_extract(id, '^[A-Za-z]+'),
#          mt_aa = str_extract(id, '[A-Za-z]+$')) %>%
#   rename(pos_am = Pos_ref)  %>%
#   rename(id_old = id)  %>%
#   rename(id_am = id_ref)  %>%
#   mutate(pos_eve = as.numeric(Pos) + 353) %>%
#   mutate(id_eve = paste0(wt_aa, pos_eve, mt_aa))
# 
# pdz3_data <- process_protein_data(
#   consurf_file = './data/scores/3_consurf/psd95_consurf/new_cleaned_file.txt',
#   am_file = './data/scores/0_am/alphamissense_psd95.tsv',
#   esm1b_file = './data/scores/1_esm1b/DLG4_P78352_esm1b.csv',
#   popeve_file = './data/scores/2_popeve/DLG4_NP_001356.1.csv'
# )

# process_and_merge_protein_data_pdz <- function(protein_data, merged_df) {
#   
#   # Merge with consurf data, dropping unnecessary columns
#   merged_df_consurf_ddg <- merged_df %>%
#     merge(protein_data$consurf, by.x = "pos_am", by.y = "position", all.x = TRUE) %>%
#     dplyr::rename(consurf = score) 
#   #select(-wt_aa, -atom, -confidence_interval, -msa_data, -`RESIDUE VARIETY`)
#   
#   # Merge with AM data and rename columns
#   merged_df_ddg <- merged_df_consurf_ddg %>%
#     merge(protein_data$am %>% dplyr::select(id, score, logits), by.x = "id_am", by.y = "id") %>%
#     dplyr::rename(AM = score, AM_logit = logits)
#   
#   # Merge with ESM1b data
#   merged_df_ddg <- merged_df_ddg %>%
#     merge(protein_data$esm1b %>% dplyr::select(id, score),by.x = "id_am", by.y = "id", all.x = FALSE) %>%
#     dplyr::rename(ESM1b = score)
#   
#   # Merge with popEVE and ESM1v data
#   merged_df_ddg <- merged_df_ddg %>%
#     merge(protein_data$popeve, by.x = "id_eve", by.y = "id", all.x = FALSE) %>%
#     merge(protein_data$esm1v, by.x = "id_eve", by.y = "id", all.x = FALSE)
#   
#   merged_df_ddg <- merged_df_ddg %>%
#     merge(protein_data$eve,  by.x = "id_eve", by.y = "id", all.x = FALSE)  
#   
#   return(merged_df_ddg)
# }

# pdz3_final_df_ddg <-process_and_merge_protein_data_pdz(pdz3_data, pdz3_ddg)
# dim(pdz3_final_df_ddg) # 1587 43
# 
# write.csv(pdz3_final_df_ddg, './data/cleaned/pdz3_ddg_cleaned.csv', row.names = FALSE)

#pdz3_final_df_ddg <- read_csv('./data/cleaned_ddg/pdz3_ddg_cleaned.csv')


# Check for NA values in each column
#colSums(is.na(pdz3_final_df))
#dim(pdz3_final_df) # 95 43

########################################
##### PSD95-PDZ3 fitness Data ##########
########################################

#domain_fit <- read_excel('./data/paper_supplements/domains_andre/41586_2022_4586_MOESM9_ESM.xlsx', sheet = "TableS6")
# pdz3_fit <- domain_fit %>% filter(protein == "PSD95-PDZ3")
# pdz3_wt <- pdz3_fit %>% filter(WT == TRUE) %>% pull(aa_seq) %>% .[1]
# #dim(pdz3_fit) #15229     9
# 
# pdz3_fit <- pdz3_fit %>% filter(Nham_aa == 1)
# pdz3_fit <- pdz3_fit %>%
#   mutate(id = sapply(aa_seq, function(x) compare_sequences(pdz3_wt, x))) %>%
#   mutate (Pos_real = as.integer(str_extract(id, '\\d+')))
# 
# #dim(pdz3_fit) #2721 9
# # colSums(is.na(pdz3_fit)) # no NA!
# 
# pdz3_fit <- pdz3_fit %>%
#   dplyr::select(id, pca_type, Pos_real, fitness, sigma, growthrate, growthrate_sigma) %>%
#   mutate(wt_aa = str_extract(id, '^[A-Za-z]+'),
#          mt_aa = str_extract(id, '[A-Za-z]+$')) %>%
#   rename(Pos = Pos_real) %>%
#   rename(id_old = id)  %>%
#   mutate(pos_eve = as.numeric(Pos) + 353) %>%  #"EVE/ESM1v/popEVE sequence starts at index 354 in the reference sequence."
#   mutate(id_eve = paste0(wt_aa, pos_eve, mt_aa)) %>%
#   mutate(pos_am = as.numeric(Pos) + 310) %>%  #"AM/ESM1b sequence starts at index 311 in the reference sequence."
#   mutate(id_am = paste0(wt_aa, pos_am, mt_aa))
# 
# dim(pdz3_fit) #2721   12
# 
# pdz3_final_fit <- process_and_merge_protein_data_pdz(pdz3_data, pdz3_fit)
# dim(pdz3_final_fit) # 2721 26
# # 
# write.csv(pdz3_final_fit, './data/cleaned/pdz3_fit_cleaned.csv', row.names = FALSE)










