########## DATA PREPROCESSING ##########

########################################
### GRB2-SH3 DDG check alignment #####
########################################
# sh3_ddg <- domain_ddg %>% filter(protein == "GRB2-SH3") %>% filter(mut_order==1) 
# 
# sh3_ddg_prep <- sh3_ddg %>%
#   mutate (wt_aa = substr(id, 1, 1),)
# 
# # Process each dataset and create sequences
# sh3_ddg_unique <- ali_process_sequence(sh3_ddg_prep, "Pos_ref", "wt_aa")
# sh3_am_unique <- ali_process_sequence(sh3_data$am, "position", "wt_aa")
# sh3_esm1v_unique <- ali_process_sequence(sh3_data$esm1v, "position", "wt_aa", "id")
# sh3_esm1b_unique <- ali_process_sequence(sh3_data$esm1b, "position", "wt_aa", "id")
# sh3_popeve_unique <- ali_process_sequence(sh3_data$popeve, "position", "wt_aa", "id")
# sh3_eve_unique <- ali_process_sequence(sh3_data$eve, "position", "wt_aa", "id")
# 
# # Find overlaps for each sequence
# ali_find_sequence_overlap(sh3_am_unique, sh3_ddg_unique ,"AM sequence")
# #"AM sequence starts at index 159 in the reference sequence."
# ali_find_sequence_overlap(sh3_esm1v_unique, sh3_ddg_unique,"ESM1v sequence")
# #"ESM1v sequence starts at index 159 in the reference sequence."
# ali_find_sequence_overlap(sh3_esm1b_unique, sh3_ddg_unique,"ESM1b sequence")
# #"ESM1b sequence starts at index 159 in the reference sequence."
# ali_find_sequence_overlap(sh3_popeve_unique, sh3_ddg_unique,"popEVE sequence")
# #"popEVE sequence starts at index 159 in the reference sequence."
# ali_find_sequence_overlap(sh3_eve_unique, sh3_ddg_unique,"EVE sequence")
# #"EVE sequence starts at index 159 in the reference sequence."

# ########################################
# # GRB2-SH3 fitness check alignment ###
# ########################################
# sh3_fit <- domain_fit %>% filter(protein == "GRB2-SH3")
# sh3_wt <- sh3_fit %>% filter(WT == TRUE) %>% pull(aa_seq) %>% .[1]
# dim(sh3_fit) #96809     9
# 
# sh3_fit <- sh3_fit %>% filter(Nham_aa == 1)
# sh3_fit <- sh3_fit %>%
#   mutate(id = sapply(aa_seq, function(x) compare_sequences(sh3_wt, x))) 
# 
# sh3_fit_prep <- sh3_fit %>%
#   mutate (Pos_real = as.integer(str_extract(id, '\\d+')),
#           wt_aa = substr(id, 1, 1),) 
# 
# sh3_fit_unique <- ali_process_sequence(sh3_fit_prep, "Pos_real", "wt_aa")
# sh3_fit_unique
# 
# # Find overlaps for each sequence
# ali_find_sequence_overlap(sh3_am_unique, sh3_fit_unique ,"AM sequence")
# #"AM sequence starts at index 159 in the reference sequence."
# ali_find_sequence_overlap(sh3_esm1v_unique, sh3_fit_unique,"ESM1v sequence")
# #"ESM1v sequence starts at index 159 in the reference sequence."
# ali_find_sequence_overlap(sh3_esm1b_unique, sh3_fit_unique,"ESM1b sequence")
# #"ESM1b sequence starts at index 159 in the reference sequence."
# ali_find_sequence_overlap(sh3_popeve_unique, sh3_fit_unique,"popEVE sequence")
# #"popEVE sequence starts at index 159 in the reference sequence."
# ali_find_sequence_overlap(sh3_eve_unique, sh3_fit_unique,"EVE sequence")
# #"EVE sequence starts at index 159 in the reference sequence."

### PSD95-sh3 Data

# Read in ddg
# domain_ddg <- read_excel("./data/paper_supplements/domains_andre/41586_2022_4586_MOESM10_ESM.xlsx", sheet = 2)
# 
# sh3_ddg <- domain_ddg %>% filter(protein == "GRB2-SH3") %>% filter(mut_order==1) 
# 
# sh3_ddg <- sh3_ddg %>%
#   mutate(wt_aa = str_extract(id, '^[A-Za-z]+'),
#          mt_aa = str_extract(id, '[A-Za-z]+$')) %>%
#   rename(old_id = id) %>% 
#   mutate(Pos_real = as.numeric(Pos) + 158) %>% 
#   mutate(id = paste0(wt_aa, Pos_real, mt_aa))
# 
# sh3_ddg
# 
# sh3_data <- process_protein_data(
#   consurf_file = './data/scores/3_consurf/grb2_consurf/new_cleaned_file.txt',
#   am_file = './data/scores/0_am/alphamissense_grb2.tsv',
#   esm1b_file = './data/scores/1_esm1b/GRB2_P62993_esm1b.csv',
#   popeve_file = './data/scores/2_popeve/GRB2_NP_002077.1.csv'
# )
# 
# 
# sh3_final_df <-process_and_merge_protein_data(sh3_data, sh3_ddg)
# 
# write.csv(sh3_final_df, './data/cleaned/sh3_ddg_cleaned.csv', row.names = FALSE)

sh3_final_df_ddg <- read_csv('./data/cleaned_ddg/sh3_ddg_cleaned.csv')


# Check for NA values in each column
#colSums(is.na(sh3_final_df))

#dim(sh3_final_df) # 1056   43

########################################
##### GRB2-SH3 fitness Data ##########
########################################
# sh3_fit <- domain_fit %>% filter(protein == "GRB2-SH3")
# sh3_wt <- sh3_fit %>% filter(WT == TRUE) %>% pull(aa_seq) %>% .[1]
# dim(sh3_fit) #96809     9
# 
# sh3_fit <- sh3_fit %>% filter(Nham_aa == 1)
# sh3_fit <- sh3_fit %>%
#   mutate(id = sapply(aa_seq, function(x) compare_sequences(sh3_wt, x))) %>%
#   mutate (Pos_real = as.integer(str_extract(id, '\\d+')))
# 
# sh3_fit <- sh3_fit %>%
#   mutate(wt_aa = str_extract(id, '^[A-Za-z]+'),
#          mt_aa = str_extract(id, '[A-Za-z]+$')) %>%
#   rename(old_id = id) %>% 
#   mutate(Pos_real = as.numeric(Pos_real) + 158) %>%  #"EVE/AM/ALL sequence starts at index 159 in the reference sequence."
#   mutate(id = paste0(wt_aa, Pos_real, mt_aa))
# 
# sh3_final_fit <- process_and_merge_protein_data(sh3_data, sh3_fit)
# dim(sh3_final_fit) # 1689   27
# 
# write.csv(sh3_final_fit, './data/cleaned/sh3_fit_cleaned.csv', row.names = FALSE)

sh3_final_df_fit <- read_csv('./data/cleaned_ddg/sh3_fit_cleaned.csv')

########################################
######### SH3 RSA Data #################
########################################
# Use the function to read your file
sh3_rsa_file <- "./data/sasa/2vwf_rsa.txt"
sh3_rsa_df <- read_sasa_file(sh3_rsa_file)
sh3_rsa_df <- sh3_rsa_df[sh3_rsa_df$Chain == "A",]


