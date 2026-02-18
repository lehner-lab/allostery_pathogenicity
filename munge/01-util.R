########## DATA PREPROCESSING ##########

#' FUNCTION TO PROCESS VEP SCORES
#' @param input extracted vep score files
#' @param eve_file only needed to get separate eve scores for KRAS so this is conditional
#'  
#' @return merged data with all vep scores for one protein
process_protein_data <- function(consurf_file, am_file, esm1b_file, popeve_file, eve_file = NULL) {
  # Load and process consurf file
  consurf_df <- read_delim(consurf_file, delim = "\t", trim_ws = TRUE) %>%
    dplyr::rename(position = POS, wt_aa = SEQ, atom = ATOM, score = SCORE, color = COLOR, 
           confidence_interval = `CONFIDENCE INTERVAL`, msa_data = `MSA DATA`) %>%
    mutate(across(c(position, score), as.numeric))
  
  # Load and process AM file
  epsilon <- 1e-10  # Small value to avoid Inf
  am_df <- read_delim(am_file, delim = '\t', col_names = FALSE) %>%
    setNames(c('uniprot', 'id', 'score', 'patho')) %>%
    mutate(
      position = as.integer(str_extract(id, '\\d+')),
      wt_aa = substr(id, 1, 1),
      score_clamped = pmin(pmax(score, epsilon), 1 - epsilon),
      logits = log(score_clamped / (1 - score_clamped))
    )
  
  # Load ESM1b file
  esm1b_df <- read_csv(esm1b_file) %>%
    setNames(c('id', 'score', 'position'))
  
  # Load and process popEVE data
  popeve_df <- read_csv(popeve_file) %>%
    dplyr::select(id = mutant, popEVE, ESM1v) %>%
    as_tibble()
  
  popeve_selected <- popeve_df[, c("id", "popEVE")]
  esm1v_df <- popeve_df[, c("id", "ESM1v")]
  
  # Conditionally load EVE data either from eve_file or from popeve data
  if (!is.null(eve_file)) {
    # Load and process the provided EVE file
    eve_df <- read_csv(eve_file) %>%
      dplyr::select(wt_aa, position, mt_aa, EVE_scores_ASM)  %>% 
      mutate(id = paste0(wt_aa, position, mt_aa))
  } else {
    # Create EVE data from the popEVE data (assuming it's constructed similarly)
    eve_df <- read_csv(popeve_file) %>%
      dplyr::select(id=mutant, EVE) %>%
      as_tibble()
  }
  
  # Return the processed data
  return(list(
    consurf = consurf_df,
    am = am_df,
    esm1b = esm1b_df,
    popeve = popeve_selected,
    esm1v = esm1v_df,
    eve = eve_df
  ))
}  

#' FUNCTION TO MERGE VEP SCORES WITH DDG / FITNESS VALUES
#' @param input merged vep data from above
#' @param exp_data cleaned ddg or fitness df
#'  
#' @return merged data with all vep scores and exp mutational effect measurements for one protein
process_and_merge_protein_data <- function(protein_data, merged_df) {
  
  # Merge with consurf data, dropping unnecessary columns
  merged_df_consurf_ddg <- merged_df %>%
    merge(protein_data$consurf, by.x = "Pos_real", by.y = "position", all.x = TRUE) %>%
    dplyr::rename(consurf = score) 
  #select(-wt_aa, -atom, -confidence_interval, -msa_data, -`RESIDUE VARIETY`)
  
  # Merge with AM data and rename columns
  merged_df_ddg <- merged_df_consurf_ddg %>%
    merge(protein_data$am %>% dplyr::select(id, score, logits), by = "id") %>%
    dplyr::rename(AM = score, AM_logit = logits)
  
  # Merge with ESM1b data
  merged_df_ddg <- merged_df_ddg %>%
    merge(protein_data$esm1b %>% dplyr::select(id, score), by = "id", all.x = FALSE) %>%
    dplyr::rename(ESM1b = score)
  
  # Merge with popEVE and ESM1v data
  merged_df_ddg <- merged_df_ddg %>%
    merge(protein_data$popeve, by = "id", all.x = FALSE) %>%
    merge(protein_data$esm1v, by = "id", all.x = FALSE)
  
  # Merge with EVE data
  if ("EVE_scores_ASM" %in% colnames(protein_data$eve)) {
    # Dedicated EVE file (e.g. KRAS) has raw ASM scores — rename to EVE
    merged_df_ddg <- merged_df_ddg %>%
      merge(protein_data$eve, by = "id", all.x = FALSE) %>%
      dplyr::rename(EVE = EVE_scores_ASM)
  } else {
    # EVE column extracted from popEVE file (SRC, PDZ3, SH3)
    merged_df_ddg <- merged_df_ddg %>%
      merge(protein_data$eve, by = "id", all.x = FALSE)
  }
  
  return(merged_df_ddg)
}


#' Parse FreeSASA RSA output into a tidy data frame.
#'
#' Reads the fixed-width RSA output from FreeSASA (lines starting with "RES"),
#' converts 3-letter residue codes to 1-letter, and fills gaps in the residue
#' numbering with NA rows so the result aligns with experimental data.
#'
#' @param file Path to FreeSASA .rsa output file.
#' @return Data frame with columns: Residue, Chain, Num, All_atoms_ABS, All_atoms_REL.
read_sasa_file <- function(file) {
  # Uses AA_CODES from lib/globals.R
  aa_codes <- AA_CODES
  
  # Read the file
  lines <- readLines(file)
  
  # Filter lines that contain residue information (start with 'RES ')
  residue_lines <- grep("^RES ", lines, value = TRUE)
  
  # Create an empty data frame to store the results
  sasa_df <- data.frame(
    Residue = character(), Chain = character(), Num = integer(),
    All_atoms_ABS = numeric(), All_atoms_REL = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Iterate over each residue line and extract the data
  for (line in residue_lines) {
    fields <- strsplit(line, "\\s+")[[1]]
    
    # Extract residue data
    res_name <- fields[2]
    chain <- fields[3]
    num <- as.integer(fields[4])
    all_atoms_abs <- as.numeric(fields[5])
    all_atoms_rel <- as.numeric(fields[6])
    
    # Translate the three-letter residue code to one-letter
    aa_code <- ifelse(res_name %in% names(aa_codes), aa_codes[[res_name]], NA)
    
    # Store the extracted data into the data frame
    sasa_df <- sasa_df %>% 
      add_row(
        Residue = aa_code, Chain = chain, Num = num,
        All_atoms_ABS = all_atoms_abs, All_atoms_REL = all_atoms_rel
      )
  }
  
  # Find the minimum and maximum residue numbers
  min_num <- min(sasa_df$Num, na.rm = TRUE)
  max_num <- max(sasa_df$Num, na.rm = TRUE)
  
  # Create a complete sequence of residue numbers
  full_res_nums <- seq(min_num, max_num)
  
  # Merge the complete sequence with the original data frame
  # This will insert NA rows for any missing residues
  sasa_df_full <- full_res_nums %>%
    as.data.frame() %>%
    dplyr::rename(Num = 1) %>%
    left_join(sasa_df, by = "Num") %>%
    mutate(Chain = coalesce(Chain, dplyr::first(sasa_df$Chain))) # Keep chain consistent
  
  return(sasa_df_full)
}

#' Compare WT and mutant protein sequences to identify substitutions.
#'
#' Returns variant IDs in the format "WtAA{position}MutAA" (e.g. "A12V").
#'
#' @param wt     Wild-type protein sequence (single string).
#' @param mut    Mutant protein sequence (single string).
#' @param offset Integer added to 1-based sequence index to obtain the
#'               reference position numbering. Use offset=0 (default) for
#'               most proteins; use offset=1 for KRAS where the DMS data
#'               is 0-indexed relative to the sequence.
#' @return Comma-separated string of variant IDs, or "" if identical.
compare_sequences <- function(wt, mut, offset = 0) {
  differences <- c()
  for (i in 1:min(nchar(wt), nchar(mut))) {
    wt_aa <- substr(wt, i, i)
    mut_aa <- substr(mut, i, i)
    if (wt_aa != mut_aa) {
      differences <- c(differences, paste0(wt_aa, i + offset, mut_aa))
    }
  }
  if (length(differences) == 0) "" else paste(differences, collapse = ", ")
}

# Backwards-compatible wrapper for KRAS (offset = 1)
compare_sequences_kras <- function(wt, mut) {
  compare_sequences(wt, mut, offset = 1)
}

######################### ALIGNMENT CHECK ##################################
#' FUNCTION TO EXTRACT PROTEIN SEQUENCE FROM FILES
#' @param input datafile
#'  
#' @return protein sequence by taking the unique ids
ali_process_sequence <- function(data, position_col, wt_aa_col, id_col = NULL) {
  if (!is.null(id_col)) {
    data <- data %>%
      mutate(position = as.integer(str_extract(!!sym(id_col), '\\d+')),
             wt_aa = substr(!!sym(id_col), 1, 1))
  }
  
  data %>%
    distinct(!!sym(position_col), .keep_all = TRUE) %>%
    filter(!is.na(!!sym(position_col))) %>%
    arrange(as.numeric(!!sym(position_col))) %>%  # Sort by numeric position
    summarise(sequence = paste(!!sym(wt_aa_col), collapse = "")) %>%
    pull(sequence)
}

#' FUNCTION TO FIND INDEX MATCH-UP
#' @param input am_seq, wt_exp_seq
#'  
#' @return index the shorter sequence is found
ali_find_sequence_overlap <- function(full_sequence, sub_sequence, sequence_name) {
  start_index <- str_locate(full_sequence, sub_sequence)
  if (!is.na(start_index[1, 1])) {
    print(paste(sequence_name, "starts at index", start_index[1, 1], "in the reference sequence."))
  } else {
    print(paste("No overlap found for", sequence_name, "in the reference sequence."))
  }
}
