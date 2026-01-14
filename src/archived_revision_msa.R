


```{r}
library(Biostrings)
library(tidyverse)

# 1. Read the MSA
msa <- readAAStringSet("~/Downloads/psd95-pdz3_phmmer_aligned.afa")
msa_mat <- as.matrix(msa)
n_seq <- nrow(msa_mat)
n_pos <- ncol(msa_mat)

# 2. Prepare Data for Fast Math
# Convert AA to integers (0 for gaps, 1-20 for residues)
num_mat <- matrix(as.integer(factor(msa_mat)), nrow = n_seq)
is_gap <- msa_mat == "-" | msa_mat == "."
num_mat[is_gap] <- 0
presence <- 1 - is_gap  # 1 where there is a residue

# 3. BLOCK-PROCESSING WEIGHT CALCULATION (Memory Efficient)
calc_neff_blocks <- function(num_mat, presence, is_gap, threshold = 0.8, block_size = 500) {
  n_seq <- nrow(num_mat)
  n_pos <- ncol(num_mat)
  m_i <- rep(0, n_seq)
  
  # Pre-calculate column indices for each residue to speed up matching
  res_indices <- lapply(1:20, function(r) which(num_mat == r, arr.ind = TRUE))
  
  num_blocks <- ceiling(n_seq / block_size)
  
  print(paste("Processing", n_seq, "sequences in", num_blocks, "blocks..."))
  
  for (b in 1:num_blocks) {
    start <- (b - 1) * block_size + 1
    end <- min(b * block_size, n_seq)
    idx_range <- start:end
    
    # Subset current block
    block_num <- num_mat[idx_range, , drop = FALSE]
    block_gap <- is_gap[idx_range, , drop = FALSE]
    
    # 1. Calculate denominator: non-gap positions in either sequence
    # Matrix: [block_size x n_seq]
    # Length = Total - positions where BOTH are gaps
    both_gaps <- block_gap %*% t(is_gap)
    non_gap_counts <- n_pos - both_gaps
    
    # 2. Calculate numerator: matches
    matches <- matrix(0, nrow = length(idx_range), ncol = n_seq)
    for (r in 1:20) {
      res_binary_block <- (block_num == r)
      res_binary_full <- (num_mat == r)
      matches <- matches + (res_binary_block %*% t(res_binary_full))
    }
    
    # 3. Count neighbors above threshold
    identity_mat <- matches / non_gap_counts
    m_i[idx_range] <- rowSums(identity_mat >= threshold, na.rm = TRUE)
    
    if (b %% 10 == 0) {
      print(paste("Completed block", b, "of", num_blocks, "-", round(100 * b / num_blocks, 1), "%"))
    }
    # Explicit garbage collection to keep RAM clean
    gc()
  }
  
  return(1 / m_i)
}

# Run the memory-safe version
# Adjust block_size (lower = less RAM, higher = faster)
sequence_weights <- calc_neff_blocks(num_mat, presence, is_gap, threshold = 0.8, block_size = 500)

# 4. Compute Positional Neff
print("Finalizing positional Neff calculation...")
# Use a loop for the final sum to avoid one last massive matrix multiplication
pos_neff <- numeric(n_pos)
for(j in 1:n_pos) {
  residue_mask <- !is_gap[, j]
  pos_neff[j] <- sum(sequence_weights[residue_mask])
}

# 5. Extract for Query
target_pattern <- "B9EGL1/351"
target_index <- grep(target_pattern, names(msa))
query_residue_indices <- which(!(msa_mat[target_index, ] %in% c("-", ".")))
final_neff <- pos_neff[query_residue_indices]

# 6. Result
msa_depth_df <- data.frame(
  res_num = (1:length(final_neff)) + 310,
  neff = final_neff
)

head(msa_depth_df)
```


```{r}
library(Biostrings)
library(tidyverse)

# 1. Read the aligned FASTA file
msa <- readAAStringSet("~/Downloads/psd95-pdz3_phmmer_aligned.afa")

# 2. Convert to matrix
msa_matrix <- as.matrix(msa)

# 3. Identify the Query Sequence
# found by querying for homo sapiens on web page, then align based on uniprot id 
target_pattern <- "B9EGL1/351"
target_index <- grep(target_pattern, names(msa))

if (length(target_index) == 0) stop("Query sequence not found!")
query_seq <- msa_matrix[target_index, ]

# 4. Define gaps and valid amino acids
gaps <- c("-", ".")
# Standard 20 amino acids (optional: restricts analysis to standard AA)
std_aa <- c("A","R","N","D","C","Q","E","G","H","I",
            "L","K","M","F","P","S","T","W","Y","V")

# 5. Compute Consensus Matrix (Counts per residue per column)
cm <- consensusMatrix(msa)

# 6. Calculate Metrics per Column

# --- Metric A: Raw Coverage (Depth) ---
# Sum of all characters that are NOT gaps
valid_row_indices <- !rownames(cm) %in% gaps
raw_coverage <- colSums(cm[valid_row_indices, ])

# --- Metric B: Neff (Shannon Entropy based) ---
# 1. Subset matrix to only standard Amino Acids (removes gaps and ambiguous chars like X, B, Z)
cm_aa <- cm[rownames(cm) %in% std_aa, ]

# 2. Calculate frequencies (p) for each column
# Add a tiny pseudocount (e.g., 1e-6) to avoid log(0) if desired, 
# or just ignore 0s as 0*log(0) = 0 in entropy limits.
col_sums_aa <- colSums(cm_aa)
# Avoid division by zero for gap-only columns
col_sums_aa[col_sums_aa == 0] <- 1 
probs <- t(t(cm_aa) / col_sums_aa)

# 3. Calculate Shannon Entropy (H = -sum(p * log2(p)))
# We treat 0 * log(0) as 0
log_probs <- log2(probs)
log_probs[is.infinite(log_probs)] <- 0
shannon_entropy <- -colSums(probs * log_probs)

# 4. Calculate Neff (Perplexity = 2^Entropy)
# This represents the "Effective number of amino acids" at that position
neff_aa <- 2^shannon_entropy

# 7. FILTER: Keep only columns where YOUR QUERY has a residue
valid_residue_indices <- which(!query_seq %in% gaps)

# Extract values for the query positions
final_raw_depth <- raw_coverage[valid_residue_indices]
final_neff <- neff_aa[valid_residue_indices]

# 8. Verification Step
real_query_seq <- as.character(msa[[target_index]])
real_query_seq <- gsub("[-.]", "", real_query_seq) 
expected_length <- nchar(real_query_seq)

if(length(final_raw_depth) != expected_length) {
  warning(paste("Mismatch! Depth length:", length(final_raw_depth), 
                "vs Query length:", expected_length))
} else {
  print(paste("Success! Dimensions match query length of:", expected_length))
}

# 9. Create final dataframe
msa_depth_df <- data.frame(
  position = 1:length(final_raw_depth),
  raw_depth = final_raw_depth,
  neff_aa = final_neff
)

# Adjust position to start at 159 (158 + 1)
msa_depth_df$position <- msa_depth_df$position + 310

# Preview
head(msa_depth_df)
```


```{r}
pdz3_final_df_ddg <- read_csv('/Users/xl7/Documents/0.Projects/01.protein-seq-evo-v1/data/cleaned_ddg/pdz3_ddg_cleaned_mochi_refit.csv')
pdz3_final_df_ddgb <- pdz3_final_df_ddg %>% filter(trait_name == "Binding") %>% 
  dplyr::rename(b_ddg_pred = ddg,
                b_ddg_pred_sd = std_ddg)

pdz3_final_df_ddgf <- pdz3_final_df_ddg %>% filter(trait_name == "Folding") %>% 
  dplyr::rename(f_ddg_pred = ddg,
                f_ddg_pred_sd = std_ddg)

nrow(pdz3_final_df_ddgb) #1567
nrow(pdz3_final_df_ddgf)
```

```{r}
test_ddg <- fread("/Users/xl7/Documents/0.Projects/01.protein-seq-evo-v1/data/decay_pdb/PSD95/ThermoMPNN_inference_AF-P78352-F1-model_v4.csv")
test_esm <- fread("/Users/xl7/Documents/0.Projects/00.large_supplements/ESM1v_proteome_wide/all_ESM1v/P78352.csv")

colnames(test_esm)[2] <- "ESM1v"

test_ddg[, new_position := position + 1]
test_ddg[, variant := paste0(wildtype, new_position, mutation)]

test_df <- merge(test_ddg, test_esm, by = "variant")

test_df %>% filter(new_position == 311)
pdz3_final_df_ddgb %>% filter (pos_am == 1)
pdz3_final_df_ddgb$pos_match <- pdz3_final_df_ddgb$pos_am + 310
pdz3_final_df_ddgb$variant_match <- paste0(pdz3_final_df_ddgb$wt_aa, 
                                           pdz3_final_df_ddgb$pos_match, 
                                           pdz3_final_df_ddgb$mt_aa)

pdz3_final_df_ddgf$pos_match <- pdz3_final_df_ddgf$pos_am + 310
pdz3_final_df_ddgf$variant_match <- paste0(pdz3_final_df_ddgf$wt_aa, 
                                           pdz3_final_df_ddgf$pos_match, 
                                           pdz3_final_df_ddgf$mt_aa)

pdz3_final_df_ddgf <- pdz3_final_df_ddgf %>% dplyr::select(variant_match, f_ddg_pred, f_ddg_pred_sd, Pos_ref)
pdz3_final_df_ddgb <- pdz3_final_df_ddgb %>% dplyr::select(variant_match, b_ddg_pred, b_ddg_pred_sd)

test_merged_df <- merge(pdz3_final_df_ddgb, test_df, by.x="variant_match", by.y="variant")
test_merged_df <- merge(pdz3_final_df_ddgf, test_merged_df, by="variant_match")
nrow(test_merged_df) #1567
head(test_merged_df)
```
```{r}
paper_anno <- read.csv("/Users/xl7/Documents/0.Projects/01.protein-seq-evo-v1/data/cleaned_ddg/pdz3_ddg_cleaned.csv")
head(paper_anno)
table(paper_anno$orthosteric)
ortho_list <- unique(paper_anno %>% filter(orthosteric == TRUE) %>% pull(pos_am))
##setdiff(unique(pdz3_out %>% filter(scHAmin_ligand < 4.3) %>% pull(Pos)), unique(merged_df_pdz3 %>% filter(orthosteric == TRUE) %>% pull (Pos_ref.x)))
ortho_list <- append(ortho_list, c(322, 323, 324 ,325, 326, 327, 328, 331, 339, 372, 376, 379, 380))

fil_test_merged_df <- test_merged_df %>%
  filter(!Pos_ref %in% ortho_list)

# Fit a loess model using the filtered data
loess_fit_comp <- loess(ESM1v ~ ddG_pred, data = fil_test_merged_df, span = 0.7, family = "symmetric")

# Predict fitted values for ALL data points using the loess model trained on fil_gck_df
test_merged_df$fitted_comp <- predict(loess_fit_comp, newdata = test_merged_df)

# Calculate residuals for ALL points
test_merged_df$residuals_comp <- test_merged_df$ESM1v - test_merged_df$fitted_comp
sum(is.na(test_merged_df$residuals_comp)) #6
head(test_merged_df)
```

```{r, fig.width=4, fig.height=4.5}
position_residuals <- test_merged_df %>%
  # Ensure you have the residuals calculated. 
  # If not, fit the model here:
  # mutate(residual = resid(loess(ESM1v ~ f_ddg_pred, span=0.7))) %>%
  group_by(position) %>%
  summarize(
    # We take the median magnitude of residuals at this site
    median_residual = median(residuals_comp, na.rm = TRUE),
    median_ddgb = median(b_ddg_pred, na.rm = TRUE)
  )

# 2. Join with MSA Depth
plot_data <- position_residuals %>%
  inner_join(msa_depth_df, by = "position")

# 3. Correlation Test
# This is the number you will quote in your rebuttal.
cor_test <- cor.test(plot_data$neff_aa, plot_data$median_residual, method = "spearman")
print(cor_test)

# 4. Create the Scatter Plot
positions_to_highlight <- c(312, 314, 318, 329, 330, 336, 338, 341, 347, 357, 359, 360, 362, 363, 367, 375, 386, 388) 

# 2. Prepare data
plot_data <- plot_data %>%
  mutate(
    # Label logic: ONLY label if the position is in your list
    label_text = ifelse(position %in% positions_to_highlight, as.character(position), ""),
    
    # Optional: Create a column to make these points look different (e.g., bigger/bolder)
    is_highlighted = position %in% positions_to_highlight
  )

# 4. Create the Scatter Plot
p1 <- ggplot(plot_data, aes(x = neff_aa, y = median_residual,  color = median_ddgb)) +
  geom_point(aes(color = median_ddgb, alpha = 0.7)) +
  scale_color_gradient2(
    low = "blue", 
    mid = "grey", 
    high = "red", 
    midpoint = 0
  ) +
  theme_classic() +
  labs(
    x = "Neff",
    y = "Median ESM-1v-ddGf residual",
    title = "PSD95-PDZ3",
    subtitle = paste("Spearman's rho =", round(cor_test$estimate, 3), 
                     "| P-value =", format.pval(cor_test$p.value))
  ) +
  theme(
    text = element_text(size = 12),
    legend.position = "bottom"
  ) +
  geom_text_repel(
    aes(label = label_text),
    size = 4,                # Slightly larger text for emphasis
    fontface = "bold",       # Bold text
    box.padding = 0.5,
    max.overlaps = Inf,      # Force it to show all your specific labels
    point.padding = 0.3,
    min.segment.length = 0   # Always draw a line to the point
  ) 
p1
```

```{r}
# 1. Read the aligned FASTA file

msa <- readAAStringSet("~/Downloads/grb2-sh3_phmmer_aligned.afa")



# 2. Convert to matrix

msa_matrix <- as.matrix(msa)

# 3. Identify the Query Sequence

target_pattern <- "Q6ICN0/159"
target_index <- grep(target_pattern, names(msa))
query_seq <- msa_matrix[target_index, ]

# 4. Define gaps
gaps <- c("-", ".")

# 5. Calculate Raw Coverage (Total non-gap residues per column)

cm <- consensusMatrix(msa)
valid_amino_acids <- !rownames(cm) %in% gaps
raw_coverage <- colSums(cm[valid_amino_acids, ])



# 6. FILTER: Keep only columns where YOUR QUERY has a residue
valid_residue_indices <- which(!query_seq %in% gaps)
residue_depths <- raw_coverage[valid_residue_indices]

# 7. Verification Step (Optional but recommended)
# Get the actual sequence string without gaps to check length

real_query_seq <- as.character(msa[[target_index]])
real_query_seq <- gsub("[-.]", "", real_query_seq) # Remove gaps
expected_length <- nchar(real_query_seq)



if(length(residue_depths) != expected_length) {
  warning(paste("Mismatch! Depth length:", length(residue_depths), 
                "vs Query length:", expected_length))
} else {
  print(paste("Success! Dimensions match query length of:", expected_length))
}

# 8. Create final dataframe

msa_depth_df <- data.frame(
  position = 1:length(residue_depths),
  msa_depth = residue_depths
)


msa_depth_df$position <- msa_depth_df$position + 158
msa_depth_df
```
```{r}
sh3_final_df_ddgb <- read_tsv('/Users/xl7/Documents/0.Projects/01.protein-seq-evo-v1/data/paper_supplements/domains_andre/GRB2-SH3/Abundance/mochi_project/task_1/weights/weights_Binding.txt')
sh3_final_df_ddgb <- sh3_final_df_ddgb %>% dplyr::select(id_ref, Pos_ref, `mean_kcal/mol`, `std_kcal/mol`)
dim(sh3_final_df_ddgb) #757   4

sh3_final_df_ddgb <- sh3_final_df_ddgb %>%
  mutate(wt_aa = str_extract(id_ref, '^[A-Za-z]+'),
         mt_aa = str_extract(id_ref, '[A-Za-z]+$')) %>%
  dplyr::rename(old_id = id_ref) %>%
  dplyr::rename(old_Pos = Pos_ref) %>%
  dplyr::rename(b_ddg_pred = `mean_kcal/mol`) %>%
  dplyr::rename(b_ddg_pred_sd = `std_kcal/mol`) %>%
  mutate(Pos_ref = as.numeric(old_Pos) + 158) %>%
  mutate(id_ref = paste0(wt_aa, Pos_ref, mt_aa))

sh3_final_df_ddgb <- sh3_final_df_ddgb %>% filter(old_id != "WT")
head(sh3_final_df_ddgb)
```

```{r}
test_ddg <- fread("/Users/xl7/Documents/0.Projects/01.protein-seq-evo-v1/data/decay_pdb/GRB2/ThermoMPNN_inference_AF-P62993-F1-model_v4.csv")
test_esm <- fread("/Users/xl7/Documents/0.Projects/00.large_supplements/ESM1v_proteome_wide/all_ESM1v/P62993.csv")

colnames(test_esm)[2] <- "ESM1v"

test_ddg[, new_position := position + 1]
test_ddg[, variant := paste0(wildtype, new_position, mutation)]

test_df <- merge(test_ddg, test_esm, by = "variant")

pred_seq <-  ali_process_sequence(test_df, "new_position", "wildtype")
exp_seq <- ali_process_sequence(sh3_final_df_ddgb, "old_Pos", "wt_aa")
ali_find_sequence_overlap(pred_seq, exp_seq ,"sequence")
#"sequence starts at index 159 in the reference sequence."
nrow(sh3_final_df_ddgb) #756
test_merged_df <- merge(sh3_final_df_ddgb, test_df, by.x="id_ref", by.y="variant")
nrow(test_merged_df) #756
head(test_merged_df)
```

```{r}
sh3_final_df_ddgf <- sh3_final_df_ddgf %>% dplyr::select(id_ref, f_ddg_pred, f_ddg_pred_sd)
test_merged_df <- test_merged_df %>% inner_join(sh3_final_df_ddgf,by = "id_ref")
nrow(test_merged_df) #756
head(test_merged_df)
```

```{r}
paper_anno <- read.csv("/Users/xl7/Documents/0.Projects/01.protein-seq-evo-v1/data/cleaned_ddg/sh3_ddg_cleaned.csv")
table(paper_anno$orthosteric)
head(paper_anno)
ortho_list <- unique(paper_anno %>% filter(orthosteric == TRUE) %>% pull(Pos_real))

ortho_list<- append(ortho_list, c(170 ,192 ,195 ,204 ,209))
length(ortho_list)

fil_test_merged_df <- test_merged_df %>%
  filter(!Pos_ref %in% ortho_list)

# Fit a loess model using the filtered data
loess_fit_comp <- loess(ESM1v ~ ddG_pred, data = fil_test_merged_df, span = 0.7, family = "symmetric")

# Predict fitted values for ALL data points using the loess model trained on fil_gck_df
test_merged_df$fitted_comp <- predict(loess_fit_comp, newdata = test_merged_df)

# Calculate residuals for ALL points
test_merged_df$residuals_comp <- test_merged_df$ESM1v - test_merged_df$fitted_comp
sum(is.na(test_merged_df$residuals_comp)) #0
head(test_merged_df)
```

```{r, fig.width=8, fig.height=4.5}
position_residuals <- test_merged_df %>%
  # Ensure you have the residuals calculated. 
  # If not, fit the model here:
  # mutate(residual = resid(loess(ESM1v ~ f_ddg_pred, span=0.7))) %>%
  group_by(position) %>%
  summarize(
    # We take the median magnitude of residuals at this site
    median_residual = median(residuals_comp, na.rm = TRUE),
    median_ddgb = median(b_ddg_pred, na.rm = TRUE)
  )

# 2. Join with MSA Depth
plot_data <- position_residuals %>%
  inner_join(msa_depth_df, by = "position")

# 3. Correlation Test
# This is the number you will quote in your rebuttal.
cor_test <- cor.test(plot_data$neff_aa, plot_data$median_residual, method = "spearman")
print(cor_test)

positions_to_highlight <- c(173, 183, 203) 

# 2. Prepare data
plot_data <- plot_data %>%
  mutate(
    # Label logic: ONLY label if the position is in your list
    label_text = ifelse(position %in% positions_to_highlight, as.character(position), ""),
    
    # Optional: Create a column to make these points look different (e.g., bigger/bolder)
    is_highlighted = position %in% positions_to_highlight
  ) cvffcdfvd

# 4. Create the Scatter Plot
p2 <- ggplot(plot_data, aes(x = neff_aa, y = median_residual,  color = median_ddgb)) +
  geom_point(aes(color = median_ddgb, alpha = 0.7)) +
  scale_color_gradient2(
    low = "blue", 
    mid = "grey", 
    high = "red", 
    midpoint = 0
  ) +
  theme_classic() +
  labs(
    x = "Neff",
    y = "Median ESM-1v-ddGf residual",
    title = "GRB2-SH3",
    subtitle = paste("Spearman's rho =", round(cor_test$estimate, 3), 
                     "| P-value =", format.pval(cor_test$p.value))
  ) +
  theme(
    text = element_text(size = 12),
    legend.position = "bottom"
  ) +
  geom_text_repel(
    aes(label = label_text),
    size = 4,                # Slightly larger text for emphasis
    fontface = "bold",       # Bold text
    box.padding = 0.5,
    max.overlaps = Inf,      # Force it to show all your specific labels
    point.padding = 0.3,
    min.segment.length = 0   # Always draw a line to the point
  ) 

plot_grid(p1,p2,ncol=2)
```
```{r}

```