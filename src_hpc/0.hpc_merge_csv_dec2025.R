
suppressPackageStartupMessages({
  library(data.table)
  library(future)
  library(future.apply)
})

## Prevent data.table from spawning its own threads
setDTthreads(1)

suppressPackageStartupMessages({
  library(data.table)
  library(future)
  library(future.apply)
})

## Prevent data.table from spawning its own threads
setDTthreads(1)

## ---------- CONFIGURE PATHS ----------
BASE <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/00.human_proteome_unzipped"

DIR_CLINVAR <- file.path(BASE, "results_unzipped/meta_clinvar_all_csv_dec2025")
DIR_ESM1V   <- file.path(BASE, "all_ESM1v")
DIR_AM      <- file.path(BASE, "all_AM")

OUT_DIR     <- file.path(BASE, "results_unzipped/merged_out_dec2025")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## ---------- COLUMNS TO KEEP ----------
CLINVAR_COLS_KEEP <- c(
  "variant",
  "Model",
  "Dataset",
  "ddG_pred",
  "new_position",
  "wildtype",
  "mutation",
  "exposure_ss",
  "exposure_plddt",
  "exposure_rasa",
  "spot_disorder",
  "chain",
  "uniprot",
  "clinvar_variation_id",
  "clinvar_clinical_significance",
  "clinvar_review_status"
)

ESM1V_COLS_KEEP <- c(
  "variant",
  "ESM-1v"
)

AM_COLS_KEEP <- c(
  "protein_variant",  # will be renamed to 'variant'
  "am_pathogenicity",
  "am_class"
)

## ---------- SAFE FILE READER ----------
read_subset <- function(path, cols_keep, rename = NULL) {
  if (!file.exists(path)) return(NULL)
  
  header <- fread(path, nrows = 0)
  cols_present <- intersect(cols_keep, names(header))
  if (length(cols_present) == 0L) return(NULL)
  
  dt <- fread(path, select = cols_present)
  
  if (!is.null(rename)) {
    old <- intersect(names(rename), names(dt))
    if (length(old) > 0L) {
      setnames(dt, old = old, new = rename[old])
    }
  }
  dt
}

## ---------- MERGE ONE UNIPROT (clinvar + esm1v + am) ----------
merge_one_uniprot <- function(uniprot_id) {
  message("Processing ", uniprot_id)
  
  # ---- ClinVar ----
  clinvar_path <- file.path(DIR_CLINVAR, paste0(uniprot_id, "_clinvar.csv"))
  clinvar_dt   <- read_subset(clinvar_path, CLINVAR_COLS_KEEP)
  
  # ---- all_ESM1v ----
  esm1v_path <- file.path(DIR_ESM1V, paste0(uniprot_id, ".csv"))
  esm1v_dt   <- read_subset(esm1v_path, ESM1V_COLS_KEEP)
  if (is.null(esm1v_dt)) {
    message("  No ESM1v file for ", uniprot_id)
  } else if ("ESM-1v" %in% names(esm1v_dt)) {
    setnames(esm1v_dt, "ESM-1v", "ESM1v_all_variants")
  }
  
  # ---- all_AM ----
  am_path <- file.path(DIR_AM, paste0(uniprot_id, ".csv"))
  am_dt   <- read_subset(
    am_path,
    AM_COLS_KEEP,
    rename = c("protein_variant" = "variant")
  )
  
  # If all three are missing, nothing to do
  if (is.null(clinvar_dt) && is.null(esm1v_dt) && is.null(am_dt)) {
    message("  No data for ", uniprot_id, ", skipping.")
    return(NULL)
  }
  
  ## pick a starting table (first non-NULL)
  base_dt <- clinvar_dt
  if (is.null(base_dt)) base_dt <- esm1v_dt
  if (is.null(base_dt)) base_dt <- am_dt
  
  setDT(base_dt)
  
  if (!is.null(esm1v_dt) && !identical(base_dt, esm1v_dt)) {
    base_dt <- merge(base_dt, esm1v_dt, by = "variant", all = TRUE)
  }
  
  if (!is.null(am_dt) && !identical(base_dt, am_dt)) {
    base_dt <- merge(base_dt, am_dt, by = "variant", all = TRUE)
  }
  
  ## Ensure AM columns exist even if AM file didn't
  if (!"am_pathogenicity" %in% names(base_dt)) {
    base_dt[, am_pathogenicity := NA_real_]
  }
  if (!"am_class" %in% names(base_dt)) {
    base_dt[, am_class := NA_character_]
  }
  
  ## Ensure ESM1v column exists even if ESM1v file didn't
  if (!"ESM1v_all_variants" %in% names(base_dt)) {
    base_dt[, ESM1v_all_variants := NA_real_]
  }
  
  out_path <- file.path(OUT_DIR, paste0(uniprot_id, "_merged.csv"))
  fwrite(base_dt, out_path)
  
  message("  Wrote ", out_path, " (", nrow(base_dt), " rows)")
  invisible(out_path)
}

## ---------- DISCOVER UNIPROT IDS FROM ALL THREE DIRS ----------
get_ids_from_dir <- function(dir, clinvar = FALSE) {
  if (!dir.exists(dir)) return(character(0))
  files <- list.files(dir, pattern = "\\.csv$", full.names = FALSE)
  if (clinvar) {
    sub("_clinvar\\.csv$", "", files)
  } else {
    sub("\\.csv$", "", files)
  }
}

ids_clinvar <- get_ids_from_dir(DIR_CLINVAR, clinvar = TRUE)
ids_esm1v   <- get_ids_from_dir(DIR_ESM1V,   clinvar = FALSE)
ids_am      <- get_ids_from_dir(DIR_AM,      clinvar = FALSE)

uniprot_ids <- sort(unique(c(ids_clinvar, ids_esm1v, ids_am)))
message("Found ", length(uniprot_ids), " UniProt IDs across ClinVar / ESM1v / AM")

## ---------- PARALLEL RUN (6 cores → 5 workers) ----------
plan(multisession, workers = 5)

res <- future_lapply(
  uniprot_ids,
  function(uid) {
    tryCatch(
      merge_one_uniprot(uid),
      error = function(e) {
        message("ERROR for ", uid, ": ", conditionMessage(e))
        NULL
      }
    )
  }
)

##check ESM1v or AM score in merged is the same as those in individual files 

## ---------- CONFIGURE PATHS (same as before) ----------
BASE <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/00.human_proteome_unzipped"

DIR_ESM1V   <- file.path(BASE, "all_ESM1v")
DIR_AM      <- file.path(BASE, "all_AM")
DIR_MERGED  <- file.path(BASE, "results_unzipped/merged_out_dec2025")

## ---------- CHECK ONE UNIPROT ----------
check_uniprot_scores <- function(uniprot_id, tol = 1e-8) {
  message("Checking ", uniprot_id, " ...")
  
  merged_path <- file.path(DIR_MERGED, paste0(uniprot_id, "_merged.csv"))
  if (!file.exists(merged_path)) {
    stop("Merged file not found for ", uniprot_id)
  }
  
  merged_dt <- fread(merged_path)
  
  ## ---------------- ESM1v CHECK ----------------
  esm1v_path <- file.path(DIR_ESM1V, paste0(uniprot_id, ".csv"))
  esm1v_ok <- NA
  esm1v_mismatches <- data.table()
  
  if (file.exists(esm1v_path)) {
    esm_dt <- fread(esm1v_path, select = c("variant", "ESM-1v"))
    setnames(esm_dt, "ESM-1v", "ESM1v_original")
    
    # join on variant
    esm_merged <- merge(
      merged_dt[, .(variant, ESM1v_all_variants)],
      esm_dt,
      by = "variant",
      all.x = FALSE,
      all.y = FALSE
    )
    
    # numeric comparison with tolerance
    esm_diff <- abs(esm_merged$ESM1v_all_variants - esm_merged$ESM1v_original)
    esm1v_ok <- all(esm_diff <= tol | (is.na(esm_merged$ESM1v_all_variants) & is.na(esm_merged$ESM1v_original)))
    
    if (!esm1v_ok) {
      esm1v_mismatches <- esm_merged[esm_diff > tol | xor(is.na(ESM1v_all_variants), is.na(ESM1v_original))]
    }
    
    message("  ESM1v check: ", if (esm1v_ok) "OK" else "MISMATCHES")
  } else {
    message("  No ESM1v file for ", uniprot_id, " (skipping ESM1v check)")
  }
  
  ## ---------------- AM CHECK ----------------
  am_path <- file.path(DIR_AM, paste0(uniprot_id, ".csv"))
  am_ok <- NA
  am_mismatches <- data.table()
  
  if (file.exists(am_path)) {
    am_dt <- fread(am_path, select = c("protein_variant", "am_pathogenicity", "am_class"))
    setnames(am_dt, "protein_variant", "variant")
    setnames(am_dt, c("am_pathogenicity", "am_class"),
             c("am_path_original", "am_class_original"))
    
    am_merged <- merge(
      merged_dt[, .(variant, am_pathogenicity, am_class)],
      am_dt,
      by = "variant",
      all.x = FALSE,
      all.y = FALSE
    )
    
    path_diff <- abs(am_merged$am_pathogenicity - am_merged$am_path_original)
    path_ok <- all(path_diff <= tol | (is.na(am_merged$am_pathogenicity) & is.na(am_merged$am_path_original)))
    class_ok <- all(am_merged$am_class == am_merged$am_class_original | (is.na(am_merged$am_class) & is.na(am_merged$am_class_original)))
    
    am_ok <- path_ok && class_ok
    
    if (!am_ok) {
      am_mismatches <- am_merged[
        path_diff > tol |
          xor(is.na(am_pathogenicity), is.na(am_path_original)) |
          am_class != am_class_original
      ]
    }
    
    message("  AM check: ", if (am_ok) "OK" else "MISMATCHES")
  } else {
    message("  No AM file for ", uniprot_id, " (skipping AM check)")
  }
  
  invisible(list(
    esm1v_ok = esm1v_ok,
    esm1v_mismatches = esm1v_mismatches,
    am_ok = am_ok,
    am_mismatches = am_mismatches
  ))
}

res_A0A075B6I1 <- check_uniprot_scores("A0A075B6I1")
res_A0A075B6I1$esm1v_ok
res_A0A075B6I1$am_ok
res_A0A075B6I1$esm1v_mismatches
res_A0A075B6I1$am_mismatches








