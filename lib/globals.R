# ── Project-Wide Constants ────────────────────────────────────────────────────
# Shared across all analysis scripts. Loaded automatically via load.project().

# ── ProjectTemplate configuration ──
add.config(apply.override = FALSE)
add.config(apply.override = TRUE)

# ── Directory paths ──
BASE_DIR <- "/Users/xl7/Documents/0.Projects/01.protein-seq-evo-v1"
SUPPLEMENTS_DIR <- "/Users/xl7/Documents/0.Projects/00.large_supplements"

# ── Numerical constants ──
EPSILON <- 1e-10           # Small value to avoid Inf in logit transforms
PDB_OUTLIER_BFACTOR <- 999 # Sentinel value for unmapped residues in PDB B-factor coloring

# ── Amino acid code lookup (3-letter to 1-letter) ──
AA_CODES <- c(
  "ALA" = "A", "CYS" = "C", "ASP" = "D", "GLU" = "E", "PHE" = "F",
  "GLY" = "G", "HIS" = "H", "ILE" = "I", "LYS" = "K", "LEU" = "L",
  "MET" = "M", "ASN" = "N", "PRO" = "P", "GLN" = "Q", "ARG" = "R",
  "SER" = "S", "THR" = "T", "VAL" = "V", "TRP" = "W", "TYR" = "Y"
)

# ── Color palettes ──

# Site classification (orthosteric / allosteric / other)
SITE_COLORS <- c(
  "orthosteric"     = "#E69F00",
  "allosteric"      = "#4DBBD5",
  "other"           = "#999999",
  "non-orthosteric" = "#8DAA9D"
)

# ClinVar classification
CLINVAR_COLORS <- c(
  "Benign"     = "grey70",
  "Pathogenic" = "firebrick3",
  "Other"      = "grey30"
)

# ggMarginal defaults (used across all scatter plots)
MARGINAL_DEFAULTS <- list(
  type = "density",
  margins = "both",
  groupColour = FALSE,
  groupFill = FALSE,
  size = 10,
  colour = "grey",
  fill = "lightgrey"
)
