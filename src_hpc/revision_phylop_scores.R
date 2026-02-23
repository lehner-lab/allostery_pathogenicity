
#infile <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/00.human_proteome_inputs/cadd/whole_genome_SNVs_inclAnno.tsv.gz"

#system2(
#  "pigz",
#  args = c("-d", "-p", parallel::detectCores(), infile)
#)

library(data.table)

cadd <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/00.human_proteome_inputs/cadd/whole_genome_SNVs_inclAnno.tsv.gz"

# 0) Point this to a GRCh38 GTF on cluster (GENCODE recommended)
gtf_path <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/00.human_proteome_inputs/gencode.v45.annotation.gtf.gz"

# Pull the first "gene" feature line for gene_name "GCK"
cmd <- sprintf(
  "zcat %s | awk -F'\t' 'BEGIN{OFS=\"\t\"} $0!~/^##/ && $3==\"gene\" && $0 ~ /gene_name \"GCK\"/ {print; exit}'",
  shQuote(gtf_path)
)

gene_line <- system(cmd, intern = TRUE)
gene_line
stopifnot(length(gene_line) == 1)

fields <- strsplit(gene_line, "\t")[[1]]
gtf_chr   <- fields[1]
gtf_start <- as.integer(fields[4])
gtf_end   <- as.integer(fields[5])

# CADD contigs are "1..22" (no chr prefix), so normalize
chr <- sub("^chr", "", gtf_chr)

message("GCK (from GTF): chr=", chr, " start=", gtf_start, " end=", gtf_end)

# GCK coords (already computed)
chr <- "7"
gtf_start <- "44143213"
gtf_end   <- "44198170"

cadd_chr7 <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/00.human_proteome_inputs/cadd/whole_genome_SNVs_inclAnno.chr7.tsv.gz"

gtf_start <- 44143213
gtf_end   <- 44198170

# tabix whole_genome_SNVs_inclAnno.tsv.gz 7:44143213-44198170 > cadd_GCK.tsv

cadd_full <- "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/00.human_proteome_inputs/cadd/whole_genome_SNVs_inclAnno.tsv.gz"

hdr_line <- system(
  paste0("zcat ", shQuote(cadd_full), " 2>/dev/null | grep -m 1 -v '^##'"),
  intern = TRUE
)

stopifnot(length(hdr_line) == 1)

colnames_cadd <- strsplit(hdr_line, "\t")[[1]]
length(colnames_cadd)

dt <- fread(
  "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/00.human_proteome_inputs/cadd/cadd_GCK.tsv",
  sep = "\t",
  header = FALSE
)

stopifnot(ncol(dt) == length(colnames_cadd))
setnames(dt, colnames_cadd)

head(dt)
unique(dt$GeneName)

nrow(dt)
dt <- dt %>% filter (GeneName == "GCK")
nrow(dt)

# choose PhyloP priPhyloP / mamPhyloP / verPhyloP
phy_col <- "verPhyloP"   # or "mamPhyloP" / "verPhyloP"

dt[, hgvsp_like := fifelse(
  !is.na(oAA) & !is.na(nAA) & !is.na(protPos),
  paste0(oAA, protPos, nAA),
  NA_character_
)]

unique(dt$hgvsp_like)

fil_dt <- dt %>% filter (!is.na(hgvsp_like))
fil_dt

unique(fil_dt$hgvsp_like)

fil_dt <- fil_dt[!grepl("\\*", hgvsp_like)]
fil_dt

nrow(fil_dt) #4040

aa_change_phylop <- fil_dt[
  !is.na(hgvsp_like) & !is.na(get(phy_col)),
  .(
    n_snv = .N,
    phylop_mean = mean(get(phy_col), na.rm = TRUE),
    phylop_median = median(get(phy_col), na.rm = TRUE),
    phylop_max = max(get(phy_col), na.rm = TRUE)
  ),
  by = .(GeneName, FeatureID, hgvsp_like, protPos)
][order(GeneName, FeatureID, hgvsp_like, protPos)]

nrow(aa_change_phylop) #3179
unique(aa_change_phylop$protPos)
aa_change_phylop

fwrite(
  aa_change_phylop,
  file = "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/git_allostery_pathogenicity/data/revision/phyloP/gck_verPhyloP_hg38.csv"
)

#################################################################################

cmd <- sprintf(
  "zcat %s | awk -F'\t' 'BEGIN{OFS=\"\t\"} $0!~/^##/ && $3==\"gene\" && $0 ~ /gene_name \"PDK1\"/ {print; exit}'",
  shQuote(gtf_path)
)

gene_line <- system(cmd, intern = TRUE)
gene_line
stopifnot(length(gene_line) == 1)

fields <- strsplit(gene_line, "\t")[[1]]
gtf_chr   <- fields[1]
gtf_start <- as.integer(fields[4])
gtf_end   <- as.integer(fields[5])

# CADD contigs are "1..22" (no chr prefix), so normalize
chr <- sub("^chr", "", gtf_chr)

message("PDK1 (from GTF): chr=", chr, " start=", gtf_start, " end=", gtf_end)
# PDK1 (from GTF): chr=2 start=172555373 end=172608669

# tabix whole_genome_SNVs_inclAnno.tsv.gz 2:172555373-172608669 > cadd_PDK1.tsv

hdr_line <- system(
  paste0("zcat ", shQuote(cadd_full), " 2>/dev/null | grep -m 1 -v '^##'"),
  intern = TRUE
)

stopifnot(length(hdr_line) == 1)

colnames_cadd <- strsplit(hdr_line, "\t")[[1]]
length(colnames_cadd)

dt <- fread(
  "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/00.human_proteome_inputs/cadd/cadd_PDK1.tsv",
  sep = "\t",
  header = FALSE
)

stopifnot(ncol(dt) == length(colnames_cadd))
setnames(dt, colnames_cadd)

head(dt)

nrow(dt)
dt <- dt %>% filter (GeneName == "PDK1")
nrow(dt)

# choose PhyloP priPhyloP / mamPhyloP / verPhyloP
phy_col <- "verPhyloP"   # or "mamPhyloP" / "verPhyloP"

dt[, hgvsp_like := fifelse(
  !is.na(oAA) & !is.na(nAA) & !is.na(protPos),
  paste0(oAA, protPos, nAA),
  NA_character_
)]

unique(dt$hgvsp_like)

fil_dt <- dt %>% filter (!is.na(hgvsp_like))
fil_dt

unique(fil_dt$hgvsp_like)

fil_dt <- fil_dt[!grepl("\\*", hgvsp_like)]
fil_dt

nrow(fil_dt) #3758

aa_change_phylop <- fil_dt[
  !is.na(hgvsp_like) & !is.na(get(phy_col)),
  .(
    n_snv = .N,
    phylop_mean = mean(get(phy_col), na.rm = TRUE),
    phylop_median = median(get(phy_col), na.rm = TRUE),
    phylop_max = max(get(phy_col), na.rm = TRUE)
  ),
  by = .(GeneName, FeatureID, hgvsp_like, protPos)
][order(GeneName, FeatureID, hgvsp_like, protPos)]

nrow(aa_change_phylop) #2992
unique(aa_change_phylop$protPos)
aa_change_phylop

fwrite(
  aa_change_phylop,
  file = "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/git_allostery_pathogenicity/data/revision/phyloP/PDK1_verPhyloP_hg38.csv"
)

#################################################################################

cmd <- sprintf(
  "zcat %s | awk -F'\t' 'BEGIN{OFS=\"\t\"} $0!~/^##/ && $3==\"gene\" && $0 ~ /gene_name \"CASP1\"/ {print; exit}'",
  shQuote(gtf_path)
)

gene_line <- system(cmd, intern = TRUE)
gene_line
stopifnot(length(gene_line) == 1)

fields <- strsplit(gene_line, "\t")[[1]]
gtf_chr   <- fields[1]
gtf_start <- as.integer(fields[4])
gtf_end   <- as.integer(fields[5])

# CADD contigs are "1..22" (no chr prefix), so normalize
chr <- sub("^chr", "", gtf_chr)

message("CASP1 (from GTF): chr=", chr, " start=", gtf_start, " end=", gtf_end)
# CASP1 (from GTF): chr=11 start=105025397 end=105035250

# tabix whole_genome_SNVs_inclAnno.tsv.gz 11:105025397-105035250 > cadd_CASP1.tsv

dt <- fread(
  "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/00.human_proteome_inputs/cadd/cadd_CASP1.tsv",
  sep = "\t",
  header = FALSE
)

stopifnot(ncol(dt) == length(colnames_cadd))
setnames(dt, colnames_cadd)

head(dt)

nrow(dt)
dt <- dt %>% filter (GeneName == "CASP1")
nrow(dt)

# choose PhyloP priPhyloP / mamPhyloP / verPhyloP
phy_col <- "verPhyloP"   # or "mamPhyloP" / "verPhyloP"

dt[, hgvsp_like := fifelse(
  !is.na(oAA) & !is.na(nAA) & !is.na(protPos),
  paste0(oAA, protPos, nAA),
  NA_character_
)]

unique(dt$hgvsp_like)

fil_dt <- dt %>% filter (!is.na(hgvsp_like))
fil_dt

unique(fil_dt$hgvsp_like)

fil_dt <- fil_dt[!grepl("\\*", hgvsp_like)]
fil_dt

nrow(fil_dt) #3485

aa_change_phylop <- fil_dt[
  !is.na(hgvsp_like) & !is.na(get(phy_col)),
  .(
    n_snv = .N,
    phylop_mean = mean(get(phy_col), na.rm = TRUE),
    phylop_median = median(get(phy_col), na.rm = TRUE),
    phylop_max = max(get(phy_col), na.rm = TRUE)
  ),
  by = .(GeneName, FeatureID, hgvsp_like, protPos)
][order(GeneName, FeatureID, hgvsp_like, protPos)]

nrow(aa_change_phylop) #2780
unique(aa_change_phylop$protPos)
aa_change_phylop

fwrite(
  aa_change_phylop,
  file = "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/git_allostery_pathogenicity/data/revision/phyloP/CASP1_verPhyloP_hg38.csv"
)

#################################################################################

cmd <- sprintf(
  "zcat %s | awk -F'\t' 'BEGIN{OFS=\"\t\"} $0!~/^##/ && $3==\"gene\" && $0 ~ /gene_name \"DLG4\"/ {print; exit}'",
  shQuote(gtf_path)
)

gene_line <- system(cmd, intern = TRUE)
gene_line
stopifnot(length(gene_line) == 1)

fields <- strsplit(gene_line, "\t")[[1]]
gtf_chr   <- fields[1]
gtf_start <- as.integer(fields[4])
gtf_end   <- as.integer(fields[5])

# CADD contigs are "1..22" (no chr prefix), so normalize
chr <- sub("^chr", "", gtf_chr)

message("DLG4 (from GTF): chr=", chr, " start=", gtf_start, " end=", gtf_end)
# DLG4 (from GTF): chr=17 start=7187187 end=7219841

# tabix whole_genome_SNVs_inclAnno.tsv.gz 17:7187187-7219841 > cadd_DLG4.tsv

dt <- fread(
  "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/00.human_proteome_inputs/cadd/cadd_DLG4.tsv",
  sep = "\t",
  header = FALSE
)

stopifnot(ncol(dt) == length(colnames_cadd))
setnames(dt, colnames_cadd)

head(dt)

nrow(dt)
dt <- dt %>% filter (GeneName == "DLG4")
nrow(dt)

# choose PhyloP priPhyloP / mamPhyloP / verPhyloP
phy_col <- "verPhyloP"   # or "mamPhyloP" / "verPhyloP"

dt[, hgvsp_like := fifelse(
  !is.na(oAA) & !is.na(nAA) & !is.na(protPos),
  paste0(oAA, protPos, nAA),
  NA_character_
)]

unique(dt$hgvsp_like)

fil_dt <- dt %>% filter (!is.na(hgvsp_like))
fil_dt

unique(fil_dt$hgvsp_like)

fil_dt <- fil_dt[!grepl("\\*", hgvsp_like)]
fil_dt

nrow(fil_dt) #6290

aa_change_phylop <- fil_dt[
  !is.na(hgvsp_like) & !is.na(get(phy_col)),
  .(
    n_snv = .N,
    phylop_mean = mean(get(phy_col), na.rm = TRUE),
    phylop_median = median(get(phy_col), na.rm = TRUE),
    phylop_max = max(get(phy_col), na.rm = TRUE)
  ),
  by = .(GeneName, FeatureID, hgvsp_like, protPos)
][order(GeneName, FeatureID, hgvsp_like, protPos)]

nrow(aa_change_phylop) #5038
unique(aa_change_phylop$protPos)
aa_change_phylop

fwrite(
  aa_change_phylop,
  file = "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/git_allostery_pathogenicity/data/revision/phyloP/DLG4_verPhyloP_hg38.csv"
)

#################################################################################

cmd <- sprintf(
  "zcat %s | awk -F'\t' 'BEGIN{OFS=\"\t\"} $0!~/^##/ && $3==\"gene\" && $0 ~ /gene_name \"GRB2\"/ {print; exit}'",
  shQuote(gtf_path)
)

gene_line <- system(cmd, intern = TRUE)
gene_line
stopifnot(length(gene_line) == 1)

fields <- strsplit(gene_line, "\t")[[1]]
gtf_chr   <- fields[1]
gtf_start <- as.integer(fields[4])
gtf_end   <- as.integer(fields[5])

# CADD contigs are "1..22" (no chr prefix), so normalize
chr <- sub("^chr", "", gtf_chr)

message("GRB2 (from GTF): chr=", chr, " start=", gtf_start, " end=", gtf_end)
# GRB2 (from GTF): chr=17 start=75318076 end=75405709

# tabix whole_genome_SNVs_inclAnno.tsv.gz 17:75318076-75405709 > cadd_GRB2.tsv

dt <- fread(
  "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/00.human_proteome_inputs/cadd/cadd_GRB2.tsv",
  sep = "\t",
  header = FALSE
)

stopifnot(ncol(dt) == length(colnames_cadd))
setnames(dt, colnames_cadd)

head(dt)

nrow(dt)
dt <- dt %>% filter (GeneName == "GRB2")
nrow(dt)

# choose PhyloP priPhyloP / mamPhyloP / verPhyloP
phy_col <- "verPhyloP"   # or "mamPhyloP" / "verPhyloP"

dt[, hgvsp_like := fifelse(
  !is.na(oAA) & !is.na(nAA) & !is.na(protPos),
  paste0(oAA, protPos, nAA),
  NA_character_
)]

unique(dt$hgvsp_like)

fil_dt <- dt %>% filter (!is.na(hgvsp_like))
fil_dt

unique(fil_dt$hgvsp_like)

fil_dt <- fil_dt[!grepl("\\*", hgvsp_like)]
fil_dt

nrow(fil_dt) #1871

aa_change_phylop <- fil_dt[
  !is.na(hgvsp_like) & !is.na(get(phy_col)),
  .(
    n_snv = .N,
    phylop_mean = mean(get(phy_col), na.rm = TRUE),
    phylop_median = median(get(phy_col), na.rm = TRUE),
    phylop_max = max(get(phy_col), na.rm = TRUE)
  ),
  by = .(GeneName, FeatureID, hgvsp_like, protPos)
][order(GeneName, FeatureID, hgvsp_like, protPos)]

nrow(aa_change_phylop) #1501
unique(aa_change_phylop$protPos)
aa_change_phylop

fwrite(
  aa_change_phylop,
  file = "/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/git_allostery_pathogenicity/data/revision/phyloP/GRB2_verPhyloP_hg38.csv"
)

#################################################################################











