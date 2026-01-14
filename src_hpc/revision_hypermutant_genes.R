---
title: "Untitled"
output: html_document
date: "2026-01-11"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
# Load necessary libraries
library(lme4)
library(data.table)
library(future.apply)

# Set up parallelization
future::plan(multisession)
```

```{r}
df <- read.csv("/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/00.human_proteome_results/all_var_lm_r2.csv")
head(df)

head(df %>% arrange(desc(adj_r2)))

protein_atlas_anno <- read.csv("/lustre/scratch126/gengen/projects_v2/alpha-allostery-global/allostery_pathogenicity/00.human_proteome_inputs/protein_atlas_meta/idmapping_gene_name.tsv")
head(protein_atlas_anno)
nrow(protein_atlas_anno)

nrow(df)

df <- merge(df, protein_atlas_anno, by.x = "protein_id", by.y = "uniprot")
df <- df %>% arrange(desc(adj_r2))
fil_df <- df %>% filter (source_file %in% c("T-cell_receptor_genes", "immunoglobulin_genes"))
fil_df <- fil_df %>% arrange(desc(adj_r2))
head(fil_df)
```

```{r}

```

```{r}

```












