############################################################
# Agilent G4112F loader: from GEO "RAW" TXT/TXT.GZ -> Ensembl
# Requirements:
#   - Directories contain the extracted .txt or .txt.gz files
#   - Columns present: accessions, gMedianSignal, gBGMedianSignal, GeneName
#   - (Optional) common_genes: vector of Ensembl IDs to filter to
############################################################

suppressPackageStartupMessages({
  library(limma)
  library(tidyverse)
  library(biomaRt)
})

##### USER: set your directories here #####
AGILENT_DIR_1 <- "~/Desktop/work_repo/data/GSE20680_RAW/"
AGILENT_DIR_2 <- "~/Desktop/work_repo/data/GSE20681_RAW/"  
# Optional: common_genes should already exist from your main script
# e.g., common_genes <- <vector of Ensembl IDs>

############################################################
# Helper: gunzip any .gz files in a directory
############################################################
gunzip_in_dir <- function(path) {
  gz_files <- list.files(path, pattern = "\\.gz$", full.names = TRUE, recursive = TRUE)
  if (length(gz_files) > 0) {
    message("Decompressing ", length(gz_files), " .gz files in ", path)
    system(paste("gunzip", paste(shQuote(gz_files), collapse = " ")))
  }
}

############################################################
# Helper: collapse duplicate Ensembl IDs by mean
############################################################
collapse_by_ensembl <- function(expr, ensembl_vec) {
  stopifnot(nrow(expr) == length(ensembl_vec))
  df <- cbind.data.frame(ensembl = ensembl_vec, expr)
  df %>%
    filter(!is.na(ensembl) & ensembl != "") %>%
    group_by(ensembl) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    as.data.frame() %>%
    tibble::column_to_rownames("ensembl")
}

############################################################
# Helper: map Agilent genes -> Ensembl via SYMBOL
############################################################
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

map_agilent_to_ensembl <- function(genes_df) {
  # Try to find a symbol column
  sym_col <- intersect(c("GeneName", "GENE_SYMBOL", "gGeneName", "Symbol"),
                       colnames(genes_df)) %>% .[1]
  if (is.na(sym_col)) {
    stop("No obvious gene symbol column found in genes_df.")
  }
  
  symbols <- genes_df[[sym_col]] %>% as.character() %>% trimws()
  symbols[symbols == ""] <- NA
  
  map_tbl <- getBM(
    attributes = c("hgnc_symbol", "ensembl_gene_id"),
    filters    = "hgnc_symbol",
    values     = unique(na.omit(symbols)),
    mart       = mart
  )
  
  ensembl_vec <- map_tbl$ensembl_gene_id[match(symbols, map_tbl$hgnc_symbol)]
  return(ensembl_vec)
}

############################################################
# Core: load + normalize one Agilent GEO "RAW" directory
############################################################
load_agilent_dir <- function(path, common_genes = NULL) {
  message("\n=== Processing Agilent directory: ", path, " ===")
  
  # 1) Decompress .gz if present
  gunzip_in_dir(path)
  
  # 2) Collect .txt files
  txt_files <- list.files(path, pattern = "\\.txt$", full.names = TRUE, recursive = TRUE)
  if (length(txt_files) == 0) {
    stop("No .txt files found in ", path)
  }
  txt_files <- txt_files[-which(str_detect(txt_files, "GPL"))]
  
  # 3) Read as *generic* tab-delimited, specifying intensity columns
  RG <- read.maimages(
    txt_files,
    source  = "agilent",
    green.only = TRUE
  )
  
  # 4) Background correction + between-array normalization (one-color)
  RG_corrected <- backgroundCorrect(RG, method = "normexp", offset = 50)
  E_norm <- normalizeBetweenArrays(RG_corrected$E, method = "quantile")
  rownames(E_norm) <- RG$genes$accessions
  
  # 5) Map probes -> Ensembl and collapse duplicates
  ensembl_ids <- map_agilent_to_ensembl(RG$genes)
  expr_ens    <- collapse_by_ensembl(E_norm, ensembl_ids)
  
  message(sprintf("Mapped %d unique Ensembl IDs from %d probes",
                  nrow(expr_ens), length(ensembl_ids)))
  
  # 6) Optionally restrict to common_genes
  if (!is.null(common_genes)) {
    keep <- intersect(rownames(expr_ens), common_genes)
    expr_ens <- expr_ens[keep, , drop = FALSE]
    message(sprintf("Filtered to %d common_genes", length(keep)))
  }
  
  return(expr_ens)  # genes (Ensembl) x samples
}

############################################################
# Run for your two datasets
############################################################

# If you already have `common_genes` defined from your Affy/RNA-seq/RPKM pipeline:
# expr_ag1 <- load_agilent_dir(AGILENT_DIR_1, common_genes)
# expr_ag2 <- load_agilent_dir(AGILENT_DIR_2, common_genes)

# If you're just loading the Agilent data for now (no filter yet):
expr_ag1 <- load_agilent_dir(AGILENT_DIR_1)
expr_ag2 <- load_agilent_dir(AGILENT_DIR_2)
expr_ag1 <- log2(expr_ag1 + 1)
expr_ag2 <- log2(expr_ag2 + 1)

# expr_ag1 / expr_ag2 are now log2-scale, quantile-normalized,
# Ensembl-row expression matrices ready to plug into EXPR / Y lists.

gx1 <- getGEO(GEO = "GSE20680")
px1 <- pData(gx1$GSE20680_series_matrix.txt.gz)
identical(rownames(px1), str_extract(colnames(expr_ag1), "GSM[0-9]+"))
# [1] TRUE
colnames(expr_ag1) <- str_extract(colnames(expr_ag1), "GSM[0-9]+")
CAD4 <- px1$`disease state:ch1`
expr_ag1 <- expr_ag1[, CAD4 %in% c("Control (0)", "Case (2)")]
CAD4 <- CAD4[CAD4 %in% c("Control (0)", "Case (2)")]
length(CAD4)
# [1] 139
CAD4 <- factor(ifelse(CAD4 == "Control (0)", "control", "case"), levels = classes)
table(CAD4)
# CAD4
# control    case 
# 52         87 

gx2 <- getGEO(GEO = "GSE20681")
px2 <- pData(gx2$GSE20681_series_matrix.txt.gz)
identical(rownames(px2), str_extract(colnames(expr_ag2), "GSM[0-9]+"))
# [1] TRUE
colnames(expr_ag2) <- str_extract(colnames(expr_ag2), "GSM[0-9]+")
CAD5 <- px2$`disease state:ch1`
table(CAD5)
# CAD5
# Case (1) Control (0) 
# 99          99 
CAD5 <- factor(ifelse(CAD5 == "Control (0)", "control", "case"), levels = classes)
table(CAD5)
# CAD5
# control    case 
# 99      99 