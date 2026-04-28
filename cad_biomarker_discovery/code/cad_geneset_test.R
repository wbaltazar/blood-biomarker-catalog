############################################################
# Gene-set testing
# 
# What this script does
# 1) Loads five cohorts (GSE12288 microarray; GSE202625 RNA‑seq counts; 
#    GSE180081 RPKM Excel; GSE20680 microarray; and GSE20681 microarray), 
#    maps to Ensembl IDs, and puts RNA‑seq datasets on a common scale 
#    (voom logCPM for counts; TPM->log for RPKM).
# 2) Intersects genes across cohorts and the common genesets from our original datasets.
# 3) Evaluates Leave‑One‑Dataset‑Out (LODO) generalization with optional
#    training‑distribution matching: anchor‑control z‑scoring (+ optional fSVA).
# 4) Compares gene set performance across datasets: TRAIT, STATE, CAD, and ALL-driven.
# 
#
# Notes
# • Replace file paths under "USER PATHS".
# • Requires Bioconductor/CRAN packages listed below.
# • fSVA is optional but recommended; set USE_FSVA <- TRUE to enable.
############################################################

########## Packages ##########
suppressPackageStartupMessages({
  library(tidyverse)
  library(GEOquery)
  library(affy)
  library(hgu133a.db)
  library(jetset)
  library(limma)
  library(edgeR)
  library(readxl)
  library(biomaRt)
  library(AnnotationDbi)
  library(caret)
  library(pROC)
  library(VennDiagram)
  library(mltools)
  library(sva)            # sva + fsva
})

set.seed(1234)

########## USER PATHS ##########
# Raw files you already have locally
CEL_DIR        <- "~/Downloads/GSE12288_RAW/"                  # CEL files for GSE12288
GSE2_COUNTS    <- "~/Downloads/GSE202625_Count_Matrix.txt.gz"  # counts (Geneid + samples)
GSE3_XLSX      <- "~/Downloads/GSE180081_CathDx_HG38_realign_RPKMs.xlsx"  # RPKM matrix

# Gene set files (by SYMBOL)
TRAIT_CSV      <- "~/Desktop/work_repo/newtest2025/trait_gene_list.csv"   # column x
STATE_CSV      <- "~/Desktop/work_repo/newtest2025/state_gene_list.csv"   # column x
CAD_TABLE      <- "~/Desktop/work_repo/newtest2025/filtered_gtex_gwas_table.csv" # column Symbol.of.blood.RNA

# Toggles
USE_FSVA       <- TRUE   # frozen SVA on training data (recommended)
INNER_CV       <- TRUE   # self performance via internal 70/30 split per dataset

############################################################
# Functions for loading datasets and caret ----
############################################################

msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n")

# caret helpers
classes <- c("control","case")
ctrl_caret <- trainControl(method = "repeatedcv", number = 5, repeats = 10,
                           classProbs = TRUE, summaryFunction = twoClassSummary)

fit_glmnet <- function(X_samp_by_feat, y, alpha){
  
  # Get a sensible lambda sequence from glmnet itself
  x_mat <- as.matrix(X_samp_by_feat)
  y_num <- ifelse(y == "case", 1, 0)
  
  glmnet_fit <- glmnet::glmnet(
    x = x_mat,
    y = y_num,
    family = "binomial",
    alpha = alpha
  )
  lambda_seq <- glmnet_fit$lambda
  
  tuneGrid <- expand.grid(
    alpha = alpha,
    lambda = lambda_seq
  )
  
  caret::train(
    x = X_samp_by_feat, y = y,
    method = "glmnet",
    family = "binomial",
    trControl = ctrl_caret,
    metric = "ROC",
    tuneGrid = tuneGrid
  )
}

auc_vec <- function(model, X, y){
  pr <- predict(model, X, type = "prob")[,"case"]
  pROC::auc(y, pr) %>% as.numeric()
}

# auc_from_prob <- function(y_true, prob_case) {
#   as.numeric(
#     pROC::auc(
#       response = y_true,
#       predictor = prob_case,
#       levels = c("control", "case"),
#       direction = "<"
#     )
#   )
# }

# Anchor control z‑scoring using TRAIN controls only
anchor_scale <- function(X_train, y_train, X_test){
  stopifnot(all(levels(y_train) == classes))
  ctrl_idx <- which(y_train == "control")
  mu  <- rowMeans(X_train[, ctrl_idx, drop=FALSE])
  sdv <- apply( X_train[, ctrl_idx, drop=FALSE], 1, sd ); sdv[sdv==0] <- 1
  Ztr <- sweep(sweep(X_train, 1, mu, "-"), 1, sdv, "/")
  Zte <- sweep(sweep(X_test,  1, mu, "-"), 1, sdv, "/")
  list(train = Ztr, test = Zte)
}

# fSVA wrapper (frozen SVA on training; project to test)
apply_fsva <- function(X_train, y_train, X_test){
  mod  <- model.matrix(~ y_train)
  mod0 <- model.matrix(~ 1, data = data.frame(y_train))
  svobj <- sva(as.matrix(X_train), mod = mod, mod0 = mod0)
  tr  <- fsva(dbdat = as.matrix(X_train), mod = mod, sv = svobj, newdat = as.matrix(X_train))$new
  te  <- fsva(dbdat = as.matrix(X_train), mod = mod, sv = svobj, newdat = as.matrix(X_test))$new
  list(train = tr, test = te)
}

# fSVA wrapper (frozen SVA on training; project to test)
apply_fsva <- function(X_train, y_train, X_test){
  mod  <- model.matrix(~ y_train)
  mod0 <- model.matrix(~ 1, data = data.frame(y_train))
  svobj <- sva(as.matrix(X_train), mod = mod, mod0 = mod0)
  tr  <- fsva(dbdat = as.matrix(X_train), mod = mod, sv = svobj, newdat = as.matrix(X_train))$new
  te  <- fsva(dbdat = as.matrix(X_train), mod = mod, sv = svobj, newdat = as.matrix(X_test))$new
  list(train = tr, test = te)
}

# Collapse by Ensembl (mean of duplicates)
collapse_by_ensembl <- function(expr, ensembl_vec){
  stopifnot(nrow(expr) == length(ensembl_vec))
  df <- cbind.data.frame(ensembl = ensembl_vec, expr)
  df %>% dplyr::filter(!is.na(ensembl) & ensembl != "") %>%
    dplyr::group_by(ensembl) %>%
    dplyr::summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    as.data.frame() %>% tibble::column_to_rownames("ensembl")
}

############################################################
# Gene ID mapping helpers (biomaRt)
############################################################

mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
attrs <- c("ucsc", "ensembl_transcript_id_version", "ensembl_gene_id",
           "affy_hg_u133a_2", "hgnc_symbol", "external_gene_name", "gene_biotype")

map_affy_to_ensembl <- function(affy_ids){
  tbl <- getBM(attributes = attrs, filters = "affy_hg_u133a_2", values = affy_ids, mart = mart)
  
  # Jetset best probeset per Ensembl gene
  keep <- tbl %>% dplyr::select(affy_hg_u133a_2, ensembl_gene_id) %>% distinct()
  j <- jetset::jmap(chip = "hgu133a", ensembl = keep$ensembl_gene_id)
  j <- j[!is.na(j)]
  best <- tibble::tibble(ensembl = names(j), affy = j)
  best
}

map_symbol_to_ensembl <- function(symbols){
  tbl <- getBM(attributes = attrs, filters = "external_gene_name", values = unique(symbols), mart = mart)
  tbl %>% dplyr::select(external_gene_name, ensembl_gene_id) %>% distinct()
}

map_ucsc_to_ensembl <- function(ucsc_ids){
  tbl <- getBM(attributes = attrs, filters = "ucsc", values = unique(ucsc_ids), mart = mart)
  tbl %>% dplyr::select(ucsc, external_gene_name, ensembl_gene_id) %>% distinct()
}

############################################################
# Load datasets
############################################################

### GSE12288 (Affymetrix HG-U133A; use pooled cohort only) ###
msg("Loading GSE12288...")
gse1 <- getGEO("GSE12288")[[1]]
ph1  <- pData(gse1)
# RMA from CELs for pooled cohort (matched men)
idx_pooled <- ph1$`cohort:ch1` == "pooled"
affy_raw   <- ReadAffy(celfile.path = CEL_DIR)
rma_expr   <- affy::rma(affy_raw[, idx_pooled]) %>% exprs()  # probes x samples
colnames(rma_expr) <- gsub("\\.CEL\\.gz", "", colnames(rma_expr))
stopifnot(identical(colnames(rma_expr), ph1$geo_accession[idx_pooled]))

CAD1 <- factor(ph1[idx_pooled, "patient group:ch1"], levels = c("control","case"))

# Map to Ensembl using jetset best probeset
map1 <- map_affy_to_ensembl(rownames(rma_expr))
expr1 <- rma_expr[rownames(rma_expr) %in% map1$affy, ]
expr1 <- collapse_by_ensembl(expr = expr1, ensembl_vec = map1$ensembl[match(rownames(expr1), map1$affy)])

### GSE202625 (RNA-seq counts; TMM -> voom logCPM) ###
msg("Loading GSE202625...")
gse2 <- getGEO("GSE202625")[[1]]
ph2  <- pData(gse2) %>%
  transmute(patient = stringr::str_extract(title, "EP....."),
            sex = `Sex:ch1`, age = `age:ch1`,
            group = `group:ch1`, batch = `sequencing batch:ch1`)

counts2 <- read.delim(GSE2_COUNTS)
stopifnot(colnames(counts2)[1] %in% c("Geneid","gene_id","symbol"))
rownames(counts2) <- counts2[[1]]
counts2 <- counts2[,-1]
counts2 <- as.matrix(counts2)

# Ensure sample order matches pheno
stopifnot(identical(colnames(counts2), ph2$patient))
CAD2 <- factor(ifelse(stringr::str_detect(ph2$group, "^no"), "control","case"), levels = classes)

# edgeR TMM + filter lowly expressed + voom
dge2 <- DGEList(counts = counts2)
dge2 <- calcNormFactors(dge2, method = "TMM")
keep2 <- filterByExpr(dge2, group = CAD2)
dge2 <- dge2[keep2, , keep.lib.sizes = FALSE]
voom2 <- limma::voom(dge2, plot = FALSE)$E  # log2 CPM

# Map SYMBOL -> Ensembl, then collapse
map2  <- map_symbol_to_ensembl(rownames(voom2))
expr2 <- voom2[rownames(voom2) %in% map2$external_gene_name, ]
expr2 <- collapse_by_ensembl(expr = expr2, ensembl_vec = map2$ensembl_gene_id[match(rownames(expr2), map2$external_gene_name)])

### GSE180081 (Excel RPKM; convert to TPM then log) ###
msg("Loading GSE180081...")
gse3 <- getGEO("GSE180081")[[1]]
ph3  <- pData(gse3)
CAD3 <- factor(ifelse(stringr::str_detect(ph3$`cad status:ch1`, "LOW"), "control","case"), levels = classes)

raw3 <- readxl::read_xlsx(GSE3_XLSX, skip = 4)
# Columns
ucsc_ids <- raw3$`UCSC ID`
# Expression columns are the "CDX" columns (like in your original script)
expr_cols <- which(stringr::str_detect(colnames(raw3), "CDX"))
rpkm3 <- as.data.frame(raw3[, expr_cols])
colnames(rpkm3) <- gsub(".raw.", "", colnames(rpkm3))
stopifnot(identical(colnames(rpkm3), ph3$title))

# Convert RPKM -> TPM per sample, then log2(TPM+1)
tpm3 <- apply(rpkm3, 2, function(x){ x / sum(x, na.rm=TRUE) * 1e6 }) %>% as.data.frame()
logtpm3 <- log2(tpm3 + 1)

# Map UCSC -> Ensembl, then collapse
map3  <- map_ucsc_to_ensembl(ucsc_ids)
idx3  <- match(ucsc_ids, map3$ucsc)
ens3  <- map3$ensembl_gene_id[idx3]
expr3 <- collapse_by_ensembl(expr = as.data.frame(logtpm3), ensembl_vec = ens3)


############################################################
# CardioDx datasets
############################################################

## Load in Agilent data
source("~/Desktop/work_repo/github/cad_biomarker_discovery/code/load_agilent_data.R")
# reads in expr_ag1 and expr_ag2

############################################################
# Harmonize: intersect genes across cohorts
############################################################

common_genes <- Reduce(intersect, list(rownames(expr1),
                                       rownames(expr2), 
                                       rownames(expr3), 
                                       rownames(expr_ag1), 
                                       rownames(expr_ag2)))
length(common_genes)
# [1] 6169

## Load in common genes from original study
original_6099 <- read.table("~/Desktop/work_repo/github/cross_study_analysis/output/common_symbols6099.txt")
original_6099_ensmbl <- ens_to_sym[ens_to_sym$symbol%in% original_6099$x, "ensembl"]
common_genes <- intersect(common_genes, original_6099_ensmbl$ensembl)

msg("Common Ensembl genes: %d", length(common_genes))
# [15:09:45]  Common Ensembl genes: 3662 

expr1 <- expr1[common_genes, ]
expr2 <- expr2[common_genes, ]
expr3 <- expr3[common_genes, ]
expr_ag1 <- expr_ag1[common_genes,]
expr_ag2 <- expr_ag2[common_genes,]

# Put into lists
EXPR <- list(GSE1 = expr1, GSE2 = expr2, GSE3 = expr3, GSE4 = expr_ag1, GSE5 = expr_ag2)  # genes x samples
Y    <- list(GSE1 = CAD1,  GSE2 = CAD2,  GSE3 = CAD3, GSE4 = CAD4, GSE5 = CAD5)
BATCH<- list(GSE1 = factor("GSE1", levels=c("GSE1","GSE2","GSE3", "GSE4", "GSE5")),
             GSE2 = factor("GSE2", levels=c("GSE1","GSE2","GSE3", "GSE4", "GSE5")),
             GSE3 = factor("GSE3", levels=c("GSE1","GSE2","GSE3", "GSE4", "GSE5")),
             GSE4 = factor("GSE4", levels=c("GSE1","GSE2","GSE3", "GSE4", "GSE5")),
             GSE5 = factor("GSE5", levels=c("GSE1","GSE2","GSE3", "GSE4", "GSE5")))

############################################################
# Load gene sets (SYMBOL -> Ensembl -> intersect)
############################################################

# Symbol lookup for readability (ENSEMBL → SYMBOL where available)
ens_to_sym <- tibble::tibble(
  ensembl = common_genes,
  symbol  = map2$external_gene_name[match(common_genes, map2$ensembl_gene_id)]
)

common_symbols <- map2$external_gene_name[match(common_genes, map2$ensembl_gene_id)]  # helpful for mapping
sym2ens <- function(symbol_vec){
  mm <- map2$ensembl_gene_id[match(symbol_vec, map2$external_gene_name)]
  mm[!is.na(mm) & mm %in% common_genes]
}

trait_syms <- read.csv(TRAIT_CSV)$x
state_syms <- read.csv(STATE_CSV)$x
cad_syms   <- read.csv(CAD_TABLE)$Symbol.of.blood.RNA %>% unique()
cad_trait_syms   <- read.csv(CAD_TABLE) %>%
  dplyr::filter(Trait.gene. == TRUE) %>%
  dplyr::select(Symbol.of.blood.RNA) %>% 
  unique()

TRAIT <- sym2ens(trait_syms)
STATE <- sym2ens(state_syms)
CAD   <- sym2ens(cad_syms)
CAD_TRAIT <- sym2ens(cad_trait_syms$Symbol.of.blood.RNA)
ALL   <- common_genes

GENESETS <- list(trait = TRAIT, state = STATE, cad = CAD, cad_trait = CAD_TRAIT, all = ALL)
GENESETS <- lapply(GENESETS, function(s) unique(s[!is.na(s)]))

msg("Gene set sizes on common genes: trait=%d, state=%d, cad=%d, cad trait=%d, all=%d",
    length(GENESETS$trait), length(GENESETS$state), length(GENESETS$cad), length(GENESETS$cad_trait), length(GENESETS$all))

# [23:37:22]  Gene set sizes on common genes: trait=813, state=798, cad=105, cad trait=24, all=3662 

############################################################
# Determine dataset specific geneset
############################################################

# Ridge per dataset → genes ranked by |coef|

gene_ranks <- list()  # will hold one tibble per dataset

for (ds in names(EXPR)) {
  msg("Training ridge on %s ...", ds)
  
  X_ds <- EXPR[[ds]]         # genes x samples
  y_ds <- Y[[ds]]            # factor(c("control","case"))
  
  # Train ridge (alpha = 0) using existing helper
  fit_ds <- fit_glmnet(t(X_ds), y = y_ds, alpha = 0)
  
  # Extract coefficients at best lambda
  lam_best <- fit_ds$bestTune$lambda
  coefs    <- as.matrix(coef(fit_ds$finalModel, s = lam_best))
  
  # coefs has an "(Intercept)" row; drop it and align to gene order
  coefs_df <- tibble::tibble(
    term     = rownames(coefs),
    coef     = as.numeric(coefs[,1])
  ) %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::rename(ensembl = term) %>%
    dplyr::mutate(abs_coef = abs(coef)) %>%
    dplyr::arrange(dplyr::desc(abs_coef)) %>%
    dplyr::mutate(rank = dplyr::row_number()) %>%
    dplyr::left_join(ens_to_sym, by = "ensembl") %>%
    dplyr::select(rank, ensembl, symbol, coef, abs_coef)
  
  gene_ranks[[ds]] <- coefs_df
  
  msg("%s: top 5 genes by |beta|:\n%s",
      ds,
      paste(utils::capture.output(print(head(coefs_df, 5))), collapse = "\n"))
}

# Access examples:
# gene_ranks$GSE1$ensembl  # ordered ENSEMBL IDs for GSE12288
# gene_ranks$GSE2          # full tibble with symbols/coefficients
# write.csv(gene_ranks$GSE3, "GSE3_ridge_ranked_genes.csv", row.names = FALSE)

lapply(gene_ranks, head)
lapply(gene_ranks, dim)
# [1] 3662    5 for all of them.
lapply(gene_ranks, function(x){
  coefs <- x$coef
  print(sum(coefs == 0))
})

# $GSE1
# [1] 0
# 
# $GSE2
# [1] 0
# 
# $GSE3
# [1] 0
# 
# $GSE4
# [1] 0
# 
# $GSE5
# [1] 0

############################################################
# Build size-matched, non-overlapping empirical sets (per dataset)
# Rule: take the ranked list for each dataset, skipping any gene in TRAIT,
# until we collect k = |TRAIT|, |STATE|, or |CAD| genes.
############################################################

# Trait gene sets
k <- length(GENESETS$trait)
msg("Target size k = |TRAIT| = %d", k)

EMPIRICAL_TOPK_NONTRAIT <- list()
SKIPPED_TRAIT <- list()

for (ds in names(EXPR)) {
  ranked <- gene_ranks[[ds]]$ensembl  # ranked by |beta| descending (ENSEMBL IDs)
  stopifnot(!is.null(ranked), length(ranked) > 0)
  
  selected <- character(0)
  skipped <- character(0)
  skips <- 0L
  
  for (g in ranked) {
    if (g %in% GENESETS$trait || g %in% selected) {
      skipped <- c(skipped, g)
      skips <- skips + 1L
      next
    }
    selected <- c(selected, g)
    if (length(selected) == k) break
  }
  
  if (length(selected) < k) {
    msg("%s: only %d non-TRAIT genes available (needed %d). Using truncated set.",
        ds, length(selected), k)
  }
  
  SKIPPED_TRAIT[[ds]] <- skipped
  EMPIRICAL_TOPK_NONTRAIT[[ds]] <- selected
  
  # For convenience (named vectors in global env matching your proposed names)
  assign(paste0(ds, "GENES"), selected, inherits = TRUE)
  
  msg("%s: built non-overlapping set of size %d (skipped %d overlaps/dups).",
      ds, length(selected), skips)
}

# [16:07:06]  GSE1: built non-overlapping set of size 813 (skipped 191 overlaps/dups). 
# [16:07:06]  GSE2: built non-overlapping set of size 813 (skipped 211 overlaps/dups). 
# [16:07:06]  GSE3: built non-overlapping set of size 813 (skipped 207 overlaps/dups). 
# [16:07:06]  GSE4: built non-overlapping set of size 813 (skipped 187 overlaps/dups). 
# [16:07:06]  GSE5: built non-overlapping set of size 813 (skipped 176 overlaps/dups). 

mean(c(191, 211, 207, 187, 176))
# [1] 194.4
(194.4 / k) * 100
# [1] 23.91144 percent of trait genes included in each list on average
813/3662 * 813
# 180.494
phyper(q = 194, m = 813, n = (3662 - 813), k = 813, lower.tail = F)
# 0.09072412

# Quick peek (top 10) and optional exports
for (ds in names(EMPIRICAL_TOPK_NONTRAIT)) {
  head_syms <- map2$external_gene_name[match(EMPIRICAL_TOPK_NONTRAIT[[ds]][1:min(10, length(EMPIRICAL_TOPK_NONTRAIT[[ds]]))],
                                             map2$ensembl_gene_id)]
  msg("%s top non-TRAIT genes (first 10): %s", ds, paste(head_syms, collapse = ", "))
}

# [16:08:39]  GSE1 top non-TRAIT genes (first 10): UBAP2L, CEP250, MUS81, PEX14, CHEK1, LIME1, AFG3L2, KCTD2, AKAP8, GSS 
# [16:08:39]  GSE2 top non-TRAIT genes (first 10): ARIH2, TEX261, DHDDS, TERF2, BLCAP, XPNPEP1, COPS5, SEC24C, PREP, MBD1 
# [16:08:39]  GSE3 top non-TRAIT genes (first 10): ITGA10, FZD3, SH3PXD2A, FAM3C, MTRF1L, CDC42EP4, POLR3F, EIF5B, CUL4B, MBTPS2 
# [16:08:39]  GSE4 top non-TRAIT genes (first 10): TWSG1, LIG3, RUFY3, PRKAR2A, PTBP1, ZNF329, SOS1, UBE2V2, CDC14A, NEK1 
# [16:08:39]  GSE5 top non-TRAIT genes (first 10): PPFIBP1, IPO8, TP53, SPRY1, VAMP1, PDE4D, ABL1, DDX31, ARL15, TWSG1 

# Examples:
length(GSE1GENES)  # vector of ENSEMBL IDs (size ~ |TRAIT|)
length(GSE2GENES)
length(GSE3GENES)
length(GSE4GENES)
length(GSE5GENES)

############################################################
# Within-dataset 70/30 evaluation for TRAIT vs GSE*GENES
# Model: Elastic net (alpha = 0.5) with caret
# Metrics: AUROC, Sensitivity, Specificity, MCC
# Requirement: exactly k = length(GENESETS$trait) features for EVERY run
############################################################

# Collect the four feature spaces into a list
feature_sets <- list(
  TRAIT     = GENESETS$trait,
  GSE1GENES = GSE1GENES,
  GSE2GENES = GSE2GENES,
  GSE3GENES = GSE3GENES,
  GSE4GENES = GSE4GENES,
  GSE5GENES = GSE5GENES
)

venn.diagram(x = feature_sets[-1], filename = "~/Desktop/work_repo/github/cad_biomarker_discovery/output/NoTraitSetOverlap.png", 
             disable.logging = T, imagetype = "png", height = 4500, width = 4500, units = "px",
             category.names = c("Set1", "Set2", "Set3", "Set4", "Set5"))
venn.diagram(x = SKIPPED_TRAIT, filename = "~/Desktop/work_repo/github/cad_biomarker_discovery/output/TraitCADRemovedSetOverlap.png", 
             disable.logging = T, imagetype = "png", height = 4500, width = 4500, units = "px",
             category.names = c("Set1", "Set2", "Set3", "Set4", "Set5"))
print(ens_to_sym[ens_to_sym$ensembl %in% Reduce(intersect, feature_sets[-1]), "symbol", drop = TRUE])
print(ens_to_sym[ens_to_sym$ensembl %in% Reduce(intersect, SKIPPED_TRAIT), "symbol", drop = TRUE])

# Validation using different machine learning algorithms
############################################################
# Compare gene sets using 3 algorithms:
# 1. Random Forest (rf)
# 2. Naive Bayes (nb)
# 3. PLS-DA (plsda)
# 4. Elastic net (en)
#
# Output: AUROC for each (dataset × feature_set × algorithm)
# Using 10x repeat, 5-fold CV on 100% of each dataset.
############################################################

set.seed(1234)
k <- length(GENESETS$trait)
algorithms <- c("rf", "nb", "pls", "glmnet")

ctrl_cv <- trainControl(
  method = "repeatedcv",
  repeats = 10,
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE
)

model_results <- list()

for (ds in names(EXPR)) {
  message("=== Cross-validated AUC on ", ds, " ===")
  
  X <- EXPR[[ds]]
  y <- Y[[ds]]
  
  for (fs_name in names(feature_sets)) {
    feats <- feature_sets[[fs_name]]
    
    genes_avail <- intersect(feats, rownames(X))
    if (length(genes_avail) != k)
      stop(sprintf("%s on %s: feature count %d != %d",
                   fs_name, ds, length(genes_avail), k))
    
    X_df <- t(X[genes_avail, , drop = FALSE])  # samples x features
    
    for (alg in algorithms) {
      
      message("Running ", alg, " on ", ds, " using ", fs_name)
      
      if (alg == "pls") {
        fit <- caret::train(
          x = X_df, y = y,
          method = "pls",
          preProcess = c("center", "scale"),
          trControl = ctrl_cv,
          metric = "ROC",
          tuneGrid = expand.grid(ncomp = 1:10)
        )
      } else {
        fit <- caret::train(
          x = X_df, y = y,
          method = alg,
          trControl = ctrl_cv,
          metric = "ROC"
        )
      }
      
      auc_cv <- max(fit$results$ROC, na.rm = TRUE)
      
      model_results[[length(model_results) + 1]] <- tibble::tibble(
        dataset     = ds,
        feature_set = fs_name,
        algorithm   = alg,
        genes_used  = k,
        AUC_5fold10xCV = round(auc_cv, 3)
      )
    }
  }
}

cv_auc_table <- bind_rows(model_results) %>%
  arrange(dataset, algorithm, desc(AUC_5fold10xCV))

print(cv_auc_table)

# Export
write.csv(cv_auc_table, "~/Desktop/work_repo/github/cad_biomarker_discovery/output/Trait_Performance.csv", row.names = FALSE)

# CAD gene sets
k <- length(GENESETS$cad)
msg("Target size k = |CAD| = %d", k)

EMPIRICAL_TOPK_NONCAD <- list()
SKIPPED_CAD <- list()

for (ds in names(EXPR)) {
  ranked <- gene_ranks[[ds]]$ensembl  # ranked by |beta| descending (ENSEMBL IDs)
  stopifnot(!is.null(ranked), length(ranked) > 0)
  
  selected <- character(0)
  skipped <- character(0)
  skips <- 0L
  
  for (g in ranked) {
    if (g %in% GENESETS$cad || g %in% selected) {
      skipped <- c(skipped, g)
      skips <- skips + 1L
      next
    }
    selected <- c(selected, g)
    if (length(selected) == k) break
  }
  
  if (length(selected) < k) {
    msg("%s: only %d non-CAD genes available (needed %d). Using truncated set.",
        ds, length(selected), k)
  }
  
  SKIPPED_CAD[[ds]] <- skipped
  EMPIRICAL_TOPK_NONCAD[[ds]] <- selected
  
  # For convenience (named vectors in global env matching your proposed names)
  assign(paste0(ds, "GENES"), selected, inherits = TRUE)
  
  msg("%s: built non-overlapping set of size %d (skipped %d overlaps/dups).",
      ds, length(selected), skips)
}

# [22:08:46]  GSE1: built non-overlapping set of size 105 (skipped 5 overlaps/dups). 
# [22:08:46]  GSE2: built non-overlapping set of size 105 (skipped 5 overlaps/dups). 
# [22:08:46]  GSE3: built non-overlapping set of size 105 (skipped 5 overlaps/dups). 
# [22:08:46]  GSE4: built non-overlapping set of size 105 (skipped 5 overlaps/dups). 
# [22:08:46]  GSE5: built non-overlapping set of size 105 (skipped 1 overlaps/dups). 

mean(c(5, 5, 5, 5, 1))
# [1] 4.2
(4.2 / k) * 100
# [1] 4 percent of CAD genes included in each list on average
105/3662 * 105
# 3.01065
phyper(q = 4, m = 105, n = (3662 - 105), k = 105, lower.tail = F)
# 0.1815355

# Quick peek (top 10) and optional exports
for (ds in names(EMPIRICAL_TOPK_NONCAD)) {
  head_syms <- map2$external_gene_name[match(EMPIRICAL_TOPK_NONCAD[[ds]][1:min(10, length(EMPIRICAL_TOPK_NONCAD[[ds]]))],
                                             map2$ensembl_gene_id)]
  msg("%s top non-CAD genes (first 10): %s", ds, paste(head_syms, collapse = ", "))
}

# [22:12:07]  GSE1 top non-CAD genes (first 10): UBAP2L, CEP250, HBB, MUS81, PEX14, CHEK1, LIME1, AFG3L2, KCTD2, AKAP8 
# [22:12:07]  GSE2 top non-CAD genes (first 10): ARIH2, TEX261, TERF2, BLCAP, XPNPEP1, COPS5, SEC24C, ARHGEF7, PREP, MBD1 
# [22:12:07]  GSE3 top non-CAD genes (first 10): COL4A3, ITGA10, FABP5, FZD3, FAM3C, MTRF1L, CDC42EP4, POLR3F, EIF5B, CUL4B 
# [22:12:07]  GSE4 top non-CAD genes (first 10): TWSG1, NCAM1, COBLL1, LIG3, RUFY3, PRKAR2A, PTBP1, ZNF329, CDK5R1, SOS1 
# [22:12:07]  GSE5 top non-CAD genes (first 10): PPFIBP1, IPO8, TP53, SPRY1, VAMP1, PDE4D, ABL1, DDX31, ARL15, TWSG1 

# Examples:
length(GSE1GENES)
length(GSE2GENES)
length(GSE3GENES)
length(GSE4GENES)
length(GSE5GENES)

############################################################
# Within-dataset 70/30 evaluation for TRAIT vs GSE*GENES
# Model: Elastic net (alpha = 0.5) with caret
# Metrics: AUROC, Sensitivity, Specificity, MCC
# Requirement: exactly k = length(GENESETS$trait) features for EVERY run
############################################################

# Collect the four feature spaces into a list
feature_sets <- list(
  CAD     = GENESETS$cad,
  GSE1GENES = GSE1GENES,
  GSE2GENES = GSE2GENES,
  GSE3GENES = GSE3GENES,
  GSE4GENES = GSE4GENES,
  GSE5GENES = GSE5GENES
)

venn.diagram(x = feature_sets[-1], filename = "~/Desktop/work_repo/github/cad_biomarker_discovery/output/NoCadSetOverlap.png", 
             disable.logging = T, imagetype = "png", height = 4500, width = 4500, units = "px",
             category.names = c("Set1", "Set2", "Set3", "Set4", "Set5"))
venn.diagram(x = SKIPPED_CAD, filename = "~/Desktop/work_repo/github/cad_biomarker_discovery/output/CadCADRemovedSetOverlap.png", 
             disable.logging = T, imagetype = "png", height = 4500, width = 4500, units = "px",
             category.names = c("Set1", "Set2", "Set3", "Set4", "Set5"))
print(ens_to_sym[ens_to_sym$ensembl %in% Reduce(intersect, feature_sets[-1]), "symbol", drop = TRUE]) # empty
print(ens_to_sym[ens_to_sym$ensembl %in% Reduce(intersect, SKIPPED_CAD), "symbol", drop = TRUE]) # empty

# Validation using different machine learning algorithms
############################################################
# Compare gene sets using 3 algorithms:
# 1. Random Forest (rf)
# 2. Naive Bayes (nb)
# 3. PLS-DA (plsda)
# 4. Elastic net (en)
#
# Output: AUROC for each (dataset × feature_set × algorithm)
# Using 10x repeat, 5-fold CV on 100% of each dataset.
############################################################

set.seed(1234)
k <- length(GENESETS$cad)
algorithms <- c("rf", "nb", "pls", "glmnet")

ctrl_cv <- trainControl(
  method = "repeatedcv",
  repeats = 10,
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE
)

model_results <- list()

for (ds in names(EXPR)) {
  message("=== Cross-validated AUC on ", ds, " ===")
  
  X <- EXPR[[ds]]
  y <- Y[[ds]]
  
  for (fs_name in names(feature_sets)) {
    feats <- feature_sets[[fs_name]]
    
    genes_avail <- intersect(feats, rownames(X))
    if (length(genes_avail) != k)
      stop(sprintf("%s on %s: feature count %d != %d",
                   fs_name, ds, length(genes_avail), k))
    
    X_df <- t(X[genes_avail, , drop = FALSE])  # samples x features
    
    for (alg in algorithms) {
      
      message("Running ", alg, " on ", ds, " using ", fs_name)
      
      if (alg == "pls") {
        fit <- caret::train(
          x = X_df, y = y,
          method = "pls",
          preProcess = c("center", "scale"),
          trControl = ctrl_cv,
          metric = "ROC",
          tuneGrid = expand.grid(ncomp = 1:10)
        )
      } else {
        fit <- caret::train(
          x = X_df, y = y,
          method = alg,
          trControl = ctrl_cv,
          metric = "ROC"
        )
      }
      
      auc_cv <- max(fit$results$ROC, na.rm = TRUE)
      
      model_results[[length(model_results) + 1]] <- tibble::tibble(
        dataset     = ds,
        feature_set = fs_name,
        algorithm   = alg,
        genes_used  = k,
        AUC_5fold10xCV = round(auc_cv, 3)
      )
    }
  }
}

cv_auc_table <- bind_rows(model_results) %>%
  arrange(dataset, algorithm, desc(AUC_5fold10xCV))

print(cv_auc_table)

# Export
write.csv(cv_auc_table, "~/Desktop/work_repo/github/cad_biomarker_discovery/output/CAD_Performance.csv", row.names = FALSE)