## Date: Jan 18 2025

## Building a tool for querying your rsIDs!
## Set input directory to where GTEx variants and rsID lookup table are
input_dir <- "~/Desktop/work_repo/data/"
setwd(input_dir)
## Output to the RShiny application directory (not on Github)
output_dir <- "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/eqtls/"

library(tidyverse)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(nanoparquet)
gtex <- read_parquet("./GTEx_Analysis_v10_eQTL_updated/Whole_Blood.v10.eQTLs.signif_pairs.parquet")
gtex$gene_id <- gsub(pattern = "(.*)(\\..*)", replacement = "\\1", x = gtex$gene_id)
gtex$Symbol <- mapIds(EnsDb.Hsapiens.v86, keys = gtex$gene_id, column = "SYMBOL", keytype = "GENEID")

variant_ids <- gtex$variant_id
length(variant_ids)
# [1] 2985690
symbols <- read.table("~/Desktop/work_repo/github/cross_study_analysis/output/common_symbols6099.txt")
symbols <- symbols$x
sum(symbols %in% gtex$Symbol)
# [1] 4452


## Get the rsIDs for the eQTLs so that they can be matched with GWAS traits. SKIP if you have run before!
library(data.table)
id_table <- fread(input = "./GTEx_Analysis_2021-02-11_v10_WholeGenomeSeq_953Indiv.lookup_table.txt.gz", nrows = 2)
id_table
# variant_id    chr   pos    ref    alt num_alt_per_site rs_id_dbSNP155_GRCh38p13  variant_id_b37
# <char> <char> <int> <char> <char>            <int>                   <char>          <char>
#   1: chr1_13550_G_A_b38   chr1 13550      G      A                1              rs554008981 1_13550_G_A_b37
# 2: chr1_14436_G_A_b38   chr1 14436      G      A                1               rs28507908 1_14436_G_A_b37
id_table <- fread(input = "./GTEx_Analysis_2021-02-11_v10_WholeGenomeSeq_953Indiv.lookup_table.txt.gz",
      select = c("variant_id", "rs_id_dbSNP155_GRCh38p13"))

id_table <- id_table %>%
  dplyr::filter(variant_id %in% variant_ids)

# Save the results for later!!
write.table(id_table, "./filtered_rsIDs.txt")

# Start HERE if you already have the id_table
id_table <- read.table("./filtered_rsIDs.txt")
dim(id_table)
# [1] 1500832       2
length(unique(gtex$variant_id))
# [1] 1500832

gtex <- merge(gtex, id_table, by = "variant_id")
dim(gtex)
# [1] 2985690      14

# Next, get the associated traits for the rsIDs.
associations <- read_tsv("./gwas_catalog_v1.0-associations_e113_r2024-11-20.tsv")
names(associations)
# [1] "DATE ADDED TO CATALOG"      "PUBMEDID"                   "FIRST AUTHOR"              
# [4] "DATE"                       "JOURNAL"                    "LINK"                      
# [7] "STUDY"                      "DISEASE/TRAIT"              "INITIAL SAMPLE SIZE"       
# [10] "REPLICATION SAMPLE SIZE"    "REGION"                     "CHR_ID"                    
# [13] "CHR_POS"                    "REPORTED GENE(S)"           "MAPPED_GENE"               
# [16] "UPSTREAM_GENE_ID"           "DOWNSTREAM_GENE_ID"         "SNP_GENE_IDS"              
# [19] "UPSTREAM_GENE_DISTANCE"     "DOWNSTREAM_GENE_DISTANCE"   "STRONGEST SNP-RISK ALLELE" 
# [22] "SNPS"                       "MERGED"                     "SNP_ID_CURRENT"            
# [25] "CONTEXT"                    "INTERGENIC"                 "RISK ALLELE FREQUENCY"     
# [28] "P-VALUE"                    "PVALUE_MLOG"                "P-VALUE (TEXT)"            
# [31] "OR or BETA"                 "95% CI (TEXT)"              "PLATFORM [SNPS PASSING QC]"
# [34] "CNV"


associations <- associations %>% 
  dplyr::select(`DISEASE/TRAIT`, `REPORTED GENE(S)`, MAPPED_GENE, SNPS, CONTEXT, `OR or BETA`,
                `P-VALUE`, PVALUE_MLOG,)
dim(associations)
# [1] 692444      8
head(associations)
# A tibble: 6 × 8
# `DISEASE/TRAIT`           `REPORTED GENE(S)` MAPPED_GENE        SNPS       CONTEXT            `OR or BETA` `P-VALUE` PVALUE_MLOG
# <chr>                     <chr>              <chr>              <chr>      <chr>                     <dbl>     <dbl>       <dbl>
#   1 Inflammatory skin disease LCE1E              LCE1F - LCE1E      rs77199844 intergenic_variant         1.23   2  e-17        16.7
# 2 Inflammatory skin disease TNFAIP3            TNFAIP3            rs643177   intron_variant             1.27   9  e-16        15.0
# 3 Inflammatory skin disease LCE3E              CRCT1 - LCE3E      rs10888499 intergenic_variant         1.49   5  e-25        24.3
# 4 Inflammatory skin disease CTB-11I22.1        LINC01932          rs10515778 intron_variant             1.29   2  e-14        13.7
# 5 Inflammatory skin disease RP1-91G5.3, LCE5A  FLG-AS1 - LCE5A    rs471144   intergenic_variant         1.54   2  e-12        11.7
# 6 Inflammatory skin disease UBLCP1             LINC01932 - UBLCP1 rs11135056 intergenic_variant         1.45   1  e-25        25  

colnames(associations)[4] <- "rs_id_dbSNP155_GRCh38p13"
gtex <- merge(gtex, associations, by = "rs_id_dbSNP155_GRCh38p13", all.x = TRUE)
dim(gtex)
# [1] 3298488      21
colnames(gtex)
# [1] "rs_id_dbSNP155_GRCh38p13" "variant_id"               "gene_id"                  "tss_distance"            
# [5] "af"                       "ma_samples"               "ma_count"                 "pval_nominal"            
# [9] "slope"                    "slope_se"                 "pval_nominal_threshold"   "min_pval_nominal"        
# [13] "pval_beta"                "Symbol"                   "DISEASE/TRAIT"            "REPORTED GENE(S)"        
# [17] "MAPPED_GENE"              "CONTEXT"                  "OR or BETA"               "P-VALUE"                 
# [21] "PVALUE_MLOG"
gtex <- gtex[,c("Symbol", "gene_id", "rs_id_dbSNP155_GRCh38p13", "variant_id", "af", "slope", 
                "pval_nominal", "pval_beta",
                "DISEASE/TRAIT", "MAPPED_GENE", "REPORTED GENE(S)", "OR or BETA", "P-VALUE", "CONTEXT")] # Re-orders the columns
colnames(gtex) <- c("Symbol of blood RNA", "Ensembl ID", "rsID of eQTL", "GTEx variant ID","Alt allele frequency in GTEx WB",
                    "eQTL slope", "eQTL nominal p_val", "eQTL beta p_value",
                    "GWAS Trait", "GWAS Catalog Mapping", "Reported Mapping", "OR or BETA", "GWAS p_value",
                    "Genomic Context") # Re-names the columns

## Attach stable_polymorphic, flexible, and housekeeping gene study statistics
stable <- read.csv("~/Desktop/work_repo/github/cross_study_analysis/output/characteristic_scores.csv")
dynamic <- read.csv("~/Desktop/work_repo/github/cross_study_analysis/output/flexible_gene_scores.csv")
stable <- stable %>% dplyr::select(Symbol, Study_counts)
colnames(stable) <- c("Symbol of blood RNA","Trait 4+ filter studies")
dynamic <- dynamic %>% dplyr::select(Symbol, P_value_study_count)
colnames(dynamic) <- c("Symbol of blood RNA","Studies below 0.05 p_value")

gtex <- merge(gtex, stable, all = T, by = "Symbol of blood RNA", sort = F)
gtex <- merge(gtex, dynamic, all = T, by = "Symbol of blood RNA", sort = F)
dim(gtex)
# [1] 3311615      16

## Filter out by common symbols. Common symbols file is in cross_study_analysis repo
gtex <- gtex[gtex$`Symbol of blood RNA` %in% symbols,]
dim(gtex)
# [1] 1058955     16

num_eqtls_gene <- c()
for (gene in symbols) {
  num_eqtls_gene <- c(num_eqtls_gene, length(unique(gtex$`rsID of eQTL`[gtex$`Symbol of blood RNA` == gene])))
}
num_traits_gene <- c()
for (gene in symbols) {
  num_traits_gene <- c(num_traits_gene, length(unique(gtex$`GWAS Trait`[gtex$`Symbol of blood RNA` == gene])))
}

median <- read.csv("~/Desktop/work_repo/github/cross_study_analysis/output/median_stability_statistics_all.csv")
median <- median[median$Symbol %in% symbols, ]
# Percentiles!
expression <- median[order(median$Average.Expression),]
expression$percentile <- (1:nrow(expression))/nrow(expression) * 100
colnames(expression)[c(10,11)] <- c("Symbol of blood RNA", "Median Average Normalized Expression Percentile")
expression <- expression[,c("Symbol of blood RNA", "Median Average Normalized Expression Percentile")]
repeatable <- median[order(median$repeatability),]
repeatable$percentile <- (1:nrow(repeatable))/nrow(repeatable) * 100
colnames(repeatable)[c(10,11)] <- c("Symbol of blood RNA", "Median Repeatability Percentile")
repeatable <- repeatable[,c("Symbol of blood RNA", "Median Repeatability Percentile")]
genvar <- median[order(median$gen.variance),]
genvar$percentile <- (1:nrow(genvar))/nrow(genvar) * 100
colnames(genvar)[c(10,11)] <- c("Symbol of blood RNA", "Median Genetic Variance Percentile")
genvar <- genvar[,c("Symbol of blood RNA","Median Genetic Variance Percentile")]

append <- data.frame(`Symbol of blood RNA` = symbols, 
                     `Number of eQTL associations in GTEx` = num_eqtls_gene,
                     `Number of trait associations in GWAS` = num_traits_gene)
colnames(append) <- c("Symbol of blood RNA", "Number of eQTL associations in GTEx", "Number of trait associations in GWAS")
append <- merge(append, expression, by = "Symbol of blood RNA")
append <- merge(append, repeatable, by = "Symbol of blood RNA")
append <- merge(append, genvar, by = "Symbol of blood RNA")

gtex <- merge(gtex, append, by = "Symbol of blood RNA", all.x = T, sort = F)
dim(gtex)
# [1] 1058955      21
for (i in 1:ncol(gtex)) {
  if (class(gtex[,i]) == 'numeric') {
    gtex[,i] <- signif(gtex[,i], digits = 2)
  } else next
}
trait <- read.table("~/Desktop/work_repo/github/cross_study_analysis/output/trait_gene_list.txt")[[1]]
state <- read.table("~/Desktop/work_repo/github/cross_study_analysis/output/state_gene_list.txt")[[1]]
gtex$`Trait gene?` <- gtex$`Symbol of blood RNA` %in% trait
gtex$`State gene?` <- gtex$`Symbol of blood RNA` %in% state
 
library(nanoparquet)
nanoparquet::write_parquet(gtex, file = paste(output_dir, "gtex_variants_and_phenotypes_whole_blood.parquet", sep = ""))

nrow(gtex[is.na(gtex$`rsID of eQTL`),])
# [1] 1647 -- genes without eQTL-GWAS link
length(unique(gtex[!is.na(gtex$`rsID of eQTL`),"Symbol of blood RNA"]))
# [1] 4452
length(unique(gtex$`GWAS Trait`))
# [1] 14952