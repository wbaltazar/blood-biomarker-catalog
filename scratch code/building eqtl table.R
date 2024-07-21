## Building a tool for querying your rsIDs!
setwd("~/Desktop/work_repo/")

library(tidyverse)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
gtex <- read_tsv("./GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz")
gtex$gene_id <- gsub(pattern = "(.*)(\\..*)", replacement = "\\1", x = gtex$gene_id)
gtex$Symbol <- mapIds(EnsDb.Hsapiens.v86, keys = gtex$gene_id, column = "SYMBOL", keytype = "GENEID")

variant_ids <- gtex$variant_id
length(variant_ids)
# [1] 2414653

# library(data.table)
# id_table <- fread(input = "./GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz",
#       select = c("variant_id", "rs_id_dbSNP151_GRCh38p7"))
# 
# id_table <- id_table %>% 
#   dplyr::filter(variant_id %in% variant_ids)
# 
# # Save the results for later!!
# write.table(id_table, "./filtered_rsIDs.txt")
id_table <- read.table("./filtered_rsIDs.txt")
head(id_table)
dim(id_table)
length(unique(gtex$variant_id))
# [1] 1277338

gtex <- merge(gtex, id_table, by = "variant_id")
dim(gtex)
# [1] 2414653      14
head(gtex)
# variant_id         gene_id tss_distance ma_samples ma_count       maf pval_nominal     slope  slope_se
# 1 chr1_100000723_G_A_b38 ENSG00000156876      -132232        153      172 0.1283580  1.55009e-04  0.147698 0.0387947
# 2   chr1_1000018_G_A_b38 ENSG00000187608        -1120         40       40 0.0298507  9.47022e-14 -0.501117 0.0657039
# 3 chr1_100002416_C_T_b38 ENSG00000156876      -130539         64       67 0.0500000  1.54726e-08  0.334512 0.0583289
# 4 chr1_100003060_C_T_b38 ENSG00000122477      -175213         62       67 0.0500000  4.36290e-06  0.401704 0.0866499
# 5 chr1_100003083_G_T_b38 ENSG00000156876      -129872         62       65 0.0485075  5.92769e-08  0.324030 0.0590180
# 6  chr1_10000449_G_A_b38 ENSG00000162444         3243         26       27 0.0201493  2.91247e-05 -0.261548 0.0620871
# pval_nominal_threshold min_pval_nominal   pval_beta Symbol rsID rs_id_dbSNP151_GRCh38p7
# 1            0.000246184      1.54726e-08 5.63557e-05  SASS6    1              rs11166389
# 2            0.000134488      1.41503e-28 2.86633e-23  ISG15    1             rs146254088
# 3            0.000246184      1.54726e-08 5.63557e-05  SASS6    1              rs12128170
# 4            0.000190134      7.36133e-09 2.28188e-05 LRRC39    1                rs547241
# 5            0.000246184      1.54726e-08 5.63557e-05  SASS6    1              rs11801439
# 6            0.000104450      1.61451e-06 5.01701e-03   RBP7    1             rs141577932

# Nice!

# Next, get the associated traits for the rsIDs.
associations <- read_tsv("./gwas_catalog_v1.0-associations_e111_r2024-04-22.tsv")
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
  dplyr::select(`DISEASE/TRAIT`, `REPORTED GENE(S)`, MAPPED_GENE, SNPS, CONTEXT, INTERGENIC, `OR or BETA`,
                `P-VALUE`, PVALUE_MLOG,)
dim(associations)
# [1] 607094      9
head(associations)
# # A tibble: 6 × 9
# `DISEASE/TRAIT`       `REPORTED GENE(S)` MAPPED_GENE SNPS  CONTEXT INTERGENIC `OR or BETA` `P-VALUE` PVALUE_MLOG
# <chr>                 <chr>              <chr>       <chr> <chr>        <dbl>        <dbl>     <dbl>       <dbl>
#   1 Mean corpuscular vol… NR                 PPIF - ZCC… rs10… interg…          1       0.0141     3e-12        11.5
# 2 Mean corpuscular vol… NR                 PTEN - MED… rs80… interg…          1       0.0192     4e-23        22.4
# 3 Mean corpuscular vol… NR                 FAS         rs38… non_co…          0       0.0153     3e-17        16.5
# 4 Mean corpuscular vol… NR                 TBC1D12     rs11… intron…          0       0.0138     5e-14        13.3
# 5 Mean corpuscular vol… NR                 NKX2-3 - S… rs10… interg…          1       0.0300     1e-47        47  
# 6 Mean corpuscular vol… NR                 NT5C2       rs12… intron…          0       0.0290     4e-48        47.4

colnames(associations)[4] <- "rs_id_dbSNP151_GRCh38p7"
gtex <- merge(gtex, associations, by = "rs_id_dbSNP151_GRCh38p7", all.x = TRUE)
dim(gtex)
# [1] 2653785      22
colnames(gtex)
# [1] "rs_id_dbSNP151_GRCh38p7" "variant_id"              "gene_id"                 "tss_distance"           
# [5] "ma_samples"              "ma_count"                "maf"                     "pval_nominal"           
# [9] "slope"                   "slope_se"                "pval_nominal_threshold"  "min_pval_nominal"       
# [13] "pval_beta"               "Symbol"                  "DISEASE/TRAIT"           "REPORTED GENE(S)"       
# [17] "MAPPED_GENE"             "CONTEXT"                 "INTERGENIC"              "OR or BETA"             
# [21] "P-VALUE"                 "PVALUE_MLOG"  
gtex <- gtex[,c("Symbol", "gene_id", "rs_id_dbSNP151_GRCh38p7", "variant_id", "maf", "slope", 
                "pval_nominal", "pval_beta",
                "DISEASE/TRAIT", "MAPPED_GENE", "REPORTED GENE(S)", "OR or BETA", "P-VALUE", "CONTEXT",
                "INTERGENIC")] # 15 columns
colnames(gtex) <- c("Symbol of blood RNA", "Ensembl ID", "rsID of eQTL", "GTEx variant ID","MAF in GTEx WB",
                    "eQTL slope", "eQTL nominal p_val", "eQTL beta p_value",
                    "GWAS Trait", "GWAS Catalog Mapping", "Reported Mapping", "OR or BETA", "GWAS p_value",
                    "Genomic Context","Intergenic?")
stable <- read.csv("./Box organization/1results/RShiny-application/data/stable_gene_filter_count.csv")
dynamic <- read.csv("./Box organization/1results/RShiny-application/data/deg_counts.csv")
house <- read.csv("./Box organization/1results/RShiny-application/data/housekeeping_scores.csv")
stable <- stable %>% dplyr::select(Symbol, Total)
colnames(stable) <- c("Symbol of blood RNA","Stable-polymorphic score")
dynamic <- dynamic %>% dplyr::select(Symbol, Count)
colnames(dynamic) <- c("Symbol of blood RNA","Flexibility score")
house <- house %>% dplyr::select(Symbol, Total)
colnames(house) <- c("Symbol of blood RNA","Housekeeping score")

gtex <- merge(gtex, stable, all.x = T, by = "Symbol of blood RNA", sort = F)
gtex <- merge(gtex, dynamic, all.x = T, by = "Symbol of blood RNA", sort = F)
gtex <- merge(gtex, house, all.x = T, by = "Symbol of blood RNA", sort = F)
dim(gtex)
# [1] 2653785      18
symbols <- read.table("./Box organization/1results/Symbols_common_to_all_datasets.txt")
symbols <- symbols$x
gtex <- gtex[gtex$`Symbol of blood RNA` %in% symbols,]
dim(gtex)
# [1] 1352617      18

write_rds(gtex, "./Box organization/1results/RShiny-application/data/eqtls/gtex_variants_and_phenotypes_whole_blood_28062024.rds")
