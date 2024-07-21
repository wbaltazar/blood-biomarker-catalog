library(gwasrapidd)

## Compare the GWAS Catalog API to gwasrapidd

query <- paste("https://www.ebi.ac.uk/gwas/rest/api/associations?gene=", "MYOM2", sep = "")
response <- GET(query)
data <- fromJSON(content(response, "text"), flatten = TRUE)
data.frame(data$`_embedded`$associations)
dim(data.frame(data$`_embedded`$associations))
GET("https://www.ebi.ac.uk/gwas/rest/api/associations/14347/efoTraits")

rapid <- gwasrapidd::get_variants(gene_name = "MYOM2")
genomic_contexts <- rapid@genomic_contexts %>% 
  filter(gene_name == "MYOM2")
variant_info <- rapid@variants %>% 
  filter(variant_id %in% genomic_contexts$variant_id)
var_query <- variant_info$variant_id

rapid2 <- gwasrapidd::get_associations(variant_id = var_query)
colnames(rapid2@associations)
asc_query <- rapid2@associations$association_id

rapid3 <- gwasrapidd::get_traits(association_id = asc_query)
rapid3

gwas_query <- function(gene) {
  # Initial query
  cat("Obtaining variant IDs...\n")
  rapid <- gwasrapidd::get_variants(gene_name = gene)
  genomic_contexts <- rapid@genomic_contexts %>% 
    filter(gene_name == gene)
  variant_info <- rapid@variants %>% 
    filter(variant_id %in% genomic_contexts$variant_id)
  var_query <- variant_info$variant_id
  whatvars <- unique(rapid@variants$variant_id)
  
  # Obtain associations
  cat("Obtaining associations...\n")
  rapid2 <- gwasrapidd::get_associations(variant_id = var_query)
  colnames(rapid2@associations)
  whatass <- unique(rapid2@risk_alleles$variant_id)
  
  # Remove extraneous IDs and create next query
  not_conserved_ids <- unique(c(whatass[which(!(whatass == whatvars))], whatvars[which(!(whatass == whatvars))]))
  not_conserved_ass <- rapid2@risk_alleles[rapid2@risk_alleles$variant_id %in% not_conserved_ids,1]
  asc_query <- rapid2@associations$association_id[which(!(rapid2@associations$association_id %in% not_conserved_ass))]
  
  # Obtain traits
  cat("Obtaining traits...\n")
  rapid3 <- gwasrapidd::get_traits(association_id = asc_query)
  
  # Constructing table
  results <- list(Variants = rapid, Associations = rapid2, Traits = rapid3)
  return(results)
}

viz <- gwas_query("MYOM2")

## FInal check: downloaded associations from the catalog
gca <- read_tsv("~/Downloads/gwas_catalog_v1.0-associations_e111_r2024-04-22.tsv")
colnames(gca)
# [1] "DATE ADDED TO CATALOG"      "PUBMEDID"                   "FIRST AUTHOR"               "DATE"                      
# [5] "JOURNAL"                    "LINK"                       "STUDY"                      "DISEASE/TRAIT"             
# [9] "INITIAL SAMPLE SIZE"        "REPLICATION SAMPLE SIZE"    "REGION"                     "CHR_ID"                    
# [13] "CHR_POS"                    "REPORTED GENE(S)"           "MAPPED_GENE"                "UPSTREAM_GENE_ID"          
# [17] "DOWNSTREAM_GENE_ID"         "SNP_GENE_IDS"               "UPSTREAM_GENE_DISTANCE"     "DOWNSTREAM_GENE_DISTANCE"  
# [21] "STRONGEST SNP-RISK ALLELE"  "SNPS"                       "MERGED"                     "SNP_ID_CURRENT"            
# [25] "CONTEXT"                    "INTERGENIC"                 "RISK ALLELE FREQUENCY"      "P-VALUE"                   
# [29] "PVALUE_MLOG"                "P-VALUE (TEXT)"             "OR or BETA"                 "95% CI (TEXT)"             
# [33] "PLATFORM [SNPS PASSING QC]" "CNV"        

myom2 <- gca %>% filter(MAPPED_GENE == "MYOM2")
View(myom2)

gtex <- read_tsv("./data/eqtls/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz")
gtex$gene_id <- gsub(pattern = "(.*)(\\..*)", replacement = "\\1", x = gtex$gene_id)
gtex$Symbol <- mapIds(EnsDb.Hsapiens.v86, keys = gtex$gene_id, column = "SYMBOL", keytype = "GENEID")
gmyom2 <- gtex %>% filter(Symbol == "MYOM2")
View(gmyom2)
## Curious if these data are located within the other?
myom2_eqtls_chrloc <- str_extract(gmyom2$variant_id, "(?<=_)[0-9]+(?=_)")
gca_in_gtex <- gca %>% filter(CHR_POS %in% myom2_eqtls_chrloc, CHR_ID == 8)
View(gca_in_gtex)

## And what do we see in gwasrapidd?
my_results <- gwas_query("MYOM2")
View(my_results)
