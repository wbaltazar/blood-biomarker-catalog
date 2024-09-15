## In order for this code to work, all variation files should be run in their respective study_data directories.
## Date: Jul 20 2024

## Load in summary data
library(tidyverse)
library(AnnotationDbi)
library(hgu133plus2.db)
library(EnsDb.Hsapiens.v86)
library(illuminaHumanv3.db)
library(biomaRt)
library(VennDiagram)
library(cowplot)
library(clusterProfiler)

## INPUT: github directory where all files are stored.
input_dir <- "~/Desktop/work_repo/github/"
## OUTPUT: cross-study analysis output folder. We will output CSV files for the genes.
output_dir <- "~/Desktop/work_repo/github/cross_study_analysis/output/"

# load in all statistics ----

## Note that the Symbol columns were all adjusted. When I saved these, if the name was NA, I
## substituted the symbol name with the key used to search. Allowing NAs makes a symbol-oriented
## analysis easier.
gomez <- read_csv(file = paste(input_dir, "Gomez_study_data/output/Gomez_variation.csv", sep = ""))
gomez$Symbol <- mapIds(x = EnsDb.Hsapiens.v86, keys = gomez$...1, column = "SYMBOL", keytype = "GENEID")
gosch <- read_csv(file = paste(input_dir, "Gosch_study_data/output/Gosch_variation.csv", sep = ""))
meaburn1 <- read_csv(file = paste(input_dir, "Meaburn_study_data/output/Meaburn_day1_variation.csv", sep = ""))
meaburn1$Symbol <- mapIds(x = hgu133plus2.db, keys = meaburn1$...1, column = "SYMBOL", keytype = "PROBEID")
meaburn2 <- read_csv(file = paste(input_dir, "Meaburn_study_data/output/Meaburn_day2_variation.csv", sep = ""))
meaburn2$Symbol <- mapIds(x = hgu133plus2.db, keys = meaburn2$...1, column = "SYMBOL", keytype = "PROBEID")
obermoser1 <- read_csv(file = paste(input_dir, "Obermoser_study_data/output/Obermoser_p1v_variation.csv", sep = ""))
obermoser1$Symbol <- mapIds(x = illuminaHumanv3.db, keys = obermoser1$...1, column = "SYMBOL", keytype = "PROBEID")
obermoser2 <- read_csv(file = paste(input_dir, "Obermoser_study_data/output/Obermoser_p1f_variation.csv", sep = ""))
obermoser2$Symbol <- mapIds(x = illuminaHumanv3.db, keys = obermoser2$...1, column = "SYMBOL", keytype = "PROBEID")
obermoser3 <- read_csv(file = paste(input_dir, "Obermoser_study_data/output/Obermoser_p2v_variation.csv", sep = ""))
obermoser3$Symbol <- mapIds(x = illuminaHumanv3.db, keys = obermoser3$...1, column = "SYMBOL", keytype = "PROBEID")
obermoser4 <- read_csv(file = paste(input_dir, "Obermoser_study_data/output/Obermoser_p2f_variation.csv", sep = ""))
obermoser4$Symbol <- mapIds(x = illuminaHumanv3.db, keys = obermoser4$...1, column = "SYMBOL", keytype = "PROBEID")
dusek <- read_csv(file = paste(input_dir, "Dusek_study_data/output/Dusek_variation.csv", sep = ""))
dusek$Symbol <- mapIds(x = hgu133plus2.db, keys = dusek$...1, column = "SYMBOL", keytype = "PROBEID")
rusch <- read_csv(file = paste(input_dir, "Rusch_study_data/output/Rusch_variation.csv", sep = ""))
rusch$Symbol <- mapIds(x = hgu133plus2.db, keys = rusch$...1, column = "SYMBOL", keytype = "PROBEID")
larocca <- read_csv(file = paste(input_dir, "LaRocca_study_data/output/LaRocca_variation.csv", sep = ""))
karlovich1 <- read_csv(file = paste(input_dir, "Karlovich_study_data/output/Karlovich_batch1_variation.csv", sep = ""))
karlovich1$Symbol <- mapIds(x = hgu133plus2.db, keys = karlovich1$...1, column = "SYMBOL", keytype = "PROBEID")
karlovich2 <- read_csv(file = paste(input_dir, "Karlovich_study_data/output/Karlovich_batch2_variation.csv", sep = ""))
karlovich2$Symbol <- mapIds(x = hgu133plus2.db, keys = karlovich2$...1, column = "SYMBOL", keytype = "PROBEID")

data <- list(gomez, gosch, meaburn1, meaburn2, obermoser1, obermoser2, obermoser3, obermoser4,
             dusek, rusch, larocca)
names(data) <- c("gomez", "gosch", "meaburn1", "meaburn2", "obermoser1", "obermoser2", "obermoser3", "obermoser4", 
                 "dusek", "rusch", "larocca")

# filtering ----
## by Total Variation low, average expression above median, gen.variance low ----

rs <- lapply(data, function (x) {
  x %>% 
    dplyr::filter(`Total Variation (SD)` <= quantile(x[["Total Variation (SD)"]])[2]) %>% 
    dplyr::filter(`Average Expression` >= quantile(x[["Average Expression"]])[3]) %>% 
    dplyr::filter(gen.variance <= quantile(x[["gen.variance"]])[2]) %>% 
    dplyr::select(`Symbol`) %>% 
    as.list() %>% unlist() %>% unname()
})

# comparison of houskeeping gene filters ----
## reformat so gene symbols are similar ----
### Replaces all '.' with '-' unless they are the final . proceeded by a number (which come about
# as a consequence of make.names(gene_symbols, unique = T)). Then, removes the . at the end.
rs <- lapply(rs, function(x) {str_replace_all(x, "\\.(?!\\d$)", "-") %>% str_remove_all(pattern = "\\..")})
## use symbols only once ----
rs <- lapply(rs, unique)

names(rs) <- paste(names(rs), "rs", sep = ".")

all_filters <- c(rs)
towrite <- lapply(all_filters, as.data.frame)
writexl::write_xlsx(towrite, path = paste(output_dir, "housekeeping_filter_summary.xlsx", sep = ""))

## How many unique filters are there?
length(all_filters)
# [1] 13

## Get results ----
## Create pool of genes
symbol_find <- function(x) {
  index <- grep("Symbol", names(x))
  return(x[,index])
}
all_symbols <- unique(unname(unlist(lapply(data, symbol_find))))
all_symbols <- str_replace_all(all_symbols, "\\.(?!\\d$)", "-") %>% 
  str_remove_all(pattern = "\\..") %>% 
  unique()
length(unique(all_symbols))
# 28913
all_symbols <- na.omit(all_symbols)
length(all_symbols)
# [1] 28912

## Test run the counting algorithm
housekeeping_results <- data.frame(Symbol = all_symbols, Score = 0, Filters = 0)
all_symbols[1] # [1] "TSPAN6"
sum(unlist(lapply(all_filters, function(x) {all_symbols["RAN3BP"] %in% x})))
# [1] 7
paste(na.omit(names(unlist(lapply(all_filters, function(x) {all_symbols["RAN3BP"] %in% x})))[unlist(lapply(all_filters, function(x) {all_symbols["RAN3BP"] %in% x}))]), collapse = " ")
# "meaburn1.rs meaburn2.rs obermoser2.rs obermoser3.rs obermoser4.rs dusek.rs rusch.rs karlovich1.rs karlovich2.rs"

### Run ----
for (i in 1:nrow(housekeeping_results)) {
  housekeeping_results[i,2] <- sum(unlist(lapply(all_filters, function(x) {housekeeping_results$Symbol[i] %in% x})))
  housekeeping_results[i,3] <- paste(na.omit(names(unlist(lapply(all_filters, function(x) {housekeeping_results$Symbol[i] %in% x})))[unlist(lapply(all_filters, function(x) {housekeeping_results$Symbol[i] %in% x}))]), collapse = " ")

}
head(housekeeping_results)

### Bin which filters the gene passed ----
housekeeping_results[,4:16] <- 0
names(housekeeping_results) <- c(names(housekeeping_results)[1:3], names(data))
head(housekeeping_results,1)
studies <- names(data)
for (study in studies) {
  for (gene in 1:nrow(housekeeping_results)) {
    housekeeping_results[gene,study] <- str_count(housekeeping_results[gene,"Filters"], study)
  }
}
housekeeping_results <- housekeeping_results %>% 
  mutate(Studies = (gomez > 0) + (meaburn1 + meaburn2 > 0) + (gosch > 0) + (obermoser1 + obermoser2 + obermoser3 + obermoser4 > 0) +
           (dusek > 0) + (rusch > 0) + (larocca > 0) + (karlovich1 + karlovich2 > 0))

### Save results -----
write.csv(housekeeping_results, file = paste(output_dir, "housekeeping_scores.csv", sep = ""))

## Table 2 ----
housekeeping_results <- read.csv(paste(output_dir, "housekeeping_scores.csv", sep = ""))
supp_table <- housekeeping_results[order(housekeeping_results$Studies, decreasing = T),]
head(supp_table)
# X  Symbol Score                                                                                     Filters gomez gosch meaburn1
# 379   379  RANBP3     8 gosch.rs meaburn2.rs obermoser3.rs obermoser4.rs dusek.rs rusch.rs larocca.rs karlovich1.rs     0     1        0
# 5572 5572   PARP6     7               gomez.rs gosch.rs meaburn1.rs obermoser4.rs rusch.rs larocca.rs karlovich1.rs     1     1        1
# 212   212  LYPLA2     7              gomez.rs gosch.rs meaburn1.rs meaburn2.rs obermoser4.rs rusch.rs karlovich1.rs     1     1        1
# 223   223  RNF216     6                             gomez.rs meaburn2.rs dusek.rs rusch.rs larocca.rs karlovich2.rs     1     0        0
# 1131 1131   REXO1     6                        gosch.rs meaburn2.rs obermoser4.rs rusch.rs larocca.rs karlovich1.rs     0     1        0
# 1724 1724 POLDIP3     6                             gomez.rs gosch.rs meaburn1.rs obermoser3.rs rusch.rs larocca.rs     1     1        1
# meaburn2 obermoser1 obermoser2 obermoser3 obermoser4 dusek rusch larocca karlovich1 karlovich2 Studies
# 379         1          0          0          1          1     1     1       1          1          0       7
# 5572        0          0          0          0          1     0     1       1          1          0       7
# 212         1          0          0          0          1     0     1       0          1          0       6
# 223         1          0          0          0          0     1     1       1          0          1       6
# 1131        1          0          0          0          1     0     1       1          1          0       6
# 1724        0          0          0          1          0     0     1       1          0          0       6
dim(housekeeping_results[housekeeping_results$Studies >= 6,])
# [1] 17 18
supp_table <- supp_table[supp_table$Studies >= 6, "Symbol"]
library(biomaRt)
mart <- useDataset(dataset = "hsapiens_gene_ensembl", mart = useMart("ENSEMBL_MART_ENSEMBL"))
supp_table <- getBM(attributes = c("hgnc_symbol","description"),
                    filters = "hgnc_symbol",
                    values = supp_table,
                    mart = mart)
supp_table
#    hgnc_symbol                                                                                   description
# 1        DHDDS              dehydrodolichyl diphosphate synthase subunit [Source:HGNC Symbol;Acc:HGNC:20603]
# 2       EFTUD2      elongation factor Tu GTP binding domain containing 2 [Source:HGNC Symbol;Acc:HGNC:30858]
# 3        GPAA1           glycosylphosphatidylinositol anchor attachment 1 [Source:HGNC Symbol;Acc:HGNC:4446]
# 4         ING4                       inhibitor of growth family member 4 [Source:HGNC Symbol;Acc:HGNC:19423]
# 5       LYPLA2                                        lysophospholipase 2 [Source:HGNC Symbol;Acc:HGNC:6738]
# 6         MBD1                        methyl-CpG binding domain protein 1 [Source:HGNC Symbol;Acc:HGNC:6916]
# 7       MRPS25                       mitochondrial ribosomal protein S25 [Source:HGNC Symbol;Acc:HGNC:14511]
# 8     NFATC2IP nuclear factor of activated T cells 2 interacting protein [Source:HGNC Symbol;Acc:HGNC:25906]
# 9        PARP6               poly(ADP-ribose) polymerase family member 6 [Source:HGNC Symbol;Acc:HGNC:26921]
# 10     POLDIP3                DNA polymerase delta interacting protein 3 [Source:HGNC Symbol;Acc:HGNC:23782]
# 11      RANBP3                                      RAN binding protein 3 [Source:HGNC Symbol;Acc:HGNC:9850]
# 12       REXO1                                 RNA exonuclease 1 homolog [Source:HGNC Symbol;Acc:HGNC:24616]
# 13       RHOT2                              ras homolog family member T2 [Source:HGNC Symbol;Acc:HGNC:21169]
# 14      RNF216                                   ring finger protein 216 [Source:HGNC Symbol;Acc:HGNC:21698]
# 15     SH3KBP1            SH3 domain containing kinase binding protein 1 [Source:HGNC Symbol;Acc:HGNC:13867]
# 16      STRADA                               STE20 related adaptor alpha [Source:HGNC Symbol;Acc:HGNC:30172]
# 17      WRNIP1                        WRN helicase interacting protein 1 [Source:HGNC Symbol;Acc:HGNC:20876]
write.csv(supp_table, paste(output_dir, "housekeeping_6_or_more_studies.csv", sep = ""))