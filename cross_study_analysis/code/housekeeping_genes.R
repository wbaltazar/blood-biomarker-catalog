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
# [1] 11

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
# "meaburn1.rs meaburn2.rs obermoser2.rs obermoser3.rs obermoser4.rs dusek.rs rusch.rs"

### Run ----
for (i in 1:nrow(housekeeping_results)) {
  housekeeping_results[i,2] <- sum(unlist(lapply(all_filters, function(x) {housekeeping_results$Symbol[i] %in% x})))
  housekeeping_results[i,3] <- paste(na.omit(names(unlist(lapply(all_filters, function(x) {housekeeping_results$Symbol[i] %in% x})))[unlist(lapply(all_filters, function(x) {housekeeping_results$Symbol[i] %in% x}))]), collapse = " ")

}
head(housekeeping_results)

### Bin which filters the gene passed ----
housekeeping_results[,4:14] <- 0
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
           (dusek > 0) + (rusch > 0) + (larocca > 0))

### Save results -----
write.csv(housekeeping_results, file = paste(output_dir, "housekeeping_scores.csv", sep = ""))
print(head(housekeeping_results[order(housekeeping_results$Studies, decreasing = T),c(1,2,length(names(housekeeping_results)))], 10))
#       Symbols Score Studies
# 379    RANBP3     7       6
# 1724  POLDIP3     6       6
# 4167   WRNIP1     6       6
# 5572    PARP6     6       6
# 6504  SH3KBP1     7       6
# 9767 NFATC2IP     6       6
# 212    LYPLA2     6       5
# 223    RNF216     5       5
# 569  TBC1D22A     6       5
# 791      DHX8     5       5