## Grouping genes according to their differential expression across studies
## In order for this code to work, all variation files should be run in their respective study_data directories.
## Date: Jan 17 2025

## Load in summary data
library(tidyverse)
library(AnnotationDbi)
library(hgu133plus2.db)
library(EnsDb.Hsapiens.v86)
library(illuminaHumanv3.db)
library(biomaRt)
library(VennDiagram)
library(cowplot)

## INPUT: github directory where all files are stored.
input_dir <- "~/Desktop/work_repo/github/"
## OUTPUT: cross-study analysis output folder. We will output CSV files for the genes.
output_dir <- "~/Desktop/work_repo/github/cross_study_analysis/output/"

# load in all statistics ----

## Adjust the symbols columns to facilitate cross-study analysis
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
             dusek, rusch, larocca, karlovich1, karlovich2)
names(data) <- c("gomez", "gosch", "meaburn1", "meaburn2", "obermoser1", "obermoser2", "obermoser3", "obermoser4", 
                 "dusek", "rusch", "larocca", "karlovich1", "karlovich2")

# THRESHOLDS ----
## by SD ----

sd <- lapply(data, function (x) {
  x %>% 
    dplyr::filter(`Within Variation (SD)` <= quantile(x[["Within Variation (SD)"]], p = 0.25)) %>% 
    dplyr::filter(`Total Variation (SD)` >= quantile(x[["Total Variation (SD)"]], p = 0.75)) %>% 
    dplyr::select(`Symbol`) %>% 
    as.list() %>% unlist() %>% unname()
})

## by varianceParition ----

vp <- lapply(data, function (x) {
  x %>% 
    dplyr::filter(time <= quantile(x[["time"]], p = 0.25)) %>% 
    dplyr::filter(subject >= quantile(x[["subject"]], p = 0.75)) %>% 
    dplyr::select(`Symbol`) %>% 
    as.list() %>% unlist() %>% unname()
})

## by repeatability & gen.variance ----

rp <- lapply(data, function (x) {
  x %>% 
    dplyr::filter(repeatability >= quantile(x[["repeatability"]], p = 0.75)) %>% 
    dplyr::filter(gen.variance >= quantile(x[["gen.variance"]], p = 0.75)) %>% 
    dplyr::select(`Symbol`) %>% 
    as.list() %>% unlist() %>% unname()
})

# Trait genes across studies ----
## reformat so gene symbols are similar ----
### Replaces all '.' with '-' unless they are the final . proceeded by a number (which come about
# as a consequence of make.names(gene_symbols, unique = T)). Then, removes the . at the end.
sd <- lapply(sd, function(x) {str_replace_all(x, "\\.(?!\\d$)", "-") %>% str_remove_all(pattern = "\\..")})
vp <- lapply(vp, function(x) {str_replace_all(x, "\\.(?!\\d$)", "-") %>% str_remove_all(pattern = "\\..")})
rp <- lapply(rp, function(x) {str_replace_all(x, "\\.(?!\\d$)", "-") %>% str_remove_all(pattern = "\\..")})
## use symbols only once ----
sd <- lapply(sd, unique)
vp <- lapply(vp, unique)
rp <- lapply(rp, unique)

names(sd) <- paste(names(sd), "sd", sep = ".")
names(vp) <- paste(names(vp), "vp", sep = ".")
names(rp) <- paste(names(rp), "rp", sep = ".")

all_filters <- c(sd, vp, rp)
towrite <- lapply(all_filters, as.data.frame)
writexl::write_xlsx(towrite, path = paste(output_dir, "trait_thresholds.xlsx", sep = ""))

## How many unique filters are there?
length(all_filters)
# [1] 39

## Get results ----
## Create pool of all genes
symbol_find <- function(x) {
  index <- grep("Symbol", names(x))
  return(x[,index])
}
all_symbols <- unname(unlist(lapply(data, symbol_find)))
# Use dash for notation instead of periods, EXCEPT when it denotes a version number. Then, remove all suspected version numbers.
all_symbols <- all_symbols %>% str_replace_all("\\.(?!\\d$)", "-") 
length(unlist(unique(all_symbols)))
# [1] 24100
all_symbols <- all_symbols %>% str_remove_all(pattern = "\\..") %>% unique()
length(unlist(all_symbols))
# [1] 23883

## Test run the counting algorithm
stable_results <- data.frame(Symbol = all_symbols, Score = 0, Filters = 0)
all_symbols[1] # [1] "TSPAN6"
sum(unlist(lapply(all_filters, function(x) {all_symbols[1] %in% x})))
# [1] 1
paste(na.omit(names(unlist(lapply(all_filters, function(x) {all_symbols[1] %in% x})))[unlist(lapply(all_filters, function(x) {all_symbols[1] %in% x}))]), collapse = " ")
# "gosch.rp"

### Run ----
for (i in 1:nrow(stable_results)) {
  stable_results[i,2] <- sum(unlist(lapply(all_filters, function(x) {stable_results$Symbol[i] %in% x})))
  stable_results[i,3] <- paste(na.omit(names(unlist(lapply(all_filters, function(x) {stable_results$Symbol[i] %in% x})))[unlist(lapply(all_filters, function(x) {stable_results$Symbol[i] %in% x}))]), collapse = " ")
}

### Studies the gene passed ----
studs <- c("gomez", "meaburn", "gosch", "obermoser", "dusek", "rusch", "larocca", "karlovich")
stable_results[,4:11] <- 0
names(stable_results)[4:11] <- studs
for (i in 1:nrow(stable_results)) {
  for (j in 4:11) {
    stable_results[i,j] <- sum(grepl(names(stable_results)[j], unlist(str_split(stable_results[i,3], pattern = " "))))
  }
}
stable_results <- stable_results %>% 
  mutate(Study_counts = (gomez > 0) + (meaburn > 0) + (gosch > 0) + (obermoser > 0) +
           (dusek > 0) + (rusch > 0) + (larocca > 0) + (karlovich > 0))
head(stable_results)

symbol_not_found <- function(x) {
  index <- grep("Symbol", names(x))
  ref <- x[[index]]
  ref <- ref %>% str_replace_all("\\.(?!\\d$)", "-") %>% str_remove_all(pattern = "\\..") %>% unique()
  return(setdiff(all_symbols, ref))
}

## If a symbol is not found in a dataset, replace its value with NA
na_index <- lapply(data, symbol_not_found)
names(na_index) <- paste(names(na_index), ".na", sep = "")
for_merge <- data.frame(Symbol = all_symbols, Studies_not_found_in = 0)
i = 1
for_merge$Symbol[i]
# [1] "TSPAN6"
paste(na.omit(names(unlist(lapply(na_index, function(x) {for_merge$Symbol[i] %in% x})))[unlist(lapply(na_index, function(x) {for_merge$Symbol[i] %in% x}))]), collapse = " ")
# [1] "meaburn1.na meaburn2.na dusek.na rusch.na larocca.na karlovich1.na karlovich2.na"
for (i in 1:nrow(for_merge)) {
  for_merge[i,2] <- paste(na.omit(names(unlist(lapply(na_index, function(x) {for_merge$Symbol[i] %in% x})))[unlist(lapply(na_index, function(x) {for_merge$Symbol[i] %in% x}))]), collapse = " ")
}
stable_results <- merge(stable_results, for_merge, by = "Symbol", sort = F)
dim(stable_results)
# [1] 23883    13
length(unique(stable_results$Symbol))
# [1] 23883

### Save results -----
write.csv(stable_results, file = paste(output_dir, "trait_scores.csv", sep = ""))

## Median statistic function
summarize_statistics <- function(gene) {
  stats <- vector(mode = "numeric", length = 8)
  names(stats) <- c("Within Variation (SD)", "Total Variation (SD)", "Rs", "subject", "time", "repeatability", "gen.variance", "Average Expression")
  n <- 11
  for (i in data) {
    if (nrow(i[which(i$Symbol == gene),]) > 0) {
      stats <- stats + colMeans(i[which(i$Symbol == gene), c("Within Variation (SD)", "Total Variation (SD)", "Rs", "subject", "time", "repeatability", "gen.variance", "Average Expression")])
    } else {
      n <- n - 1
    }
  }
  return((unlist(stats / n)))
}
# Test
summarize_statistics("ERAP2")
# Within Variation (SD)  Total Variation (SD)                    Rs               subject                  time 
#           0.279485086           0.987803585           0.406376948           0.824516204           0.001276759 
# repeatability          gen.variance    Average Expression 
#   0.917736021           0.862993756           9.191705602  

## Initialize data frame
gene_summaries <- data.frame(`Within Variation (SD)` = 0, 
                             `Total Variation (SD)` = 0,
                             Rs = 0,
                             subject = 0,
                             time = 0,
                             repeatability = 0,
                             gen.variance = 0,
                             `Average Expression` = 0)
## Summarize each gene from the list of all genes ----
for (i in 1:length(stable_results$Symbol)) {
  gene_summaries[i,] <- summarize_statistics(stable_results$Symbol[i])
}
head(gene_summaries)
gene_summaries$Symbol <- stable_results$Symbol

## Save results ----
write.csv(gene_summaries, file = paste(output_dir, "median_stability_statistics_all.csv", sep = ""))

# COMMON SYMBOLS ----
# Genes found in all datasets
common_symbols <- for_merge[which(for_merge$Studies_not_found_in == ""),"Symbol"]
length(common_symbols)
# [1] 6099
write.table(common_symbols, file = paste(output_dir, "common_symbols6099.txt", sep = ""))

## Supplementary Table 1 ----
stable_results <- read.csv(file = paste(output_dir, "trait_scores.csv", sep = ""))
stable_results <- stable_results %>% dplyr::filter(Study_counts == 8)
dim(stable_results)
# [1] 112 14
supp_table <- stable_results$Symbol
library(biomaRt)
mart <- useDataset(dataset = "hsapiens_gene_ensembl", mart = useMart("ENSEMBL_MART_ENSEMBL"))
supp_table <- getBM(attributes = c("hgnc_symbol","description"),
                    filters = "hgnc_symbol",
                    values = supp_table,
                    mart = mart)
write.csv(supp_table, file = paste(output_dir, "supplementary_table_1.csv", sep = ""))