## Grouping genes according to their differential expression across studies
## In order for this code to work, all limma_QC files should be run in their respective study_data directories.
## Date: Jul 19 2024

## INPUT: github directory where all files are stored.
input_dir <- "~/Desktop/work_repo/github/"
## OUTPUT: cross-study analysis output folder. We will output CSV files for the genes.
output_dir <- "~/Desktop/work_repo/github/cross_study_analysis/output/"

# Load libraries ----
library(tidyverse)
library(AnnotationDbi)
library(hgu133plus2.db)
library(EnsDb.Hsapiens.v86)
library(clusterProfiler)
library(org.Hs.eg.db)

# Load in all DE time analyses ----
gomez <- read.csv(file = paste(input_dir, "Gomez_study_data/output/limma_time_main.csv", sep = ""))
names(gomez)[1] <- "EnsemblID"

meaburn1 <- read.csv(file = paste(input_dir, "Meaburn_study_data/output/day1_Time_limma.csv", sep = ""))
names(meaburn1)[1] <- "ProbeID"
meaburn2 <- read.csv(file = paste(input_dir, "Meaburn_study_data/output/day2_Time_limma.csv", sep = ""))
names(meaburn2)[1] <- "ProbeID"

gosch <- list()
files <- list.files(paste(input_dir, "Gosch_study_data/output/", sep = ""), full.names = TRUE)
files <- files[grep("limma", files, perl = T)]
files <- files[-3] # Removes the limma_F.csv file
files <- files[-grep("Male", files)] # Removes interaction term files
names <- str_extract(files, "limma_([^.]*)")
names <- str_replace(names, "limma_", "")
for (i in 1:length(files)) {
  gosch[[i]] <- read.csv(files[i])
  names(gosch)[i] <- names[i]
}

obermoser1 <- list()
files <- list.files(paste(input_dir, "Obermoser_study_data/output/training_vein_limma/", sep = ""), full.names = TRUE)
files <- files[grep("ALL", files, perl = T)]
names <- str_extract(files, "Base([^.]*)")
for (i in 1:length(files)) {
  obermoser1[[i]] <- read.csv(files[i])
  names(obermoser1)[i] <- names[i]
}

obermoser2 <- list()
files <- list.files(paste(input_dir, "Obermoser_study_data/output/test_finger_limma/", sep = ""), full.names = TRUE)
files <- files[grep("All", files, perl = T)]
names <- str_extract(files, "Base([^.]*)")
for (i in 1:length(files)) {
  obermoser2[[i]] <- read.csv(files[i])
  names(obermoser2)[i] <- names[i]
}

obermoser3 <- list()
files <- list.files(paste(input_dir, "Obermoser_study_data/output/test_vein_limma/", sep = ""), full.names = TRUE)
files <- files[grep("All", files, perl = T)]
names <- str_extract(files, "Base([^.]*)")
for (i in 1:length(files)) {
  obermoser3[[i]] <- read.csv(files[i])
  names(obermoser3)[i] <- names[i]
}

obermoser4 <- list()
files <- list.files(paste(input_dir, "Obermoser_study_data/output/training_finger_limma/", sep = ""), full.names = TRUE)
files <- files[grep("All", files, perl = T)]
names <- str_extract(files, "Base([^.]*)")
for (i in 1:length(files)) {
  obermoser4[[i]] <- read.csv(files[i])
  names(obermoser4)[i] <- names[i]
}

dusek <- read.csv(file = paste(input_dir, "Dusek_study_data/output/limma_time.csv", sep = ""))
names(dusek) <- str_remove(names(dusek), "time\\.")
names(dusek)[c(1,8)] <- c("ProbeID", "Symbol")
dusek <- dusek[,-9]

rusch <- read.csv(file = paste(input_dir, "Rusch_study_data/output/limma_Time.csv", sep = ""))
names(rusch)[1] <- "ProbeID"
rusch$Symbol <- mapIds(hgu133plus2.db, rusch$ProbeID, column = "SYMBOL", keytype = "PROBEID")

larocca <- read.csv(file = paste(input_dir, "LaRocca_study_data/output/limma_Time.csv", sep = ""))
names(larocca)[1] <- "Symbol"

karlovich1 <- list()
karlovich1[[1]] <- read.csv(file = paste(input_dir, "Karlovich_study_data/output/batch1_limma/limma_Day14.csv", sep = ""))
names(karlovich1[[1]])[1] <- "ProbeID"
karlovich1[[2]] <- read.csv(file = paste(input_dir, "Karlovich_study_data/output/batch1_limma/limma_Day28.csv", sep = ""))
names(karlovich1[[2]])[1] <- "ProbeID"
names(karlovich1) <- c("karlovich_14d", "karlovich_28d")

karlovich2 <- read.csv(file = paste(input_dir, "Karlovich_study_data/output/batch2_limma/limma_Time.csv", sep = ""))
names(karlovich2)[1] <- "ProbeID"

data <- list(gomez, meaburn1, meaburn2, gosch, obermoser1, obermoser2, obermoser3, obermoser4, dusek, rusch, larocca,
             karlovich1, karlovich2)

# Function to recursively flatten a nested list
flatten_list <- function(lst) {
  flat_list <- list()
  
  # Recursive helper function
  flatten_helper <- function(x) {
    if (!is.data.frame(x)) {
      for (item in x) {
        flatten_helper(item)
      }
    } else {
      flat_list <<- append(flat_list, list(x))
    }
  }
  
  # Apply the helper function to the input list
  flatten_helper(lst)
  
  return(flat_list)
}

data <- flatten_list(data)
names(data) <- c("gomez", "meaburn1", "meaburn2", paste("gosch",names(gosch),sep=""), paste("obermoser1",names(obermoser1),sep=""), 
                 paste("obermoser2",names(obermoser2),sep=""), paste("obermoser3",names(obermoser3),sep=""), 
                 paste("obermoser4",names(obermoser4),sep=""), "dusek", "rusch", "larocca", paste(names(karlovich1)),
                 "karlovich2")
length(data)
# [1] 43

## The gosch data does not have a "symbol" column. Let's add it in.
names(data$goschEighteenHrMain)[1] <- "Symbol"
names(data$goschFifteenHrMain)[1] <- "Symbol"
names(data$goschNineHrMain)[1] <- "Symbol"
names(data$goschSixHrMain)[1] <- "Symbol"
names(data$goschThreeHrMain)[1] <- "Symbol"
names(data$goschTwelveHrMain)[1] <- "Symbol"
names(data$goschTwentyoneHrMain)[1] <- "Symbol"

# Extract all comparisons with p < 0.05 ----
p_filter <- function(x) {
  x %>%  dplyr::filter(P.Value < 0.05)
}

pfilt <- lapply(data, p_filter)

adjp_filter <- function(x) {
  x %>%  dplyr::filter(adj.P.Val < 0.05)
}

adjpfilt <- lapply(data, adjp_filter)

symbol_find <- function(x) {
  index <- grep("Symbol", names(x))
  return(x[,index])
}

## To have a common gene naming system, here are the str replacement tools we use:
## str_replace_all(x, "\\.(?!\\d$)", "-") replace dots with dashes only if they are not followed by a digit.
## RATIONALE: Genes like HLA.DRB4 should be changed to HLA-DRB4, but not genes like GSTM1.1. A name like
## GSTM1.1 appears in probe-based platforms where a transcript is targeted by multiple different probes.
## Certain transcripts (like those identified in RNA-seq) will lose their specificity using this method.
## When we map the gene symbol to the Probe ID, we use make.names with unique = T, which appends a ".X" for
## the Xth repeat name of that symbol. To correct for that, we use this function:
## str_remove_all(pattern = "\\..")

all_symbols <- lapply(data, symbol_find)
all_symbols <- lapply(all_symbols, function(x) {str_replace_all(x, "\\.(?!\\d$)", "-") %>% str_remove_all(pattern = "\\..")})
all_symbols <- lapply(all_symbols, unique)

pfilt_symbols <- lapply(pfilt, symbol_find)
pfilt_symbols <- lapply(pfilt_symbols, function(x) {str_replace_all(x, "\\.(?!\\d$)", "-") %>% str_remove_all(pattern = "\\..")})
pfilt_symbols <- lapply(pfilt_symbols, unique)

adjpfilt_symbols <- lapply(adjpfilt, symbol_find)
adjpfilt_symbols <- lapply(adjpfilt_symbols, function(x) {str_replace_all(x, "\\.(?!\\d$)", "-") %>% str_remove_all(pattern = "\\..")})
adjpfilt_symbols <- lapply(adjpfilt_symbols, unique)

all_symbols <- unique(unname(unlist(all_symbols)))
length(all_symbols)
# [1] 28974

## Test concept
all_symbols[1] # [1] "NEK8"
sum(unlist(lapply(pfilt_symbols, function(x) {all_symbols[1] %in% x})))
# [1] 9
names(unlist(lapply(pfilt_symbols, function(x) {all_symbols[1] %in% x})))[unlist(lapply(pfilt_symbols, function(x) {all_symbols[1] %in% x}))]
# [1] "gomez"                            "goschFifteenHrMain"               "goschNineHrMain"                 
# [4] "goschSixHrMain"                   "goschTwelveHrMain"                "obermoser4BaselineVSfifteenAll"  
# [7] "obermoser4BaselineVSthirtySixAll" "obermoser4BaselineVStwelveAll"    "rusch"
deg_results <- data.frame(Symbol = all_symbols, `P_counts` = 0, `P_study_names` = 0, `AdjP_counts` = 0,
                          `AdjP_study_names` = 0)
head(deg_results)
#        Symbol P_counts P_study_names AdjP_counts AdjP_study_names
# 1        NEK8        0             0           0                0
# 2         IPP        0             0           0                0
# 3  CDKN2AIPNL        0             0           0                0
# 4 RP11-346C20        0             0           0                0
# 5       MAPK7        0             0           0                0
# 6  CTC-523E23        0             0           0                0
for (i in 1:5) {
  deg_results[i,2] <- sum(unlist(lapply(pfilt_symbols, function(x) {deg_results$Symbol[i] %in% x})))
  deg_results[i,3] <- paste(names(unlist(lapply(pfilt_symbols, function(x) {all_symbols[i] %in% x})))[unlist(lapply(pfilt_symbols, function(x) {all_symbols[i] %in% x}))], collapse = " ")
  deg_results[i,4] <- sum(unlist(lapply(adjpfilt_symbols, function(x) {deg_results$Symbol[i] %in% x})))
  deg_results[i,5] <- paste(names(unlist(lapply(adjpfilt_symbols, function(x) {all_symbols[i] %in% x})))[unlist(lapply(adjpfilt_symbols, function(x) {all_symbols[i] %in% x}))], collapse = " ")
}
head(deg_results)
# Looks right.
for (i in 1:nrow(deg_results)) {
  deg_results[i,2] <- sum(unlist(lapply(pfilt_symbols, function(x) {deg_results$Symbol[i] %in% x})))
  deg_results[i,3] <- paste(names(unlist(lapply(pfilt_symbols, function(x) {all_symbols[i] %in% x})))[unlist(lapply(pfilt_symbols, function(x) {all_symbols[i] %in% x}))], collapse = " ")
  deg_results[i,4] <- sum(unlist(lapply(adjpfilt_symbols, function(x) {deg_results$Symbol[i] %in% x})))
  deg_results[i,5] <- paste(names(unlist(lapply(adjpfilt_symbols, function(x) {all_symbols[i] %in% x})))[unlist(lapply(adjpfilt_symbols, function(x) {all_symbols[i] %in% x}))], collapse = " ")
}
studs <- c("gomez", "meaburn", "gosch", "obermoser", "dusek", "rusch", "larocca", "karlovich")
deg_results[,(length(deg_results)+1):(length(deg_results) + 8)] <- 0
names(deg_results)[6:13] <- studs
deg_results[,14] <- 0
names(deg_results)[14] <- "P_value_study_count"
deg_results[,15] <- 0
names(deg_results)[15] <- "AdjP_value_study_count"
for (i in 1:nrow(deg_results)) {
  for (j in 6:13) {
    deg_results[i,j] <- sum(grepl(names(deg_results)[j], unlist(str_split(deg_results[i,3], pattern = " "))))
  }
}
for (i in 1:nrow(deg_results)) {
  deg_results[i,14] <- as.numeric(grepl("gomez", deg_results$P_study_names[i])) +
    as.numeric(grepl("meaburn", deg_results$P_study_names[i])) +
    as.numeric(grepl("gosch", deg_results$P_study_names[i])) +
    as.numeric(grepl("obermoser", deg_results$P_study_names[i])) +
    as.numeric(grepl("dusek", deg_results$P_study_names[i])) +
    as.numeric(grepl("rusch", deg_results$P_study_names[i])) +
    as.numeric(grepl("larocca", deg_results$P_study_names[i])) +
    as.numeric(grepl("karlovich", deg_results$P_study_names[i]))
  deg_results[i,15] <- as.numeric(grepl("gomez", deg_results$AdjP_study_names[i])) +
    as.numeric(grepl("meaburn", deg_results$AdjP_study_names[i])) +
    as.numeric(grepl("gosch", deg_results$AdjP_study_names[i])) +
    as.numeric(grepl("obermoser", deg_results$AdjP_study_names[i])) +
    as.numeric(grepl("dusek", deg_results$AdjP_study_names[i])) +
    as.numeric(grepl("rusch", deg_results$AdjP_study_names[i])) +
    as.numeric(grepl("larocca", deg_results$AdjP_study_names[i])) +
    as.numeric(grepl("karlovich", deg_results$P_study_names[i])) 
}
head(deg_results)
write.csv(deg_results, file = paste(output_dir, "flexible_gene_scores.csv", sep = ""))

## Supplementary table S2 ----
deg_results <- read.csv(paste(output_dir, "flexible_gene_scores.csv", sep = ""))
head(deg_results[order(deg_results$P_value_study_count, decreasing = T), c("Symbol", "P_value_study_count")],20)
#       Symbol P_value_study_count
# 613    VMA21                   7
# 9      LRIG1                   6
# 14   CCDC88A                   6
# 25      <NA>                   6
# 54       AGL                   6
# 62     USP49                   6
# 83    DIS3L2                   6
# 88    SMURF2                   6
# 106    FSD1L                   6
# 110  PLEKHA3                   6
# 165 SLC25A16                   6
# 179  C9orf85                   6
# 219    DCHS1                   6
# 267     RNF7                   6
# 286      TDG                   6
# 307      FN1                   6
# 351     HPS4                   6
# 367   ZNF273                   6
# 370     BIVM                   6
# 397     MZF1                   6
dim(deg_results[deg_results$P_value_study_count >= 6,])
# [1] 123  16
supp_table <- deg_results[deg_results$P_value_study_count >= 6,"Symbol"]
library(biomaRt)
mart <- useDataset(dataset = "hsapiens_gene_ensembl", mart = useMart("ENSEMBL_MART_ENSEMBL"))
supp_table <- getBM(attributes = c("hgnc_symbol","description"),
                    filters = "hgnc_symbol",
                    values = supp_table,
                    mart = mart)
write.csv(supp_table, file = paste(output_dir, "supplementary_table_2.csv", sep = ""))