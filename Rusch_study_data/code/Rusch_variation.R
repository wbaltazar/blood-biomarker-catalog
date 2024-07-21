## Date: Jul 12 2024

## You should run "limma_QC_Rusch.R" before running this script to remove poor quality arrays from input directory.

## Calculate statistics for stability in Rusch study
## INPUT: Expr and pheno data
## OUTPUT: .csv file of heritability, standard deviation, and variancePartition statistics

### SUMMARY: Military-age males experiencing insomnia had blood drawn before and after 3 months.

## Load libraries ----

library(AnnotationDbi)
library(limma)
library(oligo)
library(tidyverse)
library(cowplot)
library(GEOquery)
library(variancePartition)
library(heritability)
library(pheatmap)
library(affy)
library(affycoretools)
library(hgu133plus2.db)

## INPUT: QC'd .CEL files from GSE81761
input_dir <- "~/Desktop/work_repo/data/GSE81761_RAW/"
## OUTPUT: tables and graphics
output_dir <- "~/Desktop/work_repo/github/Rusch_study_data/output/"

# Load in data ----
data <- ReadAffy(celfile.path = input_dir)
data <- affy::rma(data)
norm_expr <- exprs(data)
colnames(norm_expr) <- str_extract(colnames(norm_expr), pattern = "GSM\\d+")
dim(norm_expr)
# [1] 54675    38

# Phenotype data
gset <- getGEO("GSE81761", GSEMatrix =TRUE, AnnotGPL=TRUE)
pheno_data_org <- pData(gset$GSE81761_series_matrix.txt.gz)
pheno_data <- pheno_data_org %>% 
  dplyr::select(title = title, 
                geo_accession = geo_accession,
                source = source_name_ch1,
                age = `age:ch1`,
                ethnicity = `ethnicity:ch1`,
                race = `race:ch1`,
                sex = `Sex:ch1`,
                time = `timepoint:ch1`)
pheno_data <- pheno_data[pheno_data$geo_accession %in% colnames(norm_expr),]
dim(pheno_data)
# [1] 38  8
pheno_data$subject <- str_extract(string = pheno_data$source, pattern = "Subject \\d+")
pheno_data$time[pheno_data$time != "Baseline"] <- rep(x = "Followup", times = length(pheno_data$time[pheno_data$time != "Baseline"]))
table(pheno_data$subject)
# Subject 10 Subject 108 Subject 117  Subject 37  Subject 38  Subject 46  Subject 49  Subject 53  Subject 63  Subject 66  Subject 68 
# 2           1           2           2           2           2           2           2           2           2           2 
# Subject 69  Subject 72  Subject 73  Subject 76  Subject 80   Subject 9  Subject 94  Subject 95  Subject 98 
# 2           2           2           1           2           2           2           2           2 

# Variation analysis ----

## varianceParitition ----
names(pheno_data)
# [1] "title"         "geo_accession" "source"        "age"           "ethnicity"     "race"         
# [7] "sex"           "time"          "subject"  
formula <-  ~ as.numeric(age) + (1|subject) + (1|time) + (1|ethnicity) + (1|race) + (1|time:ethnicity) + (1|time:race)
param <- SnowParam(4, "SOCK", progressbar = TRUE)
varPart <- fitExtractVarPartModel(exprObj = norm_expr, formula = formula, data = pheno_data, BPPARAM = param)
vp <- sortCols(varPart)
pdf(paste(output_dir, "vp_violin_plot.pdf", sep = ""))
plotVarPart(vp)
dev.off()
# Canonical Correlation Analysis
form <- ~ subject + time + race + ethnicity + age
C <- canCorPairs(form, pheno_data)
pdf(paste(output_dir, "cca.pdf", sep = ""))
plotCorrMatrix(C)
dev.off()

## standard deviation ----
avgexpr <- rowMeans(norm_expr)
sd_total <- apply(norm_expr, 1, sd)

colnames(norm_expr) <- pheno_data$subject
person <- as.factor(colnames(norm_expr))
within_person_sd <- list()
for (i in 1:length(levels(person))){
  if (length(grep(paste0("\\b",levels(person)[[i]],"\\b",collapse = ""), colnames(norm_expr))) != 1) {
    col_index <- grep(paste0("\\b",levels(person)[[i]],"\\b",collapse = ""), colnames(norm_expr))
    within_person_sd[[i]] <- apply(norm_expr[, col_index], 1, sd)
  } else {
    col_index <- grep(paste0("\\b",levels(person)[[i]],"\\b",collapse = ""), colnames(norm_expr))
    within_person_sd[[i]] <- rep(NA, times = length(norm_expr[, col_index]))
  }
}
within_person_sd2 <- as.data.frame(matrix(unlist(within_person_sd), nrow = length(within_person_sd[[i]])))
row.names(within_person_sd2) <- names(within_person_sd[[i]])
names(within_person_sd2) <- levels(person)
within_person_sd2
sd_within <- rowMeans(within_person_sd2, na.rm = TRUE)
head(sd_within)
# 1007_s_at   1053_at    117_at    121_at 1255_g_at   1294_at 
# 0.09887649 0.12058775 0.14011709 0.11913603 0.10232581 0.12182798 
data_var <- as.data.frame(matrix(cbind(sd_within, sd_total), nrow = length(sd_within)))
rownames(data_var) <- rownames(norm_expr)
names(data_var) <- c("Within Variation (SD)", "Total Variation (SD)")
data_var$Rs <- data_var$`Within Variation (SD)` / data_var$`Total Variation (SD)`
data_var$`Average Expression` <- avgexpr
data_var$temp <- rownames(data_var)
vp$temp <- rownames(vp)
data_var <- merge(data_var, vp, by = "temp", sort = F)
head(data_var)

## heritability ----
prepared_covariates <- pheno_data[,c(5,6,8)] # Does not include age, continuous variable
results <- vector("list", nrow(norm_expr))
pheno_data$subject <- factor(pheno_data$subject)

for (i in 1:nrow(norm_expr)) {
  results[[i]] <- repeatability(data.vector = norm_expr[i,], 
                                geno.vector = pheno_data$subject,
                                covariates.frame = prepared_covariates)
} 
broad_sense_heritability <- do.call(rbind, lapply(results, function(x) {
  as.data.frame(t(unlist(x)))
}))
head(broad_sense_heritability)
broad_sense_heritability <- broad_sense_heritability[,1:3]
broad_sense_heritability$temp <- rownames(norm_expr)
data_var <- merge(data_var, broad_sense_heritability, by = "temp")
rownames(data_var) <- data_var$temp
data_var <- data_var[,-1]
genenames <- rownames(data_var)
genenames <- mapIds(hgu133plus2.db, genenames, "SYMBOL", "PROBEID")
genenames <- make.names(ifelse(is.na(genenames), names(genenames), unname(genenames)), unique = T)
data_var$Symbol <- genenames

## Heatmap ----
pdf(paste(output_dir, "statistic_correlations.pdf", sep = ""))
pheatmap(cor(data_var[,-length(names(data_var))]), display_numbers = T, fontsize_number = 4)
dev.off()
## Save ----
write.csv(data_var, paste(output_dir, "Rusch_variation.csv", sep = ""), row.names = T)