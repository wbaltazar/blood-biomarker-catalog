## Date: Jul 11 2024

## You should run "pre_QC_Dusek.R" before running this script to remove poor quality arrays from input directory.

## Calculate statistics for stability in Dusek study
## INPUT: Expr and pheno data
## OUTPUT: .csv file of heritability, standard deviation, and variancePartition statistics

### SUMMARY: Healthy adults took 8 weeks eliciting the relaxation response. They were inexperienced
# and meditated with the help of a 20 minute guided exercise. Blood was drawn before and after the
# intervention period and hybridized to Affymetrix arrays for comparison.

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

## INPUT: pre-QC'd .CEL files
input_dir <- "~/Desktop/work_repo/data/GSE10041_RAW/"
## OUTPUT: tables and graphics
output_dir <- "~/Desktop/work_repo/github/Dusek_study_data/output/"


# Load in data ----
data <- ReadAffy(celfile.path = input_dir)
data <- affy::rma(data)
norm_expr <- exprs(data)
colnames(norm_expr) <- gsub(pattern = "(.CEL.gz)", replacement = "", x = colnames(norm_expr))
dim(norm_expr)

# Phenotype data
gset <- getGEO("GSE10041", GSEMatrix =TRUE, AnnotGPL=TRUE)
pheno_data_org <- pData(gset$GSE10041_series_matrix.txt.gz)
pheno_data <- pheno_data_org %>% 
  dplyr::select(title = title,
                geo_accession = geo_accession)
pheno_data$time <- str_extract(pheno_data_org$title, "N.")
pheno_data$id <- str_extract(pheno_data_org$title, "-.*") %>% 
  str_replace("-","")
dim(pheno_data)
# [1] 72 4
length(grep("M", pheno_data$title)) # [1] 25 This is the correct number in the original and validation groups.
pheno_data <- pheno_data[-grep("M", pheno_data$title),]
dim(pheno_data)
# [1] 47  4
pheno_data <- pheno_data[pheno_data$geo_accession %in% colnames(norm_expr),]
dim(pheno_data)
# [1] 44  4
identical(pheno_data$geo_accession, colnames(norm_expr))
# Annotate sex
rps4y_vals <- norm_expr["201909_at",]
pheno_data$sex <- NA
for (i in 1:length(rps4y_vals)) {
  if (rps4y_vals[i] > 8) {
    pheno_data$sex[i] <- "Male"
  } else {
    pheno_data$sex[i] <- "Female"
  }
}
table(pheno_data$sex, pheno_data$id)
#        01 02 03 06 07 11 12 14 23 32 33 34 35 36 39 40 41 46 47 49 50 51 59
# Female  0  0  0  0  2  2  2  2  0  2  2  2  2  0  0  0  0  2  0  0  0  2  2
# Male    2  2  2  2  0  0  0  0  2  0  0  0  0  1  2  2  2  0  2  2  1  0  0

# Variation analysis ----

## varianceParitition ----
names(pheno_data)
# [1] "title"         "geo_accession" "time"          "id"            "sex"   
names(pheno_data)[4] <- "subject"
formula <- ~ (1|subject) + (1|time) + (1|sex)  + (1|time:sex)
param <- SnowParam(4, "SOCK", progressbar = TRUE)
varPart <- fitExtractVarPartModel(exprObj = norm_expr, formula = formula, data = pheno_data, BPPARAM = param)
vp <- sortCols(varPart)
pdf(paste(output_dir, "vp_violin_plot.pdf", sep = ""))
plotVarPart(vp)
dev.off()
# Canonical Correlation Analysis
form <- ~ subject + time + sex
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
# 0.1307575 0.1623192 0.2389176 0.1488751 0.1056966 0.2072298
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
prepared_covariates <- pheno_data[,c("time","sex")] 
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
pheatmap(cor(data_var[,-length(names(data_var))]), display_numbers = T)
dev.off()
## Save ----
write.csv(data_var, paste(output_dir, "Dusek_variation.csv", sep = ""), row.names = T)