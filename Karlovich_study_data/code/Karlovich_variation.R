# Date: January 17 2025

## Calculate statistics for stability in Karlovich study
## INPUT: Expr and pheno data
## OUTPUT: .csv file of heritability, standard deviation, and variancePartition statistics

### SUMMARY: Blood was drawn from healthy volunteers over 2 weeks and 4 weeks, then at 3 months and 6 months.

# The dataset authors state that the first month of samples and the 3/6 month samples were
# processed in separate batches. We separate the statistics according to day of sampling.

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

## INPUT: .CEL files from GSE16028
input_dir <- "~/Desktop/work_repo/data/"
## OUTPUT: tables and graphics
output_dir <- "~/Desktop/work_repo/github/Karlovich_study_data/output/"

# Load in metadata and create batches ----
gset <- getGEO("GSE16028", GSEMatrix =TRUE, AnnotGPL=TRUE)
pheno_data_org <- pData(gset$GSE16028_series_matrix.txt.gz)
# select necessary columns
pheno_data <- pheno_data_org %>% 
  dplyr::select(title, 
                geo_accession,
                subject = `patient identifier:ch1`,
                age = `age (years):ch1`,
                sex = `sex:ch1`,
                time = `day:ch1`)
head(pheno_data)
#                                              title geo_accession subject age  sex    time
# GSM401059   Blood sample from patient 106 at Day 1     GSM401059     106  40 Male   Day 1
# GSM401060  Blood sample from patient 106 at Day 14     GSM401060     106  40 Male  Day 14
# GSM401061  Blood sample from patient 106 at Day 28     GSM401061     106  40 Male  Day 28
# GSM401062  Blood sample from patient 106 at Day 90     GSM401062     106  40 Male  Day 90
# GSM401063 Blood sample from patient 106 at Day 180     GSM401063     106  40 Male Day 180
# GSM401064   Blood sample from patient 108 at Day 1     GSM401064     108  24 Male   Day 1

p1 <- pheno_data[-grep("Day 90|Day 180", pheno_data$time), ]
dim(p1)
# [1] 66  6
p2 <- pheno_data[grep("Day 90|Day 180", pheno_data$time), ]
dim(p2)
# [1] 43  6

# Analysis of day 1 data ----
files <- list.files(paste(input_dir, "GSE16028_RAW", sep = ""), full.names = T)
files <- files[str_extract(files, "GSM\\d+") %in% p1$geo_accession]
length(files)
# [1] 66
data <- ReadAffy(filenames = files)
calls <- mas5calls.AffyBatch(data)
probe_pval <- assayData(calls)[["se.exprs"]] ## Returns p-values
minSamples <- min(colSums(table(p1$subject, p1$time)))
expressed <- rowSums(probe_pval < 0.05) >= minSamples

data <- affy::rma(data)
norm_expr <- exprs(data)
norm_expr <- norm_expr[expressed,]
dim(norm_expr)
# [1] 17635    66
table(p1$sex, p1$subject)
#        106 108 113 114 117 121 123 131 135 138 147 148 154 159 160 163 165 172 173 174 179 180
# Female   0   0   0   0   0   0   0   0   0   0   3   3   3   3   3   3   3   3   3   3   3   3
# Male     3   3   3   3   3   3   3   3   3   3   0   0   0   0   0   0   0   0   0   0   0   0

## variancePartition 1 ----
names(p1)
# "title"         "geo_accession" "subject"       "age"           "sex"           "time"     
formula <- ~ (1|subject) + as.numeric(age) + (1|time) + (1|sex) + (1|time:sex) + as.numeric(age):time
param <- SnowParam(4, "SOCK", progressbar = TRUE)
varPart <- fitExtractVarPartModel(exprObj = norm_expr, formula = formula, data = p1, BPPARAM = param)
vp <- sortCols(varPart)
pdf(paste(output_dir, "batch1_vp_violin_plot.pdf", sep = ""))
plotVarPart(vp)
dev.off()

## standard deviation ----
avgexpr <- rowMeans(norm_expr)
sd_total <- apply(norm_expr, 1, sd)

colnames(norm_expr) <- p1$subject
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
# 1053_at     117_at     121_at    1294_at    1316_at  1405_i_at 
# 0.17039357 0.22884449 0.09070506 0.17053974 0.13979569 0.23540320 
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
prepared_covariates <- p1[,c("time","sex")] 
results <- vector("list", nrow(norm_expr))
p1$subject <- factor(p1$subject)

for (i in 1:nrow(norm_expr)) {
  results[[i]] <- repeatability(data.vector = norm_expr[i,], 
                                geno.vector = p1$subject,
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
pdf(paste(output_dir, "batch1_statistic_heatmap.pdf", sep = ""))
pheatmap(cor(data_var[,-length(names(data_var))]), display_numbers = T, fontsize_number = 6)
dev.off()
## Save ----
write.csv(data_var, paste(output_dir, "Karlovich_batch1_variation.csv", sep = ""), row.names = T)

# Analysis of batch 2 data ----
files <- list.files(paste(input_dir, "GSE16028_RAW", sep = ""), full.names = T)
files <- files[str_extract(files, "GSM\\d+") %in% p2$geo_accession]
length(files)
# [1] 43
files <- files[-grep("GSM401157|GSM401087|GSM401148", files)]
length(files)
# [1] 40
p2 <- p2[-grep("GSM401087|GSM401148|GSM401157", p2$geo_accession),]
dim(p2)
# [1] 40  6
data <- ReadAffy(filenames = files)
calls <- mas5calls.AffyBatch(data)
probe_pval <- assayData(calls)[["se.exprs"]] ## Returns p-values
minSamples <- min(colSums(table(p2$subject, p2$time)))
print(minSamples)
# [1] 20
expressed <- rowSums(probe_pval < 0.05) >= minSamples

data <- affy::rma(data)
norm_expr <- exprs(data)
norm_expr <- norm_expr[expressed,]
dim(norm_expr)
# [1] 16352    40
colnames(norm_expr) <- gsub(pattern = "(.CEL.gz)", replacement = "", x = colnames(norm_expr)) # greps file extension and replaces it
dim(norm_expr)
df <- data.frame(norm_expr)
p2 <- p2[p2$geo_accession %in% colnames(df), ]
identical(colnames(df), p2$geo_accession)
# [1] TRUE
table(p2$sex, p2$subject)
#        106 108 113 114 117 121 123 131 135 138 147 148 154 159 160 163 165 172 173 179 180
# Female   0   0   0   0   0   0   0   0   0   0   2   2   2   2   2   2   2   1   2   2   2
# Male     2   2   2   2   2   1   2   2   2   2   0   0   0   0   0   0   0   0   0   0   0

## varianceParitition 2 ----
names(p2)
# [1] "title"         "geo_accession" "subject"       "age"           "sex"           "time"    
param <- SnowParam(4, "SOCK", progressbar = TRUE)
varPart <- fitExtractVarPartModel(exprObj = norm_expr, formula = formula, data = p2, BPPARAM = param)
vp <- sortCols(varPart)
pdf(paste(output_dir, "batch2_vp_violin_plot.pdf", sep = ""))
plotVarPart(vp) + theme(axis.text.x = element_text(size = 11))
dev.off()

## standard deviation ----
avgexpr <- rowMeans(norm_expr)
sd_total <- apply(norm_expr, 1, sd)

colnames(norm_expr) <- p2$subject
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
# 1007_s_at    1053_at     117_at     121_at  1255_g_at    1294_at 
# 0.10340392 0.06814007 0.08893897 0.14513779 0.06512726 0.04513723 
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
prepared_covariates <- p2[,c("time","sex")] 
results <- vector("list", nrow(norm_expr))
p2$subject <- factor(p2$subject)

for (i in 1:nrow(norm_expr)) {
  results[[i]] <- repeatability(data.vector = norm_expr[i,], 
                                geno.vector = p2$subject,
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
pdf(paste(output_dir, "batch2_statistic_heatmap.pdf", sep = ""))
pheatmap(cor(data_var[,-length(names(data_var))]), display_numbers = T, fontsize_number = 7)
dev.off()
## Save ----
write.csv(data_var, paste(output_dir, "Karlovich_batch2_variation.csv", sep = ""), row.names = T)