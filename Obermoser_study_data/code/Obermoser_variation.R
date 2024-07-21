## Date: Jul 19 2024

## Calculate statistics for stability in Obermoser studies
## INPUT: Expr and pheno data
## OUTPUT: .csv file of heritability, standard deviation, and variancePartition statistics

### SUMMARY: Do vaccine injections induce changes in gene expression?
# Data from multiple healthy patients sampled before, baseline, and after vaccination (w/ saline controls)
# We include the type of vaccination in this analysis to potentially explain the variability

# The data was collected in three cohorts from either the vein or the finger. Additionally, the data collected
# across cohorts measures different time intervals, including an hour-by-hour measurement.

setwd("~/Desktop/work_repo")

library(AnnotationDbi)
library(beadarray)
library(illuminaio)
library(illuminaHumanv3.db)
library(lumi)
library(limma)
library(oligo)
library(tidyverse)
library(cowplot)
library(GEOquery)
library(variancePartition)
library(heritability)
library(pheatmap)

## INPUT: Raw file folder from GSE30101
input_dir <- "~/Desktop/work_repo/data/GSE30101_RAW/"
## OUTPUT: variation tables and graphics
output_dir <- "~/Desktop/work_repo/github/Obermoser_study_data/output/"

# Load in the phenotype data ----

gset <- getGEO("GSE30101", AnnotGPL = TRUE, GSEMatrix = TRUE)
pheno_data_org <- pData(gset$GSE30101_series_matrix.txt.gz)

pheno_data <- pheno_data_org %>% 
  dplyr::select(title, 
                geo_accession, 
                subject = `subject id:ch1`,
                age = `age:ch1`, 
                time = `day:ch1`,
                treat = `vaccine:ch1`,
                sex = `gender:ch1`,
                race = `race:ch1`, 
                ethnicity = `ethnicity:ch1`,
                array = description.1,
                sample = source_name_ch1,
                group = `sample set:ch1`)
pheno_data <- pheno_data[pheno_data$sample == "Whole Blood",]
p1v <- pheno_data[grep("Training_Set_Vein", pheno_data$group),]
p1f <- pheno_data[grep("Test_Set_Finger", pheno_data$group),]
p2v <- pheno_data[grep("Test_Set_Vein", pheno_data$group),]
p2f <- pheno_data[grep("Training_Set_Finger", pheno_data$group),]
# p3v <- pheno_data[grep("Validation_Vein", pheno_data$group),]

## Load in probe annotation file ---
BGX <- readBGX(file = "./data/GSE30101_RAW/GPL6947_HumanHT-12_V3_0_R1_11283641_A.bgx.gz")
## Raw data file names ----
files <- list.files("./data/GSE30101_RAW/" , pattern = "GSM.*\\.txt.gz$", full.names = TRUE)

# Analysis 1: Training Set Vein (p1v) ----
## Load in files and organize ----
p1v <- p1v[-grep("Saline", p1v$treat),]
dim(p1v)
# 161 12
temp <- files[str_extract(files, "GSM......") %in% p1v$geo_accession]
length(temp) # [1] 160
x <- lapply(
  temp,
  function(x) read.delim(file = x))
x <- do.call(cbind, x)
dim(x)
# [1] 48803   483
rownames(x) <- x[,1]
x <- as.data.frame(x)
detectionpvalues <- x[,grep("Detection",colnames(x))]
x <- x[,grep("AVG",colnames(x))]
data <- new('EListRaw')
data@.Data[[1]] <- 'illumina'
data@.Data[[3]] <- BGX
data@.Data[[4]] <- x
data@.Data[[5]] <- NULL
data$E <- x
data$genes <- BGX
data$other$Detection <- detectionpvalues
data_norm <- neqc(data)
controls <- data_norm$genes$controls$Probe_Id
expressed <- rowSums(data_norm$other$Detection < 0.05) >= 3
data_norm$E <- data_norm$E[expressed,]
data_norm$E <- data_norm$E[-which(rownames(data_norm$E) %in% controls),]
dim(data_norm$E)
# [1] 31003   161
norm_expr <- data_norm$E

## varianceParitition ----
formula <- ~ as.numeric(age) + (1|subject) + (1|time) + (1|sex) + (1|race) + (1|treat) + (1|time:treat) + (1|time:sex) + (1|time:race) + time:as.numeric(age)
param <- SnowParam(4, "SOCK", progressbar = TRUE)
varPart <- fitExtractVarPartModel(exprObj = norm_expr, formula = formula, data = p1v, BPPARAM = param)
vp <- sortCols(varPart)
pdf(paste(output_dir, "p1v_vp_violin_plot.pdf", sep = ""))
plotVarPart(vp)
dev.off()
# Canonical Correlation Analysis
form <- ~ subject + time + sex + treat + race
C <- canCorPairs(form, p1v)
pdf(paste(output_dir, "p1v_cca.pdf", sep = ""))
plotCorrMatrix(C)
dev.off()

## standard deviation ----
avgexpr <- rowMeans(norm_expr)
sd_total <- apply(norm_expr, 1, sd)

colnames(norm_expr) <- p1v$subject
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
# ILMN_2209417 ILMN_1765401 ILMN_1652171 ILMN_1652170 ILMN_1792710 ILMN_1703558 
# 0.1275598    0.3500975    0.1903379    0.3849829    0.3832176    0.1726510 
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
prepared_covariates <- p1v[,c("time","treat","sex","race","ethnicity")] # Continuous variable (age) not included.
prepared_covariates$ethnicity[is.na(prepared_covariates$ethnicity)] <- "Not reported"
results <- vector("list", nrow(norm_expr))
p1v$id <- factor(p1v$subject)

for (i in 1:nrow(norm_expr)) {
  results[[i]] <- repeatability(data.vector = norm_expr[i,], 
                                geno.vector = p1v$id,
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
genenames <- mapIds(illuminaHumanv3.db, genenames, "SYMBOL", "PROBEID")
genenames <- make.names(ifelse(is.na(genenames), names(genenames), unname(genenames)), unique = T)
data_var$Symbol <- genenames

## Heatmap ----
pdf(paste(output_dir, "p1v_statistic_correlations.pdf", sep = ""))
pheatmap(cor(data_var[,-length(names(data_var))]), display_numbers = T, fontsize_number = 4)
dev.off()
## Save ----
write.csv(data_var, paste(output_dir, "Obermoser_p1v_variation.csv", sep = ""), row.names = T)

# Analysis 2: Test Set Finger (p1f) ----
## Load in files and organize ----
temp <- files[str_extract(files, "GSM......") %in% p1f$geo_accession]
length(temp) # [1] 49
x <- lapply(
  temp,
  function(x) read.delim(file = x))
x <- do.call(cbind, x)
dim(x)
# [1] 48803   147
rownames(x) <- x[,1]
x <- as.data.frame(x)
detectionpvalues <- x[,grep("Detection",colnames(x))]
x <- x[,grep("AVG",colnames(x))]
data <- new('EListRaw')
data@.Data[[1]] <- 'illumina'
data@.Data[[3]] <- BGX
data@.Data[[4]] <- x
data@.Data[[5]] <- NULL
data$E <- x
data$genes <- BGX
data$other$Detection <- detectionpvalues
data_norm <- neqc(data)
controls <- data_norm$genes$controls$Probe_Id
expressed <- rowSums(data_norm$other$Detection < 0.05) >= 3
data_norm$E <- data_norm$E[expressed,]
data_norm$E <- data_norm$E[-which(rownames(data_norm$E) %in% controls),]
dim(data_norm$E)
# [1] 22413    49
norm_expr <- data_norm$E

## varianceParitition ----
formula <- ~ as.numeric(age) + (1|subject) + (1|time) + (1|sex) + (1|race) + (1|treat) + (1|time:treat) + (1|time:sex) + (1|time:race) + time:as.numeric(age)
param <- SnowParam(4, "SOCK", progressbar = TRUE)
varPart <- fitExtractVarPartModel(exprObj = norm_expr, formula = formula, data = p1f, BPPARAM = param)
vp <- sortCols(varPart)
pdf(paste(output_dir, "p1f_vp_violin_plot.pdf", sep = ""))
plotVarPart(vp)
dev.off()
# Canonical Correlation Analysis
form <- ~ subject + time + sex + treat + race
C <- canCorPairs(form, p1f)
pdf(paste(output_dir, "p1f_cca.pdf", sep = ""))
plotCorrMatrix(C)
dev.off()

## standard deviation ----
avgexpr <- rowMeans(norm_expr)
sd_total <- apply(norm_expr, 1, sd)

colnames(norm_expr) <- p1f$subject
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
# ILMN_2209417 ILMN_1765401 ILMN_1652171 ILMN_1652170 ILMN_1792710 ILMN_1703558 
# 0.1275598    0.3500975    0.1903379    0.3849829    0.3832176    0.1726510 
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
prepared_covariates <- p1f[,c("time","treat","sex","race","ethnicity")] # Continuous variable (age) not included.
prepared_covariates$ethnicity[is.na(prepared_covariates$ethnicity)] <- "Not reported"
results <- vector("list", nrow(norm_expr))
p1f$id <- factor(p1f$subject)

for (i in 1:nrow(norm_expr)) {
  results[[i]] <- repeatability(data.vector = norm_expr[i,], 
                                geno.vector = p1f$id,
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
genenames <- mapIds(illuminaHumanv3.db, genenames, "SYMBOL", "PROBEID")
genenames <- make.names(ifelse(is.na(genenames), names(genenames), unname(genenames)), unique = T)
data_var$Symbol <- genenames

## Heatmap ----
pdf(paste(output_dir, "p1f_statistic_correlations.pdf", sep = ""))
pheatmap(cor(data_var[,-length(names(data_var))]), display_numbers = T, fontsize_number = 4)
dev.off()
## Save ----
write.csv(data_var, paste(output_dir, "Obermoser_p1f_variation.csv", sep = ""), row.names = T)

# Analysis 3: Test Set Vein (p2v) ----
## Load in files and organize ----
temp <- files[str_extract(files, "GSM......") %in% p2v$geo_accession]
length(temp) # [1] 143
x <- lapply(
  temp,
  function(x) read.delim(file = x))
x <- do.call(cbind, x)
dim(x)
# [1] 48803   429
rownames(x) <- x[,1]
x <- as.data.frame(x)
detectionpvalues <- x[,grep("Detection",colnames(x))]
x <- x[,grep("AVG",colnames(x))]
data <- new('EListRaw')
data@.Data[[1]] <- 'illumina'
data@.Data[[3]] <- BGX
data@.Data[[4]] <- x
data@.Data[[5]] <- NULL
data$E <- x
data$genes <- BGX
data$other$Detection <- detectionpvalues
data_norm <- neqc(data)
controls <- data_norm$genes$controls$Probe_Id
expressed <- rowSums(data_norm$other$Detection < 0.05) >= 3
data_norm$E <- data_norm$E[expressed,]
data_norm$E <- data_norm$E[-which(rownames(data_norm$E) %in% controls),]
dim(data_norm$E)
# [1] 27431   143
norm_expr <- data_norm$E

## varianceParitition ----
formula <- ~ as.numeric(age) + (1|subject) + (1|time) + (1|sex) + (1|race) + (1|treat) + (1|time:treat) + (1|time:sex) + (1|time:race) + time:as.numeric(age)
param <- SnowParam(4, "SOCK", progressbar = TRUE)
varPart <- fitExtractVarPartModel(exprObj = norm_expr, formula = formula, data = p2v, BPPARAM = param)
vp <- sortCols(varPart)
pdf(paste(output_dir, "p2v_vp_violin_plot.pdf", sep = ""))
plotVarPart(vp)
dev.off()
# Canonical Correlation Analysis
form <- ~ subject + time + sex + treat + race
C <- canCorPairs(form, p2v)
pdf(paste(output_dir, "p2v_cca.pdf", sep = ""))
plotCorrMatrix(C)
dev.off()

## standard deviation ----
avgexpr <- rowMeans(norm_expr)
sd_total <- apply(norm_expr, 1, sd)

colnames(norm_expr) <- p2v$subject
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
# ILMN_2209417 ILMN_1765401 ILMN_1652171 ILMN_1652170 ILMN_1792710 ILMN_1703558 
# 0.1275598    0.3500975    0.1903379    0.3849829    0.3832176    0.1726510 
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
prepared_covariates <- p2v[,c("time","treat","sex","race","ethnicity")] # Continuous variable (age) not included.
prepared_covariates$ethnicity[is.na(prepared_covariates$ethnicity)] <- "Not reported"
results <- vector("list", nrow(norm_expr))
p2v$id <- factor(p2v$subject)

for (i in 1:nrow(norm_expr)) {
  results[[i]] <- repeatability(data.vector = norm_expr[i,], 
                                geno.vector = p2v$id,
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
genenames <- mapIds(illuminaHumanv3.db, genenames, "SYMBOL", "PROBEID")
genenames <- make.names(ifelse(is.na(genenames), names(genenames), unname(genenames)), unique = T)
data_var$Symbol <- genenames

## Heatmap ----
pdf(paste(output_dir, "p2v_statistic_correlations.pdf", sep = ""))
pheatmap(cor(data_var[,-length(names(data_var))]), display_numbers = T, fontsize_number = 4)
dev.off()
## Save ----
write.csv(data_var, paste(output_dir, "Obermoser_p2v_variation.csv", sep = ""), row.names = T)

# Analysis 4: Training Set Finger (p2f) ----
## Load in files and organize ----
temp <- files[str_extract(files, "GSM......") %in% p2f$geo_accession]
length(temp) # [1] 185
x <- lapply(
  temp,
  function(x) read.delim(file = x))
x <- do.call(cbind, x)
dim(x)
# [1] 48803   417
rownames(x) <- x[,1]
x <- as.data.frame(x)
detectionpvalues <- x[,grep("Detection",colnames(x))]
x <- x[,grep("AVG",colnames(x))]
data <- new('EListRaw')
data@.Data[[1]] <- 'illumina'
data@.Data[[3]] <- BGX
data@.Data[[4]] <- x
data@.Data[[5]] <- NULL
data$E <- x
data$genes <- BGX
data$other$Detection <- detectionpvalues
data_norm <- neqc(data)
controls <- data_norm$genes$controls$Probe_Id
expressed <- rowSums(data_norm$other$Detection < 0.05) >= 3
data_norm$E <- data_norm$E[expressed,]
data_norm$E <- data_norm$E[-which(rownames(data_norm$E) %in% controls),]
dim(data_norm$E)
# [1] 23484   185
norm_expr <- data_norm$E

## varianceParitition ----
formula <- ~ as.numeric(age) + (1|subject) + (1|time) + (1|sex) + (1|race) + (1|treat) + (1|time:treat) + (1|time:sex) + (1|time:race) + time:as.numeric(age)
param <- SnowParam(4, "SOCK", progressbar = TRUE)
varPart <- fitExtractVarPartModel(exprObj = norm_expr, formula = formula, data = p2f, BPPARAM = param)
vp <- sortCols(varPart)
pdf(paste(output_dir, "p2f_vp_violin_plot.pdf", sep = ""))
plotVarPart(vp)
dev.off()
# Canonical Correlation Analysis
form <- ~ subject + time + sex + treat + race
C <- canCorPairs(form, p2f)
pdf(paste(output_dir, "p2f_cca.pdf", sep = ""))
plotCorrMatrix(C)
dev.off()

## standard deviation ----
avgexpr <- rowMeans(norm_expr)
sd_total <- apply(norm_expr, 1, sd)

colnames(norm_expr) <- p2f$subject
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
# ILMN_2209417 ILMN_1765401 ILMN_1652171 ILMN_1652170 ILMN_1792710 ILMN_1703558 
# 0.1275598    0.3500975    0.1903379    0.3849829    0.3832176    0.1726510 
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
prepared_covariates <- p2f[,c("time","treat","sex","race","ethnicity")] # Continuous variable (age) not included.
prepared_covariates$ethnicity[is.na(prepared_covariates$ethnicity)] <- "Not reported"
results <- vector("list", nrow(norm_expr))
p2f$id <- factor(p2f$subject)

for (i in 1:nrow(norm_expr)) {
  results[[i]] <- repeatability(data.vector = norm_expr[i,], 
                                geno.vector = p2f$id,
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
genenames <- mapIds(illuminaHumanv3.db, genenames, "SYMBOL", "PROBEID")
genenames <- make.names(ifelse(is.na(genenames), names(genenames), unname(genenames)), unique = T)
data_var$Symbol <- genenames

## Heatmap ----
pdf(paste(output_dir, "p2f_statistic_correlations.pdf", sep = ""))
pheatmap(cor(data_var[,-length(names(data_var))]), display_numbers = T, fontsize_number = 4)
dev.off()
## Save ----
write.csv(data_var, paste(output_dir, "Obermoser_p2f_variation.csv", sep = ""), row.names = T)