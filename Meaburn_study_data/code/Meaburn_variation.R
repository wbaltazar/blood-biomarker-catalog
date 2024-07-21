## Date: July 12 2024

## Calculate statistics for stability in Meaburn study
## INPUT: Expr and pheno data
## OUTPUT: .csv file of heritability, standard deviation, and variancePartition statistics

### SUMMARY: How stable are gene expression measurements across 4 hours? Are these repeatable 10 months later?
# Data from sets of healthy, 12-year-old twins were collected 10 AM and 2 PM on the same day. This
# experiment was repeated 10 months later to compare with the results from the first day.

# While not explicitly stated by the experimenters, we infer that the day 1 and day 2 samples were
# processed in separate batches due to differences in the control probe expression levels b/w the two
# groups. So, for consistency's sake, we separate the statistics according to day of sampling.

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

## INPUT: .CEL files from GSE14844
input_dir <- "~/Desktop/work_repo/data/GSE14844_RAW/"
## OUTPUT: tables and graphics
output_dir <- "~/Desktop/work_repo/github/Meaburn_study_data/output/"

## Separate the samples by day ----
gset <- getGEO("GSE14844", GSEMatrix =TRUE, AnnotGPL=TRUE)
pheno_data_org <- pData(gset$GSE14844_series_matrix.txt.gz)
feature_data <- fData(gset$GSE14844_series_matrix.txt.gz)
names(pheno_data_org)

# filter out traits we want
pheno_data <- pheno_data_org %>% 
  dplyr::select(title = title,
                geo_accession = geo_accession,
                sex = `gender:ch1`)
pheno_data$time = str_extract(pheno_data_org$title, "Day.*")
pheno_data$id = str_extract(pheno_data$title, "TD.....")
pheno_data$time[12] <- "Daytwo_10am"
dim(pheno_data)
# [1] 36  5
table(pheno_data$sex, pheno_data$time)
#        Dayone_10am Dayone_2pm Daytwo_10am Daytwo_2pm
# Female           8          8           6          8
# Male             2          1           1          2
## Suffices for analysis?
p1 <- pheno_data[grep("Dayone", pheno_data$time),]
p2 <- pheno_data[grep("Daytwo", pheno_data$time),]

# Analysis of day 1 data ----
celfiles <- list.files(input_dir, full.names = T)
c1 <- celfiles[which(str_extract(celfiles, "GSM\\d{6}") %in% p1$geo_accession)]
length(c1)
data <- ReadAffy(filenames = c1)
data <- affy::rma(data)
norm_expr <- exprs(data)
dim(norm_expr)
table(p1$sex, p1$id)
#        TD19901 TD19902 TD23461 TD23462 TD24111 TD24112 TD32681 TD32682 TD51281 TD51282
# Female       2       2       2       2       2       2       0       0       2       2
# Male         0       0       0       0       0       0       2       1       0       0


## varianceParitition ----
names(p1)
# [1] "title"         "geo_accession" "sex"           "time"          "id"    
names(p1)[5] <- "subject"
formula <- ~ (1|subject) + (1|time) + (1|sex) + (1|time:sex)
param <- SnowParam(4, "SOCK", progressbar = TRUE)
varPart <- fitExtractVarPartModel(exprObj = norm_expr, formula = formula, data = p1, BPPARAM = param)
vp <- sortCols(varPart)
pdf(paste(output_dir, "day1_vp_violin_plot.pdf", sep = ""))
plotVarPart(vp)
dev.off()
# Canonical Correlation Analysis
form <- ~ subject + time + sex
C <- canCorPairs(form, p1)
pdf(paste(output_dir, "day1_cca.pdf", sep = ""))
plotCorrMatrix(C)
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
# 1007_s_at    1053_at     117_at     121_at  1255_g_at    1294_at 
# 0.11101551 0.07088556 0.12074206 0.13273315 0.07930451 0.16932864 
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
pdf(paste(output_dir, "day1_statistic_correlations.pdf", sep = ""))
pheatmap(cor(data_var[,-length(names(data_var))]), display_numbers = T, fontsize_number = 4)
dev.off()
## Save ----
write.csv(data_var, paste(output_dir, "Meaburn_day1_variation.csv", sep = ""), row.names = T)

# Analysis of day 2 data ----
c2 <- celfiles[which(str_extract(celfiles, "GSM\\d{6}") %in% p2$geo_accession)]
length(c2) # [1] 17
data <- ReadAffy(filenames = c2)
data <- affy::rma(data)
norm_expr <- exprs(data)
dim(norm_expr)
table(p2$sex, p2$id)
#        TD19901 TD19902 TD23461 TD23462 TD24111 TD24112 TD32681 TD32682 TD51281 TD51282
# Female       1       2       2       1       2       2       0       0       2       2
# Male         0       0       0       0       0       0       2       1       0       0

## varianceParitition ----
names(p2)
# [1] "title"         "geo_accession" "sex"           "time"          "id"    
names(p2)[5] <- "subject"
param <- SnowParam(4, "SOCK", progressbar = TRUE)
varPart <- fitExtractVarPartModel(exprObj = norm_expr, formula = formula, data = p2, BPPARAM = param)
vp <- sortCols(varPart)
pdf(paste(output_dir, "day2_vp_violin_plot.pdf", sep = ""))
plotVarPart(vp)
dev.off()
# Canonical Correlation Analysis
form <- ~ subject + time + sex
C <- canCorPairs(form, p2)
pdf(paste(output_dir, "day2_cca.pdf", sep = ""))
plotCorrMatrix(C)
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
pdf(paste(output_dir, "day2_statistic_correlations.pdf", sep = ""))
pheatmap(cor(data_var[,-length(names(data_var))]), display_numbers = T, fontsize_number = 4)
dev.off()
## Save ----
write.csv(data_var, paste(output_dir, "Meaburn_day2_variation.csv", sep = ""), row.names = T)