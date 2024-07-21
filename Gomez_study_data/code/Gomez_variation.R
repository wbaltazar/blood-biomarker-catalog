## Date: Jul 11 2024

## Calculate statistics for stability in Gomez study
## INPUT: Raw count file from GSE239282
## OUTPUT: .csv file of heritability, standard deviation, average expression, and variancePartition statistics

### SUMMARY: How stable are gene expression measurements after a 50 minute classical concert?

## Load libraries ----

library(edgeR)
library(limma)
library(AnnotationDbi)
library(Homo.sapiens)
library(cowplot)
library(GEOquery)
library(tidyverse)
library(variancePartition)
library(heritability)
library(pheatmap)

## INPUT: Raw count file from GSE239282 on Gene Expression Omnibus ----
input_dir <- "~/Desktop/work_repo/data/GSE239282_raw_counts.txt.gz"
## OUTPUT: tables and graphics ----
output_dir <- "~/Desktop/work_repo/github/Gomez_study_data/output/"

# Load in the data ----
gset <- getGEO(GEO = "GSE239282")
p <- phenoData(gset$GSE239282_series_matrix.txt.gz)
pheno_data_raw <- p@data
pheno_data <- pheno_data_raw %>% 
  dplyr::select(title = title,
                geo_accession = geo_accession,
                age = `age:ch1`,
                condition = `group:ch1`,
                participant = `subject id:ch1`,
                sex = `gender:ch1`)
pheno_data$time <- as.numeric(grepl("Timepoint-2", pheno_data$title))
table(pheno_data$time)
#  0  1 
# 30 30 
library(edgeR)
library(limma)
setwd("~/Desktop/work_repo")
x1 <- read.delim("./data/GSE239282_raw_counts.txt.gz")
#### x1 saved for remembering what the data looks like LOL
x <- x1
dim(x)
# [1] 60664    60
x <- DGEList(x)
geneid <- rownames(x) 
library(EnsDb.Hsapiens.v86)
genes <- mapIds(EnsDb.Hsapiens.v86, keys=geneid, column=c("SYMBOL"), 
                keytype="GENEID")
x$genes <- genes
detach("package:EnsDb.Hsapiens.v86", unload = TRUE)
pheno_data$match <- gsub("\\[|\\]", "",str_extract(pheno_data$title, "\\[(.*?)\\]"))
index <- match(substring(text = colnames(x$counts),first = 2,last = nchar(x = colnames(x$counts))), pheno_data$match)
pheno_data <- pheno_data[index,]
pheno_data <- remove_missing(pheno_data)
index <- match(pheno_data$match, substring(text = colnames(x$counts),first = 2,last = nchar(x = colnames(x$counts))))
x$counts <- x$counts[,index]
x$samples <- x$samples[index,]
identical(colnames(x$counts), rownames(x$samples))
# [1] TRUE
x$samples$group <- pheno_data$condition
x$samples$time <- pheno_data$time
x$samples$sex <- pheno_data$sex
x$samples$age <- pheno_data$age
x$samples$participant <- pheno_data$participant
x$samples$geo <- pheno_data$geo_accession
index <- which(x$samples$group == "ACD")
x$samples <- x$samples[-index,]
x$counts <- x$counts[,-index]
dim(x$counts)
# [1] 60664    28
keep.exprs <- filterByExpr(x, group = x$samples$time) 
x$counts <- x$counts[keep.exprs,] 
dim(x)
# [1] 15416    28
x <- calcNormFactors(x, method = "TMM")
id <- factor(x$samples$participant)
time <- factor(x$samples$time)
sex <- factor(x$samples$sex)
TS <- paste(sex,time, sep = ".")
TS <- factor(TS)
design <- model.matrix(~0+TS)
colnames(design) <- make.names(levels(TS))
v <- voom(x, design, plot = FALSE)
norm_expr <- v$E
pheno_data <- pheno_data[which(pheno_data$condition != "ACD"),]
identical(pheno_data$match, substring(text = colnames(norm_expr),first = 2,last = nchar(x = colnames(norm_expr))))
# [1] TRUE
colnames(norm_expr) <- rownames(pheno_data)

# Variation analysis ----

## variancePartition ----
names(pheno_data)
# [1] "title"         "geo_accession" "age"           "condition"     "participant"   "sex"          
# [7] "match"         "time"  
names(pheno_data)[5] <- "subject"
pheno_data$time <- ifelse(pheno_data$time, "Followup", "Baseline")
## (1|) notation models variables as a random variable rather than a fixed effect. 
# The variancePartition vignette recommends modeling all categorical variables as random variables almost always.
# Interactions of continuous and categorical variables are not random for simplicity.
formula <- ~ as.numeric(age) + (1|subject) + (1|time) + (1|sex) + (1|sex:time) + as.numeric(age):time + as.numeric(age):sex
param <- SnowParam(4, "SOCK", progressbar = TRUE)
varPart <- fitExtractVarPartModel(exprObj = norm_expr, formula = formula, data = pheno_data, BPPARAM = param)
vp <- sortCols(varPart)
pdf(paste(output_dir, "vp_violin_plot.pdf", sep = ""))
plotVarPart(vp)
dev.off()
# Canonical Correlation Analysis
form <- ~ subject + time + as.numeric(age) + sex
C <- canCorPairs(form, pheno_data)
pdf(paste(output_dir, "cca.pdf", sep = ""))
plotCorrMatrix(C)
dev.off()

## standard deviation and average expression ----
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
# ENSG00000186827 ENSG00000186891 ENSG00000160072 ENSG00000041988 ENSG00000279928 ENSG00000228037 
# 0.2133203       0.1724155       0.2875533       0.1133667       0.3892533       0.2761947 
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
prepared_covariates <- pheno_data[,c(6,7)]
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
library(EnsDb.Hsapiens.v86)
genenames <- mapIds(EnsDb.Hsapiens.v86, genenames, "SYMBOL", "GENEID")
genenames <- make.names(ifelse(is.na(genenames), names(genenames), unname(genenames)), unique = T)
data_var$Symbol <- genenames

## Heatmap of all statistical values ----
pdf(paste(output_dir, "statistic_correlations.pdf", sep = ""))
pheatmap(cor(data_var[,-length(names(data_var))]), display_numbers = T, fontsize_number = 4)
dev.off()
## Save ----
write.csv(data_var, paste(output_dir, "Gomez_variation.csv"), row.names = T)
