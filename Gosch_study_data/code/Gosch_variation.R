## Date: Jul 11 2024

## Calculate statistics for stability in Gosch study
## INPUT: Raw count file from GSE223613
## OUTPUT: .csv file of heritability, standard deviation, average expression, and variancePartition statistics

### SUMMARY: What does a person's gene expression signature in blood say about the time of day?

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
input_dir <- "~/Desktop/work_repo/data/GSE223613_counts.txt"
## OUTPUT: tables and graphics ----
output_dir <- "~/Desktop/work_repo/github/Gosch_study_data/output/"

# Load in the data ----
gset <- getGEO(GEO = "GSE223613")
p <- phenoData(gset$GSE223613_series_matrix.txt.gz)
pheno_data_raw <- p@data
pheno_data <- pheno_data_raw %>% 
  dplyr::select(title = title,
                geo_accession = geo_accession,
                age = `age:ch1`,
                participant = `participant:ch1`,
                time = `sampling time:ch1`,
                sex = `Sex:ch1`,
                sample = description)
library(edgeR)
library(limma)
setwd("~/Desktop/work_repo")
x1 <- read.delim("./data/GSE223613_counts.txt")
#### x1 saved in case code must be re-run
x <- x1
dim(x)
# [1] 33121    81
gene_list <- x[,1]
rownames(x) <- make.names(gene_list, unique = TRUE)
x <- x[,-1]
match(pheno_data$sample, sub(".","",sub("\\.","-",names(x))))
names(x) <- pheno_data$title
x <- DGEList(x)
x[["samples"]]
x$samples$group <- as.factor(pheno_data$participant)
x$samples$time <- as.factor(pheno_data$time)
x$samples$sex <- as.factor(pheno_data$sex)

## Update May 14, 2024: Remove Patient J at 23 hours. Library size: 401784615
which(x$samples$lib.size == 401784615)
# [1] 60
colnames(x$counts)[60] == rownames(x$samples)[60]
# [1] TRUE
x$counts <- x$counts[,-60]
x$samples <- x$samples[-60,]
dim(x$counts)
# [1] 33121    79
dim(x$samples)
# [1] 79  5
keep.exprs <- filterByExpr(x, group=x$samples$group) 
x <- x[keep.exprs,, keep.lib.sizes=FALSE] 
dim(x)
# [1] 16389    79
x <- calcNormFactors(x, method = "TMM")
id <- factor(x$samples$group)
time <- factor(x$samples$time)
levels(time) <- c("02.00h", "05.00h", "08.00h", "11.00h", "14.00h", "17.00h", "20.00h", "23.00h")
sex <- factor(x$samples$sex)
TS <- paste(sex,time, sep = ".")
TS <- factor(TS)
design <- model.matrix(~0+TS)
colnames(design) <- make.names(levels(TS))
v <- voom(x, design, plot = FALSE)
norm_expr <- v$E
pheno_data <- pheno_data[-60,]
identical(pheno_data$title, colnames(norm_expr))
# [1] TRUE

# Variation analysis ----

## varianceParitition ----
names(pheno_data)
# [1] "title"         "geo_accession" "age"           "participant"   "time"          "sex"          
# [7] "sample"  
## (1|) notation models variables as a random variable rather than a fixed effect. 
# The variancePartition vignette recommends modeling all categorical variables as random variables almost always.
# Interactions of continuous and categorical variables are not random for simplicity.
names(pheno_data)[4] <- "subject"
formula <- ~ as.numeric(age) + (1|subject) + (1|time) + (1|sex) + as.numeric(age):sex + as.numeric(age):time + (1|sex:time)
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
# DDX11L1       WASH7P    MIR6859.1 LOC101927589    LOC729737 LOC100996442 
# 0.7878850    0.3280794    0.2687201    0.3911288    0.2454477    0.6268199 
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
prepared_covariates <- pheno_data[,c(-1,-2,-4,-7)]
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
data_var$Symbol <- rownames(data_var)

## Heatmap ----
pdf(paste(output_dir, "statistic_correlations.pdf", sep = ""))
pheatmap(cor(data_var[,-length(names(data_var))]), display_numbers = T)
dev.off()
## Save ----
write.csv(data_var, paste(output_dir, "Gosch_variation.csv", sep = ""), row.names = T)