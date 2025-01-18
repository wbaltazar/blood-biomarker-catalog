## Date: Jan 17 2025

## Calculate statistics for stability in LaRocca study
## INPUT: Expr and pheno data
## OUTPUT: .csv file of heritability, standard deviation, and variancePartition statistics

### SUMMARY: How stable are gene expression measurements after 16 weeks of cardiovascular training?

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

## INPUT: Raw count file from GSE206505 on Gene Expression Omnibus ----
input_dir <- "~/Desktop/work_repo/data/GSE206505_VO2-counts-matrix.xlsx"
## OUTPUT: tables and graphics ----
output_dir <- "~/Desktop/work_repo/github/LaRocca_study_data/output/"

# Load in the data ----
gset <- getGEO(GEO = "GSE206505")
p <- phenoData(gset$GSE206505_series_matrix.txt.gz)
pheno_data_raw <- p@data
pheno_data <- pheno_data_raw %>% 
  dplyr::select(title = title,
                geo_accession = geo_accession,
                time = `pre or post:ch1`,
                sex = `Sex:ch1`,
                response = `r or nr:ch1`)
dim(pheno_data)
# [1] 60  5
pheno_data$participant <- rep(1:30, each = 2)
pheno_data$participant <- as.factor(pheno_data$participant)
x1 <- readxl::read_xlsx("./data/GSE206505_VO2-counts-matrix.xlsx")
#### x1 saved in case of re-run
x <- x1
dim(x)
# [1] 27665    61
x <- data.frame(x)
rownames(x) <- make.names(x[,1], unique = TRUE)
x <- x[,-1]
# Sample 21, paired with 22, has no values. "Subject 11". Remove sample 21, but not sample 22.
x <- x[, -21]
x <- DGEList(x)
geneid <- rownames(x) 
library(EnsDb.Hsapiens.v86)
genes <- mapIds(EnsDb.Hsapiens.v86, keys=geneid, column=c("GENEID"), 
                keytype="SYMBOL") # Get ENSEMBL IDs.
x$genes <- genes
detach("package:EnsDb.Hsapiens.v86", unload = TRUE) ## allows dplyr functions to run properly.
identical(str_replace(colnames(x$counts), "X", ""), pheno_data[-grep(21, pheno_data$title),"title"])
# [1] TRUE
pheno_data <- pheno_data[-grep(21, pheno_data$title),]
colnames(x$counts) <- pheno_data$geo_accession
x$samples$group <- pheno_data$response
x$samples$time <- pheno_data$time
x$samples$sex <- pheno_data$sex
x$samples$participant <- pheno_data$participant
x$samples$geo <- pheno_data$geo_accession
keep.exprs <- filterByExpr(x, group=x$samples$time) 
x$counts <- x$counts[keep.exprs,] 
dim(x)
# [1] 15555    59
x <- calcNormFactors(x, method = "TMM")
id <- factor(x$samples$participant)
time <- factor(x$samples$time)
response <- factor(x$samples$group)
TS <- paste(response,time, sep = ".")
TS <- factor(TS)
design <- model.matrix(~0+TS)
colnames(design) <- make.names(levels(TS))
v <- voom(x, design, plot = FALSE)
norm_expr <- v$E
identical(rownames(pheno_data), colnames(norm_expr))
# [1] TRUE

# Variation analysis ----

## varianceParitition ----
names(pheno_data)
# [1] "title"         "geo_accession" "time"          "sex"           "response"      "participant"  
names(pheno_data)[6] <- "subject"
formula <- ~ (1|time) + (1|subject) + (1|response) + (1|time:response)
param <- SnowParam(4, "SOCK", progressbar = TRUE)
varPart <- fitExtractVarPartModel(exprObj = norm_expr, formula = formula, data = pheno_data, BPPARAM = param)
vp <- sortCols(varPart)
pdf(paste(output_dir, "vp_violin_plot.pdf", sep = ""))
plotVarPart(vp)
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
# A1BG  A1BG.AS1       A2M   A2M.AS1      AAAS      AACS 
# 0.5353166 0.3085213 0.6668374 0.5656673 0.2658260 0.2002768 
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
prepared_covariates <- pheno_data[,c(3,5)]
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
data_var$Symbol <- rownames(data_var)

## Heatmap ----
pdf(paste(output_dir, "statistic_heatmap.pdf", sep = ""))
pheatmap(cor(data_var[,-length(names(data_var))]), display_numbers = T, fontsize_number = 7)
dev.off()
## Save ----
write.csv(data_var, paste(output_dir, "LaRocca_variation.csv", sep = ""), row.names = T)