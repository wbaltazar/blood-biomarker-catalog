## Date: Jul 11 2024

# Removal of poor quality samples in Dusek et al. data, GSE10041

# Load libraries
library(GEOquery)
library(affy)
library(affycoretools)
library(oligo)
library(tidyverse)
library(affyPLM)
library(PCAtools)

## Input: raw CEL files from GSE10041. This directory will be edited
input_dir <- "~/Desktop/work_repo/data/GSE10041_RAW/"
## Output: QC plots
output_plot <- "~/Desktop/work_repo/github/Dusek_study_data/output/"

## Load in metadata
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
pheno_data <- pheno_data[-grep("M", pheno_data$title),] # M samples are not used for longitudinal analyses.
dim(pheno_data)
# [1] 47  4

files <- list.files(input_dir)
length(files)
# [1] 72
files <- files[!(files %in% paste(pheno_data$geo_accession, ".CEL.gz", sep = ""))]
length(files)
# [1] 25
file.remove(paste(input_dir, files, sep = ""))
raw <- ReadAffy(celfile.path = input_dir)
obj <- affy::rma(raw)
norm_expr <- exprs(obj)
colnames(norm_expr) <- str_replace(colnames(norm_expr), ".CEL.gz", replacement = "")
identical(colnames(norm_expr), rownames(pheno_data))
# [1] TRUE
dim(norm_expr)
# [1] 54675    47
## Since sex is not in the metadata, we annotated it according to RPS4Y value.
rps4y_vals <- norm_expr["201909_at",]
pheno_data$sex <- NA
for (i in 1:length(rps4y_vals)) {
  if (rps4y_vals[i] > 9) {
    pheno_data$sex[i] <- "Male"
  } else {
    pheno_data$sex[i] <- "Female"
  }
}
dim(pheno_data)
# [1] 47  5
table(pheno_data$id, pheno_data$sex) # patient 30 has mixed sex samples. Mislabeled.
pheno_data$geo_accession[grep("30", pheno_data$id)] # [1] "GSM253667" "GSM253695" will be removed from analyses.

## Principal component analysis
p <- pca(as.matrix(norm_expr), metadata = pheno_data, center = T, scale = T)
loads <- rownames(p$loadings)
loads <- mapIds(hgu133plus2.db, loads, "SYMBOL", "PROBEID")
loads <- make.names(ifelse(is.na(loads), names(loads), unname(loads)), unique = T)
rownames(p$loadings) <- loads
pdf(file = paste(output_plot, "QC_pca.pdf"))
biplot(p, showLoadings = T, colby = "time", shape = "sex", encircle = T, lab = pheno_data$id, legendPosition = "right")
dev.off()
## We see that patient 50's time 2 sample is removed from the central cluster.

## Check for visual abnormalities
pheno_data[grep(50, pheno_data$id),]
#           title geo_accession time id  sex
# GSM253685 N1-50     GSM253685   N1 50 Male
# GSM253709 N2-50     GSM253709   N2 50 Male
plm <- fitPLM(raw)
which(pheno_data$id == 50) # [1] 23 47
pdf(file = paste(output_plot, "subject_50_pseudoimages.pdf"))
par(mfrow = c(1,2))
image(plm, which = 23) # some spotting
image(plm, which = 47) # very intense
dev.off()
pdf(file = paste(output_plot, "NUSE.pdf"))
NUSE(plm) # Patient 50 array time 2 array stands out. Should be removed.
dev.off()

## Samples to remove from analysis: GSM253709, GSM253667, GSM253695
files <- list.files(input_dir)
files <- files[grep("GSM253709|GSM253667|GSM253695", files)]
print(files)
# [1] "GSM253667.CEL.gz" "GSM253695.CEL.gz" "GSM253709.CEL.gz"
file.remove(paste(input_dir, files, sep = ""))