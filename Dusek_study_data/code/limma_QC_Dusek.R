## Date: Jan 17 2025

## This script performs limma differential expression analysis of data from the following study:
# Dusek, Jeffery A., Hasan H. Otu, Ann L. Wohlhueter, Manoj Bhasin, Luiz F. Zerbini, Marie G. Joseph, 
# Herbert Benson, and Towia A. Libermann. “Genomic Counter-Stress Changes Induced by the Relaxation Response.” 
# PLoS ONE 3, no. 7 (July 2, 2008): e2576. https://doi.org/10.1371/journal.pone.0002576.

# Load libraries ----

library(oligo)
library(tidyverse)
library(limma)
library(affy)
library(affycoretools)
library(hgu133plus2.db)
library(GEOquery)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(cowplot)
library(PCAtools)
library(EnhancedVolcano)
library(RColorBrewer)

## INPUT: Raw, unzipped .CEL files
input_dir <- "~/Desktop/work_repo/data/GSE10041_RAW/"
## OUTPUT: limma results, plots, figures
output_dir <- "~/Desktop/work_repo/github/Dusek_study_data/output/"

## Load in metadata from GEO accession: GSE10041 ----
gset <- getGEO("GSE10041", GSEMatrix =TRUE, AnnotGPL=TRUE)
pheno_data_org <- pData(gset$GSE10041_series_matrix.txt.gz)
# feature_data <- fData(gset$GSE10041_series_matrix.txt.gz)
pheno_data <- pheno_data_org %>% 
  dplyr::select(title = title,
                geo_accession = geo_accession)
pheno_data$time <- str_extract(pheno_data_org$title, "N.")
pheno_data$id <- str_extract(pheno_data_org$title, "-.*") %>% 
  str_replace("-","")
dim(pheno_data)
# [1] 72 4
length(grep("M", pheno_data$title)) # [1] 25 This is the correct number in the original and validation groups.
# The M samples were sampled from experienced meditators at a single timepoint. This is not needed for longitudinal analyses.
pheno_data <- pheno_data[-grep("M", pheno_data$title),]
dim(pheno_data)
# [1] 47  4

## Read in the expression data ----
files <- list.files(input_dir, full.names = T)
length(files)
# [1] 72
files <- files[grep(paste(pheno_data$geo_accession, collapse = "|"), files)]
length(files)
# [1] 47
raw <- ReadAffy(filenames = files)
obj <- affy::rma(raw) # Normalizes the data via the RMA method.
norm_expr <- exprs(obj)
colnames(norm_expr) <- str_replace(colnames(norm_expr), ".CEL.gz", replacement = "")
identical(colnames(norm_expr), rownames(pheno_data))
# [1] TRUE
dim(norm_expr)
# [1] 54675    47
## Since sex was not included in the metadata, we will infer it using RPS4Y1 transcript values
rps4y1_vals <- norm_expr["201909_at",]
pheno_data$sex <- NA
for (i in 1:length(rps4y1_vals)) {
  if (rps4y1_vals[i] > 8) {
    pheno_data$sex[i] <- "Male"
  } else {
    pheno_data$sex[i] <- "Female"
  }
}
dim(pheno_data)
# [1] 47  5
table(pheno_data$id, pheno_data$sex) # Patient 30's samples have mixed sex samples. Most likely due to a labeling error.
pheno_data[grep("30", pheno_data$id),] # [1] "GSM253667" "GSM253695" will be removed from analyses.
#           title geo_accession time id    sex
# GSM253667 N1-30     GSM253667   N1 30 Female
# GSM253695 N2-30     GSM253695   N2 30   Male
grep("30", pheno_data$id)
# [1]  5 33
pdf(file = paste(output_dir, "sample_rps4y1_values.pdf", sep = ""))
rps4y1_vals %>% data.frame() %>% bind_cols(pheno_data$id) %>% 
  ggplot(aes(x = 1:length(rps4y1_vals), y = .)) + 
  geom_text(aes(label = `...2`, color = `...2`)) +
  labs(x = 'index', y = 'RPS4Y1 normalized expression') +
  geom_hline(yintercept = 8, linetype = 'dashed')
dev.off()

## Quality control ----
## AffyPLM build probe linear models that are useful for assessing data quality.
library(affyPLM)
## Check for visual abnormalities
plm <- fitPLM(raw, background = F, normalize = F)
pdf(file = paste(output_dir, "NUSE_RLE_preQC.pdf", sep = ""))
par(mfrow = c(1,2))
RLE(plm, ylim = c(-4,4))
title("RLE")
NUSE(plm, ylim = c(0.75,1.8)) # Patient 50 array time 2 array stands out.
title("NUSE")
dev.off()
pheno_data[nrow(pheno_data),]
# title geo_accession time id  sex
# GSM253709 N2-50     GSM253709   N2 50 Male

# Last sample has artifacts
which(pheno_data$id == 50) # [1] 23 47
pdf(file = paste(output_dir, "subject_50_pseudoimages.pdf", sep = ""))
par(mfrow = c(1,2))
image(plm, which = 23) # some spotting
image(plm, which = 47) # very intense
dev.off()

pdf(file = paste(output_dir, "GSM253709_pseudo.pdf", sep = ""))
image(plm, which = 47)
dev.off()

## Principal component analysis of filtered, normalized expression
p <- pca(as.matrix(norm_expr), metadata = pheno_data, center = T, scale = T)
loads <- rownames(p$loadings)
loads <- mapIds(hgu133plus2.db, loads, "SYMBOL", "PROBEID")
loads <- make.names(ifelse(is.na(loads), names(loads), unname(loads)), unique = T)
rownames(p$loadings) <- loads
pdf(file = paste(output_dir, "filter_preQC_pca.pdf", sep = ""))
biplot(p, showLoadings = T, colby = "time", shape = "sex", encircle = T, lab = pheno_data$id, legendPosition = "right")
dev.off()

# Filter unexpressed genes
calls <- mas5calls.AffyBatch(raw)
probe_pval <- assayData(calls)[["se.exprs"]] ## Returns p-values
dim(probe_pval)
# [1] 54675    47
minSamples <- min(colSums(table(pheno_data$id, pheno_data$time)))
print(minSamples)
# [1] 23
expressed <- rowSums(probe_pval < 0.05) >= minSamples
norm_expr <- norm_expr[expressed,]
dim(norm_expr)
# [1] 15938    47
## Principal component analysis of unfiltered, normalized expression
norm_expr <- exprs(obj)
colnames(norm_expr) <- str_replace(colnames(norm_expr), ".CEL.gz", replacement = "")
p <- pca(as.matrix(norm_expr), metadata = pheno_data, center = T, scale = T)
loads <- rownames(p$loadings)
loads <- mapIds(hgu133plus2.db, loads, "SYMBOL", "PROBEID")
loads <- make.names(ifelse(is.na(loads), names(loads), unname(loads)), unique = T)
rownames(p$loadings) <- loads
pdf(file = paste(output_dir, "NOfilter_preQC_pca.pdf", sep = ""))
biplot(p, showLoadings = T, colby = "time", shape = "sex", encircle = T, lab = pheno_data$id, legendPosition = "right")
dev.off()

## In this analysis, we removed GSM253709, GSM253667, GSM253695. Full explanations can be found in Supplementary
## Data File 1 of our paper.

# Import the CEL data for limma analysis ----
files <- files[-grep("GSM253709|GSM253667|GSM253695", files)]
pheno_data <- pheno_data[-grep("GSM253709|GSM253667|GSM253695", pheno_data$geo_accession),]
data <- ReadAffy(filenames = files)

# Get expression calls using MAS 5.0
calls <- mas5calls.AffyBatch(data)
probe_pval <- assayData(calls)[["se.exprs"]] ## Returns p-values
dim(probe_pval)
# [1] 54675    44
minSamples <- min(table(pheno_data$time))
print(minSamples)
# [1] 22
expressed <- rowSums(probe_pval < 0.05) >= minSamples

## Normalize the data into an eSet ----
data <- affy::rma(data)

## Get expression estimates
normalized_expression <- exprs(data)
dim(normalized_expression)
# [1] 54675    44
## Remove unexpressed probes
normalized_expression <- normalized_expression[expressed,]
dim(normalized_expression)
# [1] 15963    44

# Generate informative plots / data
pdf(file = paste(output_dir, "features_hist.pdf", sep = ""))
hist(normalized_expression)
dev.off()
pdf(file = paste(output_dir, "samples_box.pdf", sep = ""))
boxplot(normalized_expression)
dev.off()
quantile(normalized_expression)
# 0%       25%       50%       75%      100% 
# 2.513918  5.922759  7.061899  8.206119 14.765905 

# What are the highly expressed RNAs?
avg_gene_exp <- rowMeans(x = normalized_expression, na.rm = TRUE)
head(avg_gene_exp) # creates vector of each gene's mean
avg_gene_exp <- avg_gene_exp[order(avg_gene_exp, decreasing = TRUE)] # order by decreasing order
keytypes(hgu133plus2.db) # see keys for affymetrix chip
mapIds(x = hgu133plus2.db, keys = head(names(avg_gene_exp)), keytype = "PROBEID", column = "SYMBOL")

# 214414_x_at 211745_x_at 211696_x_at 209458_x_at 204018_x_at 211699_x_at 
# "HBA1"      "HBA1"       "HBB"      "HBA1"      "HBA1"      "HBA1" 

# rename expression data columns
colnames(normalized_expression)
colnames(normalized_expression) <- gsub(pattern = "(.CEL.gz)", replacement = "", x = colnames(normalized_expression))
identical(pheno_data$geo_accession, colnames(normalized_expression))
# [1] TRUE
df <- data.frame(normalized_expression)

# Annotate sex
rps4y1_vals <- norm_expr["201909_at",]
pheno_data$sex <- NA
for (i in 1:length(rps4y1_vals)) {
  if (rps4y1_vals[i] > 8) {
    pheno_data$sex[i] <- "Male"
  } else {
    pheno_data$sex[i] <- "Female"
  }
}

table(pheno_data$id, pheno_data$sex)
#    Female Male
# 01      0    2
# 02      0    2
# 03      0    2
# 06      0    2
# 07      2    0
# 11      2    0
# 12      2    0
# 14      2    0
# 23      0    2
# 32      2    0
# 33      2    0
# 34      2    0
# 35      2    0
# 36      0    1
# 39      0    2
# 40      0    2
# 41      0    2
# 46      2    0
# 47      0    2
# 49      0    2
# 50      0    1
# 51      2    0
# 59      2    0

# Run MDS and PCA ----
## PCAtools
p <- pca(df, metadata = pheno_data, center = T, scale = T)

scree <- screeplot(p, components = getComponents(p, 1:20))
ptime <- biplot(p, 
                showLoadings = F,
                lab = NULL, 
                colby = "time", colkey = c("N1" = "green", "N2" = "orange2"), 
                legendPosition = 'right', legendLabSize = 12, legendIconSize = 3,
                title = 'Dusek 1 month',
                subtitle = 'PC1 versus PC2')
psex <- biplot(p, 
               showLoadings = F,
               lab = NULL, 
               colby = "sex", 
               legendPosition = 'right', legendLabSize = 12, legendIconSize = 3)
psub <- biplot(p, 
               showLoadings = F,
               lab = NULL, 
               colby = "id",
               legendPosition = 'right', legendLabSize = 9, legendIconSize = 3)
peigencor <- eigencorplot(p,
                          components = getComponents(p, 1:10),
                          metavars = c("time", 'sex'),
                          cexCorval = 1.0,
                          fontCorval = 2,
                          posLab = 'all',
                          rotLabX = 45,
                          scale = TRUE,
                          colCorval = "burlywood1",
                          main = "PC trait correlates",
                          cexMain = 1.5,
                          plotRsquared = FALSE,
                          corFUN = 'pearson',
                          corUSE = 'pairwise.complete.obs',
                          signifSymbols = c('****', '***', '**', '*', ''),
                          signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                          returnPlot = F)
toprow <- plot_grid(ptime, psex, psub, align = 'h', nrow = 1)
botrow <- plot_grid(scree, peigencor, align = 'h', nrow = 1)
pdf(file = paste(output_dir, "pca.pdf", sep = ""), height = 12, width = 16)
plot_grid(toprow, botrow, ncol = 1, align = 'v')
dev.off()

# limma analysis ----

# Remove spike-in control probes (AFFX hybridization and polyA probes)
control_index <- grep("lys|Lys|thr|Thr|dap|Dap|phe|Phe|Trp|bioB|BioB|BioC|bioC|bioD|BioD|cre|Cre", x = row.names(df))
rownames(df)[control_index]
df <- df[-control_index,]
dim(df)
# [1] 15945    44
rownames(df)[grep("AFFX", x = rownames(df))]
# [1] "AFFX-HSAC07/X00351_3_at"     "AFFX-HSAC07/X00351_5_at"     "AFFX-HSAC07/X00351_M_at"     "AFFX-hum_alu_at"            
# [5] "AFFX-HUMGAPDH/M33197_3_at"   "AFFX-HUMGAPDH/M33197_5_at"   "AFFX-HUMGAPDH/M33197_M_at"   "AFFX-HUMISGF3A/M97935_3_at" 
# [9] "AFFX-HUMISGF3A/M97935_MB_at" "AFFX-HUMRGE/M10098_3_at"     "AFFX-HUMRGE/M10098_5_at"     "AFFX-HUMRGE/M10098_M_at"    
# [13] "AFFX-M27830_5_at"    

id <- factor(pheno_data$id)
time <- factor(pheno_data$time)
sex <- factor(pheno_data$sex)
TS <- factor(paste(time,sex, sep = "."))
design <- model.matrix(~0+TS, df)
colnames(design) <- make.names(levels(TS))
cor <- duplicateCorrelation(df, design, block=id)
cor$consensus.correlation
# [1] 0.3508609

fit <- lmFit(object=df, design=design, block=id, correlation=cor$consensus.correlation)

# contrasts
contrasts <- makeContrasts(Time = (N2.Female + N2.Male)/2 - (N1.Male + N1.Female)/2,
                                   Time_sex = (N2.Male - N1.Male) - (N2.Female - N1.Female),
                                   levels = colnames(design))


fit2 <- contrasts.fit(fit, contrasts)
x <- eBayes(fit2)

# Download results and visualizations ----
y <- topTable(x, number = Inf, adjust.method = "fdr", sort.by = "F")
y$Symbol <- mapIds(x = hgu133plus2.db, keys = row.names(y), column = "SYMBOL", keytype = "PROBEID")
write.csv(y, file = paste(output_dir, "limma_F.csv", sep = ""), row.names = T)

dt <- decideTests(x)
summary(dt)
#         Time Time_sex
# Down       0        0
# NotSig 15945    15945
# Up         0        0
dt <- decideTests(x, adjust.method = "none")
summary(dt)
#         Time Time_sex
# Down      87      470
# NotSig 15684    14708
# Up       174      767
tables <- list(time = topTable(x, coef = 1, number = Inf, adjust.method = "fdr"),
               interaction = topTable(x, coef = 2, number = Inf, adjust.method = "fdr"))
tables <- lapply(tables, function(x) {
  Symbols <- mapIds(x = hgu133plus2.db, keys = row.names(x), column = "SYMBOL", keytype = "PROBEID")
  x <- cbind(x, Symbols)
})
for (i in 1:length(tables)) {
  write.csv(tables[i], file = paste(output_dir, "limma_",names(tables)[i],".csv", sep = ""), row.names = T)
}

## Generate p-value histograms----
a <- ggplot(data = tables$time, aes(x = P.Value)) +
  geom_histogram(fill = "lightblue", color = "black") +
  ggtitle(paste("P-value histogram time")) +
  theme_cowplot() +
  theme(title = element_text(size = 7.5))
b <- ggplot(data = tables$interaction, aes(x = P.Value)) +
  geom_histogram(fill = "lightblue", color = "black") +
  ggtitle(paste("P-value histogram interaction")) +
  theme_cowplot() +
  theme(title = element_text(size = 7.5))
pdf(file = paste(output_dir, "p_value_histograms.pdf", sep = ""))
plot_grid(a,b, labels = c('A', 'B'))
dev.off()
a <- ggplot(data = tables$time, aes(x = adj.P.Val)) +
  geom_histogram(fill = "lightblue", color = "black") +
  ggtitle(paste("Adj. P-value histogram time")) +
  theme_cowplot() +
  theme(title = element_text(size = 7.5))
b <- ggplot(data = tables$interaction, aes(x = adj.P.Val)) +
  geom_histogram(fill = "lightblue", color = "black") +
  ggtitle(paste("Adj. P-value histogram interaction")) +
  theme_cowplot() +
  theme(title = element_text(size = 7.5))
pdf(file = paste(output_dir, "adj_p_value_histograms.pdf", sep = ""))
plot_grid(a,b, labels = c('A', 'B'))
dev.off()

## EnhancedVolcano plots----
volcano1 <- EnhancedVolcano(tables$time,
                            lab = tables$time$Symbols,
                            x = "logFC",
                            y = "P.Value",
                            title = "Meditiation Study: 3 month comparison", titleLabSize = 8,
                            subtitle = "time main effect", subtitleLabSize = 7,
                            caption = "",
                            col = c("black","lightblue","royalblue","blue"),
                            legendPosition = 'bottom', legendLabSize = 6.5, legendIconSize = 3,
                            pCutoff = 1e-2,
                            FCcutoff = 1,
                            labSize = (4.0),
                            colAlpha = 1,
                            pointSize = 1.0,
                            drawConnectors = TRUE,
                            max.overlaps = 50)
volcano2 <- EnhancedVolcano(tables$interaction,
                            lab = tables$interaction$Symbols,
                            x = "logFC",
                            y = "P.Value",
                            title = "", 
                            subtitle = "time interaction effect (male-female)", subtitleLabSize = 7,
                            caption = "",
                            col = c("black","lightblue","royalblue","blue"),
                            legendPosition = 'bottom', legendLabSize = 6.5, legendIconSize = 3,
                            pCutoff = 1e-2,
                            FCcutoff = 1,
                            labSize = (4.0),
                            colAlpha = 1,
                            pointSize = 1.0,
                            drawConnectors = TRUE,
                            max.overlaps = 50)
pdf(paste(output_dir, "volcanos.pdf", sep = ""))
plot_grid(volcano1, volcano2, labels = c('A', 'B'))
dev.off()