## Date: January 17 2025

# The following code contains a standard limma analysis for the following study:
# Rusch, Heather L., Jeffrey Robinson, Sijung Yun, Nicole D. Osier, Christiana Martin, Chris R. Brewin,
# and Jessica M. Gill. “Gene Expression Differences in PTSD Are Uniquely Related to the Intrusion Symptom
# Cluster: A Transcriptome-Wide Analysis in Military Service Members.” Brain, Behavior, and Immunity 80
# (August 1, 2019): 904–8. https://doi.org/10.1016/j.bbi.2019.04.039.


# Load libraries ----

library(oligo)
library(tidyverse)
library(limma)
library(affy)
library(affyPLM)
library(affycoretools)
library(hgu133plus2.db)
library(GEOquery)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(cowplot)
library(PCAtools)
library(EnhancedVolcano)
library(RColorBrewer)

## INPUT: Raw .CEL files from GSE81761
input_dir <- "~/Desktop/work_repo/data/GSE81761_RAW/"
## OUTPUT: limma results, plots, figures
output_dir <- "~/Desktop/work_repo/github/Rusch_study_data/output/"

# get study metadata ----
gset <- getGEO("GSE81761", GSEMatrix =TRUE, AnnotGPL=TRUE)
pheno_data_org <- pData(gset$GSE81761_series_matrix.txt.gz)
feature_data <- fData(gset$GSE81761_series_matrix.txt.gz)

# filter out unneccessary columns
head(pheno_data_org)
dim(pheno_data_org)
# [1] 109  48
pheno_data_org <- pheno_data_org %>% 
  filter(`case/control:ch1` == "No PTSD")
dim(pheno_data_org)
# [1] 47 48

pheno_data <- pheno_data_org %>% 
  dplyr::select(title = title, 
                geo_accession = geo_accession,
                source = source_name_ch1,
                age = `age:ch1`,
                ethnicity = `ethnicity:ch1`,
                race = `race:ch1`,
                sex = `Sex:ch1`,
                time = `timepoint:ch1`)
dim(pheno_data)
# [1] 47  8
pheno_data$patient_ID <- str_extract(string = pheno_data$source, pattern = "Subject \\d+")
dim(pheno_data)
# [1] 47  9
table(pheno_data$patient_ID)
 # Subject 10 Subject 108 Subject 117  Subject 16   Subject 2  Subject 25  Subject 37  Subject 38  Subject 46  Subject 49  Subject 50 
 #          2           2           2           1           1           1           2           2           2           2           1 
 # Subject 53  Subject 61  Subject 63  Subject 66  Subject 68  Subject 69  Subject 71  Subject 72  Subject 73  Subject 76  Subject 80 
 #          2           1           2           2           2           2           1           2           2           2           2 
 # Subject 81   Subject 9  Subject 94  Subject 95  Subject 98 
 #          1           2           2           2           2 

## Remove samples where only one timepoint is available
single_timepoint <- pheno_data$geo_accession[which(pheno_data$patient_ID %in% c("Subject 16", "Subject 2", "Subject 25", "Subject 50", "Subject 61", "Subject 71", "Subject 81"))]
files <- list.files(input_dir, full.names = T)
files <- files[-grep(paste(single_timepoint, collapse = "|"), files)]

# Import the CEL data ----
files <- files[grep(paste(pheno_data$geo_accession, collapse = "|"), files)]
length(files) # [1] 40
data <- ReadAffy(filenames = list.files(input_dir, full.names = T)[grep(paste(files, collapse = "|"), list.files(input_dir, full.names = T))])

# Raw data QC ----
# affyPLM: generates probe linear model and allows generation of pseudoimages
Pset <- fitPLM(data, background = F, normalize = F)
pdf(file = paste(output_dir, "RLE_and_NUSE.pdf", sep = ""))
par(mfrow = c(1,2))
RLE(Pset,main = "RLE for Rusch raw data")
NUSE(Pset, main = "NUSE for Rusch raw data")
dev.off()
## The second to last NUSE plot is highly elevated, RLE is low. 10th sample has a high RLE.
pheno_data <- pheno_data[pheno_data$geo_accession %in% str_extract(colnames(exprs(data)), pattern = "GSM\\d+"),]
pheno_data[c(10,39),]
# title geo_accession                   source age    ethnicity                             race  sex      time
# GSM2175202 md186    GSM2175202  Subject 38, timepoint 2  22 Not Hispanic                        Caucasian Male Follow Up
# GSM2175261 md521    GSM2175261 Subject 108, timepoint 2  36 Not Hispanic Native Hawaiian/Pacific Islander Male Follow-Up
# patient_ID
# GSM2175202  Subject 38
# GSM2175261 Subject 108
pdf(file = paste(output_dir, "GSM2175202_GSM2175261_pseudo.pdf", sep = ""))
par(mfrow = c(1,2))
image(Pset, which = 10)
image(Pset, which = 39)
dev.off()

# Background correct and normalize the data into an eSet object
raw <- data
data <- affy::rma(data)

# Get expression estimates
unfiltered_normalized_expression <- exprs(data)
dim(unfiltered_normalized_expression)
# [1] 54675    40

## QC PCAs and pseudoimages ----
colnames(unfiltered_normalized_expression) <- str_extract(colnames(unfiltered_normalized_expression), pattern = "GSM\\d+")
identical(rownames(pheno_data), colnames(unfiltered_normalized_expression)) 
# [1] TRUE
pheno_data$time <- ifelse(pheno_data$time == "Baseline", "Baseline", "Followup")
p <-  pca(unfiltered_normalized_expression, metadata = pheno_data, center = T, scale = T)
pdf(file = paste(output_dir, "pre_qc_pca_plot.pdf", sep = ""))
biplot(p, 
       colby = "patient_ID", shape = "time",
       legendPosition = "right",
       title = "Rusch et al.", subtitle = "lowly-expressed genes not filtered")
dev.off()
## Remove patient 108, timepoint 2 before further analysis.

## Re-normalize the data without this sample.
files <- files[-grep("GSM2175261", files)]
raw <- ReadAffy(filenames = files)
data <- affy::rma(raw)
unfiltered_normalized_expression <- exprs(data)
dim(unfiltered_normalized_expression)
# [1] 54675    39

## Get expression calls using MAS 5.0 ----
calls <- mas5calls.AffyBatch(raw)
probe_pval <- assayData(calls)[["se.exprs"]] ## Returns p-values
dim(probe_pval)
# [1] 54675    39
minSamples <- min(colSums(table(pheno_data$patient_ID, pheno_data$time)))
print(minSamples)
# [1] 20
expressed <- rowSums(probe_pval < 0.05) >= minSamples
normalized_expression <- unfiltered_normalized_expression[expressed,]
dim(normalized_expression)
# [1] 21133    39

colnames(unfiltered_normalized_expression) <- str_extract(colnames(unfiltered_normalized_expression), pattern = "GSM\\d+")
colnames(normalized_expression) <- str_extract(colnames(normalized_expression), pattern = "GSM\\d+")
pheno_data <- pheno_data[pheno_data$geo_accession %in% colnames(unfiltered_normalized_expression),]
dim(pheno_data)
# [1] 39  9

p_unf <-  pca(unfiltered_normalized_expression, metadata = pheno_data, center = T, scale = T)
p <-  pca(normalized_expression, metadata = pheno_data, center = T, scale = T)
pdf(file = paste(output_dir, "unfiltered_qc_pca.pdf", sep = ""))
biplot(p_unf, 
       colby = "patient_ID", shape = "time",
       legendPosition = "right",
       title = "Rusch et al. data QC",
       subtitle = "RMA normalized, unfiltered")
dev.off()
pdf(file = paste(output_dir, "filtered_qc_pca.pdf", sep = ""))
biplot(p, 
       colby = "patient_ID", shape = "time",
       legendPosition = "right", 
       title = "Rusch et al. data QC",
       subtitle = "RMA normalized, filtered")
dev.off()

grep("GSM2175187|GSM2175214", sampleNames(Pset@phenoData))
# [1]  6 18

pheno_data[pheno_data$geo_accession %in% c("GSM2175214", "GSM2175187"),]
#            title geo_accession                  source age    ethnicity      race  sex     time patient_ID
# GSM2175187 md152    GSM2175187 Subject 10, timepoint 2  49 Not Hispanic Caucasian Male Followup Subject 10
# GSM2175214 md245    GSM2175214 Subject 76, timepoint 1  46 Not Hispanic Caucasian Male Baseline Subject 76

pdf(file = paste(output_dir, "GSM2175187_GSM2175214_pseudo.pdf", sep = ""))
par(mfrow = c(1,2))
image(Pset, which = 6)
image(Pset, which = 18)
dev.off()
## Clear strip pattern on the array, plus the variation in PCA, we removed this sample as well.

# For analysis ----
files <- files[-grep("GSM2175214", files)]
data <- ReadAffy(filenames = files)
pheno_data <- pheno_data[pheno_data$geo_accession %in% str_extract(colnames(exprs(data)), pattern = "GSM\\d+"),]
dim(pheno_data)
# [1] 38  9
pheno_data$time[pheno_data$time != "Baseline"] <- rep(x = "Followup", times = length(pheno_data$time[pheno_data$time != "Baseline"]))
table(pheno_data$patient_ID, pheno_data$time)
#             Baseline Followup
# Subject 10         1        1
# Subject 108        1        0
# Subject 117        1        1
# Subject 37         1        1
# Subject 38         1        1
# Subject 46         1        1
# Subject 49         1        1
# Subject 53         1        1
# Subject 63         1        1
# Subject 66         1        1
# Subject 68         1        1
# Subject 69         1        1
# Subject 72         1        1
# Subject 73         1        1
# Subject 76         0        1
# Subject 80         1        1
# Subject 9          1        1
# Subject 94         1        1
# Subject 95         1        1
# Subject 98         1        1


## Get expression calls using MAS 5.0 ----
calls <- mas5calls.AffyBatch(data)
probe_pval <- assayData(calls)[["se.exprs"]] ## Returns p-values
dim(probe_pval)
# [1] 54675    38
minSamples <- min(colSums(table(pheno_data$patient_ID, pheno_data$time)))
print(minSamples)
# [1] 19
expressed <- rowSums(probe_pval < 0.05) >= minSamples

## Normalize the data into an eSet ----
data <- affy::rma(data)

# Get expression estimates
unfiltered_normalized_expression <- exprs(data)
dim(unfiltered_normalized_expression)
# [1] 54675    38
normalized_expression <- unfiltered_normalized_expression[expressed,]
dim(normalized_expression)
# [1] 21374    38

colnames(unfiltered_normalized_expression) <- str_extract(colnames(unfiltered_normalized_expression), pattern = "GSM\\d+")
colnames(normalized_expression) <- str_extract(colnames(normalized_expression), pattern = "GSM\\d+")
pheno_data <- pheno_data[pheno_data$geo_accession %in% colnames(normalized_expression),]
dim(pheno_data)
# [1] 38  9
p <-  pca(unfiltered_normalized_expression, metadata = pheno_data, center = T, scale = T)
pdf(file = paste(output_dir, "unfiltered_pca.pdf", sep = ""))
biplot(p, 
       colby = "patient_ID", 
       legendPosition = "right",
       title = "Rusch et al. data pre-QC",
       subtitle = "lowly-expressed genes not filtered")
dev.off()
pheno_data[pheno_data$geo_accession == "GSM2175187",]
p <-  pca(normalized_expression, metadata = pheno_data, center = T, scale = T)
pdf(file = paste(output_dir, "filtered_pca.pdf", sep = ""))
biplot(p, 
       colby = "patient_ID", 
       legendPosition = "right",
       title = "Rusch et al. data pre-QC",
       subtitle = "lowly-expressed genes filtered")
dev.off()
pheno_data[pheno_data$geo_accession == "GSM2175187",]
# GSM2175187 Subject 10, timepoint 2
which(pheno_data$geo_accession == "GSM2175187") # [1] 6
pdf(file = paste(output_dir, "GSM2175187_pseudo.pdf", sep = ""))
image(Pset, which = 6)
dev.off()

# Pre-processing ----
df <- data.frame(normalized_expression)
identical(colnames(df), pheno_data$geo_accession)
# [1] TRUE
dim(df)
# [1] 21374    38
dim(pheno_data)
# [1] 38  9

# QC ----

# Generate informative plots / data
pdf(file = paste(output_dir, "features_hist.pdf", sep =""))
hist(normalized_expression)
dev.off()
pdf(file = paste(output_dir, "samples_box.pdf", sep = ""))
boxplot(normalized_expression)
dev.off()
quantile(normalized_expression)
# 0%       25%       50%       75%      100% 
# 1.974886  6.052111  7.206498  8.371418 15.494129 

# What are the highly expressed RNAs?
avg_gene_exp <- rowMeans(x = normalized_expression, na.rm = T)
head(avg_gene_exp) # creates vector of each gene's mean
avg_gene_exp <- avg_gene_exp[order(avg_gene_exp, decreasing = T)] # order by decreasing order
keytypes(hgu133plus2.db) # see keys for affymetrix chip
mapIds(x = hgu133plus2.db, keys = head(names(avg_gene_exp)), keytype = "PROBEID", column = "SYMBOL")

# 214414_x_at 211696_x_at 211699_x_at 204018_x_at 217414_x_at 211745_x_at 
# "HBA1"       "HBB"      "HBA1"      "HBA1"      "HBA2"      "HBA1" 

# PCA and MDS plots ----
library(PCAtools)
p <- pca(df, metadata = pheno_data, center = T, scale = T)

scree <- screeplot(p, components = getComponents(p, 1:20))
ptime <- biplot(p, 
                showLoadings = F,
                lab = NULL, 
                colby = "time", colkey = c("Baseline" = "blue", "Followup" = "red"), 
                legendPosition = 'right', legendLabSize = 12, legendIconSize = 3,
                title = 'Rusch 2 months',
                subtitle = 'PC1 versus PC2')
prac <- biplot(p, 
               showLoadings = F,
               lab = NULL, 
               colby = "race", 
               legendPosition = 'right', legendLabSize = 7, legendIconSize = 3)
psub <- biplot(p, 
               showLoadings = F,
               lab = NULL, 
               colby = "patient_ID",
               legendPosition = 'right', legendLabSize = 9, legendIconSize = 3)
peigencor <- eigencorplot(p,
                          components = getComponents(p, 1:10),
                          metavars = c("time", 'ethnicity', 'race'),
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
toprow <- plot_grid(ptime, prac, psub, align = 'h', nrow = 1)
botrow <- plot_grid(scree, peigencor, align = 'h', nrow = 1)
pdf(file = paste(output_dir, "pca.pdf", sep = ""), height = 12, width = 16)
plot_grid(toprow, botrow, ncol = 1, align = 'v')
dev.off()

# limma analysis ----

# Remove spike-in probes (AFFX hybridization and polyA probes)
control_index <- grep("lys|Lys|thr|Thr|dap|Dap|phe|Phe|Trp|bioB|BioB|BioC|bioC|bioD|BioD|cre|Cre", x = row.names(df))
rownames(df)[control_index]
df <- df[-control_index,]
dim(df)
# [1] 21332    38
rownames(df)[grep("AFFX", x = rownames(df))]

id <- factor(pheno_data$patient_ID)
time <- factor(pheno_data$time)
ethnicity <- factor(pheno_data$ethnicity)
race <- factor(pheno_data$race)
TS <- paste(time,ethnicity,race, sep = ".")
TS <- factor(TS)
design <- model.matrix(~0+TS, df)
colnames(design) <- make.names(levels(TS))
cor <- duplicateCorrelation(df, design, block=id)
cor$consensus.correlation
# [1] 0.5508815

fit <- lmFit(object=df, design=design, block=id, correlation=cor$consensus.correlation)

# contrasts
my_contrasts <- makeContrasts(Time = (Followup.Hispanic.Caucasian + Followup.Hispanic.Other.Unknown + 
                                         Followup.Not.Hispanic.Black + Followup.Not.Hispanic.Caucasian + 
                                         Followup.Not.Hispanic.Native.American)/5 - 
                                        (Baseline.Hispanic.Caucasian + Baseline.Hispanic.Other.Unknown +
                                        Baseline.Not.Hispanic.Black + Baseline.Not.Hispanic.Caucasian + 
                                        Baseline.Not.Hispanic.Native.American)/5,
                              ITimeHispanic = ((Followup.Hispanic.Caucasian + Followup.Hispanic.Other.Unknown)/2 -
                                                (Baseline.Hispanic.Caucasian + Baseline.Hispanic.Other.Unknown)/2) -
                                ((Followup.Not.Hispanic.Black + Followup.Not.Hispanic.Caucasian + 
                                   Followup.Not.Hispanic.Native.American)/3 - 
                                   (Baseline.Not.Hispanic.Black + Baseline.Not.Hispanic.Caucasian + 
                                      Baseline.Not.Hispanic.Native.American)/3),
                              TimeWhite = ((Followup.Hispanic.Caucasian + Followup.Not.Hispanic.Caucasian)/2 -
                                (Baseline.Hispanic.Caucasian + Baseline.Not.Hispanic.Caucasian)/2),
                              TimeBlack = (Followup.Not.Hispanic.Black - Baseline.Not.Hispanic.Black),
                              TimeNativeAmerican = (Followup.Not.Hispanic.Native.American - Baseline.Not.Hispanic.Native.American),
                              levels = colnames(design))
fit2 <- contrasts.fit(fit, my_contrasts)
x <- eBayes(fit2)

y <- topTable(x, number = Inf, adjust.method = "fdr", sort.by = "F")
y$Symbol <- mapIds(x = hgu133plus2.db, keys = row.names(y), column = "SYMBOL", keytype = "PROBEID")

dt <- decideTests(x)
summary(dt)
#         Time ITimeHispanic TimeWhite TimeBlack TimeNativeAmerican
# Down     199             0         0         0                  0
# NotSig 18332         21332     21332     21332              21332
# Up      2801             0         0         0                  0
dt <- decideTests(x, adjust.method = 'none')
summary(dt)
#         Time ITimeHispanic TimeWhite TimeBlack TimeNativeAmerican
# Down     983           299       640       135                324
# NotSig 14309         19487     20213     20983              20829
# Up      6040          1546       479       214                179

write.csv(y, paste(output_dir, "limma_F.csv", sep = ""), row.names=T)
names(y)
tables <- list()
for (i in 1:5) {
  tables[[i]] <- topTable(x, coef = i, number = Inf, adjust.method = "fdr")
  tables[[i]]$Symbol<-mapIds(x = hgu133plus2.db, keys = row.names(tables[[i]]), column = "SYMBOL", keytype = "PROBEID")
  write.csv(tables[[i]], file = paste(output_dir, "limma_", names(y)[i],".csv", sep = ""))
}
## Plot p-value histograms ----
for (i in 1:5) {
  a <- ggplot(data = tables[[i]], aes(x = P.Value)) +
    geom_histogram(fill = "lightblue", color = "black") +
    theme_cowplot() +
    ggtitle("", subtitle = names(y)[i]) +
    theme(title = element_text(size = 7.5))
  b <- ggplot(data = tables[[i]], aes(x = adj.P.Val)) +
    geom_histogram(fill = "lightblue", color = "black") +
    theme_cowplot() +
    ggtitle("", subtitle = names(y)[i]) +
    theme(title = element_text(size = 7.5))
  pdf(file = paste(output_dir, names(y)[i], "_p_value_histograms.pdf", sep = ""), width = 10)
  plot_grid(a, b, labels = c('A. P-value histogram', 'B. BH-adjusted'), nrow = 1)
  dev.off()
}

## Create volcano plots ----
library(EnhancedVolcano)
volcanos <- list()
for (i in 1:5) {
  volcanos[[i]] <- EnhancedVolcano(tables[[1]],
                              lab = tables[[1]][,7],
                              x = "logFC",
                              y = "P.Value",
                              title = ifelse(i == 1, "Baseline vs. 3mo", ""),
                              subtitle = names(y)[i],
                              caption = "baseline vs. followup", captionLabSize = 5,
                              col = c("black","lightblue","royalblue","blue"),
                              legendPosition = 'bottom', legendLabSize = 6.5, legendIconSize = 3,
                              pCutoff = 1e-2,
                              FCcutoff = 1,
                              labSize = (4.0),
                              colAlpha = 1,
                              pointSize = 1.0,
                              drawConnectors = TRUE,
                              max.overlaps = 20)
}
pdf(paste(output_dir, "volcanos.pdf", sep = ""), width = 12, height = 10)
plot_grid(plotlist = volcanos, nrow = 2, align = 'h')
dev.off()