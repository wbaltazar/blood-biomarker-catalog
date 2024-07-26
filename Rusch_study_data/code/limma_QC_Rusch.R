## Date: July 11 2024

# The following code contains a standard limma analysis and some quality control calculations we generated.

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
files <- list.files(input_dir)
files <- files[grep(paste(single_timepoint, collapse = "|"), files)]
length(files) # [1] 7
file.remove(paste(input_dir, files, sep = ""))

# Import the CEL data ----
files <- list.files(input_dir)
files <- files[-grep(paste(pheno_data$geo_accession, collapse = "|"), files)]
length(files) # [1] 62
file.remove(paste(input_dir, files, sep = ""))

data <- ReadAffy(celfile.path = input_dir)

# Raw data QC ----
# affyPLM: generates probe linear model and allows generation of pseudoimages
Pset <- fitPLM(data) # input non-normalized data
pdf(file = paste(output_dir, "RLE_and_NUSE.pdf", sep = ""))
par(mfrow = c(1,2))
RLE(Pset,main = "RLE for Rusch raw data")
NUSE(Pset, main = "NUSE for Rusch raw data")
dev.off()
## The second to last NUSE plot is highly elevated. RLE looks normal.
## Filter out any samples in the meta data that are not in the imported data.
pheno_data <- pheno_data[pheno_data$geo_accession %in% str_extract(colnames(exprs(data)), pattern = "GSM\\d+"),]
pdf(file = paste(output_dir, "GSM2175261_pseudo.pdf", sep = ""))
image(Pset, which = 39)
dev.off()
pheno_data[pheno_data$geo_accession == "GSM2175261",] 
# GSM2175261 Subject 108, time point 2. Anomalies on right of array. This sample will probably stand out in a principal component plot.

# Background correct and normalize the data into an eSet object
data <- affy::rma(data)

# Get expression estimates
normalized_expression <- exprs(data)
normalized_expression[1:5,1:5]
dim(normalized_expression)
# [1] 54675    40

## Raw QC PCAs and pseudoimages ----
colnames(normalized_expression) <- str_extract(colnames(normalized_expression), pattern = "GSM\\d+")
identical(rownames(pheno_data), colnames(normalized_expression)) # [1] TRUE
p <-  pca(normalized_expression, metadata = pheno_data, center = T, scale = T)
pdf(file = paste(output_dir, "pre_qc_pca_plot.pdf", sep = ""))
biplot(p, 
       colby = "patient_ID", 
       legendPosition = "right",
       title = "Rusch et al. data pre-QC")
dev.off()
## Remove patient 108 before further analysis.
## Renormalize the data without this sample.
files <- list.files(input_dir)
files <- files[grep("GSM2175261", files)]
file.remove(paste(input_dir, files, sep = ""))

data <- ReadAffy(celfile.path = input_dir)
data <- affy::rma(data)
normalized_expression <- exprs(data)
dim(normalized_expression)
# [1] 54675    39
colnames(normalized_expression) <- str_extract(colnames(normalized_expression), pattern = "GSM\\d+")
pheno_data <- pheno_data[pheno_data$geo_accession %in% colnames(normalized_expression),]
p <-  pca(normalized_expression, metadata = pheno_data, center = T, scale = T)
pdf(file = paste(output_dir, "pre_qc_pca_plot2.pdf", sep = ""))
biplot(p, 
       colby = "patient_ID", 
       legendPosition = "right",
       title = "Rusch et al. data pre-QC")
dev.off()
pheno_data[pheno_data$geo_accession == "GSM2175214",]
# GSM2175214 Subject 76, timepoint 1 is driving most of the variation in this data.
which(pheno_data$geo_accession == "GSM2175214") # [1] 18
data <- ReadAffy(celfile.path = input_dir)
Pset <- fitPLM(data) # input non-normalized data
pdf(file = paste(output_dir, "RLE_and_NUSE2.pdf", sep = ""))
par(mfrow = c(1,2))
RLE(Pset,main = "RLE for Rusch pre-QC data")
NUSE(Pset, main = "NUSE for Rusch pre-QC data")
dev.off()
## Has a slightly elevated NUSE.
pdf(file = paste(output_dir, "GSM2175214_pseudo.pdf", sep = ""))
image(Pset, which = 18)
dev.off()
## Clear strip pattern on the array, plus the variation in PCA, we removed this sample as well.

files <- list.files(input_dir)
files <- files[grep("GSM2175214", files)]
file.remove(paste(input_dir, files, sep = ""))

## IF YOU HAVE RUN THIS FILE BEFORE, START HERE AFTER RUNNING LIBRARIES AND getGEO. ----
## Change column names to include only GEO sample ID and exclude pheno_data rows that don't include colnames(df).

data <- ReadAffy(celfile.path = input_dir)
data <- affy::rma(data)
normalized_expression <- exprs(data)
dim(normalized_expression)
# [1] 54675    38
colnames(normalized_expression) <- str_extract(colnames(normalized_expression), pattern = "GSM\\d+")
pheno_data <- pheno_data[pheno_data$geo_accession %in% colnames(normalized_expression),]
dim(pheno_data)
# [1] 38  9
p <-  pca(normalized_expression, metadata = pheno_data, center = T, scale = T)
pdf(file = paste(output_dir, "post_qc_pca.pdf", sep = ""))
biplot(p, 
       colby = "patient_ID", 
       legendPosition = "right",
       title = "Rusch et al. data pre-QC")
dev.off()
pheno_data[pheno_data$geo_accession == "GSM2175187",]
# GSM2175187 Subject 10, timepoint 2
which(pheno_data$geo_accession == "GSM2175187") # [1] 6
data <- ReadAffy(celfile.path = input_dir)
Pset <- fitPLM(data) # input non-normalized data
pdf(file = paste(output_dir, "RLE_and_NUSE3.pdf", sep = ""))
par(mfrow = c(1,2))
RLE(Pset,main = "RLE for Rusch post-QC data")
NUSE(Pset, main = "NUSE for Rusch post-QC data")
dev.off()
pdf(file = paste(output_dir, "GSM2175187_pseudo.pdf", sep = ""))
image(Pset, which = 18)
dev.off()

## We removed this patient's baseline sample. What does their follow-up sample look like?
image(Pset, which = which(pheno_data$patient_ID == "Subject 108"))
pheno_data[pheno_data$patient_ID == "Subject 108",]
# title geo_accession                   source age    ethnicity                             race  sex     time  patient_ID
# GSM2175248 md441    GSM2175248 Subject 108, timepoint 1  36 Not Hispanic Native Hawaiian/Pacific Islander Male Baseline Subject 108

## Pseudoimage looks normal, no concerning PCA deviations. This will be our working dataset.

# Pre-processing ----
pheno_data$time[pheno_data$time != "Baseline"] <- rep(x = "Followup", times = length(pheno_data$time[pheno_data$time != "Baseline"]))
table(pheno_data$patient_ID)
# Subject 10 Subject 108 Subject 117  Subject 37  Subject 38  Subject 46  Subject 49  Subject 53  Subject 63  Subject 66  Subject 68 
# 2           1           2           2           2           2           2           2           2           2           2 
# Subject 69  Subject 72  Subject 73  Subject 76  Subject 80   Subject 9  Subject 94  Subject 95  Subject 98 
# 2           2           2           1           2           2           2           2           2 

df <- data.frame(normalized_expression)
identical(colnames(df), pheno_data$geo_accession)
# [1] TRUE
dim(df)
# [1] 54675    38
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
# 1.974886  4.140058  5.610911  6.978722 15.494129

# What are the highly expressed RNAs?
avg_gene_exp <- rowMeans(x = normalized_expression, na.rm = T)
head(avg_gene_exp) # creates vector of each gene's mean
avg_gene_exp <- avg_gene_exp[order(avg_gene_exp, decreasing = T)] # order by decreasing order
keytypes(hgu133plus2.db) # see keys for affymetrix chip
mapIds(x = hgu133plus2.db, keys = head(names(avg_gene_exp)), keytype = "PROBEID", column = "SYMBOL")

# 214414_x_at 211696_x_at 211699_x_at 204018_x_at 217414_x_at 211745_x_at 
# "HBA1"       "HBB"      "HBA1"      "HBA1"      "HBA2"      "HBA1" 

## Control probe plots
source("~/Desktop/work_repo/github/scratch code/Affy_control_probe_check.R")
a <- polyA_check(df, pheno_data, "time")[[1]]
b <- hybrid_check(df, pheno_data, "time")[[1]]
c <- positive_check(df, pheno_data, "time")
d <- polyA_check(df, pheno_data, "title")[[1]]
e <- hybrid_check(df, pheno_data, "title")[[1]]
f <- positive_check(df, pheno_data, "title")
pdf(file = paste(output_dir, "control_probes.pdf", sep = ""), height = 12, width = 12)
plot_grid(a,b,c,d,e,f, nrow = 3)
dev.off()

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

# MDS plot
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colors <- sample(col_vector, length(unique(pheno_data$patient_ID)))
new_colors <- c()
for (i in 1:length(pheno_data$patient_ID)) {
  new_colors[i] <- switch(pheno_data$patient_ID[i],
                          "Subject 10" = colors[1], "Subject 108" = colors[2], "Subject 117" = colors[3], "Subject 37" = colors[4], 
                          "Subject 38" = colors[5], "Subject 46" = colors[6], "Subject 49" = colors[7], "Subject 53" = colors[8],
                          "Subject 63" = colors[9], "Subject 66" = colors[10], "Subject 69" = colors[11], "Subject 72" = colors[12],
                          "Subject 73" = colors[13], "Subject 76" = colors[14], "Subject 80" = colors[15], "Subject 9" = colors[16],
                          "Subject 94" = colors[17], "Subject 95" = colors[18], "Subject 98" = colors[19], "Subject 68" = colors[20]
  )
}
pdf(file = paste(output_dir, "mds.pdf", sep = ""))
limma::plotMDS(x = df, 
               cex = 1, 
               var.explained = T, 
               col = new_colors,
               labels = pheno_data$patient_ID)
title("MDS Plot")
dev.off()

# limma analysis ----

# Remove spike-in probes (AFFX hybridization and polyA probes)
control_index <- grep("lys|Lys|thr|Thr|dap|Dap|phe|Phe|Trp|bioB|BioB|BioC|bioC|bioD|BioD|cre|Cre", x = row.names(df))
rownames(df)[control_index]
df <- df[-control_index,]
dim(df)
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
# [1] 0.4531557

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
# Down       1             0         0         0                  0
# NotSig 54628         54630     54629     54630              54627
# Up         1             0         1         0                  3
dt <- decideTests(x, adjust.method = 'none')
summary(dt)
#         Time ITimeHispanic TimeWhite TimeBlack TimeNativeAmerican
# Down    7116          3545      2118       697               1118
# NotSig 41171         49488     51518     53152              52943
# Up      6343          1597       994       781                569

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
  pdf(file = paste(output_dir, names(y)[i], "p_value_histograms.pdf", sep = ""))
  plot_grid(a, b, labels = c('A. P-value histogram', 'B. BH-adjusted'))
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
                              title = ("Baseline vs. 3mo"),
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
                              max.overlaps = 50)
}
pdf(paste(output_dir, "volcanos.pdf", sep = ""), width = 12, height = 10)
plot_grid(plotlist = volcanos, nrow = 2, align = 'h')
dev.off()