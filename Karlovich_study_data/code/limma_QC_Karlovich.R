# Date: October 25 2024

# Karlovich, Chris, Guillemette Duchateau-Nguyen, Andrea Johnson, Patricia McLoughlin, Mercidita Navarro, Carole 
# Fleurbaey, Lori Steiner, et al. “A Longitudinal Study of Gene Expression in Healthy Individuals.” BMC Medical 
# Genomics 2 (June 7, 2009): 33. https://doi.org/10.1186/1755-8794-2-33.


## Note that the first timepoints of this study were hybridized in one batch, while the last two timepoints were
## hybridized in another batch. We do not compare these groups of times. See paper.

# Load libraries -----

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

## INPUT: GSE16028 raw data from GEO
input_dir <- "~/Desktop/work_repo/data/"
## OUTPUT: limma results, QC, and graphics
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

# Read in raw data for batch 1 ----
files <- list.files(paste(input_dir, "GSE16028_RAW", sep = ""), full.names = T)
files <- files[str_extract(files, "GSM\\d+") %in% p1$geo_accession]
length(files)
# [1] 66
raw <- ReadAffy(filenames = files)

## Raw data QC ----
# affyPLM: generates probe linear model and allows generation of pseudoimages
library(affyPLM)
Pset <- fitPLM(raw, background = F, normalize = F) #input non-normalized data
pdf(file = paste(output_dir, "RLE_and_NUSE_batch1.pdf", sep = ""))
par(mfrow = c(1,2))
RLE(Pset,main="RLE for Karlovich et al. batch 1", ylim = c(-1.5,1.5))
NUSE(Pset,main="NUSE for Karlovich et al. batch 1", ylim = c(0.9, 1.3))
dev.off()
# 11 and 29 are the most elevated NUSE values 48 had the highest RLE.
colnames(raw)[c(11,29,48)]
p1[c("GSM401075","GSM401105","GSM401136"),]
#                                             title geo_accession subject age    sex   time
# GSM401075 Blood sample from patient 114 at Day 14     GSM401075     114  39   Male Day 14
# GSM401105 Blood sample from patient 138 at Day 14     GSM401105     138  62   Male Day 14
# GSM401136 Blood sample from patient 163 at Day 28     GSM401136     163  64 Female Day 28
pdf(file = paste(output_dir, "pseudo_GSM401075.pdf", sep = ""))
image(Pset, which = 11)
dev.off()
pdf(file = paste(output_dir, "pseudo_GSM401105.pdf", sep = ""))
image(Pset, which = 29)
dev.off() # Strong visual artifact
pdf(file = paste(output_dir, "pseudo_GSM401136.pdf", sep = ""))
image(Pset, which = 48) # Strong visual artifact
dev.off()

## Get expression calls using MAS 5.0 ----
calls <- mas5calls.AffyBatch(raw)
probe_pval <- assayData(calls)[["se.exprs"]] ## Returns p-values
dim(probe_pval)
# [1] 54675    66
minSamples <- min(colSums(table(p1$subject, p1$time)))
print(minSamples)
# [1] 22
expressed <- rowSums(probe_pval < 0.05) >= minSamples

## Normalize the data into an eSet ----
data <- affy::rma(raw)
unfiltered_normalized_expression <- exprs(data)
dim(unfiltered_normalized_expression)
# [1] 54675    66
normalized_expression <- unfiltered_normalized_expression[expressed,]
dim(normalized_expression)
# [1] 17635    66

# Generate informative plots / data
pdf(file = paste(output_dir, "features_hist_batch1.pdf", sep = ""))
hist(normalized_expression)
dev.off()
pdf(file = paste(output_dir, "samples_box_batch1.pdf", sep = ""))
boxplot(normalized_expression)
dev.off()
quantile(normalized_expression)
# 0%       25%       50%       75%      100% 
# 1.523538  4.410875  5.636163  6.831971 14.276604 

# What are the highly expressed RNAs?
avg_gene_exp <- rowMeans(x = normalized_expression, na.rm = TRUE)
head(avg_gene_exp) # creates vector of each gene's mean
avg_gene_exp <- avg_gene_exp[order(avg_gene_exp, decreasing = TRUE)] # order by decreasing order
keytypes(hgu133plus2.db) # see keys for affymetrix chip
mapIds(x = hgu133plus2.db, keys = head(names(avg_gene_exp)), keytype = "PROBEID", column = "SYMBOL")

# 211696_x_at 214414_x_at 209116_x_at 217232_x_at 211745_x_at 204018_x_at 
# "HBB"      "HBA1"       "HBB"       "HBB"      "HBA1"      "HBA1" 

# rename expression data columns
colnames(normalized_expression)
colnames(normalized_expression) <- gsub(pattern = "(.CEL.gz)", replacement = "", x = colnames(normalized_expression)) # greps file extension and replaces it
colnames(unfiltered_normalized_expression) <- gsub(pattern = "(.CEL.gz)", replacement = "", x = colnames(normalized_expression))
dim(normalized_expression)
df <- data.frame(normalized_expression)
identical(colnames(df), p1$geo_accession)
# [1] TRUE

# Principal Component Analysis ----
p <- pca(df, metadata = p1, center = T, scale = T)
p_un <- pca(unfiltered_normalized_expression, metadata = p1, center = T, scale = T)

scree <- screeplot(p, components = getComponents(p, 1:20))

## Compare filtered to unfiltered. Do the samples with poor PLM statistics stand out?
pdf(file = paste(output_dir, "unfiltered_pca_batch1.pdf", sep = ""))
biplot(p_un, 
       showLoadings = F,
       lab = NULL, 
       colby = "subject", shape = "time",
       legendPosition = 'right', legendLabSize = 12, legendIconSize = 3,
       title = 'Karlovich batch 1',
       subtitle = 'lowly-expressed genes unfiltered')
dev.off()
pdf(file = paste(output_dir, "filtered_pca_batch1.pdf", sep = ""))
biplot(p, 
       showLoadings = F,
       lab = NULL, 
       colby = "subject", shape = "time",
       legendPosition = 'right', legendLabSize = 9, legendIconSize = 3,
       title = 'Karlovich batch 1',
       subtitle = 'lowly-expressed genes filtered')
dev.off()

## Principal component analysis on filtered, normalized expression data
ptime <- biplot(p, 
                showLoadings = F,
                lab = NULL, 
                colby = "time", 
                legendPosition = 'right', legendLabSize = 12, legendIconSize = 3,
                title = 'Karlovich batch 1',
                subtitle = 'lowly-expressed genes filtered')
psex <- biplot(p, 
               showLoadings = F,
               lab = NULL, 
               colby = "sex", colkey = c('turquoise', 'orange'),
               legendPosition = 'right', legendLabSize = 12, legendIconSize = 3)
psub <- biplot(p, 
               showLoadings = F,
               lab = NULL, 
               colby = "subject",
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
pdf(file = paste(output_dir, "pca_batch1.pdf", sep = ""), height = 12, width = 16)
plot_grid(toprow, botrow, ncol = 1, align = 'v')
dev.off()

# MDS plot
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colors <- sample(col_vector, 22)
new_colors <- c()
for (i in 1:length(p1$subject)) {
  new_colors[i] <- switch(p1$subject[i],
                          "106" = colors[1], "108" = colors[2], "113" = colors[3], "114" = colors[4], "117" = colors[5], 
                          "121" = colors[6], "123" = colors[7], "131" = colors[8], "135" = colors[9], "138" = colors[10], 
                          "147" = colors[11], "148" = colors[12], "173" = colors[13], "174" = colors[14], "179" = colors[15],
                          "180" = colors[16], "154" = colors[17], "159" = colors[18], "160" = colors[19], "163" = colors[20],
                          "165" = colors[21], "172" = colors[22],
  )
}
pdf(file = paste(output_dir, "mds_batch1.pdf", sep = ""))
limma::plotMDS(x = df, 
               cex = 1, 
               var.explained = T, 
               col = new_colors,
               labels = p1$subject)
title("MDS Plot (Batch 1)")
dev.off()

# limma analysis ----


# Remove spike-in probes (AFFX hybridization and polyA probes)
control_index <- grep("lys|Lys|thr|Thr|dap|Dap|phe|Phe|Trp|bioB|BioB|BioC|bioC|bioD|BioD|cre|Cre", x = row.names(df))
rownames(df)[control_index]
# already removed
rownames(df)[grep("AFFX", x = rownames(df))]
dim(df)

# [1] "AFFX-HSAC07/X00351_3_at"     "AFFX-HSAC07/X00351_5_at"     "AFFX-HSAC07/X00351_M_at"    
# [4] "AFFX-hum_alu_at"             "AFFX-HUMGAPDH/M33197_3_at"   "AFFX-HUMGAPDH/M33197_5_at"  
# [7] "AFFX-HUMGAPDH/M33197_M_at"   "AFFX-HUMISGF3A/M97935_3_at"  "AFFX-HUMISGF3A/M97935_5_at" 
# [10] "AFFX-HUMISGF3A/M97935_MA_at" "AFFX-HUMISGF3A/M97935_MB_at" "AFFX-HUMRGE/M10098_3_at"    
# [13] "AFFX-HUMRGE/M10098_5_at"     "AFFX-HUMRGE/M10098_M_at"     "AFFX-M27830_3_at"           
# [16] "AFFX-M27830_5_at"            "AFFX-M27830_M_at"     

# combine all factors into unique identifiers, using ID as indices
sex <- factor(p1$sex)
time <- factor(p1$time)
id <- factor(p1$subject)
TS <- paste(sex, time, sep=".")
TS <- factor(TS)

design <- model.matrix(~0+TS, df) # 0 indicates y-intercept

colnames(design) <- make.names(levels(TS))
cor <- duplicateCorrelation(df, design, block=id) # here, ID identifies which samples correlate
cor$consensus.correlation
# [1] 0.168133

fit <- lmFit(object=df, design=design, block=id, correlation=cor$consensus.correlation)

#To find the genes changing across time regardless of sex
my_contrasts <- makeContrasts(Day14 = (Female.Day.14 + Male.Day.14)/2 - (Female.Day.1 + Male.Day.1)/2,
                                   Day28 = (Female.Day.28 + Male.Day.28)/2 - (Female.Day.1 + Male.Day.1)/2,
                                   IMaleDay14 = (Male.Day.14 - Male.Day.1)/2 - (Female.Day.14 - Female.Day.1)/2,
                                   IMaleDay28 = (Male.Day.14 - Male.Day.1)/2 - (Female.Day.14 - Female.Day.1)/2,
                                   levels = colnames(design))
                                   

fit2 <- contrasts.fit(fit, my_contrasts)
x <- eBayes(fit2)

# Download results and visualizations ----
dir.create(paste(output_dir, "batch1_limma", sep = ""))
batch_dir <- paste(output_dir, "batch1_limma/", sep = "")
y <- topTable(x, number = Inf, adjust.method = "fdr", sort.by = "F")
y$Symbol <- mapIds(x = hgu133plus2.db, keys = row.names(y), column = "SYMBOL", keytype = "PROBEID")
write.csv(y, file = paste(batch_dir, "limma_F.csv", sep = ""), row.names = T)

dt <- decideTests(x)
summary(dt)
#         Time Time_sex
# Down       0        0
# NotSig 54630    54630
# Up         0        0
dt <- decideTests(x, adjust.method = "none")
summary(dt)
#         Time Time_sex
# Down     874     3399
# NotSig 52928    49655
# Up       828     1576
names(y)
# [1] "Day14"      "Day28"      "IMaleDay14" "IMaleDay28" "AveExpr"    "F"          "P.Value"    "adj.P.Val"  "Symbol"    
tables <- list()
for (i in 1:4) {
  tables[[i]] <- topTable(x, number = Inf, coef = i, adjust.method = "fdr")
  names(tables)[i] <- names(y)[i]
  tables[[i]]$Symbol <- mapIds(x = hgu133plus2.db, keys = row.names(tables[[i]]), column = "SYMBOL", keytype = "PROBEID")
  write.csv(tables[[i]], file = paste(batch_dir, "limma_", names(y)[i], ".csv", sep = ""), row.names = T)
}

## Generate p-value histograms----
histlist <- list()
for (i in 1:4) {
  histlist[[i]] <- ggplot(data = tables[[i]], aes(x = P.Value)) +
    geom_histogram(fill = "lightblue", color = "black") +
    ggtitle(names(tables)[i]) +
    theme_cowplot() +
    theme(title = element_text(size = 7.5))
}
pdf(file = paste(batch_dir, "p_value_histograms.pdf", sep = ""))
plot_grid(plotlist = histlist, nrow = 2)
dev.off()

histlist <- list()
for (i in 1:4) {
  histlist[[i]] <- ggplot(data = tables[[i]], aes(x = adj.P.Val)) +
    geom_histogram(fill = "lightblue", color = "black") +
    ggtitle(names(tables)[i]) +
    theme_cowplot() +
    theme(title = element_text(size = 7.5))
}
pdf(file = paste(batch_dir, "adj_p_value_histograms.pdf", sep = ""))
plot_grid(plotlist = histlist, nrow = 2)
dev.off()

## EnhancedVolcano plots----
volcanos <- list()
for (i in 1:4) {
  volcanos[[i]] <- EnhancedVolcano(tables[[i]],
                                   lab = tables[[i]]$Symbol,
                                   x = "logFC",
                                   y = "P.Value",
                                   title = "", titleLabSize = 8,
                                   subtitle = names(tables)[i], subtitleLabSize = 12,
                                   caption = "",
                                   col = c("black","lightblue","royalblue","blue"),
                                   legendPosition = 'none',
                                   pCutoff = 1e-2,
                                   FCcutoff = 1,
                                   labSize = (4.0),
                                   colAlpha = 1,
                                   pointSize = 1.0,
                                   drawConnectors = TRUE,
                                   max.overlaps = 50)
}
pdf(paste(batch_dir, "volcanos.pdf", sep = ""))
plot_grid(plotlist = volcanos, nrow = 2)
dev.off()

# Read in raw data for batch 2 ----
files <- list.files(paste(input_dir, "GSE16028_RAW", sep = ""), full.names = T)
files <- files[str_extract(files, "GSM\\d+") %in% p2$geo_accession]
length(files)
# [1] 43
data <- ReadAffy(filenames = files)

## Raw data QC ----
# affyPLM: generates probe linear model and allows generation of pseudoimages
library(affyPLM)
Pset <- fitPLM(data, background = F, normalize = F) #input non-normalized data
pdf(file = paste(output_dir, "RLE_and_NUSE_batch2.pdf", sep = ""))
par(mfrow = c(1,2))
RLE(Pset,main="RLE for Karlovich et al. batch 2")
NUSE(Pset,main="NUSE for Karlovich et al. batch 2")
dev.off()
# Most elevated RLE was 5
# Elevated NUSE: 5, 11, 12, 20, 21, 36, 39
inspect <- colnames(data)[c(5, 11, 12, 20, 21, 36, 39)]
p2[c(str_replace(inspect, "\\.CEL\\.gz", "")),]
#                                              title geo_accession subject age    sex    time
# GSM401072  Blood sample from patient 113 at Day 90     GSM401072     113  28   Male  Day 90
# GSM401087  Blood sample from patient 121 at Day 90     GSM401087     121  62   Male  Day 90
# GSM401088 Blood sample from patient 121 at Day 180     GSM401088     121  62   Male Day 180
# GSM401108 Blood sample from patient 138 at Day 180     GSM401108     138  62   Male Day 180
# GSM401112  Blood sample from patient 147 at Day 90     GSM401112     147  31 Female  Day 90
# GSM401148 Blood sample from patient 172 at Day 180     GSM401148     172  60 Female Day 180
# GSM401157  Blood sample from patient 174 at Day 90     GSM401157     174  56 Female  Day 90
index <- c(5, 11, 12, 20, 21, 36, 39)
for (i in index) {
  pdf(file = paste(output_dir, str_replace(inspect, "\\.CEL\\.gz", "_pseudo.pdf")[which(index == i)], sep = ""))
  image(Pset, which = i)
  dev.off()
}

## Normalize the data ----
raw <- data
data <- affy::rma(data)
unfiltered_normalized_expression <- exprs(data)
dim(unfiltered_normalized_expression)
# [1] 54675    43

## Get expression calls using MAS 5.0 ----
calls <- mas5calls.AffyBatch(raw)
probe_pval <- assayData(calls)[["se.exprs"]] ## Returns p-values
dim(probe_pval)
# [1] 54675    43
minSamples <- min(colSums(table(p2$subject, p2$time)))
print(minSamples)
# [1] 21
expressed <- rowSums(probe_pval < 0.05) >= minSamples
normalized_expression <- unfiltered_normalized_expression[expressed,]
dim(normalized_expression)
# [1] 16384    43

colnames(normalized_expression) <- gsub(pattern = "(.CEL.gz)", replacement = "", x = colnames(normalized_expression))
colnames(unfiltered_normalized_expression) <- gsub(pattern = "(.CEL.gz)", replacement = "", x = colnames(unfiltered_normalized_expression))
identical(colnames(normalized_expression), p2$geo_accession)
# [1] TRUE
p <- pca(normalized_expression, metadata = p2, center = T, scale = T)
p_unf <- pca(unfiltered_normalized_expression, metadata = p2, center = T, scale = T)
pdf(file = paste(output_dir, "filtered_pca_batch2.pdf", sep = ""))
biplot(p, 
       showLoadings = F,
       lab = NULL, 
       colby = "subject", shape = "time",
       legendPosition = 'right', legendLabSize = 12, legendIconSize = 3,
       title = 'Karlovich batch 2',
       subtitle = 'lowly-expressed genes filtered')
dev.off()
pdf(file = paste(output_dir, "unfiltered_pca_batch2.pdf", sep = ""))
biplot(p_unf, 
       showLoadings = F,
       lab = NULL, 
       colby = "subject", shape = "time",
       legendPosition = 'right', legendLabSize = 12, legendIconSize = 3,
       title = 'Karlovich batch 2',
       subtitle = 'lowly-expressed genes not filtered')
dev.off()

## Samples to remove: GSM401087, GSM401148, GSM401157.
p2 <- p2[-grep("GSM401087|GSM401148|GSM401157", p2$geo_accession),]
files <- list.files(paste(input_dir, "GSE16028_RAW", sep = ""), full.names = T)
files <- files[str_extract(files, "GSM\\d+") %in% p2$geo_accession]
length(files)
# [1] 40
data <- ReadAffy(filenames = files)
calls <- mas5calls.AffyBatch(data)
probe_pval <- assayData(calls)[["se.exprs"]] ## Returns p-values
dim(probe_pval)
# [1] 54675    40
minSamples <- min(colSums(table(p2$subject, p2$time)))
print(minSamples)
# [1] 20
expressed <- rowSums(probe_pval < 0.05) >= minSamples
data <- affy::rma(data)
unfiltered_normalized_expression <- exprs(data)
normalized_expression <- unfiltered_normalized_expression[expressed,]
dim(normalized_expression)
# [1] 16352    40

# Generate informative plots / data
pdf(file = paste(output_dir, "features_hist_batch2.pdf", sep = ""))
hist(normalized_expression)
dev.off()
pdf(file = paste(output_dir, "samples_box_batch2.pdf", sep = ""))
boxplot(normalized_expression)
dev.off()
quantile(normalized_expression)
# 0%       25%       50%       75%      100% 
# 2.193795  5.202867  6.430439  7.594720 14.671843 

# What are the highly expressed RNAs?
avg_gene_exp <- rowMeans(x = normalized_expression, na.rm = TRUE)
head(avg_gene_exp) # creates vector of each gene's mean
avg_gene_exp <- avg_gene_exp[order(avg_gene_exp, decreasing = TRUE)] # order by decreasing order
keytypes(hgu133plus2.db) # see keys for affymetrix chip
mapIds(x = hgu133plus2.db, keys = head(names(avg_gene_exp)), keytype = "PROBEID", column = "SYMBOL")

# 211696_x_at 214414_x_at 209116_x_at 217232_x_at 211745_x_at 204018_x_at 
# "HBB"      "HBA1"       "HBB"       "HBB"      "HBA1"      "HBA1" 

# rename expression data columns
colnames(normalized_expression)
colnames(normalized_expression) <- gsub(pattern = "(.CEL.gz)", replacement = "", x = colnames(normalized_expression)) # greps file extension and replaces it
dim(normalized_expression)
df <- data.frame(normalized_expression)
identical(colnames(df), p2$geo_accession)
# [1] TRUE

# Principal Component Analysis on QC'd, filtered data ----
p <- pca(df, metadata = p2, center = T, scale = T)
scree <- screeplot(p, components = getComponents(p, 1:20))
ptime <- biplot(p, 
                showLoadings = F,
                lab = NULL, 
                colby = "time", 
                legendPosition = 'right', legendLabSize = 12, legendIconSize = 3,
                title = 'Karlovich batch 2',
                subtitle = 'lowly-expressed genes filtered')
psex <- biplot(p, 
               showLoadings = F,
               lab = NULL, 
               colby = "sex", colkey = c('turquoise', 'orange'),
               legendPosition = 'right', legendLabSize = 12, legendIconSize = 3)
psub <- biplot(p, 
               showLoadings = F,
               lab = NULL, 
               colby = "subject",
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
pdf(file = paste(output_dir, "pca_batch2.pdf", sep = ""), height = 12, width = 16)
plot_grid(toprow, botrow, ncol = 1, align = 'v')
dev.off()
# Day 90 subject 121 (GSM401087) and Day 180 subject 172 (GSM401148) both had elevated NUSEs and
# drive PC1 and PC2. Removing both of these samples later.

# MDS plot
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colors <- sample(col_vector, 22)
new_colors <- c()
for (i in 1:length(p2$subject)) {
  new_colors[i] <- switch(p2$subject[i],
                          "106" = colors[1], "108" = colors[2], "113" = colors[3], "114" = colors[4], "117" = colors[5], 
                          "121" = colors[6], "123" = colors[7], "131" = colors[8], "135" = colors[9], "138" = colors[10], 
                          "147" = colors[11], "148" = colors[12], "173" = colors[13], "174" = colors[14], "179" = colors[15],
                          "180" = colors[16], "154" = colors[17], "159" = colors[18], "160" = colors[19], "163" = colors[20],
                          "165" = colors[21], "172" = colors[22],
  )
}
pdf(file = paste(output_dir, "mds_batch2.pdf", sep = ""))
limma::plotMDS(x = df, 
               cex = 1, 
               var.explained = T, 
               col = new_colors,
               labels = paste(p2$subject, p2$time))
title("MDS Plot (batch 2)")
dev.off()

# limma analysis ----
# Remove spike-in probes (AFFX hybridization and polyA probes)
control_index <- grep("lys|Lys|thr|Thr|dap|Dap|phe|Phe|Trp|bioB|BioB|BioC|bioC|bioD|BioD|cre|Cre", x = row.names(df))
rownames(df)[control_index]
dim(df)
rownames(df)[grep("AFFX", x = rownames(df))]

# [1] "AFFX-HSAC07/X00351_3_at"     "AFFX-HSAC07/X00351_5_at"     "AFFX-HSAC07/X00351_M_at"    
# [4] "AFFX-hum_alu_at"             "AFFX-HUMGAPDH/M33197_3_at"   "AFFX-HUMGAPDH/M33197_5_at"  
# [7] "AFFX-HUMGAPDH/M33197_M_at"   "AFFX-HUMISGF3A/M97935_3_at"  "AFFX-HUMISGF3A/M97935_5_at" 
# [10] "AFFX-HUMISGF3A/M97935_MA_at" "AFFX-HUMISGF3A/M97935_MB_at" "AFFX-HUMRGE/M10098_3_at"    
# [13] "AFFX-HUMRGE/M10098_5_at"     "AFFX-HUMRGE/M10098_M_at"     "AFFX-M27830_3_at"           
# [16] "AFFX-M27830_5_at"            "AFFX-M27830_M_at"     

# combine all factors into unique identifiers, using ID as indices
sex <- factor(p2$sex)
time <- factor(p2$time)
id <- factor(p2$subject)
TS <- paste(sex, time, sep=".")
TS <- factor(TS)

design <- model.matrix(~0+TS, df) # 0 indicates y-intercept

colnames(design) <- make.names(levels(TS))
cor <- duplicateCorrelation(df, design, block=id) # here, ID identifies which samples correlate
cor$consensus.correlation
# [1] 0.2329959

fit <- lmFit(object=df, design=design, block=id, correlation=cor$consensus.correlation)

#To find the genes changing across time regardless of sex
my_contrasts <- makeContrasts(Time = (Female.Day.180 + Male.Day.180)/2 - (Female.Day.90 + Male.Day.90)/2,
                              IMaleTime = (Male.Day.180 - Male.Day.90)/2 - (Female.Day.180 - Female.Day.90)/2,
                              levels = colnames(design))


fit2 <- contrasts.fit(fit, my_contrasts)
x <- eBayes(fit2)

# Download results and visualizations ----
dir.create(paste(output_dir, "batch2_limma", sep = ""))
batch_dir <- paste(output_dir, "batch2_limma/", sep = "")
y <- topTable(x, number = Inf, adjust.method = "fdr", sort.by = "F")
y$Symbol <- mapIds(x = hgu133plus2.db, keys = row.names(y), column = "SYMBOL", keytype = "PROBEID")
write.csv(y, file = paste(batch_dir, "limma_F.csv", sep = ""), row.names = T)

dt <- decideTests(x)
summary(dt)
#         Time IMaleTime
# Down     459         0
# NotSig 15858     16352
# Up        35         0
dt <- decideTests(x, adjust.method = "none")
summary(dt)
#         Time IMaleTime
# Down    1957       335
# NotSig 13526     15913
# Up       869       104
names(y)
# [1] "Time"      "IMaleTime" "AveExpr"   "F"         "P.Value"   "adj.P.Val" "Symbol"
tables <- list()
for (i in 1:2) {
  tables[[i]] <- topTable(x, number = Inf, coef = i, adjust.method = "fdr")
  names(tables)[i] <- names(y)[i]
  tables[[i]]$Symbol <- mapIds(x = hgu133plus2.db, keys = row.names(tables[[i]]), column = "SYMBOL", keytype = "PROBEID")
  write.csv(tables[[i]], file = paste(batch_dir, "limma_", names(y)[i], ".csv", sep = ""), row.names = T)
}

## Generate p-value histograms----
histlist <- list()
for (i in 1:2) {
  histlist[[i]] <- ggplot(data = tables[[i]], aes(x = P.Value)) +
    geom_histogram(fill = "lightblue", color = "black") +
    ggtitle(names(tables)[i]) +
    theme_cowplot() +
    theme(title = element_text(size = 7.5))
}
pdf(file = paste(batch_dir, "p_value_histograms.pdf", sep = ""))
plot_grid(plotlist = histlist, nrow = 1)
dev.off()

histlist <- list()
for (i in 1:2) {
  histlist[[i]] <- ggplot(data = tables[[i]], aes(x = adj.P.Val)) +
    geom_histogram(fill = "lightblue", color = "black") +
    ggtitle(names(tables)[i]) +
    theme_cowplot() +
    theme(title = element_text(size = 7.5))
}
pdf(file = paste(batch_dir, "adj_p_value_histograms.pdf", sep = ""))
plot_grid(plotlist = histlist, nrow = 1)
dev.off()

## EnhancedVolcano plots----
volcanos <- list()
for (i in 1:2) {
  volcanos[[i]] <- EnhancedVolcano(tables[[i]],
                                   lab = tables[[i]]$Symbol,
                                   x = "logFC",
                                   y = "P.Value",
                                   title = "", titleLabSize = 8,
                                   subtitle = names(tables)[i], subtitleLabSize = 12,
                                   caption = "",
                                   col = c("black","lightblue","royalblue","blue"),
                                   legendPosition = 'none',
                                   pCutoff = 1e-2,
                                   FCcutoff = 1,
                                   labSize = (4.0),
                                   colAlpha = 1,
                                   pointSize = 1.0,
                                   drawConnectors = TRUE,
                                   max.overlaps = 50)
}
pdf(paste(batch_dir, "volcanos.pdf", sep = ""))
plot_grid(plotlist = volcanos, nrow = 1)
dev.off()