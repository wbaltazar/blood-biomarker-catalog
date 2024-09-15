## Date: July 11 2024
## You should run "pre_QC_Dusek.R" before running this script to remove poor quality arrays from input directory.

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

## INPUT: pre-QC'd .CEL files
input_dir <- "~/Desktop/work_repo/data/GSE10041_RAW/"
## OUTPUT: limma results, plots, figures
output_dir <- "~/Desktop/work_repo/github/Dusek_study_data/output/"

# Import the CEL data ----
data <- ReadAffy(celfile.path = input_dir)

## Raw data QC ----
# affyPLM: generates probe linear model and allows generation of pseudoimages
library(affyPLM)
Pset <- fitPLM(data) #input non-normalized data
pdf(file = paste(output_dir, "RLE_and_NUSE.pdf", sep = ""))
par(mfrow = c(1,2))
RLE(Pset,main="RLE for Dusek et al.")
NUSE(Pset,main="NUSE for Dusek et al.")
dev.off()

## Normalize the data into an eSet ----
data <- affy::rma(data)

# Get expression estimates
normalized_expression <- exprs(data)
dim(normalized_expression)
# [1] 54675    44

# Generate informative plots / data
pdf(file = paste(output_dir, "features_hist.pdf", sep = ""))
hist(normalized_expression)
dev.off()
pdf(file = paste(output_dir, "samples_box.pdf", sep = ""))
boxplot(normalized_expression)
dev.off()
quantile(normalized_expression)
# 0%       25%       50%       75%      100% 
# 2.416043  4.457443  5.600040  6.714875 14.765905 

# What are the highly expressed RNAs?
avg_gene_exp <- rowMeans(x = normalized_expression, na.rm = TRUE)
head(avg_gene_exp) # creates vector of each gene's mean
avg_gene_exp <- avg_gene_exp[order(avg_gene_exp, decreasing = TRUE)] # order by decreasing order
keytypes(hgu133plus2.db) # see keys for affymetrix chip
mapIds(x = hgu133plus2.db, keys = head(names(avg_gene_exp)), keytype = "PROBEID", column = "SYMBOL")

# 214414_x_at 211745_x_at 211696_x_at 209458_x_at 204018_x_at 211699_x_at 
# "HBA1"      "HBA1"       "HBB"      "HBA1"      "HBA1"      "HBA1" 

# Load in metadata ----
gset <- getGEO("GSE10041", GSEMatrix =TRUE, AnnotGPL=TRUE)
pheno_data_org <- pData(gset$GSE10041_series_matrix.txt.gz)
feature_data <- fData(gset$GSE10041_series_matrix.txt.gz)
dim(pheno_data_org)
# [1] 72 33
# filter out traits we want
pheno_data <- pheno_data_org %>% 
  dplyr::select(title = title,
                geo_accession = geo_accession)
pheno_data$time <- str_extract(pheno_data_org$title, "N.")
pheno_data$id <- str_extract(pheno_data_org$title, "-.*") %>% 
  str_replace("-","")
dim(pheno_data)
# [1] 72 4
length(grep("M", pheno_data$title)) # [1] 25 This is the correct number in the original and validation groups.
pheno_data <- pheno_data[-grep("M", pheno_data$title),]
dim(pheno_data)
# [1] 47  4

table(pheno_data$time, pheno_data$id)
#    01 02 03 06 07 11 12 14 23 30 32 33 34 35 36 39 40 41 46 47 49 50 51 59
# N1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  0  1  1  1  1  1  1  1  1  1
# N2  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1
sum(na.omit(pheno_data$time) == "N1") #[1] 23
sum(na.omit(pheno_data$time) == "N2") #[1] 24
# In this study, the researchers collected 19 N1 data points and 20 N2 datapoints for their analysis (one
# of the N1 arrays had too low of an intensity.) They also had a validation cohort of 5 additional individuals.
# The pheno data indicates we are one sample short. This is explained in GEO: there were 4 validation samples.
# However, there isn't a way to discern which cohort is which using the pheno_data.

# rename expression data columns
colnames(normalized_expression)
colnames(normalized_expression) <- gsub(pattern = "(.CEL.gz)", replacement = "", x = colnames(normalized_expression))
dim(normalized_expression)
df <- data.frame(normalized_expression)

pheno_data <- pheno_data[pheno_data$geo_accession %in% colnames(df),]
dim(pheno_data)
# [1] 44  4

# remove unnecessary columns
df <- df[,colnames(df) %in% pheno_data$geo_accession]
normalized_expression <- normalized_expression[,colnames(normalized_expression) %in% pheno_data$geo_accession]
dim(df)
# [1] 54675    44

# make sure that the order of samples in the pheno data matches the order in exp data
index <- match(pheno_data$geo_accession, colnames(df))
pheno_data <- pheno_data[index,]
dim(pheno_data)
# [1] 44  4
identical(pheno_data$geo_accession, colnames(df))
# [1] TRUE

# Annotate sex
rps4y_vals <- normalized_expression["201909_at",]
plot(rps4y_vals)
pheno_data$sex <- NA
for (i in 1:length(rps4y_vals)) {
  if (rps4y_vals[i] > 8) {
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

# MDS plot
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colors <- sample(col_vector, 23)
new_colors <- c()
for (i in 1:length(pheno_data$id)) {
  new_colors[i] <- switch(pheno_data$id[i],
                          "01" = colors[1], "02" = colors[2], "03" = colors[3], "06" = colors[4], "07" = colors[5], "59" = colors[6],
                          "11" = colors[7], "12" = colors[8], "14" = colors[9], "23" = colors[10], "32" = colors[11], "33" = colors[12],
                          "34" = colors[13], "35" = colors[14], "36" = colors[15], "39" = colors[16], "40" = colors[17], "41" = colors[18],
                          "46" = colors[19], "47" = colors[20], "49" = colors[21], "50" = colors[22], "51" = colors[23],
  )
}
pdf(file = paste(output_dir, "mds.pdf", sep = ""))
limma::plotMDS(x = df, 
               cex = 1, 
               var.explained = T, 
               col = new_colors,
               labels = paste("P", pheno_data$id, sep = "-"))
title("MDS Plot")
dev.off()


# limma analysis ----

# Remove spike-in probes (AFFX hybridization and polyA probes)
control_index <- grep("lys|Lys|thr|Thr|dap|Dap|phe|Phe|Trp|bioB|BioB|BioC|bioC|bioD|BioD|cre|Cre", x = row.names(df))
rownames(df)[control_index]
df <- df[-control_index,]
dim(df)
# [1] 54630    42
rownames(df)[grep("AFFX", x = rownames(df))]

id <- factor(pheno_data$id)
time <- factor(pheno_data$time)
sex <- factor(pheno_data$sex)

TS <- factor(paste(time,sex, sep = "."))
design <- model.matrix(~0+TS, df)
colnames(design) <- make.names(levels(TS))
cor <- duplicateCorrelation(df, design, block=id)
cor$consensus.correlation
# [1] 0.2657351

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
# NotSig 54630    54630
# Up         0        0
dt <- decideTests(x, adjust.method = "none")
summary(dt)
#         Time Time_sex
# Down     874     3399
# NotSig 52928    49655
# Up       828     1576
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