## Date: Jan 17 2025

## This script performs limma differential expression analysis of data from the following study:
# Gómez-Carballa, Alberto, Laura Navarro, Jacobo Pardo-Seco, Xabier Bello, Sara Pischedda, Sandra Viz-Lasheras, 
# Alba Camino-Mera, et al. “Music Compensates for Altered Gene Expression in Age-Related Cognitive Disorders.” 
# Scientific Reports 13, no. 1 (December 2, 2023): 21259. https://doi.org/10.1038/s41598-023-48094-5.

## Subjects had blood drawn into PAXgene Tubes before and after a 50 minute concert. Data on
## healthy controls was collected.

## Load libraries ----
library(edgeR)
library(limma)
library(AnnotationDbi)
library(Homo.sapiens)
library(cowplot)
library(PCAtools) 
library(GEOquery)
library(tidyverse)
library(ggalt)

## INPUT: Raw count file from GSE239282 on Gene Expression Omnibus ----
input_dir <- "~/Desktop/work_repo/data/GSE239282_raw_counts.txt.gz"
## OUTPUT: Limma top table, quality control plots, volcano plots, and p-value histograms ----
output_dir <- "~/Desktop/work_repo/github/Gomez_study_data/output/"

# Step 1: pull sample information ----
gset <- getGEO(GEO = "GSE239282")
p <- phenoData(gset$GSE239282_series_matrix.txt.gz)
pheno_data_raw <- p@data
dim(pheno_data_raw)
# [1] 60 55
colnames(pheno_data_raw)
# [1] "title"                   "geo_accession"           "status"                  "submission_date"         "last_update_date"        "type"                   
# [7] "channel_count"           "source_name_ch1"         "organism_ch1"            "characteristics_ch1"     "characteristics_ch1.1"   "characteristics_ch1.2"  
# [13] "characteristics_ch1.3"   "characteristics_ch1.4"   "characteristics_ch1.5"   "molecule_ch1"            "extract_protocol_ch1"    "extract_protocol_ch1.1" 
# [19] "taxid_ch1"               "description"             "description.1"           "description.2"           "data_processing"         "data_processing.1"      
# [25] "data_processing.2"       "data_processing.3"       "data_processing.4"       "data_processing.5"       "data_processing.6"       "data_processing.7"      
# [31] "data_processing.8"       "platform_id"             "contact_name"            "contact_email"           "contact_department"      "contact_institute"      
# [37] "contact_address"         "contact_city"            "contact_state"           "contact_zip/postal_code" "contact_country"         "data_row_count"         
# [43] "instrument_model"        "library_selection"       "library_source"          "library_strategy"        "relation"                "relation.1"             
# [49] "supplementary_file_1"    "age:ch1"                 "gender:ch1"              "group:ch1"               "subject id:ch1"          "timepoint:ch1"          
# [55] "tissue:ch1"   

table(pheno_data_raw$`timepoint:ch1`)
# 1  2 
# 6 54 This column is mislabeled, and contains faulty sample labels.
length(grep("Timepoint-2", pheno_data_raw$title[which(pheno_data_raw$`group:ch1` == "Control")]))
# [1] 14 Correct number of healthy controls, as stated by study authors.
pheno_data <- pheno_data_raw %>% 
  dplyr::select(title = title,
         geo_accession = geo_accession,
         age = `age:ch1`,
         condition = `group:ch1`,
         participant = `subject id:ch1`,
         sex = `gender:ch1`)
dim(pheno_data)
# [1] 60  6

pheno_data$time <- as.numeric(grepl("Timepoint-2", pheno_data$title))
table(pheno_data$time)
# 0  1 
# 30 30 

# Step 2: Pull in the count data from GEO ----
x1 <- read.delim(input_dir)
x <- x1 # x will be used for analysis, x1 for raw counts
dim(x)
# [1] 60664    60

# Create a DGEList object
x <- DGEList(x)
geneid <- rownames(x) 
library(EnsDb.Hsapiens.v86)
genes <- mapIds(EnsDb.Hsapiens.v86, keys=geneid, column=c("SYMBOL"), 
                keytype="GENEID")
# Warning message:
  # Unable to map 3533 of 60664 requested IDs. 
x$genes <- genes
detach("package:EnsDb.Hsapiens.v86", unload = TRUE) ## stops masking dplyr functions

head(colnames(x$counts))
# [1] "X22_P_1001_1_M" "X22_P_1001_2_M" "X22_P_1003_1_t" "X22_P_1003_2_t" "X22_P_1004_1_M" "X22_P_1004_2_M"
pheno_data$match <- gsub("\\[|\\]", "",str_extract(pheno_data$title, "\\[(.*?)\\]"))
index <- match(substring(text = colnames(x$counts),first = 2,last = nchar(x = colnames(x$counts))), pheno_data$match)
pheno_data <- pheno_data[index,]
pheno_data <- remove_missing(pheno_data)
head(colnames(x$counts))
# [1] "X22_P_1001_1_M" "X22_P_1001_2_M" "X22_P_1003_1_t" "X22_P_1003_2_t" "X22_P_1004_1_M" "X22_P_1004_2_M"
head(pheno_data$participant)
# [1] "22-P-1001-M" "22-P-1001-M" "22-P-1003"   "22-P-1003"   "22-P-1004-M" "22-P-1004-M"
index <- match(pheno_data$match, substring(text = colnames(x$counts),first = 2,last = nchar(x = colnames(x$counts))))
x$counts <- x$counts[,index]
dim(x$counts)
# [1] 60664    60
x$samples <- x$samples[index,]
dim(x$samples)
# [1] 60  3
identical(colnames(x$counts), rownames(x$samples))
# [1] TRUE

### Add group information to the metadata of our DGEList.
x$samples$group <- pheno_data$condition
x$samples$time <- pheno_data$time
x$samples$sex <- pheno_data$sex
x$samples$age <- pheno_data$age
x$samples$participant <- pheno_data$participant
x$samples$geo <- pheno_data$geo_accession

x$samples$time <- as.factor(x$samples$time)

### Filter out age-related cognitive disorders (ACD)
index <- which(x$samples$group == "ACD")
x$samples <- x$samples[-index,]
x$counts <- x$counts[,-index]
dim(x$counts)
# [1] 60664    28
identical(colnames(x$counts), rownames(x$samples))
# [1] TRUE
table(x$samples$participant, x$samples$time)
#              0 1
# 22-P-1001-M 1 1
# 22-P-1004-M 1 1
# 22-P-1013   1 1
# 22-P-1020   1 1
# 22-P-1021   1 1
# 22-P-1022   1 1
# 22-P-1056   1 1
# 22-P-1060   1 1
# 22-P-1062   1 1
# 22-P-1065   1 1
# 22-P-1069   1 1
# 22-P-1079   1 1
# 22-P-1103   1 1
# 22-P-1120   1 1

# Step 3: Filter out lowly expressed genes ----
## For the pre-filtering plot, save unfiltered log CPM values
lcpm <- cpm(x, log = TRUE)
## Filter by expr
keep.exprs <- filterByExpr(x, group=x$samples$time) 
x$counts <- x$counts[keep.exprs,] 
dim(x)
# [1] 15416    28

# Plot the difference.
L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
# Set the log-cpm cutoff
lcpm.cutoff <- log2(10/M + 2/L)

pdf(file = paste(output_dir, "filtered_counts.pdf", sep = ""))
### Code to generate plot for filtered genes
## Load the RColorBrewer library
library(RColorBrewer)
# Calculate the number of samples (nsamples)
nsamples <- ncol(x)
# Generate color palette
col <- brewer.pal(nsamples, "Paired")
# Set the plot layout
par(mfrow=c(1,2))
# Plot density for raw data
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
# Loop through samples and plot density
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

# Transform the new, filtered data to log scale
lcpm <- cpm(x, log=TRUE)
# Plot density for filtered data
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
# Loop through samples and plot density
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
## Save
dev.off()

# Step 3: Normalize ----
pdf(file = paste(output_dir, "lib_sizes.pdf", sep = ""))
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors
# [1] 0.8868048 0.9796594 0.8275745 0.8016391 0.9803654 0.9830938 1.0379586 1.0887652 1.1860290 1.0659665 1.0450978
# [12] 1.1100882 0.9473309 1.0234163 1.0105564 1.0110818 1.0124836 1.0631334 1.0512861 1.1113776 1.0323291 1.0080275
# [23] 1.0098306 1.1854714 0.9462405 0.9909154 0.8616843 0.8661332
range(x$samples$norm.factors)
# [1] 0.8016391 1.1860290
par(mfrow=c(1,2))
barplot(x$samples$lib.size, main = "library size")
barplot(lcpm, main = "lcpm lib sizes")
dev.off()

## Step 4: Unsupervised clustering analysis with PCAtools for QC ----
p <- pca(lcpm, metadata = x$samples, center = T, scale = T) 
## Change ENSG IDs to gene symbols.
p$xvars <- make.names(unname(mapIds(x = Homo.sapiens, keys = rownames(p$loadings), column = "SYMBOL", keytype = "ENSEMBL")), unique = T)
rownames(p$loadings) <- make.names(unname(mapIds(x = Homo.sapiens, keys = rownames(p$loadings), column = "SYMBOL", keytype = "ENSEMBL")), unique = T)

scree <- screeplot(p)
ptime <- biplot(p, 
                showLoadings = F,
                lab = NULL, 
                colby = "time", colkey = c("0" = "darkblue", "1" = "lightblue"),
                legendPosition = 'right', legendLabSize = 12, legendIconSize = 3,
                title = 'Gomez baseline vs. 50 minutes',
                subtitle = 'PC1 versus PC2')
psex <- biplot(p, 
               showLoadings = F,
               lab = NULL, 
               colby = "sex", 
               legendPosition = 'right', legendLabSize = 12, legendIconSize = 3)
psub <- biplot(p, 
               showLoadings = F,
               lab = NULL, 
               colby = "participant", 
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

# Step 5: voom transform the data for limma ----
## First, we produce a design matrix.

id <- factor(x$samples$participant)
time <- factor(x$samples$time)
sex <- factor(x$samples$sex)
TS <- paste(sex,time, sep = ".")
TS <- factor(TS)

design <- model.matrix(~0+TS)
colnames(design) <- make.names(levels(TS))

pdf(file = paste(output_dir, "mean_var_trend.pdf", sep = ""))
v <- voom(x, design, plot = TRUE)
dev.off()

# Step 6: Build contrasts, run limma ----
contrast_matrix <- makeContrasts(Time = (Female.1 + Male.1)/2 - (Female.0 + Male.0)/2,
                            TimeMale = (Male.1 - Male.0) - (Female.1 - Female.0),
                            levels = colnames(design))

cor <- duplicateCorrelation(v$E, design, block= id)
cor$consensus.correlation
# [1] 0.4306578

fit <- lmFit(object = v, design = design, block = id, correlation = cor$consensus.correlation)
fit2 <- contrasts.fit(fit, contrasts = contrast_matrix)
e <- eBayes(fit2)
pdf(file = paste(output_dir, "sa_plot.pdf", sep = ""))
plotSA(e, main = "Gomez: standardized residuals vs average expression")
dev.off()

# Step 7: View results ----
y <- topTable(e, number = Inf, sort.by = "F", adjust.method = "fdr")

yt <- topTable(e, coef = 1, number = Inf, adjust.method = "fdr")
ys <- topTable(e, coef = 2, number = Inf, adjust.method = "fdr")

dt <- decideTests(e, adjust.method = "none")
summary(dt)
#         Time TimeMale
# Down     305      269
# NotSig 14740    14915
# Up       371      232
dt <- decideTests(e, adjust.method = "BH")
summary(dt)
#         Time TimeMale
# Down       0        0
# NotSig 15416    15416
# Up         0        0

## Fixing symbols
library(EnsDb.Hsapiens.v86)
names(y)[1] <- "Symbol"
names(yt)[1] <- "Symbol"
names(ys)[1] <- "Symbol"
y$Symbol <- mapIds(EnsDb.Hsapiens.v86, keys = rownames(y), column = "SYMBOL", keytype = "GENEID")
yt$Symbol <- mapIds(EnsDb.Hsapiens.v86, keys = rownames(yt), column = "SYMBOL", keytype = "GENEID")
ys$Symbol <- mapIds(EnsDb.Hsapiens.v86, keys = rownames(ys), column = "SYMBOL", keytype = "GENEID")
# Warning message:
#   Unable to map 236 of 15416 requested IDs. 

write.csv(y, paste(output_dir, "limma_F.csv",sep=""), row.names=T)
write.csv(yt, paste(output_dir, "limma_time_main.csv",sep=""), row.names=T)
write.csv(ys, paste(output_dir, "limma_time_sex_interaction.csv",sep=""), row.names=T)

## Adjusted and unadjusted p-value histograms ----
library(cowplot)
a <- ggplot(data = yt, aes(x = adj.P.Val)) +
  geom_histogram(fill = "lightblue", color = "black") +
  ggtitle(paste("Adj. p histogram", names(y)[2])) +
  theme_cowplot() +
  theme(title = element_text(size = 10))
b <- ggplot(data = ys, aes(x = adj.P.Val)) +
  geom_histogram(fill = "lightblue", color = "black") +
  ggtitle(paste("Adj. p histogram", names(y)[3])) +
  theme_cowplot() +
  theme(title = element_text(size = 10))
c <- ggplot(data = yt, aes(x = P.Value)) +
  geom_histogram(fill = "lightblue", color = "black") +
  ggtitle(paste("P-value histogram", names(y)[2])) +
  theme_cowplot() +
  theme(title = element_text(size = 10))
d <-  ggplot(data = ys, aes(x = P.Value)) +
  geom_histogram(fill = "lightblue", color = "black") +
  ggtitle(paste("P-value histogram", names(y)[3])) +
  theme_cowplot() +
  theme(title = element_text(size = 10))
pdf(file = paste(output_dir, "p_value_histograms.pdf", sep = ""))
plot_grid(a,b,c,d, labels = LETTERS[1:4], nrow = 2, align = 'hv')
dev.off()


## Creating volcano plots ----
library(EnhancedVolcano)
# Time one and two
volcano1 <- EnhancedVolcano(yt,
                            lab = yt$Symbol,
                            x = "logFC",
                            y = "P.Value",
                            title = ("Gomez: 50-minute samples"), titleLabSize = 12,
                            subtitle = "Main effect of time", caption = NULL,
                            labSize = (4.0),
                            col = c("black","lightblue","royalblue","blue"),
                            colAlpha = 1,
                            legendLabSize = 8, legendIconSize = 3,
                            pointSize = 1.0,
                            boxedLabels = FALSE,
                            FCcutoff = 1,
                            pCutoff = 0.05,
                            drawConnectors = TRUE)
volcano1
volcano2 <- EnhancedVolcano(ys,
                            lab = ys$Symbol,
                            x = "logFC",
                            y = "P.Value",
                            title = (""), titleLabSize = 12,
                            subtitle = "Male to female interaction term", caption = NULL,
                            labSize = (4.0),
                            col = c("black","lightblue","royalblue","blue"),
                            colAlpha = 1,
                            legendLabSize = 8, legendIconSize = 3,
                            pointSize = 1.0,
                            boxedLabels = FALSE,
                            FCcutoff = 1,
                            pCutoff = 0.05,
                            drawConnectors = TRUE)
volcano2
pdf(file = paste(output_dir, "volcano_plots.pdf", sep = ""), width = 10)
plot_grid(volcano1, volcano2, align = "h", nrow = 1)
dev.off()