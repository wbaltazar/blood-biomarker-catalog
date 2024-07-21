## Date: Jul 11 2024

### SUMMARY: Do responders/non-responders to exercise induced changes in VO2 max have different
# transcriptional profiles? Data from 30 women sampled at baseline and 4 months after random
# assignment of exercise regimens of varying intensities. Subjects between robust responders (RR)
# and non-responders (NR) were matched as reasonably as possible.

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
library(EnhancedVolcano)
library(RColorBrewer)

## INPUT: Raw count file from GSE206505 on Gene Expression Omnibus ----
input_dir <- "~/Desktop/work_repo/data/GSE206505_VO2-counts-matrix.xlsx"
## OUTPUT: Limma top table, quality control plots, volcano plots, and p-value histograms ----
output_dir <- "~/Desktop/work_repo/github/LaRocca_study_data/output/"

# Step 1: pull sample information ----
gset <- getGEO(GEO = "GSE206505")
p <- phenoData(gset$GSE206505_series_matrix.txt.gz)
pheno_data_raw <- p@data
dim(pheno_data_raw)
# [1] 60 42
colnames(pheno_data_raw)
# [1] "title"                   "geo_accession"           "status"                  "submission_date"        
# [5] "last_update_date"        "type"                    "channel_count"           "source_name_ch1"        
# [9] "organism_ch1"            "characteristics_ch1"     "characteristics_ch1.1"   "characteristics_ch1.2"  
# [13] "characteristics_ch1.3"   "molecule_ch1"            "extract_protocol_ch1"    "extract_protocol_ch1.1" 
# [17] "taxid_ch1"               "data_processing"         "data_processing.1"       "data_processing.2"      
# [21] "data_processing.3"       "platform_id"             "contact_name"            "contact_email"          
# [25] "contact_institute"       "contact_address"         "contact_city"            "contact_state"          
# [29] "contact_zip/postal_code" "contact_country"         "data_row_count"          "instrument_model"       
# [33] "library_selection"       "library_source"          "library_strategy"        "relation"               
# [37] "relation.1"              "supplementary_file_1"    "cell type:ch1"           "pre or post:ch1"        
# [41] "r or nr:ch1"             "Sex:ch1" 
## We have reason to believe that the source name is improperly labeled in this data. According to paper
## the blood was collected in PAXgene tubes, making PBMC analysis impossible
pheno_data <- pheno_data_raw %>% 
  dplyr::select(title = title,
         geo_accession = geo_accession,
         # age = `age:ch1`, may need to search for this info. Ages were between 30-45
         # participant = `subject id:ch1`, this is not recorded. the table suggests that the info is organized by participant.
         time = `pre or post:ch1`,
         sex = `Sex:ch1`,
         response = `r or nr:ch1`)
dim(pheno_data)
# [1] 60  5
pheno_data$participant <- rep(1:30, each = 2)
pheno_data$participant <- as.factor(pheno_data$participant)
pheno_data$participant

# Step 2: Pull in the counts ----
x1 <- readxl::read_xlsx(input_dir)
x <- x1 # x stores raw counts
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
detach("package:EnsDb.Hsapiens.v86", unload = TRUE) ## unmask dplyr functions

## Ensure that the metadata is aligned with the count data
identical(str_replace(colnames(x$counts), "X", ""), pheno_data[-grep(21, pheno_data$title),"title"])
# [1] TRUE
pheno_data <- pheno_data[-grep(21, pheno_data$title),]
colnames(x$counts) <- pheno_data$geo_accession

### Add group information to the samples table.
x$samples$group <- pheno_data$response
x$samples$time <- pheno_data$time
x$samples$sex <- pheno_data$sex
x$samples$participant <- pheno_data$participant
x$samples$geo <- pheno_data$geo_accession


# Step 3: Filter out lowly expressed genes ----
lcpm <- cpm(x, log = TRUE) # Stored for pre-filter plots
keep.exprs <- filterByExpr(x, group=x$samples$time) 
x$counts <- x$counts[keep.exprs,] 
dim(x)
# [1] 15555    59

L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
# Set the log-cpm cutoff
lcpm.cutoff <- log2(10/M + 2/L)

pdf(file = paste(output_dir, "filtered_counts.pdf"))
### Code to generate plot for filtered genes
# Calculate the number of samples (nsamples)
nsamples <- ncol(x)
# Set sample names
samplenames <- x$samples$participant
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

# Now use the filtered data.
lcpm <- cpm(x, log = TRUE)
# Plot density for filtered data
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
# Loop through samples and plot density
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
dev.off()

# Step 3: Normalize ----
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors
# [1] 2.2287677 1.6190875 1.7846828 1.9006497 1.6567835 1.3588008 1.3682596 0.6970378 1.0299705 1.5197897 2.0550928
# [12] 1.5925446 0.9181631 1.2960379 1.1868339 1.5669682 1.8425490 1.7785557 1.4320234 0.2585464 0.5640930 1.0981034
# [23] 0.4371921 0.4952125 0.3877064 0.3344057 0.6296602 0.4160134 0.9540915 0.5443197 1.1778069 0.5424408 1.5275817
# [34] 0.5923576 1.3911888 0.5604255 0.6028132 0.6163293 0.5681378 0.2823831 0.4266784 0.7275542 1.2544331 1.2061025
# [45] 1.2507532 1.5668099 1.5945944 0.9664580 1.6642729 1.5475576 0.6543195 1.7358335 1.1764426 0.9732007 1.2972947
# [56] 1.4328381 1.4951922 1.7746299 1.8576358
range(x$samples$norm.factors)
# [1] 0.2585464 2.2287677

pdf(file = paste(output_dir, "library_sizes.pdf"))
par(mfrow = c(1,2))
barplot(x$samples$lib.size, main = "library size")
barplot(colSums(lcpm), col=col, main = "lcpm lib sizes")
dev.off()

# Step 4: MDS & PCA ----
# PCAtools
library(PCAtools) #response, time, participant
p <- pca(x$counts, metadata = data.frame(x$samples, row.names = colnames(x$counts)), center = T, scale = T)

scree <- screeplot(p, components = getComponents(p, 1:20))
ptime <- biplot(p, 
                showLoadings = F,
                lab = NULL, 
                colby = "time", colkey = c("pre" = "lightblue2", "post" = "darkblue"),
                legendPosition = 'right', legendLabSize = 12, legendIconSize = 3,
                title = 'LaRocca baseline vs. 4 months',
                subtitle = 'PC1 versus PC2')
pres <- biplot(p, 
               showLoadings = F,
               lab = NULL, 
               colby = "group", subtitle = "group = responder status", 
               legendPosition = 'right', legendLabSize = 12, legendIconSize = 3)
psub <- biplot(p, 
               showLoadings = F,
               lab = NULL, 
               colby = "participant",
               legendPosition = 'right', legendLabSize = 9, legendIconSize = 3)
peigencor <- eigencorplot(p,
                          components = getComponents(p, 1:10),
                          metavars = c("time", 'group'),
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
toprow <- plot_grid(ptime, pres, psub, align = 'h', nrow = 1)
botrow <- plot_grid(scree, peigencor, align = 'h', nrow = 1)
pdf(file = paste(output_dir, "pca.pdf", sep = ""), height = 12, width = 16)
plot_grid(toprow, botrow, ncol = 1, align = 'v')
dev.off()

# MDS plot
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colors <- sample(col_vector, 30)
colors <- rep(colors, each = 2)
colors <- colors[-20] # Account for missing sample
pdf(file = paste(output_dir, "mds.pdf", sep = ""))
limma::plotMDS(x = lcpm, 
               cex = 1, 
               var.explained = T, 
               col = colors,
               labels = paste("P-",x$samples$participant, sep = ""))
title("MDS Plot")
dev.off()

# Step 5: voom transform the data for limma ----
## First, we produce a design matrix.

id <- factor(x$samples$participant)
time <- factor(x$samples$time)
response <- factor(x$samples$group)
TS <- paste(response,time, sep = ".")
TS <- factor(TS)

design <- model.matrix(~0+TS)
colnames(design) <- make.names(levels(TS))
pdf(file = paste(output_dir, "mean_variance.pdf", sep = ""))
v <- voom(x, design, plot = TRUE)
dev.off()
cor <- duplicateCorrelation(v$E, design, block=id)
cor$consensus.correlation
# [1] 0.2582235

# Step 6: Build contrasts, run limma ----
my_contrasts <- makeContrasts(Time = (R.post + NR.post)/2 - (R.pre + NR.pre)/2,
                              Time_response = (R.post - R.pre) - (NR.post - NR.pre),
                              Response = (R.post + R.pre)/2 - (NR.post + NR.pre)/2,
                              levels = levels(TS))

fit <- lmFit(object = v, design = design, block = id, correlation = cor$consensus.correlation)
fit2 <- contrasts.fit(fit, contrasts = my_contrasts)
e <- eBayes(fit2)
pdf(file = paste(output_dir, "sa_plot.pdf", sep = ""))
plotSA(e, main = "LaRocca: standardized residuals vs average expression")
dev.off()

# Step 7: View results ----
y <- topTable(e, number = Inf, sort.by = "F", adjust.method = "fdr")
dt <- decideTests(e)
summary(dt)
#         Time Time_response Response
# Down       0             0     1764
# NotSig 15554         15555    12925
# Up         1             0      866
dt <- decideTests(e, adjust.method = "none")
summary(dt)
#         Time Time_response Response
# Down     204           191     3299
# NotSig 15138         15273     9187
# Up       213            91     3069

library(EnsDb.Hsapiens.v86)
y$ProbeID <- mapIds(EnsDb.Hsapiens.v86, keys = rownames(y), column = "GENEID", keytype = "SYMBOL")
# Warning message:
# Unable to map 2251 of 15555 requested IDs. 
names(y)[1] <- "GeneID"
write.csv(y, paste(output_dir, "limma_F.csv",sep=""), row.names=TRUE)

tables <- list()
names(y)
# [1] "GeneID"       "Time"          "Time_response" "Response"      "AveExpr"       "F"             "P.Value"      
# [8] "adj.P.Val"    
for (i in 1:3) {
  tables[[i]] <- topTable(e, coef = i, number = Inf, adjust.method = "fdr")
  tables[[i]]$GeneID <- mapIds(EnsDb.Hsapiens.v86, keys = rownames(tables[[i]]), column = "GENEID", keytype = "SYMBOL")
  write.csv(tables[[i]], file = paste(output_dir, "limma_", names(y)[i + 1], ".csv", sep = ""))
}
detach("package:EnsDb.Hsapiens.v86", unload = TRUE)

## Plot p-value histograms ----
histlist <- list()
for (i in 1:3) {
  histlist[[i]] <- ggplot(data = tables[[i]], aes(x = adj.P.Val)) +
    geom_histogram(fill = "lightblue", color = "black") +
    ggtitle(paste("Adj. p histogram", names(y)[i+1])) +
    theme_cowplot() +
    theme(title = element_text(size = 7.5))
}
histlist[[4]] <- ggplot(data = y, aes(x = adj.P.Val)) +
  geom_histogram(fill = "lightblue", color = "black") +
  ggtitle("Adj. p histogram of F-statistic") +
  theme_cowplot() +
  theme(title = element_text(size = 7.5))
for (i in 1:3) {
  histlist[[i + 4]] <- ggplot(data = tables[[i]], aes(x = P.Value)) +
    geom_histogram(fill = "lightblue", color = "black") +
    ggtitle(paste("P-value histogram", names(y)[i + 1])) +
    theme_cowplot() +
    theme(title = element_text(size = 7.5))
}
histlist[[8]] <- ggplot(data = y, aes(x = P.Value)) +
  geom_histogram(fill = "lightblue", color = "black") +
  ggtitle("P-value histogram of time F-statistic") +
  theme_cowplot() +
  theme(title = element_text(size = 7.5))
pdf(file = paste(output_dir, "p_value_histograms.pdf", sep = ""))
plot_grid(plotlist = histlist, nrow = 2, ncol = 4, labels = LETTERS[1:8], align = 'hv')
dev.off()

## Enhanced volcano plots ----
volcano1 <- EnhancedVolcano(tables[[1]],
                            lab = rownames(tables[[1]]),
                            x = "logFC",
                            y = "P.Value",
                            title = ("LaRocca: 4-month volcano plot of p-values"), titleLabSize = 8,
                            subtitle = "",
                            caption = "",
                            labSize = (4.0),
                            legendLabSize = 5, legendIconSize = 2,
                            col = c("black","lightblue","royalblue","blue"),
                            colAlpha = 1,
                            pointSize = 1.0, 
                            pCutoff = 1e-2,
                            boxedLabels = FALSE,
                            drawConnectors = TRUE)
volcano2 <- EnhancedVolcano(tables[[2]],
                            lab = rownames(tables[[2]]),
                            x = "logFC",
                            y = "P.Value",
                            title = ("LaRocca: 4-month interaction with response"), titleLabSize = 8,
                            subtitle = "",
                            caption = "",
                            labSize = (4.0),
                            legendLabSize = 5, legendIconSize = 2, 
                            col = c("black","lightblue","royalblue","blue"),
                            colAlpha = 1,
                            pointSize = 1.0,
                            pCutoff = 1e-2,
                            boxedLabels = FALSE,
                            drawConnectors = TRUE)
volcano3 <- EnhancedVolcano(tables[[3]],
                            lab = rownames(tables[[3]]),
                            x = "logFC",
                            y = "P.Value",
                            title = ("LaRocca: NR v. RR Response volcano"), titleLabSize = 8,
                            subtitle = "",
                            caption = "",
                            labSize = (4.0),
                            legendLabSize = 5, legendIconSize = 2,
                            col = c("black","lightblue","royalblue","blue"),
                            colAlpha = 1,
                            pointSize = 1.0,
                            pCutoff = 1e-2,
                            boxedLabels = FALSE,
                            drawConnectors = TRUE)
pdf(paste(output_dir, "volcano_plots.pdf", sep = ""), width = 10)
plot_grid(volcano1, volcano2, volcano3, align = 'h', labels = LETTERS[1:3], nrow = 1)
dev.off()