## Date: Jul 11 2024

### SUMMARY: Traces of Time paper, Gosch et al 2023

# Updated May 13 2024: Removing sample from Patient J at 23 hours.

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
library(RColorBrewer)
library(EnhancedVolcano)

## INPUT: Raw count file from GSE223613 on Gene Expression Omnibus ----
input_dir <- "~/Desktop/work_repo/data/GSE223613_counts.txt"
## OUTPUT: Limma top table, quality control plots, volcano plots, and p-value histograms ----
output_dir <- "~/Desktop/work_repo/github/Gosch_study_data/output/"

# Step 1: pull sample information ----
library(GEOquery)
gset <- getGEO(GEO = "GSE223613")
p <- phenoData(gset$GSE223613_series_matrix.txt.gz)
pheno_data_raw <- p@data
dim(pheno_data_raw)
# [1] 80 49
colnames(pheno_data_raw)
# [1] "title"                   "geo_accession"           "status"                  "submission_date"        
# [5] "last_update_date"        "type"                    "channel_count"           "source_name_ch1"        
# [9] "organism_ch1"            "characteristics_ch1"     "characteristics_ch1.1"   "characteristics_ch1.2"  
# [13] "characteristics_ch1.3"   "characteristics_ch1.4"   "treatment_protocol_ch1"  "growth_protocol_ch1"    
# [17] "molecule_ch1"            "extract_protocol_ch1"    "extract_protocol_ch1.1"  "extract_protocol_ch1.2" 
# [21] "taxid_ch1"               "description"             "data_processing"         "data_processing.1"      
# [25] "data_processing.2"       "data_processing.3"       "data_processing.4"       "data_processing.5"      
# [29] "platform_id"             "contact_name"            "contact_department"      "contact_institute"      
# [33] "contact_address"         "contact_city"            "contact_zip/postal_code" "contact_country"        
# [37] "data_row_count"          "instrument_model"        "library_selection"       "library_source"         
# [41] "library_strategy"        "relation"                "relation.1"              "supplementary_file_1"   
# [45] "age:ch1"                 "participant:ch1"         "sampling time:ch1"       "Sex:ch1"                
# [49] "tissue:ch1" 
pheno_data <- pheno_data_raw %>% 
  dplyr::select(title = title,
         geo_accession = geo_accession,
         age = `age:ch1`,
         participant = `participant:ch1`,
         time = `sampling time:ch1`,
         sex = `Sex:ch1`,
         sample = description)
dim(pheno_data)
# [1] 80  7

# Step 2: Pull in the counts ----
x1 <- read.delim(input_dir)
x <- x1 # x will be used for analysis, x1 for raw counts
dim(x)
# [1] 33121    81
gene_list <- x[,1]
rownames(x) <- make.names(gene_list, unique = TRUE)
x <- x[,-1]
identical(pheno_data$sample, sub(".","",sub("\\.","-",names(x))))
# [1] TRUE
names(x) <- pheno_data$title
x <- DGEList(x)
x[["samples"]]
x$samples$group <- as.factor(pheno_data$participant)
x$samples$time <- as.factor(pheno_data$time)
x$samples$sex <- as.factor(pheno_data$sex)

## Remove Patient J at 23 hours. Library size: 401784615
pdf(file = paste(output_dir, "lib_size_before.pdf", sep = ""))
barplot(x$samples$lib.size, names=substr(colnames(x$counts), start = 19, stop = 28), las=2)
title("Barplot of library sizes")
dev.off()
which(x$samples$lib.size == 401784615)
# [1] 60
pheno_data[60,]
#                                  title geo_accession age participant   time    sex                                sample
# GSM6965634 Blood, ParticipantJ, 23.00h    GSM6965634  22           J 23.00h Female 21Aug27_F10-L1_S96_L004.featureCounts
colnames(x$counts)[60] == rownames(x$samples)[60]
# [1] TRUE
x$counts <- x$counts[,-60]
x$samples <- x$samples[-60,]
dim(x$counts)
# [1] 33121    79
dim(x$samples)
# [1] 79  5

# Step 3: Filter out lowly expressed genes ----
# Create an object for plots before filter
lcpm <- cpm(x, log = TRUE)

## from RNA-seq analysis is easy as 1-2-3...
keep.exprs <- filterByExpr(x, group=x$samples$group) 
x <- x[keep.exprs,, keep.lib.sizes=FALSE] 
dim(x)
# [1] 16389    79

L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
# Set the log-cpm cutoff
lcpm.cutoff <- log2(10/M + 2/L)

pdf(file = paste(output_dir, "filtered_counts.pdf", sep = ""))
### Code to generate plot for filtered genes
# Calculate the number of samples (nsamples)
nsamples <- ncol(x)
# Set sample names
samplenames <- x$samples$group
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
# Plot density for filtered data
lcpm <- cpm(x, log = TRUE)
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
# [1] 1.0184021 0.9149087 1.0024654 1.0841515 0.9866103 1.0357074 0.9969307 1.0172417 0.7927552 1.1115282 0.9833401 0.8450658
# [13] 0.9236970 1.0382787 0.8727475 1.0104439 0.9968386 1.0047867 0.7887432 1.0278945 0.9788101 0.8686775 0.9271797 1.0502645
# [25] 0.9031842 1.0144149 0.9983534 1.0569251 0.7785320 1.0203860 1.0104843 0.8862879 0.9717633 1.0830962 0.9617000 1.0427094
# [37] 1.0102041 1.0247231 0.8608458 1.0732166 1.0037541 0.9664944 1.0214881 1.1030054 0.9748720 1.0323586 0.9944404 1.0469698
# [49] 0.8879460 1.1054414 1.0308126 0.9220745 1.0271236 1.1214018 1.0217899 1.0473623 1.0150967 1.0298751 0.9319958 1.0477921
# [61] 0.9718627 1.0464417 1.1751787 1.0365402 1.0510682 1.0241392 1.0539804 0.9896065 1.1094076 1.0419033 0.9901191 1.0463826
# [73] 1.1444592 0.9919014 1.0378734 1.0653931 1.0556027 1.0153664 1.0939273
range(x$samples$norm.factors)
# [1] 0.778532 1.175179

## Create a barplot of library sizes
pdf(file = paste(output_dir, "lib_size_after.pdf", sep = ""))
par(mfrow = c(1,2))
barplot(x$samples$lib.size, names=substr(colnames(x$counts), start = 19, stop = 28), las=2)
title("Barplot of library sizes")
barplot(colSums(lcpm), names=substr(colnames(x$counts), start = 19, stop = 28), las=2)
title("lcpm library sizes")
dev.off()

# Step 4: PCA and MDS plots ----
p <- pca(x$counts, metadata = x$samples, center = T, scale = T)

scree <- screeplot(p, components = getComponents(p, 1:20))
ptime <- biplot(p, 
                showLoadings = F,
                lab = NULL, 
                colby = "time",
                legendPosition = 'right', legendLabSize = 12, legendIconSize = 3,
                title = 'Gosch 24 hours',
                subtitle = 'PC1 versus PC2')
psex <- biplot(p, 
               showLoadings = F,
               lab = NULL, 
               colby = "sex", 
               legendPosition = 'right', legendLabSize = 12, legendIconSize = 3)
psub <- biplot(p, 
               showLoadings = F,
               lab = NULL, 
               colby = "group", subtitle = "group = participant that samples belonged to",
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
colors <- brewer.pal(n = 10, name = "Paired")
colors <- rep(colors, times = 8)
colors <- colors[-60] # Since patient sample was removed
pdf(file = paste(output_dir, "mds.pdf", sep = ""))
limma::plotMDS(x = lcpm, 
               cex = 1, 
               var.explained = T, 
               col = colors,
               labels = x$samples$participant)
title("MDS Plot")
dev.off()

# Step 5: voom transform the data for limma ----
## First, we produce a design matrix.

id <- factor(x$samples$group)
time <- factor(x$samples$time)
levels(time) <- c("02.00h", "05.00h", "08.00h", "11.00h", "14.00h", "17.00h", "20.00h", "23.00h")
sex <- factor(x$samples$sex)
TS <- paste(sex,time, sep = ".")
TS <- factor(TS)

design <- model.matrix(~0+TS)
colnames(design) <- make.names(levels(TS))

pdf(file = paste(output_dir, "mean_variance.pdf", sep = ""))
v <- voom(x, design, plot = TRUE)
dev.off()

# Step 6: Build contrasts, run limma ----

constrasts_twohrs_all <- makeContrasts(# Main Effect of Time
                            ThreeHrMain = (Female.05.00h + Male.05.00h)/2 - (Female.02.00h + Male.02.00h)/2,
                            SixHrMain = (Female.08.00h + Male.08.00h)/2 - (Female.02.00h + Male.02.00h)/2,
                            NineHrMain = (Female.11.00h + Male.11.00h)/2 - (Female.02.00h + Male.02.00h)/2,
                            TwelveHrMain = (Female.14.00h + Male.14.00h)/2 - (Female.02.00h + Male.02.00h)/2,
                            FifteenHrMain = (Female.17.00h + Male.17.00h)/2 - (Female.02.00h + Male.02.00h)/2,
                            EighteenHrMain = (Female.20.00h + Male.20.00h)/2 - (Female.02.00h + Male.02.00h)/2,
                            TwentyoneHrMain = (Female.23.00h + Male.23.00h)/2 - (Female.02.00h + Male.02.00h)/2,
                            # Time-sex interaction effects
                            ThreeHrMale = (Male.05.00h - Male.02.00h) - (Female.05.00h - Female.02.00h),
                            SixHrMale = (Male.08.00h - Male.02.00h) - (Female.08.00h - Female.02.00h),
                            NineHrMale = (Male.11.00h - Male.02.00h) - (Female.11.00h - Female.02.00h),
                            TwelveHrMale = (Male.14.00h - Male.02.00h) - (Female.14.00h - Female.02.00h),
                            FifteenHrMale = (Male.17.00h - Male.02.00h) - (Female.17.00h - Female.02.00h),
                            EighteenHrMale = (Male.20.00h - Male.02.00h) - (Female.20.00h - Female.02.00h),
                            TwentyoneHrMale = (Male.23.00h - Male.02.00h) - (Female.23.00h - Female.02.00h),
                            levels = colnames(design))
cor <- duplicateCorrelation(v$E, design, block=id)
cor$consensus.correlation
# [1] 0.6236048

fit <- lmFit(object = v, design = design, block = id, correlation = cor$consensus.correlation)
fit2 <- contrasts.fit(fit, contrasts = constrasts_twohrs_all)
e <- eBayes(fit2)
pdf(file = paste(output_dir, "sa_plot.pdf", sep = ""))
plotSA(e, main="Residual standard deviation versus average log expression")
dev.off()

# Step 8: View results ----
y <- topTable(e, number = Inf, adjust.method = "fdr", sort.by = "F")
dt <- decideTests(e)
summary(dt)
#        ThreeHrMain SixHrMain NineHrMain TwelveHrMain FifteenHrMain EighteenHrMain TwentyoneHrMain ThreeHrMale SixHrMale NineHrMale
# Down             0      1332       3601         3250          1192            299               6           0         0          0
# NotSig       16389     13482       8968         9603         13506          15746           16382       16389     16389      16389
# Up               0      1575       3820         3536          1691            344               1           0         0          0
#        TwelveHrMale FifteenHrMale EighteenHrMale TwentyoneHrMale
# Down              0             0              0               0
# NotSig        16389         16389          16389           16389
# Up                0             0              0               0
dt <- decideTests(e, adjust.method = "none")
summary(dt)
#        ThreeHrMain SixHrMain NineHrMain TwelveHrMain FifteenHrMain EighteenHrMain TwentyoneHrMain ThreeHrMale SixHrMale NineHrMale
# Down           349      2437       4225         3927          2405           1654             873         272       304        384
# NotSig       15711     11159       7906         8380         11177          12822           14583       15717     15814      15616
# Up             329      2793       4258         4082          2807           1913             933         400       271        389
#        TwelveHrMale FifteenHrMale EighteenHrMale TwentyoneHrMale
# Down            672           236            238             255
# NotSig        15144         15928          15911           15953
# Up              573           225            240             181

write.csv(y, paste(output_dir, "limma_F.csv", sep = ""), row.names=T)

tables <- list()
for (i in 1:14) {
  tables[[i]] <- topTable(e, coef = i, number = Inf, adjust.method = "fdr")
  write.csv(tables[[i]], file = paste(output_dir,"limma_",names(y)[i],".csv", sep = ""))
}

## P-value histograms ----
histlist <- list()
for (i in 1:7) {
  histlist[[i]] <- ggplot(data = tables[[i]], aes(x = adj.P.Val)) +
    geom_histogram(fill = "lightblue", color = "black") +
    ggtitle(paste("Adj. p histogram", names(y)[i])) +
    theme_cowplot() +
    theme(title = element_text(size = 7.5))
}
histlist[[8]] <- ggplot(data = y, aes(x = adj.P.Val)) +
  geom_histogram(fill = "lightblue", color = "black") +
  ggtitle(paste("Adj. p histogram F-statistic")) +
  theme_cowplot() +
  theme(title = element_text(size = 7.5))
pdf(file = paste(output_dir, "adjusted_p_hist_main.pdf", sep = ""))
plot_grid(plotlist = histlist)
dev.off()

histlist <- list()
for (i in 1:7) {
  histlist[[i]] <- ggplot(data = tables[[i]], aes(x = P.Value)) +
    geom_histogram(fill = "lightblue", color = "black") +
    ggtitle(paste("P-value histogram", names(y)[i])) +
    theme_cowplot() +
    theme(title = element_text(size = 7.5))
}
histlist[[8]] <- ggplot(data = y, aes(x = P.Value)) +
  geom_histogram(fill = "lightblue", color = "black") +
  ggtitle(paste("P-value overall F-statistic")) +
  theme_cowplot() +
  theme(title = element_text(size = 7.5))
pdf(file = paste(output_dir, "p_hist_main.pdf", sep = ""))
plot_grid(plotlist = histlist)
dev.off()

histlist <- list()
for (i in 1:7) {
  histlist[[i]] <- ggplot(data = tables[[i + 7]], aes(x = P.Value)) +
    geom_histogram(fill = "lightblue", color = "black") +
    ggtitle(paste("P-value histogram", names(y)[i + 7])) +
    theme_cowplot() +
    theme(title = element_text(size = 7.5))
}
pdf(file = paste(output_dir, "p_hist_time_sex_int.pdf", sep = ""))
plot_grid(plotlist = histlist, ncol = 2)
dev.off()

histlist <- list()
for (i in 1:7) {
  histlist[[i]] <- ggplot(data = tables[[i + 7]], aes(x = adj.P.Val)) +
    geom_histogram(fill = "lightblue", color = "black") +
    ggtitle(paste("Adj. P-value histogram", names(y)[i + 7])) +
    theme_cowplot() +
    theme(title = element_text(size = 7.5))
}
pdf(file = paste(output_dir, "adj_p_hist_time_sex_int.pdf", sep = ""))
plot_grid(plotlist = histlist, ncol = 2)
dev.off()

## Volcano plots ----
volcanos <- list()
for (i in 1:7) {
  volcanos[[i]] <- EnhancedVolcano(toptable = tables[[i]], 
                              lab = rownames(tables[[i]]), 
                              x = "logFC",
                              y = "P.Value",
                              title = paste(names(y)[i], "p-value volcano"), titleLabSize = 8.5,
                              subtitle = "main effect of time", subtitleLabSize = 5.5,
                              caption = "",
                              labSize = (4.0),
                              colAlpha = 1,
                              pointSize = 1.0, 
                              pCutoff = 1e-2,
                              boxedLabels = FALSE,
                              drawConnectors = TRUE)
}
toprow <- plot_grid(volcanos[[1]], volcanos[[2]], volcanos[[3]], volcanos[[4]], nrow = 1, labels = LETTERS[1:4])
botrow <-  plot_grid(volcanos[[5]], volcanos[[6]], volcanos[[7]], nrow = 1, labels = LETTERS[5:7])
pdf(file = paste(output_dir, "main_volcanos.pdf", sep = ""), width = 12, height = 10)
plot_grid(toprow, botrow, nrow = 2)
dev.off()

volcanos <- list()
for (i in 1:7) {
  volcanos[[i]] <- EnhancedVolcano(toptable = tables[[i + 7]], 
                                   lab = rownames(tables[[i + 7]]), 
                                   x = "logFC",
                                   y = "P.Value",
                                   title = paste(names(y)[i], "p-value volcano"), titleLabSize = 8.5,
                                   subtitle = "time-sex interaction, male - female", subtitleLabSize = 5.5,
                                   caption = "",
                                   labSize = (4.0),
                                   colAlpha = 1,
                                   pointSize = 1.0, 
                                   pCutoff = 1e-2,
                                   boxedLabels = FALSE,
                                   drawConnectors = TRUE)
}
toprow <- plot_grid(volcanos[[1]], volcanos[[2]], volcanos[[3]], volcanos[[4]], nrow = 1, labels = LETTERS[1:4])
botrow <-  plot_grid(volcanos[[5]], volcanos[[6]], volcanos[[7]], nrow = 1, labels = LETTERS[5:7])
pdf(file = paste(output_dir, "time_sex_volcanos.pdf", sep = ""), width = 12, height = 10)
plot_grid(toprow, botrow, nrow = 2)
dev.off()

### Supplementary Figure S3-A ----
volcanos <- list()
for (i in 1:7) {
  volcanos[[i]] <- EnhancedVolcano(toptable = tables[[i]], 
                                   lab = "", 
                                   x = "logFC",
                                   y = "P.Value",
                                   title = "",
                                   subtitle = "",
                                   caption = "", legendPosition = 'none',
                                   labSize = (4.0),
                                   colAlpha = 1,
                                   pointSize = 1.0, 
                                   pCutoff = 1e-2,
                                   xlim = c(-4,4),
                                   ylim = c(0,25),
                                   boxedLabels = FALSE)
}
toprow <- plot_grid(volcanos[[1]], volcanos[[2]], volcanos[[3]], volcanos[[4]], nrow = 1)
botrow <-  plot_grid(plotlist = list(volcanos[[5]], volcanos[[6]], volcanos[[7]], NULL), nrow = 1)
pdf(file = paste(output_dir, "figureS3A.pdf", sep = ""), width = 12, height = 10)
plot_grid(toprow, botrow, nrow = 2, align = 'hv')
dev.off()

### Supplementary Figure S3-B ----
tables <- tables[1:7]
degs <- lapply(tables, function(x){
  x %>% filter(P.Value < 0.05) %>% nrow() %>% return()
})
degstb <- tibble(time = factor(c(3,6,9,12,15,18,21), levels = c(3,6,9,12,15,18,21), ordered = T), nums = unlist(degs))
pdf(file = paste(output_dir, "figureS3B.pdf", sep = ""), width = 8, height = 8)
ggplot(degstb, aes(x = time, y = nums)) + 
  geom_col(fill = "black", color = "black") +
  labs(title = "", subtitle = "",
       x = "hours from baseline", y = "DEGs (p < 0.05)") +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
dev.off()

## Diurnal only for the F-statistic ----
diurnal_contrasts <- makeContrasts(# Main Effect of Time ONLY
  ThreeHrMain = (Female.05.00h + Male.05.00h)/2 - (Female.02.00h + Male.02.00h)/2,
  SixHrMain = (Female.08.00h + Male.08.00h)/2 - (Female.02.00h + Male.02.00h)/2,
  NineHrMain = (Female.11.00h + Male.11.00h)/2 - (Female.02.00h + Male.02.00h)/2,
  TwelveHrMain = (Female.14.00h + Male.14.00h)/2 - (Female.02.00h + Male.02.00h)/2,
  FifteenHrMain = (Female.17.00h + Male.17.00h)/2 - (Female.02.00h + Male.02.00h)/2,
  EighteenHrMain = (Female.20.00h + Male.20.00h)/2 - (Female.02.00h + Male.02.00h)/2,
  TwentyoneHrMain = (Female.23.00h + Male.23.00h)/2 - (Female.02.00h + Male.02.00h)/2,
  levels = colnames(design))
fit2 <- contrasts.fit(fit, contrasts = diurnal_contrasts)
e <- eBayes(fit2)
y <- topTable(e, number = Inf, adjust.method = "fdr", sort.by = "F")
write.csv(y, paste(output_dir, "diurnal_genes_limma_F.csv", sep = ""), row.names=T)