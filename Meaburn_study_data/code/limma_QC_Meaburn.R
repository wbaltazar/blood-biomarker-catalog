## Date: Oct 21 2024

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
input_dir <- "~/Desktop/work_repo/data/GSE14844_RAW/"
## OUTPUT: limma results, plots, figures
output_dir <- "~/Desktop/work_repo/github/Meaburn_study_data/output/"

## Get study metadata and separate by day ----
gset <- getGEO("GSE14844", GSEMatrix =TRUE, AnnotGPL=TRUE)
pheno_data_org <- pData(gset$GSE14844_series_matrix.txt.gz)
feature_data <- fData(gset$GSE14844_series_matrix.txt.gz)
names(pheno_data_org)

# filter out traits we want
pheno_data <- pheno_data_org %>% 
  dplyr::select(title = title,
                geo_accession = geo_accession,
                sex = `gender:ch1`)
pheno_data$time = str_extract(pheno_data_org$title, "Day.*")
pheno_data$id = str_extract(pheno_data$title, "TD.....")
pheno_data$time[12] <- "Daytwo_10am"
dim(pheno_data)
# [1] 36  5
table(pheno_data$sex, pheno_data$time)
#        Dayone_10am Dayone_2pm Daytwo_10am Daytwo_2pm
# Female           8          8           6          8
# Male             2          1           1          2
## Suffices for analysis?
p1 <- pheno_data[grep("Dayone", pheno_data$time),]
p2 <- pheno_data[grep("Daytwo", pheno_data$time),]

## We can test if the data were batched using the control probes:
raw <- ReadAffy(celfile.path = input_dir)
data <- affy::rma(raw)
norm_expr <- exprs(data)
colnames(norm_expr) <- str_extract(colnames(norm_expr), "GSM\\d+")
identical(rownames(pheno_data), colnames(norm_expr)) # [1] TRUE
source("~/Desktop/work_repo/github/misc/Affy_control_probe_check.R")
a <- polyA_check(norm_expr, pheno_data, "time")[[1]]
b <- hybrid_check(norm_expr, pheno_data, "time")[[1]]
c <- positive_check(norm_expr, pheno_data, "time")
pdf(file = paste(output_dir, "control_probes_batched.pdf", sep = ""), height = 7, width = 12)
plot_grid(a,b,c, nrow = 1)
dev.off()
## Poly A an hybridization show concerning batch effects, justifying our separation of the days.

# Analysis of day 1 data ----
celfiles <- list.files(input_dir, full.names = T)
c1 <- celfiles[which(str_extract(celfiles, "GSM\\d+") %in% p1$geo_accession)]
length(c1) # [1] 19
raw <- ReadAffy(filenames = c1)

## Raw data QC----
## According to the researchers, poor or abnormal arrays were already filtered out using the same methods
## we use (Mearburn et al. 2009). This is simply to double-check array quality.
library(affyPLM)
Pset <- fitPLM(raw, background = F, normalize = F) #input non-normalized data
pdf(file = paste(output_dir, "RLE_and_NUSE_day_1.pdf", sep = ""))
par(mfrow = c(1,2))
RLE(Pset,main = "RLE for Meaburn day 1 raw data")
NUSE(Pset, main = "NUSE for Meaburn day 1 raw data")
dev.off()
# These arrays have good statistics.

## Get expression calls using MAS 5.0 ----
calls <- mas5calls.AffyBatch(raw)
probe_pval <- assayData(calls)[["se.exprs"]] ## Returns p-values
dim(probe_pval)
# [1] 54675    19
minSamples <- min(colSums(table(p1$id, p1$time)))
print(minSamples)
# [1] 9
expressed <- rowSums(probe_pval < 0.05) >= minSamples

## Normalize the data into an eSet ----
data <- affy::rma(raw)

# Get expression estimates
unfiltered_normalized_expression <- exprs(data)
## Remove unexpressed probes
normalized_expression <- unfiltered_normalized_expression[expressed,]
dim(unfiltered_normalized_expression)
# [1] 54675    19
dim(normalized_expression)
# [1] 25419    19

normalized_expression[1:5,1:5]
#           GSM371373.CEL.gz GSM371374.CEL.gz GSM371376.CEL.gz GSM371378.CEL.gz GSM371383.CEL.gz
# 1007_s_at         7.376076         7.171558         7.129033         7.127949         6.941682
# 1053_at           7.195994         7.161228         7.194788         7.208823         7.268647
# 117_at            9.833055         9.674456         9.631341         9.478013         9.385483
# 121_at            7.799649         7.734322         7.524356         7.258983         7.480098
# 1255_g_at         2.731875         2.619617         2.555998         2.635680         2.636198

colnames(normalized_expression) <- gsub(pattern = "(.CEL.gz)", replacement = "", x = colnames(normalized_expression))
colnames(unfiltered_normalized_expression) <- gsub(pattern = "(.CEL.gz)", replacement = "", x = colnames(unfiltered_normalized_expression))
df <- data.frame(normalized_expression)
identical(colnames(df), rownames(p1))
# [1] TRUE
table(p1$sex, p1$id)
#        TD19901 TD19902 TD23461 TD23462 TD24111 TD24112 TD32681 TD32682 TD51281 TD51282
# Female       2       2       2       2       2       2       0       0       2       2
# Male         0       0       0       0       0       0       2       1       0       0
table(p1$time, p1$id)
#             TD19901 TD19902 TD23461 TD23462 TD24111 TD24112 TD32681 TD32682 TD51281 TD51282
# Dayone_10am       1       1       1       1       1       1       1       1       1       1
# Dayone_2pm        1       1       1       1       1       1       1       0       1       1


## Normalized data QC----
### Plots and quantiles
pdf(file = paste(output_dir, "features_hist_day1.pdf", sep = ""))
hist(normalized_expression)
dev.off()
pdf(file = paste(output_dir, "samples_box_day1.pdf", sep = ""))
boxplot(normalized_expression)
dev.off()
quantile(normalized_expression)
# 0%       25%       50%       75%      100% 
# 1.881568  5.579088  6.916074  8.236908 14.629243 


# What are the highly expressed RNAs?
avg_gene_exp <- rowMeans(x = normalized_expression, na.rm = T)
head(avg_gene_exp) # creates vector of each gene's mean
avg_gene_exp <- avg_gene_exp[order(avg_gene_exp, decreasing = T)] # order by decreasing order
mapIds(x = hgu133plus2.db, keys = head(names(avg_gene_exp)), keytype = "PROBEID", column = "SYMBOL")

# AFFX-hum_alu_at     217414_x_at     204018_x_at     214414_x_at     211745_x_at     209458_x_at 
#              NA          "HBA2"          "HBA1"          "HBA1"          "HBA1"          "HBA1" 

## MDS and PCA ----
### PCA tools ----
library(PCAtools)
p <- pca(df, metadata = p1, center = T, scale = T)
p_unfilt <- pca(unfiltered_normalized_expression, metadata = p1, center = T, scale = T)
pdf(file = paste(output_dir, "day1_unfiltered_pca.pdf", sep = ""))
biplot(p_unfilt, 
       showLoadings = F,
       lab = p1$id, 
       colby = "time", colkey = c("Dayone_10am" = "blue", "Dayone_2pm" = "red"), 
       legendPosition = 'right', legendLabSize = 12, legendIconSize = 3,
       title = 'Meaburn 4 hours',
       subtitle = 'day 1 - lowly epxressed genes not filtered')
dev.off()


scree <- screeplot(p)
ptime <- biplot(p, 
                showLoadings = F,
                lab = NULL, 
                colby = "time", colkey = c("Dayone_10am" = "blue", "Dayone_2pm" = "red"), 
                legendPosition = 'right', legendLabSize = 12, legendIconSize = 3,
                title = 'Meaburn 4 hours',
                subtitle = 'day 1 - lowly epxressed genes filtered')
psex <- biplot(p, 
               showLoadings = F,
               lab = NULL, 
               colby = "sex", 
               legendPosition = 'right', legendLabSize = 9, legendIconSize = 3)
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
pdf(file = paste(output_dir, "pca_day_one.pdf", sep = ""), height = 12, width = 16)
plot_grid(toprow, botrow, ncol = 1, align = 'v')
dev.off()


### MDS ----
qual_col_pals <- brewer.pal.info
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colors <- sample(col_vector, length(unique(p1$id)))
new_colors <- c()
for (i in 1:length(p1$id)) {
  new_colors[i] <- switch(p1$id[i],
                          "TD19901" = colors[1], "TD19902" = colors[2], "TD23461" = colors[3], "TD23462" = colors[4], 
                          "TD24111" = colors[5], "TD24112" = colors[6], "TD32681" = colors[7], "TD32682" = colors[8],
                          "TD51281" = colors[9], "TD51282" = colors[10]
  )
}
pdf(file = paste(output_dir, "mds_day_one.pdf", sep = ""))
limma::plotMDS(x = df, 
               cex = 1, 
               var.explained = T, 
               col = new_colors,
               labels = paste(p1$id, p1$sex, sep = ", "))
title("MDS Plot")
dev.off()


## limma analysis----

### Remove control probes ----
control_index <- grep("lys|Lys|thr|Thr|dap|Dap|phe|Phe|Trp|bioB|BioB|BioC|bioC|bioD|BioD|cre|Cre", x = row.names(df))
rownames(df)[control_index]
df <- df[-control_index,]
dim(df)
# [1] 26545    19
rownames(df)[grep("AFFX", x = rownames(df))]

### Fit LM ----
id <- factor(p1$id)
time <- factor(p1$time)
sex <- factor(p1$sex)

TS <- factor(paste(time,sex, sep = "."))
design <- model.matrix(~0+TS, df)
colnames(design) <- make.names(levels(TS))
cor <- duplicateCorrelation(df, design, block=id)
cor$consensus.correlation
# 0.393599

fit <- lmFit(object=df, design=design, block=id, correlation=cor$consensus.correlation)

#### Contrasts ----
my_contrasts <- makeContrasts(Time = (Dayone_2pm.Female + Dayone_2pm.Male)/2 -
                                (Dayone_10am.Female + Dayone_10am.Male)/2,
                              TimeMale = (Dayone_2pm.Male - Dayone_10am.Male) -
                                (Dayone_2pm.Female - Dayone_10am.Female),
                              levels = colnames(design))

fit2 <- contrasts.fit(fit, my_contrasts)
x <- eBayes(fit2)

## Collect results ----
y <- topTable(x, number = Inf, adjust.method = "fdr", sort.by = "F")
y$Symbol <- mapIds(x = hgu133plus2.db, keys = row.names(y), column = "SYMBOL", keytype = "PROBEID")

dt <- decideTests(x)
summary(dt)
#         Time TimeMale
# Down       0        0
# NotSig 26545    26545
# Up         0        0
dt <- decideTests(x, adjust.method = 'none')
summary(dt)
#         Time TimeMale
# Down     680      733
# NotSig 24352    24657
# Up      1513     1155

write.csv(y, file = paste(output_dir, "day1_limma_F.csv", sep = ""))

names(y)
tables <- list()
for (i in 1:2) {
  tables[[i]] <- topTable(x, coef = i, number = Inf, adjust.method = "fdr")
  names(tables)[i] <- names(y)[i]
  tables[[i]]$Symbol <- mapIds(x = hgu133plus2.db, keys = row.names(tables[[i]]), column = "SYMBOL", keytype = "PROBEID")
  write.csv(tables[[i]], file = paste(output_dir, "day1_", names(y)[i], "_limma.csv", sep = ""), row.names = T)
}

## Plot p-value histograms ----
for (i in 1:2) {
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
  pdf(file = paste(output_dir, names(y)[i], "day_1_p_value_histograms.pdf", sep = ""))
  plot_grid(a, b, labels = c('A. P-value histogram', 'B. BH-adjusted'))
  dev.off()
}


### Enhanced Volcano ----
a <- EnhancedVolcano(toptable = tables[[1]], 
                     x = "logFC",
                     y = "P.Value",
                     lab = tables[[1]]$Symbol,
                     title = "4-hour comparison",
                     subtitle ="Time main effect",
                     caption = "Meaburn Day 1",
                     col = c("black","lightblue","royalblue","blue"),
                     legendPosition = 'top', legendLabSize = 6.5, legendIconSize = 3,
                     pCutoff = 1e-2,
                     FCcutoff = 1,
                     labSize = (4.0),
                     colAlpha = 1,
                     pointSize = 1.0,
                     drawConnectors = TRUE)
b <- EnhancedVolcano(toptable = tables[[2]], 
                     x = "logFC",
                     y = "P.Value",
                     lab = tables[[2]]$Symbol,
                     title = "4-hour comparison",
                     subtitle ="Time interaction (male - female)",
                     caption = "Meaburn Day 1",
                     col = c("black","lightblue","royalblue","blue"),
                     legendPosition = 'top', legendLabSize = 6.5, legendIconSize = 3,
                     pCutoff = 1e-2,
                     FCcutoff = 1,
                     labSize = (4.0),
                     colAlpha = 1,
                     pointSize = 1.0,
                     drawConnectors = TRUE)
pdf(file = paste(output_dir, "day_1_volcanos.pdf", sep = ""), width = 9)
plot_grid(a, b, nrow=1)
dev.off()

# Analysis of day 2 data ----
c2 <- celfiles[which(str_extract(celfiles, "GSM\\d{6}") %in% p2$geo_accession)]
length(c2) # [1] 17
raw <- ReadAffy(filenames = c2)

## Raw data QC----
## According to the researchers, poor or abnormal arrays were already filtered out using the same methods
## we use (Mearburn et al. 2009). This is simply to double-check array quality.
library(affyPLM)
Pset <- fitPLM(raw, background = F, normalize = F) #input non-normalized data
pdf(file = paste(output_dir, "RLE_and_NUSE_day_2.pdf", sep = ""))
par(mfrow = c(1,2))
RLE(Pset,main = "RLE for Meaburn day 2 raw data")
NUSE(Pset, main = "NUSE for Meaburn day 2 raw data")
dev.off()
# These arrays have good statistics.

## Get expression calls using MAS 5.0 ----
calls <- mas5calls.AffyBatch(raw)
probe_pval <- assayData(calls)[["se.exprs"]] ## Returns p-values
dim(probe_pval)
# [1] 54675    17
minSamples <- min(colSums(table(p2$id, p2$time)))
print(minSamples)
# [1] 7
expressed <- rowSums(probe_pval < 0.05) >= minSamples

## Normalize the data into an eSet ----
data <- affy::rma(raw)

# Get expression estimates
unfiltered_normalized_expression <- exprs(data)
normalized_expression <- unfiltered_normalized_expression[expressed,]
dim(normalized_expression)
# [1] 26474    17

normalized_expression[1:5,1:5]
#          GSM371375.CEL.gz GSM371377.CEL.gz GSM371379.CEL.gz GSM371380.CEL.gz GSM371381.CEL.gz
# 1007_s_at         7.408323         7.290897         7.013294         6.694431         6.882215
# 1053_at           7.028865         7.160074         7.468729         7.249186         7.172744
# 117_at            9.305568         9.376824         9.305617        10.009949        10.174146
# 121_at            7.689531         7.639373         7.268682         7.098646         7.354284
# 1294_at           8.158172         8.307702         8.385276         8.202762         8.054674

colnames(normalized_expression) <- gsub(pattern = "(.CEL.gz)", replacement = "", x = colnames(normalized_expression))
colnames(unfiltered_normalized_expression) <- gsub(pattern = "(.CEL.gz)", replacement = "", x = colnames(unfiltered_normalized_expression))
df <- data.frame(normalized_expression)
identical(colnames(df), rownames(p2))
# [1] TRUE
table(p2$sex, p2$id)
#        TD19901 TD19902 TD23461 TD23462 TD24111 TD24112 TD32681 TD32682 TD51281 TD51282
# Female       1       2       2       1       2       2       0       0       2       2
# Male         0       0       0       0       0       0       2       1       0       0
table(p2$time, p2$id)
#             TD19901 TD19902 TD23461 TD23462 TD24111 TD24112 TD32681 TD32682 TD51281 TD51282
# Daytwo_10am       0       1       1       0       1       1       1       0       1       1
# Daytwo_2pm        1       1       1       1       1       1       1       1       1       1


## Normalized data QC----
### PCA tools ----
library(PCAtools)
p <- pca(df, metadata = p2, center = T, scale = T)
p_unfilt <- pca(unfiltered_normalized_expression, metadata = p2, center = T, scale = T)
pdf(file = paste(output_dir, "day2_unfiltered_pca.pdf", sep = ""))
biplot(p_unfilt, 
       showLoadings = F,
       lab = p2$id, 
       colby = "time", colkey = c("Daytwo_10am" = "blue", "Daytwo_2pm" = "red"), 
       legendPosition = 'right', legendLabSize = 12, legendIconSize = 3,
       title = 'Meaburn 4 hours',
       subtitle = 'day 2 - lowly epxressed genes not filtered')
dev.off()
# TD32681 and TD32682, 2pm samples are separated from most other samples on PCA.
p2[grep("TD32681|TD32682", p2$id),]
#                               title geo_accession  sex        time      id
# GSM371396 Blood_TD32681_Daytwo_10am     GSM371396 Male Daytwo_10am TD32681
# GSM371398  Blood_TD32681_Daytwo_2pm     GSM371398 Male  Daytwo_2pm TD32681
# GSM371400  Blood_TD32682_Daytwo_2pm     GSM371400 Male  Daytwo_2pm TD32682
grep("GSM371398|GSM371400",colnames(Pset@chip.coefs))
# [1] 12 13
pdf(paste(output_dir, "GSM371398_GSM371400_pseudo.pdf", sep = ""), width = 12, height = 9)
par(mfrow = c(1,2))
image(Pset, which = 12)
image(Pset, which = 13)
dev.off()


scree <- screeplot(p)
ptime <- biplot(p, 
                showLoadings = F,
                lab = NULL, 
                colby = "time", colkey = c("Daytwo_10am" = "blue", "Daytwo_2pm" = "red"), 
                legendPosition = 'right', legendLabSize = 12, legendIconSize = 3,
                title = 'Meaburn 4 hours',
                subtitle = 'day 2 - lowly expressed probes filtered')
psex <- biplot(p, 
               showLoadings = F,
               lab = NULL, 
               colby = "sex", 
               legendPosition = 'right', legendLabSize = 9, legendIconSize = 3)
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
pdf(file = paste(output_dir, "pca_day_two.pdf", sep = ""), height = 12, width = 16)
plot_grid(toprow, botrow, ncol = 1, align = 'v')
dev.off()


### MDS ----
qual_col_pals <- brewer.pal.info
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colors <- sample(col_vector, length(unique(p2$id)))
new_colors <- c()
for (i in 1:length(p2$id)) {
  new_colors[i] <- switch(p1$id[i],
                          "TD19901" = colors[1], "TD19902" = colors[2], "TD23461" = colors[3], "TD23462" = colors[4], 
                          "TD24111" = colors[5], "TD24112" = colors[6], "TD32681" = colors[7], "TD32682" = colors[8],
                          "TD51281" = colors[9], "TD51282" = colors[10]
  )
}
pdf(file = paste(output_dir, "mds_day_two.pdf", sep = ""))
limma::plotMDS(x = df, 
               cex = 1, 
               var.explained = T, 
               col = new_colors,
               labels = paste(p2$id, p2$sex, sep = ", "))
title("MDS Plot")
dev.off()

### Plots and quantiles
pdf(file = paste(output_dir, "features_hist_day2.pdf", sep = ""))
hist(normalized_expression)
dev.off()
pdf(file = paste(output_dir, "samples_box_day2.pdf", sep = ""))
boxplot(normalized_expression)
dev.off()
quantile(normalized_expression)
# 0%       25%       50%       75%      100% 
# 1.870306  5.594583  6.951321  8.279332 14.721409 

# What are the highly expressed RNAs?
avg_gene_exp <- rowMeans(x = normalized_expression, na.rm = T)
head(avg_gene_exp) # creates vector of each gene's mean
avg_gene_exp <- avg_gene_exp[order(avg_gene_exp, decreasing = T)] # order by decreasing order
mapIds(x = hgu133plus2.db, keys = head(names(avg_gene_exp)), keytype = "PROBEID", column = "SYMBOL")

# AFFX-hum_alu_at     204892_x_at     202917_s_at     206559_x_at     208834_x_at     208825_x_at 
#              NA        "EEF1A1"        "S100A8"        "EEF1A1"        "RPL23A"        "RPL23A" 


## limma analysis----

### Remove control probes ----
control_index <- grep("lys|Lys|thr|Thr|dap|Dap|phe|Phe|Trp|bioB|BioB|BioC|bioC|bioD|BioD|cre|Cre", x = row.names(df))
rownames(df)[control_index]
df <- df[-control_index,]
dim(df)
# [1] 26432    17
rownames(df)[grep("AFFX", x = rownames(df))]

### Fit LM ----
id <- factor(p2$id)
time <- factor(p2$time)
sex <- factor(p2$sex)

TS <- factor(paste(time,sex, sep = "."))
design <- model.matrix(~0+TS, df)
colnames(design) <- make.names(levels(TS))
cor <- duplicateCorrelation(df, design, block=id)
cor$consensus.correlation
# 0.5628484

fit <- lmFit(object=df, design=design, block=id, correlation=cor$consensus.correlation)

#### Contrasts ----
my_contrasts <- makeContrasts(Time = (Daytwo_2pm.Female + Daytwo_2pm.Male)/2 -
                                (Daytwo_10am.Female + Daytwo_10am.Male)/2,
                              TimeMale = (Daytwo_2pm.Male - Daytwo_10am.Male) -
                                (Daytwo_2pm.Female - Daytwo_10am.Female),
                              levels = colnames(design))

fit2 <- contrasts.fit(fit, my_contrasts)
x <- eBayes(fit2)

## Collect results ----
y <- topTable(x, number = Inf, adjust.method = "fdr", sort.by = "F")
y$Symbol <- mapIds(x = hgu133plus2.db, keys = row.names(y), column = "SYMBOL", keytype = "PROBEID")

dt <- decideTests(x)
summary(dt)
#         Time TimeMale
# Down     107       97
# NotSig 25766    26014
# Up       559      321
dt <- decideTests(x, adjust.method = 'none')
summary(dt)
#         Time TimeMale
# Down    1127     1016
# NotSig 21743    22149
# Up      3562     3267

write.csv(y, file = paste(output_dir, "day2_limma_F.csv", sep = ""))

names(y)
tables <- list()
for (i in 1:2) {
  tables[[i]] <- topTable(x, coef = i, number = Inf, adjust.method = "fdr")
  names(tables)[i] <- names(y)[i]
  tables[[i]]$Symbol <- mapIds(x = hgu133plus2.db, keys = row.names(tables[[i]]), column = "SYMBOL", keytype = "PROBEID")
  write.csv(tables[[i]], file = paste(output_dir, "day2_", names(y)[i], "_limma.csv", sep = ""), row.names = T)
}

## Plot p-value histograms ----
for (i in 1:2) {
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
  pdf(file = paste(output_dir, names(y)[i], "day_2_p_value_histograms.pdf", sep = ""))
  plot_grid(a, b, labels = c('A. P-value histogram', 'B. BH-adjusted'))
  dev.off()
}


### Enhanced Volcano ----
a <- EnhancedVolcano(toptable = tables[[1]], 
                     x = "logFC",
                     y = "P.Value",
                     lab = tables[[1]]$Symbol,
                     title = "4-hour comparison",
                     subtitle ="Time main effect",
                     caption = "Meaburn Day 2",
                     col = c("black","lightblue","royalblue","blue"),
                     legendPosition = 'top', legendLabSize = 6.5, legendIconSize = 3,
                     pCutoff = 1e-2,
                     FCcutoff = 1,
                     labSize = (4.0),
                     colAlpha = 1,
                     pointSize = 1.0,
                     drawConnectors = TRUE)
b <- EnhancedVolcano(toptable = tables[[2]], 
                     x = "logFC",
                     y = "P.Value",
                     lab = tables[[2]]$Symbol,
                     title = "4-hour comparison",
                     subtitle ="Time interaction (male - female)",
                     caption = "Meaburn Day 2",
                     col = c("black","lightblue","royalblue","blue"),
                     legendPosition = 'top', legendLabSize = 6.5, legendIconSize = 3,
                     pCutoff = 1e-2,
                     FCcutoff = 1,
                     labSize = (4.0),
                     colAlpha = 1,
                     pointSize = 1.0,
                     drawConnectors = TRUE)
pdf(file = paste(output_dir, "day_2_volcanos.pdf", sep = ""), width = 9)
plot_grid(a, b, nrow=1)
dev.off()