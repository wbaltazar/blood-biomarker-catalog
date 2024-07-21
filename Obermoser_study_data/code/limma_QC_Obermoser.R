### SUMMARY: Do vaccine injections induce changes in gene expression?
# Data from multiple healthy patients sampled before, baseline, and after vaccination (w/ saline controls)
# We include the type of vaccination in this analysis to potentially explain the variability

# The data was collected in three cohorts from either the vein or the finger. Additionally, the data collected
# across cohorts measures different time intervals, including an hour-by-hour measurement.

# Load libraries ----

library(AnnotationDbi)
library(beadarray)
library(illuminaio)
library(illuminaHumanv3.db)
library(lumi)
library(limma)
library(oligo)
library(tidyverse)
library(cowplot)
library(PCAtools)
library(GEOquery)
library(RColorBrewer)
library(EnhancedVolcano)

## INPUT: Raw file folder from GSE30101
input_dir <- "~/Desktop/work_repo/data/GSE30101_RAW/"
## OUTPUT: QC plots, limma tables, volcano plots and p-value histograms
output_dir <- "~/Desktop/work_repo/github/Obermoser_study_data/output/"

# Load in metadata ----
gset <- getGEO("GSE30101", AnnotGPL = TRUE, GSEMatrix = TRUE)
pheno_data_org <- pData(gset$GSE30101_series_matrix.txt.gz)
feature_data <- fData(gset$GSE30101_series_matrix.txt.gz)

colnames(pheno_data_org)
pheno_data_org$data_processing[1]
pheno_data_org$description.1
pheno_data <- pheno_data_org %>% 
  dplyr::select(title, 
                geo_accession, 
                subject = `subject id:ch1`,
                age = `age:ch1`, 
                time = `day:ch1`,
                treat = `vaccine:ch1`,
                sex = `gender:ch1`,
                race = `race:ch1`, 
                ethnicity = `ethnicity:ch1`,
                array = description.1,
                sample = source_name_ch1,
                group = `sample set:ch1`)
dim(pheno_data)
# [1] 693  12
pheno_data <- pheno_data[pheno_data$sample == "Whole Blood",]
dim(pheno_data)
# [1] 621  12

## Divide the data into appropriate pools ----
## KEY: p = pheno_data, 1,2, or 3 = cohort, f,v = finger or vein, as specified in GEO
## We do not use cohort 3: no saline controls were included.
unique(pheno_data$group)
# [1] "Training_Set_Vein"   "Test_Set_Vein"       "Validation_Vein"     "Training_Set_Finger" "Test_Set_Finger"  
p1v <- pheno_data[grep("Training_Set_Vein", pheno_data$group),]
p1f <- pheno_data[grep("Test_Set_Finger", pheno_data$group),]
p2v <- pheno_data[grep("Test_Set_Vein", pheno_data$group),]
p2f <- pheno_data[grep("Training_Set_Finger", pheno_data$group),]
## Note that p1f actually belongs to cohort 2:
p1f_patients <- paste(p1f$subject, p1f$age, p1f$sex, p1f$race, sep = "-") %>% unique()
p2f_patients <- paste(p2f$subject, p2f$age, p2f$sex, p2f$race, sep = "-") %>% unique()
p2v_patients <- paste(p2v$subject, p2v$age, p2v$sex, p2v$race, sep = "-") %>% unique()
p1v_patients <- paste(p1v$subject, p1v$age, p1v$sex, p1v$race, sep = "-") %>% unique()
all(p2f_patients %in% p2v_patients) # [1] TRUE
all(p2f_patients %in% p1f_patients) # [1] TRUE
all(p1v_patients %in% p1f_patients) # [1] FALSE

## Open probe annotation file ----
BGX <- readBGX(file = paste(input_dir,"GPL6947_HumanHT-12_V3_0_R1_11283641_A.bgx.gz", sep = ""))
str(BGX)
## Raw data file names ----
# View cohort clustering: justifies the separate analyses of batches ----
files <- list.files(input_dir, pattern = "GSM.*\\.txt.gz$", full.names = TRUE)
length(files) # [1] 621
files <- files[str_extract(files, "GSM......") %in% pheno_data$geo_accession]
length(files) # [1] 621
x <- lapply(
  files,
  function(x) read.delim(file = x))
x <- do.call(cbind, x)
x[1:10,1:10]
dim(x)
# [1] 48803  1863
rownames(x) <- x[,1]
x <- as.data.frame(x)
detectionpvalues <- x[,grep("Detection",colnames(x))]
x <- x[,grep("AVG",colnames(x))]
x[1:5,1:5]

data <- new('EListRaw')
data@.Data[[1]] <- 'illumina'
# data@.Data[[2]] <- targetinfo
data@.Data[[3]] <- BGX
data@.Data[[4]] <- x
data@.Data[[5]] <- NULL
data$E <- x
data$genes <- BGX
data$other$Detection <- detectionpvalues
dim(data$E)
# [1] 48803   621
data_norm <- neqc(data)
controls <- data_norm$genes$controls$Probe_Id
expressed <- rowSums(data_norm$other$Detection < 0.05) >= 3
data_norm$E <- data_norm$E[expressed,]
data_norm$E <- data_norm$E[-which(rownames(data_norm$E) %in% controls),]
dim(data_norm$E)

## Prinicpal Component for batches ----
library(PCAtools)
identical(substr(pheno_data$array, start = 1, stop = 10), str_extract(colnames(data$E),"\\d........."))
colnames(data_norm$E) <- rownames(pheno_data)
# Organize the pheno_data in an interpretable way
pheno_data <- pheno_data %>% 
  mutate(`Blood source` = str_extract(pheno_data$group, "Vein|Finger")) %>% 
  mutate(`Cohort` = factor(ifelse(pheno_data$group %in% c("Training_Set_Vein","Test_Set_Finger"), 1, ifelse(pheno_data$group %in% c("Test_Set_Vein","Training_Set_Finger"), 2, 3)))) %>%
  mutate(`Treat2` = ifelse(pheno_data$treat %in% c("Flu","FLUZONE", "Influenza"), "Fluzone", ifelse(pheno_data$treat %in% c("PNEUM","Pneumovax"), "Pneumovax", "Saline"))) %>% 
  mutate(`CohortTreat` = paste(`Cohort`, `Treat2`, sep = ', '), `CohortBloodTreat` = paste(`Cohort`, `Blood source`, `Treat2`, sep = ', '))
p <- pca(data_norm$E, metadata = pheno_data, center = T, scale = T)
pdf(file = paste(output_dir, "batches_pca.pdf"))
biplot(p, 
       colby = "group", 
       legendPosition = "right",
       lab = "", 
       encircle = T,
       title = "Obermoser: normalized data files colored by cohort labels", titleLabSize = 14,
       subtitle = "Note: in final analysis, cohort 3 is excluded. Additionally, Training_Set_Vein samples with 'Saline' are removed.",
       subtitleLabSize = 7, legendLabSize = 8)
dev.off()
pdf(file = paste(output_dir, "batches_samplesite_pca.pdf"))
biplot(p, 
       colby = "Blood source", 
       legendPosition = "right",
       lab = "", 
       encircle = T,
       title = "Obermoser: normalized data files colored by blood source", titleLabSize = 14,
       subtitleLabSize = 7)
dev.off()
pdf(file = paste(output_dir, "batches_treatments_pca.pdf"))
biplot(p, colby = "group", 
       legendPosition = "right", 
       lab = "", shape = "treat", 
       shapekey = c("saline" = 8, "Saline" = 9, "NS" = 10),
       title = "Obermoser normalized samples, saline",
       subtitle = "Bottom left shows Saline. It's unusual that these samples cluster near pneumovax samples.",
       subtitleLabSize = 7)
dev.off()

# Analysis 1: Training_Set_Vein (p1v) ----
## This seems to be the data that the researchers primarily analyze in the Obermoser paper.
## Reading in illumina data: https://www.biostars.org/p/403439/#403446

## Step 1: Read in raw expression data ----
p1v[1:5,1:5]
#                               title geo_accession subject age time
# GSM744835 BAYVAC002T0 [4835507041F]     GSM744835       2  40   -7
# GSM744836 BAYVAC002T2 [4835507052A]     GSM744836       2  40    0
# GSM744837 BAYVAC002T3 [4835507049K]     GSM744837       2  40    1
# GSM744838 BAYVAC002T4 [4853887047C]     GSM744838       2  40    3
# GSM744839 BAYVAC002T5 [4835507042G]     GSM744839       2  40    7
levels(factor(p1v$time))
# [1] "-7" "0"  "1"  "10" "14" "21" "28" "3"  "7" 
## Note that time units here are in days, and 0 represents the day vaccination was administered.
dim(p1v)
# [1] 214  12

temp <- files[str_extract(files, "GSM\\d+") %in% p1v$geo_accession]
x <- lapply(
  temp,
  function(x) read.delim(file = x))
x <- do.call(cbind, x)
x[1:5,1:5]
#              X X4835507041_F.AVG_Signal X4835507041_F.Detection_Pval          X.1 X4835507052_A.AVG_Signal
# 1 ILMN_1885639                   2.5759                       0.3531 ILMN_1885639                  10.3486
# 2 ILMN_2209417                 -10.2824                       0.8748 ILMN_2209417                  -7.2999
# 3 ILMN_1765401                  55.5706                       0.0013 ILMN_1765401                  38.5928
# 4 ILMN_1652171                  11.9083                       0.0856 ILMN_1652171                  -1.8236
# 5 ILMN_1652170                  16.5814                       0.0356 ILMN_1652170                   5.1670
dim(x)
# [1] 48803   642
rownames(x) <- x[,1]
x <- as.data.frame(x)
detectionpvalues <- x[,grep("Detection",colnames(x))]
x <- x[,grep("AVG",colnames(x))]
x[1:2,1:5]
# X4835507041_F.AVG_Signal X4835507052_A.AVG_Signal X4835507049_K.AVG_Signal X4853887047_C.AVG_Signal
# ILMN_1885639                   2.5759                  10.3486                  -0.4804                  -3.0539
# ILMN_2209417                 -10.2824                  -7.2999                 -11.3488                   4.5483
# X4835507042_G.AVG_Signal
# ILMN_1885639                   6.9538
# ILMN_2209417                  -1.2148

### Create an EList object ----
data <- new('EListRaw')
data@.Data[[1]] <- 'illumina'
data@.Data[[3]] <- BGX
data@.Data[[4]] <- x
data@.Data[[5]] <- NULL
data$E <- x
data$genes <- BGX
data$other$Detection <- detectionpvalues

### Unsupervised clustering analysis: determine samples for analysis ----
### This cohort is interesting because there appear to be some sort of replicates in the analysis. However,
### it has been difficult to tell which of these sample groups represent the actual samples used in analysis.
### We use unsupervised clustering to justify the usage of 'saline' samples over 'Saline' samples.
identical(substr(p1v$array, start = 1, stop = 10),str_extract(colnames(data$E),"\\d........."))
names(data$E) <- rownames(p1v)
p <- pca(data$E, metadata = p1v, center = T, scale = T)
pdf(file = paste(output_dir, "training_vein_pre_QC_pca.pdf"))
biplot(p, colby = "treat", encircle = T, legendPosition = "right",
       title = "Clustering all samples", subtitle = "Saline, rather than saline, cluster with Pneumovax.")
dev.off()

### Hierarchical clustering
corp <- data$E[,which(p1v$treat %in% c("saline", "Saline"))]
colnames(corp) <- paste(p1v[which(p1v$treat %in% c("saline", "Saline")), "subject"], p1v[which(p1v$treat %in% c("saline", "Saline")), "treat"])
library(pheatmap)
pdf(file = paste(output_dir, "training_vein_pre_QC_hierarchical_clusters.pdf"))
pheatmap(cor(corp), cluster_cols = F, cluster_rows = F)
dev.off()
# This makes it clear that the saline samples are more correlated with one another than the Saline samples are.
# This most fits our expectations for a saline treatment, so we use it for analysis.

## Repeat step 1 w/o unknown samples ----
p1v <- p1v[-grep("Saline",p1v$treat), ]
temp <- files[str_extract(files, "GSM\\d+") %in% p1v$geo_accession]
x <- lapply(
  temp,
  function(x) read.delim(file = x))
x <- do.call(cbind, x)
dim(x)
# [1] 48803   483
rownames(x) <- x[,1]
x <- as.data.frame(x)
detectionpvalues <- x[,grep("Detection",colnames(x))]
x <- x[,grep("AVG",colnames(x))]
data <- new('EListRaw')
data@.Data[[1]] <- 'illumina'
data@.Data[[3]] <- BGX
data@.Data[[4]] <- x
data@.Data[[5]] <- NULL
data$E <- x
data$genes <- BGX
data$other$Detection <- detectionpvalues

### Background correction and normalization ----
data_norm <- neqc(data)
table(data_norm$genes$controls$Reporter_Group_Name)
# biotin            cy3_hyb       housekeeping           labeling low_stringency_hyb 
# 2                  6                  7                  2                  8 
# negative 
# 759 


## Remove unexpressed probes and controls
controls <- data_norm$genes$controls$Probe_Id
expressed <- rowSums(data_norm$other$Detection < 0.05) >= 3
data_norm$E <- data_norm$E[expressed,]
data_norm$E <- data_norm$E[-which(rownames(data_norm$E) %in% controls),]
dim(data_norm$E)
# [1] 30998   161

## Step 2: PCA and MDS QC analysis ----
identical(substr(p1v$array, start = 1, stop = 10), str_extract(colnames(data_norm$E),"\\d........."))
# [1] TRUE
colnames(data_norm$E) <- rownames(p1v)
### PCAtools ----
p <- pca(data_norm$E, metadata = p1v, center = T, scale = T)

scree <- screeplot(p, getComponents(p, c(1:10)))
ptime <- biplot(p, 
                showLoadings = F,
                lab = NULL, 
                colby = "time", 
                legendPosition = 'right', legendLabSize = 12, legendIconSize = 3,
                title = 'Obermoser QC',
                subtitle = 'p1v samples')
psex <- biplot(p, 
               showLoadings = F,
               lab = NULL, 
               colby = "sex", 
               legendPosition = 'right', legendLabSize = 9, legendIconSize = 3)
psub <- biplot(p, 
               showLoadings = F,
               lab = NULL, 
               colby = "subject",
               legendPosition = 'right', legendLabSize = 9, legendIconSize = 3)
ptre <-  biplot(p, 
                showLoadings = F,
                lab = NULL, 
                colby = "treat",
                legendPosition = 'right', legendLabSize = 9, legendIconSize = 3)
peigencor <- eigencorplot(p,
                          components = getComponents(p, 1:10),
                          metavars = c("time", 'sex', 'race', 'ethnicity', 'treat'),
                          cexCorval = 0.5,
                          fontCorval = 1, 
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
botrow <- plot_grid(ptre, scree, peigencor, align = 'h', nrow = 1)
pdf(file = paste(output_dir, "pca_p1v_QC.pdf", sep = ""), height = 12, width = 16)
plot_grid(toprow, botrow, ncol = 1, align = 'v')
dev.off()

### MDS ----
### MDS can be good for getting really nice separation and pulling out outliers.
qual_col_pals <- brewer.pal.info
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colors <- sample(col_vector, length(unique(p1v$subject)))
new_colors <- c()
for (i in 1:length(p1v$subject)) {
  new_colors[i] <- switch(p1v$subject[i],
                          "2" = colors[1], "9" = colors[2], "10" = colors[3], "12" = colors[4], 
                          "15" = colors[5], "16" = colors[6], "20" = colors[7], "21" = colors[8],
                          "22" = colors[9], "23" = colors[10], "24" = colors[11], "25" = colors[12],
                          "1" = colors[13], "3" = colors[14], "5" = colors[15], "6" = colors[16],
                          "13" = colors[17], "18" = colors[18]
  )
}
pdf(file = paste(output_dir, "mds_p1v_qc.pdf", sep = ""))
limma::plotMDS(x = data_norm$E, 
               cex = 1, 
               var.explained = T, 
               col = new_colors,
               labels = paste(p1v$subject, p1v$sex, p1v$treat, sep = ", "))
title("MDS Plot- training set vein")
dev.off()

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
## Step 4: limma analysis ----
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

id <- factor(p1v$subject)
time <- factor(p1v$time, levels = c("-7", "0","1", "3","7","10","14","21","28"), ordered = T)
sex <- factor(p1v$sex)
vax <- factor(p1v$treat)
eth <- factor(p1v$ethnicity)
race <- factor(p1v$race)

df <- data.frame(data_norm$E)

TS <- factor(paste(sex,time,vax,race,eth, sep = "."))
design <- model.matrix(~0+TS, df)
colnames(design) <- make.names(levels(TS))
cor <- duplicateCorrelation(df, design, block=id)
cor$consensus.correlation
# [1] 0.1135356

fit <- lmFit(object=df, design=design, block=id, correlation=cor$consensus.correlation)

# contrasts
levels(time)
# [1] "-7" "0"  "1"  "3"  "7"  "10" "14" "21" "28"
my_contrasts <- makeContrasts(# Time contrasts
                              BaseVSzeroALL = (F.0.Flu.African.American.NA + F.0.Flu.Caucasian.Not.Hispanic.or.Latino + F.0.Pneumovax.Asian.NA + F.0.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.0.saline.Caucasian.Hispanic + F.0.saline.Caucasian.Not.Hispanic.or.Latino + M.0.Flu.Asian.NA + M.0.Flu.Caucasian.Not.Hispanic.or.Latino + M.0.Pneumovax.Asian.NA + M.0.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.0.saline.Asian.NA + M.0.saline.Caucasian.Not.Hispanic.or.Latino)/12-
                                (F..7.Flu.African.American.NA + F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Asian.NA + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Asian.NA + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Asian.NA + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Asian.NA + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/12,
                              BaseVSoneALL = (F.1.Flu.African.American.NA + F.1.Flu.Caucasian.Not.Hispanic.or.Latino + F.1.Pneumovax.Asian.NA + F.1.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.1.saline.Caucasian.Hispanic + F.1.saline.Caucasian.Not.Hispanic.or.Latino + M.1.Flu.Asian.NA + M.1.Flu.Caucasian.Not.Hispanic.or.Latino + M.1.Pneumovax.Asian.NA + M.1.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.1.saline.Asian.NA + M.1.saline.Caucasian.Not.Hispanic.or.Latino)/12-
                                (F..7.Flu.African.American.NA + F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Asian.NA + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Asian.NA + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Asian.NA + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Asian.NA + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/12,
                              BaseVSthreeALL = (F.3.Flu.African.American.NA + F.3.Flu.Caucasian.Not.Hispanic.or.Latino + F.3.Pneumovax.Asian.NA + F.3.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.3.saline.Caucasian.Hispanic + F.3.saline.Caucasian.Not.Hispanic.or.Latino + M.3.Flu.Asian.NA + M.3.Flu.Caucasian.Not.Hispanic.or.Latino + M.3.Pneumovax.Asian.NA + M.3.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.3.saline.Asian.NA + M.3.saline.Caucasian.Not.Hispanic.or.Latino)/12-
                                (F..7.Flu.African.American.NA + F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Asian.NA + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Asian.NA + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Asian.NA + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Asian.NA + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/12,
                              BaseVSsevenALL = (F.7.Flu.African.American.NA + F.7.Flu.Caucasian.Not.Hispanic.or.Latino + F.7.Pneumovax.Asian.NA + F.7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.7.saline.Caucasian.Hispanic + F.7.saline.Caucasian.Not.Hispanic.or.Latino + M.7.Flu.Asian.NA + M.7.Flu.Caucasian.Not.Hispanic.or.Latino + M.7.Pneumovax.Asian.NA + M.7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.7.saline.Asian.NA + M.7.saline.Caucasian.Not.Hispanic.or.Latino)/12-
                                (F..7.Flu.African.American.NA + F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Asian.NA + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Asian.NA + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Asian.NA + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Asian.NA + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/12,
                              BaseVStenALL = (F.10.Flu.African.American.NA + F.10.Flu.Caucasian.Not.Hispanic.or.Latino + F.10.Pneumovax.Asian.NA + F.10.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.10.saline.Caucasian.Not.Hispanic.or.Latino + M.10.Flu.Asian.NA + M.10.Flu.Caucasian.Not.Hispanic.or.Latino + M.10.Pneumovax.Asian.NA + M.10.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.10.saline.Asian.NA + M.10.saline.Caucasian.Not.Hispanic.or.Latino)/11-
                                (F..7.Flu.African.American.NA + F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Asian.NA + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Asian.NA + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Asian.NA + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Asian.NA + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/12,
                              BaseVSfourteenALL = (F.14.Flu.African.American.NA + F.14.Flu.Caucasian.Not.Hispanic.or.Latino + F.14.Pneumovax.Asian.NA + F.14.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.14.saline.Caucasian.Hispanic + F.14.saline.Caucasian.Not.Hispanic.or.Latino + M.14.Flu.Asian.NA + M.14.Flu.Caucasian.Not.Hispanic.or.Latino + M.14.Pneumovax.Asian.NA + M.14.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.14.saline.Asian.NA + M.14.saline.Caucasian.Not.Hispanic.or.Latino)/12-
                                (F..7.Flu.African.American.NA + F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Asian.NA + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Asian.NA + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Asian.NA + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Asian.NA + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/12,
                              BaseVStwentyOneALL = (F.21.Flu.African.American.NA + F.21.Flu.Caucasian.Not.Hispanic.or.Latino + F.21.Pneumovax.Asian.NA + F.21.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.21.saline.Caucasian.Hispanic + F.21.saline.Caucasian.Not.Hispanic.or.Latino + M.21.Flu.Asian.NA + M.21.Flu.Caucasian.Not.Hispanic.or.Latino + M.21.Pneumovax.Asian.NA + M.21.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.21.saline.Asian.NA + M.21.saline.Caucasian.Not.Hispanic.or.Latino)/12-
                                (F..7.Flu.African.American.NA + F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Asian.NA + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Asian.NA + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Asian.NA + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Asian.NA + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/12,
                              BaseVStwentyEightALL = (F.28.Flu.African.American.NA + F.28.Flu.Caucasian.Not.Hispanic.or.Latino + F.28.Pneumovax.Asian.NA + F.28.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.28.saline.Caucasian.Hispanic + F.28.saline.Caucasian.Not.Hispanic.or.Latino + M.28.Flu.Asian.NA + M.28.Flu.Caucasian.Not.Hispanic.or.Latino + M.28.Pneumovax.Asian.NA + M.28.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.28.saline.Asian.NA + M.28.saline.Caucasian.Not.Hispanic.or.Latino)/12-
                                (F..7.Flu.African.American.NA + F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Asian.NA + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Asian.NA + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Asian.NA + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Asian.NA + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/12,
                              # Saline
                              SalineBaseVSzero = (F.0.saline.Caucasian.Hispanic + F.0.saline.Caucasian.Not.Hispanic.or.Latino + M.0.saline.Asian.NA + M.0.saline.Caucasian.Not.Hispanic.or.Latino)/4-
  (F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Asian.NA + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/4,
                              SalineBaseVSone = (F.1.saline.Caucasian.Hispanic + F.1.saline.Caucasian.Not.Hispanic.or.Latino + M.1.saline.Asian.NA + M.1.saline.Caucasian.Not.Hispanic.or.Latino)/4-
  (F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Asian.NA + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/4,
                              SalineBaseVSthree = (F.3.saline.Caucasian.Hispanic + F.3.saline.Caucasian.Not.Hispanic.or.Latino + M.3.saline.Asian.NA + M.3.saline.Caucasian.Not.Hispanic.or.Latino)/4-
  (F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Asian.NA + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/4,
                              SalineBaseVSseven = (F.7.saline.Caucasian.Hispanic + F.7.saline.Caucasian.Not.Hispanic.or.Latino + M.7.saline.Asian.NA + M.7.saline.Caucasian.Not.Hispanic.or.Latino)/4-
  (F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Asian.NA + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/4,
                              SalineBaseVSten = (F.10.saline.Caucasian.Not.Hispanic.or.Latino + M.10.saline.Asian.NA + M.10.saline.Caucasian.Not.Hispanic.or.Latino)/3-
  (F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Asian.NA + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/4,
                              SalineBaseVSfourteen = (F.14.saline.Caucasian.Hispanic + F.14.saline.Caucasian.Not.Hispanic.or.Latino + M.14.saline.Asian.NA + M.14.saline.Caucasian.Not.Hispanic.or.Latino)/4-
  (F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Asian.NA + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/4,
                              SalineBaseVStwentyOne = (F.21.saline.Caucasian.Hispanic + F.21.saline.Caucasian.Not.Hispanic.or.Latino + M.21.saline.Asian.NA + M.21.saline.Caucasian.Not.Hispanic.or.Latino)/4-
  (F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Asian.NA + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/4,
                              SalineBaseVStwentyEight = (F.28.saline.Caucasian.Hispanic + F.28.saline.Caucasian.Not.Hispanic.or.Latino + M.28.saline.Asian.NA + M.28.saline.Caucasian.Not.Hispanic.or.Latino)/4-
  (F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Asian.NA + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/4,
                              # Pneumovax
                              PneumBaseVSzero = (F.0.Pneumovax.Asian.NA + F.0.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.0.Pneumovax.Asian.NA + M.0.Pneumovax.Caucasian.Not.Hispanic.or.Latino)/4-
  (F..7.Pneumovax.Asian.NA + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Asian.NA + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino)/4,
                              PneumBaseVSone = (F.1.Pneumovax.Asian.NA + F.1.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.1.Pneumovax.Asian.NA + M.1.Pneumovax.Caucasian.Not.Hispanic.or.Latino)/4-
  (F..7.Pneumovax.Asian.NA + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Asian.NA + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino)/4,
                              PneumBaseVSthree = (F.3.Pneumovax.Asian.NA + F.3.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.3.Pneumovax.Asian.NA + M.3.Pneumovax.Caucasian.Not.Hispanic.or.Latino)/4-
  (F..7.Pneumovax.Asian.NA + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Asian.NA + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino)/4,
                              PneumBaseVSseven = (F.7.Pneumovax.Asian.NA + F.7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.7.Pneumovax.Asian.NA + M.7.Pneumovax.Caucasian.Not.Hispanic.or.Latino)/4-
  (F..7.Pneumovax.Asian.NA + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Asian.NA + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino)/4,
                              PneumBaseVSten = (F.10.Pneumovax.Asian.NA + F.10.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.10.Pneumovax.Asian.NA + M.10.Pneumovax.Caucasian.Not.Hispanic.or.Latino)/4-
  (F..7.Pneumovax.Asian.NA + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Asian.NA + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino)/4,
                              PneumBaseVSfourteen = (F.14.Pneumovax.Asian.NA + F.14.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.14.Pneumovax.Asian.NA + M.14.Pneumovax.Caucasian.Not.Hispanic.or.Latino)/4-
  (F..7.Pneumovax.Asian.NA + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Asian.NA + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino)/4,
                              PneumBaseVStwentyOne = (F.21.Pneumovax.Asian.NA + F.21.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.21.Pneumovax.Asian.NA + M.21.Pneumovax.Caucasian.Not.Hispanic.or.Latino)/4-
  (F..7.Pneumovax.Asian.NA + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Asian.NA + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino)/4,
                              PneumBaseVStwentyEight = (F.28.Pneumovax.Asian.NA + F.28.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.28.Pneumovax.Asian.NA + M.28.Pneumovax.Caucasian.Not.Hispanic.or.Latino)/4-
  (F..7.Pneumovax.Asian.NA + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Asian.NA + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino)/4,
                              # Flu
                              FluBaseVSzero = (F.0.Flu.African.American.NA + F.0.Flu.Caucasian.Not.Hispanic.or.Latino + M.0.Flu.Asian.NA + M.0.Flu.Caucasian.Not.Hispanic.or.Latino)/4-
  (F..7.Flu.African.American.NA + F..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Asian.NA + M..7.Flu.Caucasian.Not.Hispanic.or.Latino)/4,
                              FluBaseVSone = (F.1.Flu.African.American.NA + F.1.Flu.Caucasian.Not.Hispanic.or.Latino + M.1.Flu.Asian.NA + M.1.Flu.Caucasian.Not.Hispanic.or.Latino)/4-
  (F..7.Flu.African.American.NA + F..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Asian.NA + M..7.Flu.Caucasian.Not.Hispanic.or.Latino)/4,
                              FluBaseVSthree = (F.3.Flu.African.American.NA + F.3.Flu.Caucasian.Not.Hispanic.or.Latino + M.3.Flu.Asian.NA + M.3.Flu.Caucasian.Not.Hispanic.or.Latino)/4-
  (F..7.Flu.African.American.NA + F..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Asian.NA + M..7.Flu.Caucasian.Not.Hispanic.or.Latino)/4,
                              FluBaseVSseven = (F.7.Flu.African.American.NA + F.7.Flu.Caucasian.Not.Hispanic.or.Latino + M.7.Flu.Asian.NA + M.7.Flu.Caucasian.Not.Hispanic.or.Latino)/4-
  (F..7.Flu.African.American.NA + F..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Asian.NA + M..7.Flu.Caucasian.Not.Hispanic.or.Latino)/4,
                              FluBaseVSten = (F.10.Flu.African.American.NA + F.10.Flu.Caucasian.Not.Hispanic.or.Latino + M.10.Flu.Asian.NA + M.10.Flu.Caucasian.Not.Hispanic.or.Latino)/4-
  (F..7.Flu.African.American.NA + F..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Asian.NA + M..7.Flu.Caucasian.Not.Hispanic.or.Latino)/4,
                              FluBaseVSfourteen = (F.14.Flu.African.American.NA + F.14.Flu.Caucasian.Not.Hispanic.or.Latino + M.14.Flu.Asian.NA + M.14.Flu.Caucasian.Not.Hispanic.or.Latino)/4-
  (F..7.Flu.African.American.NA + F..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Asian.NA + M..7.Flu.Caucasian.Not.Hispanic.or.Latino)/4,
                              FluBaseVStwentyOne = (F.21.Flu.African.American.NA + F.21.Flu.Caucasian.Not.Hispanic.or.Latino + M.21.Flu.Asian.NA + M.21.Flu.Caucasian.Not.Hispanic.or.Latino)/4-
  (F..7.Flu.African.American.NA + F..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Asian.NA + M..7.Flu.Caucasian.Not.Hispanic.or.Latino)/4,
                              FluBaseVStwentyEight = (F.28.Flu.African.American.NA + F.28.Flu.Caucasian.Not.Hispanic.or.Latino + M.28.Flu.Asian.NA + M.28.Flu.Caucasian.Not.Hispanic.or.Latino)/4-
  (F..7.Flu.African.American.NA + F..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Asian.NA + M..7.Flu.Caucasian.Not.Hispanic.or.Latino)/4,
                              # Time Male - Female
                              IMaleBaseVSzero = ((M.0.Flu.Asian.NA + M.0.Flu.Caucasian.Not.Hispanic.or.Latino + M.0.Pneumovax.Asian.NA + M.0.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.0.saline.Asian.NA + M.0.saline.Caucasian.Not.Hispanic.or.Latino)/6- 
    (M..7.Flu.Asian.NA + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Asian.NA + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Asian.NA + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/6) -
    ((F.0.Flu.African.American.NA + F.0.Flu.Caucasian.Not.Hispanic.or.Latino + F.0.Pneumovax.Asian.NA + F.0.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.0.saline.Caucasian.Hispanic + F.0.saline.Caucasian.Not.Hispanic.or.Latino)/6-
       (F..7.Flu.African.American.NA + F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Asian.NA + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino)/6),
                              IMaleBaseVSone = ((M.1.Flu.Asian.NA + M.1.Flu.Caucasian.Not.Hispanic.or.Latino + M.1.Pneumovax.Asian.NA + M.1.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.1.saline.Asian.NA + M.1.saline.Caucasian.Not.Hispanic.or.Latino)/6-
    (M..7.Flu.Asian.NA + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Asian.NA + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Asian.NA + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/6) -
    ((F.1.Flu.African.American.NA + F.1.Flu.Caucasian.Not.Hispanic.or.Latino + F.1.Pneumovax.Asian.NA + F.1.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.1.saline.Caucasian.Hispanic + F.1.saline.Caucasian.Not.Hispanic.or.Latino)/6-
       (F..7.Flu.African.American.NA + F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Asian.NA + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino)/6),
                              IMaleBaseVSthree = ((M.3.Flu.Asian.NA + M.3.Flu.Caucasian.Not.Hispanic.or.Latino + M.3.Pneumovax.Asian.NA + M.3.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.3.saline.Asian.NA + M.3.saline.Caucasian.Not.Hispanic.or.Latino)/6-
    (M..7.Flu.Asian.NA + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Asian.NA + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Asian.NA + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/6) -
    ((F.3.Flu.African.American.NA + F.3.Flu.Caucasian.Not.Hispanic.or.Latino + F.3.Pneumovax.Asian.NA + F.3.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.3.saline.Caucasian.Hispanic + F.3.saline.Caucasian.Not.Hispanic.or.Latino)/6-
       (F..7.Flu.African.American.NA + F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Asian.NA + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino)/6),
                              IMaleBaseVSseven = ((M.7.Flu.Asian.NA + M.7.Flu.Caucasian.Not.Hispanic.or.Latino + M.7.Pneumovax.Asian.NA + M.7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.7.saline.Asian.NA + M.7.saline.Caucasian.Not.Hispanic.or.Latino)/6-
    (M..7.Flu.Asian.NA + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Asian.NA + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Asian.NA + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/6) -
    ((F.7.Flu.African.American.NA + F.7.Flu.Caucasian.Not.Hispanic.or.Latino + F.7.Pneumovax.Asian.NA + F.7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.7.saline.Caucasian.Hispanic + F.7.saline.Caucasian.Not.Hispanic.or.Latino)/6-
       (F..7.Flu.African.American.NA + F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Asian.NA + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino)/6),
                              IMaleBaseVSten = ((M.10.Flu.Asian.NA + M.10.Flu.Caucasian.Not.Hispanic.or.Latino + M.10.Pneumovax.Asian.NA + M.10.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.10.saline.Asian.NA + M.10.saline.Caucasian.Not.Hispanic.or.Latino)/6-
    (M..7.Flu.Asian.NA + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Asian.NA + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Asian.NA + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/6) -
    ((F.10.Flu.African.American.NA + F.10.Flu.Caucasian.Not.Hispanic.or.Latino + F.10.Pneumovax.Asian.NA + F.10.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.10.saline.Caucasian.Hispanic + F.10.saline.Caucasian.Not.Hispanic.or.Latino)/6-
       (F..7.Flu.African.American.NA + F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Asian.NA + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino)/6),
                              IMaleBaseVSfourteen = ((M.14.Flu.Asian.NA + M.14.Flu.Caucasian.Not.Hispanic.or.Latino + M.14.Pneumovax.Asian.NA + M.14.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.14.saline.Asian.NA + M.14.saline.Caucasian.Not.Hispanic.or.Latino)/6-
    (M..7.Flu.Asian.NA + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Asian.NA + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Asian.NA + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/6) -
    ((F.14.Flu.African.American.NA + F.14.Flu.Caucasian.Not.Hispanic.or.Latino + F.14.Pneumovax.Asian.NA + F.14.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.14.saline.Caucasian.Hispanic + F.14.saline.Caucasian.Not.Hispanic.or.Latino)/6-
       (F..7.Flu.African.American.NA + F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Asian.NA + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino)/6),
                              IMaleBaseVStwentyOne = ((M.21.Flu.Asian.NA + M.21.Flu.Caucasian.Not.Hispanic.or.Latino + M.21.Pneumovax.Asian.NA + M.21.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.21.saline.Asian.NA + M.21.saline.Caucasian.Not.Hispanic.or.Latino)/6-
    (M..7.Flu.Asian.NA + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Asian.NA + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Asian.NA + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/6) -
    ((F.21.Flu.African.American.NA + F.21.Flu.Caucasian.Not.Hispanic.or.Latino + F.21.Pneumovax.Asian.NA + F.21.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.21.saline.Caucasian.Hispanic + F.21.saline.Caucasian.Not.Hispanic.or.Latino)/6-
       (F..7.Flu.African.American.NA + F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Asian.NA + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino)/6),
                              IMaleBaseVStwentyEight = ((M.28.Flu.Asian.NA + M.28.Flu.Caucasian.Not.Hispanic.or.Latino + M.28.Pneumovax.Asian.NA + M.28.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.28.saline.Asian.NA + M.28.saline.Caucasian.Not.Hispanic.or.Latino)/6-
    (M..7.Flu.Asian.NA + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Asian.NA + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Asian.NA + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/6) -
    ((F.28.Flu.African.American.NA + F.28.Flu.Caucasian.Not.Hispanic.or.Latino + F.28.Pneumovax.Asian.NA + F.28.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.28.saline.Caucasian.Hispanic + F.28.saline.Caucasian.Not.Hispanic.or.Latino)/6-
       (F..7.Flu.African.American.NA + F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Asian.NA + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino)/6),
                              # Caucasian only
  CaucBaseVSzero = (F.0.Flu.Caucasian.Not.Hispanic.or.Latino + F.0.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.0.saline.Caucasian.Hispanic + F.0.saline.Caucasian.Not.Hispanic.or.Latino + M.0.Flu.Caucasian.Not.Hispanic.or.Latino + M.0.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.0.saline.Caucasian.Not.Hispanic.or.Latino)/7-
    (F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/7,
  CaucBaseVSone = (F.1.Flu.Caucasian.Not.Hispanic.or.Latino + F.1.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.1.saline.Caucasian.Hispanic + F.1.saline.Caucasian.Not.Hispanic.or.Latino + M.1.Flu.Caucasian.Not.Hispanic.or.Latino + M.1.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.1.saline.Caucasian.Not.Hispanic.or.Latino)/7-
    (F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/7,
    CaucBaseVSthree = (F.3.Flu.Caucasian.Not.Hispanic.or.Latino + F.3.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.3.saline.Caucasian.Hispanic + F.3.saline.Caucasian.Not.Hispanic.or.Latino + M.3.Flu.Caucasian.Not.Hispanic.or.Latino + M.3.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.3.saline.Caucasian.Not.Hispanic.or.Latino)/7-
      (F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/7,
    CaucBaseVSseven = (F.7.Flu.Caucasian.Not.Hispanic.or.Latino + F.7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.7.saline.Caucasian.Hispanic + F.7.saline.Caucasian.Not.Hispanic.or.Latino + M.7.Flu.Caucasian.Not.Hispanic.or.Latino + M.7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.7.saline.Caucasian.Not.Hispanic.or.Latino)/7-
                       (F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/7,
  CaucBaseVSten = (F.10.Flu.Caucasian.Not.Hispanic.or.Latino + F.10.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.10.saline.Caucasian.Hispanic + F.10.saline.Caucasian.Not.Hispanic.or.Latino + M.10.Flu.Caucasian.Not.Hispanic.or.Latino + M.10.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.10.saline.Caucasian.Not.Hispanic.or.Latino)/7-
    (F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/7,
  CaucBaseVSfourteen = (F.14.Flu.Caucasian.Not.Hispanic.or.Latino + F.14.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.14.saline.Caucasian.Hispanic + F.14.saline.Caucasian.Not.Hispanic.or.Latino + M.14.Flu.Caucasian.Not.Hispanic.or.Latino + M.14.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.14.saline.Caucasian.Not.Hispanic.or.Latino)/7-
    (F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/7,
  CaucBaseVStwentyOne = (F.21.Flu.Caucasian.Not.Hispanic.or.Latino + F.21.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.21.saline.Caucasian.Hispanic + F.21.saline.Caucasian.Not.Hispanic.or.Latino + M.21.Flu.Caucasian.Not.Hispanic.or.Latino + M.21.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.21.saline.Caucasian.Not.Hispanic.or.Latino)/7-
    (F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/7,
  CaucBaseVStwentyEight = (F.28.Flu.Caucasian.Not.Hispanic.or.Latino + F.28.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.28.saline.Caucasian.Hispanic + F.28.saline.Caucasian.Not.Hispanic.or.Latino + M.28.Flu.Caucasian.Not.Hispanic.or.Latino + M.28.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.28.saline.Caucasian.Not.Hispanic.or.Latino)/7-
    (F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Hispanic + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/7,
                              # Asian only
  AsianBaseVSzero = (F.0.Pneumovax.Asian.NA + M.0.Flu.Asian.NA + M.0.Pneumovax.Asian.NA + M.0.saline.Asian.NA)/4-
    (F..7.Pneumovax.Asian.NA + M..7.Flu.Asian.NA + M..7.Pneumovax.Asian.NA + M..7.saline.Asian.NA)/4,
  AsianBaseVSone = (F.1.Pneumovax.Asian.NA + M.1.Flu.Asian.NA + M.1.Pneumovax.Asian.NA + M.1.saline.Asian.NA)/4-
    (F..7.Pneumovax.Asian.NA + M..7.Flu.Asian.NA + M..7.Pneumovax.Asian.NA + M..7.saline.Asian.NA)/4,
  AsianBaseVSthree = (F.3.Pneumovax.Asian.NA + M.3.Flu.Asian.NA + M.3.Pneumovax.Asian.NA + M.3.saline.Asian.NA)/4-
    (F..7.Pneumovax.Asian.NA + M..7.Flu.Asian.NA + M..7.Pneumovax.Asian.NA + M..7.saline.Asian.NA)/4,
  AsianBaseVSseven = (F.7.Pneumovax.Asian.NA + M.7.Flu.Asian.NA + M.7.Pneumovax.Asian.NA + M.7.saline.Asian.NA)/4-
    (F..7.Pneumovax.Asian.NA + M..7.Flu.Asian.NA + M..7.Pneumovax.Asian.NA + M..7.saline.Asian.NA)/4,
  AsianBaseVSten = (F.10.Pneumovax.Asian.NA + M.10.Flu.Asian.NA + M.10.Pneumovax.Asian.NA + M.10.saline.Asian.NA)/4-
    (F..7.Pneumovax.Asian.NA + M..7.Flu.Asian.NA + M..7.Pneumovax.Asian.NA + M..7.saline.Asian.NA)/4,
  AsianBaseVSfourteen = (F.14.Pneumovax.Asian.NA + M.14.Flu.Asian.NA + M.14.Pneumovax.Asian.NA + M.14.saline.Asian.NA)/4-
    (F..7.Pneumovax.Asian.NA + M..7.Flu.Asian.NA + M..7.Pneumovax.Asian.NA + M..7.saline.Asian.NA)/4,
  AsianBaseVStwentyOne = (F.21.Pneumovax.Asian.NA + M.21.Flu.Asian.NA + M.21.Pneumovax.Asian.NA + M.21.saline.Asian.NA)/4-
    (F..7.Pneumovax.Asian.NA + M..7.Flu.Asian.NA + M..7.Pneumovax.Asian.NA + M..7.saline.Asian.NA)/4,
  AsianBaseVStwentyEight = (F.28.Pneumovax.Asian.NA + M.28.Flu.Asian.NA + M.28.Pneumovax.Asian.NA + M.28.saline.Asian.NA)/4-
    (F..7.Pneumovax.Asian.NA + M..7.Flu.Asian.NA + M..7.Pneumovax.Asian.NA + M..7.saline.Asian.NA)/4,
                              # African American only
  AfricanAmerBaseVSzero = (F.0.Flu.African.American.NA) - (F..7.Flu.African.American.NA),
  AfricanAmerBaseVSone = (F.1.Flu.African.American.NA) - (F..7.Flu.African.American.NA),
  AfricanAmerBaseVSthree = (F.3.Flu.African.American.NA) - (F..7.Flu.African.American.NA),
  AfricanAmerBaseVSseven = (F.7.Flu.African.American.NA) - (F..7.Flu.African.American.NA),
  AfricanAmerBaseVSten = (F.10.Flu.African.American.NA) - (F..7.Flu.African.American.NA),
  AfricanAmerBaseVSfourteen = (F.14.Flu.African.American.NA) - (F..7.Flu.African.American.NA),
  AfricanAmerBaseVStwentyOne = (F.21.Flu.African.American.NA) - (F..7.Flu.African.American.NA),
  AfricanAmerBaseVStwentyEight = (F.28.Flu.African.American.NA) - (F..7.Flu.African.American.NA),
                              # Time Hispanic - Non-hispanic
  IHispanicBaseVSzero = ((F.0.saline.Caucasian.Hispanic) - (F..7.saline.Caucasian.Hispanic)) - 
    ((F.0.Flu.Caucasian.Not.Hispanic.or.Latino + F.0.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.0.saline.Caucasian.Not.Hispanic.or.Latino + M.0.Flu.Caucasian.Not.Hispanic.or.Latino + M.0.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.0.saline.Caucasian.Not.Hispanic.or.Latino)/6 -
        (F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/6),
  IHispanicBaseVSone = ((F.1.saline.Caucasian.Hispanic) - (F..7.saline.Caucasian.Hispanic)) - 
    ((F.1.Flu.Caucasian.Not.Hispanic.or.Latino + F.1.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.1.saline.Caucasian.Not.Hispanic.or.Latino + M.1.Flu.Caucasian.Not.Hispanic.or.Latino + M.1.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.1.saline.Caucasian.Not.Hispanic.or.Latino)/6 -
        (F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/6),
  IHispanicBaseVSthree = ((F.3.saline.Caucasian.Hispanic) - (F..7.saline.Caucasian.Hispanic)) - 
    ((F.3.Flu.Caucasian.Not.Hispanic.or.Latino + F.3.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.3.saline.Caucasian.Not.Hispanic.or.Latino + M.3.Flu.Caucasian.Not.Hispanic.or.Latino + M.3.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.3.saline.Caucasian.Not.Hispanic.or.Latino)/6 -
        (F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/6),
  IHispanicBaseVSseven = ((F.7.saline.Caucasian.Hispanic) - (F..7.saline.Caucasian.Hispanic)) - 
    ((F.7.Flu.Caucasian.Not.Hispanic.or.Latino + F.7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.7.saline.Caucasian.Not.Hispanic.or.Latino + M.7.Flu.Caucasian.Not.Hispanic.or.Latino + M.7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.7.saline.Caucasian.Not.Hispanic.or.Latino)/6 -
        (F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/6),
  IHispanicBaseVSten = ((F.10.saline.Caucasian.Hispanic) - (F..7.saline.Caucasian.Hispanic)) - 
    ((F.10.Flu.Caucasian.Not.Hispanic.or.Latino + F.10.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.10.saline.Caucasian.Not.Hispanic.or.Latino + M.10.Flu.Caucasian.Not.Hispanic.or.Latino + M.10.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.10.saline.Caucasian.Not.Hispanic.or.Latino)/6 -
        (F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/6),
  IHispanicBaseVSfourteen = ((F.14.saline.Caucasian.Hispanic) - (F..7.saline.Caucasian.Hispanic)) - 
    ((F.14.Flu.Caucasian.Not.Hispanic.or.Latino + F.14.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.14.saline.Caucasian.Not.Hispanic.or.Latino + M.14.Flu.Caucasian.Not.Hispanic.or.Latino + M.14.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.14.saline.Caucasian.Not.Hispanic.or.Latino)/6 -
        (F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/6),
  IHispanicBaseVStwentyOne = ((F.21.saline.Caucasian.Hispanic) - (F..7.saline.Caucasian.Hispanic)) - 
    ((F.21.Flu.Caucasian.Not.Hispanic.or.Latino + F.21.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.21.saline.Caucasian.Not.Hispanic.or.Latino + M.21.Flu.Caucasian.Not.Hispanic.or.Latino + M.21.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.21.saline.Caucasian.Not.Hispanic.or.Latino)/6 -
        (F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/6),
  IHispanicBaseVStwentyEight = ((F.28.saline.Caucasian.Hispanic) - (F..7.saline.Caucasian.Hispanic)) - 
    ((F.28.Flu.Caucasian.Not.Hispanic.or.Latino + F.28.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F.28.saline.Caucasian.Not.Hispanic.or.Latino + M.28.Flu.Caucasian.Not.Hispanic.or.Latino + M.28.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M.28.saline.Caucasian.Not.Hispanic.or.Latino)/6 -
        (F..7.Flu.Caucasian.Not.Hispanic.or.Latino + F..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + F..7.saline.Caucasian.Not.Hispanic.or.Latino + M..7.Flu.Caucasian.Not.Hispanic.or.Latino + M..7.Pneumovax.Caucasian.Not.Hispanic.or.Latino + M..7.saline.Caucasian.Not.Hispanic.or.Latino)/6),
                              levels = colnames(design))

fit2 <- contrasts.fit(fit, my_contrasts)
x <- eBayes(fit2)

### Results tables ----
y <- topTable(x, number = Inf, adjust.method = "fdr")
y$Symbol <- mapIds(x = illuminaHumanv3.db, keys = row.names(y), column = "SYMBOL", keytype = "PROBEID")
dir.create(paste(output_dir, "training_vein_limma", sep = ""))
write.csv(y, file = paste(output_dir, "training_vein_limma/limma_F.csv", sep = ""), row.names = TRUE)

dt <- decideTests(x)
summary(dt)
#        BaseVSzeroALL BaseVSoneALL BaseVSthreeALL BaseVSsevenALL BaseVStenALL BaseVSfourteenALL BaseVStwentyOneALL BaseVStwentyEightALL
# Down               0          371              0            961            0                 0                  1                    3   
# NotSig         31003        30415          31003          29735        31001             31003              31001                30998
# Up                 0          217              0            307            2                 0                  1                    2 
### many more
dt <- decideTests(x, adjust.method = 'none')
summary(dt)
#        BaseVSzeroALL BaseVSoneALL BaseVSthreeALL BaseVSsevenALL BaseVStenALL BaseVSfourteenALL BaseVStwentyOneALL BaseVStwentyEightALL 
# Down             734         2554            845           2608         1011               895               1338                 1043 
# NotSig         29621        26515          29457          26419        28988             29393              28587                29238 
# Up               648         1934            701           1976         1004               715               1078                  722 
### many more

names(y) # 72 columns
tables <- list()
for (i in 1:72) {
  tables[[i]] <- topTable(x, coef = i, number = Inf, adjust.method = "fdr")
  tables[[i]]$Symbol <- mapIds(illuminaHumanv3.db, keys = rownames(tables[[i]]), column = "SYMBOL", keytype = "PROBEID")
  names(tables)[[i]] <- names(y)[i]
  write.csv(tables[[i]], file = paste(output_dir, "training_vein_limma/",
                                      names(y)[i],".csv", sep = ""), row.names = TRUE)
}

### Plot p-value histograms for timepoints ----
histlist <- list()
for (i in 1:8) {
  histlist[[i]] <- ggplot(data = tables[[i]], aes(x = adj.P.Val)) +
    geom_histogram(fill = "lightblue", color = "black") +
    ggtitle(paste("Adj. p histogram", names(y)[i])) +
    theme_cowplot() +
    theme(title = element_text(size = 7.5))
  names(histlist)[i] <- names(y)[i]
}
pdf(file = paste(output_dir, "training_vein_limma/p_value_histogram.pdf", sep = ""), width = 12)
plot_grid(plotlist = histlist[grep("^Base", names(histlist))], align = 'hv',  nrow = 2)
dev.off()
for (i in 1:8) {
  histlist[[i]] <- ggplot(data = tables[[i]], aes(x = P.Value)) +
    geom_histogram(fill = "lightblue", color = "black") +
    ggtitle(paste("P-value histogram", names(y)[i])) +
    theme_cowplot() +
    theme(title = element_text(size = 7.5))
  names(histlist)[i] <- names(y)[i]
}
pdf(file = paste(output_dir, "training_vein_limma/adj_p_value_histogram.pdf", sep = ""), width = 12)
plot_grid(plotlist = histlist[grep("^Base", names(histlist))], align = 'hv',  nrow = 2)
dev.off()

### Create volcano plots ----
volcanos <- list()
for (i in 1:8) {
  volcanos[[i]] <- EnhancedVolcano(tables[[i]],
                                   lab = tables[[i]]$Symbol,
                                   x = "logFC",
                                   y = "P.Value",
                                   title = (paste("Training_Vein:", names(y)[i])), titleLabSize = 6,
                                   subtitle = "",
                                   labSize = (4.0),
                                   colAlpha = 1,
                                   pointSize = 1.0, 
                                   pCutoff = 1e-2,
                                   boxedLabels = FALSE,
                                   drawConnectors = TRUE,
                                   legendLabSize = 4,
                                   legendIconSize = 2,
                                   axisLabSize = 5,
                                   caption = "")
  names(volcanos)[i] <- names(y)[i]
}
pdf(file = paste(output_dir, "training_vein_limma/volcanos.pdf", sep = ""), width = 12, height = 9)
plot_grid(plotlist = volcanos[grep("^Base", names(volcanos))], align = 'hv', nrow = 2)
dev.off()

# Analysis 2: Test_Set_Finger (p1f) ----
## Step 1: Read in raw expression data ----
p1f[1:5,1:5]
#                                    title geo_accession subject age time
# GSM745399  009-282-018-TP2 [5447320035B]     GSM745399      18  23    0
# GSM745400  009-282-010-TP1 [5447320035G]     GSM745400      10  29   -7
# GSM745401  009-282-024-TP1 [5447320036A]     GSM745401      24  39   -7
# GSM745402 009-282-010-TP13 [5447320036B]     GSM745402      10  29    7
# GSM745403  009-282-018-TP1 [5447320036C]     GSM745403      18  23   -7
levels(factor(p1f$time))
# [1] "-7" "0"  "7" 
## It appears the time units are in days, and 0 represents the day vaccination was administered.
dim(p1f)
# [1] 49 12
table(p1f$time, p1f$subject)
#    1 10 11 13 14 15 18 20 21 22 24 26 3 4 6 7 8
# -7 1  1  1  1  1  1  1  0  1  1  1  1 1 1 1 1 1
# 0  1  1  1  1  1  1  1  1  1  1  1  1 1 1 1 1 0
# 7  1  1  1  1  1  1  1  1  1  1  1  1 1 1 1 1 1
table(p1f$sex, p1f$subject) ## Sanity check.
#   1 10 11 13 14 15 18 20 21 22 24 26 3 4 6 7 8
# F 3  0  0  3  0  0  3  0  0  3  0  3 3 3 3 3 0
# M 0  3  3  0  3  3  0  2  3  0  3  0 0 0 0 0 2
table(p1f$treat, p1f$subject)
#         1 10 11 13 14 15 18 20 21 22 24 26 3 4 6 7 8
# FLUZONE 3  3  0  0  0  0  3  0  0  3  3  0 0 0 0 0 2
# NS      0  0  0  3  3  0  0  2  0  0  0  0 3 0 3 0 0
# PNEUM   0  0  3  0  0  3  0  0  3  0  0  3 0 3 0 3 0

## Load in expression data
temp <- files[str_extract(files, "GSM......") %in% p1f$geo_accession]
length(temp) # [1] 49
x <- lapply(
  temp,
  function(x) read.delim(file = x))
x <- do.call(cbind, x)
x[1:5,1:5]
#              X X5447320035_B.AVG_Signal X5447320035_B.Detection_Pval          X.1 X5447320035_G.AVG_Signal
# 1 ILMN_1885639                   2.4490                       0.3294 ILMN_1885639                  -2.4503
# 2 ILMN_2209417                  -7.6996                       0.7945 ILMN_2209417                 -15.4126
# 3 ILMN_1765401                  26.4028                       0.0184 ILMN_1765401                  -0.8125
# 4 ILMN_1652171                 -12.6385                       0.9341 ILMN_1652171                 -12.7959
# 5 ILMN_1652170                   0.9654                       0.3992 ILMN_1652170                  -0.4988
dim(x)
# [1] 48803   147
rownames(x) <- x[,1]
x <- as.data.frame(x)
detectionpvalues <- x[,grep("Detection",colnames(x))]
x <- x[,grep("AVG",colnames(x))]
x[1:2,1:5]
#              X5447320035_B.AVG_Signal X5447320035_G.AVG_Signal X5447320036_A.AVG_Signal X5447320036_B.AVG_Signal
# ILMN_1885639                   2.4490                  -2.4503                  -5.5602                  -3.5167
# ILMN_2209417                  -7.6996                 -15.4126                  -6.5398                 -10.1059
#              X5447320036_C.AVG_Signal
# ILMN_1885639                  -0.9949
# ILMN_2209417                 -19.1151

## Step 2: Create an EList object ----
data <- new('EListRaw')
data@.Data[[1]] <- 'illumina'
# data@.Data[[2]] <- targetinfo
data@.Data[[3]] <- BGX
data@.Data[[4]] <- x
data@.Data[[5]] <- NULL
data$E <- x
# data$targets <- targetinfo
data$genes <- BGX
data$other$Detection <- detectionpvalues

### Background correction and normalization ----
data_norm <- neqc(data)
table(data_norm$genes$controls$Reporter_Group_Name)
# biotin            cy3_hyb       housekeeping           labeling low_stringency_hyb 
# 2                  6                  7                  2                  8 
# negative 
# 759 

## Remove unexpressed probes and controls
controls <- data_norm$genes$controls$Probe_Id
expressed <- rowSums(data_norm$other$Detection < 0.05) >= 3
data_norm$E <- data_norm$E[expressed,]
data_norm$E <- data_norm$E[-which(rownames(data_norm$E) %in% controls),]
dim(data_norm$E)
# [1] 22413    49

## Step 3: PCA and MDS analysis ----
identical(substr(p1f$array, start = 1, stop = 10), str_extract(colnames(data_norm$E),"\\d........."))
# [1] TRUE
### PCAtools ----
p <- pca(data_norm$E, metadata = p1f, center = T, scale = T)

scree <- screeplot(p, getComponents(p, c(1:10)))
ptime <- biplot(p, 
                showLoadings = F,
                lab = NULL, 
                colby = "time", 
                legendPosition = 'right', legendLabSize = 12, legendIconSize = 3,
                title = 'Obermoser QC',
                subtitle = 'p1f samples')
psex <- biplot(p, 
               showLoadings = F,
               lab = NULL, 
               colby = "sex", 
               legendPosition = 'right', legendLabSize = 9, legendIconSize = 3)
psub <- biplot(p, 
               showLoadings = F,
               lab = NULL, 
               colby = "subject",
               legendPosition = 'right', legendLabSize = 9, legendIconSize = 3)
ptre <-  biplot(p, 
                showLoadings = F,
                lab = NULL, 
                colby = "treat",
                legendPosition = 'right', legendLabSize = 9, legendIconSize = 3)
peigencor <- eigencorplot(p,
                          components = getComponents(p, 1:10),
                          metavars = c("time", 'sex', 'race', 'ethnicity', 'treat'),
                          cexCorval = 0.5,
                          fontCorval = 1, 
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
botrow <- plot_grid(ptre, scree, peigencor, align = 'h', nrow = 1)
pdf(file = paste(output_dir, "pca_p1f_QC.pdf", sep = ""), height = 12, width = 16)
plot_grid(toprow, botrow, ncol = 1, align = 'v')
dev.off()

### MDS ----
### MDS can be good for getting really nice separation and pulling out outliers.
qual_col_pals <- brewer.pal.info
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colors <- sample(col_vector, length(unique(p1f$subject)))
new_colors <- c()
for (i in 1:length(p1f$subject)) {
  new_colors[i] <- switch(p1f$subject[i],
                          "18" = colors[1], "10" = colors[2], "24" = colors[3], "1" = colors[4], 
                          "22" = colors[5], "8" = colors[6], "6" = colors[7], "13" = colors[8],
                          "3" = colors[9], "14" = colors[10], "20" = colors[11], "7" = colors[12],
                          "4" = colors[13], "15" = colors[14], "21" = colors[15], "26" = colors[16],
                          "11" = colors[17]
  )
}
pdf(file = paste(output_dir, "mds_p1f_qc.pdf", sep = ""))
limma::plotMDS(x = data_norm$E, 
               cex = 1, 
               var.explained = T, 
               col = new_colors,
               labels = paste(p1f$subject, p1f$sex, p1f$treat, sep = ", "))
title("MDS Plot - test set finger")
dev.off()

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
## Step 4: limma analysis ----
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

id <- factor(p1f$subject)
time <- factor(p1f$time)
sex <- factor(p1f$sex)
vax <- factor(p1f$treat)
eth <- factor(p1f$ethnicity)
race <- factor(p1f$race)

df <- data.frame(data_norm$E)

TS <- factor(paste(sex,time,vax,race,eth, sep = "."))
design <- model.matrix(~0+TS, df)
colnames(design) <- make.names(levels(TS))
cor <- duplicateCorrelation(df, design, block=id)
cor$consensus.correlation
# [1] 0.1753894

fit <- lmFit(object=df, design=design, block=id, correlation=cor$consensus.correlation)

# contrasts
my_contrasts <- makeContrasts(# Time main effect
                              BaselineVSzeroAll = (F.0.FLUZONE.Asian.NA+F.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino+F.0.NS.African.American.NA+F.0.NS.Caucasian.Hispanic+F.0.NS.Caucasian.Non.Hispanic.or.Latino+F.0.PNEUM.Caucasian.Non.Hispanic.or.Latino+M.0.FLUZONE.Asian.NA+M.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino+M.0.NS.Caucasian.Non.Hispanic.or.Latino+M.0.PNEUM.African.American.NA+M.0.PNEUM.Caucasian.Hispanic+M.0.PNEUM.Caucasian.Non.Hispanic.or.Latino)/12 -
                                (F..7.FLUZONE.Asian.NA+F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino+F..7.NS.African.American.NA+F..7.NS.Caucasian.Hispanic+F..7.NS.Caucasian.Non.Hispanic.or.Latino+F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino+M..7.FLUZONE.Asian.NA+M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino+M..7.NS.Caucasian.Non.Hispanic.or.Latino+M..7.PNEUM.African.American.NA+M..7.PNEUM.Caucasian.Hispanic+M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino)/12,
                              BaselineVSsevenAll = (F.7.FLUZONE.Asian.NA+F.7.FLUZONE.Caucasian.Non.Hispanic.or.Latino+F.7.NS.African.American.NA+F.7.NS.Caucasian.Hispanic+F.7.NS.Caucasian.Non.Hispanic.or.Latino+F.7.PNEUM.Caucasian.Non.Hispanic.or.Latino+M.7.FLUZONE.Asian.NA+M.7.FLUZONE.Caucasian.Non.Hispanic.or.Latino+M.7.NS.Caucasian.Non.Hispanic.or.Latino+M.7.PNEUM.African.American.NA+M.7.PNEUM.Caucasian.Hispanic+M.7.PNEUM.Caucasian.Non.Hispanic.or.Latino)/12 -
                                (F..7.FLUZONE.Asian.NA+F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino+F..7.NS.African.American.NA+F..7.NS.Caucasian.Hispanic+F..7.NS.Caucasian.Non.Hispanic.or.Latino+F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino+M..7.FLUZONE.Asian.NA+M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino+M..7.NS.Caucasian.Non.Hispanic.or.Latino+M..7.PNEUM.African.American.NA+M..7.PNEUM.Caucasian.Hispanic+M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino)/12,
                              # Saline times
                              SalineBaseVSZero = (F.0.NS.African.American.NA + F.0.NS.Caucasian.Hispanic + F.0.NS.Caucasian.Non.Hispanic.or.Latino + M.0.NS.Caucasian.Non.Hispanic.or.Latino)/4 -
                                (F..7.NS.African.American.NA + F..7.NS.Caucasian.Hispanic + F..7.NS.Caucasian.Non.Hispanic.or.Latino + M..7.NS.Caucasian.Non.Hispanic.or.Latino)/4,
                              SalineBaseVSseven = (F.7.NS.African.American.NA + F.7.NS.Caucasian.Hispanic + F.7.NS.Caucasian.Non.Hispanic.or.Latino + M.7.NS.Caucasian.Non.Hispanic.or.Latino)/4 - 
                                (F..7.NS.African.American.NA + F..7.NS.Caucasian.Hispanic + F..7.NS.Caucasian.Non.Hispanic.or.Latino + M..7.NS.Caucasian.Non.Hispanic.or.Latino)/4,
                              # Pneum times
                              PNEUMBaseVSZero = (F.0.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.0.PNEUM.African.American.NA + M.0.PNEUM.Caucasian.Hispanic + M.0.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4 -
                                (F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.African.American.NA + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4,
                              PNEUMBaseVSseven = (F.7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.7.PNEUM.African.American.NA + M.7.PNEUM.Caucasian.Hispanic + M.7.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4 -
                                (F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.African.American.NA + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4,
                              # Flu times
                              FLUBaseVSZero = (F.0.FLUZONE.Asian.NA + F.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.0.FLUZONE.Asian.NA + M.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4 -
                                (F..7.FLUZONE.Asian.NA + F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Asian.NA + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4,
                              FLUBaseVSseven = (F.7.FLUZONE.Asian.NA + F.7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.7.FLUZONE.Asian.NA + M.7.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4 -
                                (F..7.FLUZONE.Asian.NA + F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Asian.NA + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4,
                              # Time-sex interaction, Male - Female
                              IMaleBaseVSZero = ((M.0.FLUZONE.Asian.NA + M.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.0.NS.Caucasian.Non.Hispanic.or.Latino + M.0.PNEUM.African.American.NA + M.0.PNEUM.Caucasian.Hispanic + M.0.PNEUM.Caucasian.Non.Hispanic.or.Latino)/6-
                                                   (M..7.FLUZONE.Asian.NA + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.NS.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.African.American.NA + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino)/6) -
                                ((F.0.FLUZONE.Asian.NA + F.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.0.NS.African.American.NA + F.0.NS.Caucasian.Hispanic + F.0.NS.Caucasian.Non.Hispanic.or.Latino + F.0.PNEUM.Caucasian.Non.Hispanic.or.Latino)/6-
                                   (F..7.FLUZONE.Asian.NA + F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.NS.African.American.NA + F..7.NS.Caucasian.Hispanic + F..7.NS.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino)/6),
                              IMaleBaseVSseven = ((M.7.FLUZONE.Asian.NA + M.7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.7.NS.Caucasian.Non.Hispanic.or.Latino + M.7.PNEUM.African.American.NA + M.7.PNEUM.Caucasian.Hispanic + M.7.PNEUM.Caucasian.Non.Hispanic.or.Latino)/6-
                                                   (M..7.FLUZONE.Asian.NA + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.NS.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.African.American.NA + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino)/6) -
                                ((F.7.FLUZONE.Asian.NA + F.7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.7.NS.African.American.NA + F.7.NS.Caucasian.Hispanic + F.7.NS.Caucasian.Non.Hispanic.or.Latino + F.7.PNEUM.Caucasian.Non.Hispanic.or.Latino)/6-
                                   (F..7.FLUZONE.Asian.NA + F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.NS.African.American.NA + F..7.NS.Caucasian.Hispanic + F..7.NS.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino)/6),
                              # Caucasian only
                              CaucBaseVSZero = (F.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.0.NS.Caucasian.Hispanic + F.0.NS.Caucasian.Non.Hispanic.or.Latino + F.0.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.0.NS.Caucasian.Non.Hispanic.or.Latino + M.0.PNEUM.Caucasian.Hispanic + M.0.PNEUM.Caucasian.Non.Hispanic.or.Latino)/8-
                                (F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.NS.Caucasian.Hispanic + F..7.NS.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.NS.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino)/8,
                              CaucBaseVSseven = (F.7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.7.NS.Caucasian.Hispanic + F.7.NS.Caucasian.Non.Hispanic.or.Latino + F.7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.7.NS.Caucasian.Non.Hispanic.or.Latino + M.7.PNEUM.Caucasian.Hispanic + M.7.PNEUM.Caucasian.Non.Hispanic.or.Latino)/8-
                                (F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.NS.Caucasian.Hispanic + F..7.NS.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.NS.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino)/8,
                              # Asian only
                              AsianBaseVSZero = (F.0.FLUZONE.Asian.NA + M.0.FLUZONE.Asian.NA)/2 - (F..7.FLUZONE.Asian.NA + M..7.FLUZONE.Asian.NA)/2,
                              AsianBaseVSseven = (F.7.FLUZONE.Asian.NA + M.7.FLUZONE.Asian.NA)/2 - (F..7.FLUZONE.Asian.NA + M..7.FLUZONE.Asian.NA)/2,
                              # African American only
                              AfricanAmerBaseVSZero = (F.0.NS.African.American.NA + M.0.PNEUM.African.American.NA)/2 - (F..7.NS.African.American.NA + M..7.PNEUM.African.American.NA)/2,
                              AfricanAmerBaseVSseven = (F.7.NS.African.American.NA + M.7.PNEUM.African.American.NA)/2 - (F..7.NS.African.American.NA + M..7.PNEUM.African.American.NA)/2,
                              # Time-ethnicity interaction, Hispanic - non hispanic
                              IHispanicBaseVSZero = ((F.0.NS.Caucasian.Hispanic + M.0.PNEUM.Caucasian.Hispanic)/2 - (F..7.NS.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Hispanic)/2) -
                                ((F.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.0.NS.Caucasian.Non.Hispanic.or.Latino + F.0.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.0.NS.Caucasian.Non.Hispanic.or.Latino + M.0.PNEUM.Caucasian.Non.Hispanic.or.Latino)/6-
                                   (F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.NS.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.NS.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino)/6),
                              IHispanicBaseVSseven = ((F.7.NS.Caucasian.Hispanic + M.7.PNEUM.Caucasian.Hispanic)/2 - (F..7.NS.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Hispanic)/2) -
                                ((F.7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.7.NS.Caucasian.Non.Hispanic.or.Latino + F.7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.7.NS.Caucasian.Non.Hispanic.or.Latino + M.7.PNEUM.Caucasian.Non.Hispanic.or.Latino)/6-
                                   (F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.NS.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.NS.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino)/6),
                                levels = colnames(design))

fit2 <- contrasts.fit(fit, my_contrasts)
x <- eBayes(fit2)

### Results tables ----

y <- topTable(x, number = Inf, adjust.method = "fdr")
y$Symbol <- mapIds(x = illuminaHumanv3.db, keys = row.names(y), column = "SYMBOL", keytype = "PROBEID")
dir.create(paste(output_dir, "test_finger_limma", sep = ""))
write.csv(y, file = paste(output_dir, "test_finger_limma/limma_F.csv", sep = ""), row.names = TRUE)

dt <- decideTests(x)
summary(dt)
#       BaselineVSzeroAll BaselineVSsevenAll 
# Down                   0                  3
# NotSig             22413              22386
# Up                     0                 24
dt <- decideTests(x, adjust.method = 'none')
summary(dt)
#        BaselineVSzeroAll BaselineVSsevenAll 
# Down                 312               1345 
# NotSig             21731              19883 
# Up                   370               1185 

names(y)
tables <- list()
for (i in 1:18) {
  tables[[i]] <- topTable(x, coef = i, number = Inf, adjust.method = "fdr")
  tables[[i]]$Symbol <- mapIds(x = illuminaHumanv3.db, keys = row.names(tables[[i]]), column = "SYMBOL", keytype = "PROBEID")
  write.csv(tables[[i]], file = paste(output_dir,"test_finger_limma/",names(y)[i],".csv",sep = ""), row.names = TRUE)
}

### Plot p-value histograms ----
histlist <- list()
for (i in 1:2) {
  histlist[[i]] <- ggplot(data = tables[[i]], aes(x = adj.P.Val)) +
    geom_histogram(fill = "lightblue", color = "black") +
    ggtitle(paste("Adj. p histogram", names(y)[i])) +
    theme_cowplot() +
    theme(title = element_text(size = 7.5))
}
for (i in 1:2) {
  histlist[[i + 2]] <- ggplot(data = tables[[i]], aes(x = P.Value)) +
    geom_histogram(fill = "lightblue", color = "black") +
    ggtitle(paste("P-value histogram", names(y)[i])) +
    theme_cowplot() +
    theme(title = element_text(size = 7.5))
}

toprow <- plot_grid(plotlist = histlist[1:2], align = 'hv', nrow = 1)
botrow <- plot_grid(plotlist = histlist[3:4], align = 'hv', nrow = 1)
pdf(file = paste(output_dir, "test_finger_limma/p_value_histograms.pdf", sep = ""))
plot_grid(toprow, botrow, ncol = 1)
dev.off()

### Create volcano plots ----
volcanos <- list()
for (i in 1:2) {
  volcanos[[i]] <- EnhancedVolcano(tables[[i]],
                                   lab = tables[[i]]$Symbol,
                                   x = "logFC",
                                   y = "P.Value",
                                   title = (paste("Test_Finger:", names(y)[i])), titleLabSize = 10,
                                   legendLabSize = 7, legendIconSize = 4,
                                   subtitle = "",
                                   labSize = (4.0),
                                   colAlpha = 1,
                                   pointSize = 1.0, 
                                   pCutoff = 1e-2,
                                   boxedLabels = FALSE,
                                   drawConnectors = TRUE)
}
pdf(file = paste(output_dir, "test_finger_limma/volcanos.pdf", sep = ""), width = 9)
plot_grid(plotlist = volcanos, align = 'hv', nrow = 1)
dev.off()

# Analysis 3: Test_Set_Vein (p2v) ----
## Step 1: Read in raw expression data ----
p2v[1:5,1:5]
#                                    title geo_accession subject age time
# GSM745049  009-282-001-TP1 [5356583035D]     GSM745049       1  49   -7
# GSM745050  009-282-001-TP2 [5356583041E]     GSM745050       1  49    0
# GSM745051  009-282-001-TP7 [5356583031A]     GSM745051       1  49  0.5
# GSM745052  009-282-001-TP9 [5356583049D]     GSM745052       1  49    1
# GSM745053 009-282-001-TP12 [5356583040D]     GSM745053       1  49    3
levels(factor(p2v$time))
# [1] "-7"  "0"   "0.5" "1"   "10"  "28"  "3"   "7"  
p2v$time <- factor(p2v$time, levels = c("-7","0","0.5","1","3","7","10","28"), ordered = T)
## This is consistent with supplementary data from Obermoser et al.
dim(p2v)
# [1] 143  12
table(p2v$time, p2v$subject)
#     1 10 11 12 13 14 15 18 20 21 22 24 26 3 4 6 7 8
# -7  1  1  1  1  1  1  1  1  1  1  1  1  1 1 1 1 1 1
# 0   1  1  1  1  1  1  1  1  1  1  1  1  1 1 1 1 1 1
# 0.5 1  1  1  1  1  1  1  1  1  1  1  1  1 1 1 1 1 1
# 1   1  1  1  1  1  1  1  1  1  1  1  1  1 1 1 1 1 1
# 3   1  1  1  1  1  1  1  1  1  1  1  1  1 1 1 1 1 1
# 7   1  1  1  1  1  1  1  1  1  1  1  1  1 1 1 1 1 1
# 10  1  1  1  1  1  1  1  1  1  1  1  1  1 1 1 1 1 1
# 28  1  1  0  1  1  1  1  1  1  1  1  1  1 1 1 1 1 1
table(p2v$sex, p2v$subject) ## Sanity check.
#   1 10 11 12 13 14 15 18 20 21 22 24 26 3 4 6 7 8
# F 8  0  0  0  8  0  0  8  0  0  8  0  8 8 8 8 8 0
# M 0  8  7  8  0  8  8  0  8  8  0  8  0 0 0 0 0 8
table(p2v$treat, p2v$subject)
#         1 10 11 12 13 14 15 18 20 21 22 24 26 3 4 6 7 8
# FLUZONE 8  8  0  0  0  0  0  8  0  0  8  8  0 0 0 0 0 8
# PNEUM   0  0  7  0  0  0  8  0  0  8  0  0  8 0 8 0 8 0
# Saline  0  0  0  8  8  8  0  0  8  0  0  0  0 8 0 8 0 0

## Read in files
temp <- files[str_extract(files, "GSM......") %in% p2v$geo_accession]
length(temp) # [1] 143
x <- lapply(
  temp,
  function(x) read.delim(file = x))
x <- do.call(cbind, x)
x[1:4,1:4]
#              X X5356583035_D.AVG_Signal X5356583035_D.Detection_Pval          X.1
# 1 ILMN_1885639                   5.6323                       0.1700 ILMN_1885639
# 2 ILMN_2209417                  -9.6776                       0.8972 ILMN_2209417
# 3 ILMN_1765401                 102.6921                       0.0013 ILMN_1765401
# 4 ILMN_1652171                   2.8823                       0.2622 ILMN_1652171
dim(x)
# [1] 48803   429
rownames(x) <- x[,1]
x <- as.data.frame(x)
detectionpvalues <- x[,grep("Detection",colnames(x))]
x <- x[,grep("AVG",colnames(x))]
x[1:2,1:4]
#              X5356583035_D.AVG_Signal X5356583041_E.AVG_Signal X5356583031_A.AVG_Signal X5356583049_D.AVG_Signal
# ILMN_1885639                   5.6323                  -5.9032                  -4.7824                   2.0096
# ILMN_2209417                  -9.6776                  -0.5463                  -2.5845                   9.7934

## Step 2: Create an EList object ----
data <- new('EListRaw')
data@.Data[[1]] <- 'illumina'
# data@.Data[[2]] <- targetinfo
data@.Data[[3]] <- BGX
data@.Data[[4]] <- x
data@.Data[[5]] <- NULL
data$E <- x
# data$targets <- targetinfo
data$genes <- BGX
data$other$Detection <- detectionpvalues

### Background correction and normalization ----
data_norm <- neqc(data)
table(data_norm$genes$controls$Reporter_Group_Name)
# biotin            cy3_hyb       housekeeping           labeling low_stringency_hyb 
# 2                  6                  7                  2                  8 
# negative 
# 759 


## Remove unexpressed probes and controls
controls <- data_norm$genes$controls$Probe_Id
expressed <- rowSums(data_norm$other$Detection < 0.05) >= 3
data_norm$E <- data_norm$E[expressed,]
data_norm$E <- data_norm$E[-which(rownames(data_norm$E) %in% controls),]
dim(data_norm$E)
# [1] 27531   143

## Step 3: PCA and MDS analysis ----
# are pheno and matrix aligned?
identical(substr(p2v$array, start = 1, stop = 10), str_extract(colnames(data_norm$E),"\\d........."))
# [1] TRUE
colnames(data_norm$E) <- rownames(p2v)
### PCAtools ----
p <- pca(data_norm$E, metadata = p2v, center = T, scale = T)

scree <- screeplot(p, getComponents(p, c(1:10)))
ptime <- biplot(p, 
                showLoadings = F,
                lab = NULL, 
                colby = "time", 
                legendPosition = 'right', legendLabSize = 12, legendIconSize = 3,
                title = 'Obermoser QC',
                subtitle = 'p2v samples')
psex <- biplot(p, 
               showLoadings = F,
               lab = NULL, 
               colby = "sex", 
               legendPosition = 'right', legendLabSize = 9, legendIconSize = 3)
psub <- biplot(p, 
               showLoadings = F,
               lab = NULL, 
               colby = "subject",
               legendPosition = 'right', legendLabSize = 9, legendIconSize = 3)
ptre <-  biplot(p, 
                showLoadings = F,
                lab = NULL, 
                colby = "treat",
                legendPosition = 'right', legendLabSize = 9, legendIconSize = 3)
peigencor <- eigencorplot(p,
                          components = getComponents(p, 1:10),
                          metavars = c("time", 'sex', 'race', 'ethnicity', 'treat'),
                          cexCorval = 0.5,
                          fontCorval = 1, 
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
botrow <- plot_grid(ptre, scree, peigencor, align = 'h', nrow = 1)
pdf(file = paste(output_dir, "pca_p2v_QC.pdf", sep = ""), height = 12, width = 16)
plot_grid(toprow, botrow, ncol = 1, align = 'v')
dev.off()
rownames(p$loadings) <- make.names(mapIds(x = illuminaHumanv3.db, keys = rownames(p$loadings), column = "SYMBOL", keytype = "PROBEID"), unique = T)
pdf(file = paste(output_dir, "ethnicity_p2v_pca.pdf", sep = ""))
biplot(p, 
       showLoadings = T,
       lab = NULL, 
       colby = "ethnicity", 
       legendPosition = 'right', legendLabSize = 10, legendIconSize = 3,
       title = 'Ethnicity of samples',
       subtitle = 'p2v samples')
dev.off()

### MDS ----
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colors <- sample(col_vector, length(unique(p2v$subject)))
new_colors <- c()
for (i in 1:length(p2v$subject)) {
  new_colors[i] <- switch(p2v$subject[i],
                          "18" = colors[1], "10" = colors[2], "24" = colors[3], "1" = colors[4], 
                          "22" = colors[5], "8" = colors[6], "6" = colors[7], "13" = colors[8],
                          "3" = colors[9], "14" = colors[10], "20" = colors[11], "7" = colors[12],
                          "4" = colors[13], "15" = colors[14], "21" = colors[15], "26" = colors[16],
                          "11" = colors[17], "12" = colors[18]
  )
}
pdf(file = paste(output_dir, "mds_p2v_qc.pdf", sep = ""))
limma::plotMDS(x = data_norm$E, 
               cex = 1, 
               var.explained = T, 
               col = new_colors,
               labels = paste(p2v$subject, p2v$sex, p2v$treat, sep = ", "))
title("MDS Plot - test set vein")
dev.off()

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
## Step 4: limma analysis ----
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

id <- factor(p2v$subject)
time <- factor(p2v$time)
sex <- factor(p2v$sex)
vax <- factor(p2v$treat)
eth <- factor(p2v$ethnicity)
race <- factor(p2v$race)

df <- data.frame(data_norm$E)

TS <- factor(paste(sex,time,vax,race,eth, sep = "."))
design <- model.matrix(~0+TS, df)
colnames(design) <- make.names(levels(TS))
cor <- duplicateCorrelation(df, design, block=id)
cor$consensus.correlation
# [1] 0.1599962

fit <- lmFit(object=df, design=design, block=id, correlation=cor$consensus.correlation)

# contrasts
my_contrasts <- makeContrasts(# Time effects
  BaselineVSzeroAll = (F.0.FLUZONE.Asian.NA + F.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.0.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.0.Saline.African.American.NA + F.0.Saline.Caucasian.Hispanic + F.0.Saline.Caucasian.Non.Hispanic.or.Latino + M.0.FLUZONE.Asian.NA + M.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.0.PNEUM.African.American.NA + M.0.PNEUM.Caucasian.Hispanic + M.0.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.0.Saline.Caucasian.Non.Hispanic.or.Latino)/12 -
    (F..7.FLUZONE.Asian.NA + F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.African.American.NA + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Asian.NA + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.African.American.NA + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/12,
  BaselineVSpoint5All = (F.0.5.FLUZONE.Asian.NA + F.0.5.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.0.5.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.0.5.Saline.African.American.NA + F.0.5.Saline.Caucasian.Hispanic + F.0.5.Saline.Caucasian.Non.Hispanic.or.Latino + M.0.5.FLUZONE.Asian.NA + M.0.5.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.0.5.PNEUM.African.American.NA + M.0.5.PNEUM.Caucasian.Hispanic + M.0.5.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.0.5.Saline.Caucasian.Non.Hispanic.or.Latino)/12 -
    (F..7.FLUZONE.Asian.NA + F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.African.American.NA + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Asian.NA + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.African.American.NA + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/12,
  BaselineVSoneAll = (F.1.FLUZONE.Asian.NA + F.1.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.1.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.1.Saline.African.American.NA + F.1.Saline.Caucasian.Hispanic + F.1.Saline.Caucasian.Non.Hispanic.or.Latino + M.1.FLUZONE.Asian.NA + M.1.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.1.PNEUM.African.American.NA + M.1.PNEUM.Caucasian.Hispanic + M.1.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.1.Saline.Caucasian.Non.Hispanic.or.Latino)/12 -
    (F..7.FLUZONE.Asian.NA + F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.African.American.NA + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Asian.NA + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.African.American.NA + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/12,
  BaselineVSthreeAll = (F.3.FLUZONE.Asian.NA + F.3.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.3.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.3.Saline.African.American.NA + F.3.Saline.Caucasian.Hispanic + F.3.Saline.Caucasian.Non.Hispanic.or.Latino + M.3.FLUZONE.Asian.NA + M.3.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.3.PNEUM.African.American.NA + M.3.PNEUM.Caucasian.Hispanic + M.3.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.3.Saline.Caucasian.Non.Hispanic.or.Latino)/12 -
    (F..7.FLUZONE.Asian.NA + F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.African.American.NA + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Asian.NA + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.African.American.NA + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/12,
  BaselineVSsevenAll = (F.7.FLUZONE.Asian.NA + F.7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.7.Saline.African.American.NA + F.7.Saline.Caucasian.Hispanic + F.7.Saline.Caucasian.Non.Hispanic.or.Latino + M.7.FLUZONE.Asian.NA + M.7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.7.PNEUM.African.American.NA + M.7.PNEUM.Caucasian.Hispanic + M.7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.7.Saline.Caucasian.Non.Hispanic.or.Latino)/12 -
    (F..7.FLUZONE.Asian.NA + F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.African.American.NA + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Asian.NA + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.African.American.NA + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/12,
  BaselineVStenAll = (F.10.FLUZONE.Asian.NA + F.10.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.10.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.10.Saline.African.American.NA + F.10.Saline.Caucasian.Hispanic + F.10.Saline.Caucasian.Non.Hispanic.or.Latino + M.10.FLUZONE.Asian.NA + M.10.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.10.PNEUM.Caucasian.Hispanic + M.10.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.10.Saline.Caucasian.Non.Hispanic.or.Latino + M.10.PNEUM.African.American.NA)/12 -
    (F..7.FLUZONE.Asian.NA + F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.African.American.NA + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Asian.NA + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.African.American.NA + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/12,
  BaselineVStwentyEightAll = (F.28.FLUZONE.Asian.NA + F.28.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.28.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.28.Saline.African.American.NA + F.28.Saline.Caucasian.Hispanic + F.28.Saline.Caucasian.Non.Hispanic.or.Latino + M.28.FLUZONE.Asian.NA + M.28.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.28.PNEUM.Caucasian.Hispanic + M.28.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.28.Saline.Caucasian.Non.Hispanic.or.Latino)/11 -
    (F..7.FLUZONE.Asian.NA + F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.African.American.NA + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Asian.NA + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.African.American.NA + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/12,
  # Saline times
  SalineBaseVSzero = (F.0.Saline.African.American.NA + F.0.Saline.Caucasian.Hispanic + F.0.Saline.Caucasian.Non.Hispanic.or.Latino + M.0.Saline.Caucasian.Non.Hispanic.or.Latino)/4 -
    (F..7.Saline.African.American.NA + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/4,
  SalineBaseVSpoint5 = (F.0.5.Saline.African.American.NA + F.0.5.Saline.Caucasian.Hispanic + F.0.5.Saline.Caucasian.Non.Hispanic.or.Latino + M.0.5.Saline.Caucasian.Non.Hispanic.or.Latino)/4 -
    (F..7.Saline.African.American.NA + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/4,
  SalineBaseVSone = (F.1.Saline.African.American.NA + F.1.Saline.Caucasian.Hispanic + F.1.Saline.Caucasian.Non.Hispanic.or.Latino + M.1.Saline.Caucasian.Non.Hispanic.or.Latino)/4 -
    (F..7.Saline.African.American.NA + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/4,
  SalineBaseVSthree = (F.3.Saline.African.American.NA + F.3.Saline.Caucasian.Hispanic + F.3.Saline.Caucasian.Non.Hispanic.or.Latino + M.3.Saline.Caucasian.Non.Hispanic.or.Latino)/4 -
    (F..7.Saline.African.American.NA + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/4,
  SalineBaseVSseven = (F.7.Saline.African.American.NA + F.7.Saline.Caucasian.Hispanic + F.7.Saline.Caucasian.Non.Hispanic.or.Latino + M.7.Saline.Caucasian.Non.Hispanic.or.Latino)/4 -
    (F..7.Saline.African.American.NA + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/4,
  SalineBaseVSten = (F.10.Saline.African.American.NA + F.10.Saline.Caucasian.Hispanic + F.10.Saline.Caucasian.Non.Hispanic.or.Latino + M.10.Saline.Caucasian.Non.Hispanic.or.Latino)/4 -
    (F..7.Saline.African.American.NA + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/4,
  SalineBaseVStwentyEight = (F.28.Saline.African.American.NA + F.28.Saline.Caucasian.Hispanic + F.28.Saline.Caucasian.Non.Hispanic.or.Latino + M.28.Saline.Caucasian.Non.Hispanic.or.Latino)/4 -
    (F..7.Saline.African.American.NA + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/4,
  # Pneum times
  PNEUMBaseVSzero = (F.0.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.0.PNEUM.African.American.NA + M.0.PNEUM.Caucasian.Hispanic + M.0.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4 -
    (F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.African.American.NA + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4,
  PNEUMBaseVSpoint5 = (F.0.5.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.0.5.PNEUM.African.American.NA + M.0.5.PNEUM.Caucasian.Hispanic + M.0.5.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4 -
    (F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.African.American.NA + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4,
  PNEUMBaseVSone = (F.1.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.1.PNEUM.African.American.NA + M.1.PNEUM.Caucasian.Hispanic + M.1.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4 -
    (F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.African.American.NA + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4,
  PNEUMBaseVSthree = (F.3.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.3.PNEUM.African.American.NA + M.3.PNEUM.Caucasian.Hispanic + M.3.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4 -
    (F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.African.American.NA + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4,
  PNEUMBaseVSseven = (F.7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.7.PNEUM.African.American.NA + M.7.PNEUM.Caucasian.Hispanic + M.7.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4 -
    (F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.African.American.NA + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4,
  PNEUMBaseVSten = (F.10.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.10.PNEUM.Caucasian.Hispanic + M.10.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.10.PNEUM.African.American.NA)/4 -
    (F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.African.American.NA + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4,
  PNEUMBaseVStwentyEight = (F.28.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.28.PNEUM.Caucasian.Hispanic + M.28.PNEUM.Caucasian.Non.Hispanic.or.Latino)/3 -
    (F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.African.American.NA + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4,
  # Flu times
  FLUBaseVSzero = (F.0.FLUZONE.Asian.NA + F.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.0.FLUZONE.Asian.NA + M.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4 -
    (F..7.FLUZONE.Asian.NA + F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Asian.NA + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4,
  FLUBaseVSpoint5 = (F.0.5.FLUZONE.Asian.NA + F.0.5.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.0.5.FLUZONE.Asian.NA + M.0.5.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4 -
    (F..7.FLUZONE.Asian.NA + F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Asian.NA + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4,
  FLUBaseVSone = (F.1.FLUZONE.Asian.NA + F.1.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.1.FLUZONE.Asian.NA + M.1.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4 -
    (F..7.FLUZONE.Asian.NA + F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Asian.NA + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4,
  FLUBaseVSthree = (F.3.FLUZONE.Asian.NA + F.3.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.3.FLUZONE.Asian.NA + M.3.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4 -
    (F..7.FLUZONE.Asian.NA + F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Asian.NA + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4,
  FLUBaseVSseven = (F.7.FLUZONE.Asian.NA + F.7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.7.FLUZONE.Asian.NA + M.7.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4 -
    (F..7.FLUZONE.Asian.NA + F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Asian.NA + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4,
  FLUBaseVSten = (F.10.FLUZONE.Asian.NA + F.10.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.10.FLUZONE.Asian.NA + M.10.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4 -
    (F..7.FLUZONE.Asian.NA + F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Asian.NA + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4,
  FLUBaseVStwentyEight = (F.28.FLUZONE.Asian.NA + F.28.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.28.FLUZONE.Asian.NA + M.28.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4 -
    (F..7.FLUZONE.Asian.NA + F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Asian.NA + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4,
  # Time-sex interaction Male - Female
  IMaleBaseVSzero = ((M.0.FLUZONE.Asian.NA + M.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.0.PNEUM.African.American.NA + M.0.PNEUM.Caucasian.Hispanic + M.0.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.0.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
                       (M..7.FLUZONE.Asian.NA + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.African.American.NA + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/6) -
    ((F.0.FLUZONE.Asian.NA + F.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.0.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.0.Saline.African.American.NA + F.0.Saline.Caucasian.Hispanic + F.0.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..7.FLUZONE.Asian.NA + F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.African.American.NA + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino)/6),
  IMaleBaseVSpoint5 = ((M.0.5.FLUZONE.Asian.NA + M.0.5.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.0.5.PNEUM.African.American.NA + M.0.5.PNEUM.Caucasian.Hispanic + M.0.5.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.0.5.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
                       (M..7.FLUZONE.Asian.NA + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.African.American.NA + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/6) -
    ((F.0.5.FLUZONE.Asian.NA + F.0.5.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.0.5.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.0.5.Saline.African.American.NA + F.0.5.Saline.Caucasian.Hispanic + F.0.5.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..7.FLUZONE.Asian.NA + F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.African.American.NA + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino)/6),
  IMaleBaseVSone = ((M.1.FLUZONE.Asian.NA + M.1.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.1.PNEUM.African.American.NA + M.1.PNEUM.Caucasian.Hispanic + M.1.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.1.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
                         (M..7.FLUZONE.Asian.NA + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.African.American.NA + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/6) -
    ((F.1.FLUZONE.Asian.NA + F.1.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.1.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.1.Saline.African.American.NA + F.1.Saline.Caucasian.Hispanic + F.1.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..7.FLUZONE.Asian.NA + F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.African.American.NA + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino)/6),
  IMaleBaseVSthree = ((M.3.FLUZONE.Asian.NA + M.3.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.3.PNEUM.African.American.NA + M.3.PNEUM.Caucasian.Hispanic + M.3.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.3.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
                      (M..7.FLUZONE.Asian.NA + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.African.American.NA + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/6) -
    ((F.3.FLUZONE.Asian.NA + F.3.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.3.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.3.Saline.African.American.NA + F.3.Saline.Caucasian.Hispanic + F.3.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..7.FLUZONE.Asian.NA + F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.African.American.NA + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino)/6),
  IMaleBaseVSseven = ((M.7.FLUZONE.Asian.NA + M.7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.7.PNEUM.African.American.NA + M.7.PNEUM.Caucasian.Hispanic + M.7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.7.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
                        (M..7.FLUZONE.Asian.NA + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.African.American.NA + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/6) -
    ((F.7.FLUZONE.Asian.NA + F.7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.7.Saline.African.American.NA + F.7.Saline.Caucasian.Hispanic + F.7.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..7.FLUZONE.Asian.NA + F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.African.American.NA + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino)/6),
  IMaleBaseVSten = ((M.10.FLUZONE.Asian.NA + M.10.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.10.PNEUM.African.American.NA + M.10.PNEUM.Caucasian.Hispanic + M.10.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.10.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
                        (M..7.FLUZONE.Asian.NA + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.African.American.NA + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/6) -
    ((F.10.FLUZONE.Asian.NA + F.10.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.10.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.10.Saline.African.American.NA + F.10.Saline.Caucasian.Hispanic + F.10.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..7.FLUZONE.Asian.NA + F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.African.American.NA + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino)/6),
  IMaleBaseVStwentyEight = ((M.28.FLUZONE.Asian.NA + M.28.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.28.PNEUM.Caucasian.Hispanic + M.28.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.28.Saline.Caucasian.Non.Hispanic.or.Latino)/5-
                      (M..7.FLUZONE.Asian.NA + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.African.American.NA + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/6) -
    ((F.28.FLUZONE.Asian.NA + F.28.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.28.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.28.Saline.African.American.NA + F.28.Saline.Caucasian.Hispanic + F.28.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..7.FLUZONE.Asian.NA + F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.African.American.NA + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino)/6),
  # Caucasian only
  CaucasianBaseVSzero = (F.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.0.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.0.Saline.Caucasian.Hispanic + F.0.Saline.Caucasian.Non.Hispanic.or.Latino + M.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.0.PNEUM.Caucasian.Hispanic + M.0.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.0.Saline.Caucasian.Non.Hispanic.or.Latino)/8-
    (F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/8,
  CaucasianBaseVSpoint5 = (F.0.5.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.0.5.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.0.5.Saline.Caucasian.Hispanic + F.0.5.Saline.Caucasian.Non.Hispanic.or.Latino + M.0.5.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.0.5.PNEUM.Caucasian.Hispanic + M.0.5.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.0.5.Saline.Caucasian.Non.Hispanic.or.Latino)/8-
    (F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/8,
  CaucasianBaseVSone = (F.1.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.1.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.1.Saline.Caucasian.Hispanic + F.1.Saline.Caucasian.Non.Hispanic.or.Latino + M.1.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.1.PNEUM.Caucasian.Hispanic + M.1.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.1.Saline.Caucasian.Non.Hispanic.or.Latino)/8-
    (F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/8,
  CaucasianBaseVSthree = (F.3.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.3.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.3.Saline.Caucasian.Hispanic + F.3.Saline.Caucasian.Non.Hispanic.or.Latino + M.3.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.3.PNEUM.Caucasian.Hispanic + M.3.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.3.Saline.Caucasian.Non.Hispanic.or.Latino)/8-
    (F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/8,
  CaucasianBaseVSseven = (F.7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.7.Saline.Caucasian.Hispanic + F.7.Saline.Caucasian.Non.Hispanic.or.Latino + M.7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.7.PNEUM.Caucasian.Hispanic + M.7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.7.Saline.Caucasian.Non.Hispanic.or.Latino)/8-
    (F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/8,
  CaucasianBaseVSten = (F.10.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.10.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.10.Saline.Caucasian.Hispanic + F.10.Saline.Caucasian.Non.Hispanic.or.Latino + M.10.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.10.PNEUM.Caucasian.Hispanic + M.10.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.10.Saline.Caucasian.Non.Hispanic.or.Latino)/8-
    (F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/8,
  CaucasianBaseVStwentyEight = (F.28.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.28.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.28.Saline.Caucasian.Hispanic + F.28.Saline.Caucasian.Non.Hispanic.or.Latino + M.28.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.28.PNEUM.Caucasian.Hispanic + M.28.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.28.Saline.Caucasian.Non.Hispanic.or.Latino)/8-
    (F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.Caucasian.Hispanic + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/8,
  # Asian only
  AsianBaseVSzero = (F.0.FLUZONE.Asian.NA + M.0.FLUZONE.Asian.NA)/2 - (F..7.FLUZONE.Asian.NA + M..7.FLUZONE.Asian.NA)/2,
  AsianBaseVSpoint5 = (F.0.5.FLUZONE.Asian.NA + M.0.5.FLUZONE.Asian.NA)/2 - (F..7.FLUZONE.Asian.NA + M..7.FLUZONE.Asian.NA)/2,
  AsianBaseVSone = (F.1.FLUZONE.Asian.NA + M.1.FLUZONE.Asian.NA)/2 - (F..7.FLUZONE.Asian.NA + M..7.FLUZONE.Asian.NA)/2,
  AsianBaseVSthree = (F.3.FLUZONE.Asian.NA + M.3.FLUZONE.Asian.NA)/2 - (F..7.FLUZONE.Asian.NA + M..7.FLUZONE.Asian.NA)/2,
  AsianBaseVSseven = (F.7.FLUZONE.Asian.NA + M.7.FLUZONE.Asian.NA)/2 - (F..7.FLUZONE.Asian.NA + M..7.FLUZONE.Asian.NA)/2,
  AsianBaseVSten = (F.10.FLUZONE.Asian.NA + M.10.FLUZONE.Asian.NA)/2 - (F..7.FLUZONE.Asian.NA + M..7.FLUZONE.Asian.NA)/2,
  AsianBaseVStwentyEight = (F.28.FLUZONE.Asian.NA + M.28.FLUZONE.Asian.NA)/2 - (F..7.FLUZONE.Asian.NA + M..7.FLUZONE.Asian.NA)/2,
  # African American only
  AfricanAmerBaseVSzero = (F.0.Saline.African.American.NA + M.0.PNEUM.African.American.NA)/2 - (F..7.Saline.African.American.NA + M..7.PNEUM.African.American.NA)/2,
  AfricanAmerBaseVSpoint5 = (F.0.5.Saline.African.American.NA + M.0.5.PNEUM.African.American.NA)/2 - (F..7.Saline.African.American.NA + M..7.PNEUM.African.American.NA)/2,
  AfricanAmerBaseVSone = (F.1.Saline.African.American.NA + M.1.PNEUM.African.American.NA)/2 - (F..7.Saline.African.American.NA + M..7.PNEUM.African.American.NA)/2,
  AfricanAmerBaseVSthree = (F.3.Saline.African.American.NA + M.3.PNEUM.African.American.NA)/2 - (F..7.Saline.African.American.NA + M..7.PNEUM.African.American.NA)/2,
  AfricanAmerBaseVSseven = (F.7.Saline.African.American.NA + M.7.PNEUM.African.American.NA)/2 - (F..7.Saline.African.American.NA + M..7.PNEUM.African.American.NA)/2,
  AfricanAmerBaseVSten = (F.10.Saline.African.American.NA + M.10.PNEUM.African.American.NA)/2 - (F..7.Saline.African.American.NA + M..7.PNEUM.African.American.NA)/2,
  AfricanAmerBaseVStwentyEight = (F.28.Saline.African.American.NA) - (F..7.Saline.African.American.NA + M..7.PNEUM.African.American.NA)/2,
  # Ethnicity-time interaction Hispanic - not Hispanic
  IHispanicBaseVSzero = ((F.0.Saline.Caucasian.Hispanic + M.0.PNEUM.Caucasian.Hispanic)/2 - (F..7.Saline.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Hispanic)/2) -
    ((F.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.0.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.0.Saline.Caucasian.Non.Hispanic.or.Latino + M.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.0.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.0.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/6),
  IHispanicBaseVSpoint5 = ((F.0.5.Saline.Caucasian.Hispanic + M.0.5.PNEUM.Caucasian.Hispanic)/2 - (F..7.Saline.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Hispanic)/2) -
    ((F.0.5.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.0.5.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.0.5.Saline.Caucasian.Non.Hispanic.or.Latino + M.0.5.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.0.5.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.0.5.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/6),
  IHispanicBaseVSone = ((F.1.Saline.Caucasian.Hispanic + M.1.PNEUM.Caucasian.Hispanic)/2 - (F..7.Saline.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Hispanic)/2) -
    ((F.1.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.1.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.1.Saline.Caucasian.Non.Hispanic.or.Latino + M.1.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.1.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.1.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/6),
  IHispanicBaseVSthree = ((F.3.Saline.Caucasian.Hispanic + M.3.PNEUM.Caucasian.Hispanic)/2 - (F..7.Saline.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Hispanic)/2) -
    ((F.3.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.3.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.3.Saline.Caucasian.Non.Hispanic.or.Latino + M.3.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.3.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.3.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/6),
  IHispanicBaseVSseven = ((F.7.Saline.Caucasian.Hispanic + M.7.PNEUM.Caucasian.Hispanic)/2 - (F..7.Saline.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Hispanic)/2) -
    ((F.7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.7.Saline.Caucasian.Non.Hispanic.or.Latino + M.7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.7.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/6),
  IHispanicBaseVSten = ((F.10.Saline.Caucasian.Hispanic + M.10.PNEUM.Caucasian.Hispanic)/2 - (F..7.Saline.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Hispanic)/2) -
    ((F.10.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.10.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.10.Saline.Caucasian.Non.Hispanic.or.Latino + M.10.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.10.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.10.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/6),
  IHispanicBaseVStwentyEight = ((F.28.Saline.Caucasian.Hispanic + M.28.PNEUM.Caucasian.Hispanic)/2 - (F..7.Saline.Caucasian.Hispanic + M..7.PNEUM.Caucasian.Hispanic)/2) -
    ((F.28.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.28.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.28.Saline.Caucasian.Non.Hispanic.or.Latino + M.28.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.28.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.28.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..7.Saline.Caucasian.Non.Hispanic.or.Latino + M..7.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..7.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..7.Saline.Caucasian.Non.Hispanic.or.Latino)/6),
  levels = colnames(design))

fit2 <- contrasts.fit(fit, my_contrasts)
x <- eBayes(fit2)

### Results tables ----

y <- topTable(x, number = Inf, adjust.method = "fdr")
y$Symbol <- mapIds(x = illuminaHumanv3.db, keys = row.names(y), column = "SYMBOL", keytype = "PROBEID")
dir.create(paste(output_dir, "test_vein_limma", sep = ""))
write.csv(y, file = ,paste(output_dir, "test_vein_limma/limma_F.csv", sep = ""), row.names = TRUE)

dt <- decideTests(x)
summary(dt)
#        BaselineVSzeroAll BaselineVSpoint5All BaselineVSoneAll BaselineVSthreeAll BaselineVSsevenAll BaselineVStenAll BaselineVStwentyEightAll
# Down                   0                   4                0                  0                  6                0                        0
# NotSig             27530               27525            27530              27531              27450            27531                    27531
# Up                     1                   2                1                  0                 75                0                        0
dt <- decideTests(x, adjust.method = 'none')
summary(dt)
#        BaselineVSzeroAll BaselineVSpoint5All BaselineVSoneAll BaselineVSthreeAll BaselineVSsevenAll BaselineVStenAll BaselineVStwentyEightAll
# Down                 711                 713              463                832               1075             1124                      945
# NotSig             26355               25965            26581              26014              25513            25373                    25753
# Up                   465                 853              487                685                943             1034                      833

names(y)
tables <- list()
for (i in 1:63) {
  tables[[i]] <- topTable(x, coef = i, number = Inf, adjust.method = "fdr")
  tables[[i]]$Symbol <- mapIds(x = illuminaHumanv3.db, keys = row.names(tables[[i]]), column = "SYMBOL", keytype = "PROBEID")
  write.csv(tables[[i]], file = paste(output_dir, "test_vein_limma/",
                                      names(y)[i],".csv",sep = ""), row.names = TRUE)
}

### Plot p-value histograms ----
histlist <- list()
for (i in 1:7) {
  histlist[[i]] <- ggplot(data = tables[[i]], aes(x = adj.P.Val)) +
    geom_histogram(fill = "lightblue", color = "black") +
    ggtitle(paste("Adj. p histogram", names(y)[i])) +
    theme_cowplot() +
    theme(title = element_text(size = 7.5))
  names(histlist)[i] <- names(y)[i]
}
pdf(file = paste(output_dir, "test_vein_limma/adj_p_value_histograms.pdf", sep = ""), width = 12)
plot_grid(plotlist = histlist, align = 'hv', nrow = 2)
dev.off()
for (i in 1:7) {
  histlist[[i]] <- ggplot(data = tables[[i]], aes(x = P.Value)) +
    geom_histogram(fill = "lightblue", color = "black") +
    ggtitle(paste("P-value histogram", names(y)[i])) +
    theme_cowplot() +
    theme(title = element_text(size = 7.5))
  names(histlist)[i] <- names(y)[i]
}
pdf(file = paste(output_dir, "test_vein_limma/p_value_histograms.pdf", sep = ""), width = 12)
plot_grid(plotlist = histlist, align = 'hv', nrow = 2)
dev.off()

### Create volcano plots ----
volcanos <- list()
for (i in 1:7) {
  volcanos[[i]] <- EnhancedVolcano(tables[[i]],
                                   lab = tables[[i]]$Symbol,
                                   x = "logFC",
                                   y = "P.Value",
                                   title = (paste("Test_Vein:", names(y)[i])),
                                   titleLabSize = 8,
                                   subtitle = "",
                                   labSize = (4.0),
                                   colAlpha = 1,
                                   pointSize = 1.0, 
                                   pCutoff = 1e-2,
                                   boxedLabels = FALSE,
                                   drawConnectors = TRUE,
                                   legendLabSize = 5, legendIconSize = 2,
                                   axisLabSize = 5,
                                   caption = "")
  names(volcanos)[i] <- names(y)[i]
}
pdf(file = paste(output_dir, "test_vein_limma/volcanos.pdf", sep = ""), width = 12, height = 10)
plot_grid(plotlist = volcanos, align = 'hv', nrow = 2)
dev.off()

# Analysis 4: Training_Set_Finger (p2f) ----
## Step 1: Read in raw expression data ----
p2f[1:5,1:5]
#                                   title geo_accession subject age time
# GSM745214 009-282-001-TP1 [5307141014E]     GSM745214       1  49 -168
# GSM745215 009-282-001-TP2 [5307141024E]     GSM745215       1  49    0
# GSM745216 009-282-001-TP3 [5307141024K]     GSM745216       1  49  1.5
# GSM745217 009-282-001-TP4 [5307141024F]     GSM745217       1  49    3
# GSM745218 009-282-001-TP5 [5322244022H]     GSM745218       1  49    6
levels(factor(p2f$time))
# [1] "-168" "0"    "1.5"  "12"   "15"   "24"   "3"    "36"   "48"   "6"    "9"  
p2f$time <- factor(p2f$time, levels = c("-168","0","1.5","3","6","9","12","15","24","36","48"), ordered = T)
## This is recorded in units of hours.
dim(p2f)
# [1] 185  12
table(p2f$time, p2f$subject)
#      1 10 11 13 14 15 18 20 21 22 24 26 3 4 6 7 8
# -168 1  1  1  1  1  1  1  0  1  1  1  1 1 1 1 1 1 # one missing sample
# 0    1  1  1  1  1  1  1  1  1  1  1  1 1 1 1 1 1
# 1.5  1  1  1  1  1  1  1  1  1  1  1  1 1 1 1 1 1
# 3    1  1  1  1  1  1  1  1  1  1  1  1 1 1 1 1 1
# 6    1  1  1  1  1  1  1  1  1  1  1  1 1 1 1 1 1
# 9    1  0  1  1  1  1  1  1  1  1  1  1 1 1 1 1 1 # one missing sample
# 12   1  1  1  1  1  1  1  1  1  1  1  1 1 1 1 1 1
# 15   1  1  1  1  1  1  1  1  1  1  1  1 1 1 1 1 1
# 24   1  1  1  1  1  1  1  1  1  1  1  1 1 1 1 1 1
# 36   1  1  1  1  1  1  1  1  1  1  1  1 1 1 1 1 1
# 48   1  1  1  1  1  1  1  1  1  1  1  1 1 1 1 1 1
p2f[which(p2f$subject == 10),] # Non-hispanic, caucasian male treated w/ FLUZONE
p2f[which(p2f$subject == 20),] # Non-hispanic, caucasian male treated w/ Saline
table(p2f$sex, p2f$subject) ## Sanity check.
#    1 10 11 13 14 15 18 20 21 22 24 26  3  4  6  7  8
# F 11  0  0 11  0  0 11  0  0 11  0 11 11 11 11 11  0
# M  0 10 11  0 11 11  0 10 11  0 11  0  0  0  0  0 11
table(p2f$treat, p2f$subject)
#          1 10 11 13 14 15 18 20 21 22 24 26  3  4  6  7  8
# FLUZONE 11 10  0  0  0  0 11  0  0 11 11  0  0  0  0  0 11
# PNEUM    0  0 11  0  0 11  0  0 11  0  0 11  0 11  0 11  0
# Saline   0  0  0 11 11  0  0 10  0  0  0  0 11  0 11  0  0

## Load in files
temp <- files[str_extract(files, "GSM......") %in% p2f$geo_accession]
length(temp) # [1] 185
x <- lapply(
  temp,
  function(x) read.delim(file = x))
x <- do.call(cbind, x)
x[1:4,1:4]
#              X X5307141014_E.AVG_Signal X5307141014_E.Detection_Pval          X.1
# 1 ILMN_1885639                  -0.1356                       0.3544 ILMN_1885639
# 2 ILMN_2209417                  18.5837                       0.1054 ILMN_2209417
# 3 ILMN_1765401                  44.3429                       0.0237 ILMN_1765401
# 4 ILMN_1652171                  -5.8566                       0.5455 ILMN_1652171
dim(x)
# [1] 48803   555
rownames(x) <- x[,1]
x <- as.data.frame(x)
detectionpvalues <- x[,grep("Detection",colnames(x))]
x <- x[,grep("AVG",colnames(x))]
x[1:2,1:4]
#              X5307141014_E.AVG_Signal X5307141024_E.AVG_Signal X5307141024_K.AVG_Signal X5307141024_F.AVG_Signal
# ILMN_1885639                  -0.1356                   2.5434                  12.4263                   1.9173
# ILMN_2209417                  18.5837                   8.9706                  -4.8892                   0.9122

## Step 2: Create an EList object ----
data <- new('EListRaw')
data@.Data[[1]] <- 'illumina'
data@.Data[[3]] <- BGX
data@.Data[[4]] <- x
data@.Data[[5]] <- NULL
data$E <- x
data$genes <- BGX
data$other$Detection <- detectionpvalues

### Background correction and normalization ----
data_norm <- neqc(data)

## Remove unexpressed probes and controls
controls <- data_norm$genes$controls$Probe_Id
expressed <- rowSums(data_norm$other$Detection < 0.05) >= 3
data_norm$E <- data_norm$E[expressed,]
data_norm$E <- data_norm$E[-which(rownames(data_norm$E) %in% controls),]
dim(data_norm$E)
# [1] 23484   185

## Step 3: PCA and MDS analysis ----
identical(substr(p2f$array, start = 1, stop = 10), str_extract(colnames(data_norm$E),"\\d........."))
# [1] TRUE
colnames(data_norm$E) <- rownames(p2f)
### PCAtools ----
p <- pca(data_norm$E, metadata = p2f, center = T, scale = T)

scree <- screeplot(p, getComponents(p, c(1:10)))
ptime <- biplot(p, 
                showLoadings = F,
                lab = NULL, 
                colby = "time", 
                legendPosition = 'right', legendLabSize = 12, legendIconSize = 3,
                title = 'Obermoser QC',
                subtitle = 'p2f samples')
psex <- biplot(p, 
               showLoadings = F,
               lab = NULL, 
               colby = "sex", 
               legendPosition = 'right', legendLabSize = 9, legendIconSize = 3)
psub <- biplot(p, 
               showLoadings = F,
               lab = NULL, 
               colby = "subject",
               legendPosition = 'right', legendLabSize = 9, legendIconSize = 3)
ptre <-  biplot(p, 
                showLoadings = F,
                lab = NULL, 
                colby = "treat",
                legendPosition = 'right', legendLabSize = 9, legendIconSize = 3)
peigencor <- eigencorplot(p,
                          components = getComponents(p, 1:10),
                          metavars = c("time", 'sex', 'race', 'ethnicity', 'treat'),
                          cexCorval = 0.5,
                          fontCorval = 1, 
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
botrow <- plot_grid(ptre, scree, peigencor, align = 'h', nrow = 1)
pdf(file = paste(output_dir, "pca_p2f_QC.pdf", sep = ""), height = 12, width = 16)
plot_grid(toprow, botrow, ncol = 1, align = 'v')
dev.off()
rownames(p$loadings) <- make.names(mapIds(x = illuminaHumanv3.db, keys = rownames(p$loadings), column = "SYMBOL", keytype = "PROBEID"), unique = T)
pdf(file = paste(output_dir, "ethnicity_p2f_pca.pdf", sep = ""))
biplot(p, 
       showLoadings = T,
       lab = NULL, 
       colby = "ethnicity", 
       legendPosition = 'right', legendLabSize = 10, legendIconSize = 3,
       title = 'Ethnicity of samples',
       subtitle = 'p2f samples')
dev.off()

### MDS ----
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colors <- sample(col_vector, length(unique(p2f$subject)))
new_colors <- c()
for (i in 1:length(p2f$subject)) {
  new_colors[i] <- switch(p2f$subject[i],
                          "18" = colors[1], "10" = colors[2], "24" = colors[3], "1" = colors[4], 
                          "22" = colors[5], "8" = colors[6], "6" = colors[7], "13" = colors[8],
                          "3" = colors[9], "14" = colors[10], "20" = colors[11], "7" = colors[12],
                          "4" = colors[13], "15" = colors[14], "21" = colors[15], "26" = colors[16],
                          "11" = colors[17], "12" = colors[18]
  )
}
pdf(file = paste(output_dir, "mds_p2f_qc.pdf", sep = ""))
limma::plotMDS(x = data_norm$E, 
               cex = 1, 
               var.explained = T, 
               col = new_colors,
               labels = paste(p2v$subject, p2v$sex, p2v$treat, sep = ", "))
title("MDS Plot - training set finger")
dev.off()

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#
## Step 4: limma analysis ----
#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

id <- factor(p2f$subject)
time <- factor(p2f$time)
sex <- factor(p2f$sex)
vax <- factor(p2f$treat)
eth <- factor(p2f$ethnicity)
race <- factor(p2f$race)

df <- data.frame(data_norm$E)

TS <- factor(paste(sex,time,vax,race,eth, sep = "."))
design <- model.matrix(~0+TS, df)
colnames(design) <- make.names(levels(TS))
cor <- duplicateCorrelation(df, design, block=id)
cor$consensus.correlation
# [1] 0.2495819

fit <- lmFit(object=df, design=design, block=id, correlation=cor$consensus.correlation)

# contrasts
# Levels: -168 < 0 < 1.5 < 3 < 6 < 9 < 12 < 15 < 24 < 36 < 48
my_contrasts <- makeContrasts(# Time effects
  BaselineVSzeroAll = (F.0.FLUZONE.Asian.NA + F.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.0.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.0.Saline.African.American.NA + F.0.Saline.Caucasian.Hispanic + F.0.Saline.Caucasian.Non.Hispanic.or.Latino + M.0.FLUZONE.Asian.NA + M.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.0.PNEUM.African.American.NA + M.0.PNEUM.Caucasian.Hispanic + M.0.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.0.Saline.Caucasian.Non.Hispanic.or.Latino)/12-
    (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/12,
  BaselineVSone.5All = (F.1.5.FLUZONE.Asian.NA + F.1.5.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.1.5.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.1.5.Saline.African.American.NA + F.1.5.Saline.Caucasian.Hispanic + F.1.5.Saline.Caucasian.Non.Hispanic.or.Latino + M.1.5.FLUZONE.Asian.NA + M.1.5.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.1.5.PNEUM.African.American.NA + M.1.5.PNEUM.Caucasian.Hispanic + M.1.5.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.1.5.Saline.Caucasian.Non.Hispanic.or.Latino)/12-
    (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/12,
  BaselineVSthreeAll = (F.3.FLUZONE.Asian.NA + F.3.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.3.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.3.Saline.African.American.NA + F.3.Saline.Caucasian.Hispanic + F.3.Saline.Caucasian.Non.Hispanic.or.Latino + M.3.FLUZONE.Asian.NA + M.3.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.3.PNEUM.African.American.NA + M.3.PNEUM.Caucasian.Hispanic + M.3.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.3.Saline.Caucasian.Non.Hispanic.or.Latino)/12-
    (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/12,
  BaselineVSsixAll = (F.6.FLUZONE.Asian.NA + F.6.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.6.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.6.Saline.African.American.NA + F.6.Saline.Caucasian.Hispanic + F.6.Saline.Caucasian.Non.Hispanic.or.Latino + M.6.FLUZONE.Asian.NA + M.6.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.6.PNEUM.African.American.NA + M.6.PNEUM.Caucasian.Hispanic + M.6.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.6.Saline.Caucasian.Non.Hispanic.or.Latino)/12-
    (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/12,
  BaselineVSnineAll = (F.9.FLUZONE.Asian.NA + F.9.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.9.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.9.Saline.African.American.NA + F.9.Saline.Caucasian.Hispanic + F.9.Saline.Caucasian.Non.Hispanic.or.Latino + M.9.FLUZONE.Asian.NA + M.9.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.9.PNEUM.African.American.NA + M.9.PNEUM.Caucasian.Hispanic + M.9.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.9.Saline.Caucasian.Non.Hispanic.or.Latino)/12-
    (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/12,
  BaselineVStwelveAll = (F.12.FLUZONE.Asian.NA + F.12.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.12.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.12.Saline.African.American.NA + F.12.Saline.Caucasian.Hispanic + F.12.Saline.Caucasian.Non.Hispanic.or.Latino + M.12.FLUZONE.Asian.NA + M.12.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.12.PNEUM.African.American.NA + M.12.PNEUM.Caucasian.Hispanic + M.12.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.12.Saline.Caucasian.Non.Hispanic.or.Latino)/12-
    (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/12,
  BaselineVSfifteenAll = (F.15.FLUZONE.Asian.NA + F.15.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.15.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.15.Saline.African.American.NA + F.15.Saline.Caucasian.Hispanic + F.15.Saline.Caucasian.Non.Hispanic.or.Latino + M.15.FLUZONE.Asian.NA + M.15.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.15.PNEUM.African.American.NA + M.15.PNEUM.Caucasian.Hispanic + M.15.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.15.Saline.Caucasian.Non.Hispanic.or.Latino)/12-
    (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/12,
  BaselineVStwentyFourAll = (F.24.FLUZONE.Asian.NA + F.24.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.24.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.24.Saline.African.American.NA + F.24.Saline.Caucasian.Hispanic + F.24.Saline.Caucasian.Non.Hispanic.or.Latino + M.24.FLUZONE.Asian.NA + M.24.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.24.PNEUM.African.American.NA + M.24.PNEUM.Caucasian.Hispanic + M.24.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.24.Saline.Caucasian.Non.Hispanic.or.Latino)/12-
    (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/12,
  BaselineVSthirtySixAll = (F.36.FLUZONE.Asian.NA + F.36.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.36.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.36.Saline.African.American.NA + F.36.Saline.Caucasian.Hispanic + F.36.Saline.Caucasian.Non.Hispanic.or.Latino + M.36.FLUZONE.Asian.NA + M.36.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.36.PNEUM.African.American.NA + M.36.PNEUM.Caucasian.Hispanic + M.36.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.36.Saline.Caucasian.Non.Hispanic.or.Latino)/12-
    (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/12,
  BaselineVSfortyEightAll = (F.48.FLUZONE.Asian.NA + F.48.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.48.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.48.Saline.African.American.NA + F.48.Saline.Caucasian.Hispanic + F.48.Saline.Caucasian.Non.Hispanic.or.Latino + M.48.FLUZONE.Asian.NA + M.48.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.48.PNEUM.African.American.NA + M.48.PNEUM.Caucasian.Hispanic + M.48.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.48.Saline.Caucasian.Non.Hispanic.or.Latino)/12-
    (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/12,
  # Saline times
  SalineBaselineVSzero = (F.0.Saline.African.American.NA + F.0.Saline.Caucasian.Hispanic + F.0.Saline.Caucasian.Non.Hispanic.or.Latino + M.0.Saline.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/4,
  SalineBaselineVSone.5 = (F.1.5.Saline.African.American.NA + F.1.5.Saline.Caucasian.Hispanic + F.1.5.Saline.Caucasian.Non.Hispanic.or.Latino + M.1.5.Saline.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/4,
  SalineBaselineVSthree = (F.3.Saline.African.American.NA + F.3.Saline.Caucasian.Hispanic + F.3.Saline.Caucasian.Non.Hispanic.or.Latino + M.3.Saline.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/4,
  SalineBaselineVSsix = (F.6.Saline.African.American.NA + F.6.Saline.Caucasian.Hispanic + F.6.Saline.Caucasian.Non.Hispanic.or.Latino + M.6.Saline.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/4,
  SalineBaselineVSnine = (F.9.Saline.African.American.NA + F.9.Saline.Caucasian.Hispanic + F.9.Saline.Caucasian.Non.Hispanic.or.Latino + M.9.Saline.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/4,
  SalineBaselineVStwelve = (F.12.Saline.African.American.NA + F.12.Saline.Caucasian.Hispanic + F.12.Saline.Caucasian.Non.Hispanic.or.Latino + M.12.Saline.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/4,
  SalineBaselineVSfifteen = (F.15.Saline.African.American.NA + F.15.Saline.Caucasian.Hispanic + F.15.Saline.Caucasian.Non.Hispanic.or.Latino + M.15.Saline.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/4,
  SalineBaselineVStwentyFour = (F.24.Saline.African.American.NA + F.24.Saline.Caucasian.Hispanic + F.24.Saline.Caucasian.Non.Hispanic.or.Latino + M.24.Saline.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/4,
  SalineBaselineVSthirtySix = (F.36.Saline.African.American.NA + F.36.Saline.Caucasian.Hispanic + F.36.Saline.Caucasian.Non.Hispanic.or.Latino + M.36.Saline.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/4,
  SalineBaselineVSfortyEight = (F.48.Saline.African.American.NA + F.48.Saline.Caucasian.Hispanic + F.48.Saline.Caucasian.Non.Hispanic.or.Latino + M.48.Saline.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/4,
  # PNEUM times
  PNEUMBaselineVSzero = (F.0.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.0.PNEUM.African.American.NA + M.0.PNEUM.Caucasian.Hispanic + M.0.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4,
  PNEUMBaselineVSone.5 = (F.1.5.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.1.5.PNEUM.African.American.NA + M.1.5.PNEUM.Caucasian.Hispanic + M.1.5.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4,
  PNEUMBaselineVSthree = (F.3.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.3.PNEUM.African.American.NA + M.3.PNEUM.Caucasian.Hispanic + M.3.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4,
  PNEUMBaselineVSsix = (F.6.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.6.PNEUM.African.American.NA + M.6.PNEUM.Caucasian.Hispanic + M.6.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4,
  PNEUMBaselineVSnine = (F.9.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.9.PNEUM.African.American.NA + M.9.PNEUM.Caucasian.Hispanic + M.9.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4,
  PNEUMBaselineVStwelve = (F.12.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.12.PNEUM.African.American.NA + M.12.PNEUM.Caucasian.Hispanic + M.12.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4,
  PNEUMBaselineVSfifteen = (F.15.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.15.PNEUM.African.American.NA + M.15.PNEUM.Caucasian.Hispanic + M.15.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4,
  PNEUMBaselineVStwentyFour = (F.24.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.24.PNEUM.African.American.NA + M.24.PNEUM.Caucasian.Hispanic + M.24.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4,
  PNEUMBaselineVSthirtySix = (F.36.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.36.PNEUM.African.American.NA + M.36.PNEUM.Caucasian.Hispanic + M.36.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4,
  PNEUMBaselineVSfortyEight = (F.48.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.48.PNEUM.African.American.NA + M.48.PNEUM.Caucasian.Hispanic + M.48.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino)/4,
  # Flu times
  FLUBaselineVSzero = (F.0.FLUZONE.Asian.NA + F.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.0.FLUZONE.Asian.NA + M.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4,
  FLUBaselineVSone.5 = (F.1.5.FLUZONE.Asian.NA + F.1.5.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.1.5.FLUZONE.Asian.NA + M.1.5.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4,
  FLUBaselineVSthree = (F.3.FLUZONE.Asian.NA + F.3.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.3.FLUZONE.Asian.NA + M.3.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4,
  FLUBaselineVSsix = (F.6.FLUZONE.Asian.NA + F.6.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.6.FLUZONE.Asian.NA + M.6.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4,
  FLUBaselineVSnine = (F.9.FLUZONE.Asian.NA + F.9.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.9.FLUZONE.Asian.NA + M.9.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4,
  FLUBaselineVStwelve = (F.12.FLUZONE.Asian.NA + F.12.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.12.FLUZONE.Asian.NA + M.12.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4,
  FLUBaselineVSfifteen = (F.15.FLUZONE.Asian.NA + F.15.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.15.FLUZONE.Asian.NA + M.15.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4,
  FLUBaselineVStwentyFour = (F.24.FLUZONE.Asian.NA + F.24.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.24.FLUZONE.Asian.NA + M.24.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4,
  FLUBaselineVSthirtySix = (F.36.FLUZONE.Asian.NA + F.36.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.36.FLUZONE.Asian.NA + M.36.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4,
  FLUBaselineVSfortyEight = (F.48.FLUZONE.Asian.NA + F.48.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.48.FLUZONE.Asian.NA + M.48.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4-
    (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino)/4,
  # Time-sex interaction Male - Female
  IMaleBaselineVSzero = ((M.0.FLUZONE.Asian.NA + M.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.0.PNEUM.African.American.NA + M.0.PNEUM.Caucasian.Hispanic + M.0.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.0.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
                           (M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6) -
    ((F.0.FLUZONE.Asian.NA + F.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.0.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.0.Saline.African.American.NA + F.0.Saline.Caucasian.Hispanic + F.0.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6),
  IMaleBaselineVSone.5 = ((M.1.5.FLUZONE.Asian.NA + M.1.5.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.1.5.PNEUM.African.American.NA + M.1.5.PNEUM.Caucasian.Hispanic + M.1.5.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.1.5.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
                            (M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6) -
    ((F.1.5.FLUZONE.Asian.NA + F.1.5.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.1.5.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.1.5.Saline.African.American.NA + F.1.5.Saline.Caucasian.Hispanic + F.1.5.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6),
  IMaleBaselineVSthree = ((M.3.FLUZONE.Asian.NA + M.3.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.3.PNEUM.African.American.NA + M.3.PNEUM.Caucasian.Hispanic + M.3.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.3.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
                            (M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6) -
    ((F.3.FLUZONE.Asian.NA + F.3.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.3.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.3.Saline.African.American.NA + F.3.Saline.Caucasian.Hispanic + F.3.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6),
  IMaleBaselineVSsix = ((M.6.FLUZONE.Asian.NA + M.6.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.6.PNEUM.African.American.NA + M.6.PNEUM.Caucasian.Hispanic + M.6.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.6.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
                          (M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6) -
    ((F.6.FLUZONE.Asian.NA + F.6.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.6.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.6.Saline.African.American.NA + F.6.Saline.Caucasian.Hispanic + F.6.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6),
  IMaleBaselineVSnine = ((M.9.FLUZONE.Asian.NA + M.9.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.9.PNEUM.African.American.NA + M.9.PNEUM.Caucasian.Hispanic + M.9.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.9.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
                           (M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6) -
    ((F.9.FLUZONE.Asian.NA + F.9.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.9.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.9.Saline.African.American.NA + F.9.Saline.Caucasian.Hispanic + F.9.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6),
  IMaleBaselineVStwelve = ((M.12.FLUZONE.Asian.NA + M.12.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.12.PNEUM.African.American.NA + M.12.PNEUM.Caucasian.Hispanic + M.12.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.12.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
                             (M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6) -
    ((F.12.FLUZONE.Asian.NA + F.12.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.12.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.12.Saline.African.American.NA + F.12.Saline.Caucasian.Hispanic + F.12.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6),
  IMaleBaselineVSfifteen = ((M.15.FLUZONE.Asian.NA + M.15.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.15.PNEUM.African.American.NA + M.15.PNEUM.Caucasian.Hispanic + M.15.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.15.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
                              (M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6) -
    ((F.15.FLUZONE.Asian.NA + F.15.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.15.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.15.Saline.African.American.NA + F.15.Saline.Caucasian.Hispanic + F.15.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6),
  IMaleBaselineVStwentyFour = ((M.24.FLUZONE.Asian.NA + M.24.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.24.PNEUM.African.American.NA + M.24.PNEUM.Caucasian.Hispanic + M.24.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.24.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
                                 (M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6) -
    ((F.24.FLUZONE.Asian.NA + F.24.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.24.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.24.Saline.African.American.NA + F.24.Saline.Caucasian.Hispanic + F.24.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6),
  IMaleBaselineVSthirtySix = ((M.36.FLUZONE.Asian.NA + M.36.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.36.PNEUM.African.American.NA + M.36.PNEUM.Caucasian.Hispanic + M.36.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.36.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
                                (M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6) -
    ((F.36.FLUZONE.Asian.NA + F.36.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.36.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.36.Saline.African.American.NA + F.36.Saline.Caucasian.Hispanic + F.36.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6),
  IMaleBaselineVSfortyEight = ((M.48.FLUZONE.Asian.NA + M.48.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.48.PNEUM.African.American.NA + M.48.PNEUM.Caucasian.Hispanic + M.48.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.48.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
                                 (M..168.FLUZONE.Asian.NA + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.African.American.NA + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6) -
    ((F.48.FLUZONE.Asian.NA + F.48.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.48.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.48.Saline.African.American.NA + F.48.Saline.Caucasian.Hispanic + F.48.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..168.FLUZONE.Asian.NA + F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.African.American.NA + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6),
  # Caucasian only
  CaucasianBaseVSzero = (F.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.0.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.0.Saline.Caucasian.Hispanic + F.0.Saline.Caucasian.Non.Hispanic.or.Latino + M.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.0.PNEUM.Caucasian.Hispanic + M.0.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.0.Saline.Caucasian.Non.Hispanic.or.Latino)/8- 
    (F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/8,
  CaucasianBaseVSone.5 = (F.1.5.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.1.5.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.1.5.Saline.Caucasian.Hispanic + F.1.5.Saline.Caucasian.Non.Hispanic.or.Latino + M.1.5.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.1.5.PNEUM.Caucasian.Hispanic + M.1.5.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.1.5.Saline.Caucasian.Non.Hispanic.or.Latino)/8- 
    (F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/8,
  CaucasianBaseVSthree = (F.3.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.3.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.3.Saline.Caucasian.Hispanic + F.3.Saline.Caucasian.Non.Hispanic.or.Latino + M.3.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.3.PNEUM.Caucasian.Hispanic + M.3.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.3.Saline.Caucasian.Non.Hispanic.or.Latino)/8- 
    (F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/8,
  CaucasianBaseVSsix = (F.6.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.6.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.6.Saline.Caucasian.Hispanic + F.6.Saline.Caucasian.Non.Hispanic.or.Latino + M.6.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.6.PNEUM.Caucasian.Hispanic + M.6.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.6.Saline.Caucasian.Non.Hispanic.or.Latino)/8- 
    (F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/8,
  CaucasianBaseVSnine = (F.9.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.9.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.9.Saline.Caucasian.Hispanic + F.9.Saline.Caucasian.Non.Hispanic.or.Latino + M.9.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.9.PNEUM.Caucasian.Hispanic + M.9.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.9.Saline.Caucasian.Non.Hispanic.or.Latino)/8- 
    (F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/8,
  CaucasianBaseVStwelve = (F.12.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.12.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.12.Saline.Caucasian.Hispanic + F.12.Saline.Caucasian.Non.Hispanic.or.Latino + M.12.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.12.PNEUM.Caucasian.Hispanic + M.12.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.12.Saline.Caucasian.Non.Hispanic.or.Latino)/8- 
    (F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/8,
  CaucasianBaseVSfifteen = (F.15.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.15.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.15.Saline.Caucasian.Hispanic + F.15.Saline.Caucasian.Non.Hispanic.or.Latino + M.15.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.15.PNEUM.Caucasian.Hispanic + M.15.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.15.Saline.Caucasian.Non.Hispanic.or.Latino)/8- 
    (F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/8,
  CaucasianBaseVStwentyFour = (F.24.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.24.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.24.Saline.Caucasian.Hispanic + F.24.Saline.Caucasian.Non.Hispanic.or.Latino + M.24.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.24.PNEUM.Caucasian.Hispanic + M.24.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.24.Saline.Caucasian.Non.Hispanic.or.Latino)/8- 
    (F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/8,
  CaucasianBaseVSthirtySix = (F.36.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.36.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.36.Saline.Caucasian.Hispanic + F.36.Saline.Caucasian.Non.Hispanic.or.Latino + M.36.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.36.PNEUM.Caucasian.Hispanic + M.36.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.36.Saline.Caucasian.Non.Hispanic.or.Latino)/8- 
    (F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/8,
  CaucasianBaseVSfortyEight = (F.48.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.48.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.48.Saline.Caucasian.Hispanic + F.48.Saline.Caucasian.Non.Hispanic.or.Latino + M.48.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.48.PNEUM.Caucasian.Hispanic + M.48.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.48.Saline.Caucasian.Non.Hispanic.or.Latino)/8- 
    (F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.Caucasian.Hispanic + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/8,
  # Asian only
  AsianBaseVSzero = (F.0.FLUZONE.Asian.NA + M.0.FLUZONE.Asian.NA)/2 - (F..168.FLUZONE.Asian.NA + M..168.FLUZONE.Asian.NA)/2,
  AsianBaseVSone.5 = (F.1.5.FLUZONE.Asian.NA + M.1.5.FLUZONE.Asian.NA)/2 - (F..168.FLUZONE.Asian.NA + M..168.FLUZONE.Asian.NA)/2,
  AsianBaseVSthree = (F.3.FLUZONE.Asian.NA + M.3.FLUZONE.Asian.NA)/2 - (F..168.FLUZONE.Asian.NA + M..168.FLUZONE.Asian.NA)/2,
  AsianBaseVSsix = (F.6.FLUZONE.Asian.NA + M.6.FLUZONE.Asian.NA)/2 - (F..168.FLUZONE.Asian.NA + M..168.FLUZONE.Asian.NA)/2,
  AsianBaseVSnine = (F.9.FLUZONE.Asian.NA + M.9.FLUZONE.Asian.NA)/2 - (F..168.FLUZONE.Asian.NA + M..168.FLUZONE.Asian.NA)/2,
  AsianBaseVStwelve = (F.12.FLUZONE.Asian.NA + M.12.FLUZONE.Asian.NA)/2 - (F..168.FLUZONE.Asian.NA + M..168.FLUZONE.Asian.NA)/2,
  AsianBaseVSfifteen = (F.15.FLUZONE.Asian.NA + M.15.FLUZONE.Asian.NA)/2 - (F..168.FLUZONE.Asian.NA + M..168.FLUZONE.Asian.NA)/2,
  AsianBaseVStwentyFour = (F.24.FLUZONE.Asian.NA + M.24.FLUZONE.Asian.NA)/2 - (F..168.FLUZONE.Asian.NA + M..168.FLUZONE.Asian.NA)/2,
  AsianBaseVSthirtySix = (F.36.FLUZONE.Asian.NA + M.36.FLUZONE.Asian.NA)/2 - (F..168.FLUZONE.Asian.NA + M..168.FLUZONE.Asian.NA)/2,
  AsianBaseVSfortyEight = (F.48.FLUZONE.Asian.NA + M.48.FLUZONE.Asian.NA)/2 - (F..168.FLUZONE.Asian.NA + M..168.FLUZONE.Asian.NA)/2,
  # African American only
  AfricanAmerBaseVSzero = (F.0.Saline.African.American.NA + M.0.PNEUM.African.American.NA)/2 - (F..168.Saline.African.American.NA + M..168.PNEUM.African.American.NA)/2,
  AfricanAmerBaseVSone.5 = (F.1.5.Saline.African.American.NA + M.1.5.PNEUM.African.American.NA)/2 - (F..168.Saline.African.American.NA + M..168.PNEUM.African.American.NA)/2,
  AfricanAmerBaseVSthree = (F.3.Saline.African.American.NA + M.3.PNEUM.African.American.NA)/2 - (F..168.Saline.African.American.NA + M..168.PNEUM.African.American.NA)/2,
  AfricanAmerBaseVSsix = (F.6.Saline.African.American.NA + M.6.PNEUM.African.American.NA)/2 - (F..168.Saline.African.American.NA + M..168.PNEUM.African.American.NA)/2,
  AfricanAmerBaseVSnine = (F.9.Saline.African.American.NA + M.9.PNEUM.African.American.NA)/2 - (F..168.Saline.African.American.NA + M..168.PNEUM.African.American.NA)/2,
  AfricanAmerBaseVStwelve = (F.12.Saline.African.American.NA + M.12.PNEUM.African.American.NA)/2 - (F..168.Saline.African.American.NA + M..168.PNEUM.African.American.NA)/2,
  AfricanAmerBaseVSfifteen = (F.15.Saline.African.American.NA + M.15.PNEUM.African.American.NA)/2 - (F..168.Saline.African.American.NA + M..168.PNEUM.African.American.NA)/2,
  AfricanAmerBaseVStwentyFour = (F.24.Saline.African.American.NA + M.24.PNEUM.African.American.NA)/2 - (F..168.Saline.African.American.NA + M..168.PNEUM.African.American.NA)/2,
  AfricanAmerBaseVSthirtySix = (F.36.Saline.African.American.NA + M.36.PNEUM.African.American.NA)/2 - (F..168.Saline.African.American.NA + M..168.PNEUM.African.American.NA)/2,
  AfricanAmerBaseVSfortyEight = (F.48.Saline.African.American.NA + M.48.PNEUM.African.American.NA)/2 - (F..168.Saline.African.American.NA + M..168.PNEUM.African.American.NA)/2,
  # Ethnicity-time interaction Hispanic - non-hispanic
  IHispanicBaseVSzero = ((F.0.Saline.Caucasian.Hispanic + M.0.PNEUM.Caucasian.Hispanic)/2-(F..168.Saline.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Hispanic)/2) -
    ((F.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.0.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.0.Saline.Caucasian.Non.Hispanic.or.Latino + M.0.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.0.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.0.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6) ,
  IHispanicBaseVSone.5 = ((F.1.5.Saline.Caucasian.Hispanic + M.1.5.PNEUM.Caucasian.Hispanic)/2-(F..168.Saline.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Hispanic)/2) -
    ((F.1.5.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.1.5.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.1.5.Saline.Caucasian.Non.Hispanic.or.Latino + M.1.5.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.1.5.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.1.5.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6) ,
  IHispanicBaseVSthree = ((F.3.Saline.Caucasian.Hispanic + M.3.PNEUM.Caucasian.Hispanic)/2-(F..168.Saline.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Hispanic)/2) -
    ((F.3.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.3.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.3.Saline.Caucasian.Non.Hispanic.or.Latino + M.3.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.3.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.3.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6) ,
  IHispanicBaseVSsix = ((F.6.Saline.Caucasian.Hispanic + M.6.PNEUM.Caucasian.Hispanic)/2-(F..168.Saline.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Hispanic)/2) -
    ((F.6.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.6.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.6.Saline.Caucasian.Non.Hispanic.or.Latino + M.6.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.6.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.6.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6) ,
  IHispanicBaseVSnine = ((F.9.Saline.Caucasian.Hispanic + M.9.PNEUM.Caucasian.Hispanic)/2-(F..168.Saline.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Hispanic)/2) -
    ((F.9.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.9.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.9.Saline.Caucasian.Non.Hispanic.or.Latino + M.9.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.9.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.9.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6) ,
  IHispanicBaseVStwelve = ((F.12.Saline.Caucasian.Hispanic + M.12.PNEUM.Caucasian.Hispanic)/2-(F..168.Saline.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Hispanic)/2) -
    ((F.12.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.12.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.12.Saline.Caucasian.Non.Hispanic.or.Latino + M.12.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.12.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.12.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6) ,
  IHispanicBaseVSfifteen = ((F.15.Saline.Caucasian.Hispanic + M.15.PNEUM.Caucasian.Hispanic)/2-(F..168.Saline.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Hispanic)/2) -
    ((F.15.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.15.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.15.Saline.Caucasian.Non.Hispanic.or.Latino + M.15.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.15.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.15.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6) ,
  IHispanicBaseVStwentyFour = ((F.24.Saline.Caucasian.Hispanic + M.24.PNEUM.Caucasian.Hispanic)/2-(F..168.Saline.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Hispanic)/2) -
    ((F.24.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.24.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.24.Saline.Caucasian.Non.Hispanic.or.Latino + M.24.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.24.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.24.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6) ,
  IHispanicBaseVSthirtySix = ((F.36.Saline.Caucasian.Hispanic + M.36.PNEUM.Caucasian.Hispanic)/2-(F..168.Saline.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Hispanic)/2) -
    ((F.36.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.36.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.36.Saline.Caucasian.Non.Hispanic.or.Latino + M.36.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.36.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.36.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6) ,
  IHispanicBaseVSfortyEight = ((F.48.Saline.Caucasian.Hispanic + M.48.PNEUM.Caucasian.Hispanic)/2-(F..168.Saline.Caucasian.Hispanic + M..168.PNEUM.Caucasian.Hispanic)/2) -
    ((F.48.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F.48.PNEUM.Caucasian.Non.Hispanic.or.Latino + F.48.Saline.Caucasian.Non.Hispanic.or.Latino + M.48.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M.48.PNEUM.Caucasian.Non.Hispanic.or.Latino + M.48.Saline.Caucasian.Non.Hispanic.or.Latino)/6-
       (F..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + F..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + F..168.Saline.Caucasian.Non.Hispanic.or.Latino + M..168.FLUZONE.Caucasian.Non.Hispanic.or.Latino + M..168.PNEUM.Caucasian.Non.Hispanic.or.Latino + M..168.Saline.Caucasian.Non.Hispanic.or.Latino)/6) ,
  levels = colnames(design))

fit2 <- contrasts.fit(fit, my_contrasts)
x <- eBayes(fit2)

### Results tables ----

y <- topTable(x, number = Inf, adjust.method = "fdr")
y$Symbol <- mapIds(x = illuminaHumanv3.db, keys = row.names(y), column = "SYMBOL", keytype = "PROBEID")
dir.create(paste(output_dir, "training_finger_limma", sep = ""))
write.csv(y, file = paste(output_dir, "training_finger_limma/limma_F.csv", sep = ""), row.names = TRUE)

dt <- decideTests(x)
summary(dt)
#        BaselineVSzeroAll BaselineVSone.5All BaselineVSthreeAll BaselineVSsixAll
# Down                   5                 10                  8               54
# NotSig             23478              23446              23467            23308
# Up                     1                 28                  9              122
#        BaselineVSnineAll BaselineVStwelveAll BaselineVSfifteenAll BaselineVStwentyFourAll
# Down                 324                 395                 1291                      20
# NotSig             22614               22624                21166                   23429
# Up                   546                 465                 1027                      35
#        BaselineVSthirtySixAll BaselineVSfortyEightAll 
# Down                      126                      28               
# NotSig                  23169                   23399               
# Up                        189                      57
dt <- decideTests(x, adjust.method = 'none')
summary(dt)
#        BaselineVSzeroAll BaselineVSone.5All BaselineVSthreeAll BaselineVSsixAll
# Down                 710               1276               1176             1424
# NotSig             21723              20607              20648            20123
# Up                  1051               1601               1660             1937
#        BaselineVSnineAll BaselineVStwelveAll BaselineVSfifteenAll BaselineVStwentyFourAll
# Down                2110                2046                 2920                    1144
# NotSig             18893               19009                17473                   20876
# Up                  2481                2429                 3091                    1464
#        BaselineVSthirtySixAll BaselineVSfortyEightAll 
# Down                     1511                    1286             
# NotSig                  19877                   20540              
# Up                       2096                    1658

names(y)
tables <- list()
for (i in 1:90) {
  tables[[i]] <- topTable(x, coef = i, number = Inf, adjust.method = "fdr")
  tables[[i]]$Symbol <- mapIds(x = illuminaHumanv3.db, keys = row.names(tables[[i]]), column = "SYMBOL", keytype = "PROBEID")
  write.csv(tables[[i]], file = paste(output_dir, "training_finger_limma/",
                                      names(y)[i],".csv",sep = ""), row.names = TRUE)
}

### Plot p-value histograms ----
histlist <- list()
for (i in 1:10) {
  histlist[[i]] <- ggplot(data = tables[[i]], aes(x = adj.P.Val)) +
    geom_histogram(fill = "lightblue", color = "black") +
    ggtitle(paste("Adj. p histogram", names(y)[i])) +
    theme_cowplot() +
    theme(title = element_text(size = 7.5))
  names(histlist)[i] <- names(y)[i]
}
pdf(file = paste(output_dir, "training_finger_limma/adj_p_value_histograms.pdf", sep = ""), width = 12)
plot_grid(plotlist = histlist, align = 'hv', nrow = 2)
dev.off()
for (i in 1:10) {
  histlist[[i]] <- ggplot(data = tables[[i]], aes(x = P.Value)) +
    geom_histogram(fill = "lightblue", color = "black") +
    ggtitle(paste("P-value histogram", names(y)[i])) +
    theme_cowplot() +
    theme(title = element_text(size = 7.5))
  names(histlist)[i] <- names(y)[i]
}
pdf(file = paste(output_dir, "training_finger_limma/p_value_histograms.pdf", sep = ""), width = 12)
plot_grid(plotlist = histlist, align = 'hv', nrow = 2)
dev.off()

### Create volcano plots ----
volcanos <- list()
for (i in 1:10) {
  volcanos[[i]] <- EnhancedVolcano(tables[[i]],
                                   lab = tables[[i]]$Symbol,
                                   x = "logFC",
                                   y = "P.Value",
                                   title = "Training_Finger", titleLabSize = 6,
                                   subtitle = names(y)[i], subtitleLabSize = 5,
                                   labSize = (4.0),
                                   colAlpha = 1,
                                   pointSize = 1.0, 
                                   pCutoff = 1e-2,
                                   boxedLabels = FALSE,
                                   drawConnectors = TRUE,
                                   legendLabSize = 4,
                                   legendIconSize = 2,
                                   axisLabSize = 5,
                                   caption = "")
  names(volcanos)[i] <- names(y)[i]
}
pdf(file = paste(output_dir, "training_finger_limma/volcanos.pdf", sep = ""), height = 10, width = 14)
plot_grid(plotlist = volcanos, align = 'h', nrow = 2)
dev.off()