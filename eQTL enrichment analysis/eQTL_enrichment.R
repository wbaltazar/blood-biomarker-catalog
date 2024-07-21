# Date: June 11, 2024

## Input: .csv files of genes with stability filters, median stability statistics, and DE counts
## Output: GSEA plots and results

setwd("./Desktop/work_repo/Box organization")
library(tidyverse)
library(clusterProfiler)

## Tables
stable_statistics <- read.csv("./type 2 gene analysis/10062024 analysis/median_stability_statistics.csv")
head(stable_statistics)
head(stable_statistics[order(stable_statistics$repeatability, decreasing = T), "Symbol"], 20)
# [1] "LOC102724528" "LOC644172"    "LOC647859"    "LOC101929545" "X8-Mar"       "LOC100233156" "LOC101927752"
# [8] "LOC102723839" "CTC-490G23"   "LOC101929354" "CTD-2006K23"  "LOC100507460" "AC131056"     "LOC101928192"
# [15] "OVOS"         "LOC102723750" "FLJ38379"     "GPS2P1"       "FAM46C"       "FLJ20444" 

## Data from source("./type 2 gene analysis/06062024 analysis/new_finding_type_2.R")
symbols_list <- lapply(data, symbol_find)
symbols_list <- lapply(symbols_list, function(x) {str_replace_all(x$Symbol, "\\.(?!\\d$)", "-") %>% str_remove_all(pattern = "\\..")})
common_symbols <- c()
for (i in 1:(length(symbols_list) - 1)) {
  if (i == 1) {
    common_symbols <- intersect(symbols_list[[i]], symbols_list[[i + 1]])
  } else {
    common_symbols <- intersect(common_symbols, symbols_list[[i + 1]])
  }
}
length(common_symbols)
# [1] 9474
write.table(common_symbols, "./1results/Symbols_common_to_all_datasets.txt")
sum(is.na(common_symbols))
# [1] 0
common_symbol_statistics <- stable_statistics[stable_statistics$Symbol %in% common_symbols, ]
dim(common_symbol_statistics)
# [1] 9474   10
write.csv(common_symbol_statistics, "./eQTL enrichment analysis/median_stability_for_ubiquitous_genes.csv")
head(common_symbol_statistics[order(common_symbol_statistics$repeatability, decreasing = T), "Symbol"], 20)
# [1] "MYOM2"    "TMEM176A" "FAM118A"  "RPS26"    "HLA-DRB5" "FOLR3"    "GZMH"     "CFD"      "GSTM1"    "TUBB2A"  
# [11] "CLEC12A"  "PAM"      "HLA-DQB1" "GNLY"     "SMIM1"    "KRT1"     "TMEM158"  "ZNF683"   "TCL1A"    "CENPK"   

degs <- read.csv("./DEG analysis/11062024 analysis/deg_counts.csv")
common_symbol_de <- degs[degs$Symbol %in% common_symbols, ]
head(common_symbol_de[order(common_symbol_de$P.value.Studies, decreasing = T), "Symbol"], 20)
# [1] "USP49"     "VMA21"     "TTC38"     "LRIG1"     "LINC00921" "CCDC88A"   "DCAF15"    "AGL"       "DIS3L2"   
# [10] "SMURF2"    "HMGB3"     "CA13"      "FSD1L"     "PLEKHA3"   "SLC38A2"   "HDAC9"     "RBM18"     "TCAIM"    
# [19] "SLC25A16"  "C9orf85"  
write.csv(common_symbol_de,"./eQTL enrichment analysis/deg_counts_for_ubiquitous_genes.csv", row.names = F)

## Gene Set Enrichment Analysis ---- 
## (tutorial: https://www.youtube.com/watch?v=Mi6u4r0lJvo)
### Create gene set for top 500 GTEx eQTLs ----
egenes <- read_tsv("../GTEx_Analysis_v8_eQTL/Whole_Blood.v8.egenes.txt.gz")
length(unique(egenes$gene_name))
# [1] 20299
range(egenes$pval_nominal)
# [1] 7.72767e-278  1.47722e-02
## Use the q-value, or the FDR-controlled, beta-distributed p-value of the permutated norminal p-values
egenes <- egenes[order(egenes$qval, decreasing = F),]
gtex_eqtls <- egenes[1:500, "gene_name"]
write.csv(gtex_eqtls, "./eQTL enrichment analysis/GTEx eQTLs used for GSEA.csv")

### Median repeatabilty enrichment ----
common_symbol_statistics <- common_symbol_statistics[order(common_symbol_statistics$repeatability, decreasing = T),]
repeatability <- common_symbol_statistics$repeatability
names(repeatability) <- common_symbol_statistics$Symbol
eqtlgse <- clusterProfiler::GSEA(geneList = repeatability, 
                                 TERM2GENE = data.frame(term = "GTEx_eQTLs", gene = gtex_eqtls$gene_name),
                                 eps = 1e-300, pvalueCutoff = 1,
                                 minGSSize = 10, maxGSSize = 500)
gseaplot(eqtlgse, geneSetID = 1, title = "GTEx eQTL enrichment by repeatability")
as.data.frame(eqtlgse)
#                    ID Description setSize enrichmentScore      NES       pvalue     p.adjust qvalue rank
# GTEx_eQTLs GTEx_eQTLs  GTEx_eQTLs     262       0.5596715 2.738177 4.102507e-40 4.102507e-40     NA 1686
# leading_edge
# GTEx_eQTLs tags=58%, list=18%, signal=49%

### Median subject_vp enrichment ----
common_symbol_statistics <- common_symbol_statistics[order(common_symbol_statistics$subject, decreasing = T),]
subject <- common_symbol_statistics$subject
names(subject) <- common_symbol_statistics$Symbol
eqtlgse <- clusterProfiler::GSEA(geneList = subject, 
                                 TERM2GENE = data.frame(term = "GTEx_eQTLs", gene = gtex_eqtls$gene_name),
                                 eps = 1e-300, pvalueCutoff = 1,
                                 minGSSize = 10, maxGSSize = 500)
gseaplot(eqtlgse, geneSetID = 1, title = "GTEx eQTL enrichment by subject varPart")
as.data.frame(eqtlgse)
#                    ID Description setSize enrichmentScore      NES       pvalue     p.adjust qvalue rank
# GTEx_eQTLs GTEx_eQTLs  GTEx_eQTLs     262       0.5675835 2.800669 1.227463e-41 1.227463e-41     NA 1979
# leading_edge
# GTEx_eQTLs tags=62%, list=21%, signal=50%

### Median Within Var enrichment ----
common_symbol_statistics <- common_symbol_statistics[order(common_symbol_statistics$Within.Variation..SD., decreasing = T),]
within <- common_symbol_statistics$Within.Variation..SD.
names(within) <- common_symbol_statistics$Symbol
eqtlgse <- clusterProfiler::GSEA(geneList = within, 
                                 TERM2GENE = data.frame(term = "GTEx_eQTLs", gene = gtex_eqtls$gene_name),
                                 eps = 1e-300, pvalueCutoff = 1,
                                 minGSSize = 10, maxGSSize = 500)
gseaplot(eqtlgse, geneSetID = 1, title = "GTEx eQTL enrichment by within SD")
as.data.frame(eqtlgse)
# ID Description setSize enrichmentScore      NES pvalue p.adjust qvalue rank
# GTEx_eQTLs GTEx_eQTLs  GTEx_eQTLs     262       0.1170396 1.081992  0.365    0.365     NA 4946
# leading_edge
# GTEx_eQTLs tags=56%, list=52%, signal=28%

### Median Total Var enrichment ----
common_symbol_statistics <- common_symbol_statistics[order(common_symbol_statistics$Total.Variation..SD., decreasing = T),]
total <- common_symbol_statistics$Total.Variation..SD.
names(total) <- common_symbol_statistics$Symbol
eqtlgse <- clusterProfiler::GSEA(geneList = total, 
                                 TERM2GENE = data.frame(term = "GTEx_eQTLs", gene = gtex_eqtls$gene_name),
                                 eps = 1e-300, pvalueCutoff = 1,
                                 minGSSize = 10, maxGSSize = 500)
gseaplot(eqtlgse, geneSetID = 1, title = "GTEx eQTL enrichment by total SD")
as.data.frame(eqtlgse)
# ID Description setSize enrichmentScore      NES       pvalue     p.adjust qvalue rank
# GTEx_eQTLs GTEx_eQTLs  GTEx_eQTLs     262       0.4490156 3.241552 1.365526e-26 1.365526e-26     NA 2281
# leading_edge
# GTEx_eQTLs tags=51%, list=24%, signal=40%

### Median Rs enrichment ----
common_symbol_statistics <- common_symbol_statistics[order(common_symbol_statistics$Rs, decreasing = T),]
rs <- common_symbol_statistics$Rs
names(rs) <- common_symbol_statistics$Symbol
eqtlgse <- clusterProfiler::GSEA(geneList = rs, 
                                 TERM2GENE = data.frame(term = "GTEx_eQTLs", gene = gtex_eqtls$gene_name),
                                 eps = 1e-300, pvalueCutoff = 1,
                                 minGSSize = 10, maxGSSize = 500)
gseaplot(eqtlgse, geneSetID = 1, title = "GTEx eQTL enrichment by Rs")
as.data.frame(eqtlgse)
# ID Description setSize enrichmentScore       NES       pvalue     p.adjust qvalue rank
# GTEx_eQTLs GTEx_eQTLs  GTEx_eQTLs     262      -0.3485371 -7.025066 2.629541e-32 2.629541e-32     NA 3023
# leading_edge
# GTEx_eQTLs tags=73%, list=32%, signal=51%

### Median gen.variance enrichment ----
common_symbol_statistics <- common_symbol_statistics[order(common_symbol_statistics$gen.variance, decreasing = T),]
genvar <- common_symbol_statistics$gen.variance
names(genvar) <- common_symbol_statistics$Symbol
eqtlgse <- clusterProfiler::GSEA(geneList = genvar, 
                                 TERM2GENE = data.frame(term = "GTEx_eQTLs", gene = gtex_eqtls$gene_name),
                                 eps = 1e-300, pvalueCutoff = 1,
                                 minGSSize = 10, maxGSSize = 500)
gseaplot(eqtlgse, geneSetID = 1, title = "GTEx eQTL enrichment by gen.variance")
as.data.frame(eqtlgse)
# ID Description setSize enrichmentScore      NES       pvalue     p.adjust qvalue rank
# GTEx_eQTLs GTEx_eQTLs  GTEx_eQTLs     262       0.7751781 1.787651 1.157492e-17 1.157492e-17     NA 1101
# leading_edge
# GTEx_eQTLs tags=42%, list=12%, signal=39%

### Median Average Expression enrichment ----
common_symbol_statistics <- common_symbol_statistics[order(common_symbol_statistics$Average.Expression, decreasing = T),]
avexpr <- common_symbol_statistics$Average.Expression
names(avexpr) <- common_symbol_statistics$Symbol
eqtlgse <- clusterProfiler::GSEA(geneList = avexpr, 
                                 TERM2GENE = data.frame(term = "GTEx_eQTLs", gene = gtex_eqtls$gene_name),
                                 eps = 1e-300, pvalueCutoff = 1,
                                 minGSSize = 10, maxGSSize = 500)
gseaplot(eqtlgse, geneSetID = 1, title = "GTEx eQTL enrichment by Avg. Expression")
as.data.frame(eqtlgse)
#                    ID Description setSize enrichmentScore      NES    pvalue  p.adjust qvalue rank
# GTEx_eQTLs GTEx_eQTLs  GTEx_eQTLs     262       0.1552173 1.317816 0.0969031 0.0969031     NA 4764
# leading_edge
# GTEx_eQTLs tags=58%, list=50%, signal=30%

### Stability filter enrichment ----
stable_counts_common <- read.csv("./type 2 gene analysis/06062024 analysis/stable_gene_filter_count.csv")
dim(stable_counts_common)
# [1] 28841    16
stable_counts_common <- stable_counts_common[stable_counts_common$Symbols %in% common_symbols,]
dim(stable_counts_common)
# [1] 9511   16
stable_counts_common <- stable_counts_common[order(stable_counts_common$Score, decreasing = T),]
stable <- stable_counts_common$Score
names(stable) <- stable_counts_common$Symbols
eqtlgse <- clusterProfiler::GSEA(geneList = stable, 
                                 TERM2GENE = data.frame(term = "GTEx_eQTLs", gene = gtex_eqtls$gene_name),
                                 eps = 1e-300, pvalueCutoff = 1,
                                 minGSSize = 10, maxGSSize = 500)
gseaplot(eqtlgse, geneSetID = 1, title = "GTEx eQTL enrichment by stability filter")
as.data.frame(eqtlgse)
#                    ID Description setSize enrichmentScore      NES       pvalue     p.adjust qvalue rank
# GTEx_eQTLs GTEx_eQTLs  GTEx_eQTLs     262       0.6366862 1.904505 1.252482e-32 1.252482e-32     NA 1667
# leading_edge
# GTEx_eQTLs tags=57%, list=18%, signal=49%


### Median time vp enrichment ----
common_symbol_statistics <- common_symbol_statistics[order(common_symbol_statistics$time, decreasing = T),]
time <- common_symbol_statistics$time
names(time) <- common_symbol_statistics$Symbol
eqtlgse <- clusterProfiler::GSEA(geneList = time, 
                                 TERM2GENE = data.frame(term = "GTEx_eQTLs", gene = gtex_eqtls$gene_name),
                                 eps = 1e-300, pvalueCutoff = 1,
                                 minGSSize = 10, maxGSSize = 500)
gseaplot(eqtlgse, geneSetID = 1, title = "GTEx eQTL enrichment by time varPart")
as.data.frame(eqtlgse)
#                    ID Description setSize enrichmentScore      NES    pvalue  p.adjust qvalue rank
# GTEx_eQTLs GTEx_eQTLs  GTEx_eQTLs     262       0.1552173 1.317816 0.0969031 0.0969031     NA 4764
# leading_edge
# GTEx_eQTLs tags=58%, list=50%, signal=30%

### DEG enrichment ----
common_symbol_de <- common_symbol_de[order(common_symbol_de$P.value.Studies, decreasing = T),]
deg_score <- common_symbol_de$P.value.Studies
names(deg_score) <- common_symbol_statistics$Symbol
eqtlgse <- clusterProfiler::GSEA(geneList = deg_score, 
                                 TERM2GENE = data.frame(term = "GTEx_eQTLs", gene = gtex_eqtls$gene_name),
                                 eps = 1e-300, pvalueCutoff = 1,
                                 minGSSize = 10, maxGSSize = 500)
gseaplot(eqtlgse, geneSetID = 1, title = "GTEx eQTL enrichment by study called DE (nom. p)")
as.data.frame(eqtlgse)
# ID Description setSize enrichmentScore      NES    pvalue  p.adjust qvalue rank
# GTEx_eQTLs GTEx_eQTLs  GTEx_eQTLs     262       0.1552173 1.317816 0.0969031 0.0969031     NA 4764
# leading_edge
# GTEx_eQTLs tags=58%, list=50%, signal=30%

common_symbol_de <- common_symbol_de[order(common_symbol_de$Adjusted.P.value.Studies, decreasing = T),]
deg_score <- common_symbol_de$Adjusted.P.value.Studies
names(deg_score) <- common_symbol_statistics$Symbol
eqtlgse <- clusterProfiler::GSEA(geneList = deg_score, 
                                 TERM2GENE = data.frame(term = "GTEx_eQTLs", gene = gtex_eqtls$gene_name),
                                 eps = 1e-300, pvalueCutoff = 1,
                                 minGSSize = 10, maxGSSize = 500)
gseaplot(eqtlgse, geneSetID = 1, title = "GTEx eQTL enrichment by study called DE (adj. p)")
as.data.frame(eqtlgse)
# ID Description setSize enrichmentScore      NES    pvalue  p.adjust qvalue rank
# GTEx_eQTLs GTEx_eQTLs  GTEx_eQTLs     262       0.1552173 1.317816 0.0969031 0.0969031     NA 4764
# leading_edge
# GTEx_eQTLs tags=58%, list=50%, signal=30%

## Housekeeping gene enrichment ----
hkg <- read.csv("./Box organization/housekeeping gene/housekeeping_scores.csv")
hkg <- hkg[order(hkg$Score, decreasing = T),]
hkg <- hkg[which(hkg$Symbols %in% common_symbols),]
dim(hkg)
# [1] 9511   15
hkgs <- hkg$Score
names(hkgs) <- hkg$Symbols
eqtlgse <- clusterProfiler::GSEA(geneList = hkgs, 
                                 TERM2GENE = data.frame(term = "GTEx_eQTLs", gene = gtex_eqtls$gene_name),
                                 eps = 1e-300, pvalueCutoff = 1,
                                 minGSSize = 10, maxGSSize = 500)
gseaplot(eqtlgse, geneSetID = 1, title = "GTEx eQTL enrichment by housekeeping gene score")
as.data.frame(eqtlgse)
#                    ID Description setSize enrichmentScore      NES    pvalue  p.adjust qvalue  rank
# GTEx_eQTLs GTEx_eQTLs  GTEx_eQTLs     415        0.352386 1.056264 0.1498501 0.1498501     NA 15234
# leading_edge
# GTEx_eQTLs tags=62%, list=53%, signal=30%

# Figure 4B ----
fourb <- tibble(Statistic = c("Repeatability", "Subject variancePartition", "Within variation", "Total variation",
                     "Rs statistic", "Genetic variance", "Average expression", "Stability score", 
                     "Time variancePartition", "Flexibility score", "Housekeeping score"), 
       `P-value` = c(4.102507e-40, 1.227463e-41, 0.365, 1.365526e-26, 2.629541e-32, 1.157492e-17, 0.0969031,
                     1.252482e-32, 0.0969031, 0.0969031, 0.1498501))
fourb <- fourb[order(fourb$`P-value`),]
ggplot(fourb, aes(x = factor(Statistic, levels = Statistic), y = -log10(`P-value`))) +
  geom_col(fill = "purple3", alpha = 0.8, color = "grey40") +
  theme_minimal_hgrid() +
  theme(axis.text.x = element_text(size = 7.5, angle = 45, vjust = 0.5)) +
  labs(x = "Statistic", y = "-log10 p-value")
