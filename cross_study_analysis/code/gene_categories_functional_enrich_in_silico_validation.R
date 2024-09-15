## Jul 21 2024
## New functional enrichment: don't do row by row, but score by score correlations.
## Also, only analyze genes in the common symbols.

## Maximizing the difference between type 2 genes (ideal diagnostic / prognostic / risk biomarkers)
## and temporally flexible / dynamic genes (ideal safety / pharmacodynamic biomarkers)
library(tidyverse)
library(clusterProfiler)
library(VennDiagram)
library(biomaRt)
library(org.Hs.eg.db)

## INPUT: CSV files for housekeeping, flexible, and stable-polymorphic genes stored in the cross-study analysis directory
input_dir <- "~/Desktop/work_repo/github/cross_study_analysis/output/"
## OUTPUT: Venn diagrams and functional enrichment figures and tables
output_dir <- "~/Desktop/work_repo/github/cross_study_analysis/output/"

# First, load in all data ----
characteristic <- read.csv(file = paste(input_dir, "/characteristic_scores.csv", sep = ""))
## Contents: the .sd, .vp, .rp analysis from July 1st, described in the Methods section
head(characteristic)
dim(characteristic)
# [1] 28913    13
flexible <- read.csv(file = paste(input_dir, "/flexible_gene_scores.csv", sep = ""))
head(flexible)
dim(flexible)
# [1] 28974    16

hkg <- read.csv(file = paste(input_dir, "/housekeeping_scores.csv", sep = ""))
dim(hkg)
# [1] 28912    18

symbols <- read.table(file = paste(input_dir, "common_symbols9474.txt", sep = ""))
symbols <- symbols$x
characteristic <- characteristic[which(characteristic$Symbol %in% symbols),]
dim(characteristic)
# [1] 9474   13
flexible <- flexible[which(flexible$Symbol %in% symbols),]
hkg <- hkg[which(hkg$Symbol %in% symbols),]

# Count the number of microarray and RNA-seq ----
characteristic <- characteristic %>% 
  mutate(RNAseq = gosch + gomez + larocca) %>% 
  mutate(Microarray = obermoser + meaburn + dusek + rusch + karlovich)

flexible <- flexible %>% 
  mutate(RNAseq = gomez + gosch + larocca) %>% 
  mutate(Microarray = meaburn + obermoser + dusek + rusch + karlovich)

hkg <- hkg %>% 
  mutate(RNAseq = gomez + gosch + larocca) %>% 
  mutate(Microarray = obermoser1 + obermoser2 + obermoser3 + obermoser4 + dusek + rusch + meaburn1 + meaburn2
         + karlovich1 + karlovich2)

# Create enrichment categories ----
stable_enrich <- characteristic %>% 
  filter(Study_counts >= 4, Microarray > 0, RNAseq > 0) %>% 
  dplyr::select(Symbol) %>% 
  unlist()
length(stable_enrich)
# [1] 3910

dynamic_enrich <- flexible %>% 
  filter(P_value_study_count >= 4, Microarray > 0, RNAseq > 0) %>% 
  dplyr::select(Symbol) %>% 
  unlist()
length(dynamic_enrich)
# [1] 2855

house_enrich <- hkg %>% 
  filter(Studies >= 4, Microarray > 0, RNAseq > 0) %>% 
  dplyr::select(Symbol) %>% 
  unlist()
length(house_enrich)
# [1] 289

## Create the Venn Diagram before removing intersections ----
venn.diagram(x = list(stable_enrich, dynamic_enrich, house_enrich), 
             category.names = c("Characteristic", "Flexible", "Housekeeping"),
             fill = c("#84cbeb", "#fffeb0", "#f6c6ad"), col = "black", alpha = c(0.8,0.8,0.8),
             main = "Gene categories", imagetype = "png", cex = 1.5,
             main.fontfamily = "Arial", main.fontface = "bold", main.cex = 1.5, fontfamily = "Arial",
             cat.fontfamily = "Arial", cat.fontface = "bold", cat.dist = c(-0.075,-0.05, -0.05), cat.cex = c(1.1,1.2,1.2),
             disable.logging = T,
             filename = paste(output_dir, "venn_diagram_overlap.png", sep = ""))

## Remove intersections and create new Venn diagram ----
house_enrich2 <- setdiff(house_enrich, stable_enrich) %>% setdiff(y = dynamic_enrich)
length(house_enrich2) # 188
stable_enrich2 <- setdiff(stable_enrich, house_enrich) %>% setdiff(y = dynamic_enrich)
length(stable_enrich2) # 2768
dynamic_enrich2 <- setdiff(dynamic_enrich, house_enrich) %>% setdiff(y = stable_enrich)
length(dynamic_enrich2) # 1668
write.table(stable_enrich2, paste(output_dir, "characteristic_gene_list.txt", sep = ""))
write.table(dynamic_enrich2, paste(output_dir, "flexible_gene_list.txt", sep = ""))
write.table(house_enrich2, paste(output_dir, "housekeeping_gene_list.txt", sep = ""))

# In silico validation of gene sets ----
### The Eisenberg et al. study identified over 3800 gene symbols they classified as housekeeping genes in 16 human
### tissues. How many genes in our enrichment category fall into this?
eisenberg <- read.table("https://www.tau.ac.il/~elieis/HKG/HK_genes.txt")
head(eisenberg)
#      V1        V2
# 1  AAAS NM_015665
# 2 AAGAB NM_024666
# 3  AAMP NM_001087
# 4  AAR2 NM_015511
# 5  AARS NM_001605
# 6 AARS2 NM_020745
## Store for future use
dir.create(paste(output_dir, "in_silico_validation", sep= ""))
write.table(eisenberg, file = paste(input_dir, "in_silico_validation/eisenberg_housekeeping_genes_August_27_2024.txt", sep = ""))
# Over-enrichment test
phyper(q = length(intersect(house_enrich2, eisenberg$V1)) - 1, # Number of genes - 1, exclusive greater than
       m = length(intersect(eisenberg$V1, symbols)),
       n = length(setdiff(symbols, eisenberg$V1)),
       k = length(house_enrich2),
       lower.tail = F)
# [1] 9.664911e-18
### Figure for paper (crop title)
venn.diagram(x = list(house_enrich2, intersect(eisenberg$V1, symbols)), 
             category.names = c("Housekeeping set", "Eisenberg et al."),
             fill = c("#f6c6ad", "#fba200"), col = "black", alpha = c(0.9, 0.9),
             main = "Previously identified housekeeping", imagetype = "png", cex = 1.3,
             main.fontfamily = "Arial", main.fontface = "bold", fontfamily = "Arial",
             cat.fontfamily = "Arial", cat.dist = c(0.1,-0.08), cat.pos = c(-90, 170), cat.cex = c(1.15,1.17),
             ext.dist = c(0.0005,0.008), ext.length = c(0.9,0.8), margin = 0.16,
             disable.logging = T, force.unique = T,
             filename = paste(output_dir, "in_silico_validation/housekeeping_eisenberg_validation.png", sep = ""))

## Dynamic and DeBoever seasonal gene enrichment ----
### DeBoever et al. performed a similar analysis to us. They created a list of genes 'differentially expressed' with
### time 'within-season' and 'between-season', or short and long-term. 
## Download online at https://doi.org/10.1016/j.ygeno.2013.11.006 Supp. Table 3 and Supp. Table 4
input_dir <- "~/Desktop/work_repo/data/"
deboever_within <- readxl::read_xlsx(paste(input_dir, "1-s2.0-S088875431300222X-mmc3.xlsx", sep = "")) ## Table S3
deboever_between <- readxl::read_xlsx(paste(input_dir, "1-s2.0-S088875431300222X-mmc4.xlsx", sep = "")) ## Table S4
deboever <- na.omit(c(deboever_between$GeneSymbol, deboever_within$GeneSymbol))
length(unique(deboever))
# [1] 1495
phyper(q = length(intersect(dynamic_enrich2, unique(deboever))) - 1,
       m = length(intersect(unique(deboever), symbols)),
       n = length(setdiff(symbols, unique(deboever))),
       k = length(dynamic_enrich2),
       lower.tail = F)
# [1] 1.083316e-06
venn.diagram(x = list(dynamic_enrich2, intersect(deboever, symbols)), 
             category.names = c("Flexible set", "DeBoever et al."),
             fill = c("#fffeb0", "yellow4"), col = "black", cex = 1.4, alpha = c(0.8,0.8),
             main = "Previously identified flexible genes", imagetype = "png", 
             main.fontfamily = "Arial", main.fontface = "bold", fontfamily = "Arial",
             cat.fontfamily = "Arial", cat.dist = c(-0.056,-0.07), cat.cex = c(1.2,1.2),
             disable.logging = T, force.unique = T,
             filename = paste(output_dir, "in_silico_validation/flexible_deboever_validation.png", sep = ""))

## Stable-polymorphic and GTEx eQTL enrichment ----
gtex <- read_tsv("~/Desktop/work_repo/data/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.egenes.txt.gz")
gtex <- gtex %>% 
  dplyr::filter(qval < 0.01)
length(gtex$gene_name)
# [1] 10071
length(intersect(gtex$gene_name, symbols))
# [1] 5333
## lintersect
lintersect <- function(x,y) {
  length(intersect(unique(x),unique(y)))
}
phyper(q = lintersect(stable_enrich2, gtex$gene_name) - 1,
       m = lintersect(gtex$gene_name, symbols),
       n = length(setdiff(symbols, gtex$gene_name)),
       k = length(stable_enrich2), 
       lower.tail = F)
# [1] 5.497379e-32
venn.diagram(x = list(stable_enrich2, intersect(gtex$gene_name, symbols)), 
             category.names = c("Characteristic", "GTEx eGenes (q < 0.01)"),
             fill = c("#84cbeb", "#005fbf"), col = "black", cex = 1.4, alpha = c(0.8,0.8),
             main = "Previously identified GTEx eQTLs", imagetype = "png", 
             main.fontfamily = "Arial", main.fontface = "bold", fontfamily = "Arial",
             cat.fontfamily = "Arial", cat.dist = c(0.035,-0.08), cat.pos = c(200,45), cat.cex = c(1.3, 1.2),
             disable.logging = T, force.unique = T,
             filename = paste(output_dir, "in_silico_validation/characteristic_GTEx_validation.png", sep = ""))


## Checking for under-enrichment ----
## Lower tail is true
## You can see expected values at https://systems.crump.ucla.edu/hypergeometric/index.php 
## Stable-polymorphic with Eisenberg
phyper(q = lintersect(stable_enrich2, eisenberg$V1), #652, inclusive less than
       m = lintersect(eisenberg$V1, symbols), #3175
       n = length(setdiff(symbols, eisenberg$V1)), #6299
       k = length(stable_enrich2), #2768
       lower.tail = T) 
# [1] 2.368478e-41, under-enriched, 927 expected
supp_figure_4 <- data.frame(expected = 927, actual = 652)
a <- ggplot(pivot_longer(supp_figure_4, cols = c('expected', 'actual')), aes(x = factor(name, levels = c('expected', 'actual')), y = value, fill = name)) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("#84cbeb","#fba200")) +
  ylim(c(0,1000)) +
  labs(x = "", y = "") +
  theme_minimal_hgrid() +
  theme(legend.position = 'none') +
  annotate(geom = 'text', label = '****', x = 2, y = 660, size = 15)
## Dynamic with Eisenberg
phyper(q = lintersect(dynamic_enrich2, eisenberg$V1), #656
       m = lintersect(eisenberg$V1, symbols), #3175
       n = length(setdiff(symbols, eisenberg$V1)), #6299
       k = length(dynamic_enrich2), #1668
       lower.tail = T)
# [1] 1 , over-enriched, 558 expected
phyper(q = lintersect(dynamic_enrich2, eisenberg$V1) - 1, #656
       m = lintersect(eisenberg$V1, symbols), #3175
       n = length(setdiff(symbols, eisenberg$V1)), #6299
       k = length(dynamic_enrich2), #1668
       lower.tail = F)
# [1] 2.347654e-08 over-enriched p-value
supp_figure_4 <- data.frame(expected = 558, actual = 656)
b <- ggplot(pivot_longer(supp_figure_4, cols = c('expected', 'actual')), aes(x = factor(name, levels = c('expected', 'actual')), y = value, fill = name)) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("#fffeb0", "#fba200")) +
  ylim(c(0,1000)) +
  labs(x = "", y = "") +
  theme_minimal_hgrid() +
  theme(legend.position = 'none') +
  annotate(geom = 'text', label = '****', x = 2, y = 660, size = 15)
pdf(file = paste(output_dir, "/in_silico_validation/figure_4c_bottom.pdf", sep = ""), width = 9, height = 6)
plot_grid(a,b, labels = c('Characteristic','Flexible'), nrow = 1)
dev.off()

## Stable-polymorphic with De Boever
phyper(q = lintersect(stable_enrich2, deboever), #247
       m = lintersect(deboever, symbols), #1019
       n = length(setdiff(symbols, deboever)), #8455
       k = length(stable_enrich2), #2768
       lower.tail = T)
# [1] 0.000100876, expected 297 overlaps
supp_figure_4 <- data.frame(expected = 297, actual = 247)
c <- ggplot(pivot_longer(supp_figure_4, cols = c('expected', 'actual')), aes(x = factor(name, levels = c('expected', 'actual')), y = value, fill = name)) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("#84cbeb","yellow4")) +
  ylim(c(0,310)) +
  labs(x = "", y = "") +
  theme_minimal_hgrid() +
  theme(legend.position = 'none') +
  annotate(geom = 'text', label = '***', x = 2, y = 255, size = 15)
## Housekeeping with De Boever
phyper(q = lintersect(house_enrich2, deboever), #26
       m = lintersect(deboever, symbols), #1019
       n = length(setdiff(symbols, deboever)), #8455
       k = length(house_enrich2), #188
       lower.tail = T)
# [1] 0.9281765 over-enriched, 20 expected
phyper(q = lintersect(house_enrich2, deboever)-1, #26
       m = lintersect(deboever, symbols), #1019
       n = length(setdiff(symbols, deboever)), #8455
       k = length(house_enrich2), #188
       lower.tail = F)
# [1] 0.1072417 p-value
supp_figure_4 <- data.frame(expected = 20, actual = 26)
d <- ggplot(pivot_longer(supp_figure_4, cols = c('expected', 'actual')), aes(x = factor(name, levels = c('expected', 'actual')), y = value, fill = name)) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("#f6c6ad","yellow4")) +
  ylim(c(0,310)) +
  labs(x = "", y = "") +
  theme_minimal_hgrid() +
  theme(legend.position = 'none')
pdf(file = paste(output_dir, "/in_silico_validation/figure_4b_bottom.pdf", sep = ""), width = 9, height = 6)
plot_grid(c,d, labels = c('Characteristic','Housekeeping'), nrow = 1)
dev.off()

## Dynamic with GTEx
phyper(q = lintersect(dynamic_enrich2, gtex$gene_name), #777
       m = lintersect(gtex$gene_name, symbols), #5333
       n = length(setdiff(symbols, gtex$gene_name)), #4141
       k = length(dynamic_enrich2), #1668
       lower.tail = T)
# [1] 1.105556e-18, 938 expected
supp_figure_4 <- data.frame(expected = 938, actual = 777)
e <- ggplot(pivot_longer(supp_figure_4, cols = c('expected', 'actual')), aes(x = factor(name, levels = c('expected', 'actual')), y = value, fill = name)) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("#fffeb0","#005fbf")) +
  ylim(c(0,1000)) +
  labs(x = "", y = "") +
  theme_minimal_hgrid() +
  theme(legend.position = 'none') +
  annotate(geom = 'text', label = '****', x = 2, y = 785, size = 15)
## Housekeeping with GTEx
phyper(q = lintersect(house_enrich2, gtex$gene_name), #89
       m = lintersect(gtex$gene_name, symbols), #5333
       n = length(setdiff(symbols, gtex$gene_name)), #4141
       k = length(house_enrich2), #188
       lower.tail = T)
# [1] 0.007860481, 105 expected
supp_figure_4 <- data.frame(expected = 105, actual = 89)
f <- ggplot(pivot_longer(supp_figure_4, cols = c('expected', 'actual')), aes(x = factor(name, levels = c('expected', 'actual')), y = value, fill = name)) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("#f6c6ad","#005fbf")) +
  ylim(c(0,1000)) +
  labs(x = "", y = "") +
  theme_minimal_hgrid() +
  theme(legend.position = 'none') +
  annotate(geom = 'text', label = '**', x = 2, y = 95, size = 15)
pdf(file = paste(output_dir, "/in_silico_validation/figure_4a_bottom.pdf", sep = ""), width = 9, height = 6)
plot_grid(e,f, labels = c('Flexible','Housekeeping'), nrow = 1)
dev.off()

# Functional enrichment analysis ----
## Get entrez IDs for enrichKEGG
categories <- list(stable_enrich2, dynamic_enrich2, house_enrich2)
names(categories) <- c("characteristic", "flexible", "housekeeping")
## biomaRt
mart <- useDataset(dataset = "hsapiens_gene_ensembl", mart = useMart("ENSEMBL_MART_ENSEMBL"))
categories_entrez <- lapply(categories, function(x) {
  getBM(attributes = c("external_gene_name","entrezgene_id"),
        filters = "external_gene_name",
        values = x,
        mart = mart)
})
## Universe of entrez IDs
entrez_universe <- getBM(attributes = c("external_gene_name","entrezgene_id"),
                         filters = "external_gene_name",
                         values = symbols,
                         mart = mart)
nrow(entrez_universe)
# [1] 9562

resultsGO <- list()
resultsKEGG <- list()

dir.create(paste(output_dir, "functional_enrichment", sep = ""))

## Perform enrichment loop ----
for (i in 1:3) {
  resGO <- enrichGO(categories[[i]],
                    OrgDb = org.Hs.eg.db,
                    keyType = "SYMBOL",
                    ont = "ALL", 
                    universe = symbols, 
                    pvalueCutoff = 1,
                    qvalueCutoff = 1,
                    pool = FALSE)
  rdf <- as.data.frame(resGO)
  rdf <- rdf[order(rdf$pvalue),]
  resultsGO[[i]] <- rdf
  write.csv(rdf, paste(output_dir, "functional_enrichment/",names(categories)[i],"_GO.csv", sep = ""))
  resKEGG <- enrichKEGG(as.character(categories_entrez[[i]]$entrezgene_id), 
                        keyType = "ncbi-geneid",
                        maxGSSize = 500, 
                        pvalueCutoff = 1,
                        qvalueCutoff = 1,
                        universe = as.character(entrez_universe$entrezgene_id))
  rkf <- as.data.frame(resKEGG)
  rkf <- rkf[order(rkf$pvalue),]
  resultsKEGG[[i]] <- rkf
  write.csv(rkf, paste(output_dir, "functional_enrichment/",names(categories)[i],"_KEGG.csv", sep = ""))
}

## Plot results ---
library(cowplot)
go_bp <- list()
go_cc <- list()
go_mf <- list()
kegg <- list()


resultsGO[[1]]$passed <- 0
resultsGO[[1]]$passed <- resultsGO[[1]]$p.adjust < 0.05
resultsGO[[2]]$passed <- 0
resultsGO[[2]]$passed <- resultsGO[[2]]$p.adjust < 0.05
resultsGO[[3]]$passed <- 0
resultsGO[[3]]$passed <- resultsGO[[3]]$p.adjust < 0.05
resultsKEGG[[1]]$passed <- 0
resultsKEGG[[1]]$passed <- resultsKEGG[[1]]$p.adjust < 0.05
resultsKEGG[[2]]$passed <- 0
resultsKEGG[[2]]$passed <- resultsKEGG[[2]]$p.adjust < 0.05
resultsKEGG[[3]]$passed <- 0
resultsKEGG[[3]]$passed <- resultsKEGG[[3]]$p.adjust < 0.05

## Fixing some of the descriptions with longer names

## Flexible genes
resultsGO[[2]]$Description[12]
# "RNA splicing, via transesterification reactions with bulged adenosine as nucleophile"
resultsGO[[2]]$Description[12] <- "RNA splicing, via transesterification reactions*"
resultsGO[[2]]$Description[38]
# "basal RNA polymerase II transcription machinery binding"
resultsGO[[2]]$Description[38] <- "basal RNA pol II transcription machinery binding"
resultsKEGG[[2]]$Description[6]
# "Parathyroid hormone synthesis, secretion and action"
resultsKEGG[[2]]$Description[6] <- "PTH synthesis, secretion and action"

## Houskeeping genes
resultsGO[[3]]$Description[1]
# "RNA splicing, via transesterification reactions with bulged adenosine as nucleophile"
resultsGO[[3]]$Description[1] <- "RNA splicing, via transesterification reactions*"
resultsKEGG[[3]]$Description[2]
# "Glycosylphosphatidylinositol (GPI)-anchor biosynthesis"
resultsKEGG[[3]]$Description[2] <- "GPI-anchor biosynthesis"

## Generate plots
colors <- list(sp = c("FALSE" = "#84cbeb", "TRUE" = "#005fbf"), fl = c("FALSE" = "#fffeb0", "TRUE" = "yellow4"), ho = c("FALSE" = "#f6c6ad", "TRUE" = "#fba200"))
for (i in 1:3) {
  lims <- c(0,max(max(-log10(resultsGO[[i]]$pvalue)),
                  max(-log10(resultsKEGG[[i]]$pvalue)))
  )
  ## GO_BP
  go_bp[[i]] <- ggplot(na.omit(head(resultsGO[[i]][which(resultsGO[[i]]$ONTOLOGY == "BP"),],10)), aes(x = factor(Description, levels = rev(Description)), y = -log10(pvalue), fill = passed)) +
    geom_col(color = "black", width = 0.5, show.legend = F) +
    scale_fill_manual(values = colors[[i]]) +
    theme_cowplot() +
    coord_flip() +
    labs(y = "", x = "",
         title =  "") +
    ylim(lims) +
    theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16))
  ## GO_CC
  go_cc[[i]] <- ggplot(na.omit(head(resultsGO[[i]][which(resultsGO[[i]]$ONTOLOGY == "CC"),],10)), aes(x = factor(Description, levels = rev(Description)), y = -log10(pvalue), fill = passed)) +
    geom_col(color = "black", width = 0.5, show.legend = F) +
    scale_fill_manual(values = colors[[i]]) +
    theme_cowplot() +
    coord_flip() +
    labs(y = "", x = "",
         title =  "") +
    ylim(lims) +
    theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16))
  ## GO_MF
  go_mf[[i]] <- ggplot(na.omit(head(resultsGO[[i]][which(resultsGO[[i]]$ONTOLOGY == "MF"),],10)), aes(x = factor(Description, levels = rev(Description)), y = -log10(pvalue), fill = passed)) +
    geom_col(color = "black", width = 0.5, show.legend = F) +
    scale_fill_manual(values = colors[[i]]) +
    theme_cowplot() +
    coord_flip() +
    labs(y = "", x = "",
         title = "") +
    ylim(lims) +
    theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16))
  ## KEGG
  kegg[[i]] <- ggplot(na.omit(head(resultsKEGG[[i]], 10)), aes(x = factor(Description, levels = rev(Description)), y = -log10(pvalue), fill = passed)) +
    geom_col(color = "black", width = 0.5, show.legend = F) +
    scale_fill_manual(values = colors[[i]]) +
    theme_cowplot() +
    coord_flip() +
    labs(y = "-log10(p-value)", x = "",
         title =  "") +
    ylim(lims) +
    theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16))
  
}

## Stable-polymorphic function
a <- plot_grid(go_bp[[1]], go_cc[[1]], go_mf[[1]], kegg[[1]], ncol = 1, align = 'v')
## FLexible function
b <- plot_grid(go_bp[[2]], go_cc[[2]], go_mf[[2]], kegg[[2]], ncol = 1, align = 'v')
## Housekeeping function
c <- plot_grid(go_bp[[3]], go_cc[[3]], go_mf[[3]], kegg[[3]], ncol = 1, align = 'v')

pdf(file = paste(output_dir, "functional_enrichment/functional_enrichment_figure.pdf", sep = ""), height = 12, width = 23)
plot_grid(a, b, c, ncol = 3, nrow = 1, labels = c('A. Characteristic genes', 'B. Flexible genes', 'C. Housekeeping genes'),
          align = 'hv')
dev.off()

# Supplementary data file 2 ----
sp <- read.table(paste(output_dir, "characteristic_gene_list.txt", sep = ""))
fl <- read.table(paste(output_dir, "flexible_gene_list.txt", sep = ""))
hk <- read.table(paste(output_dir, "housekeeping_gene_list.txt", sep = ""))
sdf4 <- list(sp, fl, hk)
names(sdf4) <- c("characteristic", "flexible", "housekeeping")
dataset <- useDataset(dataset = 'hsapiens_gene_ensembl', mart = useMart('ENSEMBL_MART_ENSEMBL'))
gene_names <- c()
sdf4 <- lapply(sdf4, function(x) {
  gene_names <- getBM(attributes = c('hgnc_symbol','description'),
                      filters = 'hgnc_symbol',
                      values = x[[1]], 
                      mart = dataset)
  x <- gene_names
})
lapply(sdf4, head)
## First list
# $characteristic
# hgnc_symbol                                                                description
# 1        ABAT          4-aminobutyrate aminotransferase [Source:HGNC Symbol;Acc:HGNC:23]
# 2       ABCA1 ATP binding cassette subfamily A member 1 [Source:HGNC Symbol;Acc:HGNC:29]
# 3       ABCA2 ATP binding cassette subfamily A member 2 [Source:HGNC Symbol;Acc:HGNC:32]
# 4       ABCA5 ATP binding cassette subfamily A member 5 [Source:HGNC Symbol;Acc:HGNC:35]
# 5       ABCA7 ATP binding cassette subfamily A member 7 [Source:HGNC Symbol;Acc:HGNC:37]
# 6       ABCB1 ATP binding cassette subfamily B member 1 [Source:HGNC Symbol;Acc:HGNC:40]
writexl::write_xlsx(x = sdf4, path = paste(output_dir, "supplementary_data_2.xlsx", sep = ""))