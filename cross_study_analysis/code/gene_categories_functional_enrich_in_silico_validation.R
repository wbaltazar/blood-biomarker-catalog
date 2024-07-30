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
stable_polymorphic <- read.csv(file = paste(input_dir, "/stable_polymorphic_scores.csv", sep = ""))
## Contents: the .sd, .vp, .rp analysis from July 1st, described in the Methods section
head(stable_polymorphic)
dim(stable_polymorphic)
# [1] 28913    12
flexible <- read.csv(file = paste(input_dir, "/flexible_gene_scores.csv", sep = ""))
head(flexible)
dim(flexible)
# [1] 27140    15

hkg <- read.csv(file = paste(input_dir, "/housekeeping_scores.csv", sep = ""))
dim(hkg)
# [1] 28912    16

symbols <- read.table(file = paste(input_dir, "common_symbols9474.txt", sep = ""))
symbols <- symbols$x
stable_polymorphic <- stable_polymorphic[which(stable_polymorphic$Symbol %in% symbols),]
dim(stable_polymorphic)
# [1] 9474   12
flexible <- flexible[which(flexible$Symbol %in% symbols),]
hkg <- hkg[which(hkg$Symbol %in% symbols),]

# Count the number of microarray and RNA-seq ----
stable_polymorphic <- stable_polymorphic %>% 
  mutate(RNAseq = gosch + gomez + larocca) %>% 
  mutate(Microarray = obermoser + meaburn + dusek + rusch)

flexible <- flexible %>% 
  mutate(RNAseq = gomez + gosch + larocca) %>% 
  mutate(Microarray = meaburn + obermoser + dusek + rusch)

hkg <- hkg %>% 
  mutate(RNAseq = gomez + gosch + larocca) %>% 
  mutate(Microarray = obermoser1 + obermoser2 + obermoser3 + obermoser4 + dusek + rusch + meaburn1 + meaburn2)

# Create enrichment categories ----
stable_enrich <- stable_polymorphic %>% 
  filter(Study_counts >= 4, Microarray > 0, RNAseq > 0) %>% 
  dplyr::select(Symbol) %>% 
  unlist()
length(stable_enrich)
# [1] 3225

dynamic_enrich <- flexible %>% 
  filter(P_value_study_count >= 4, Microarray > 0, RNAseq > 0) %>% 
  dplyr::select(Symbol) %>% 
  unlist()
length(dynamic_enrich)
# [1] 1739

house_enrich <- hkg %>% 
  filter(Studies >= 4, Microarray > 0, RNAseq > 0) %>% 
  dplyr::select(Symbol) %>% 
  unlist()
length(house_enrich)
# [1] 222

## Create the Venn Diagram before removing intersections ----
venn.diagram(x = list(stable_enrich, dynamic_enrich, house_enrich), 
             category.names = c("Stable-polymorphic", "Flexible", "Housekeeping"),
             fill = c("lightblue", "pink", "lightgreen"), col = "black",
             main = "Gene categories", imagetype = "png", 
             main.fontfamily = "Arial", main.fontface = "bold", fontfamily = "Arial",
             cat.fontfamily = "Arial", cat.dist = c(-0.075,-0.05, -0.05),
             disable.logging = T,
             filename = paste(output_dir, "venn_diagram_overlap.png", sep = ""))

## Remove intersections and create new Venn diagram ----
house_enrich2 <- setdiff(house_enrich, stable_enrich) %>% setdiff(y = dynamic_enrich)
length(house_enrich2) # 175
stable_enrich2 <- setdiff(stable_enrich, house_enrich) %>% setdiff(y = dynamic_enrich)
length(stable_enrich2) # 2688
dynamic_enrich2 <- setdiff(dynamic_enrich, house_enrich) %>% setdiff(y = stable_enrich)
length(dynamic_enrich2) # 1179
write.table(stable_enrich2, paste(output_dir, "stable_polymorphic_gene_list.txt", sep = ""))
write.table(dynamic_enrich2, paste(output_dir, "flexible_gene_list.txt", sep = ""))
write.table(house_enrich2, paste(output_dir, "housekeeping_gene_list.txt", sep = ""))

venn.diagram(x = list(stable_enrich2, dynamic_enrich2, house_enrich2), 
             category.names = c("Stable-polymorphic", "Flexible", "Housekeeping"),
             fill = c("lightblue", "pink", "lightgreen"), col = "black",
             main = "Functional enrichment categories", imagetype = "png", 
             main.fontfamily = "Arial", main.fontface = "bold", fontfamily = "Arial",
             cat.fontfamily = "Arial", cat.dist = c(-0.1, -0.06, 0.017),
             disable.logging = T,
             filename = paste(output_dir, "genes_used_for_func_enrich.png", sep = ""))

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
write.table(eisenberg, file = paste(input_dir, "in_silico_validation/eisenberg_housekeeping_genes_July_21_2024.txt", sep = ""))
# Over-enrichment test
phyper(q = length(intersect(house_enrich2, eisenberg$V1)) - 1,
       m = length(intersect(eisenberg$V1, symbols)),
       n = length(setdiff(symbols, eisenberg$V1)),
       k = length(house_enrich2),
       lower.tail = F)
# [1] 6.52438e-18
### Figure for paper
venn.diagram(x = list(house_enrich2, intersect(eisenberg$V1, symbols)), 
             category.names = c("Housekeeping set", "Eisenberg et al."),
             fill = c("grey80", "grey50"), col = "black",
             main = "Previously identified housekeeping", imagetype = "png", 
             main.fontfamily = "Arial", main.fontface = "bold", fontfamily = "Arial",
             cat.fontfamily = "Arial", cat.dist = c(0.015,-0.08), cat.pos = c(0, 90),
             disable.logging = T, hyper.test = T, lower.tail = F, total.population = length(symbols), force.unique = T,
             filename = paste(output_dir, "in_silico_validation/housekeeping_eisenberg_validation.png", sep = ""))

## Dynamic and DeBoever seasonal gene enrichment ----
### DeBoever et al. performed a similar analysis to us. They created a list of genes 'differentially expressed' with
### time 'within-season' and 'between-season', or short and long-term. 
## Download online at https://doi.org/10.1016/j.ygeno.2013.11.006 Supp. Table 3 and Supp. Table 5
input_dir <- "~/Desktop/work_repo/data/"
deboever_within <- readxl::read_xlsx(paste(input_dir, "DeBoever short-term.xlsx", sep = "")) ## Table S3
deboever_between <- readxl::read_xlsx(paste(input_dir, "DeBoever long-term.xlsx", sep = "")) ## Table S5
deboever <- c(deboever_between$`DeBoever long-term`, deboever_within$`DeBoever within-season (short-term)`)
length(unique(deboever))
# [1] 1495
phyper(q = length(intersect(dynamic_enrich2, unique(deboever))) - 1,
       m = length(intersect(unique(deboever), symbols)),
       n = length(setdiff(symbols, unique(deboever))),
       k = length(dynamic_enrich2),
       lower.tail = F)
# [1] 2.054972e-06
venn.diagram(x = list(dynamic_enrich2, intersect(deboever, symbols)), 
             category.names = c("Flexible set", "DeBoever et al."),
             fill = c("grey80", "grey50"), col = "black",
             main = "Previously identified flexible genes", imagetype = "png", 
             main.fontfamily = "Arial", main.fontface = "bold", fontfamily = "Arial",
             cat.fontfamily = "Arial", cat.dist = c(-0.05,-0.07),
             disable.logging = T, hyper.test = T, lower.tail = F, total.population = 9474, force.unique = T,
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
# [1] 1.970359e-39
venn.diagram(x = list(stable_enrich2, intersect(gtex$gene_name, symbols)), 
             category.names = c("Stable-polymorphic", "GTEx eGenes (q < 0.01)"),
             fill = c("grey80", "grey50"), col = "black",
             main = "Previously identified GTEx eQTLs", imagetype = "png", 
             main.fontfamily = "Arial", main.fontface = "bold", fontfamily = "Arial",
             cat.fontfamily = "Arial", cat.dist = c(0.03,-0.075), cat.pos = c(200,45),
             disable.logging = T, hyper.test = T, lower.tail = F, total.population = 9474, force.unique = T,
             filename = paste(output_dir, "in_silico_validation/stable_poly_GTEx_validation.png", sep = ""))


## Checking for unexpected enrichment ----
## Lower tail is false because we don't expect enrichment in any particular direction
## You can see expected values at https://systems.crump.ucla.edu/hypergeometric/index.php 
## Stable-polymorphic with Eisenberg
phyper(q = lintersect(stable_enrich2, eisenberg$V1), #608
       m = lintersect(eisenberg$V1, symbols), #3175
       n = length(setdiff(symbols, eisenberg$V1)), #6299
       k = length(stable_enrich2), #2688
       lower.tail = T) # Lower and upper tails
# [1] 1.889889e-47, under-enriched
## Stable-polymorphic with De Boever
phyper(q = lintersect(stable_enrich2, deboever), #259
       m = lintersect(deboever, symbols), #1019
       n = length(setdiff(symbols, deboever)), #8455
       k = length(stable_enrich2), #2688
       lower.tail = T)
# [1] 0.01407925
## Dynamic with GTEx
phyper(q = lintersect(dynamic_enrich2, gtex$gene_name), #545
       m = lintersect(gtex$gene_name, symbols), #5333
       n = length(setdiff(symbols, gtex$gene_name)), #4141
       k = length(dynamic_enrich2), #1179
       lower.tail = T)
# [1] 7.702749e-14
## Dynamic with Eisenberg
phyper(q = lintersect(dynamic_enrich2, eisenberg$V1), #440
       m = lintersect(eisenberg$V1, symbols), #3175
       n = length(setdiff(symbols, eisenberg$V1)), #6299
       k = length(dynamic_enrich2), #1179
       lower.tail = T)
# [1] 0.9985237
phyper(q = lintersect(dynamic_enrich2, eisenberg$V1), m = lintersect(eisenberg$V1, symbols),
       n = length(setdiff(symbols, eisenberg$V1)), k = length(dynamic_enrich2), lower.tail = F)
# [1] 0.001476319
## Housekeeping with GTEx
phyper(q = lintersect(house_enrich2, gtex$gene_name), #83
       m = lintersect(gtex$gene_name, symbols), #5333
       n = length(setdiff(symbols, gtex$gene_name)), #4141
       k = length(house_enrich2), #175
       lower.tail = T)
# [1] 0.01072431
## Housekeeping with De Boever
phyper(q = lintersect(house_enrich2, deboever), #27
       m = lintersect(deboever, symbols), #1019
       n = length(setdiff(symbols, deboever)), #8455
       k = length(house_enrich2), #175
       lower.tail = T)


# Functional enrichment analysis ----
## Get entrez IDs for enrichKEGG
categories <- list(stable_enrich2, dynamic_enrich2, house_enrich2)
names(categories) <- c("stable_polymorphic", "flexible", "housekeeping")
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
resultsGO[[2]]$Description[6]
# "RNA splicing, via transesterification reactions with bulged adenosine as nucleophile"
resultsGO[[2]]$Description[6] <- "RNA splicing, via transesterification reactions*"
resultsGO[[2]]$Description[23]
# "DNA-binding transcription activator activity, RNA polymerase II-specific"
resultsGO[[2]]$Description[23] <- "DNA-binding transcription activator activity**"

## Houskeeping genes
resultsGO[[3]]$Description[6]
# "RNA splicing, via transesterification reactions with bulged adenosine as nucleophile"
resultsGO[[3]]$Description[6] <- "RNA splicing, via transesterification reactions*"
resultsGO[[3]]$Description[73]
# "hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds, in linear amides"
resultsGO[[3]]$Description[73] <- "hydrolase activity***"
resultsKEGG[[3]]$Description[2]
# "Glycosylphosphatidylinositol (GPI)-anchor biosynthesis"
resultsKEGG[[3]]$Description[2] <- "GPI-anchor biosynthesis"

## Generate plots
for (i in 1:3) {
  lims <- c(0,max(max(-log10(resultsGO[[i]]$pvalue)),
                  max(-log10(resultsKEGG[[i]]$pvalue)))
  )
  ## GO_BP
  go_bp[[i]] <- ggplot(na.omit(head(resultsGO[[i]][which(resultsGO[[i]]$ONTOLOGY == "BP"),],10)), aes(x = factor(Description, levels = rev(Description)), y = -log10(pvalue), fill = as.character(passed))) +
    geom_col(color = "black", width = 0.5, show.legend = F) +
    scale_fill_manual(values = c("FALSE" = "grey80" , "TRUE" = "grey45")) +
    theme_cowplot() +
    coord_flip() +
    labs(y = "", x = "",
         title =  "") +
    ylim(lims)
  theme(axis.text.y = element_text(size = 12), title = element_text(size = 12))
  ## GO_CC
  go_cc[[i]] <- ggplot(na.omit(head(resultsGO[[i]][which(resultsGO[[i]]$ONTOLOGY == "CC"),],10)), aes(x = factor(Description, levels = rev(Description)), y = -log10(pvalue), fill = as.character(passed))) +
    geom_col(color = "black", width = 0.5, show.legend = F) +
    scale_fill_manual(values = c("FALSE" = "grey80" , "TRUE" = "grey45")) +
    theme_cowplot() +
    coord_flip() +
    labs(y = "", x = "",
         title =  "") +
    ylim(lims) +
    theme(axis.text.y = element_text(size = 12), title = element_text(size = 12))
  ## GO_MF
  go_mf[[i]] <- ggplot(na.omit(head(resultsGO[[i]][which(resultsGO[[i]]$ONTOLOGY == "MF"),],10)), aes(x = factor(Description, levels = rev(Description)), y = -log10(pvalue), fill = as.character(passed))) +
    geom_col(color = "black", width = 0.5, show.legend = F) +
    scale_fill_manual(values = c("FALSE" = "grey80" , "TRUE" = "grey45")) +
    theme_cowplot() +
    coord_flip() +
    labs(y = "", x = "",
         title = "") +
    ylim(lims) +
    theme(axis.text.y = element_text(size = 12), title = element_text(size = 12))
  ## KEGG
  kegg[[i]] <- ggplot(na.omit(head(resultsKEGG[[i]], 10)), aes(x = factor(Description, levels = rev(Description)), y = -log10(pvalue), fill = passed)) +
    geom_col(color = "black", width = 0.5, show.legend = F) +
    scale_fill_manual(values = c("FALSE" = "grey80" , "TRUE" = "grey45")) +
    theme_cowplot() +
    coord_flip() +
    labs(y = "-log10(p-value)", x = "",
         title =  "") +
    ylim(lims) +
    theme(axis.text.y = element_text(size = 12), title = element_text(size = 12))
  
}

## Stable-polymorphic function
a <- plot_grid(go_bp[[1]], go_cc[[1]], go_mf[[1]], kegg[[1]], ncol = 1, align = 'v')
## FLexible function
b <- plot_grid(go_bp[[2]], go_cc[[2]], go_mf[[2]], kegg[[2]], ncol = 1, align = 'v')
## Housekeeping function
c <- plot_grid(go_bp[[3]], go_cc[[3]], go_mf[[3]], kegg[[3]], ncol = 1, align = 'v')

pdf(file = paste(output_dir, "functional_enrichment/functional_enrichment_figure.pdf", sep = ""), height = 12, width = 20)
plot_grid(a, b, c, ncol = 3, nrow = 1, labels = c('A. Stable-polymorphic genes', 'B. Flexible genes', 'C. Housekeeping genes'),
          align = 'hv')
dev.off()