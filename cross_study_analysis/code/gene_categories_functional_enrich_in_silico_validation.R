## Jan 18 2025
## New functional enrichment: don't do row by row, but score by score correlations.
## Also, only analyze genes in the common symbols.

## Maximizing the difference between type 2 genes (ideal diagnostic / prognostic / risk biomarkers)
## and temporally flexible / dynamic genes (ideal safety / pharmacodynamic biomarkers)
library(tidyverse)
library(clusterProfiler)
library(VennDiagram)
library(biomaRt)
library(org.Hs.eg.db)
library(cowplot)

## INPUT: CSV files for housekeeping, flexible, and stable-polymorphic genes stored in the cross-study analysis directory
input_dir <- "~/Desktop/work_repo/github/cross_study_analysis/output/"
## OUTPUT: Venn diagrams and functional enrichment figures and tables
output_dir <- "~/Desktop/work_repo/github/cross_study_analysis/output/"

# First, load in all data ----
trait <- read.csv(file = paste(input_dir, "/trait_scores.csv", sep = ""))
## Contents: the .sd, .vp, .rp analysis from July 1st, described in the Methods section
head(trait)
dim(trait)
# [1] 23883    14
state <- read.csv(file = paste(input_dir, "/state_gene_scores.csv", sep = ""))
head(state)
dim(state)
# [1] 23883    16

symbols <- read.table(file = paste(input_dir, "common_symbols6099.txt", sep = ""))
symbols <- symbols$x
trait <- trait[which(trait$Symbol %in% symbols),]
dim(trait)
# [1] 6099   14
state <- state[which(state$Symbol %in% symbols),]

# Count the number of microarray and RNA-seq ----
trait <- trait %>% 
  mutate(RNAseq = gosch + gomez + larocca) %>% 
  mutate(Microarray = obermoser + meaburn + dusek + rusch + karlovich)

state <- state %>% 
  mutate(RNAseq = gomez + gosch + larocca) %>% 
  mutate(Microarray = meaburn + obermoser + dusek + rusch + karlovich)

# Create enrichment categories ----
trait_enrich <- trait %>% 
  filter(Study_counts >= 4, Microarray > 0, RNAseq > 0) %>% 
  dplyr::select(Symbol) %>% 
  unlist()
length(trait_enrich)
# [1] 1934

state_enrich <- state %>% 
  filter(P_value_study_count >= 4, Microarray > 0, RNAseq > 0) %>% 
  dplyr::select(Symbol) %>% 
  unlist()
length(state_enrich)
# [1] 1905

## Create the Venn Diagram before removing intersections ----
venn.diagram(x = list(trait_enrich, state_enrich), 
             category.names = c("Trait", "State"),
             fill = c("#84cbeb", "#f6c6ad"), col = "black", alpha = c(0.8,0.8),
             main = "Gene categories", imagetype = "png", cex = 1.5,
             main.fontfamily = "Arial", main.fontface = "bold", main.cex = 1.5, fontfamily = "Arial",
             cat.fontfamily = "Arial", cat.fontface = "bold", cat.dist = c(-0.075,-0.05), cat.cex = c(1.1,1.2),
             disable.logging = T,
             filename = paste(output_dir, "venn_diagram_overlap.png", sep = ""))
venn.diagram(x = list(trait_enrich, state_enrich), 
             category.names = c("", ""),
             fill = c("#84cbeb", "#f6c6ad"), col = "black", alpha = c(0.8,0.8),
             main = "Gene categories", imagetype = "png", cex = 1.5,
             main.fontfamily = "Arial", main.fontface = "bold", main.cex = 1.5, fontfamily = "Arial",
             cat.fontfamily = "Arial", cat.fontface = "bold", cat.dist = c(-0.075,-0.05), cat.cex = c(1.1,1.2),
             disable.logging = T,
             filename = paste(output_dir, "venn_diagram_overlap_nolabels.png", sep = ""))

## Remove intersections and create new Venn diagram ----
trait_enrich2 <- setdiff(trait_enrich, state_enrich)
length(trait_enrich2) # 1356
state_enrich2 <- setdiff(state_enrich, trait_enrich)
length(state_enrich2) # 1327
write.table(trait_enrich2, paste(output_dir, "trait_gene_list.txt", sep = ""))
write.table(state_enrich2, paste(output_dir, "state_gene_list.txt", sep = ""))

# In silico validation of gene sets ----
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
phyper(q = length(intersect(state_enrich2, unique(deboever))) - 1, #186 - 1
       m = length(intersect(unique(deboever), symbols)), #681
       n = length(setdiff(symbols, unique(deboever))), #5418
       k = length(state_enrich2), #1327
       lower.tail = F)
# [1] 0.0001598957 (148 expected)
venn.diagram(x = list(state_enrich2, intersect(deboever, symbols)), 
             category.names = c("State genes", "DeBoever et al."),
             fill = c("#f6c6ad", "#fba200"), col = "black", cex = 1.4, alpha = c(0.8,0.8),
             main = "Previously identified flexible genes", imagetype = "png", 
             main.fontfamily = "Arial", main.fontface = "bold", fontfamily = "Arial",
             cat.fontfamily = "Arial", cat.dist = c(-0.056,-0.07), cat.cex = c(1.2,1.2),
             disable.logging = T, force.unique = T,
             filename = paste(output_dir, "in_silico_validation/flexible_deboever_validation.png", sep = ""))
venn.diagram(x = list(state_enrich2, intersect(deboever, symbols)), 
             category.names = c("", ""),
             fill = c("#f6c6ad", "#fba200"), col = "black", cex = 1.4, alpha = c(0.8,0.8),
             main = "Previously identified flexible genes", imagetype = "png", 
             main.fontfamily = "Arial", main.fontface = "bold", fontfamily = "Arial",
             cat.fontfamily = "Arial", cat.dist = c(-0.056,-0.07), cat.cex = c(1.2,1.2),
             disable.logging = T, force.unique = T,
             filename = paste(output_dir, "in_silico_validation/flexible_deboever_validation_nolabels.png", sep = ""))

## Stable-polymorphic and GTEx eQTL enrichment ----
gtex <- read_tsv("~/Desktop/work_repo/data/GTEx_Analysis_v10_eQTL_updated/Whole_Blood.v10.eGenes.txt.gz")
gtex <- gtex %>% 
  dplyr::filter(qval < 0.01)
length(gtex$gene_name)
# [1] 11494
length(intersect(gtex$gene_name, symbols))
# [1] 3781
## lintersect
lintersect <- function(x,y) {
  length(intersect(unique(x),unique(y)))
}
phyper(q = lintersect(trait_enrich2, gtex$gene_name) - 1, #983 - 1
       m = lintersect(gtex$gene_name, symbols), #3781
       n = length(setdiff(symbols, gtex$gene_name)), #2318
       k = length(trait_enrich2), #1356
       lower.tail = F)
# [1] 2.8339e-20 (840 expected)
venn.diagram(x = list(trait_enrich2, intersect(gtex$gene_name, symbols)), 
             category.names = c("Trait genes", "GTEx eGenes (q < 0.01)"),
             fill = c("#84cbeb", "#005fbf"), col = "black", cex = 1.4, alpha = c(0.8,0.8),
             main = "Previously identified GTEx eQTLs", imagetype = "png", 
             main.fontfamily = "Arial", main.fontface = "bold", fontfamily = "Arial",
             cat.fontfamily = "Arial", cat.dist = c(0.035,-0.08), cat.pos = c(200,45), cat.cex = c(1.3, 1.2),
             disable.logging = T, force.unique = T,
             filename = paste(output_dir, "in_silico_validation/characteristic_GTEx_validation.png", sep = ""))
venn.diagram(x = list(trait_enrich2, intersect(gtex$gene_name, symbols)), 
             category.names = c("", ""),
             fill = c("#84cbeb", "#005fbf"), col = "black", cex = 1.4, alpha = c(0.8,0.8),
             main = "Previously identified GTEx eQTLs", imagetype = "png", 
             main.fontfamily = "Arial", main.fontface = "bold", fontfamily = "Arial",
             cat.fontfamily = "Arial", cat.dist = c(0.035,-0.08), cat.pos = c(200,45), cat.cex = c(1.3, 1.2),
             disable.logging = T, force.unique = T,
             filename = paste(output_dir, "in_silico_validation/characteristic_GTEx_validation_nolabels.png", sep = ""))


## Checking for under-enrichment ----
## Lower tail is true
## You can see expected values at https://systems.crump.ucla.edu/hypergeometric/index.php 
## Stable-polymorphic with DeBoever
phyper(q = lintersect(trait_enrich2, deboever), #107
       m = lintersect(deboever, symbols), #681
       n = length(setdiff(symbols, deboever)), #5418
       k = length(trait_enrich2), #1356
       lower.tail = T)
# [1] 4.541019e-06, expected 151 overlaps

## Dynamic with GTEx
phyper(q = lintersect(state_enrich2, gtex$gene_name), #741
       m = lintersect(gtex$gene_name, symbols), #3781
       n = length(setdiff(symbols, gtex$gene_name)), #2318
       k = length(state_enrich2), #1327
       lower.tail = T)
# [1] 1.251241e-07, 822 expected

# Functional enrichment analysis ----
## Get entrez IDs for enrichKEGG
categories <- list(trait_enrich2, state_enrich2)
names(categories) <- c("trait", "state")
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
# [1] 6151

resultsGO <- list()
resultsKEGG <- list()

dir.create(paste(output_dir, "functional_enrichment", sep = ""))

## Perform enrichment loop ----
for (i in 1:2) {
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
  rkf <- as.data.frame(setReadable(resKEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID"))
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
resultsKEGG[[1]]$passed <- 0
resultsKEGG[[1]]$passed <- resultsKEGG[[1]]$p.adjust < 0.05
resultsKEGG[[2]]$passed <- 0
resultsKEGG[[2]]$passed <- resultsKEGG[[2]]$p.adjust < 0.05

## Fixing some of the descriptions with longer names

## Flexible genes
grep("transferring groups", resultsGO[[2]]$Description)
resultsGO[[2]]$Description[29]
# "acyltransferase activity, transferring groups other than amino-acyl groups"
resultsGO[[2]]$Description[29] <- "acyltransferase activity*"

## Generate plots
colors <- list(sp = c("#84cbeb"), fl = c("#f6c6ad"))
for (i in 1:2) {
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
         title =  "GO Biological Process") +
    ylim(lims) +
    theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16))
  ## GO_CC
  go_cc[[i]] <- ggplot(na.omit(head(resultsGO[[i]][which(resultsGO[[i]]$ONTOLOGY == "CC"),],10)), aes(x = factor(Description, levels = rev(Description)), y = -log10(pvalue), fill = passed)) +
    geom_col(color = "black", width = 0.5, show.legend = F) +
    scale_fill_manual(values = colors[[i]]) +
    theme_cowplot() +
    coord_flip() +
    labs(y = "", x = "",
         title =  "GO Cellular Compartment") +
    ylim(lims) +
    theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16))
  ## GO_MF
  go_mf[[i]] <- ggplot(na.omit(head(resultsGO[[i]][which(resultsGO[[i]]$ONTOLOGY == "MF"),],10)), aes(x = factor(Description, levels = rev(Description)), y = -log10(pvalue), fill = passed)) +
    geom_col(color = "black", width = 0.5, show.legend = F) +
    scale_fill_manual(values = colors[[i]]) +
    theme_cowplot() +
    coord_flip() +
    labs(y = "", x = "",
         title = "GO Molecular Function") +
    ylim(lims) +
    theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16))
  ## KEGG
  kegg[[i]] <- ggplot(na.omit(head(resultsKEGG[[i]], 10)), aes(x = factor(Description, levels = rev(Description)), y = -log10(pvalue), fill = passed)) +
    geom_col(color = "black", width = 0.5, show.legend = F) +
    scale_fill_manual(values = colors[[i]]) +
    theme_cowplot() +
    coord_flip() +
    labs(y = "", x = "",
         title =  "KEGG Pathways") +
    ylim(lims) +
    theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16))
  
}

## Stable-polymorphic function
a <- plot_grid(go_bp[[1]], go_cc[[1]], go_mf[[1]], kegg[[1]], ncol = 2, align = 'v')
## FLexible function
b <- plot_grid(go_bp[[2]], go_cc[[2]], go_mf[[2]], kegg[[2]], ncol = 2, align = 'v')

pdf(file = paste(output_dir, "functional_enrichment/figure3.pdf", sep = ""), height = 10, width = 20)
a
dev.off()
pdf(file = paste(output_dir, "functional_enrichment/figure4.pdf", sep = ""), height = 10, width = 20)
b 
dev.off()

# Supplementary data file 2 ----
sp <- read.table(paste(output_dir, "trait_gene_list.txt", sep = ""))
fl <- read.table(paste(output_dir, "state_gene_list.txt", sep = ""))
sdf2 <- list(sp, fl)
names(sdf2) <- c("trait", "state")
dataset <- useDataset(dataset = 'hsapiens_gene_ensembl', mart = useMart('ENSEMBL_MART_ENSEMBL'))
gene_names <- c()
sdf2 <- lapply(sdf2, function(x) {
  gene_names <- getBM(attributes = c('hgnc_symbol','description'),
                      filters = 'hgnc_symbol',
                      values = x[[1]], 
                      mart = dataset)
  x <- gene_names
})
lapply(sdf2, head)
## First list
# $trait
# hgnc_symbol                                                                                  description
# 1       ABCA1                   ATP binding cassette subfamily A member 1 [Source:HGNC Symbol;Acc:HGNC:29]
# 2       ABCB1                   ATP binding cassette subfamily B member 1 [Source:HGNC Symbol;Acc:HGNC:40]
# 3       ABCC5                   ATP binding cassette subfamily C member 5 [Source:HGNC Symbol;Acc:HGNC:56]
# 4       ABCG1                   ATP binding cassette subfamily G member 1 [Source:HGNC Symbol;Acc:HGNC:73]
# 5      ABLIM1                                 actin binding LIM protein 1 [Source:HGNC Symbol;Acc:HGNC:78]
# 6       ACAP2 ArfGAP with coiled-coil, ankyrin repeat and PH domains 2 [Source:HGNC Symbol;Acc:HGNC:16469]
writexl::write_xlsx(x = sdf2, path = paste(output_dir, "supplementary_data_2.xlsx", sep = ""))

# Supplementary data file 3 ----
cgo <- read.csv(paste(output_dir, "/functional_enrichment/trait_GO.csv", sep = ""))
ckg <- read.csv(paste(output_dir, "/functional_enrichment/trait_KEGG.csv", sep = ""))
fgo <- read.csv(paste(output_dir, "/functional_enrichment/state_GO.csv", sep = ""))
fkg <- read.csv(paste(output_dir, "/functional_enrichment/state_KEGG.csv", sep = ""))
sdf3 <- list(cgo, ckg, fgo, fkg)
names(sdf3) <- c("Trait GO", "Trait KEGG", "State GO", "State KEGG")
lapply(sdf3, head)
writexl::write_xlsx(x = sdf3, path = paste(output_dir, "supplementary_data_3.xlsx", sep = ""))