## INPUT: the F-statistic from limma comparing single-day measurements over time
## OUTPUT: the intersection of differentially expressed genes in one day, "diurnal genes"
setwd("~/Desktop/work_repo/github/")
common_symbols <- read.table("./cross_study_analysis/output/common_symbols6099.txt")
common_symbols <- common_symbols$x

gosch <- read.csv(file = "./Gosch_study_data/output/diurnal_genes_limma_F.csv")
gosch %>% filter(X %in% common_symbols, P.Value < 0.05) %>% nrow()
# [1] 4017
diurnal_genes_gosch_F <- gosch %>% filter(X %in% common_symbols, P.Value < 0.05) %>% dplyr::select(X)
obermoser <- read.csv(file = "./Obermoser_study_data/output/diurnal_genes/diurnal_limma_F.csv")
obermoser %>% filter(Symbol %in% common_symbols, P.Value < 0.05) %>% unique() %>% nrow()
# [1] 1157
diurnal_genes_obermoser_F <- obermoser %>% filter(Symbol %in% common_symbols, P.Value < 0.05) %>% dplyr::select(Symbol)

library(VennDiagram)
venn.diagram(x = list(diurnal_genes_gosch_F$X, diurnal_genes_obermoser_F$Symbol),
             filename = "./misc/diurnal_venn_diagram.png", imagetype = "png", disable.logging = T,
             category.names = c("Gosch", "Obermoser"), cat.cex = 1.25, cat.fontfamily = "Arial",
             fontfamily = "Arial", cex = 1.5, fill = c("pink", "skyblue"), cat.pos = c(-30, 150))
diurnal_gene_intersect <- intersect(diurnal_genes_gosch_F$X, diurnal_genes_obermoser_F$Symbol)
gosch_p_values <- gosch[gosch$X %in% diurnal_gene_intersect, c("X", "P.Value")]
colnames(gosch_p_values)[1] <- "Symbol"
obermoser_p_values <- obermoser[obermoser$Symbol %in% diurnal_gene_intersect, c("Symbol", "P.Value")]
diurnal_gene_intersect <- data.frame(Symbol = diurnal_gene_intersect)
diurnal_gene_intersect <- merge(diurnal_gene_intersect, obermoser_p_values)
diurnal_gene_intersect <- merge(diurnal_gene_intersect, gosch_p_values, by = "Symbol")
colnames(diurnal_gene_intersect) <- c("Symbol", "Obermoser_F_pvalue", "Gosch_F_pvalue")
write.csv(diurnal_gene_intersect, file = "./misc/diurnal_genes.csv")