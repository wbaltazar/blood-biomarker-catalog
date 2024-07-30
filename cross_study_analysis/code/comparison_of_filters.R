## Which is the best for finding stable-polymorphic genes: .sd, .vp, or .rp?
## July 30 2024
library(tidyverse)
library(VennDiagram)

## Input: the code and filters from stable-polymorphic gene analysis
input <- "~/Desktop/work_repo/github/cross_study_analysis/"
## Output: venn diagram of the filter stable-polymorphic filter overlap
output <- "~/Desktop/work_repo/github/cross_study_analysis/output/"
source(paste(input, "code/stable_polymorphic_genes.R", sep = ""))
sd
rp
vp
symbols <- read.table(paste(input, "output/common_symbols9474.txt", sep = ""))
symbols <- symbols$x
sd_count <- data.frame(matrix(ncol = 11, nrow = 9474))
rownames(sd_count) <- symbols
colnames(sd_count) <- names(sd)
for (i in 1:nrow(sd_count)) {
  for (j in 1:ncol(sd_count)) {
    sd_count[i,j] <- rownames(sd_count)[i] %in% sd[[j]]
  }
}
head(sd_count)
sd_count <- sd_count %>% 
  mutate(Microarray = obermoser1.sd+obermoser2.sd+obermoser3.sd+obermoser4.sd+dusek.sd+rusch.sd+meaburn1.sd+meaburn2.sd,
         RNAseq = gomez.sd+gosch.sd+larocca.sd) %>% 
  mutate(Studies = gomez.sd + gosch.sd + (obermoser1.sd+obermoser2.sd+obermoser3.sd+obermoser4.sd > 0) + dusek.sd +
           rusch.sd + larocca.sd + (meaburn1.sd+meaburn2.sd > 1))
head(sd_count)
sd_symbols <- sd_count %>% 
  dplyr::filter(Microarray > 0, RNAseq > 0, Studies >= 4) %>% 
  rownames()
sd_symbols
# [1] "ERAP2"    "HLA-DQB1"

vp_count <- data.frame(matrix(ncol = 11, nrow = 9474))
rownames(vp_count) <- symbols
colnames(vp_count) <- names(vp)
for (i in 1:nrow(vp_count)) {
  for (j in 1:ncol(vp_count)) {
    vp_count[i,j] <- rownames(vp_count)[i] %in% vp[[j]]
  }
}
head(vp_count)
vp_count <- vp_count %>% 
  mutate(Microarray = obermoser1.vp+obermoser2.vp+obermoser3.vp+obermoser4.vp+dusek.vp+rusch.vp+meaburn1.vp+meaburn2.vp,
         RNAseq = gomez.vp+gosch.vp+larocca.vp) %>% 
  mutate(Studies = gomez.vp + gosch.vp + (obermoser1.vp+obermoser2.vp+obermoser3.vp+obermoser4.vp > 0) + dusek.vp +
           rusch.vp + larocca.vp + (meaburn1.vp+meaburn2.vp > 1))
head(vp_count)
vp_symbols <- vp_count %>% 
  dplyr::filter(Microarray > 0, RNAseq > 0, Studies >= 4) %>% 
  rownames()
length(vp_symbols)
# [1] 1357

rp_count <- data.frame(matrix(ncol = 11, nrow = 9474))
rownames(rp_count) <- symbols
colnames(rp_count) <- names(rp)
for (i in 1:nrow(rp_count)) {
  for (j in 1:ncol(rp_count)) {
    rp_count[i,j] <- rownames(rp_count)[i] %in% rp[[j]]
  }
}
head(rp_count)
rp_count <- rp_count %>% 
  mutate(Microarray = obermoser1.rp+obermoser2.rp+obermoser3.rp+obermoser4.rp+dusek.rp+rusch.rp+meaburn1.rp+meaburn2.rp,
         RNAseq = gomez.rp+gosch.rp+larocca.rp) %>% 
  mutate(Studies = gomez.rp + gosch.rp + (obermoser1.rp+obermoser2.rp+obermoser3.rp+obermoser4.rp > 0) + dusek.rp +
           rusch.rp + larocca.rp + (meaburn1.rp+meaburn2.rp > 1))
head(rp_count)
rp_symbols <- rp_count %>% 
  dplyr::filter(Microarray > 0, RNAseq > 0, Studies >= 4) %>% 
  rownames()
length(rp_symbols)
# [1] 1320

venn.diagram(x = list(sd_symbols, vp_symbols, rp_symbols), 
             category.names = c(".sd", ".vp", ".rp"),
             fill = c("lightblue", "royalblue3", "grey75"), col = "black",
             main = "Stable-polymorphic filters (>= 4 studies)", imagetype = "png", 
             main.fontfamily = "Arial", main.fontface = "bold", fontfamily = "Arial",
             cat.fontfamily = "Arial", cat.dist = c(-0.05,-0.05, -0.05),
             disable.logging = T,
             filename = paste(output, "comparison_of_filters.png", sep = ""))