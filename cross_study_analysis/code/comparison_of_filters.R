## Which is the best for finding stable-polymorphic genes: , .vp, or .rp?
## October 27 2024
library(tidyverse)
library(VennDiagram)

## Input: the code and filters from stable-polymorphic gene analysis
input <- "~/Desktop/work_repo/github/cross_study_analysis/"
## Output: venn diagram of the filter stable-polymorphic filter overlap
output <- "~/Desktop/work_repo/github/cross_study_analysis/output/"
sourcePartial <- function(fn, startTag='#from here', endTag='#to here') {
  lines <- scan(fn, what = character(), sep = "\n", quiet = TRUE)
  st <- grep(startTag, lines, fixed = T)
  en <- grep(endTag, lines, fixed = T)
  tc <- textConnection(lines[(st+1):(en-1)])
  source(tc)
  close(tc)
}
sourcePartial(fn = paste(input, "code/trait_genes.R", sep = ""), 
              startTag = '## Load in', 
              endTag = '# Trait genes across studies ----'
              )

symbols <- read.table(paste(input, "output/common_symbols6099.txt", sep = ""))
symbols <- symbols$x
sd_count <- data.frame(matrix(ncol = 13, nrow = 6099))
rownames(sd_count) <- symbols
colnames(sd_count) <- names(sd)
for (i in 1:nrow(sd_count)) {
  for (j in 1:ncol(sd_count)) {
    sd_count[i,j] <- rownames(sd_count)[i] %in% sd[[j]]
  }
}
head(sd_count)
sd_count <- sd_count %>% 
  mutate(Microarray = obermoser1 + obermoser2 + obermoser3 + obermoser4 + dusek + rusch + meaburn1 +
           meaburn2 + karlovich1 + karlovich2,
         RNAseq = gomez + gosch + larocca) %>% 
  mutate(Studies = gomez  + gosch  + (obermoser1 + obermoser2 + obermoser3 + obermoser4  > 0) + dusek  +
           rusch  + larocca  + (meaburn1 + meaburn2  > 0) + (karlovich1 + karlovich2 > 0))
head(sd_count)
sd_symbols <- sd_count %>% 
  dplyr::filter(Microarray > 0, RNAseq > 0, Studies >= 4) %>% 
  rownames()
sd_symbols
# [1] "ERAP2"    "HLA-DQB1" "HLA-DRB5"

vp_count <- data.frame(matrix(ncol = 13, nrow = 6099))
rownames(vp_count) <- symbols
colnames(vp_count) <- names(vp)
for (i in 1:nrow(vp_count)) {
  for (j in 1:ncol(vp_count)) {
    vp_count[i,j] <- rownames(vp_count)[i] %in% vp[[j]]
  }
}
head(vp_count)
vp_count <- vp_count %>% 
  mutate(Microarray = obermoser1  + obermoser2  + obermoser3  + obermoser4  + dusek  + rusch  + meaburn1  + meaburn2 +
           karlovich1 + karlovich2,
         RNAseq = gomez  +gosch  +larocca  ) %>% 
  mutate(Studies = gomez   + gosch   + (obermoser1  +obermoser2  +obermoser3  +obermoser4   > 0) + dusek   +
           rusch   + larocca   + (meaburn1  +meaburn2   > 0) + (karlovich1 + karlovich2 > 0))
head(vp_count)
vp_symbols <- vp_count %>% 
  dplyr::filter(Microarray > 0, RNAseq > 0, Studies >= 4) %>% 
  rownames()
length(vp_symbols)
# [1] 1260

rp_count <- data.frame(matrix(ncol = 13, nrow = 6099))
rownames(rp_count) <- symbols
colnames(rp_count) <- names(rp)
for (i in 1:nrow(rp_count)) {
  for (j in 1:ncol(rp_count)) {
    rp_count[i,j] <- rownames(rp_count)[i] %in% rp[[j]]
  }
}
head(rp_count)
rp_count <- rp_count %>% 
  mutate(Microarray = obermoser1  + obermoser2  + obermoser3  + obermoser4  + dusek  + rusch  + meaburn1  + meaburn2 +
           karlovich1 + karlovich2,
         RNAseq = gomez  +gosch  +larocca  ) %>% 
  mutate(Studies = gomez   + gosch   + (obermoser1  +obermoser2  +obermoser3  +obermoser4   > 0) + dusek   +
           rusch   + larocca   + (meaburn1  +meaburn2   > 0) + (karlovich1 + karlovich2 > 0))
head(rp_count)
rp_symbols <- rp_count %>% 
  dplyr::filter(Microarray > 0, RNAseq > 0, Studies >= 4) %>% 
  rownames()
length(rp_symbols)
# [1] 869

venn.diagram(x = list(sd_symbols, vp_symbols, rp_symbols), 
             category.names = c("Standard deviation", "variancePartition", "Repeatability"),
             fill = c("orange1", "royalblue3", "darkviolet"), col = "black",
             main = "Trait gene filters (>= 4 studies)", imagetype = "png", 
             main.fontfamily = "Arial", main.fontface = "bold", fontfamily = "Arial", cat.cex = c(1.078,1.13,1.23),
             cat.fontfamily = "Arial", cat.dist = c(-0.078,-0.09, -0.09),
             disable.logging = T,
             filename = paste(output, "comparison_of_filters.png", sep = ""))
## WITHOUT LABELS
venn.diagram(x = list(sd_symbols, vp_symbols, rp_symbols), 
             category.names = c("", "", ""),
             fill = c("orange1", "royalblue3", "darkviolet"), col = "black",
             main = "", imagetype = "png", 
             main.fontfamily = "Arial", main.fontface = "bold", fontfamily = "Arial", cat.cex = c(1.078,1.13,1.23),
             cat.fontfamily = "Arial", cat.dist = c(-0.078,-0.09, -0.09), cex = 1.8,
             disable.logging = T,
             filename = paste(output, "supplementary_figure_s1.png", sep = ""))
