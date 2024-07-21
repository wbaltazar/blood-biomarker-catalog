## PrediXcan is a model that takes eQTL data and genotypes to try and predict a model transcriptome.
## MASHR comes from the same lab. I downloaded their data to explore which genes are best predicted
## by genotype (in their original paper, ERAP2 was the best predicted, and conincidentally the most
## stable gene from our data.)

library(tidyverse)
library(AnnotationDbi)
setwd("./Desktop/work_repo/")

mashr_wb <- read.table("./data/eqtl/mashr/mashr_Whole_Blood.txt.gz", header = T)
head(mashr_wb)
#                 GENE                  RSID1                  RSID2       VALUE
# 1 ENSG00000000457.13 chr1_169894240_G_T_b38 chr1_169894240_G_T_b38  0.17362599
# 2 ENSG00000000457.13 chr1_169894240_G_T_b38 chr1_169891332_G_A_b38  0.09154263
# 3 ENSG00000000457.13 chr1_169891332_G_A_b38 chr1_169891332_G_A_b38  0.50753525
# 4 ENSG00000000460.16 chr1_169661963_G_A_b38 chr1_169661963_G_A_b38  0.21071987
# 5 ENSG00000000460.16 chr1_169661963_G_A_b38 chr1_169697456_A_T_b38 -0.03676570
# 6 ENSG00000000460.16 chr1_169697456_A_T_b38 chr1_169697456_A_T_b38  0.39746015
db <- dbFileConnect("./data/eqtl/mashr/mashr_Whole_Blood.db")


## not sure what to do after this