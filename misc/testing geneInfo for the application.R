use <- affy_data[affy_data$affy_hg_u133_plus_2 == "204141_at",]
data.frame(`Gene ID` = unlist(str_split(use$ensembl_gene_id, ",")))
print(use)
# ENSG00000137267
# "TUBB2A"
# 204141_at

## Genes not annotated
# ILMN_1837476
# X1552563_a_at
use <- anno_data$gomez[anno_data$gomez$ensembl_gene_id == attr(rownames(e1),"names")[which(rownames(e1) == "TUBB2A")],]
# works
use <- anno_data$affy[anno_data$affy$affy_hg_u133_plus_2 == "204141_at",]
# works
use <- anno_data$affy[anno_data$affy$affy_hg_u133_plus_2 == "1552563_a_at",]
# works
use <- anno_data$ilmn[anno_data$ilmn$illumina_humanwg_6_v3 == "ILMN_1837476",]
# doesn't work, no annotation
dim(use)
## so if rows = 0, return nothing
use <- anno_data$larocca[anno_data$larocca$ensembl_gene_id == attr(rownames(e11),"names")[which(rownames(e11) == "TUBB2A")],]
# works
head(rownames(e11)[which(is.na(attr(rownames(e11), "names")))])
# APOBEC3A_B
use <- anno_data$larocca[anno_data$larocca$ensembl_gene_id == attr(rownames(e11),"names")[which(rownames(e11) == "APOBEC3A_B")],]

## things to try
## Gosch - DDX11L
head(e4)
use <- anno_data$gosch[anno_data$gosch$ensembl_gene_id == attr(rownames(e4),"names")[which(rownames(e4) == "RP11.465B22.5")],]
## PROBLEM: gene name is WAS protein family homolog 8, pseudogene [Source:HGNC Symbol;Acc:HGNC:53913]
## also, the "APPRIS" column is empty after str_split is used on it. So, we add a ","
## Recognize if there is a , in the gene description
if (grepl(",", substr(use$description, start = 1, stop = str_locate(use$description, "\\[")))) {
  use$description <- gsub("\\],", "\\]/", use$description)
  skipdesc = T
}
## Replace empty vectors with empty space
test <- apply(use, 2, function(x) {trimws(unlist(strsplit(x,",")))})
test <- lapply(test, function(x) {
  if (length(x) == 0) {
    rep(", ", times = str_count(use$ensembl_gene_id, pattern = "E"))
  }
})
lapply(test, function(x) {
  if (!is.null(x)) {
    print(x)
  }
})

t(as.data.frame(apply(use, 2, function(x) {trimws(unlist(strsplit(x,",")))})))
