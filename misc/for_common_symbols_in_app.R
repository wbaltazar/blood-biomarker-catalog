## June 28 2024
## For the application, we decided to use only the gene symbols which we are confident were captured across
## all seven datasets. For us, those are the 9474 symbols that I found in my eQTl enrichment analysis. So,
## we will use that to filter all the rds files and annotation files to reduce space and speed up computation.
## formatting_expr_for_app.R should be run prior to this code
## biomaRt_all_transcripts.R should be run prior to this code

## input: RDS files generated from study_data files + common symbols file
symbols <- read.table("~/Desktop/work_repo/github/cross_study_analysis/output/common_symbols6099.txt")
symbols <- symbols$x
length(symbols)
# [1] 6099
setwd("~/Desktop/work_repo/Box organization/1results/RShiny-application/")
e1 <- read_rds("./data/expr/gomez.rds")
e2 <- read_rds("./data/expr/Meaburn1.rds")
e3 <- read_rds("./data/expr/Meaburn2.rds")
e4 <- read_rds("./data/expr/gosch.rds")
e5 <- read_rds("./data/expr/obermoser1.rds")
e6 <- read_rds("./data/expr/obermoser2.rds")
e7 <- read_rds("./data/expr/obermoser3.rds")
e8 <- read_rds("./data/expr/obermoser4.rds")
e9 <- read_rds("./data/expr/dusek.rds")
e10 <- read_rds("./data/expr/rusch.rds")
e11 <- read_rds("./data/expr/larocca.rds")
e12 <- read_rds("./data/expr/Karlovich1.rds")
e13 <- read_rds("./data/expr/Karlovich2.rds")

# Gomez ----
e1 <- e1[which(rownames(e1) %in% symbols),]
dim(e1)
# [1] 6099   28
library(EnsDb.Hsapiens.v86)
names <- mapIds(EnsDb.Hsapiens.v86, keys = rownames(e1), column = "GENEID", keytype = "SYMBOL")
any(is.na(names))
# [1] FALSE
attr(rownames(e1), "names") <- names
head(rownames(e1))
# ENSG00000157933 ENSG00000142655 ENSG00000162441 ENSG00000177000 ENSG00000069424 ENSG00000074800 
# "SKI"         "PEX14"          "LZIC"         "MTHFR"        "KCNAB2"          "ENO1" 
write_rds(e1, "./data/expr/gomez.rds")

# Meaburn ----
find <- str_replace_all(rownames(e2), "\\.(?!\\d$)", "-") %>% str_remove_all(pattern = "\\..")
length(unique(find)) # [1] 14948
index <- which(find %in% symbols)
length(index) # [1] 13658
## This makes sure we include all probes that represent our common set of genes.
e2 <- e2[index,]
dim(e2)
# [1] 13658    20
rownames(e2) <-  str_replace_all(rownames(e2), "\\.(?!\\d$)", "-")
write_rds(e2, "./data/expr/Meaburn1.rds")
find <- str_replace_all(rownames(e3), "\\.(?!\\d$)", "-") %>% str_remove_all(pattern = "\\..")
length(unique(find)) # [1] 15598
index <- which(find %in% symbols)
length(index) # [1] 13792
e3 <- e3[index,]
dim(e3)
# [1] 13792    18
rownames(e3) <- str_replace_all(rownames(e3), "\\.(?!\\d$)", "-")
write_rds(e3, "./data/expr/Meaburn2.rds")

# Gosch ----
dim(e4[str_replace_all(rownames(e4), "\\.(?!\\d$)", "-") %in% symbols,])
# [1] 6099   79
e4 <- e4[str_replace_all(rownames(e4), "\\.(?!\\d$)", "-") %in% symbols,]
rownames(e4) <- str_replace_all(rownames(e4), "\\.(?!\\d$)", "-")
names <- mapIds(EnsDb.Hsapiens.v86, keys = rownames(e4), column = "GENEID", keytype = "SYMBOL")
any(is.na(names))
# [1] FALSE
attr(rownames(e4), "names") <- names
head(rownames(e4))
# ENSG00000187608 ENSG00000078808 ENSG00000160087 ENSG00000175756 ENSG00000221978 ENSG00000242485 
# "ISG15"          "SDF4"        "UBE2J2"      "AURKAIP1"         "CCNL2"        "MRPL20" 
write_rds(e4, "./data/expr/gosch.rds")

# Obermoser ----
## Similar to what we did for Meaburn
find <- str_replace_all(rownames(e5), "\\.(?!\\d$)", "-") %>% str_remove_all(pattern = "\\..")
length(unique(find)) # [1] 16492
index <- which(find %in% symbols)
length(index)
# [1] 9566
e5 <- e5[index,]
dim(e5)
# [1] 9566  162
rownames(e5) <- str_replace_all(rownames(e5), "\\.(?!\\d$)", "-")
sum(rownames(e5) %in% symbols)
# [1] 6099
head(rownames(e5))
head(e5$Probe)
# [1] "ILMN_1812666" "ILMN_1800197" "ILMN_1765409" "ILMN_1798083" "ILMN_1732300" "ILMN_1755909"
write_rds(e5, "./data/expr/obermoser1.rds")

find <- str_replace_all(rownames(e6), "\\.(?!\\d$)", "-") %>% str_remove_all(pattern = "\\..")
length(unique(find)) # [1] 12562
index <- which(find %in% symbols)
length(index) 
# [1] 8711
e6 <- e6[index,]
dim(e6)
# [1] 8711   50
rownames(e6) <- str_replace_all(rownames(e6), "\\.(?!\\d$)", "-")
sum(rownames(e6) %in% symbols)
# [1] 6099
head(rownames(e6))
head(e6$Probe)
write_rds(e6, "./data/expr/obermoser2.rds")

find <- str_replace_all(rownames(e7), "\\.(?!\\d$)", "-") %>% str_remove_all(pattern = "\\..")
length(unique(find)) # [1] 15945
index <- which(find %in% symbols)
length(index) 
# [1] 9601
e7 <- e7[index,]
dim(e7)
# [1] 9601  144
rownames(e7) <- str_replace_all(rownames(e7), "\\.(?!\\d$)", "-")
sum(rownames(e7) %in% symbols)
# [1] 6099
head(rownames(e7))
head(e7$Probe)
write_rds(e7, "./data/expr/obermoser3.rds")

find <- str_replace_all(rownames(e8), "\\.(?!\\d$)", "-") %>% str_remove_all(pattern = "\\..")
length(unique(find)) # [1] 15299
index <- which(find %in% symbols)
length(index) 
# [1] 9297
e8 <- e8[index,]
dim(e8)
# [1] 9297   186
rownames(e8) <- str_replace_all(rownames(e8), "\\.(?!\\d$)", "-")
sum(rownames(e8) %in% symbols)
# [1] 60
head(rownames(e8))
head(e8$Probe)
write_rds(e8, "./data/expr/obermoser4.rds")

# Dusek ----
find <- str_replace_all(rownames(e9), "\\.(?!\\d$)", "-") %>% str_remove_all(pattern = "\\..")
length(unique(find)) # [1] 10912
index <- which(find %in% symbols)
length(index)
# [1] 10243
e9 <- e9[index,]
dim(e9)
# [1] 10243    45
rownames(e9) <- str_replace_all(rownames(e9), "\\.(?!\\d$)", "-")
sum(rownames(e9) %in% symbols)
# [1] 6099
write_rds(e9, "./data/expr/dusek.rds")

# Rusch ----
find <- str_replace_all(rownames(e10), "\\.(?!\\d$)", "-") %>% str_remove_all(pattern = "\\..")
length(unique(find)) # [1] 14100
index <- which(find %in% symbols)
length(index)
# [1] 12839
e10 <- e10[index,]
dim(e10)
# [1] 12839    38
rownames(e10) <- str_replace_all(rownames(e10), "\\.(?!\\d$)", "-")
sum(rownames(e10) %in% symbols)
# [1] 6099
write_rds(e10, "./data/expr/rusch.rds")

# LaRocca ----
dim(e11[str_replace_all(rownames(e11), "\\.(?!\\d$)", "-") %in% symbols,])
# [1] 6099   59
e11 <- e11[str_replace_all(rownames(e11), "\\.(?!\\d$)", "-") %in% symbols,]
rownames(e11) <- str_replace_all(rownames(e11), "\\.(?!\\d$)", "-")
names <- mapIds(EnsDb.Hsapiens.v86, keys = rownames(e11), column = "GENEID", keytype = "SYMBOL")
any(is.na(names))
# [1] FALSE
attr(rownames(e11), "names") <- names
head(rownames(e11))
write_rds(e11, "./data/expr/larocca.rds")

# Karlovich ----
find <- str_replace_all(rownames(e12), "\\.(?!\\d$)", "-") %>% str_remove_all(pattern = "\\..")
length(unique(find)) # [1] 11635
index <- which(find %in% symbols)
length(index) # [1] 11001
## This makes sure we include all probes that represent our common set of genes.
e12 <- e12[index,]
dim(e12)
# [1] 11001   67
rownames(e12) <-  str_replace_all(rownames(e12), "\\.(?!\\d$)", "-")
sum(rownames(e12) %in% symbols)
# [1] 6099
write_rds(e12, "./data/expr/Karlovich1.rds")
find <- str_replace_all(rownames(e13), "\\.(?!\\d$)", "-") %>% str_remove_all(pattern = "\\..")
length(unique(find)) # [1] 10954
index <- which(find %in% symbols)
length(index) # [1] 10450
e13 <- e13[index,]
dim(e13)
# [1] 10544    41
rownames(e13) <- str_replace_all(rownames(e13), "\\.(?!\\d$)", "-")
sum(rownames(e13) %in% symbols)
# 6099
write_rds(e13, "./data/expr/Karlovich2.rds")

## Testing to see if app.R still works. Everything works!
## Shrinking file sizes of annotations so that they load faster.
affy <- read.csv("./data/annotate/affy.csv")
affy <- affy[affy$affy_hg_u133_plus_2 %in% e2$Probe,]
dim(affy)
# [1] 13024    12
## Backed up files to /1results/annotation backup
write.csv(affy, "./data/annotate/affy.csv")

ilmn <- read.csv("./data/annotate/ilmn.csv")
probes <- union(union(e5$Probe, e6$Probe), union(e7$Probe, e8$Probe))
length(probes)
ilmn <- ilmn[ilmn$illumina_humanwg_6_v3 %in% probes,]
dim(ilmn)
# [1] 9383    12
write.csv(ilmn, "./data/annotate/ilmn.csv")

gomez <- read.csv("./data/annotate/gomez.csv")
all(attr(rownames(e1), "names") %in% gomez$ensembl_gene_id)
# [1] FALSE
notfound <- rownames(e1)[!(attr(rownames(e1), "names") %in% gomez$ensembl_gene_id)]
gomez <- gomez[gomez$ensembl_gene_id %in% attr(rownames(e1), "names"),]
dim(gomez)
# [1] 6098   10
library(biomaRt)
get <- c("ensembl_gene_id", "description", "external_gene_name", "gene_biotype", "ensembl_transcript_id",
         "transcript_appris", "external_transcript_name", "transcript_tsl", "transcript_biotype")
mart <- useDataset(dataset = "hsapiens_gene_ensembl", mart = useMart("ENSEMBL_MART_ENSEMBL"))
gomez_extra <- getBM(attributes = get, values = names(notfound), mart = mart, filters = "ensembl_gene_id")
gomez_extra <- gomez_extra %>% 
  group_by(ensembl_gene_id) %>% 
  summarise(
    description = paste(description, collapse = ", "),
    external_gene_name = paste(external_gene_name, collapse = ", "),
    gene_biotype = paste(gene_biotype, collapse = ", "),
    ensembl_transcript_id = paste(ensembl_transcript_id, collapse = ", "),
    transcript_appris = paste(transcript_appris, collapse = ", "),
    external_transcript_name = paste(external_transcript_name, collapse = ", "),
    transcript_tsl = paste(transcript_tsl, collapse = ", "),
    transcript_biotype = paste(transcript_biotype, collapse = ", ")
  )
dim(gomez_extra)
# [1] 0  9, this ENS id was updated for GRCh38.p14
gomez_extra <- getBM(attributes = get, values = "ENSG00000291237", mart = mart, filters = "ensembl_gene_id")
gomez_extra <- gomez_extra %>% 
  group_by(ensembl_gene_id) %>% 
  summarise(
    description = paste(description, collapse = ", "),
    external_gene_name = paste(external_gene_name, collapse = ", "),
    gene_biotype = paste(gene_biotype, collapse = ", "),
    ensembl_transcript_id = paste(ensembl_transcript_id, collapse = ", "),
    transcript_appris = paste(transcript_appris, collapse = ", "),
    external_transcript_name = paste(external_transcript_name, collapse = ", "),
    transcript_tsl = paste(transcript_tsl, collapse = ", "),
    transcript_biotype = paste(transcript_biotype, collapse = ", ")
  )
dim(gomez_extra)
# [1] 1 9
gomez <- rbind(gomez[,-1], gomez_extra)
write.csv(gomez, "./data/annotate/gomez.csv")

gosch <- read.csv("./data/annotate/gosch.csv")
all(attr(rownames(e4), "names") %in% gosch$ensembl_gene_id)
# [1] FALSE
notfound <- rownames(e4)[!(attr(rownames(e4), "names") %in% gosch$ensembl_gene_id)]
print(notfound)
# ENSG00000112096 
# "SOD2" same as gomez
gosch <- gosch[gosch$ensembl_gene_id %in% attr(rownames(e4), "names"),]
gosch <- rbind(gosch[,-1], gomez_extra)
write.csv(gosch, "./data/annotate/gosch.csv")

larocca <- read.csv("./data/annotate/larocca.csv")
all(attr(rownames(e11), "names") %in% larocca$ensembl_gene_id)
# [1] FALSE
notfound <- rownames(e11)[!(attr(rownames(e11), "names") %in% larocca$ensembl_gene_id)]
length(notfound)
# [1] 1
print(notfound)
# ENSG00000112096 
# "SOD2" same as gomez
larocca <- rbind(larocca[,-1], gomez_extra)
write.csv(larocca, "./data/annotate/larocca.csv")