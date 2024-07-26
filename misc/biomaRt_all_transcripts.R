## Date: Jul 24 2024

# This code uses the rds files generated from formatting_expr_for_app.R and for_common_symbols_in_app.R.
# This code pulls the gene information for the searchable symbols in the RShiny app.

# All files will be output to a folder in our Rshiny application in /data/annotate/
setwd("~/Desktop/work_repo/Box organization/1results/RShiny-application/")

## Repository of gene names
library(biomaRt)
mart <- useDataset(dataset = "hsapiens_gene_ensembl", mart = useMart("ENSEMBL_MART_ENSEMBL"))
ats <- listAttributes(mart = mart, page = "feature_page")
View(ats)
fil <- listFilters(mart = mart)
View(fil)
get <- c("ensembl_gene_id", "description", "external_gene_name", "gene_biotype", "ensembl_transcript_id",
         "transcript_appris", "external_transcript_name", "transcript_tsl", "transcript_biotype")
## Possible filters
c("ensembl_gene_id", "external_gene_name", "external_synonym", "affy_hg_u133_plus_2", "illumina_humanwg_6_v3")
## Test
getBM(attributes = get, values = "TUBB2A", 
      mart = mart, filters = "external_gene_name")

## Affymetrix ----
get_affy <- c(get, "affy_hg_u133_plus_2")
affy_anno <- getBM(attributes = get_affy, values = e2$Probe, 
                   mart = mart, filters = "affy_hg_u133_plus_2")
## Takes 5 minutes to pull
dim(affy_anno)
# [1] 224377     10
length(unique(affy_anno$affy_hg_u133_plus_2))
# [1] 43321
affy_data <- affy_anno %>% 
  group_by(affy_hg_u133_plus_2) %>% 
  summarise(
    ensembl_gene_id = paste(ensembl_gene_id, collapse = ", "),
    description = paste(description, collapse = ", "),
    external_gene_name = paste(external_gene_name, collapse = ", "),
    gene_biotype = paste(gene_biotype, collapse = ", "),
    ensembl_transcript_id = paste(ensembl_transcript_id, collapse = ", "),
    transcript_appris = paste(transcript_appris, collapse = ", "),
    external_transcript_name = paste(external_transcript_name, collapse = ", "),
    transcript_tsl = paste(transcript_tsl, collapse = ", "),
    transcript_biotype = paste(transcript_biotype, collapse = ", ")
  )
dim(affy_data)
# [1] 43321    10
head(affy_data)
write.csv(affy_data, "./data/annotate/affy.csv")


## Illumina ----
get_ilmn <- c(get, "illumina_humanwg_6_v3")
ilmn_ids <- union(union(e5$Probe, e6$Probe), union(e7$Probe, e8$Probe))
length(ilmn_ids)
# [1] 33049
ilmn_anno <- getBM(attributes = get_ilmn, values = ilmn_ids, 
                   mart = mart, filters = "illumina_humanwg_6_v3")
## Takes 3 minutes to pull
dim(ilmn_anno)
# [1] 134674     10
length(unique(ilmn_anno$illumina_humanwg_6_v3))
# [1] 25938
ilmn_data <- ilmn_anno %>% 
  group_by(illumina_humanwg_6_v3) %>% 
  summarise(
    ensembl_gene_id = paste(ensembl_gene_id, collapse = ", "),
    description = paste(description, collapse = ", "),
    external_gene_name = paste(external_gene_name, collapse = ", "),
    gene_biotype = paste(gene_biotype, collapse = ", "),
    ensembl_transcript_id = paste(ensembl_transcript_id, collapse = ", "),
    transcript_appris = paste(transcript_appris, collapse = ", "),
    external_transcript_name = paste(external_transcript_name, collapse = ", "),
    transcript_tsl = paste(transcript_tsl, collapse = ", "),
    transcript_biotype = paste(transcript_biotype, collapse = ", ")
  )
dim(ilmn_data)
# [1] 25938    10
write.csv(ilmn_data, "./data/annotate/ilmn.csv")


## Gomez ----
### By ensembl IDs ----
gomez_ens <- names(rownames(e1)[!is.na(names(rownames(e1)))])
length(gomez_ens)
# [1] 15416. Recall all entries in Gomez were collected using Ensembl IDs.
gomez_ens_anno <- getBM(attributes = get, values = gomez_ens, 
                        mart = mart, filters = "ensembl_gene_id")
## Takes 30 seconds
dim(gomez_ens_anno)
# [1] 148885      9
length(unique(gomez_ens_anno$ensembl_gene_id))
# [1] 15372
gomez_data <- gomez_ens_anno %>% 
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
dim(gomez_data)
# [1] 15372     9
write.csv(gomez_data, file = "./data/annotate/gomez.csv")

## Gosch ----
### By ensembl IDs ----
gosch_ens <- names(rownames(e4)[!is.na(names(rownames(e4)))])
length(rownames(e4)) - length(gosch_ens)
# [1] 1849
gosch_ens_anno <- getBM(attributes = get, values = gosch_ens, 
                        mart = mart, filters = "ensembl_gene_id")
gosch_data <- gosch_ens_anno %>% 
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
dim(gosch_data)
# [1] 14488     9

### By synonyms ----
gosch_syn <- unname(rownames(e4)[is.na(names(rownames(e4)))])
length(gosch_syn)
# [1] 1849
gosch_syn_anno <- getBM(attributes = c(get, "external_synonym"), values = gosch_syn, 
                        mart = mart, filters = "external_synonym")
dim(gosch_syn_anno)
# [1] 3769   10
length(unique(gosch_syn_anno$ensembl_gene_id))
# [1] 467
gosch_data2 <- gosch_syn_anno %>% 
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
dim(gosch_data2)
# [1] 467   9
gosch_data <- rbind(gosch_data, gosch_data2)
dim(gosch_data)
# [1] 14955     9
length(rownames(e4)) - dim(gosch_data)[1]
# [1] 1434 not found.
write.csv(gosch_data, "./data/annotate/gosch.csv")
## Add ensembl IDs to e4 names. chatgpt solution
# Step 1: Identify indices of NA rownames
na_indices <- which(is.na(attr(rownames(e4), "names")))

# Step 2: Extract the corresponding gene symbols from the rows with NA
na_symbols <- rownames(e4)[na_indices]

# Create a lookup vector for new IDs
new_ids <- gosch_syn_anno$ensembl_gene_id
names(new_ids) <- gosch_syn_anno$external_synonym

# Step 3: Replace NA symbols with their corresponding new IDs
for (i in seq_along(na_indices)) {
  symbol <- na_symbols[i]
  if (symbol %in% names(new_ids)) {
    new_id <- new_ids[symbol]
    attr(rownames(e4), "names")[na_indices[i]] <- new_id
  }
}
rownames(e4)[unname(rownames(e4)) %in% gosch_syn_anno$external_synonym]
head(e4)
class(e4)
write_rds(e4, "./data/expr/gosch.rds")

## LaRocca ----
### By ensembl IDs ----
larocca_ens <- names(rownames(e11)[!is.na(names(rownames(e11)))])
length(rownames(e11)) - length(larocca_ens)
# [1] 1764
larocca_ens_anno <- getBM(attributes = get, values = larocca_ens, 
                        mart = mart, filters = "ensembl_gene_id")
larocca_data <- larocca_ens_anno %>% 
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
dim(larocca_data)
# [1] 13738     9

### By synonyms ----
larocca_syn <- unname(rownames(e11)[is.na(names(rownames(e11)))])
length(larocca_syn)
# [1] 1764
larocca_syn_anno <- getBM(attributes = c(get, "external_synonym"), values = larocca_syn, 
                        mart = mart, filters = "external_synonym")
dim(larocca_syn_anno)
# [1] 1859   10
length(unique(larocca_syn_anno$ensembl_gene_id))
# [1] 202
larocca_data2 <- larocca_syn_anno %>% 
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
dim(larocca_data2)
# [1] 202   9
larocca_data <- rbind(larocca_data, larocca_data2)
dim(larocca_data)
# [1] 13940     9
length(rownames(e11)) - dim(larocca_data)[1]
# [1] 1817 not found.
write.csv(larocca_data, "./data/annotate/larocca.csv")
## Add ensembl IDs to e11 names. chatgpt solution
# Step 1: Identify indices of NA rownames
na_indices <- which(is.na(attr(rownames(e11), "names")))

# Step 2: Extract the corresponding gene symbols from the rows with NA
na_symbols <- rownames(e11)[na_indices]

# Create a lookup vector for new IDs
new_ids <- larocca_syn_anno$ensembl_gene_id
names(new_ids) <- larocca_syn_anno$external_synonym

# Step 3: Replace NA symbols with their corresponding new IDs
for (i in seq_along(na_indices)) {
  symbol <- na_symbols[i]
  if (symbol %in% names(new_ids)) {
    new_id <- new_ids[symbol]
    attr(rownames(e11), "names")[na_indices[i]] <- new_id
  }
}
rownames(e11)[unname(rownames(e11)) %in% larocca_syn_anno$external_synonym]
head(e11)
class(e11)
write_rds(e11, "./data/expr/larocca.rds")