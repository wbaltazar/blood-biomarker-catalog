biomaRt::listMarts()
mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL")
View(biomaRt::listDatasets(mart))
dataset <- biomaRt::useDataset(dataset = "hsapiens_gene_ensembl", mart = mart)
possible_attributes <- biomaRt::listAttributes(mart = dataset)
View(possible_attributes)
unique(possible_attributes$page)
View(possible_attributes[possible_attributes$page == "feature_page",])
attributes <- c("description", "transcript_gencode_basic", "transcript_appris", "external_gene_name",
                "external_gene_source", "external_transcript_name", "external_transcript_source_name",
                "gene_biotype", "transcript_biotype")
possible_filters <- biomaRt::listFilters(mart = dataset)
View(possible_filters)
filters <- c("external_gene_name", "external_transcript_name", "external_synonym", "affy_hg_u133_plus_2",
             "illumina_humanwg_6_v3", "ensembl_gene_id")
mapIds(hgu133plus2.db, "TUBB2A", "PROBEID", "SYMBOL")
# 204141_at
mapIds(illuminaHumanv3.db, "TUBB2A", "PROBEID", "SYMBOL")
# ILMN_2038775

gene_info(gene = "204141_at", filter = "affy_hg_u133_plus_2") # works
gene_info(gene = "ILMN_2038775", filter = "illumina_humanwg_6_v3") # works
head(rownames(e11))
# ENSG00000121410 ENSG00000268895 ENSG00000175899 ENSG00000245105 ENSG00000094914 ENSG00000081760 
# "A1BG"      "A1BG.AS1"           "A2M"       "A2M.AS1"          "AAAS"          "AACS" 
gene_info(gene = "ENSG00000268895", filter = "ensembl_gene_id")
head(rownames(e11)[is.na(names(rownames(e11)))])
# <NA>       <NA>       <NA>       <NA>       <NA>       <NA> 
#   "ADCK3"    "ADCK4"   "ADRBK1"   "ADRBK2" "AHCTF1P1" "ALOX15P1" 
gene_info(gene = "ADCK3", filter = "external_synonym") # gene name doesn't work. ADCK3 is another name for COQ8A
# description transcript_gencode_basic transcript_appris external_gene_name
# 1 coenzyme Q8A [Source:HGNC Symbol;Acc:HGNC:16812]            GENCODE basic                                COQ8A
# 2 coenzyme Q8A [Source:HGNC Symbol;Acc:HGNC:16812]                                                         COQ8A
# 3 coenzyme Q8A [Source:HGNC Symbol;Acc:HGNC:16812]            GENCODE basic        principal1              COQ8A
tail(rownames(e11)[is.na(names(rownames(e11)))])
# <NA>                        <NA>                        <NA>                        <NA> 
#   "X8_LINE.CR1.LINE"           "X9_LINE.L1.LINE"    "Zaphod2.hAT.Tip100.DNA"    "Zaphod3.hAT.Tip100.DNA" 
# <NA>                        <NA> 
#   "Zaphod.hAT.Tip100.DNA" "hAT.N1_Mam.hAT.Tip100.DNA" 
head(rownames(e11)[grep("^LOC", rownames(e11))])
gene_info(gene = "646214", filter = "entrezgene_trans_name") # found this out on ucsc genome browser
View(rownames(e11)[is.na(names(rownames(e11)))]) # retrotransposons / repeatmasker elements, some noncoding elements
View(rownames(e4)[is.na(names(rownames(e4)))])
gene_info(gene = "C1orf233", filter = "external_synonym")


## Was gonna use in the app, but I changed my mind.
# Return gene information from biomaRt
output$geneInfo <- renderDataTable({
  # Check for BeadChip or Microarray conditions
  if (isBeadChip()) {
    return(gene_info(gene = input$probeInput, filter = "illumina_humanwg_6_v3"))
  } else if (isMicroarray()) {
    return(gene_info(gene = input$probeInput, filter = "affy_hg_u133_plus_2"))
  }
  
  # Function to handle gene search within data frames
  search_gene <- function(df, geneName) {
    gene_row <- grep(geneName, rownames(df))
    if (length(gene_row) > 0) {
      return(rownames(df)[gene_row])
    }
    return(NULL)
  }
  
  # Handle specific pData cases
  if (input$pData == "pheno_gomez") {
    gene <- search_gene(e1, input$geneName)
  } else if (input$pData == "pheno_gosch") {
    gene <- search_gene(e4, input$geneName)
  } else if (input$pData == "pheno_larocca") {
    gene <- search_gene(e11, input$geneName)
  } else {
    gene <- NULL
  }
  
  # If gene is found, perform biomaRt query
  if (!is.na(names(gene))) {
    return(gene_info(gene = names(gene), filter = "ensembl_gene_id"))
  } else if (!is.null(gene)) {
    result <- NULL
    for (filter in c("external_gene_name", "external_synonym")) {
      result <- gene_info(gene = unname(gene), filter = filter)
      if (nrow(result) > 0) {
        break
      }
    }
    
    if (!is.null(result) && nrow(result) > 0) {
      return(result)
    } else {
      return(data.frame(Unfortunately = "no info was found."))
    }
  }
  
  # Default case if no gene is found
  return(data.frame(Unfortunately = "no info was found."))
})