## Data for RShiny-app

## File order
study_names <- c("Gomez-Carballa", "Meaburn Day 1", "Meaburn Day 2",
                 "Gosch", "Obermoser C1 Vein", "Obermoser C1 Finger",
                 "Obermoser C2 Vein", "Obermoser C2 Finger",
                 "Bhasin", "Rusch", "LaRocca", "Karlovich Batch 1", "Karlovich Batch 2")

## Variation tables
variation_tables <- list()
# files <- c("1variation.csv", "2variation.csv", "3variation.csv", "4variation.csv",
#            "5variation.csv", "6variation.csv", "7variation.csv", "8variation.csv",
#            "9variation.csv", "10variation.csv", "11variation.csv", "12variation.csv", "13variation.csv")
files <- c("1variation.parquet", "2variation.parquet", "3variation.parquet", "4variation.parquet",
           "5variation.parquet", "6variation.parquet", "7variation.parquet", "8variation.parquet",
           "9variation.parquet", "10variation.parquet", "11variation.parquet", "12variation.parquet", "13variation.parquet")
for (i in 1:length(files)) {
  variation_tables[[i]] <- read_parquet(file = paste("./data/variation/", files[i], sep = ""))
  names(variation_tables)[i] <- study_names[i]
}
names(variation_tables[[1]]) <- c("Ensembl ID", "Within Variation", "Total Variation", "Rs", "Average Expression",
                                  "Subject VP", "Age-sex VP", "Age VP", "Age-time VP", "Sex VP", "Sex-time VP", "Time VP",
                                  "Residual VP", "Repeatability", "Genetic Variance", "Residual Variance", "Symbol")
names(variation_tables[[2]]) <- c("Probe ID", "Within Variation", "Total Variation", "Rs", "Average Expression",
                                  "Subject VP", "Sex VP", "Time VP", "Time-sex VP", "Residual VP", "Repeatability",
                                  "Genetic Variance", "Residual Variance", "Symbol")
names(variation_tables[[3]]) <- c("Probe ID", "Within Variation", "Total Variation", "Rs", "Average Expression",
                                  "Subject VP", "Sex VP", "Time VP", "Time-sex VP", "Residual VP", "Repeatability",
                                  "Genetic Variance", "Residual Variance", "Symbol")
names(variation_tables[[4]]) <- c("Probe ID", "Within Variation", "Total Variation", "Rs", "Average Expression", 
                                  "Subject VP", "Age VP", "Age-male VP", "Age-11h VP", "Age-14h VP", "Age-17h VP", "Age-8h VP",
                                  "Age-20h VP", "Age-23h VP", "Age-5h VP", "Sex VP VP", "Sex-time VP", "Time VP", "Residual VP", 
                                  "Repeatability", "Genetic Variance", "Residual Variance", "Symbol2")
names(variation_tables[[5]]) <- c("Probe ID", "Within Variation", "Total Variation", "Rs", "Average Expression",
                                  "Subject VP", "Treatment VP", "Age VP", "Age-8d VP", "Age-14d VP", "Age-17d VP", "Age-28d VP", 
                                  "Age-35d VP", "Age-21d VP", "Age-3d VP", "Age-7d VP", "Sex VP", "Time-treatment VP", 
                                  "Race VP", "Time VP", "Time-race VP", "Time-sex VP", "Residual VP", "Repeatability",
                                  "Genetic Variance", "Residiual Variance", "Symbol")
names(variation_tables[[6]]) <- c("Probe ID", "Within Variation", "Total Variation", "Rs", "Average Expression",
                                  "Subject VP", "Age VP", "Age-14d VP", "Age-7d VP", "Race VP", "Time-sex VP", "Sex VP", "Time VP",
                                  "Time-race VP", "Time-treatment VP", "Treatment VP", "Residual VP", "Repeatability",
                                  "Genetic Variance", "Residual Variance", "Symbol")
names(variation_tables[[7]]) <- c("Probe ID", "Within Variation", "Total Variation", "Rs", "Average Expression",
                                  "Subject VP", "Age VP", "Age-14d VP", "Age-7.5d VP", "Age-17d VP", "Age-35d VP", "Age-10d VP",
                                  "Age-7d VP", "Age-8d VP", "Race VP", "Time-sex VP", "Sex VP", "Time VP", "Time-race VP", 
                                  "Time-treatment VP", "Treatment VP", "Residual VP", "Repeatability","Genetic Variance",
                                  "Residual Variance", "Symbol")
names(variation_tables[[8]]) <- c("Probe ID", "Within Variation", "Total Variation", "Rs", "Average Expression",
                                  "Subject VP", "Age VP", "Age-7d15h VP", "Age-7d12h VP", "Age-7d9h VP", "Age-8.5d VP", "Age-7d6h",
                                  "Age-9d VP", "Age7d1.5h VP", "Age-9d VP", "Age7d3h VP", "Age-7d VP", "Time-treatment VP", "Sex VP",
                                  "Race VP", "Time VP", "Time-race VP", "Time-sex VP", "Treatment VP", "Residual VP", "Repeatability",
                                  "Genetic Variance", "Residual Variance", "Symbol")
names(variation_tables[[9]]) <- c("Probe ID", "Within Variation", "Total Variation", "Rs", "Average Expression",
                                  "Subject VP", "Sex VP", "Time VP", "Time-sex VP",
                                  "Residual VP", "Repeatability","Genetic Variance", "Residual Variance", "Symbol")
names(variation_tables[[10]]) <- c("Probe ID", "Within Variation", "Total Variation", "Rs", "Average Expression",
                                   "Subject VP", "Age VP", "Time-race VP", "Time-ethnicity VP", "Time VP", "Race VP", "Ethnicity VP",
                                   "Residual VP", "Repeatability","Genetic Variance", "Residual Variance", "Symbol")
names(variation_tables[[11]]) <- c("Symbol", "Within Variation", "Total Variation", "Rs", "Average Expression",
                                   "Subject VP", "Response VP", "Time-response VP", "Time VP", "Residual VP", "Repeatability", 
                                   "Genetic Variance", "Residual Variance", "Symbol2")
names(variation_tables[[12]]) <- c("Probe ID", "Within Variation", "Total Variation SD", "Rs", "Average Expression", "Subject VP", 
                                   "Age VP", "Age-28d VP", "Age-14d VP", "Time-sex VP", "Sex VP", "Time VP", "Residual VP", 
                                   "Repeatability", "Genetic Variance", "Residual Variance", "Symbol")
names(variation_tables[[13]]) <- c("Probe ID", "Within Variation", "Total Variation SD", "Rs", "Average Expression", "Subject VP", 
                                   "Age-day90 VP", "Age VP", "Sex VP", "Time VP", "Time-sex VP", "Residual VP", 
                                   "Repeatability", "Genetic Variance", "Residual Variance", "Symbol")

variation_tables[[1]] <- variation_tables[[1]][,c("Ensembl ID", "Symbol", "Repeatability", "Genetic Variance", "Within Variation", 
                                                  "Total Variation", "Average Expression", "Subject VP", "Age-sex VP", "Age VP", 
                                                  "Age-time VP", "Sex VP", "Sex-time VP", "Time VP", "Residual VP", "Residual Variance", "Rs")]
variation_tables[[2]] <- variation_tables[[2]][,c("Probe ID", "Symbol", "Repeatability", "Genetic Variance", "Within Variation",
                                                  "Total Variation", "Average Expression", "Subject VP", "Sex VP", "Time VP", "Time-sex VP",
                                                  "Residual VP", "Residual Variance", "Rs")]
variation_tables[[3]] <- variation_tables[[3]][,c("Probe ID", "Symbol", "Repeatability", "Genetic Variance", "Within Variation",
                                                  "Total Variation", "Average Expression", "Subject VP", "Sex VP", "Time VP", 
                                                  "Time-sex VP", "Residual VP",  "Residual Variance", "Rs")]
variation_tables[[4]] <- variation_tables[[4]][,c("Probe ID", "Symbol2", "Repeatability", "Genetic Variance", "Within Variation", 
                                                  "Total Variation",  "Average Expression", "Subject VP", "Age VP", "Age-male VP",
                                                  "Age-11h VP", "Age-14h VP", "Age-17h VP", "Age-8h VP", "Age-20h VP", "Age-23h VP", 
                                                  "Age-5h VP", "Sex VP VP", "Sex-time VP", "Time VP", "Residual VP", 
                                                  "Residual Variance", "Rs")]
variation_tables[[5]] <- variation_tables[[5]][,c("Probe ID", "Symbol", "Repeatability", "Genetic Variance", "Within Variation",
                                                  "Total Variation", "Average Expression", "Subject VP", "Treatment VP", "Age VP",
                                                  "Age-8d VP", "Age-14d VP", "Age-17d VP", "Age-28d VP", "Age-35d VP", "Age-21d VP",
                                                  "Age-3d VP", "Age-7d VP", "Sex VP", "Time-treatment VP", "Race VP", "Time VP",
                                                  "Time-race VP", "Time-sex VP", "Residual VP",  "Residiual Variance", "Rs")]
variation_tables[[6]] <- variation_tables[[6]][,c("Probe ID", "Symbol", "Repeatability", "Genetic Variance", "Within Variation",
                                                  "Total Variation", "Average Expression", "Subject VP", "Age VP", "Age-14d VP", 
                                                  "Age-7d VP", "Race VP", "Time-sex VP", "Sex VP", "Time VP", "Time-race VP", 
                                                  "Time-treatment VP", "Treatment VP", "Residual VP",  "Residual Variance", "Rs")]
variation_tables[[7]] <- variation_tables[[7]][,c("Probe ID", "Symbol", "Repeatability","Genetic Variance", "Within Variation",
                                                  "Total Variation", "Average Expression", "Subject VP", "Age VP", "Age-14d VP",
                                                  "Age-7.5d VP", "Age-17d VP", "Age-35d VP", "Age-10d VP", "Age-7d VP", "Age-8d VP",
                                                  "Race VP", "Time-sex VP", "Sex VP", "Time VP", "Time-race VP", "Time-treatment VP",
                                                  "Treatment VP", "Residual VP", "Residual Variance", "Rs")]
variation_tables[[8]] <- variation_tables[[8]][,c("Probe ID", "Symbol", "Repeatability", "Genetic Variance", "Within Variation",
                                                  "Total Variation", "Average Expression", "Subject VP", "Age VP", "Age-7d15h VP",
                                                  "Age-7d12h VP", "Age-7d9h VP", "Age-8.5d VP", "Age-7d6h", "Age-9d VP", "Age7d1.5h VP",
                                                  "Age-9d VP", "Age7d3h VP", "Age-7d VP", "Time-treatment VP", "Sex VP", "Race VP", 
                                                  "Time VP", "Time-race VP", "Time-sex VP", "Treatment VP", "Residual VP",
                                                  "Residual Variance", "Rs")]
variation_tables[[9]] <- variation_tables[[9]][,c("Probe ID", "Symbol", "Repeatability","Genetic Variance", "Within Variation",
                                                  "Total Variation", "Average Expression", "Subject VP", "Sex VP", "Time VP", 
                                                  "Time-sex VP", "Residual VP",  "Residual Variance", "Rs")]
variation_tables[[10]] <- variation_tables[[10]][,c("Probe ID", "Symbol", "Repeatability", "Genetic Variance", "Within Variation",
                                                    "Total Variation", "Average Expression", "Subject VP", "Age VP", "Time-race VP", 
                                                    "Time-ethnicity VP", "Time VP", "Race VP", "Ethnicity VP", "Residual VP", 
                                                    "Residual Variance", "Rs")]
variation_tables[[11]] <- variation_tables[[11]][,c("Symbol", "Repeatability", "Genetic Variance", "Within Variation", "Total Variation",
                                                    "Average Expression", "Subject VP", "Response VP", "Time-response VP", "Time VP",
                                                    "Residual VP",  "Residual Variance", "Rs", "Symbol2")]
variation_tables[[12]] <- variation_tables[[12]][,c("Probe ID", "Symbol", "Repeatability", "Genetic Variance", "Within Variation",
                                                    "Total Variation SD", "Average Expression", "Subject VP", "Age VP", "Age-28d VP",
                                                    "Age-14d VP", "Time-sex VP", "Sex VP", "Time VP", "Residual VP", "Residual Variance", "Rs")]
variation_tables[[13]] <- variation_tables[[13]][,c("Probe ID", "Symbol", "Repeatability", "Genetic Variance", "Within Variation",
                                                    "Total Variation SD", "Average Expression", "Subject VP", "Age-day90 VP", "Age VP",
                                                    "Sex VP", "Time VP", "Time-sex VP", "Residual VP", "Residual Variance", "Rs")]


## Stability and dynamic data
stable <- read.csv("./data/characteristic_scores.csv")
dynamic <-  read.csv("./data/flexible_gene_scores.csv")

## Phenotype data
pheno_data <- list()
files <- list.files("./data/phenotypes/")
for (i in 1:length(files)) {
  pheno_data[[i]] <- read.csv(file = paste("./data/phenotypes/", files[i], sep = ""))
  names(pheno_data)[i] <- study_names[i]
}
for (i in 1:length(pheno_data)) {
  if ("participant" %in% colnames(pheno_data[[i]]) | "id" %in% colnames(pheno_data[[i]])) {
    colnames(pheno_data[[i]])[grep("participant|id", colnames(pheno_data[[i]]))] <- "subject"
  }
}

## eQTL data
# gtex <- read_rds("./data/eqtls/gtex_variants_and_phenotypes_whole_blood.rds")
gtex <- nanoparquet::read_parquet("./data/eqtls/gtex_variants_and_phenotypes_whole_blood.parquet")
gtexColumnNames <- c(
  '<span title="GENCODE gene name for blood RNA">Symbol of blood RNA</span>',
  '<span title="GENCODE/Ensembl gene ID for blood RNA">Ensembl ID</span>',
  '<span title="Reference SNP cluster ID in dbSNP database">rsID of eQTL</span>',
  '<span title="variant ID in the format {chr}_{pos_first_ref_base}_{ref_seq}_{alt_seq}_b38">GTEx variant ID</span>',
  '<span title="Minor allele frequency observed in the set of donors for GTEx whole blood tissue">MAF in GTEx WB</span>',
  '<span title="Regression slope, i.e. the effect size on expression level due to the variant">eQTL slope</span>',
  '<span title="eQTL nominal p-value">eQTL nominal p_val</span>',
  '<span title="Beta-approximated permutation p-value for the gene">eQTL beta p_value</span>',
  '<span title="Disease or trait examined in study submitted to GWAS Catalog">GWAS Trait</span>',
  '<span title=" Gene(s) mapped to the strongest SNP. If the SNP is located within a gene, that gene is listed. If the SNP is located within multiple genes, these genes are listed separated by commas. If the SNP is intergenic, the upstream and downstream genes are listed, separated by a hyphen.">GWAS Catalog Mapping</span>',
  '<span title="Gene(s) reported by catalog submission author">Reported Mapping</span>',
  '<span title="Reported odds ratio or beta-coefficient associated with strongest SNP risk allele, i.e. the effect size of the risk allele.">OR or BETA</span>',
  '<span title="Reported p-value for strongest SNP risk allele (linked to dbGaP Association Browser). Note that p-values are rounded to 1 significant digit (for example, a published p-value of 4.8 x 10-7 is rounded to 5 x 10-7)">GWAS p_value</span>',
  '<span title="Provides information on a variant’s predicted most severe functional effect from Ensembl">Genomic Context</span>',
  '<span title="A score ranging from 0-7. Genes that score higher are generally expressed stably over time for an individual, but expression levels vary between individuals.">Trait no. studies</span>',
  '<span title="A score ranging from 0-6. Genes that score higher were called differentially expressed over time in more studies (p < 0.05 using limma comparing baseline and follow-up gene expression), indicating varying expression levels over time.">State no. studies</span>'
)

# Gene graph
gtex_gene <- as.data.frame(table(gtex$`Symbol of blood RNA`))
gtex_gene <- gtex_gene[order(gtex_gene$Freq, decreasing = T),]
gtex_gene$index = 1:nrow(gtex_gene)

# rsID graph
gtex_rsID <- as.data.frame(table(gtex$`rsID of eQTL`))
gtex_rsID <- gtex_rsID[order(gtex_rsID$Freq, decreasing = T),]
gtex_rsID$index = 1:nrow(gtex_rsID)

# Trait graph
gtex_trait <- as.data.frame(table(gtex$`GWAS Trait`))
gtex_trait <- gtex_trait[order(gtex_trait$Freq, decreasing = T),]
gtex_trait$index = 1:nrow(gtex_trait)

# Gene Trait Graph
gtex_gene_trait <- gtex[!is.na(gtex$`GWAS Trait`),]
gtex_gene_trait <- as.data.frame(table(gtex_gene_trait$`Symbol of blood RNA`))
gtex_gene_trait <- gtex_gene_trait[order(gtex_gene_trait$Freq, decreasing = T),]
gtex_gene_trait$index = 1:nrow(gtex_gene_trait)

## ENSEMBL annotation data
anno_data <- list()
files <- list.files("./data/annotate/", full.names = T)
for (i in 1:length(files)) {
  anno_data[[i]] <- read_parquet(files[[i]])
  anno_data[[i]] <- anno_data[[i]][,-1]
}
names(anno_data) <- c("affy", "gomez", "gosch", "ilmn", "larocca")