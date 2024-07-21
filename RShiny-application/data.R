## Data for RShiny-app

## File order
study_names <- c("Gomez-Carballa", "Meaburn Day 1", "Meaburn Day 2",
                 "Gosch", "Obermoser C1 Vein", "Obermoser C1 Finger",
                 "Obermoser C2 Vein", "Obermoser C2 Finger",
                 "Bhasin", "Rusch", "LaRocca")

## Variation tables
variation_tables <- list()
files <- c("1variation.csv", "2variation.csv", "3variation.csv", "4variation.csv",
           "5variation.csv", "6variation.csv", "7variation.csv", "8variation.csv",
           "9variation.csv", "10variation.csv", "11variation.csv")
for (i in 1:length(files)) {
  variation_tables[[i]] <- read.csv(file = paste("./data/variation/", files[i], sep = ""))
  names(variation_tables)[i] <- study_names[i]
}
names(variation_tables[[1]]) <- c("Ensembl ID", "Within Variation", "Total Variation", "Rs", "Average Expression",
                                  "Subject VP", "Age VP", "Sex VP", "Time VP", "Residual VP", "Repeatability",
                                  "Genetic Variance", "Residual Variance", "Symbol")
names(variation_tables[[2]]) <- c("Probe ID", "Within Variation", "Total Variation", "Rs", "Average Expression",
                                  "Subject VP", "Sex VP", "Time VP", "Residual VP", "Repeatability",
                                  "Genetic Variance", "Residual Variance", "Symbol")
names(variation_tables[[3]]) <- c("Probe ID", "Within Variation", "Total Variation", "Rs", "Average Expression",
                                  "Subject VP", "Sex VP", "Time VP", "Residual VP", "Repeatability",
                                  "Genetic Variance", "Residual Variance", "Symbol")
names(variation_tables[[4]]) <- c("Symbol", "Within Variation", "Total Variation", "Rs", "Average Expression",
                                  "Subject VP", "Age VP", "Sex VP", "Time VP", "Residual VP", "Repeatability",
                                  "Genetic Variance", "Residual Variance", "Symbol2")
names(variation_tables[[5]]) <- c("Probe ID", "Within Variation", "Total Variation", "Rs", "Average Expression",
                                  "Treatment VP", "Subject VP", "Age VP", "Time VP", "Sex VP", "Race VP",
                                  "Residual VP", "Repeatability","Genetic Variance", "Residual Variance", "Symbol")
names(variation_tables[[6]]) <- c("Probe ID", "Within Variation", "Total Variation", "Rs", "Average Expression",
                                  "Treatment VP", "Subject VP", "Age VP", "Time VP", "Sex VP", "Race VP",
                                  "Residual VP", "Repeatability","Genetic Variance", "Residual Variance", "Symbol")
names(variation_tables[[7]]) <- c("Probe ID", "Within Variation", "Total Variation", "Rs", "Average Expression",
                                  "Treatment VP", "Subject VP", "Age VP", "Time VP", "Sex VP", "Race VP",
                                  "Residual VP", "Repeatability","Genetic Variance", "Residual Variance", "Symbol")
names(variation_tables[[8]]) <- c("Probe ID", "Within Variation", "Total Variation", "Rs", "Average Expression",
                                  "Treatment VP", "Subject VP", "Age VP", "Time VP", "Sex VP", "Race VP",
                                  "Residual VP", "Repeatability","Genetic Variance", "Residual Variance", "Symbol")
names(variation_tables[[9]]) <- c("Probe ID", "Within Variation", "Total Variation", "Rs", "Average Expression",
                                   "Subject VP", "Sex VP", "Time VP", 
                                   "Residual VP", "Repeatability","Genetic Variance", "Residual Variance", "Symbol")
names(variation_tables[[10]]) <- c("Probe ID", "Within Variation", "Total Variation", "Rs", "Average Expression",
                                   "Subject VP", "Age VP", "Time VP", "Ethnicity VP", "Race VP","Residual VP",
                                   "Repeatability","Genetic Variance", "Residual Variance", "Symbol")
names(variation_tables[[11]]) <- c("Symbol", "Within Variation", "Total Variation", "Rs", "Average Expression",
                                   "Subject VP", "Time VP", "Residual VP", "Repeatability", "Genetic Variance",
                                   "Residual Variance", "Symbol2")


## Stability and dynamic data
stable <- read.csv("./data/stable_gene_filter_count.csv")
dynamic <-  read.csv("./data/deg_counts.csv")
housekeeping <- read.csv("./data/housekeeping_scores.csv")

## Phenotype data
pheno_data <- list()
files <- list.files("./data/phenotypes/")
for (i in 1:length(files)) {
  pheno_data[[i]] <- read.csv(file = paste("./data/phenotypes/", files[i], sep = ""))
  names(pheno_data)[i] <- study_names[i]
}

## eQTL data
gtex <- read_rds("./data/eqtls/gtex_variants_and_phenotypes_whole_blood_28062024.rds")
gtex$`Stable-polymorphic score` <- replace_na(gtex$`Stable-polymorphic score`, 0)
gtex$`Flexibility score` <- replace_na(gtex$`Flexibility score`, 0)
gtex$`Housekeeping score` <- replace_na(gtex$`Housekeeping score`, 0)
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
  '<span title="Denotes whether SNP is in intergenic region (0 = no; 1 = yes)">Intergenic?</span>',
  '<span title="A score ranging from 1-22. Genes that score higher are generally expressed stably over time for an individual, but expression levels vary between individuals.">Stable-polymorphic score</span>',
  '<span title="A score ranging from 1-6. Genes that score higher were called differentially expressed over time in more studies, indicating varying expression levels over time.">Flexibility score</span>',
  '<span title="A score ranging from 1-39. Genes that score higher show little variation across time or across individuals.">Housekeeping score</span>'
)

gtexplot <- as_tibble(as.data.frame(table(gtex$`Symbol of blood RNA`)))
gtexplot2 <- as_tibble(as.data.frame(table(gtex$`GWAS Trait`)))

## ENSEMBL annotation data
anno_data <- list()
files <- list.files("./data/annotate/", full.names = T)
for (i in 1:length(files)) {
  anno_data[[i]] <- read.csv(files[[i]])
  anno_data[[i]] <- anno_data[[i]][,-1]
}
names(anno_data) <- c("affy", "gomez", "gosch", "ilmn", "larocca")