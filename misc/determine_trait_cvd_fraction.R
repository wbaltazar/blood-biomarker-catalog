library(tidyverse)
library(nanoparquet)
tab <- read.csv("~/Downloads/filtered_gtex_gwas_table (3).csv") # CVD file. For CAD, use (4) instead of (3).
sym <- tab$Symbol.of.blood.RNA
gtex <- read_parquet("~/Desktop/work_repo/Box organization/1results/RShiny-application/data/eqtls/gtex_variants_and_phenotypes_whole_blood.parquet")
cvd_spec_table <- gtex[gtex$`Symbol of blood RNA` %in% sym,]
cvd_spec_table <- cvd_spec_table %>% select(`Symbol of blood RNA`, `GWAS Trait`) %>% group_by(`Symbol of blood RNA`)
cvd_spec_table <- cvd_spec_table %>% reframe(traits = unique(`GWAS Trait`)) %>% remove_missing()
sum(is.na(cvd_spec_table$traits))
# [1] 0

cvd_traits <- read.csv("~/Desktop/listofterms.csv")
cvd_traits <- cvd_traits$Coronary.artery.disease
cvd_traits[607] <- "Coronary artery disease"
cvd_traits <- cvd_traits[1:607]

cvd_traits_only <- cvd_spec_table[cvd_spec_table$traits %in% cvd_traits,]

cvd <- data.frame(table(cvd_traits_only$`Symbol of blood RNA`)) 
all <- data.frame(table(cvd_spec_table$`Symbol of blood RNA`))
table <- merge(cvd, all - 1, by = "Var1", all.y = T) # all - 1 removes cardiovascular disease from all traits
table$Freq.x <- replace_na(table$Freq.x, replace = 0)
table$Percent <- signif(x = table$Freq.x / table$Freq.y * 100, digits = 2)
colnames(table) <- c("RNA", "CVD traits", "All traits", "Percent")
range(table$Percent)
table[which.max(table$Percent),]
#     RNA CVD traits All traits Percent
# 64 RPS6          6          9       67
write.csv(table, "~/Desktop/cvd_percent_90_genes.csv")



tab <- read.csv("~/Downloads/filtered_gtex_gwas_table (4).csv")
sym <- tab$Symbol.of.blood.RNA
gtex <- read_parquet("~/Desktop/work_repo/Box organization/1results/RShiny-application/data/eqtls/gtex_variants_and_phenotypes_whole_blood.parquet")
cad_spec_table <- gtex[gtex$`Symbol of blood RNA` %in% sym,]
cad_spec_table <- cad_spec_table %>% select(`Symbol of blood RNA`, `GWAS Trait`) %>% group_by(`Symbol of blood RNA`)
cad_spec_table <- cad_spec_table %>% reframe(traits = unique(`GWAS Trait`)) %>% remove_missing()
sum(is.na(cad_spec_table$traits))
# [1] 0
all_traits <- unique(cad_spec_table$traits)
write.csv(all_traits[str_detect(all_traits, "age|gender|lipid|smok|diabetes|cholesterol|systolic|diastolic
                      |hypertension|coronary artery|adiposity|adipose|diet|physical activity|alcohol|
                      triglyceride|glucose|calori|lifestyle|consumption|exercise|atherosclerosis|
                      homocysteine|sleep apnea|homocystinuria|left ventric|insulin")],
            file = "~/Desktop/work_repo/github/worked_example/listofterms_CAD_new.csv")
# Manually prune terms which are not risk factors for CAD from the csv file. Then, re-read the data in.
cad_traits <- read.csv("~/Desktop/work_repo/github/worked_example/listofterms_CAD.csv")
cad_traits <- cad_traits$Traits
head(cad_traits)
# [1] "BMI and adiposity (confirmatory factor analysis Factor 7)"                      
# [2] "Body fat percentage"                                                            
# [3] "Total cholesterol levels (MTAG)"                                                
# [4] "Body mass index (joint analysis main effects and physical activity interaction)"
# [5] "Visceral adipose tissue volumes to abdominal adipose tissue volumes ratio"      
# [6] "BMI in non-smokers"   
cad_traits <- c(cad_traits, "Coronary artery disease")

cad_traits_only <- cad_spec_table[cad_spec_table$traits %in% cad_traits,]
# SAP130, DOCK5 CAD only!
cad_spec_table[cad_spec_table$`Symbol of blood RNA` %in% c("SAP130", "DOCK5"),]
cad_traits_only[cad_traits_only$`Symbol of blood RNA` %in% c("SAP130", "DOCK5"),]

cad <- data.frame(table(cad_traits_only$`Symbol of blood RNA`)) 
all <- data.frame(table(cad_spec_table$`Symbol of blood RNA`))
table <- merge(cad, all, by = "Var1", all.y = T) 
table$Freq.x <- replace_na(table$Freq.x, replace = 0)
table$Percent <- signif(x = table$Freq.x / table$Freq.y * 100, digits = 2)
colnames(table) <- c("RNA", "CAD traits", "All traits", "Percent")
range(table$Percent)
table[which.max(table$Percent),]
#      RNA CAD traits All traits Percent
# 44 DOCK5          1          1     100
write.csv(table, "~/Desktop/work_repo/github/worked_example/cad_percent_179_genes.csv")
nrow(table)
# [1] 179
length(cad_traits)
# [1] 185