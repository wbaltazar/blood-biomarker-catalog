setwd("~/Desktop/work_repo/Box organization/1results/")
library(tidyverse)
sample <- read_rds("./shiny-data-new/data/eqtls/gtex_variants_and_phenotypes_whole_blood.rds")
sample %>% 
  filter(Symbol == "MYOM2", !is.na(`DISEASE/TRAIT`))
# rs_id_dbSNP151_GRCh38p7           variant_id         gene_id tss_distance ma_samples ma_count       maf pval_nominal    slope
# 1              rs11987758 chr8_2180583_G_A_b38 ENSG00000036448       135543        161      170 0.1268660  3.59994e-11 0.424966
# 2              rs17685410 chr8_2165853_A_G_b38 ENSG00000036448       120813        327      383 0.2858210  8.79864e-10 0.289132
# 3              rs36098404 chr8_2134258_G_A_b38 ENSG00000036448        89218        113      121 0.0902985  7.43704e-29 0.861834
# 4              rs36098404 chr8_2134258_G_A_b38 ENSG00000036448        89218        113      121 0.0902985  7.43704e-29 0.861834
# 5               rs4876199 chr8_2229242_T_G_b38 ENSG00000036448       184202        196      209 0.1559700  5.83200e-05 0.243290
# slope_se pval_nominal_threshold min_pval_nominal   pval_beta Symbol                            DISEASE/TRAIT
# 1 0.0630019            6.36103e-05      2.06347e-96 1.03386e-84  MYOM2 Amyotrophic lateral sclerosis (sporadic)
# 2 0.0464119            6.36103e-05      2.06347e-96 1.03386e-84  MYOM2                  Toenail selenium levels
# 3 0.0733183            6.36103e-05      2.06347e-96 1.03386e-84  MYOM2                     Opioid use cessation
# 4 0.0733183            6.36103e-05      2.06347e-96 1.03386e-84  MYOM2                     Opioid use cessation
# 5 0.0600949            6.36103e-05      2.06347e-96 1.03386e-84  MYOM2                Heschl's gyrus morphology
#                  REPORTED GENE(S)       MAPPED_GENE                            CONTEXT INTERGENIC
# 1                    LOC100286951 MYOM2 - LINC03021                     intron_variant          1
# 2 MYOM2, CSMD1, KBTBD11, ARHGEF10 MYOM2 - LINC03021                     intron_variant          1
# 3                           MYOM2             MYOM2 non_coding_transcript_exon_variant          0
# 4                          SEPT14             MYOM2 non_coding_transcript_exon_variant          0
# 5                           MYOM2 MYOM2 - LINC03021                 intergenic_variant          1
sample %>% 
  filter(Symbol == "TUBB2A", !is.na(`DISEASE/TRAIT`))
# rs_id_dbSNP151_GRCh38p7               variant_id         gene_id tss_distance ma_samples ma_count       maf pval_nominal     slope
# 1               rs13204642     chr6_3143926_T_C_b38 ENSG00000137267       -13600         71       75 0.0559702  4.74326e-08  0.502096
# 2               rs13204642     chr6_3143926_T_C_b38 ENSG00000137267       -13600         71       75 0.0559702  4.74326e-08  0.502096
# 3               rs13206421     chr6_3149037_G_A_b38 ENSG00000137267        -8489         69       73 0.0544776  1.18691e-07  0.491912
# 4              rs145692049 chr6_3149626_TTTTG_T_b38 ENSG00000137267        -7900         70       74 0.0552239  5.13220e-08  0.504141
# 5               rs17548315     chr6_3114713_G_A_b38 ENSG00000137267       -42813         19       20 0.0149254  3.67448e-05  0.697355
# 6               rs17548315     chr6_3114713_G_A_b38 ENSG00000137267       -42813         19       20 0.0149254  3.67448e-05  0.697355
# 7               rs35351985     chr6_3166734_G_A_b38 ENSG00000137267         9208        400      507 0.3783580 3.69770e-102 -0.816365
# 8               rs35351985     chr6_3166734_G_A_b38 ENSG00000137267         9208        400      507 0.3783580 3.69770e-102 -0.816365
# 9                rs3799219     chr6_3150869_A_G_b38 ENSG00000137267        -6657         70       74 0.0552239  5.13220e-08  0.504141
# 10              rs55635148   chr6_3162237_C_CCG_b38 ENSG00000137267         4711        366      465 0.3470150 2.19257e-116 -0.849248
# 11               rs6936244     chr6_3140822_C_T_b38 ENSG00000137267       -16704         96      106 0.0791045  2.31073e-08  0.462101
# 12                rs727260     chr6_3152171_G_A_b38 ENSG00000137267        -5355         70       74 0.0552239  5.13220e-08  0.504141
# 13               rs7765391     chr6_3138265_A_G_b38 ENSG00000137267       -19261        208      247 0.1843280  7.87557e-07 -0.279569
# 14               rs7765391     chr6_3138265_A_G_b38 ENSG00000137267       -19261        208      247 0.1843280  7.87557e-07 -0.279569
# 15               rs9392465     chr6_3162144_C_A_b38 ENSG00000137267         4618        395      498 0.3716420 1.70945e-143 -0.899735
# 16               rs9392465     chr6_3162144_C_A_b38 ENSG00000137267         4618        395      498 0.3716420 1.70945e-143 -0.899735
# 17               rs9392465     chr6_3162144_C_A_b38 ENSG00000137267         4618        395      498 0.3716420 1.70945e-143 -0.899735
# 18               rs9392465     chr6_3162144_C_A_b38 ENSG00000137267         4618        395      498 0.3716420 1.70945e-143 -0.899735
# 19               rs9392465     chr6_3162144_C_A_b38 ENSG00000137267         4618        395      498 0.3716420 1.70945e-143 -0.899735
# slope_se pval_nominal_threshold min_pval_nominal    pval_beta Symbol
# 1  0.0907689             8.8144e-05     1.70945e-143 1.79087e-131 TUBB2A
# 2  0.0907689             8.8144e-05     1.70945e-143 1.79087e-131 TUBB2A
# 3  0.0917703             8.8144e-05     1.70945e-143 1.79087e-131 TUBB2A
# 4  0.0913789             8.8144e-05     1.70945e-143 1.79087e-131 TUBB2A
# 5  0.1677030             8.8144e-05     1.70945e-143 1.79087e-131 TUBB2A
# 6  0.1677030             8.8144e-05     1.70945e-143 1.79087e-131 TUBB2A
# 7  0.0310011             8.8144e-05     1.70945e-143 1.79087e-131 TUBB2A
# 8  0.0310011             8.8144e-05     1.70945e-143 1.79087e-131 TUBB2A
# 9  0.0913789             8.8144e-05     1.70945e-143 1.79087e-131 TUBB2A
# 10 0.0292586             8.8144e-05     1.70945e-143 1.79087e-131 TUBB2A
# 11 0.0816024             8.8144e-05     1.70945e-143 1.79087e-131 TUBB2A
# 12 0.0913789             8.8144e-05     1.70945e-143 1.79087e-131 TUBB2A
# 13 0.0560129             8.8144e-05     1.70945e-143 1.79087e-131 TUBB2A
# 14 0.0560129             8.8144e-05     1.70945e-143 1.79087e-131 TUBB2A
# 15 0.0262396             8.8144e-05     1.70945e-143 1.79087e-131 TUBB2A
# 16 0.0262396             8.8144e-05     1.70945e-143 1.79087e-131 TUBB2A
# 17 0.0262396             8.8144e-05     1.70945e-143 1.79087e-131 TUBB2A
# 18 0.0262396             8.8144e-05     1.70945e-143 1.79087e-131 TUBB2A
# 19 0.0262396             8.8144e-05     1.70945e-143 1.79087e-131 TUBB2A
# DISEASE/TRAIT REPORTED GENE(S)       MAPPED_GENE
# 1                                                                  X-17676 levels             <NA>              BPHL
# 2                                                                  X-17676 levels             <NA>              BPHL
# 3                                         Protein quantitative trait loci (liver)               NR              BPHL
# 4                                                                  X-17676 levels             <NA>              BPHL
# 5  Depressive symptoms x independent stressful life events interaction (2df test)            RIPK1             RIPK1
# 6  Depressive symptoms x independent stressful life events interaction (1df test)            RIPK1             RIPK1
# 7                                                     Mean corpuscular hemoglobin             <NA> TUBB2A - TUBB2BP1
# 8                                                     Red cell distribution width               NR TUBB2A - TUBB2BP1
# 9                                                   Valacyclovir hydrolase levels             <NA>              BPHL
# 10                                                       Mean reticulocyte volume           TUBB2A TUBB2A - TUBB2BP1
# 11                                                                 X-17676 levels             <NA>              BPHL
# 12                                                                 X-17676 levels             <NA>              BPHL
# 13                                                                 X-17676 levels             <NA>              BPHL
# 14                                                  Valacyclovir hydrolase levels             <NA>              BPHL
# 15                                                    Red cell distribution width           TUBB2A TUBB2A - TUBB2BP1
# 16                                                    Mean corpuscular hemoglobin             <NA> TUBB2A - TUBB2BP1
# 17                                                    Red cell distribution width             <NA> TUBB2A - TUBB2BP1
# 18                                                    Mean corpuscular hemoglobin           TUBB2A TUBB2A - TUBB2BP1
# 19                                                    Red cell distribution width           TUBB2A TUBB2A - TUBB2BP1
# CONTEXT INTERGENIC
# 1                      intron_variant          0
# 2                      intron_variant          0
# 3                      intron_variant          0
# 4                      intron_variant          0
# 5                 3_prime_UTR_variant          0
# 6                 3_prime_UTR_variant          0
# 7                  intergenic_variant          1
# 8                  intergenic_variant          1
# 9                      intron_variant          0
# 10          regulatory_region_variant          1
# 11                     intron_variant          0
# 12                     intron_variant          0
# 13 non_coding_transcript_exon_variant          0
# 14 non_coding_transcript_exon_variant          0
# 15          regulatory_region_variant          1
# 16          regulatory_region_variant          1
# 17          regulatory_region_variant          1
# 18          regulatory_region_variant          1
# 19          regulatory_region_variant          1

colnames(sample)
sample <- sample[,c("Symbol", "rs_id_dbSNP151_GRCh38p7", "variant_id", "maf", "DISEASE/TRAIT", "CONTEXT", "INTERGENIC", "gene_id", "MAPPED_GENE", "REPORTED GENE(S)")]
colnames(sample) <- c("Gene Symbol", "rs_id_dbSNP151_GRCh38p7", "GTEx variant ID", "MAF", "GWAS Trait", "Genomic Context", "Intergenic?", "Ensembl ID", "GWAS Mapping", "Reported Mapping")

sample %>% 
  filter(`Gene Symbol` == "FOLR3", !is.na(`GWAS Trait`)) %>% 
  head()
# Gene Symbol rs_id_dbSNP151_GRCh38p7            GTEx variant ID       MAF                       GWAS Trait
# 1       FOLR3             rs113841147 chr11_72357234_ACGTG_A_b38 0.0992537                           Height
# 2       FOLR3              rs12421674     chr11_72345547_C_T_b38 0.1246270          Systolic blood pressure
# 3       FOLR3                rs504217     chr11_72295042_C_T_b38 0.1171640         Diastolic blood pressure
# 4       FOLR3                rs504217     chr11_72295042_C_T_b38 0.1171640         Diastolic blood pressure
# 5       FOLR3              rs57191897     chr11_72224396_A_C_b38 0.0492537 FOLR2/HAVCR2 protein level ratio
# 6       FOLR3                rs652197     chr11_72138695_C_T_b38 0.1865670              Serum folate levels
# Genomic Context Intergenic?      Ensembl ID GWAS Mapping Reported Mapping
# 1                     intron_variant           0 ENSG00000110203         CLPB             CLPB
# 2                     intron_variant           0 ENSG00000110203         CLPB             <NA>
#   3 non_coding_transcript_exon_variant           0 ENSG00000110203         CLPB             CLPB
# 4 non_coding_transcript_exon_variant           0 ENSG00000110203         CLPB             <NA>
#   5                     intron_variant           0 ENSG00000110203       INPPL1             <NA>
#   6                     intron_variant           0 ENSG00000110203        FOLR3             <NA>

stable_study_counts <- stable # Loaded in from shiny-app/data.R
colnames(stable_study_counts)[3] <- "Gene Symbol"
colnames(stable_study_counts)[2] <- "Stable Study Count"
stable_study_counts <- stable_study_counts[,c(2,3)]
sample <- merge(sample, stable_study_counts, by = "Gene Symbol", all.x = T)

dynamic_study_counts <- dynamic
colnames(dynamic_study_counts)[3] <- "Gene Symbol"
colnames(dynamic_study_counts)[2] <- "Dynamic Study Count"
dynamic_study_counts <- dynamic_study_counts[,c(2,3)]
sample <- merge(sample, dynamic_study_counts, by = "Gene Symbol", all.x = T)

write_rds(sample, "./shiny-data-new/data/eqtls/gtex_variants_and_phenotypes_whole_blood.rds")
