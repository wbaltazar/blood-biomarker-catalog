## June 28 2024
# I want to do two things: 1., confirm what the common symbols are between all 7 datasets, and 2., see if
# using the Spearman correlation coefficient significantly affects enrichment results.

# First, load in all the stability data for our genes.
setwd("~/Desktop/work_repo/Box organization/")
source("./type 2 gene analysis/06062024 analysis/new_finding_type_2.R")

# Also, run the code in eQTL enrichment, lines 18-34.
# The first time I did this, I got 9511 common symbols. Now, I only get 9474. What genes did I lose?
old_symbols <- read.csv("./eQTL enrichment analysis/median_stability_for_ubiquitous_genes.csv")
common_symbol_statistics # this is our current
missing_symbols <- setdiff(old_symbols$Symbol, common_symbol_statistics$Symbol)
# [1] "C1orf112"  "C12orf4"   "HSPB11"    "C20orf27"  "TTC26"     "SLC9A3R1"  "PRPF4B"    "FAM172A"   "C3orf14"  
# [10] "C1orf109"  "FAM104A"   "C12orf29"  "MRPS36"    "C11orf1"   "C17orf80"  "RNF165"    "SPATA5"    "SOGA1"    
# [19] "C18orf25"  "PPP1R32"   "KIAA1522"  "C9orf64"   "C19orf48"  "SPATA5L1"  "CCDC96"    "C11orf80"  "C16orf91" 
# [28] "ODF3B"     "FAM104B"   "KIAA2026"  "C22orf46"  "METTL7A"   "C19orf54"  "HHLA3"     "TP73-AS1"  "C4orf48"  
# [37] "LINC00672"
## Is it in affymetrix? illumina?
## If they don't run in the app, we cut them. Testing now... Yup, they don't work. So we cut them.

## Okay, now we go to func_enrich.R where we did the correlation stuff.
### We used ALL the data (28,000 symbols) and Pearson correlation
combined_data <- type2[,c("Symbols", "Score")]
names(combined_data)[1] <- "Symbol"
combined_data <- merge(combined_data, deg[,c("Symbol", "P.value.Studies")], by = "Symbol", sort = F)
head(combined_data)
#     Symbol Score P.value.Studies
# 1    ERAP2    22               2
# 2 HLA-DQB1    21               1
# 3   ARL17A    21               4
# 4  FAM118A    20               0
# 5     NAAA    20               2
# 6    GSTM1    19               1
combined_data$`Additive Correlation` <- 0
for (i in 1:nrow(combined_data)) {
  combined_data$`Additive Correlation`[i] <- cor(combined_data$Score[1:i], combined_data$P.value.Studies[1:i])
}
min_cor <- min(combined_data[-1,"Additive Correlation"])
print(min_cor)
# [1] -0.2659681
which.min(combined_data[-1,"Additive Correlation"])
# [1] 560, plus 1

## NEW METHOD: SPEARMAN ----
combined_data$`Spearman` <- 0
for (i in 1:nrow(combined_data)) {
  combined_data$Spearman[i] <- cor(combined_data$Score[1:i], combined_data$P.value.Studies[1:i], method = "spearman")
}
min_cor <- min(combined_data[-1,"Spearman"])
print(min_cor)
# [1] -0.3685968
which.min(combined_data[-1,"Spearman"])
# [1] 80, plus 1
## Entire distribution plot
ggplot(combined_data[-1,], aes(x = 2:nrow(combined_data), y = `Spearman`, color = `Spearman`)) +
  geom_line() + 
  scale_color_gradient(name = NULL, low = "blue", high = "red") +
  theme_minimal() + 
  xlab("Number of datapoints used") +
  geom_hline(yintercept = min_cor, color = "darkblue", linetype = "dashed")
ggplot(combined_data[2:500,], aes(x = 2:500, y = `Spearman`, color = `Spearman`)) +
  geom_line() + 
  scale_color_gradient(name = NULL, low = "blue", high = "red") +
  theme_minimal() + 
  xlab("Number of datapoints used") +
  geom_hline(yintercept = min_cor, color = "darkblue", linetype = "dashed") +
  geom_text(label = ifelse(combined_data[2:500,]$`Spearman` != min_cor, "", paste("Minimum correlation:", round(min_cor,2))), 
            nudge_y = -0.005, color = "darkblue")
ggplot(combined_data[2:10000,], aes(x = 2:10000, y = `Spearman`, color = `Spearman`)) +
  geom_line() + 
  scale_color_gradient(name = NULL, low = "blue", high = "red") +
  theme_minimal() + 
  xlab("Number of datapoints used") +
  geom_hline(yintercept = min_cor, color = "darkblue", linetype = "dashed") +
  geom_text(label = ifelse(combined_data[2:10000,]$`Spearman` != min_cor, "", paste("Minimum correlation:", round(min_cor,2))), 
            nudge_y = -0.005, color = "darkblue")


# ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= =======
## NEW METHOD: FILTERED, SPEARMAN ----
common_symbols <- read.table("./Box organization/1results/Symbols_common_to_all_datasets.txt")
common_symbols <- common_symbols$x
combined_data_filt <- combined_data[combined_data$Symbol %in% common_symbols, ]
dim(combined_data_filt)
# [1] 9474    5
## Recalculate both Pearson and Spearman correlations. What's different?
combined_data_filt$`Additive Correlation` <- 0
for (i in 1:nrow(combined_data_filt)) {
  combined_data_filt$`Additive Correlation`[i] <- cor(combined_data_filt$Score[1:i], combined_data_filt$P.value.Studies[1:i])
}
min_cor <- min(combined_data_filt[-1,"Additive Correlation"])
print(min_cor)
# [1] -0.2787637
which.min(combined_data_filt[-1,"Additive Correlation"])
# [1] 177, plus 1
combined_data_filt$`Spearman` <- 0
for (i in 1:nrow(combined_data_filt)) {
  combined_data_filt$Spearman[i] <- cor(combined_data_filt$Score[1:i], combined_data_filt$P.value.Studies[1:i], method = "spearman")
}
min_cor <- min(combined_data_filt[-1,"Spearman"])
print(min_cor)
# [1] -0.3685397
which.min(combined_data_filt[-1,"Spearman"])
# [1] 57
## Spearman distribution
ggplot(combined_data_filt[-1,], aes(x = 2:nrow(combined_data_filt), y = `Spearman`, color = `Spearman`)) +
  geom_line() + 
  scale_color_gradient(name = NULL, low = "blue", high = "red") +
  theme_minimal() + 
  xlab("Number of datapoints used") +
  geom_hline(yintercept = min_cor, color = "darkblue", linetype = "dashed")
ggplot(combined_data_filt[-1,], aes(x = 2:nrow(combined_data_filt), y = `Additive Correlation`, color = `Additive Correlation`)) +
  geom_line() + 
  scale_color_gradient(name = NULL, low = "blue", high = "red") +
  theme_minimal() + 
  xlab("Number of datapoints used") +
  geom_hline(yintercept = min_cor, color = "darkblue", linetype = "dashed")





# ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= =======
## TRYING NEW METHOD on  FLEXIBLE GENES
deg <- deg[order(deg$P.value.Studies, decreasing = T),]
combined_data2 <- deg[,c("Symbol","P.value.Studies")]
names(type2)[2] <- "Symbol"
combined_data2 <- merge(combined_data2, type2[,c("Symbol","Score")], sort = F)
head(combined_data2)
#      Symbol P.value.Studies Score
# 1     USP49               6     5
# 2    GTF3C4               6     4
# 3     VMA21               6     2
# 4     TTC38               6     7
# 5     LRIG1               5     4
# 6 LINC00921               5     1
combined_data2$`Additive Correlation` <- 0
for (i in 1:nrow(combined_data2)) {
  combined_data2$`Additive Correlation`[i] <- cor(combined_data2$Score[1:i], combined_data2$P.value.Studies[1:i])
}
min_cor <- min(na.omit(combined_data2$`Additive Correlation`))
print(min_cor)
# [1] -0.01975839
length(na.omit(combined_data2$`Additive Correlation`))
# [1] 28818, 4 less
which.min(na.omit(combined_data2$`Additive Correlation`))
# [1] [1] 2143, plus 4 because the top genes have the same score
combined_data2$`Spearman` <- 0
for (i in 1:nrow(combined_data2)) {
  combined_data2$`Spearman`[i] <- cor(combined_data2$Score[1:i], combined_data2$P.value.Studies[1:i], method = "spearman")
}
min_cor <- min(na.omit(combined_data2$`Spearman`))
print(min_cor)
# [1] -0.0173573
length(na.omit(combined_data2$`Spearman`))
# [1] 28818, 4 less
which.min(na.omit(combined_data2$`Spearman`))
# [1] 2143, crazy how it's the exact same!
ggplot(na.omit(combined_data2), aes(x = 5:nrow(combined_data2), y = `Spearman`, color = `Spearman`)) +
  geom_line() + 
  scale_color_gradient(name = NULL, low = "blue", high = "red") +
  theme_minimal() + 
  xlab("Number of datapoints used") +
  geom_hline(yintercept = min_cor, color = "darkblue", linetype = "dashed")
