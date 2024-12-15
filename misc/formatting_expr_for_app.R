## Date: Jul 24 2024
# This code writes RDS files which can be read by the Rshiny-application app.R application and data.R source code.
# This code works after all other study_data has been generated.

sourcePartial <- function(fn, startTag='#from here', endTag='#to here') {
  lines <- scan(fn, what = character(), sep = "\n", quiet = TRUE)
  st <- grep(startTag, lines, fixed = T)
  en <- grep(endTag, lines, fixed = T)
  tc <- textConnection(lines[(st+1):(en-1)])
  source(tc)
  close(tc)
}

## Gomez ----
## Load in norm_expr from Gomez_variation.R
sourcePartial("~/Desktop/work_repo/github/Gomez_study_data/code/Gomez_variation.R",
              startTag = "## Load libraries ----", 
              endTag = "# Variation analysis ----")
library(EnsDb.Hsapiens.v86)
new_names <- mapIds(EnsDb.Hsapiens.v86, keys = rownames(norm_expr), column = "SYMBOL", keytype = "GENEID")
length(new_names)
# [1] 15416
length(na.omit(new_names))
# [1] 15180
rownames(norm_expr) <- new_names
head(norm_expr)
str(norm_expr)
# num [1:15416, 1:28] 4.52 3.32 3.38 4.46 2.05 ...
# - attr(*, "dimnames")=List of 2
# ..$ : Named chr [1:15416] "TNFRSF4" "TNFRSF18" "ATAD3B" "THAP3" ...
# .. ..- attr(*, "names")= chr [1:15416] "ENSG00000186827" "ENSG00000186891" "ENSG00000160072" "ENSG00000041988" ...
# ..$ : chr [1:28] "GSM7662572" "GSM7662571" "GSM7662575" "GSM7662576" ...
## Write to RShiny application folder
write_rds(norm_expr, file = "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/expr/gomez.rds")
## Write phenotype data to application folder
write.csv(pheno_data, "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/phenotypes/1Gomez_pheno.csv")

## Meaburn ----
## Load in norm_expr files from Meaburn variationm DAY 1
## Note that the rownames are the symbols: the probe is a separate column. The gene_graph_p function utilizes this.
## Load in norm_expr from Meaburn_variation.R, day 1
sourcePartial("~/Desktop/work_repo/github/Meaburn_study_data/code/Meaburn_variation.R",
              startTag = "## Load libraries ----", 
              endTag = "## varianceParitition 1 ----")
head(norm_expr)
dim(norm_expr)
# [1] 25419    19
colnames(norm_expr) <- str_replace(colnames(norm_expr), pattern = "\\.CEL.gz", replacement = "")
probes <- unlist(rownames(norm_expr))
class(probes)
# [1] "character"
class(norm_expr)
# [1] "matrix" "array" 
norm_expr <- as.data.frame(norm_expr)
norm_expr$Probe <- probes
new_names <- mapIds(x = hgu133plus2.db, keys = probes, column = "SYMBOL", keytype = "PROBEID")
head(unname(new_names), 100)
# [1] "DDR1"      "RFC2"      "HSPA6"     "PAX8"      "UBA7"      "THRA"      "CCL5"      "CYP2E1"    "ESRRA"     "CYP2A6"   
# [11] "SCARB1"    "TTLL12"    NA          "MAPK1"     "MAPK1"     "PXK"       "PXK"       "MSANTD3"   "SLC46A1"   "ZDHHC11"  
# [21] "ATP6V1E2"  "AFG3L1P"   "PIGX"      "CORO6"     "TMEM106A"  "TMEM106A"  "ALG10"     "TTC39C"    "NEXN"      "C15orf40" 
# [31] "MFAP3"     "GIMAP1"    "GIMAP1"    "GIMAP1"    "CCDC65"    "CCDC65"    "PABIR3"    "PABIR3"    "RBBP6"     "CENPBD1P" 
# [41] "TRIOBP"    "PDE7A"     "CNOT7"     "CRYZL1"    "PRSS33"    "PRSS33"    "TIRAP"     "LEAP2"     "MSI2"      "C4orf33"  
# [51] "ZNF333"    "GAPT"      "CLEC12A"   "CLEC4F"    "NEDD1"     "ETV3"      "TM2D3"     "ZNF485"    "ACAP2"     "PTPRC"    
# [61] "LACTB"     "LACTB"     "AFG1L"     "SLAMF6"    "GPBAR1"    "HIPK1"     "ACVR1C"    "TMEM74"    "VTI1A"     "TAGAP"    
# [71] "TAGAP"     "LETM2"     "CLEC4C"    "NLRC4"     "NLRC4"     "ZNF570"    "ABCC13"    "IL12RB1"   "FAM223B"   "JAK1"     
# [81] "JAK1"      "CDC42SE2"  "CDC42SE2"  "COP1"      "STX6"      "POLR2J2"   "PMS2P5"    "HSH2D"     "TRNT1"     "HERPUD2"  
# [91] "ARSG"      "ZNF101"    "ZNF101"    "PHC3"      "TNFRSF10A" "RFFL"      "HPS4"      "C5orf22"   "FLCN"      "SH2D3C" 
sum(is.na(new_names))
# [1] 3100
for (i in 1:length(new_names)) {
  if (is.na(new_names[i])) {
    new_names[i] <- names(new_names)[i]
  } else {
    next
  }
}
sum(is.na(new_names))
# [1] 0
length(grep("^\\d", new_names))
# [1] 3052
new_names <- make.names(new_names, unique = T)
rownames(norm_expr) <- new_names
write_rds(norm_expr, file = "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/expr/Meaburn1.rds")
## Write phenotype data to application folder
write.csv(p1, "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/phenotypes/2Meaburn1_pheno.csv")

## DAY 2
sourcePartial("~/Desktop/work_repo/github/Meaburn_study_data/code/Meaburn_variation.R",
              startTag = "# Analysis of day 2 data ----", 
              endTag = "## varianceParitition 2 ----")
colnames(norm_expr) <- str_replace(colnames(norm_expr), pattern = "\\.CEL.gz", replacement = "")
probes <- unlist(rownames(norm_expr))
norm_expr <- as.data.frame(norm_expr)
dim(norm_expr)
# [1] 26474    17
norm_expr$Probe <- probes
new_names <- mapIds(x = hgu133plus2.db, keys = probes, column = "SYMBOL", keytype = "PROBEID")
head(unname(new_names), 100)
for (i in 1:length(new_names)) {
  if (is.na(new_names[i])) {
    new_names[i] <- names(new_names)[i]
  } else {
    next
  }
}
sum(is.na(new_names))
# [1] 0
new_names <- make.names(new_names, unique = T)
rownames(norm_expr) <- new_names
write_rds(norm_expr, file = "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/expr/Meaburn2.rds")
## Write phenotype data to application folder
write.csv(p2, "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/phenotypes/2Meaburn2_pheno.csv")

## Gosch ----
sourcePartial("~/Desktop/work_repo/github/Gosch_study_data/code/Gosch_variation.R",
              startTag = "## Load libraries ----", 
              endTag = "# Variation analysis ----")
norm_expr[1:5,1:5]
# Blood, ParticipantA, 08.00h Blood, ParticipantB, 08.00h Blood, ParticipantC, 08.00h
# DDX11L1                       -3.4947730                   -1.790080                  -1.8193774
# WASH7P                         3.2010353                    2.945245                   2.5220909
# MIR6859.1                      0.4121176                    1.090338                   0.6519283
# LOC101927589                   0.1297179                    1.379845                   1.1969244
# LOC729737                      5.9307926                    6.755762                   6.3663991
# Blood, ParticipantD, 08.00h Blood, ParticipantE, 08.00h
# DDX11L1                       -0.8110214                   -1.598915
# WASH7P                         2.8159829                    4.170287
# MIR6859.1                      1.8122111                    1.976770
# LOC101927589                   0.2120622                    1.777461
# LOC729737                      5.7160772                    7.483439
identical(colnames(norm_expr), pheno_data$title)
# [1] TRUE
colnames(norm_expr) <- pheno_data$geo_accession
write_rds(norm_expr, file = "./Box organization/1results/RShiny-application/data/expr/gosch.rds")
## We added ENSG IDs to the row names at the bottom of the script
## Write phenotype data to application folder
write.csv(pheno_data, "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/phenotypes/3Gosch_pheno.csv")

## Obermoser ----
## C1
sourcePartial("~/Desktop/work_repo/github/Obermoser_study_data/code/Obermoser_variation.R",
              startTag = "# across cohorts measures different time intervals, including an hour-by-hour measurement.", 
              endTag = "## varianceParitition 1 ----")
norm_expr[1:4,1:4]
#              X4835507041_F.AVG_Signal X4835507052_A.AVG_Signal X4835507049_K.AVG_Signal X4853887047_C.AVG_Signal
# ILMN_1765401                 6.433098                 6.350373                 5.234039                 7.160030
# ILMN_1652170                 5.184856                 4.757480                 5.822481                 4.850829
# ILMN_1792710                 5.365804                 4.830513                 5.307994                 5.406890
# ILMN_1703558                 8.187752                 8.537027                 8.147067                 8.438531
identical(substr(p1v$array, start = 1, stop = 10), str_extract(colnames(norm_expr), "\\d+"))
# [1] TRUE
colnames(norm_expr) <- p1v$geo_accession
dim(norm_expr)
# [1] 21699   161
newnames <- mapIds(x = illuminaHumanv3.db, keys = rownames(norm_expr), column = "SYMBOL", keytype = "PROBEID")
for (i in 1:length(newnames)) {
  if (is.na(newnames[i])) {
    newnames[i] <- names(newnames)[i]
  } else {
    next
  }
}
sum(is.na(newnames))
# [1] 0
class(norm_expr)
# [1] "matrix" "array" 
norm_expr <- as.data.frame(norm_expr)
norm_expr$Probe <- rownames(norm_expr)
rownames(norm_expr) <- make.names(newnames, unique = T)
norm_expr[1:5, 1:5]
#           GSM744835 GSM744836 GSM744837 GSM744838 GSM744839
# RRN3       6.433098  6.350373  5.234039  7.160030  6.368260
# LINC00273  5.184856  4.757480  5.822481  4.850829  5.259615
# DAPK3      5.365804  4.830513  5.307994  5.406890  4.854150
# FHL3       8.187752  8.537027  8.147067  8.438531  8.051755
# DNAJC15    7.545133  6.856175  6.924462  7.075560  7.333289
head(norm_expr$Probe)
# [1] "ILMN_1765401" "ILMN_1652170" "ILMN_1792710" "ILMN_1703558" "ILMN_1812666" "ILMN_1720799"
write_rds(norm_expr, file = "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/expr/obermoser1.rds")
## Write phenotype data to application folder
write.csv(p1v, "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/phenotypes/4obermoser1_pheno.csv")

## C2
sourcePartial("~/Desktop/work_repo/github/Obermoser_study_data/code/Obermoser_variation.R",
              startTag = "# Analysis 2: Test Set Finger (p1f) ----", 
              endTag = "## varianceParitition 2 ----")
head(colnames(norm_expr))
# X5447320035_B.AVG_Signal" "X5447320035_G.AVG_Signal" "X5447320036_A.AVG_Signal" "X5447320036_B.AVG_Signal"
# [5] "X5447320036_C.AVG_Signal" "X5447320036_J.AVG_Signal"
identical(substr(p1f$array, start = 1, stop = 10), str_extract(colnames(norm_expr), "\\d+"))
# [1] TRUE
colnames(norm_expr) <- p1f$geo_accession
dim(norm_expr)
# [1] 16054    49
newnames <- mapIds(x = illuminaHumanv3.db, keys = rownames(norm_expr), column = "SYMBOL", keytype = "PROBEID")
for (i in 1:length(newnames)) {
  if (is.na(newnames[i])) {
    newnames[i] <- names(newnames)[i]
  } else {
    next
  }
}
sum(is.na(newnames))
# [1] 0
norm_expr <- as.data.frame(norm_expr)
norm_expr$Probe <- rownames(norm_expr)
rownames(norm_expr) <- make.names(newnames, unique = T)
norm_expr[1:5, 1:5]
head(norm_expr$Probe)
dim(norm_expr)
# [1] 16054    50 (same as Obermoser_variation.R, but one more column for Probe)
write_rds(norm_expr, file = "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/expr/obermoser2.rds")
## Write phenotype data to application folder
write.csv(p1f, "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/phenotypes/4obermoser2_pheno.csv")

## C3
sourcePartial("~/Desktop/work_repo/github/Obermoser_study_data/code/Obermoser_variation.R",
              startTag = "# Analysis 3: Test Set Vein (p2v) ----", 
              endTag = "## varianceParitition 3 ----")
head(colnames(norm_expr))
# [1] "X5356583035_D.AVG_Signal" "X5356583041_E.AVG_Signal" "X5356583031_A.AVG_Signal" "X5356583049_D.AVG_Signal" "X5356583040_D.AVG_Signal"
# [6] "X5356583038_G.AVG_Signal"
identical(substr(p2v$array, start = 1, stop = 10), str_extract(colnames(norm_expr), "\\d+"))
# [1] TRUE
colnames(norm_expr) <- p2v$geo_accession
dim(norm_expr)
# [1] 21167   144
newnames <- mapIds(x = illuminaHumanv3.db, keys = rownames(norm_expr), column = "SYMBOL", keytype = "PROBEID")
for (i in 1:length(newnames)) {
  if (is.na(newnames[i])) {
    newnames[i] <- names(newnames)[i]
  } else {
    next
  }
}
sum(is.na(newnames))
# [1] 0
norm_expr <- as.data.frame(norm_expr)
norm_expr$Probe <- rownames(norm_expr)
rownames(norm_expr) <- make.names(newnames, unique = T)
norm_expr[1:5, 1:5]
head(norm_expr$Probe)
dim(norm_expr)
# [1] 21167   144
write_rds(norm_expr, file = "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/expr/obermoser3.rds")
## Write phenotype data to application folder
write.csv(p2v, "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/phenotypes/4obermoser3_pheno.csv")

## C4
sourcePartial("~/Desktop/work_repo/github/Obermoser_study_data/code/Obermoser_variation.R",
              startTag = "# Analysis 4: Training Set Finger (p2f) ----", 
              endTag = "## varianceParitition 4 ----")
head(colnames(norm_expr))
# [1] "X5307141014_E.AVG_Signal" "X5307141024_E.AVG_Signal" "X5307141024_K.AVG_Signal" "X5307141024_F.AVG_Signal" "X5322244022_H.AVG_Signal"
# [6] "X5307141024_J.AVG_Signal"
identical(substr(p2f$array, start = 1, stop = 10), str_extract(colnames(norm_expr), "\\d+"))
# [1] TRUE
colnames(norm_expr) <- p2f$geo_accession
dim(norm_expr)
# [1] 19825   185
newnames <- mapIds(x = illuminaHumanv3.db, keys = rownames(norm_expr), column = "SYMBOL", keytype = "PROBEID")
for (i in 1:length(newnames)) {
  if (is.na(newnames[i])) {
    newnames[i] <- names(newnames)[i]
  } else {
    next
  }
}
sum(is.na(newnames))
# [1] 0
norm_expr <- as.data.frame(norm_expr)
norm_expr$Probe <- rownames(norm_expr)
rownames(norm_expr) <- make.names(newnames, unique = T)
norm_expr[1:5, 1:5]
head(norm_expr$Probe)
dim(norm_expr)
# [1] 19825   186
write_rds(norm_expr, file = "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/expr/obermoser4.rds")
## Write phenotype data to application folder
write.csv(p2f, "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/phenotypes/4obermoser4_pheno.csv")

## Dusek ----
sourcePartial("~/Desktop/work_repo/github/Dusek_study_data/code/Dusek_variation.R",
              startTag = "## Load libraries ----", 
              endTag = "## varianceParitition ----")
norm_expr[1:5,1:5]
#         GSM253663 GSM253664 GSM253665 GSM253666 GSM253668
# 1053_at  6.020718  5.967918  6.221690  6.340561  5.705007
# 117_at   8.405883  8.248298  8.831063  8.202793  8.309342
# 121_at   8.494362  8.781406  8.425188  8.529056  8.447904
# 1294_at  7.761764  7.398413  7.644221  7.626453  7.577665
# 1316_at  5.800415  6.072798  5.869718  5.698908  5.951719
dim(norm_expr)
# [1] 15963    44
newnames <- mapIds(x = hgu133plus2.db, keys = rownames(norm_expr), column = "SYMBOL", keytype = "PROBEID")
for (i in 1:length(newnames)) {
  if (is.na(newnames[i])) {
    newnames[i] <- names(newnames)[i]
  } else {
    next
  }
}
sum(is.na(newnames))
# [1] 0
norm_expr <- as.data.frame(norm_expr)
norm_expr$Probe <- rownames(norm_expr)
rownames(norm_expr) <- make.names(newnames, unique = T)
norm_expr[1:5, 1:5]
head(norm_expr$Probe)
write_rds(norm_expr, file = "./Box organization/1results/RShiny-application/data/expr/dusek.rds")
## Write phenotype data to application folder
write.csv(pheno_data, "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/phenotypes/5Dusek_pheno.csv")

## Rusch ----
sourcePartial("~/Desktop/work_repo/github/Rusch_study_data/code/Rusch_variation.R",
              startTag = "## Load libraries ----", 
              endTag = "## varianceParitition ----")
dim(norm_expr)
# [1] 21374    38
norm_expr[1:5,1:3]
# GSM2175170_md109_HG-U133_Plus_2.CEL.gz GSM2175171_md112_HG-U133_Plus_2.CEL.gz GSM2175182_md138_HG-U133_Plus_2.CEL.gz
# 1053_at                               6.599784                               6.712995                               6.859752
# 117_at                                9.007288                               8.805269                               8.375141
# 121_at                                8.227604                               8.331777                               8.171151
# 1294_at                               7.829571                               7.858019                               7.998548
# 1316_at                               5.646710                               5.949423                               6.005740
newnames <- mapIds(x = hgu133plus2.db, keys = rownames(norm_expr), column = "SYMBOL", keytype = "PROBEID")
for (i in 1:length(newnames)) {
  if (is.na(newnames[i])) {
    newnames[i] <- names(newnames)[i]
  } else {
    next
  }
}
sum(is.na(newnames))
# [1] 0
norm_expr <- as.data.frame(norm_expr)
norm_expr$Probe <- rownames(norm_expr)
rownames(norm_expr) <- make.names(newnames, unique = T)
norm_expr[1:5, 1:5]
#       GSM2175170_md109_HG-U133_Plus_2.CEL.gz GSM2175171_md112_HG-U133_Plus_2.CEL.gz GSM2175182_md138_HG-U133_Plus_2.CEL.gz
# RFC2                                6.599784                               6.712995                               6.859752
# HSPA6                               9.007288                               8.805269                               8.375141
# PAX8                                8.227604                               8.331777                               8.171151
# UBA7                                7.829571                               7.858019                               7.998548
# THRA                                5.646710                               5.949423                               6.005740
#       GSM2175183_md139_HG-U133_Plus_2.CEL.gz GSM2175185_md141_HG-U133_Plus_2.CEL.gz
# RFC2                                7.145761                               6.817535
# HSPA6                               9.490100                               9.197723
# PAX8                                8.175178                               8.197186
# UBA7                                7.650091                               8.150917
# THRA                                5.540203                               5.482492
colnames(norm_expr) <- str_extract(colnames(norm_expr), pattern = "^GSM\\d+")
colnames(norm_expr)[39] <- "Probe"
head(norm_expr$Probe)
# [1] "1053_at"   "117_at"    "121_at"    "1294_at"   "1316_at"   "1405_i_at"
write_rds(norm_expr, file = "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/expr/rusch.rds")
## Write phenotype data to application folder
write.csv(pheno_data, "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/phenotypes/6Rusch_pheno.csv")

## LaRocca ----
sourcePartial("~/Desktop/work_repo/github/LaRocca_study_data/code/LaRocca_variation.R",
              startTag = "## Load libraries ----", 
              endTag = "## varianceParitition ----")
norm_expr[1:5,1:5]
# GSM6255519 GSM6255520 GSM6255521 GSM6255522 GSM6255523
# A1BG      0.98359095  1.7893584  1.1537217  1.1959563  1.3428967
# A1BG.AS1  1.00488573  1.2552097  1.1945554  0.8192853  1.2876672
# A2M       0.01016067 -0.7059013 -0.5763575  0.5191117 -0.2420658
# A2M.AS1  -1.02178222 -1.5056026 -1.1273727 -0.8992009 -1.5338319
# AAAS      3.83173682  4.4902168  4.0416200  3.7063796  4.0396847
dim(norm_expr)
# [1] 15555    59
write_rds(norm_expr, file = "./Box organization/1results/RShiny-application/data/expr/larocca.rds")

## Since all RNA-seq data sets are matrices, we can replace their row names with named lists that include ENSEMBL IDs
e11 <- read_rds("~/Desktop/work_repo/Box organization/1results/RShiny-application/data/expr/larocca.rds")
library(EnsDb.Hsapiens.v86)
new_names <- mapIds(x = EnsDb.Hsapiens.v86, keys = str_replace(rownames(e11), "\\.", "-"), column = "GENEID", keytype = "SYMBOL")
names(rownames(e11)) <- unname(new_names)
head(rownames(e11))
# ENSG00000121410 ENSG00000268895 ENSG00000175899 ENSG00000245105 ENSG00000094914 ENSG00000081760 
# "A1BG"      "A1BG.AS1"           "A2M"       "A2M.AS1"          "AAAS"          "AACS"
write_rds(e11, "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/expr/larocca.rds")

e4 <- read_rds("~/Desktop/work_repo/Box organization/1results/RShiny-application/data/expr/gosch.rds")
new_names <- mapIds(x = EnsDb.Hsapiens.v86, keys = str_replace(rownames(e4), "\\.", "-"), column = "GENEID", keytype = "SYMBOL")
names(rownames(e4)) <- unname(new_names)
head(rownames(e4))
# ENSG00000223972 ENSG00000226210 ENSG00000278267            <NA>            <NA>            <NA> 
  # "DDX11L1"        "WASH7P"     "MIR6859.1"  "LOC101927589"     "LOC729737"  "LOC100996442"
write_rds(e4, "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/expr/gosch.rds")
## Write phenotype data to application folder
write.csv(pheno_data, "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/phenotypes/7LaRocca_pheno.csv")

## Karlovich ----
## Load in norm_expr files from Karlovich variation batch 1
## Note that the rownames are the symbols: the probe is a separate column. The gene_graph_p function utilizes this.
## Load in norm_expr from Meaburn_variation.R, day 1
sourcePartial("~/Desktop/work_repo/github/Karlovich_study_data/code/Karlovich_variation.R",
              startTag = "## Load libraries ----", 
              endTag = "## varianceParitition 1 ----")
head(norm_expr)
dim(norm_expr)
# [1] 17635    66
colnames(norm_expr) <- str_replace(colnames(norm_expr), pattern = "\\.CEL.gz", replacement = "")
probes <- unlist(rownames(norm_expr))
class(probes)
# [1] "character"
class(norm_expr)
# [1] "matrix" "array" 
norm_expr <- as.data.frame(norm_expr)
norm_expr$Probe <- probes
new_names <- mapIds(x = hgu133plus2.db, keys = probes, column = "SYMBOL", keytype = "PROBEID")
head(unname(new_names), 100)
# [1] "RFC2"       "HSPA6"      "PAX8"       "UBA7"       "THRA"       "CCL5"       "ESRRA"      "SCARB1"     "TTLL12"    
# [10] "MAPK1"      "MAPK1"      "PXK"        "SLC46A1"    "ZDHHC11"    "AFG3L1P"    "PIGX"       "CORO6"      "TMEM106A"  
# [19] "NEXN"       "C15orf40"   "GIMAP1"     "GIMAP1"     "GIMAP1"     "CCDC65"     "RBBP6"      "CENPBD1P"   "PDE7A"     
# [28] "CRYZL1"     "PRSS33"     "LEAP2"      "MSI2"       "ZNF333"     "GAPT"       "CLEC12A"    "ETV3"       "TM2D3"     
# [37] "ZNF485"     "SERPINB11"  "ACAP2"      "PTPRC"      "LACTB"      "SLAMF6"     "GPBAR1"     "ACVR1C"     "TMEM74"    
# [46] "VTI1A"      "CLEC4C"     "NLRC4"      "RAI1"       NA           "IL12RB1"    "FAM223B"    "JAK1"       "JAK1"      
# [55] "CDC42SE2"   "CDC42SE2"   "COP1"       "POLR2J2"    "PMS2P5"     "HSH2D"      "HERPUD2"    "ARSG"       "ZNF101"    
# [64] "ZNF101"     "PHC3"       "HPS4"       "FLCN"       "SH2D3C"     "PPP1R3B"    "UCN3"       "ARL11"      "CARD16"    
# [73] "CARD16"     "SLC4A1"     "GOLM2"      "RHOXF1"     "KLHDC1"     "CADM2"      "ZNF563"     "CLEC4D"     "CLEC4D"    
# [82] "HELB"       "SEC62"      "SEC62"      "TLR4"       "SIGLEC10"   NA           "CACNG6"     "EXOC3L2"    "LINC01312" 
# [91] "MMEL1"      "TRIM67-AS1" "ZMYM6"      "DNAJB13"    "CNPY3"      "SCAMP1"     "SDCCAG8"    "CD300LF"    "PIP4K2B"   
# [100] "SLFN5"     
sum(is.na(new_names))
# [1] 1913
for (i in 1:length(new_names)) {
  if (is.na(new_names[i])) {
    new_names[i] <- names(new_names)[i]
  } else {
    next
  }
}
sum(is.na(new_names))
# [1] 0
length(grep("^\\d", new_names))
# [1] 1907
new_names <- make.names(new_names, unique = T)
rownames(norm_expr) <- new_names
write_rds(norm_expr, file = "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/expr/Karlovich1.rds")
## Write phenotype data to application folder
write.csv(p1, "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/phenotypes/8Karlovich1_pheno.csv")

## DAY 2
sourcePartial("~/Desktop/work_repo/github/Karlovich_study_data/code/Karlovich_variation.R",
              startTag = "# Analysis of batch 2 data ----", 
              endTag = "## varianceParitition 2 ----")
colnames(norm_expr) <- str_replace(colnames(norm_expr), pattern = "\\.CEL.gz", replacement = "")
dim(norm_expr)
# [1] 16352    40
probes <- unlist(rownames(norm_expr))
norm_expr <- as.data.frame(norm_expr)
norm_expr$Probe <- probes
new_names <- mapIds(x = hgu133plus2.db, keys = probes, column = "SYMBOL", keytype = "PROBEID")
head(unname(new_names), 100)
for (i in 1:length(new_names)) {
  if (is.na(new_names[i])) {
    new_names[i] <- names(new_names)[i]
  } else {
    next
  }
}
sum(is.na(new_names))
# [1] 0
new_names <- make.names(new_names, unique = T)
rownames(norm_expr) <- new_names
write_rds(norm_expr, file = "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/expr/Karlovich2.rds")
## Write phenotype data to application folder
write.csv(p2, "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/phenotypes/8Karlovich2_pheno.csv")