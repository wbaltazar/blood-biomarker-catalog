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

## Gomez
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

## Meaburn
## Load in norm_expr files from Meaburn variationm DAY 1
## Note that the rownames are the symbols: the probe is a separate column. The gene_graph_p function utilizes this.
## Load in norm_expr from Meaburn_variation.R, day 1
sourcePartial("~/Desktop/work_repo/github/Meaburn_study_data/code/Meaburn_variation.R",
              startTag = "## Load libraries ----", 
              endTag = "## varianceParitition 1 ----")
head(norm_expr)
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
# [1] "DDR1"     "RFC2"     "HSPA6"    "PAX8"     "GUCA1A"   "UBA7"     "THRA"     "PTPN21"   "CCL5"    
# [10] "CYP2E1"   "EPHB3"    "ESRRA"    "CYP2A6"   "SCARB1"   "TTLL12"   NA         "WFDC2"    "MAPK1"   
# [19] "MAPK1"    "ADAM32"   "SPATA17"  "PRR22"    "PRR22"    "PXK"      "PXK"      "VPS18"    "MSANTD3" 
# [28] "SLC46A1"  "SLC46A1"  "TIMD4"    "SLC39A5"  "ZDHHC11"  "ATP6V1E2" "AFG3L1P"  "CILP2"    "CILP2"   
# [37] "PIGX"     "TMEM196"  "SLC39A13" "BEST4"    "AK9"      "CORO6"    "TMEM106A" "TMEM106A" "ALG10"   
# [46] "ALG10"    "TTC39C"   "NEXN"     "C15orf40" "RAX2"     "MFAP3"    "EYA3"     "GIMAP1"   "GIMAP1"  
# [55] "GIMAP1"   "KLK8"     "CCDC65"   "CCDC65"   "PABIR3"   "PABIR3"   "CFAP53"   "CFAP53"   "ARMCX4"  
# [64] "RBBP6"    "CENPBD1P" "TRIOBP"   "TRIOBP"   "CATSPER1" "HOXD4"    "GSC"      "SP7"      "PDE7A"   
# [73] "CNOT7"    "CRYZL1"   "PRSS33"   "PRSS33"   "CBARP"    "CBARP"    "MCMDC2"   "TIRAP"    "LEAP2"   
# [82] "MSI2"     "SCIN"     "SCIN"     "CTCFL"    "C4orf33"  "C4orf33"  "C4orf33"  "ZNF333"   "TVP23C"  
# [91] "RDH10"    "RDH10"    "SRSF12"   "GARIN4"   "GARIN4"   "GAPT"     "SCUBE1"   "ERICH5"   "ERICH5"  
# [100] "CCDC185" 
sum(is.na(new_names))
# [1] 10025
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
# [1] 9973
new_names <- make.names(new_names, unique = T)
rownames(norm_expr) <- new_names
write_rds(norm_expr, file = "./Box organization/1results/RShiny-application/data/expr/Meaburn1.rds")

## DAY 2
sourcePartial("~/Desktop/work_repo/github/Meaburn_study_data/code/Meaburn_variation.R",
              startTag = "# Analysis of day 2 data ----", 
              endTag = "## varianceParitition 2 ----")
colnames(norm_expr) <- str_replace(colnames(norm_expr), pattern = "\\.CEL.gz", replacement = "")
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
write_rds(norm_expr, file = "./Box organization/1results/RShiny-application/data/expr/Meaburn2.rds")

## Gosch
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

## Obermoser
## C1
sourcePartial("~/Desktop/work_repo/github/Obermoser_study_data/code/Obermoser_variation.R",
              startTag = "# across cohorts measures different time intervals, including an hour-by-hour measurement.", 
              endTag = "## varianceParitition 1 ----")
norm_expr[1:5,1:5]
# X4835507041_F.AVG_Signal X4835507052_A.AVG_Signal X4835507049_K.AVG_Signal X4853887047_C.AVG_Signal
# ILMN_2209417                 4.480964                 4.525181                 4.509781                 4.731090
# ILMN_1765401                 6.427961                 6.345370                 5.231324                 7.154682
# ILMN_1652171                 4.971882                 4.602338                 4.752796                 4.838257
# ILMN_1652170                 5.182336                 4.756486                 5.818189                 4.849556
# ILMN_1792710                 5.362644                 4.829314                 5.304991                 5.403610
# X4835507042_G.AVG_Signal
# ILMN_2209417                 4.570871
# ILMN_1765401                 6.363239
# ILMN_1652171                 4.605877
# ILMN_1652170                 5.256817
# ILMN_1792710                 4.852865
identical(substr(p1v$array, start = 1, stop = 10), str_extract(colnames(norm_expr), "\\d+"))
# [1] TRUE
colnames(norm_expr) <- p1v$geo_accession
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
# GSM744835 GSM744836 GSM744837 GSM744838 GSM744839
# ATP4A      4.480964  4.525181  4.509781  4.731090  4.570871
# RRN3       6.427961  6.345370  5.231324  7.154682  6.363239
# CDKN2C     4.971882  4.602338  4.752796  4.838257  4.605877
# LINC00273  5.182336  4.756486  5.818189  4.849556  5.256817
# DAPK3      5.362644  4.829314  5.304991  5.403610  4.852865
head(norm_expr$Probe)
# [1] "ILMN_2209417" "ILMN_1765401" "ILMN_1652171" "ILMN_1652170" "ILMN_1792710" "ILMN_1703558"
write_rds(norm_expr, file = "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/expr/obermoser1.rds")

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
# [1] 22413    50 (same as Obermoser_variation.R, but one more column for Probe)
write_rds(norm_expr, file = "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/expr/obermoser2.rds")

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
# [1] 27431   144
write_rds(norm_expr, file = "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/expr/obermoser3.rds")

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
# [1] 23484   186 (same as Obermoser_variation.R, but one more column for Probe)
write_rds(norm_expr, file = "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/expr/obermoser4.rds")

## Dusek
sourcePartial("~/Desktop/work_repo/github/Dusek_study_data/code/Dusek_variation.R",
              startTag = "## Load libraries ----", 
              endTag = "## varianceParitition ----")
norm_expr[1:5,1:5]
# GSM253663 GSM253664 GSM253665 GSM253666 GSM253668
# 1007_s_at  6.893945  6.933046  7.038875  7.360868  7.033034
# 1053_at    6.020718  5.967918  6.221690  6.340561  5.705007
# 117_at     8.405883  8.248298  8.831063  8.202793  8.309342
# 121_at     8.494362  8.781406  8.425188  8.529056  8.447904
# 1255_g_at  3.325685  3.624963  3.340179  3.407576  3.516416
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

## Rusch
sourcePartial("~/Desktop/work_repo/github/Rusch_study_data/code/Rusch_variation.R",
              startTag = "## Load libraries ----", 
              endTag = "## varianceParitition ----")
dim(norm_expr)
# [1] 54675    38
norm_expr[1:5,1:5]
# GSM2175170 GSM2175171 GSM2175182 GSM2175183 GSM2175185
# 1007_s_at   6.723975   6.796586   7.164297   6.610909   6.818922
# 1053_at     6.599784   6.712995   6.859752   7.145761   6.817535
# 117_at      9.007288   8.805269   8.375141   9.490100   9.197723
# 121_at      8.227604   8.331777   8.171151   8.175178   8.197186
# 1255_g_at   3.011893   2.993132   2.910271   2.848268   2.735265
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
# GSM2175170 GSM2175171 GSM2175182 GSM2175183 GSM2175185
# DDR1     6.723975   6.796586   7.164297   6.610909   6.818922
# RFC2     6.599784   6.712995   6.859752   7.145761   6.817535
# HSPA6    9.007288   8.805269   8.375141   9.490100   9.197723
# PAX8     8.227604   8.331777   8.171151   8.175178   8.197186
# GUCA1A   3.011893   2.993132   2.910271   2.848268   2.735265
head(norm_expr$Probe)
# [1] "1007_s_at" "1053_at"   "117_at"    "121_at"    "1255_g_at" "1294_at"  
write_rds(norm_expr, file = "~/Desktop/work_repo/Box organization/1results/RShiny-application/data/expr/rusch.rds")

## LaRocca
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