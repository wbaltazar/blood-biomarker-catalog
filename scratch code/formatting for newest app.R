## Gomez
## Load in norm_expr from Gomez_variation.R
print(norm_expr)
old <- readRDS("~/Desktop/work_repo/Box organization/1results/shiny-data-new/data/expr/gomez.rds")
print(old)
old_names <- rownames(old)
library(EnsDb.Hsapiens.v86)
new_names <- mapIds(EnsDb.Hsapiens.v86, keys = rownames(norm_expr), column = "SYMBOL", keytype = "GENEID")
length(old_names)
length(new_names)
length(na.omit(new_names))
# [1] 15180
new_names <- na.omit(new_names)
anyDuplicated(new_names)
print(new_names[570])
# ENSG00000143248
# "RGS5"
rownames(old)[which(rownames(old) == "RGS5")]
## Rownames of old matrix is a named list?
rownames(norm_expr)
new_names <- mapIds(EnsDb.Hsapiens.v86, keys = rownames(norm_expr), column = "SYMBOL", keytype = "GENEID")
rownames(norm_expr) <- new_names
head(norm_expr)
## Write to new application folder
write_rds(norm_expr, file = "./Box organization/1results/RShiny-application/data/expr/gomez.rds")

## Meaburn
## Load in norm_expr files from Meaburn variationm DAY 1
old <- readRDS("./Box organization/1results/shiny-data-new/data/expr/meaburn.rds")
head(old)
## The rownames are the symbols: the probe is a separate column. The gene_graph_p function utilizes this.
head(norm_expr)
colnames(norm_expr) <- str_replace(colnames(norm_expr), pattern = "\\.CEL.gz", replacement = "")
probes <- unlist(rownames(norm_expr))
class(probes)
class(norm_expr)
norm_expr <- as.data.frame(norm_expr)
norm_expr$Probe <- probes
old_names <- rownames(old)
head(old_names, 100)
# [1] "DDR1"       "RFC2"       "HSPA6"      "PAX8"       "GUCA1A"     "UBA7"       "THRA"       "PTPN21"    
# [9] "CCL5"       "CYP2E1"     "EPHB3"      "ESRRA"      "CYP2A6"     "SCARB1"     "TTLL12"     "NA."       
# [17] "WFDC2"      "MAPK1"      "MAPK1.1"    "ADAM32"     "SPATA17"    "PRR22"      "PRR22.1"    "PXK"       
# [25] "PXK.1"      "VPS18"      "MSANTD3"    "SLC46A1"    "SLC46A1.1"  "TIMD4"      "SLC39A5"    "ZDHHC11"   
# [33] "ATP6V1E2"   "AFG3L1P"    "CILP2"      "CILP2.1"    "PIGX"       "TMEM196"    "SLC39A13"   "BEST4"     
# [41] "AK9"        "CORO6"      "TMEM106A"   "TMEM106A.1" "ALG10"      "ALG10.1"    "TTC39C"     "NEXN"      
# [49] "C15orf40"   "RAX2"       "MFAP3"      "EYA3"       "GIMAP1"     "GIMAP1.1"   "GIMAP1.2"   "KLK8"      
# [57] "CCDC65"     "CCDC65.1"   "PABIR3"     "PABIR3.1"   "CFAP53"     "CFAP53.1"   "ARMCX4"     "RBBP6"     
# [65] "CENPBD1P"   "TRIOBP"     "TRIOBP.1"   "CATSPER1"   "HOXD4"      "GSC"        "SP7"        "PDE7A"     
# [73] "CNOT7"      "CRYZL1"     "PRSS33"     "PRSS33.1"   "CBARP"      "CBARP.1"    "MCMDC2"     "TIRAP"     
# [81] "LEAP2"      "MSI2"       "SCIN"       "SCIN.1"     "CTCFL"      "C4orf33"    "C4orf33.1"  "C4orf33.2" 
# [89] "ZNF333"     "TVP23C"     "RDH10"      "RDH10.1"    "SRSF12"     "GARIN4"     "GARIN4.1"   "GAPT"      
# [97] "SCUBE1"     "ERICH5"     "ERICH5.1"   "CCDC185"  
length(grep("^NA", old_names))
# [1] 10214
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
grep("^\\d", new_names)
# many
new_names <- make.names(new_names, unique = T)
rownames(norm_expr) <- new_names
write_rds(norm_expr, file = "./Box organization/1results/RShiny-application/data/expr/Meaburn1.rds")

## DAY 2
colnames(norm_expr) <- str_replace(colnames(norm_expr), pattern = "\\.CEL.gz", replacement = "")
probes <- unlist(rownames(norm_expr))
class(probes)
class(norm_expr)
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
grep("^\\d", new_names)
new_names <- make.names(new_names, unique = T)
rownames(norm_expr) <- new_names
write_rds(norm_expr, file = "./Box organization/1results/RShiny-application/data/expr/Meaburn2.rds")

## Gosch
old <- readRDS("./Box organization/1results/shiny-data-new/data/expr/gosch.rds")
old[1:5,1:5]
#              GSM6965575 GSM6965576 GSM6965577 GSM6965578 GSM6965579
# DDX11L1      -3.4684370 -1.9183553 -1.8157963 -0.6944219  -1.618335
# WASH7P        3.2273713  2.8169700  2.5256720  2.9325824   4.150867
# MIR6859.1     0.4384536  0.9620631  0.6555094  1.9288106   1.957349
# LOC101927589  0.1560539  1.2515697  1.2005055  0.3286617   1.758041
# LOC729737     5.9571286  6.6274868  6.3699802  5.8326767   7.464019
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

## Obermoser
## C1
old <- readRDS("./Box organization/1results/shiny-data-new/data/expr/obermoser.rds")
old[1:5,1:5]
#           GSM744835 GSM744836 GSM744837 GSM744838 GSM744839
# ATP4A      4.481638  4.525867  4.510462  4.732031  4.571580
# RRN3       6.433098  6.350373  5.234039  7.160030  6.368260
# CDKN2C     4.973508  4.603080  4.753782  4.839487  4.606624
# LINC00273  5.184856  4.757480  5.822481  4.850829  5.259615
# DAPK3      5.365804  4.830513  5.307994  5.406890  4.854150
old$Probe[1:5]
# "ILMN_2209417" "ILMN_1765401" "ILMN_1652171" "ILMN_1652170" "ILMN_1792710"
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
write_rds(norm_expr, file = "./Box organization/1results/RShiny-application/data/expr/obermoser1.rds")

## C2
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
write_rds(norm_expr, file = "./Box organization/1results/RShiny-application/data/expr/obermoser2.rds")

## C3
head(colnames(norm_expr))
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
# [1] 27431   140 (same as Obermoser_variation.R, but one more column for Probe)
write_rds(norm_expr, file = "./Box organization/1results/RShiny-application/data/expr/obermoser3.rds")

## C4
head(colnames(norm_expr))
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
write_rds(norm_expr, file = "./Box organization/1results/RShiny-application/data/expr/obermoser4.rds")

## Dusek
old <- readRDS("./Box organization/1results/shiny-data-new/data/expr/bhasin.rds")
old[1:5,1:5]
# GSM253663 GSM253664 GSM253665 GSM253666 GSM253668
# DDR1    6.909879  6.944320  7.055973  7.381506  7.059420
# RFC2    6.023353  5.955863  6.220316  6.332936  5.713215
# HSPA6   8.432983  8.243785  8.843591  8.209773  8.324805
# PAX8    8.504396  8.805405  8.452344  8.546908  8.463327
# GUCA1A  3.345014  3.643303  3.358237  3.425642  3.520712
head(old$Probe)
# [1] "1007_s_at" "1053_at"   "117_at"    "121_at"    "1255_g_at" "1294_at"  
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
old <- readRDS("./Box organization/1results/shiny-data-new/data/expr/rusch.rds")
old[1:5,1:5]
head(old$Probe)
norm_expr[1:5,1:5]
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
write_rds(norm_expr, file = "./Box organization/1results/RShiny-application/data/expr/rusch.rds")

## LaRocca
old <- readRDS("./Box organization/1results/shiny-data-new/data/expr/larocca.rds")
old[1:5,1:5]
# GSM6255519 GSM6255520 GSM6255521 GSM6255522 GSM6255523
# A1BG      0.99875224  1.8045462  1.1688350  1.2111409  1.3580736
# A1BG-AS1  1.02004702  1.2703975  1.2096687  0.8344699  1.3028441
# A2M       0.02532196 -0.6907135 -0.5612442  0.5342963 -0.2268889
# A2M-AS1  -1.00662093 -1.4904148 -1.1122594 -0.8840163 -1.5186550
# AAAS      3.84689812  4.5054047  4.0567333  3.7215642  4.0548616
colnames(old)
dim(old)
# [1] 15553    58
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

## From Rshiny-application
class(e1) # Gosch
class(e2) # Meaburn
class(e4) # Gosch
class(e5) # Obermoser
class(e9) # Dusek
class(e11) # LaRocca
## Since all RNA-seq datasets are matrices, we can replace their rownames with named lists of ENSEMBL IDs
new_names <- mapIds(x = EnsDb.Hsapiens.v86, keys = str_replace(rownames(e11), "\\.", "-"), column = "GENEID", keytype = "SYMBOL")
names(rownames(e11)) <- unname(new_names)
head(rownames(e11))
# ENSG00000121410            <NA> ENSG00000175899            <NA> ENSG00000094914 ENSG00000081760 
# "A1BG"      "A1BG.AS1"           "A2M"       "A2M.AS1"          "AAAS"          "AACS" 

# ENSG00000121410 ENSG00000268895 ENSG00000175899 ENSG00000245105 ENSG00000094914 ENSG00000081760 
# "A1BG"      "A1BG.AS1"           "A2M"       "A2M.AS1"          "AAAS"          "AACS"
## replacing all . with -


new_names <- mapIds(x = EnsDb.Hsapiens.v86, keys = str_replace(rownames(e4), "\\.", "-"), column = "GENEID", keytype = "SYMBOL")
names(rownames(e4)) <- unname(new_names)
head(rownames(e4))
# ENSG00000223972 ENSG00000226210 ENSG00000278267            <NA>            <NA>            <NA> 
  # "DDX11L1"        "WASH7P"     "MIR6859.1"  "LOC101927589"     "LOC729737"  "LOC100996442"

write_rds(e4, "./data/expr/gosch.rds")
write_rds(e11, "./data/expr/larocca.rds")
