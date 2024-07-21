## June 24. Context: ran up to line 25 in the `formatting for newest app.R` script

## How many unique symbols are there before and after running str_replace?
length(unique(new_names))
new_names2 <- new_names %>% str_replace_all("\\.(?!\\d$)", "-") %>% str_remove_all(pattern = "\\..")
length(new_names2)
# [1] 15416
print(unname(new_names[grep("^RP", new_names)]))
length(unname(new_names[grep("^RP", new_names)]))
# [1] 1283
length(unique(unname(new_names[grep("^RP", new_names)])))
# [1] 1283
print(unname(new_names2[grep("^RP", new_names2)]))
length(unname(new_names2[grep("^RP", new_names2)]))
# [1] 1283
length(unique(unname(new_names2[grep("^RP", new_names2)])))
# [1] 1157
## We lose some names. Why is this?
anyDuplicated(new_names2[grep("^RP", new_names2)])
# [1] 4
new_names2[grep("^RP", new_names2)][4]
# "RP3-395M20"
new_names2[grep("RP3-395M20", new_names2)]
# [1] "RP3-395M20"    "RP3-395M20"    "RP3-395M20-12"
new_names[grep("RP3-395M20", new_names)]
# "RP3-395M20.9"  "RP3-395M20.8" "RP3-395M20.12" 
## In UCSC genome browser
# chr1:2,581,560-2,584,533 for RP3-395M20.9. A lncRNA with poor TSL (transcript support level)
# RP3-395M20.8 returns no results
# chr1:2,530,064-2,547,460 for RP3-395M20.12. About 50 kb upstream. TSL is NA (see more: https://genome.ucsc.edu/cgi-bin/hgc?hgsid=2303384264_auJteTAKMwTMROzlnxaA624VEv6i&db=hg38&c=chr1&l=2530063&r=2547460&o=2530063&t=2547460&g=wgEncodeGencodeBasicV38&i=ENST00000644246.1)
## Probably will be okay if these genes are annotated strangely. Are there any others we should be concerned about?
new_names2[grep("^RP", new_names2)][duplicated(new_names2[grep("^RP", new_names2)])]
# [1] "RP3-395M20"   "RP11-443B7"   "RP11-134G8"   "RP11-345P4"   "RP1-283E3"    "RP11-206L10"  "RP11-430C7"  
# [8] "RP11-373D23"  "RP11-107E5"   "RP11-222K16"  "RP11-415F23"  "RP11-154H23"  "RP11-299J3"   "RP11-185E8"  
# [15] "RP11-539L10"  "RP11-127B20"  "RP11-673E1"   "RP11-290F5"   "RP11-326I11"  "RP11-83A24"   "RP11-109G23" 
# [22] "RP11-10L12"   "RP11-1023L17" "RP11-321E2"   "RP11-423H2"   "RP11-1277A3"  "RP11-532F6"   "RP5-1112D6"  
# [29] "RP1-257A7"    "RP11-138A9"   "RP11-1275H24" "RP11-274B21"  "RP11-274B21"  "RP11-875O11"  "RP11-875O11" 
# [36] "RP11-110G21"  "RP11-473O4"   "RP5-855D21"   "RP5-855D21"   "RP11-87H9"    "RP11-531H8"   "RP11-677M14" 
# [43] "RP1-68D18"    "RP11-351I24"  "RP11-111M22"  "RP11-563J2"   "RP11-45A17"   "RP11-380G5"   "RP11-705C15" 
# [50] "RP11-996F15"  "RP11-1105G2"  "RP11-256L6"   "RP11-467L13"  "RP1-197B17"   "RP11-493L12"  "RP11-631N16" 
# [57] "RP11-338E21"  "RP13-554M15"  "RP11-278C7"   "RP11-327P2"   "RP11-596C23"  "RP11-47I22"   "RP11-930O11" 
# [64] "RP11-930O11"  "RP11-7F17"    "RP11-300J18"  "RP11-7F17"    "RP5-1021I20"  "RP11-73M18"   "RP11-76E17"  
# [71] "RP11-81A1"    "RP11-139H15"  "RP11-16E12"   "RP11-617F23"  "RP11-161M6"   "RP11-1102P22" "RP11-166B2"  
# [78] "RP11-266L9"   "RP11-1102P22" "RP11-266L9"   "RP11-166B2"   "RP11-830F9"   "RP11-196G11"  "RP11-196G11" 
# [85] "RP11-264B17"  "RP11-426C22"  "RP11-391L3"   "RP11-96D1"    "RP11-452L6"   "RP11-20I23"   "RP11-44F14"  
# [92] "RP11-44F14"   "RP11-44F14"   "RP11-44F14"   "RP11-44F14"   "RP11-1094M14" "RP11-82O19"   "RP5-890E16"  
# [99] "RP5-867C24"   "RP11-927P21"  "RP11-334C17"  "RP11-360N9"   "RP11-855A2"   "RP11-1055B8"  "RP11-45M22"  
# [106] "RP11-45M22"   "RP11-214O1"   "RP11-670E13"  "RP11-159D12"  "RP13-104F24"  "RP11-68I3"    "RP11-848P1"  
# [113] "RP11-888D10"  "RP11-879F14"  "RP11-158H5"   "RP11-405M12"  "RP11-4B16"    "RP11-162A12"  "RP11-35G9"   
# [120] "RP11-53B2"    "RP11-737O24"  "RP4-673D20"   "RP11-678G14"  "RP11-678G14"  "RP6-109B7"    "RP4-671O14"
return <- new_names2[grep("^RP", new_names2)][duplicated(new_names2[grep("^RP", new_names2)])]
poopoogenes <- names(new_names[grep(paste(return, collapse = "|"), new_names)])
write.table(poopoogenes, "~/Desktop/for_biomart.txt")
## Used web-version of biomart
poopooinfo <- read.csv("~/Downloads/mart_export (6).txt")
head(poopooinfo)
length(unique(poopooinfo$Gene.stable.ID))
# [1] 240
length(poopoogenes)
# [1] 244
unique(poopooinfo$HGNC.symbol)
# [1] ""                          "GVQW3"                     "LINC02363"                
# [4] "LINC02656"                 "RPL11P3"                   "GUSBP18"                  
# [7] "LINC02649"                 "ARHGAP27P1-BPTFP1-KPNA2P3" "MARK2P18"                 
# [10] "RPS10P7"                   "IFITM3P2"                  "CSRP1-AS1"                
# [13] "PCMTD1-DT"                 "NPIPB2"                    "LINC02939"                
# [16] "LINC01409"                 "LNCATV"                    "TNFRSF14-AS1"             
# [19] "MSANTD2-AS1"               "LINC02481"                 "LINC02288"                
# [22] "UBE2D3-AS1"                "SIMC1P1"                   "LINC02217"                
# [25] "LINC02362"                 "LINC02218"                 "ESAM-AS1"                 
# [28] "LINC01303"                 "LINC02482"                 "CHILL1"                   
# [31] "ADM-DT"                    "LINC02751"                 "EMSY-DT"                  
# [34] "CD44-AS1"                  "MIRLET7IHG"                "RPAP3-DT"                 
# [37] "LINC03108"                 "PPIAP45"                   "LINC02207"                
# [40] "LINC02289"                 "LINC01147"                 "MIDEAS-AS1"               
# [43] "ZNF710-AS1"                "LINC02352"                 "LINC02284"                
# [46] "LINC03034"                 "MPHOSPH10P1"               "BMP2K-DT"                 
# [49] "CEROX1"                    "NFATC2IP-AS1"              "LINC02175"                
# [52] "RSL1D1-DT"                 "NFE2L1-DT"                 "ZNF236-DT"                
# [55] "TWSG1-DT"                  "ATP8B1-AS1"                "SNHG30"                   
# [58] "TAF5LP1"                   "FOXP1-DT"                  "CKAP2LP1"                 
# [61] "IRF2-DT"                   "TPM4P2"                    "LINC02084"                
# [64] "HNRNPD-DT"                 "WDR5B-DT"                  "LINC03073"                
# [67] "RAB33B-AS1"                "TMEM272"  

## Anyways, I honestly don't have time to investigate this. So I'm moving on.