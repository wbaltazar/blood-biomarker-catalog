## January 5 2025
## Will Baltazar
## Our RShiny Application uses parquet files to load large tables in faster. This script tests file
## loading times in their original format and parquet format. Where time improvements were noted, files
## were rewrtiten to the application folder from this code.

setwd("~/Desktop/work_repo/Box organization/1results/RShiny-application/")
library(tidyverse)
library(nanoparquet)

# rds to parquet
e1 <- read_rds("./data/expr/gomez.rds")


e1_thing <- as.data.frame(e1)
str(e1_thing)

write_parquet(e1_thing,"./data/expr/gomez.parquet")
file.size("./data/expr/gomez.parquet") * 1e-6
file.size("./data/expr/gomez.rds") * 1e-6
system.time(read_parquet("./data/expr/gomez.parquet")) # this format isn't sufficient for application
# user  system elapsed 
# 0.006   0.002   0.008 
testfunc <- function() {
  obj <- read_parquet("./data/expr/gomez.parquet")
  obj <- as.matrix(obj)
  dimnames(obj)[[1]] <- dimnames(e1)[[1]]
}
system.time(testfunc())
# user  system elapsed 
# 0.011   0.002   0.014
system.time(read_rds("./data/expr/gomez.rds"))
# user  system elapsed 
# 0.011   0.001   0.012 they are practically the same
system.time(read_rds("./data/eqtls/gtex_variants_and_phenotypes_whole_blood.rds"))
# user  system elapsed 
# 1.645   0.094   1.817
system.time(read_parquet("./data/eqtls/gtex_variants_and_phenotypes_whole_blood.parquet"))
# user  system elapsed 
# 0.857   0.108   0.971 

thing <- read.csv(file = "./data/variation/1variation.csv")
dim(thing)
write_parquet(thing, "./data/variation/1variation.parquet")
system.time(read_parquet("./data/variation/1variation.parquet"))
# user  system elapsed 
# 0.020   0.004   0.023 this is 10 times faster!!!
system.time(read.csv(file = "./data/variation/1variation.csv"))
# user  system elapsed 
# 0.210   0.004   0.214 
file.remove("./data/variation/1variation.parquet")

files <- c("1variation.csv", "2variation.csv", "3variation.csv", "4variation.csv",
            "5variation.csv", "6variation.csv", "7variation.csv", "8variation.csv",
            "9variation.csv", "10variation.csv", "11variation.csv", "12variation.csv", "13variation.csv")
for (i in 1:13) {
  temp <- read.csv(paste("./data/variation/", files[i], sep = ""))
  write_parquet(temp, file = paste("./data/variation/", i, "variation.parquet", sep = ""))
}

files <- list.files(path = "./data/variation/", full.names = T)
print(files)
# [1] "./data/variation//10variation.csv"     "./data/variation//10variation.parquet" "./data/variation//11variation.csv"    
# [4] "./data/variation//11variation.parquet" "./data/variation//12variation.csv"     "./data/variation//12variation.parquet"
# [7] "./data/variation//13variation.csv"     "./data/variation//13variation.parquet" "./data/variation//1variation.csv"     
# [10] "./data/variation//1variation.parquet"  "./data/variation//2variation.csv"      "./data/variation//2variation.parquet" 
# [13] "./data/variation//3variation.csv"      "./data/variation//3variation.parquet"  "./data/variation//4variation.csv"     
# [16] "./data/variation//4variation.parquet"  "./data/variation//5variation.csv"      "./data/variation//5variation.parquet" 
# [19] "./data/variation//6variation.csv"      "./data/variation//6variation.parquet"  "./data/variation//7variation.csv"     
# [22] "./data/variation//7variation.parquet"  "./data/variation//8variation.csv"      "./data/variation//8variation.parquet" 
# [25] "./data/variation//9variation.csv"      "./data/variation//9variation.parquet" 

variation_tables <- list()
files1 <- c("1variation.csv", "2variation.csv", "3variation.csv", "4variation.csv",
           "5variation.csv", "6variation.csv", "7variation.csv", "8variation.csv",
           "9variation.csv", "10variation.csv", "11variation.csv", "12variation.csv", "13variation.csv")
files2 <- c("1variation.parquet", "2variation.parquet", "3variation.parquet", "4variation.parquet",
           "5variation.parquet", "6variation.parquet", "7variation.parquet", "8variation.parquet",
           "9variation.parquet", "10variation.parquet", "11variation.parquet", "12variation.parquet", "13variation.parquet")

testfunc2 <- function(files) {
  for (i in 1:length(files)) {
    variation_tables[[i]] <- read.csv(file = paste("./data/variation/", files[i], sep = ""))
    names(variation_tables)[i] <- study_names[i]
  } 
}
testfunc3 <- function(files) {
  for (i in 1:length(files)) {
    variation_tables[[i]] <- read_parquet(file = paste("./data/variation/", files[i], sep = ""))
    names(variation_tables)[i] <- study_names[i]
  } 
}

system.time(testfunc2(files1))
# user  system elapsed 
# 4.000   0.105   4.139 
system.time(testfunc3(files2))
# user  system elapsed 
# 0.248   0.038   0.286 16 times faster

filestoremove <- list.files("./data/variation/", full.names = T)
filestoremove <- filestoremove[grep('csv', filestoremove)]
file.remove(filestoremove)

system.time(read.csv("./data/characteristic_scores.csv"))
# user  system elapsed 
# 0.088   0.005   0.095 
temp <- read.csv("./data/characteristic_scores.csv")
write_parquet(temp, "./data/characteristic_scores.parquet")
system.time(read_parquet("./data/characteristic_scores.parquet"))
# user  system elapsed 
# 0.022   0.002   0.024
temp <- read.csv("./data/flexible_gene_scores.csv")
write_parquet(temp, "./data/flexible_gene_scores.parquet")

files <- list.files("./data/phenotypes/", full.names = T)
length(files)
# [1] 13
filenames <- list.files("./data/phenotypes/")
for (i in 1:13) {
  temp <- read.csv(files[i])
  write_parquet(temp, file = paste("./data/phenotypes/", str_extract(filenames[i], ".*pheno"), ".parquet", sep = ""))
}

files <- list.files("./data/phenotypes/", full.names = T)
length(files)
# [1] 26

pheno_data <- list()
files1 <- files[grep('csv', files)]
files2 <- files[grep('parquet', files)]
testfunc4 <- function(files) {
  for (i in 1:length(files)) {
    pheno_data[[i]] <- read.csv(file = paste(files[i], sep = ""))
    names(pheno_data)[i] <- study_names[i]
  }
}
testfunc5 <- function(files) {
  for (i in 1:length(files)) {
    pheno_data[[i]] <- read_parquet(file = paste(files[i], sep = ""))
    names(pheno_data)[i] <- study_names[i]
  }
}

system.time(testfunc4(files1))
# user  system elapsed 
# 0.027   0.001   0.029
system.time(testfunc5(files2))
# user  system elapsed 
# 0.073   0.002   0.075 # this is not faster
file.remove(files2)

## Annotation files
files <- list.files("./data/annotate/", full.names = T)
print(files)
# [1] "./data/annotate//affy.csv"    "./data/annotate//gomez.csv"   "./data/annotate//gosch.csv"   "./data/annotate//ilmn.csv"   
# [5] "./data/annotate//larocca.csv"
temp <- read.csv(files[[1]])
temp
write_parquet(temp, "./data/annotate/affy.parquet")
system.time(read.csv(files[[1]]))
# user  system elapsed 
# 0.218   0.010   0.228 
system.time(read_parquet("./data/annotate/affy.parquet"))
# user  system elapsed 
# 0.079   0.012   0.090 good time reduction
file.remove("./data/annotate/affy.parquet")
filenames <- list.files("./data/annotate/")
for (i in 1:5) {
  temp <- read.csv(files[i])
  write_parquet(temp, file = paste("./data/annotate/", str_replace(filenames[i], pattern = "csv", replacement = "parquet"), sep = ""))
}
files <- list.files("./data/annotate", full.names = T)
files1 <- files[grep("csv", files)]
files2 <- files[grep("parquet", files)]
testfunca <- function(files) {
  for (i in 1:length(files)) {
    anno_data[[i]] <- read.csv(files[[i]])
    anno_data[[i]] <- anno_data[[i]][,-1]
  }
}
testfuncb <- function(files) {
  for (i in 1:length(files)) {
    anno_data[[i]] <- read_parquet(files[[i]])
    anno_data[[i]] <- anno_data[[i]][,-1]
  }
}
system.time(testfunca(files1))
# user  system elapsed 
# 0.894   0.031   0.934 
system.time(testfuncb(files2))
# user  system elapsed 
# 0.264   0.036   0.338 
file.remove(files1)

system.time(read.csv("./data/median_stability_statistics_all.csv"))
# user  system elapsed 
# 0.154   0.010   0.163
med <- read.csv("./data/median_stability_statistics_all.csv")
write_parquet(med, file = "./data/median_stability_statistics_all.parquet")
system.time(read_parquet("./data/median_stability_statistics_all.parquet"))
# user  system elapsed 
# 0.019   0.003   0.022 much faster!

## For the characteristic vs. flexible genes, I will read in supplementary data file 2
## and turn it into a table, save as .parquet
supp <- readxl::read_xlsx("~/Desktop/work_repo/github/cross_study_analysis/output/supplementary_data_2.xlsx", sheet = 1)
supp2 <- readxl::read_xlsx("~/Desktop/work_repo/github/cross_study_analysis/output/supplementary_data_2.xlsx", sheet = 2)
head(supp) 
head(supp2)
supp$category <- "trait"
supp2$category <- "state"
newsupp <- rbind(supp, supp2)
table(newsupp$category)
# state trait 
# 1327  1356 
length(unique(newsupp$hgnc_symbol))
# [1] 2683
dim(newsupp)
# [1] 2683    3
## looks good, write it to parquet
write_parquet(newsupp, "./data/supplementary_data_2.parquet")