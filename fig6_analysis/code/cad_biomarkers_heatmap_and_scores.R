
# Coronary artery disease (CAD)
# blood RNA biomarker candidates
# cluster RNA's based on various RNA properties
# display results as a heatmap
# give each RNA a score based on a modified rank product across the various RNA properties

# read in table
# this was output from the app after applying "Disease or trait" filter
# using the term "Coronary artery disease"
# on the Discover RNA Biomarker Candidates page of the application

df <- read.csv("~/Downloads/filtered_gtex_gwas_table.csv")

# Get statistics about the CAD search results

# number of state genes
state<-read.csv("state_gene_list.csv")
sum(unique(df$Symbol.of.blood.RNA) %in% state$x)
# [1] 34
# number of trait genes
trait<-read.csv("trait_gene_list.csv")
sum(unique(df$Symbol.of.blood.RNA) %in% trait$x)
# [1] 48
# number of uncategorized genes
length(unique(df$Symbol.of.blood.RNA)) - sum(unique(df$Symbol.of.blood.RNA) %in% trait$x) - sum(unique(df$Symbol.of.blood.RNA) %in% state$x)
# [1] 97

snps<-df[!duplicated(df$rsID.of.eQTL),]
table(snps$Genomic.Context)

# 3_prime_UTR_variant 
# 14 
# 5_prime_UTR_variant 
# 2 
# intergenic_variant 
# 40 
# intron_variant 
# 238 
# missense_variant 
# 19 
# non_coding_transcript_exon_variant 
# 24 
# regulatory_region_variant 
# 18 
# splice_donor_5th_base_variant 
# 1 
# splice_region_variant 
# 1 
# synonymous_variant 
# 1 
# TF_binding_site_variant 
# 4

# Calculate the -log10 p values from the GWAS and eqtl p values
df$`-log10.eQTL.beta.p_value` <- -log(df$eQTL.beta.p_value, 10)
df$`-log10.GWAS.p_value` <- -log(df$GWAS.p_value, 10)

# retain the RNA features most important for prioritizing candidates
df2 <- df[,c("Symbol.of.blood.RNA", "rsID.of.eQTL",
             "-log10.eQTL.beta.p_value", "-log10.GWAS.p_value",
             "Studies.below.0.05.p_value", "Trait.4..filter.studies",
             "Median.Average.Normalized.Expression.Percentile",
             "Median.Repeatability.Percentile",
             "Median.Genetic.Variance.Percentile",
             "Number.of.trait.associations.in.GWAS")]

# rename columns with clearer titles for heatmap display 
names(df2) <- c("Symbol of Blood RNA", "CAD eQTL",
                "eQTL p_val", "GWAS p_val",
                "#State Studies", "#Trait Studies",
                "Expression Level",
                "Repeatability",
                "Genetic Var",
                "#Trait Assoc")

# round numbers for heatmap display 
df2[,3:ncol(df2)] <- apply(df2[,3:ncol(df2)], 2, function(x) round(x, 2))

# make numeric portion a matrix
mat <- as.matrix(df2[,3:ncol(df2)])
row.names(mat) <- df2$`Symbol of Blood RNA`
write.csv(mat, row.names = T, "CAD_biomarkers_heatmap_mat.csv")


# remove duplicate genes
# for each duplicate, keep gene with strongest association with eqtl
# first, sort by -log10 eqtl pval
# then get unique genes only
mat <- mat[order(mat[,"GWAS p_val"], decreasing = T),]
mat <- mat[order(mat[,"eQTL p_val"], decreasing = T),]
dim(mat)
# [1] 617   8

mat_dedup <- mat[!duplicated(row.names(mat)),]
dim(mat_dedup)
# [1] 179   8
write.csv(mat_dedup, row.names = T, "CAD_biomarkers_heatmap_mat_deduplicated.csv")


mat_state <- mat[row.names(mat) %in% state$x,]
dim(mat_state)
# [1] 88  8
mat_state_dedup <- mat_state[!duplicated(row.names(mat_state)),]
dim(mat_state_dedup)
# [1] 34  8

mat_trait <- mat[row.names(mat) %in% trait$x,]
dim(mat_trait)
# [1] 235   8
mat_trait_dedup <- mat_trait[!duplicated(row.names(mat_trait)),]
dim(mat_trait_dedup)
# [1] 48  8

# make heatmaps
clrsp <- colorRampPalette(c("white", "lemonchiffon1", "goldenrod1"))   
clrs <- clrsp(200) 
breaks1 <- seq(-2, 2, length.out = 200)

pheatmap::pheatmap(mat_trait_dedup, 
                   scale = "column", 
                   cluster_cols = T,
                   angle_col = 45, 
                   fontsize_row = 14,
                   fontsize_col = 18,
                   border_color = "grey",
                   display_numbers = mat_trait_dedup,
                   filename = "CAD_biomarkers_heatmap_48_trait_genes_with_trait_associations.png",
                   breaks = breaks1,
                   color =  clrs,
                   cellwidth = 50,
                   cellheight = 15,
                   height = 13,
                   width = 9)



# give each gene a score

ranks <- list(eQTL = rank(mat_dedup[,"eQTL p_val"],ties.method= "average"),
              GWAS = rank(mat_dedup[,"GWAS p_val"],ties.method= "average"),
              expression = rank(mat_dedup[,"Expression Level"],ties.method= "average"),
              repeatability = rank(mat_dedup[,"Repeatability"],ties.method= "average"),
              genetic = rank(mat_dedup[,"Genetic Var"],ties.method= "average"),
              trait = rank(mat_dedup[,"#Trait Studies"],ties.method= "average"),
              state = rank(-mat_dedup[,"#State Studies"],ties.method= "average"),
              spec = rank(-mat_dedup[,"#Trait Assoc"],ties.method= "average"))

str(ranks)
mat_dedup_ranks <- as.data.frame(do.call(cbind, ranks))
dim(mat_dedup_ranks)
head(mat_dedup_ranks)
mat_dedup_ranks$rank_prod <- apply(mat_dedup_ranks, 1, function(x) prod(x)^(1/ncol(mat_dedup_ranks)))
mat_dedup_ranks <- mat_dedup_ranks[order(mat_dedup_ranks$rank_prod, decreasing = T),]
head(mat_dedup_ranks,20)

write.csv(mat_dedup_ranks, row.names = T, "~/Downloads/CAD_biomarkers_scores.csv")

pdf("~/Downloads/CAD_biomarkers_scores.pdf", height = 6, width = 6)
plot(mat_dedup_ranks$rank_prod, pch=1, col = "black", 
     cex = 0.5, 
     main = "CAD Biomarker Candidate Scores",
     xlab="RNA Biomarker Candidate", 
     ylab="Rank Product")
text(x = 30, y = 133, 
     labels = paste(row.names(mat_dedup_ranks)[1:10],collapse = "\n"),
       cex = 0.7)
dev.off()

paste(row.names(mat_dedup_ranks)[1:10], collapse = ", ")


