# cluster and make heatmap of CAD blood RNA biomarker
# using temporal expression data on our app

# read in table
# this was output from the app after applying the term 
# cardiovascular disease in the Disease or trait to filer by box
# on the Discover RNA Biomarker Candidates page of the application

df <- read.csv("~/Downloads/filtered_gtex_gwas_table (7).csv")
df2 <- read.csv("~/Downloads/cad_percent_179_genes.csv") # see note at bottom for method from Will's email
df$Non.CAD.Traits <- df2$All.traits[match(df$Symbol.of.blood.RNA,df2$RNA)]-df2$CAD.traits[match(df$Symbol.of.blood.RNA,df2$RNA)]

# Calucalte the -log10 p values from the GWAS and eqtl p values
df$`-log10.eQTL.beta.p_value` <- -log(df$eQTL.beta.p_value, 10)
df$`-log10.GWAS.p_value` <- -log(df$GWAS.p_value, 10)

df2 <- df[,c("Symbol.of.blood.RNA", "rsID.of.eQTL",
             "-log10.eQTL.beta.p_value", "-log10.GWAS.p_value",
             "Studies.below.0.05.p_value", "Trait.4..filter.studies",
             "Median.Average.Normalized.Expression.Percentile",
             "Median.Repeatability.Percentile",
             "Non.CAD.Traits")]

names(df2) <- c("Symbol of Blood RNA", "CAD eQTL",
                "eQTL p_val", "GWAS p_val",
                "#State Studies", "#Trait Studies",
                "Expression Level",
                "Repeatability",
                "#Non-CAD Traits")

# it is hard to see the heatmap with this. many genes
# filter by trait and state studies first
df3 <- df[df$Trait.4..filter.studies>=3 & df$Studies.below.0.05.p_value <=3,
          c("Symbol.of.blood.RNA", "rsID.of.eQTL",
            "-log10.eQTL.beta.p_value", "-log10.GWAS.p_value",
            "Studies.below.0.05.p_value", "Trait.4..filter.studies",
            "Median.Average.Normalized.Expression.Percentile",
            "Median.Repeatability.Percentile",
            "Non.CAD.Traits")]
names(df3) <- c("Symbol of Blood RNA", "CAD eQTL",
                "eQTL p_val", "GWAS p_val",
                "#State Studies", "#Trait Studies",
                "Expression Level",
                "Repeatability",
                "#Non-CAD Traits")


# round numbers
df2[,3:ncol(df2)] <- apply(df2[,3:ncol(df2)], 2, function(x) round(x, 2))
df3[,3:ncol(df3)] <- apply(df3[,3:ncol(df3)], 2, function(x) round(x, 2))

# make numeric portion a matrix
mat <- as.matrix(df2[,3:ncol(df2)])
row.names(mat) <- df2$`Symbol of Blood RNA`
write.csv(mat, row.names = T, "~/Downloads/CAD_biomarkers_heatmap_mat.csv")

mat2 <- as.matrix(df3[,3:ncol(df3)])
row.names(mat2) <- df3$`Symbol of Blood RNA`
write.csv(mat2, row.names = T, "~/Downloads/CAD_biomarkers_heatmap_mat_filtered.csv")


# remove duplicate genes, but which one to pick for heatmap?
# pick one with strongest association with eqtl
# first, sort by -log10 eqtl pval
# then get unique genes only
mat <- mat[order(mat[,"GWAS p_val"], decreasing = T),]
mat <- mat[order(mat[,"eQTL p_val"], decreasing = T),]
dim(mat)
# [1] 617   7

mat_dedup <- mat[!duplicated(row.names(mat)),]
dim(mat_dedup)
# [1] 179   7
write.csv(mat_dedup, row.names = T, "~/Downloads/CAD_biomarkers_heatmap_mat_deduplicated.csv")


mat2 <- mat2[order(mat2[,"GWAS p_val"], decreasing = T),]
mat2 <- mat2[order(mat2[,"eQTL p_val"], decreasing = T),]
dim(mat2)
# [1] 329   7

mat2_dedup <- mat2[!duplicated(row.names(mat2)),]
dim(mat2_dedup)
# [1] 81   7
write.csv(mat2_dedup, row.names = T, "~/Downloads/CAD_biomarkers_heatmap_mat_filtered_deduplicated.csv")


# make heatmaps
clrsp <- colorRampPalette(c("white", "lemonchiffon1", "goldenrod1"))   
clrs <- clrsp(200) 
breaks1 <- seq(-2, 2, length.out = 200)

pheatmap::pheatmap(mat2_dedup, scale = "column", 
                   angle_col = 45, fontsize = 19,
                   display_numbers = mat2_dedup,
                   # cutree_rows = 10,
                   breaks = breaks1,
                   color =  clrs,
                   filename = "~/Downloads/CAD_biomarkers_heatmap_mat_filtered_deduplicated.png",
                   height = 20,
                   width = 12)

pheatmap::pheatmap(mat_dedup, scale = "column", 
                   angle_col = 45, fontsize = 19,
                   display_numbers = mat_dedup,
                   # cutree_rows = 10,
                   breaks = breaks1,
                   color =  clrs,
                   filename = "~/Downloads/CAD_biomarkers_heatmap_mat_deduplicated.png",
                   height = 40,
                   width = 20)


# give each gene a score
# based on a modified rank product 

ranks <- list(eQTL = rank(mat_dedup[,"eQTL p_val"],ties.method= "average"),
              GWAS = rank(mat_dedup[,"GWAS p_val"],ties.method= "average"),
              trait = rank(mat_dedup[,"#Trait Studies"],ties.method= "average"),
              expression = rank(mat_dedup[,"Expression Level"],ties.method= "average"),
              repeatability = rank(mat_dedup[,"Repeatability"],ties.method= "average"),
              state = rank(-mat_dedup[,"#State Studies"],ties.method= "average"),
              spec = rank(-mat_dedup[,"#Non-CAD Traits"],ties.method= "average"))
str(ranks)
mat_dedup_ranks <- as.data.frame(do.call(cbind, ranks))
dim(mat_dedup_ranks)
head(mat_dedup_ranks)
mat_dedup_ranks$rank_prod <- apply(mat_dedup_ranks, 1, function(x) prod(x)^(1/ncol(mat_dedup_ranks)))
mat_dedup_ranks <- mat_dedup_ranks[order(mat_dedup_ranks$rank_prod, decreasing = T),]
head(mat_dedup_ranks)

write.csv(mat_dedup_ranks, row.names = T, "~/Downloads/CAD_biomarkers_scores.csv")

pdf("~/Downloads/CAD_biomarkers_scores.pdf")
plot(mat_dedup_ranks$rank_prod, pch=1, col = "black", 
     cex = 0.5, 
     main = "CAD Biomarker Candidate Scores",
     xlab="RNA Biomarker Candidate", 
     ylab="Rank Product")
text(x = 30, y = 135, 
     labels = paste(row.names(mat_dedup_ranks)[1:10],collapse = "\n"),
       cex = 0.7)
dev.off()

# legend("topright", legend = row.names(mat_dedup_ranks)[1:11])
# textxy(X = c(1:10), Y = mat_dedup_ranks$rank_prod[1:10], 
#        labs = row.names(mat_dedup_ranks)[1:10],
#        cex = 0.2)
# library(ggrepel)
# pdf("~/Downloads/CAD_biomarkers_scores.pdf")
# ggplot(mat_dedup_ranks, aes(x=seq(nrow(mat_dedup_ranks)),y=rank_prod)) +
#   geom_point(alpha = 0.2) +
#   theme_classic() +
#   labs(title = "CAD Biomarker Candidate Scores",
#        x = "RNA Biomarker Candidate", 
#        y = "Rank Product") 
# # +geom_text_repel(label = c(row.names(mat_dedup_ranks)[1:10], rep("",nrow(mat_dedup_ranks)-10)), cex = 1)
# dev.off()



