############################################################
# cluster and make heatmap of CAD blood RNA biomarker
# using temporal expression data on our app
############################################################

options(stringsAsFactors = FALSE)

input_dir <- "~/Desktop/work_repo/github/biomarker_discovery/input/"
output_dir <- "~/Desktop/work_repo/github/biomarker_discovery/output/"

# read in table
# this was output from the app after applying the term 
# cardiovascular disease in the Disease or trait to filer by box
# on the Discover RNA Biomarker Candidates page of the application

df <- read.csv(paste(input_dir, "filtered_gtex_gwas_table.csv", sep = ""))

# Calucalte the -log10 p values from the GWAS and eqtl p values

df$`-log10.eQTL.beta.p_value` <- -log(df$eQTL.beta.p_value, 10)
df$`-log10.GWAS.p_value` <- -log(df$GWAS.p_value, 10)

df2 <- df[, c("Symbol.of.blood.RNA", "rsID.of.eQTL",
              "-log10.eQTL.beta.p_value", "-log10.GWAS.p_value",
              "Number.of.trait.associations.in.GWAS",
              "Studies.below.0.05.p_value", "Trait.4..filter.studies",
              "Average.Expression",
              "repeatability", "Trait.gene.", "State.gene.")]

names(df2) <- c("Symbol of Blood RNA", "CAD eQTL",
                "-log10 eQTL p_val", "-log10 GWAS p_val", "#GWAS Associations",
                "#State Studies", "#Trait Studies",
                "Expression Level",
                "Repeatability", "Is Trait", "Is State")

# make gene-level summary table
# this keeps the strongest CAD association but also retains how many CAD eQTLs support each gene

genes <- sort(unique(df2$`Symbol of Blood RNA`))

gene_summary <- data.frame(
  `Symbol of Blood RNA` = genes,
  check.names = FALSE
)

gene_summary$`n CAD eQTLs` <- NA
gene_summary$`-log10 eQTL p_val` <- NA
gene_summary$`-log10 GWAS p_val` <- NA
gene_summary$`#GWAS Associations` <- NA
gene_summary$`#State Studies` <- NA
gene_summary$`#Trait Studies` <- NA
gene_summary$`Expression Level` <- NA
gene_summary$`Repeatability` <- NA
gene_summary$`Is Trait` <- NA
gene_summary$`Is State` <- NA

for (i in seq_along(genes)) {
  
  gene_now <- genes[i]
  tmp <- df2[df2$`Symbol of Blood RNA` == gene_now, ]
  
  gene_summary$`n CAD eQTLs`[i] <- length(unique(tmp$`CAD eQTL`))
  gene_summary$`-log10 eQTL p_val`[i] <- max(tmp$`-log10 eQTL p_val`, na.rm = TRUE)
  gene_summary$`-log10 GWAS p_val`[i] <- max(tmp$`-log10 GWAS p_val`, na.rm = TRUE)
  gene_summary$`#GWAS Associations`[i] <- min(tmp$`#GWAS Associations`, na.rm = TRUE)
  gene_summary$`#State Studies`[i] <- max(tmp$`#State Studies`, na.rm = TRUE)
  gene_summary$`#Trait Studies`[i] <- max(tmp$`#Trait Studies`, na.rm = TRUE)
  gene_summary$`Expression Level`[i] <- max(tmp$`Expression Level`, na.rm = TRUE)
  gene_summary$`Repeatability`[i] <- max(tmp$`Repeatability`, na.rm = TRUE)
  gene_summary$`Is Trait`[i] <- any(tmp$`Is Trait`, na.rm = TRUE)
  gene_summary$`Is State`[i] <- any(tmp$`Is State`, na.rm = TRUE)
}

all_df <- gene_summary[, c("Symbol of Blood RNA",
                           "n CAD eQTLs",
                           "-log10 eQTL p_val",
                           "-log10 GWAS p_val",
                           "#GWAS Associations",
                           "#State Studies",
                           "#Trait Studies",
                           "Expression Level",
                           "Repeatability")]

trait_df <- gene_summary[gene_summary$`Is Trait`,
                         c("Symbol of Blood RNA",
                           "n CAD eQTLs",
                           "-log10 eQTL p_val",
                           "-log10 GWAS p_val",
                           "#GWAS Associations",
                           "Expression Level",
                           "Repeatability")]

state_df <- gene_summary[gene_summary$`Is State`,
                         c("Symbol of Blood RNA",
                           "n CAD eQTLs",
                           "-log10 eQTL p_val",
                           "-log10 GWAS p_val",
                           "#GWAS Associations",
                           "Expression Level",
                           "Repeatability")]

other_df <- gene_summary[!gene_summary$`Is Trait` & !gene_summary$`Is State`,
                         c("Symbol of Blood RNA",
                           "n CAD eQTLs",
                           "-log10 eQTL p_val",
                           "-log10 GWAS p_val",
                           "#GWAS Associations",
                           "Expression Level",
                           "Repeatability")]

# round numbers

all_df[, 2:ncol(all_df)] <- apply(all_df[, 2:ncol(all_df)], 2, function(x) round(x, 2))
trait_df[, 2:ncol(trait_df)] <- apply(trait_df[, 2:ncol(trait_df)], 2, function(x) round(x, 2))
state_df[, 2:ncol(state_df)] <- apply(state_df[, 2:ncol(state_df)], 2, function(x) round(x, 2))
other_df[, 2:ncol(other_df)] <- apply(other_df[, 2:ncol(other_df)], 2, function(x) round(x, 2))

# make numeric portion a matrix

all_df_mat <- as.matrix(all_df[, 2:ncol(all_df)])
row.names(all_df_mat) <- all_df$`Symbol of Blood RNA`

trait_df_mat <- as.matrix(trait_df[, 2:ncol(trait_df)])
row.names(trait_df_mat) <- trait_df$`Symbol of Blood RNA`

state_df_mat <- as.matrix(state_df[, 2:ncol(state_df)])
row.names(state_df_mat) <- state_df$`Symbol of Blood RNA`

other_df_mat <- as.matrix(other_df[, 2:ncol(other_df)])
row.names(other_df_mat) <- other_df$`Symbol of Blood RNA`

# score each gene
# higher scores are better for all components
# specificity is calculated as the inverse of the number of GWAS trait associations

library(ggplot2)
library(ggrepel)
library(pheatmap)
library(cowplot)
library(magick)
library(grid)

rescale01 <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

scores <- data.frame(
  GWAS = rescale01(all_df_mat[, "-log10 GWAS p_val"]),
  eQTL = rescale01(all_df_mat[, "-log10 eQTL p_val"]),
  # specificity = rescale01(-all_df_mat[, "#GWAS Associations"]),
  # expression = rescale01(-all_df_mat[, "Expression Level"]),
  # repeatability = rescale01(-all_df_mat[, "Repeatability"]),
  check.names = FALSE
)


row.names(scores) <- row.names(all_df_mat)

# add small floor so one zero component does not make the whole score zero

# score_cols <- c("GWAS", "eQTL", "specificity", "expression", "repeatability")
# score_cols <- c("GWAS", "eQTL", "expression", "repeatability")
# score_cols <- c("GWAS", "eQTL", "expression")
score_cols <- c("GWAS", "eQTL")

scores_for_plot <- scores

scores_for_plot[, score_cols] <- apply(
  scores_for_plot[, score_cols],
  2,
  function(x) pmax(x, 0.01)
)

scores_for_plot$combined_score <- apply(
  scores_for_plot[, score_cols],
  1,
  function(x) prod(x)^(1 / length(x))
)

scores_for_plot <- scores_for_plot[order(scores_for_plot$combined_score, decreasing = TRUE), ]

scores <- scores_for_plot

# write.csv(scores, row.names = TRUE, "~/Downloads/TableS4_CAD_biomarkers_scores.csv")
write.csv(scores, row.names = TRUE, paste(output_dir, "CAD_biomarkers_scores.csv", sep = ""))

# order heatmaps by combined score

all_df_mat <- all_df_mat[row.names(scores), ]

trait_df_mat <- trait_df_mat[intersect(row.names(scores), row.names(trait_df_mat)), ]
state_df_mat <- state_df_mat[intersect(row.names(scores), row.names(state_df_mat)), ]
other_df_mat <- other_df_mat[intersect(row.names(scores), row.names(other_df_mat)), ]

# make category annotation
# this version uses the matrices above, so it works even if you change the gene summary code

gene_category <- rep("Other", length(row.names(scores)))
names(gene_category) <- row.names(scores)

gene_category[names(gene_category) %in% row.names(state_df_mat)] <- "State"
gene_category[names(gene_category) %in% row.names(trait_df_mat)] <- "Trait"

state_and_trait_genes <- intersect(row.names(state_df_mat), row.names(trait_df_mat))
gene_category[names(gene_category) %in% state_and_trait_genes] <- "State and Trait"

# summarize rankings by state, trait, other category

rank_summary <- data.frame(
  gene = row.names(scores),
  rank = seq_len(nrow(scores)),
  combined_score = scores$combined_score,
  category = gene_category[row.names(scores)]
)

print(table(rank_summary$category))
print(table(rank_summary$category[rank_summary$rank <= 50]))

print(aggregate(rank ~ category, rank_summary, median))
print(aggregate(combined_score ~ category, rank_summary, median))

print(wilcox.test(rank ~ category, data = rank_summary[rank_summary$category %in% c("Trait", "State"), ]))
print(wilcox.test(rank ~ category, data = rank_summary[rank_summary$category %in% c("Trait", "Other"), ]))
print(wilcox.test(rank ~ category, data = rank_summary[rank_summary$category %in% c("State", "Other"), ]))

# make ranked candidate score plot

plot_df <- data.frame(
  gene = row.names(scores),
  combined_score = scores$combined_score
)

plot_df$rank <- seq_len(nrow(plot_df))

top_n_label <- 9
plot_df$top_gene <- plot_df$rank <= top_n_label

p_score <- ggplot(plot_df, aes(x = rank, y = combined_score)) +
  geom_point(size = 0.65, color = "black") +
  geom_point(
    data = plot_df[plot_df$top_gene, ],
    color = "firebrick",
    size = 1.3
  ) +
  geom_vline(
    xintercept = top_n_label,
    linetype = "dashed",
    color = "gray55",
    linewidth = 0.35
  ) +
  geom_text_repel(
    data = plot_df[plot_df$top_gene, ],
    aes(label = gene),
    size = 3,
    box.padding = 0.35,
    point.padding = 0.25,
    force = 1,
    min.segment.length = 0,
    segment.color = "black",
    segment.size = 0.3,
    max.overlaps = Inf,
    seed = 1
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "Biomarker Candidate Prioritization",
    x = "RNA Biomarker Candidate Rank",
    y = "Prioritization Score"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.08))) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.08))) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13),
    panel.grid.minor = element_blank()
  )

# make top 50 heatmap
# use score components so that darker color = stronger contribution to prioritization

top_n <- min(50, nrow(scores))
top_genes <- row.names(scores)[1:top_n]

# top_score_mat <- as.matrix(scores[top_genes, c(
#   "GWAS",
#   "eQTL",
#   "expression",
#   "repeatability",
#   "specificity"
# )])

# top_score_mat <- as.matrix(scores[top_genes, c(
#   "GWAS",
#   "eQTL",
#   "expression",
#   "repeatability"
# )])

# top_score_mat <- as.matrix(scores[top_genes, c(
#   "GWAS",
#   "eQTL",
#   "expression"
# )])

top_score_mat <- as.matrix(scores[top_genes, c(
  "GWAS",
  "eQTL"
)])

# colnames(top_score_mat) <- c(
#   "CAD GWAS",
#   "Blood eQTL",
#   "Expression",
#   "Repeatability",
#   "Specificity\n(low pleiotropy)"
# )

# colnames(top_score_mat) <- c(
#   "CAD GWAS",
#   "Blood eQTL",
#   "Expression",
#   "Repeatability"
# )

# colnames(top_score_mat) <- c(
#   "CAD GWAS",
#   "Blood eQTL",
#   "Expression"
# )

colnames(top_score_mat) <- c(
  "CAD GWAS",
  "Blood eQTL"
)

# make row annotation for state, trait, other

top_annotation <- data.frame(
  Category = gene_category[top_genes]
)

row.names(top_annotation) <- top_genes

# improved color palette
# darker = better biomarker candidate

clrs <- colorRampPalette(c("white", "lightsteelblue1", "steelblue"))(200)

top_score_mat_scaled <- apply(top_score_mat, 2, function(x) {
  (x - min(x)) / (max(x) - min(x))
})

p_heat <- pheatmap::pheatmap(
  top_score_mat_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = top_annotation,
  annotation_names_row = FALSE,
  annotation_legend = TRUE,
  angle_col = 45,
  fontsize = 8,
  fontsize_row = 6,
  fontsize_col = 8.5,
  border_color = "gray85",
  display_numbers = FALSE,
  color = clrs,
  breaks = seq(0.1, 1, length.out = 201),
  main = "Top 50 Prioritized Candidates",
  silent = TRUE
)

# make category score plot

p_category <- ggplot(rank_summary, aes(x = category, y = combined_score)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.12, size = 0.6, alpha = 0.5) +
  theme_bw(base_size = 11) +
  labs(
    title = "Scores by Temporal Category",
    x = "",
    y = "Prioritization Score"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    panel.grid.minor = element_blank()
  )

# read schematic panel

schematic_path <- paste(input_dir, "fig3_schematic.pdf", sep = "")

schematic_img <- magick::image_read_pdf(
  path = schematic_path,
  density = 600
)

schematic_img <- magick::image_trim(schematic_img)

p_schematic <- ggdraw() +
  draw_image(
    image = schematic_img,
    x = 0,
    y = 0,
    width = 1,
    height = 1
  )

# convert heatmap to drawable object

p_heat_grob <- ggdraw() +
  draw_grob(p_heat$gtable)

# add panel labels

panel_label_size <- 16

p_schematic_labeled <- p_schematic +
  draw_label("A", x = 0.01, y = 0.98, hjust = 0, vjust = 1,
             fontface = "bold", size = panel_label_size)

p_score_labeled <- ggdraw(p_score) +
  draw_label("B", x = 0.01, y = 0.98, hjust = 0, vjust = 1,
             fontface = "bold", size = panel_label_size)

p_heat_labeled <- p_heat_grob +
  draw_label("C", x = 0.01, y = 0.99, hjust = 0, vjust = 1,
             fontface = "bold", size = panel_label_size)

p_category_labeled <- ggdraw(p_category) +
  draw_label("D", x = 0.01, y = 0.98, hjust = 0, vjust = 1,
             fontface = "bold", size = panel_label_size)

# assemble layout

# top row: schematic (smaller) + score plot (larger)
top_row <- plot_grid(
  p_schematic_labeled,
  p_score_labeled,
  ncol = 2,
  rel_widths = c(0.7, 1.3),
  align = "hv"
)

# bottom row: heatmap (larger) + category plot (smaller)
bottom_row <- plot_grid(
  p_heat_labeled,
  p_category_labeled,
  ncol = 2,
  rel_widths = c(1.4, 0.8),
  align = "hv"
)

# combine rows
final_fig <- plot_grid(
  top_row,
  bottom_row,
  ncol = 1,
  rel_heights = c(0.85, 1.15)
)


# save final figure

ggsave(
  filename = paste(output_dir, "Figure3_CAD_biomarker_prioritization.png", sep = ""),
  plot = final_fig,
  width = 16,
  height = 10,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = paste(output_dir, "Figure3_CAD_biomarker_prioritization.pdf", sep = ""),
  plot = final_fig,
  width = 16,
  height = 10,
  bg = "white"
)


# # make heatmaps
# # values are column-scaled for visualization
# # displayed numbers are the original gene-level summary values
# 
# clrsp <- colorRampPalette(c("white", "lemonchiffon1", "goldenrod1"))   
# clrs <- clrsp(200) 
# breaks1 <- seq(-2, 2, length.out = 200)
# 
# pheatmap::pheatmap(trait_df_mat,
#                    scale = "column", 
#                    cluster_rows = TRUE,
#                    cluster_cols = TRUE,
#                    angle_col = 45,
#                    fontsize = 19,
#                    display_numbers = trait_df_mat,
#                    breaks = breaks1,
#                    color = clrs,
#                    filename = "~/Downloads/CAD_biomarkers_heatmap_trait.png",
#                    height = 20,
#                    width = 12)
# 
# pheatmap::pheatmap(state_df_mat,
#                    scale = "column", 
#                    cluster_rows = TRUE,
#                    cluster_cols = TRUE,
#                    angle_col = 45,
#                    fontsize = 19,
#                    display_numbers = state_df_mat,
#                    breaks = breaks1,
#                    color = clrs,
#                    filename = "~/Downloads/CAD_biomarkers_heatmap_state.png",
#                    height = 20,
#                    width = 12)
# 
# pheatmap::pheatmap(other_df_mat,
#                    scale = "column", 
#                    cluster_rows = TRUE,
#                    cluster_cols = TRUE,
#                    angle_col = 45,
#                    fontsize = 19,
#                    display_numbers = other_df_mat,
#                    breaks = breaks1,
#                    color = clrs,
#                    filename = "~/Downloads/CAD_biomarkers_heatmap_other.png",
#                    height = 20,
#                    width = 12)
# 
# # make a heatmap of the top-ranked candidates only
# # this may be the clearest version for the manuscript
# 
# top_n <- min(40, nrow(all_df_mat))
# top_genes <- row.names(scores)[1:top_n]
# top_df_mat <- all_df_mat[top_genes, ]
# 
# pheatmap::pheatmap(top_df_mat,
#                    scale = "column", 
#                    cluster_rows = FALSE,
#                    cluster_cols = TRUE,
#                    angle_col = 45,
#                    fontsize = 14,
#                    display_numbers = top_df_mat,
#                    breaks = breaks1,
#                    color = clrs,
#                    filename = "~/Downloads/CAD_biomarkers_heatmap_top40.png",
#                    height = 14,
#                    width = 12)
# 

