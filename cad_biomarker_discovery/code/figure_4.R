############################################################
# Visualization - Figure 4
############################################################

library(tidyverse)
library(cowplot)

traitResults <- read.csv("~/Desktop/work_repo/github/cad_biomarker_discovery/output/Trait_Performance.csv")
cadResults <- read.csv("~/Desktop/work_repo/github/cad_biomarker_discovery/output/CAD_Performance.csv")

############################################################
# Trait Genes
############################################################

traitResults <- traitResults %>% filter(dataset != str_remove(feature_set, "GENES"))
traitResults$dataset <- case_when(traitResults$dataset == "GSE1" ~ "GSE12288",
                                  traitResults$dataset == "GSE2" ~ "GSE202625",
                                  traitResults$dataset == "GSE3" ~ "GSE180081",
                                  traitResults$dataset == "GSE4" ~ "GSE20680",
                                  traitResults$dataset == "GSE5" ~ "GSE20681")
traitResults$feature_set <- case_when(traitResults$feature_set == "GSE1GENES" ~ "GSE12288-nominated genes",
                                      traitResults$feature_set == "GSE2GENES" ~ "GSE202625-nominated genes",
                                      traitResults$feature_set == "GSE3GENES" ~ "GSE180081-nominated genes",
                                      traitResults$feature_set == "GSE4GENES" ~ "GSE20680-nominated genes",
                                      traitResults$feature_set == "GSE5GENES" ~ "GSE20681-nominated genes",
                                      traitResults$feature_set == "TRAIT" ~ "Trait genes")
traitResults$algorithm <- case_when(traitResults$algorithm == "glmnet" ~ "Elastic Net",
                                    traitResults$algorithm == "nb" ~ "Naive Bayes",
                                    traitResults$algorithm == "pls" ~ "PLS-DA",
                                    traitResults$algorithm == "rf" ~ "Random Forest")

rank_df <- traitResults %>%
  group_by(dataset, algorithm) %>%
  mutate(
    rank = rank(-AUC_5fold10xCV, ties.method = "average")
  ) %>%
  ungroup()

rank_summary <- rank_df %>%
  group_by(feature_set, algorithm) %>%
  summarise(
    mean_rank = mean(rank),
    sd_rank   = sd(rank),
    .groups = "drop"
  )

fs_mean <- rank_summary %>%
  group_by(feature_set) %>%
  summarise(
    overall_mean_rank = mean(mean_rank),
    overall_sd_rank   = sd(mean_rank),
    .groups = "drop"
  )

## Plot
# pdf(file = "~/Desktop/work_repo/github/cad_biomarker_discovery/output/7B_Trait_Performance.pdf", width = 4, height = 7)
a <- ggplot() +
  # Algorithm-specific points (light, secondary)
  geom_point(
    data = rank_summary,
    aes(x = reorder(feature_set, mean_rank),
        y = mean_rank,
        color = algorithm),
    size = 2,
    alpha = 0.6,
    position = position_jitter(width = 0.15)
  ) +
  
  # Overall mean marker (primary)
  geom_point(
    data = fs_mean,
    aes(x = reorder(feature_set, overall_mean_rank),
        y = overall_mean_rank),
    color = "black",
    size = 4,
    shape = 18   # diamond
  ) +
  
  # Optional error bars for overall variability
  geom_errorbar(
    data = fs_mean,
    aes(x = reorder(feature_set, overall_mean_rank),
        ymin = overall_mean_rank - overall_sd_rank,
        ymax = overall_mean_rank + overall_sd_rank),
    width = 0.25,
    color = "black"
  ) +
  
  coord_flip() +
  labs(
    x = "Feature set",
    y = "Mean rank",
    title = "Trait Genes",
    # subtitle = "Black diamonds indicate mean rank averaged across algorithms",
    color = "Algorithm"
  ) +
  scale_y_continuous(
    breaks = seq(1, 6, 1),    # Major grid lines every 1 unit on the x-axis
    minor_breaks = NULL        # Remove minor grid lines on the x-axis
  ) +
  theme(legend.position = "none", panel.background = element_rect(fill = "white"), panel.grid = element_line(color = "grey80"))
# dev.off()

############################################################
# CAD Genes
############################################################

cadResults <- cadResults %>% filter(dataset != str_remove(feature_set, "GENES"))
cadResults$dataset <- case_when(cadResults$dataset == "GSE1" ~ "GSE12288",
                                cadResults$dataset == "GSE2" ~ "GSE202625",
                                cadResults$dataset == "GSE3" ~ "GSE180081",
                                cadResults$dataset == "GSE4" ~ "GSE20680",
                                cadResults$dataset == "GSE5" ~ "GSE20681")
cadResults$feature_set <- case_when(cadResults$feature_set == "GSE1GENES" ~ "GSE12288-nominated genes",
                                    cadResults$feature_set == "GSE2GENES" ~ "GSE202625-nominated genes",
                                    cadResults$feature_set == "GSE3GENES" ~ "GSE180081-nominated genes",
                                    cadResults$feature_set == "GSE4GENES" ~ "GSE20680-nominated genes",
                                    cadResults$feature_set == "GSE5GENES" ~ "GSE20681-nominated genes",
                                    cadResults$feature_set == "CAD" ~ "CAD-associated genes")
cadResults$algorithm <- case_when(cadResults$algorithm == "glmnet" ~ "Elastic Net",
                                  cadResults$algorithm == "nb" ~ "Naive Bayes",
                                  cadResults$algorithm == "pls" ~ "PLS-DA",
                                  cadResults$algorithm == "rf" ~ "Random Forest")

rank_df <- cadResults %>%
  group_by(dataset, algorithm) %>%
  mutate(
    rank = rank(-AUC_5fold10xCV, ties.method = "average")
  ) %>%
  ungroup()

rank_summary <- rank_df %>%
  group_by(feature_set, algorithm) %>%
  summarise(
    mean_rank = mean(rank),
    sd_rank   = sd(rank),
    .groups = "drop"
  )

fs_mean <- rank_summary %>%
  group_by(feature_set) %>%
  summarise(
    overall_mean_rank = mean(mean_rank),
    overall_sd_rank   = sd(mean_rank),
    .groups = "drop"
  )

## Plot
# pdf(file = "~/Desktop/work_repo/github/cad_biomarker_discovery/output/7C_CAD_Performance.pdf", width = 4, height = 7)
b <- ggplot() +
  # Algorithm-specific points (light, secondary)
  geom_point(
    data = rank_summary,
    aes(x = reorder(feature_set, mean_rank),
        y = mean_rank,
        color = algorithm),
    size = 2,
    alpha = 0.6,
    position = position_jitter(width = 0.15)
  ) +
  
  # Overall mean marker (primary)
  geom_point(
    data = fs_mean,
    aes(x = reorder(feature_set, overall_mean_rank),
        y = overall_mean_rank),
    color = "black",
    size = 4,
    shape = 18   # diamond
  ) +
  
  # Optional error bars for overall variability
  geom_errorbar(
    data = fs_mean,
    aes(x = reorder(feature_set, overall_mean_rank),
        ymin = overall_mean_rank - overall_sd_rank,
        ymax = overall_mean_rank + overall_sd_rank),
    width = 0.25,
    color = "black"
  ) +
  
  coord_flip() +
  labs(
    x = "Feature set",
    y = "Mean rank",
    title = "CAD-related Genes",
    # subtitle = "Black diamonds indicate mean rank averaged across algorithms",
    color = "Algorithm"
  ) +
  scale_y_continuous(
    breaks = seq(1, 6, 1),    # Major grid lines every 1 unit on the x-axis
    minor_breaks = NULL        # Remove minor grid lines on the x-axis
  ) +
  theme(legend.position = "none", panel.background = element_rect(fill = "white"), panel.grid = element_line(color = "grey80"))
# dev.off()

leg <- cowplot::get_plot_component(a + theme(legend.position = "left"), "guide-box", return_all = TRUE)
## legend is stored in box 2

pdf(file = "~/Desktop/work_repo/github/cad_biomarker_discovery/output/Figure7.pdf", width = 16, height = 4.5)
plot_grid(
  leg[[2]],
  plot_grid(a, b, nrow = 1, align = "h", axis = "l"),
  ncol = 2,
  rel_widths = c(0.1, 1)
)
dev.off()
