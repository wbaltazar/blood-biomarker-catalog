# Code and result repository for The Healthy Whole Blood Transcriptome Project

Below is a description of what each directory contains.

### (First author)\_study_data

These folders contain the raw data and code used to preprocess and perform analysis on the data.

-   code folder

    -   Contains the R scripts used to perform analysis.

    -   limma-QC R scripts perform pre-processing steps on the data as well as quality control. It also contains the code for differential expression analysis using limma.

        -   These files all take the files as they are stored on GEO (which you must unzip first) for analysis. As long as your input directory has access to these files, the code runs correctly.

    -   The variation R scripts perform the pre-processing steps as before and calculate all statistics needed for gene categorization.

-   output folder

    -   Contains all outputs from the R code. Contains various outputs, including CSV files and image files for figures.

### cross_study_analysis

-   code folder

    -   Scripts for finding

        -   trait genes (trait_genes.R)

        -   state genes (state_genes.R)

    -   gene_categories_functional_enrich_in_silico_validation.R finds genes which meet their respective categories criteria in 4 or more studies (one microarray, one RNA-seq) and creates a list of those genes in the form of a TXT file. Creates sub-directories in output for storing the outputs of in silico validation analysis and functional enrichment analysis. Performs functional enrichment analysis using clusterProfiler and produces all figures.

-   output folder

    -   Contains all outputs from the R code. PDFs, CSVs, PNGs, and TXT files. Includes supplementary data files.

### fig6_analysis

-   code folder

    -   cad_biomarkers_heatmap_and_scores.R calculates and organizes the parameters used for our prioritization of biomarker candidates for coronary artery disease (CAD) using outputs from the RShiny application.

    -   gene_categories_functional_enrich_in_silico_validation.R finds genes which meet their respective categories criteria in 4 or more studies (one microarray, one RNA-seq) and creates a list of those genes in the form of a TXT file. Creates sub-directories in output for storing the outputs of in silico validation analysis and functional enrichment analysis. Performs functional enrichment analysis using clusterProfiler and produces all figures.

-   input folder

    -   filtered_gtex_gwas_table.csv contains the filtered GTEx-GWAS pairs and calculated state and trait values for RNAs with "Coronary artery disease" as the filter. To replicate this using the application, navigate to the “Discover RNA Biomarker Candidates” page of the RShiny application and search for “Coronary artery disease” (CAD) with the “Disease/trait” radio button activated.

    -   state_gene_list.csv and trait_gene_list.csv can be downloaded from the "Download" page of the RShiny application.

-   output folder

    -   Contains all outputs from the R code. PNG and PDF for figure 6, plus the table used to generate these data in CSV format.

### RShiny-application

-   Contains all code and files needed to reproduce the current version of our RShiny application (<https://wcaar.shinyapps.io/blood-biomarker-catalog/>), but does not contain data files used due to storage issues

    -   Files can be downloaded directly from our application on the Download tab

    -   Alternatively, the code in the misc directory and study_data folders can be used to produce the necessary data and results used in the application. The folders should be organized such that the data is accessible in directories called by app.R and data.R (functions.R does not require any files).

-   misc

    -   building_eqtl_table.R contains code used to format the eQTL search table utilized by the eQTL variant search tab in the RShiny application. It requires the significant variant-gene pairs file from GTEx v10, the lookup table for the rsIDs from GTEx, and the most up-to-date version of the GWAS associations table as of December 14, 2024.

    -   biomaRt_all_transcripts.R contains code used to download gene information data used for the gene search tab in the RShiny application

    -   The diurnal_gene_analysis folder contains diurnal_genes.R and outputs. It compares the genes studied across 24 hours for two different datasets we used for this study.

    -   formatting_expr_for_app.R contains code used to create RDS files used for RShiny app spaghetti plots. It matches gene names to Probe IDs and Gene IDs for the Gene Info application function.

    -   for_common_symbols_in_app.R contains code used to create RDS files used in RShiny app spaghetti plots as well. It filters the genes in the RDS files so that only the 6,099 common symbols are used in the analysis.

    -   writing_parquet_files.R contains code that transforms files for the RShiny application into parquet files for faster loading.

-   As of February 18, 2025, this RShiny code is up to date.

## Session Info

R version 4.4.1 (2024-06-14) Platform: x86_64-apple-darwin20 Running under: macOS 15.1.1

Matrix products: default BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib LAPACK: /Library/Frameworks/R.framework/Versions/4.4-x86_64/Resources/lib/libRlapack.dylib; LAPACK version 3.12.0

locale: [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Chicago tzcode source: internal

attached base packages: [1] grid parallel stats4 stats graphics grDevices utils datasets methods base

other attached packages: [1] VennDiagram_1.7.3 futile.logger_1.4.3 variancePartition_1.34.0\
[4] BiocParallel_1.38.0 shiny_1.8.1.1 shinytoastr_2.2.0\
[7] shinyWidgets_0.8.6 shinythemes_1.2.0 shinyBS_0.61.1\
[10] rJava_1.0-11 readxl_1.4.3 png_0.1-8\
[13] PCAtools_2.16.0 oligo_1.68.2 Biostrings_2.72.1\
[16] XVector_0.44.0 oligoClasses_1.66.0 mailR_0.8\
[19] Homo.sapiens_1.3.1 TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 GO.db_3.19.1\
[22] OrganismDbi_1.46.0 hgu133plus2.db_3.13.0 hgu133a.db_3.13.0\
[25] heritability_1.4 MASS_7.3-61 glmnet_4.1-8\
[28] Matrix_1.7-0 EnsDb.Hsapiens.v86_2.99.0 ensembldb_2.28.0\
[31] AnnotationFilter_1.28.0 GenomicFeatures_1.56.0 GenomicRanges_1.56.1\
[34] GenomeInfoDb_1.40.1 edgeR_4.2.0 EnhancedVolcano_1.22.0\
[37] ggrepel_0.9.5 DT_0.33 cowplot_1.1.3\
[40] clusterProfiler_4.12.0 bslib_0.7.0 BiocManager_1.30.23\
[43] affycoretools_1.76.0 affy_1.82.0 pROC_1.18.5\
[46] caret_7.0-1 lattice_0.22-6 illuminaHumanv3.db_1.26.0\
[49] org.Hs.eg.db_3.19.1 AnnotationDbi_1.66.0 IRanges_2.38.0\
[52] S4Vectors_0.42.0 limma_3.60.3 GEOquery_2.72.0\
[55] Biobase_2.64.0 BiocGenerics_0.50.0 pheatmap_1.0.12\
[58] nanoparquet_0.3.1 lubridate_1.9.3 forcats_1.0.0\
[61] stringr_1.5.1 dplyr_1.1.4 purrr_1.0.2\
[64] readr_2.1.5 tidyr_1.3.1 tibble_3.2.1\
[67] ggplot2_3.5.1 tidyverse_2.0.0
