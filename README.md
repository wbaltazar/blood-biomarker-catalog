---
editor_options: 
  markdown: 
    wrap: 72
---

# Code and result repository for The Healthy Whole Blood Transcriptome Project

Below is a description of what each directory contains.

### (First author)\_study_data

These folders contain the raw data and code used to preprocess and
perform analysis on the data.

-   code folder

    -   Contains the R scripts used to perform analysis.

    -   limma-QC R scripts perform pre-processing steps on the data as
        well as quality control. It also contains the code for
        differential expression analysis using limma.

        -   Dusek has an extra script called pre_QC which should be run
            before the limma_QC file.

        -   These files all take the files as they are stored on GEO
            (which you must unzip first) for analysis. As long as your
            input directory has access to these files, the code runs
            correctly. Since some files are removed during quality
            control, the code may not run as expected the second time
            (you will receive non-existent file errors). You can get
            around this by running only the parts of the code needed for
            analysis or re-downloading the raw data from GEO and setting
            it as the input directory.

    -   The variation R scripts perform the pre-processing steps as
        before and calculate all statistics needed for gene
        categorization.

-   output folder

    -   Contains all outputs from the R code. Contains various outputs,
        including csv files and pdfs of figures.

### cross_study_analysis

-   code folder

    -   Scripts for finding

        -   stable-polymorphic genes (stable_polymorphic_genes.R)

        -   flexible genes (flexible_genes.R)

        -   housekeeping genes (housekeeping.R)

    -   gene_categories_functional_enrich_in_silico_validation.R finds
        genes which meet their respective categories criteria in 4 or
        more studies (one microarray, one RNA-seq) and creates a list of
        those genes in the form of a txt file. Creates sub-directories
        in output for storing the outputs of in silico validation
        analysis and functional enrichment analysis. Performs functional
        enrichment analysis using clusterProfiler and produces all
        figures.

-   output folder

    -   Contains all outputs from the R code. PDFs, CSVs, PNGs, and TXT
        files.

### misc

-   Affy_control_probe_check.R contains functions for assessing
    Affymetrix array control probe statistics

-   building_eqtl_table.R contains code used to format the eQTL search
    table utilized by the eQTL variant search tab in the RShiny
    application

-   biomaRt_all_transcripts.R contains code used to download gene
    information data used for the gene search tab in the RShiny
    application
    
-   formatting_expr_for_app.R contains code used to create RDS files 
    used for RShiny app spaghetti plots. It matches gene names to Probe
    IDs and Gene IDs for the Gene Info application function.

-   for_common_symbols_in_app.R contains code used to create RDS files 
    used in RShiny app spaghetti plots as well. It filters the genes in
    the RDS files so that only the 9474 common symbols are used in the
    analysis.

### RShiny-application

-   Contains all code and files needed to reproduce the current version
    of our RShiny application
    (<https://wcaar.shinyapps.io/blood-biomarker-catalog/>), but
    does not contain data files used due to storage issues

    -   Files can be downloaded directly from our application on the
        Download tab

    -   Alternatively, the code in the misc directory and study_data
        folders can be used to produce the necessary data and results
        used in the application

-   As of July 21, 2024, this code is out of date and will be updated
    soon.

## Session Info (and all packages used)

R version 4.4.1 (2024-06-14) Platform: x86_64-apple-darwin20 Running
under: macOS Sonoma 14.5

Matrix products: default BLAS:
/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK:
/Library/Frameworks/R.framework/Versions/4.4-x86_64/Resources/lib/libRlapack.dylib;
LAPACK version 3.12.0

locale: [1]
en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Chicago tzcode source: internal

attached base packages: [1] stats4 stats graphics grDevices utils
datasets methods base

other attached packages: [1] EnhancedVolcano_1.22.0 RColorBrewer_1.1-3
pheatmap_1.0.12 GEOquery_2.72.0 PCAtools_2.16.0\
[6] ggrepel_0.9.5 cowplot_1.1.3 lubridate_1.9.3 forcats_1.0.0
stringr_1.5.1\
[11] dplyr_1.1.4 purrr_1.0.2 readr_2.1.5 tidyr_1.3.1 tibble_3.2.1\
[16] ggplot2_3.5.1 tidyverse_2.0.0 oligo_1.68.2 Biostrings_2.72.1
GenomeInfoDb_1.40.1\
[21] XVector_0.44.0 oligoClasses_1.66.0 limma_3.60.3 lumi_2.56.0
illuminaHumanv3.db_1.26.0 [26] org.Hs.eg.db_3.19.1 illuminaio_0.46.0
beadarray_2.54.0 hexbin_1.28.3 AnnotationDbi_1.66.0\
[31] IRanges_2.38.0 S4Vectors_0.42.0 Biobase_2.64.0 BiocGenerics_0.50.0

loaded via a namespace (and not attached): [1] splines_4.4.1
BiocIO_1.14.0 bitops_1.0-7 R.oo_1.26.0 preprocessCore_1.66.0\
[6] methylumi_2.50.0 affxparser_1.76.0 XML_3.99-0.16.1 lifecycle_1.0.4
edgeR_4.2.0\
[11] lattice_0.22-6 MASS_7.3-61 base64_2.0.1 scrime_1.3.5
magrittr_2.0.3\
[16] minfi_1.50.0 yaml_2.3.8 doRNG_1.8.6 askpass_1.2.0 DBI_1.2.3\
[21] maps_3.4.2 abind_1.4-5 zlibbioc_1.50.0 quadprog_1.5-8
GenomicRanges_1.56.1\
[26] R.utils_2.12.3 RCurl_1.98-1.14 GenomeInfoDbData_1.2.12
irlba_2.3.5.1 genefilter_1.86.0\
[31] dqrng_0.4.1 annotate_1.82.0 DelayedMatrixStats_1.26.0
codetools_0.2-20 DelayedArray_0.30.1\
[36] xml2_1.3.6 tidyselect_1.2.1 farver_2.1.2 UCSC.utils_1.0.0
ScaledMatrix_1.12.0\
[41] beanplot_1.3.1 ash_1.0-15 matrixStats_1.3.0
GenomicAlignments_1.40.0 jsonlite_1.8.8\
[46] multtest_2.60.0 survival_3.7-0 iterators_1.0.14 foreach_1.5.2
tools_4.4.1\
[51] Rcpp_1.0.12 glue_1.7.0 Rttf2pt1_1.3.12 SparseArray_1.4.8
mgcv_1.9-1\
[56] DESeq2_1.44.0 MatrixGenerics_1.16.0 HDF5Array_1.32.0 withr_3.0.0
BiocManager_1.30.23\
[61] fastmap_1.2.0 rhdf5filters_1.16.0 fansi_1.0.6 openssl_2.2.0
rsvd_1.0.5\
[66] digest_0.6.36 BeadDataPackR_1.56.0 timechange_0.3.0 R6_2.5.1
colorspace_2.1-0\
[71] RSQLite_2.3.7 R.methodsS3_1.8.2 utf8_1.2.4 generics_0.1.3
data.table_1.15.4\
[76] rtracklayer_1.64.0 httr_1.4.7 htmlwidgets_1.6.4 S4Arrays_1.4.1
pkgconfig_2.0.3\
[81] gtable_0.3.5 rsconnect_1.3.1 blob_1.2.4 siggenes_1.78.0
Glimma_2.14.0\
[86] htmltools_0.5.8.1 scales_1.3.0 png_0.1-8 rstudioapi_0.16.0
tzdb_0.4.0\
[91] reshape2_1.4.4 rjson_0.2.21 nlme_3.1-165 curl_5.2.1
bumphunter_1.46.0\
[96] cachem_1.1.0 rhdf5_2.48.0 KernSmooth_2.23-24 parallel_4.4.1
extrafont_0.19\
[101] restfulr_0.0.15 pillar_1.9.0 grid_4.4.1 reshape_0.8.9 vctrs_0.6.5\
[106] BiocSingular_1.20.0 beachmat_2.20.0 ff_4.0.12 xtable_1.8-4
extrafontdb_1.0\
[111] GenomicFeatures_1.56.0 cli_3.6.3 locfit_1.5-9.10 compiler_4.4.1
Rsamtools_2.20.0\
[116] rlang_1.1.4 crayon_1.5.3 rngtools_1.5.2 labeling_0.4.3
nor1mix_1.3-3\
[121] mclust_6.1.1 affy_1.82.0 plyr_1.8.9 stringi_1.8.4
viridisLite_0.4.2\
[126] BiocParallel_1.38.0 ggalt_0.4.0 nleqslv_3.3.5 munsell_0.5.1
proj4_1.0-14\
[131] Matrix_1.7-0 hms_1.1.3 sparseMatrixStats_1.16.0 bit64_4.0.5
Rhdf5lib_1.26.0\
[136] KEGGREST_1.44.1 statmod_1.5.0 SummarizedExperiment_1.34.0
memoise_2.0.1 affyio_1.74.0\
[141] bit_4.0.5
