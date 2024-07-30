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
    (<https://wcaar.shinyapps.io/blood-biomarker-catalog/>), but does
    not contain data files used due to storage issues

    -   Files can be downloaded directly from our application on the
        Download tab

    -   Alternatively, the code in the misc directory and study_data
        folders can be used to produce the necessary data and results
        used in the application

-   As of July 21, 2024, this code is out of date and will be updated
    soon.

## Session Info

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

attached base packages: [1] grid stats4 stats graphics grDevices utils
datasets methods base

other attached packages: [1] biomaRt_2.60.0 VennDiagram_1.7.3
futile.logger_1.4.3\
[4] clusterProfiler_4.12.0 EnsDb.Hsapiens.v86_2.99.0 ggalt_0.4.0\
[7] RColorBrewer_1.1-3 EnhancedVolcano_1.22.0 PCAtools_2.16.0\
[10] ggrepel_0.9.5 lumi_2.56.0 illuminaHumanv3.db_1.26.0\
[13] illuminaio_0.46.0 beadarray_2.54.0 hexbin_1.28.3\
[16] hgu133plus2cdf_2.18.0 hgu133plus2.db_3.13.0 affycoretools_1.76.0\
[19] affy_1.82.0 oligo_1.68.2 Biostrings_2.72.1\
[22] XVector_0.44.0 oligoClasses_1.66.0 ensembldb_2.28.0\
[25] AnnotationFilter_1.28.0 pheatmap_1.0.12 heritability_1.4\
[28] MASS_7.3-61 variancePartition_1.34.0 BiocParallel_1.38.0\
[31] GEOquery_2.72.0 Homo.sapiens_1.3.1
TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 [34] org.Hs.eg.db_3.19.1
GO.db_3.19.1 OrganismDbi_1.46.0\
[37] GenomicFeatures_1.56.0 GenomicRanges_1.56.1 GenomeInfoDb_1.40.1\
[40] AnnotationDbi_1.66.0 IRanges_2.38.0 S4Vectors_0.42.0\
[43] Biobase_2.64.0 BiocGenerics_0.50.0 edgeR_4.2.0\
[46] limma_3.60.3 cowplot_1.1.3 jsonlite_1.8.8\
[49] httr_1.4.7 rJava_1.0-11 mailR_0.8\
[52] DT_0.33 lubridate_1.9.3 forcats_1.0.0\
[55] stringr_1.5.1 dplyr_1.1.4 purrr_1.0.2\
[58] readr_2.1.5 tidyr_1.3.1 tibble_3.2.1\
[61] ggplot2_3.5.1 tidyverse_2.0.0 bslib_0.7.0\
[64] shiny_1.8.1.1

loaded via a namespace (and not attached): [1] igraph_2.0.3 graph_1.82.0
minfi_1.50.0 Formula_1.2-5 maps_3.4.2\
[6] zlibbioc_1.50.0 tidyselect_1.2.1 bit_4.0.5 lattice_0.22-6
rjson_0.2.21\
[11] nor1mix_1.3-3 blob_1.2.4 rngtools_1.5.2 S4Arrays_1.4.1
pbkrtest_0.5.2\
[16] base64_2.0.1 parallel_4.4.1 dichromat_2.0-0.1 scrime_1.3.5
png_0.1-8\
[21] ggplotify_0.1.2 cli_3.6.3 GOstats_2.70.0 ProtGenerics_1.36.0
askpass_1.2.0\
[26] multtest_2.60.0 openssl_2.2.0 nleqslv_3.3.5 textshaping_0.4.0
BiocIO_1.14.0\
[31] shadowtext_0.1.3 curl_5.2.1 tidytree_0.4.6 mime_0.12
evaluate_0.24.0\
[36] stringi_1.8.4 backports_1.5.0 lmerTest_3.1-3 XML_3.99-0.16.1
httpuv_1.6.15\
[41] affxparser_1.76.0 magrittr_2.0.3 rappdirs_0.3.3 splines_4.4.1
mclust_6.1.1\
[46] doRNG_1.8.6 ggraph_2.2.1 ff_4.0.12 DBI_1.2.3 HDF5Array_1.32.0\
[51] jquerylib_0.1.4 genefilter_1.86.0 withr_3.0.0 corpcor_1.6.10
systemfonts_1.1.0\
[56] enrichplot_1.24.0 RBGL_1.80.0 GSEABase_1.66.0 tidygraph_1.3.1
formatR_1.14\
[61] rtracklayer_1.64.0 BiocManager_1.30.23 htmlwidgets_1.6.4 fs_1.6.4
labeling_0.4.3\
[66] fANCOVA_0.6-1 SparseArray_1.4.8 DESeq2_1.44.0 cellranger_1.1.0
MatrixGenerics_1.16.0\
[71] extrafont_0.19 annotate_1.82.0 VariantAnnotation_1.50.0
GGally_2.2.1 knitr_1.47\
[76] UCSC.utils_1.0.0 beanplot_1.3.1 RhpcBLASctl_0.23-42
timechange_0.3.0 foreach_1.5.2\
[81] fansi_1.0.6 Glimma_2.14.0 patchwork_1.2.0 caTools_1.18.2
ggtree_3.12.0\
[86] data.table_1.15.4 rhdf5_2.48.0 R.oo_1.26.0 irlba_2.3.5.1
extrafontdb_1.0\
[91] gridGraphics_0.5-1 lazyeval_0.2.2 yaml_2.3.8 survival_3.7-0
crayon_1.5.3\
[96] tweenr_2.0.3 later_1.3.2 Rgraphviz_2.48.0 codetools_0.2-20
base64enc_0.1-3\
[101] aod_1.3.3 KEGGREST_1.44.1 Rsamtools_2.20.0 filelock_1.0.3
foreign_0.8-86\
[106] pkgconfig_2.0.3 xml2_1.3.6 EnvStats_2.8.1 GenomicAlignments_1.40.0
aplot_0.2.3\
[111] ape_5.8 viridisLite_0.4.2 BSgenome_1.72.0 biovizBase_1.52.0
xtable_1.8-4\
[116] Category_2.70.0 hwriter_1.3.2.1 plyr_1.8.9 rbibutils_2.2.16
tools_4.4.1\
[121] htmlTable_2.4.2 AnnotationForge_1.46.0 broom_1.0.6 checkmate_2.3.1
nlme_3.1-165\
[126] lambda.r_1.2.4 HDO.db_0.99.1 dbplyr_2.5.0 crosstalk_1.2.1
assertthat_0.2.1\
[131] lme4_1.1-35.4 digest_0.6.36 numDeriv_2016.8-1.1 Matrix_1.7-0
farver_2.1.2\
[136] tzdb_0.4.0 remaCor_0.0.18 PFAM.db_3.19.1 reshape2_1.4.4
viridis_0.6.5\
[141] yulab.utils_0.1.4 rpart_4.1.23 glue_1.7.0 cachem_1.1.0
BiocFileCache_2.12.0\
[146] polyclip_1.10-6 methylumi_2.50.0 Hmisc_5.1-3 generics_0.1.3
mvtnorm_1.2-5\
[151] txdbmaker_1.0.0 statmod_1.5.0 ragg_1.3.2 ScaledMatrix_1.12.0
minqa_1.2.7\
[156] httr2_1.0.1 SummarizedExperiment_1.34.0 gson_0.1.0 ggbio_1.52.0
dqrng_0.4.1\
[161] ash_1.0-15 utf8_1.2.4 graphlayouts_1.1.1 siggenes_1.78.0
gtools_3.9.5\
[166] readxl_1.4.3 fontawesome_0.5.2 preprocessCore_1.66.0 affyio_1.74.0
gridExtra_2.3\
[171] GenomeInfoDbData_1.2.12 R.utils_2.12.3 rhdf5filters_1.16.0
RCurl_1.98-1.14 memoise_2.0.1\
[176] rmarkdown_2.27 scales_1.3.0 ReportingTools_2.44.0
R.methodsS3_1.8.2 reshape_0.8.9\
[181] rstudioapi_0.16.0 cluster_2.1.6 hms_1.1.3 munsell_0.5.1
colorspace_2.1-0\
[186] rlang_1.1.4 quadprog_1.5-8 DelayedMatrixStats_1.26.0
sparseMatrixStats_1.16.0 ggforce_0.4.2\
[191] mgcv_1.9-1 xfun_0.45 bumphunter_1.46.0 iterators_1.0.14
matrixStats_1.3.0\
[196] proj4_1.0-14 abind_1.4-5 GOSemSim_2.30.0 treeio_1.28.0
Rhdf5lib_1.26.0\
[201] futile.options_1.0.1 bitops_1.0-7 Rdpack_2.6 promises_1.3.0
scatterpie_0.2.3\
[206] RSQLite_2.3.7 qvalue_2.36.0 fgsea_1.30.0 DelayedArray_0.30.1
compiler_4.4.1\
[211] prettyunits_1.2.0 boot_1.3-30 beachmat_2.20.0 Rcpp_1.0.12
Rttf2pt1_1.3.12\
[216] BiocSingular_1.20.0 progress_1.2.3 R6_2.5.1 fastmap_1.2.0
fastmatch_1.1-4\
[221] ggstats_0.6.0 rsvd_1.0.5 nnet_7.3-19 gtable_0.3.5
KernSmooth_2.23-24\
[226] gcrma_2.76.0 htmltools_0.5.8.1 bit64_4.0.5 lifecycle_1.0.4
nloptr_2.1.0\
[231] restfulr_0.0.15 sass_0.4.9 vctrs_0.6.5 rsconnect_1.3.1
DOSE_3.30.1\
[236] ggfun_0.1.5 BeadDataPackR_1.56.0 pillar_1.9.0 gplots_3.1.3.1
locfit_1.5-9.10
