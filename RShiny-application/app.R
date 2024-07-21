# Load libraries ----
library(shiny)
library(bslib)
library(tidyverse)
library(DT)

# Load functions and objects ----
source("functions.R")
source("data.R")
e1 <- read_rds("./data/expr/gomez.rds")
e2 <- read_rds("./data/expr/Meaburn1.rds")
e3 <- read_rds("./data/expr/Meaburn2.rds")
e4 <- read_rds("./data/expr/gosch.rds")
e5 <- read_rds("./data/expr/obermoser1.rds")
e6 <- read_rds("./data/expr/obermoser2.rds")
e7 <- read_rds("./data/expr/obermoser3.rds")
e8 <- read_rds("./data/expr/obermoser4.rds")
e9 <- read_rds("./data/expr/dusek.rds")
e10 <- read_rds("./data/expr/rusch.rds")
e11 <- read_rds("./data/expr/larocca.rds")
symbolOptions <- read.table("./data/symbol.txt")
symbolOptions <- symbolOptions$x
eqtlOptionsGene <- unique(gtex$`Symbol`)
eqtlOptionsTrait <- unique(gtex$`GWAS Trait`)

tutorial <- function(text) {
  tooltip(icon("comments"), paste(text), placement = "top")
}
info <- function(text) {
  tooltip(icon("circle-info"), paste(text), placement = "top")
}

# Build ui ----
ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "litera"),
  titlePanel(title = "", windowTitle = "Blood Biomarker Catalog"),
  # Define the main navigation tabs
  navbarPage(
    title = "Blood Biomarker Catalog",
    
    ## Home page----
    tabPanel("Home",
             navlistPanel(widths = c(3,9), well = FALSE,
                          ### About----
                          tabPanel(title = "About",
                                   tags$h5(tags$em("Graphical abstract")),
                                   img(src = "biomarker-guide.png", height = "482px"),
                                   tags$h3("Purpose of this application"),
                                   tags$p("The purpose of this application is to find suitable blood RNA biomarkers for clinical research.
                               The application is formatted so that researchers can"),
                                   tags$ul(
                                     tags$li("easily find RNAs associated with their clinical trait of interest"),
                                     tags$li("view the genes, variants, and expression quantitative trait loci (eQTLs) associated with a trait"),
                                     tags$li("see differences in whole blood gene expression based on time or individual.")
                                   ),
                                   tags$p("This application features an analysis on 814 whole-blood transcriptomes collected from 141 healthy individuals sampled at multiple timepoints,
                               from 50 minutes up to 4 months. The results give key insights into the natural variation and flucuation that occurs in whole blood RNA.",
                                          tags$b("To get started, click on the \"Tutorial\" tab to the left.")),
                                   tags$p("The suggested workflow for this application focuses on cataloging transcripts by which of the seven biomarker categories it would be best suited for,
                        as defined by the NIH BEST resource (see 'Acknowledgements'). We hypothesized that transcripts that maintained relatively constant expression 
                        over time but varied between individuals are under genetic regulatory control and are resistant to environment perturbation. This suggests
                        that biomarkers with temporal stability and interindividual variation would be strong diagnostic, predictive, prognostic, and risk biomarkers."),
                                   tags$p("We also reasoned that blood transcripts that are differentially expressed over time are dynamically responding to changes in the environment. Indeed,
                               these genes tended to be enriched in homeostatic regulatory processes, which would make their transcripts ideal candidates for safety biomarkers (monitoring
                               exposure to harmful environments) and for pharmacodynamic response biomarkers."),
                                   
                                   tags$h3("About the Data"),
                                   tags$p("The datasets included in our study derive from seven independent research studies, each analyzing whole blood RNA, 
                focusing on different time intervals. These datasets provide a comprehensive overview of how gene expression varies 
                over time and between individuals under normal health conditions."),
                                   
                                   accordion(
                                     accordion_panel(title = "Gomez-Carballa et al.",
                                                     tags$h5("About this dataset"),
                                                     tags$p("This study observed changes in the whole-blood transcriptome of heatlhy and cognitively declining
                                                 individuals before and after a 50 minute concert. The researchers observed whether any changes in 
                                                 cognitive behavior marker genes occured as a result of musical intervention. They obtained blood via capillary puncture."),
                                                     tags$p("For the purposes of our study, we selected only individuals defined as healthy by the study authors."),
                                                     tags$br(),
                                                     tags$i("Example entries"),
                                                     DT::dataTableOutput("gomezTableHead")
                                     ),
                                     accordion_panel(title = "Meaburn et al.",
                                                     tags$h5("About this dataset"),
                                                     tags$p("This study collected whole blood from sets of 12-year-old twins four hours apart on the same day and 10 months later.
                                                 The purpose of this study was to test how repeatable and reliable microarray measurements are for the same individual.
                                                 We used all data from this dataset that was complete and met our quality control standards. They collected blood via venipuncture."),
                                                     tags$br(),
                                                     tags$i("Example entries"),
                                                     DT::dataTableOutput("meaburnTableHead")
                                     ),
                                     accordion_panel(title = "Gosch et al.",
                                                     tags$h5("About this dataset"),
                                                     tags$p("This study observed changes in the whole-blood transcriptome of heatlhy adults under the age of 30
                                                 every three hours for a 24-hour period. The subjects did not sleep throughout the period, but were
                                                 permitted to eat and drink normally. The researchers of this study were interested in modeling blood
                                                 gene expression for forensic purposes. All samples from this study were used for analysis. They colleceted blood via capillary puncture."),
                                                     tags$br(),
                                                     tags$i("Example entries"),
                                                     DT::dataTableOutput("goschTableHead")
                                     ),
                                     accordion_panel(title = "Obermoser et al.",
                                                     tags$br(),
                                                     tags$h5("About this dataset"),
                                                     tags$p("This study observed changes in the whole-blood transcriptome of heatlhy adults over the course of 5 weeks,
                                                 with multiple time intervals. Subjects were enrolled in two cohorts and given one of three treatments:
                                                 a flu vaccine (FLUZONE), a pneumonia vaccine (PNEUMOVAX), or a saline control. The collected blood via capillary puncture and venipuncture."),
                                                     tags$p("We included all subjects from the study who had at least one sample from each of the timepoints the researchers measured."),
                                                     tags$i("Example entries"),
                                                     DT::dataTableOutput("obermoserTableHead1"),
                                                     DT::dataTableOutput("obermoserTableHead2")
                                                     
                                     ),
                                     accordion_panel(title = "Dusek et al.",
                                                     tags$h5("About this dataset"),
                                                     tags$p("This study observed changes in the whole-blood transcriptome of heatlhy adults before and after 8 weeks.
                                                 The subjects had no prior experience in meditation and blood collected at baseline. After 8 weeks of self-guided meditation exercises
                                                 (designed to induce the relaxation response), blood was collected again for comparison. All subjects with comeplete data
                                                 were used for our study."),
                                                     tags$br(),
                                                     tags$i("Example entries"),
                                                     DT::dataTableOutput("DusekTableHead")
                                     ),
                                     accordion_panel(title = "Rusch et al.",
                                                     tags$h5("About this dataset"),
                                                     tags$p("This study observed changes in the whole-blood transcriptome of military-age men before and after 3 months.
                                                 The researchers were interested in comparing non-PTSD blood transcriptomes to those of PTSD servicemembers. We
                                                 only used data from servicemen who the researchers determined did not have PTSD and consider these individuals as healthy.
                                                            They collected blood via venipuncture."),
                                                     tags$br(),
                                                     tags$i("Example entries"),
                                                     DT::dataTableOutput("ruschTableHead"),
                                     ),
                                     accordion_panel(title = "LaRocca et al.",
                                                     tags$h5("About this dataset"),
                                                     tags$p("This study observed changes in the whole-blood transcriptome of heatlhy women before and after being assigned
                                                 a 4-month exercise regimen. The intensity and duration of the regimen was randomized. The researchers were interested
                                                 in observing what differences in the transcriptome were observed between women who experienced a change in VO2 max (reponders, R)
                                                 and those who did not respond to exercise (non-responders, NR). We used all samples from this dataset for our analysis."),
                                                     tags$br(),
                                                     tags$i("Example entries"),
                                                     DT::dataTableOutput("laroccaTableHead")
                                     )
                                   )),
                          ### Tutorial----
                          tabPanel(title = "Tutorial",
                                   tags$h3("Step-by-Step Guide to Using This Application"),
                                   tags$p("This tutorial is a step-by-step guide for investigating RNA biomarkers of interest. We use gene symbols when referring to both genes and their transcripts in whole blood."),
                                   tags$p("The application is structured so that you may start with the biological trait of interest, identify the RNAs that correspond with that trait, and filtering those RNAs according to
                               the biomarker type you are interested in (go back to the About section for a brief summary on what this means)."),
                                   tags$p("We will do a worked example from start to finish. The", tutorial("Example text"), "icon indicates optional commentary which can be viewed by hovering your mouse over the icon. At any time you use the application, we would appreciate any feedback you have and encourage you to submit it via the Submit Feedback tab. Your message will immediately
                               be sent to the developers!"),
                                   
                                   tags$h4("Step 1: Finding RNAs that correspond to your trait of interest"),
                                   tags$p("The first tab in the navigation bar is called the eQTL Variant Search. This tab consists of a synthesis of the GWAS Catalog, GTEx database, and healthy, longitudinal whole-blood transcriptome data for candidate biomarker identification."),
                                   tags$p("When you open the page, you are immediately brought to a table consisting of variant data on 9,474 unique genes. On the left side are the filtering tools which will help us find RNAs we are interested in."),
                                   tags$ul(
                                     tags$li("the 'Gene' button filters the table to contain only results for the gene in the dropdown menu. You can adjust the gene in the box by clicking through the options or searching your own symbol by typing it in.", tutorial("When you select 'Gene' as your filtering option, the sliders will not work. This is because scores are assigned to gene symbols, so all results will have the same stable-polymorphic, flexibility, and housekeeping score.")),
                                     tags$li("the 'Disease/Trait' button filters the table to contain only genes which have been associated with the selected trait in the GWAS Catalog."),
                                     tags$li("the 'No filter' button does not filter by trait or gene, only by the slider values selected below.")
                                   ),
                                   div(img(src = "tutorial1.png", height = "450px")),
                                   tags$h5("Choosing your genes of interest"),
                                   tags$p("For researchers who already have genes they are interested in studying, we suggest using the 'Gene' filtering option. Otherwise, the disease/trait filter option can help with discovering potentially new biomarkers for clinical validation. The slider options below these buttons filter each gene according to their expression profile in whole blood RNA. When the slider is set to -1, it is turned off, and the table will not be filtered by that score. These sliders will help us find candidate biomarker genes."),
                                   tags$ul(
                                     tags$li(tags$b("Stable-polymorphic score:"), "high score indicates stable expression over time, variable expression over time (see top right panel in Graphical Abstract). Ideal candidates for diagnostic, predictive, prognostic, and risk biomarkers."),
                                     tags$li(tags$b("Flexibility score:"), "high score indicates differential expression over time. Ideal candidates for safety and pharmacodynamic response biomarkers."),
                                     tags$li(tags$b("Housekeeping score:"), "high score indicates above median expression level, stable over time, and low variation across individuals. Since our samples come from healthy individuals, we recommend using these RNAs to create healthy reference ranges for a gene.", tutorial("You can view the expression level of genes in the 'Gene Symbol Search' panel, which we discuss later on in the tutorial."))
                                   ),
                                   tags$p("For this example, let's look for RNAs associated with cardiovascular disease. Ensure the 'Disease/Trait' button is ticked and search for 'cardiovascular' in the search box. To load results, click on the trait of interest."),
                                   div(img(src = "tutorial2.png", height = "300px"), img(src = "tutorial3.png", height = "300px")),
                                   tags$p("Let's find RNAs that are stable over time and don't vary between individuals. These will be housekeeping genes, and we will use our example to find healthy reference ranges for these RNAs in whole blood. We can find the houseekping score slider at the bottom of the search panel. Set the slider to 20 and click the '>='. button to filter for all genes with a housekeeping score greater than or equal to 20."),
                                   div(img(src = "tutorial4.png", height = "300px"), img(src = "tutorial5.png", height = "300px"), tutorial("Because the table is very large, you must slide all the way right to view the scores. As you can see, all housekeeping scores for these genes are greater than or equal to 20.")),
                                   tags$p("Some genes like ", tags$code("ERAP2") ," are stable over time, but variable between individuals. We can keep adjusting our filters until we find the most housekeeping-like genes possible.", tutorial("There are many ways to do this: raising the housekeeping threshold, or minimizing the stable-polymorphic / flexibility score. You can experiment with these filters and get different results each time.")),
                                   div(img(src = "tutorial6.png", height = "450px"), tutorial("Hovering over the column names shows an HTML title that explains each variable. You may need to hold your mouse there for a few seconds- if it's still not loading, try moving between columns until the app recongizes you are searching for titles.")),
                                   tags$p("We've narrowed our options down significantly. In GTEx, ", tags$code("rs10858023"), " is associated with expression levels of ", tags$code("AP4B1"),". In GWAS, ", tags$code("rs10858023"), " maps to an intron variant of ", tags$code("DCLRE1B"), ". Both the eQTL slope and GWAS association are significant."),
                                   tags$p("From here, we can use the 'Gene' filter button to see what other eQTLs and GWAS traits are associated with ", tags$code("AP4B1") ," or ", tags$code("DCLRE1B"),"."),
                                   div(img(src = "tutorial7.png", height = "300px"), img(src = "tutorial8.png", height = "300px"), tutorial("The top right of the table is a plotting feature that shows the number of entries a Gene or Trait has in the table.")),
                                   tags$p("We see that, despite being relatively stable in blood for healthy people", tutorial("AP4B1 was found to be a housekeeping gene in 16 normal tissues: see Eisenberg, Eli, and Erez Y. Levanon. “Human Housekeeping Genes, Revisited.” Trends in Genetics, Human Genetics, 29, no. 10 (October 1, 2013): 569–74. https://doi.org/10.1016/j.tig.2013.05.010.
") ,", ", tags$code("AP4B1"), " has many eQTL variants that are also associated with GWAS traits. We can explore ", tags$code("AP4B1"), "in depth using the", tags$b("Gene Symbol Search tab.")),
                                   
                                   tags$h4("Step 2: Analyzing your RNAs of interest in healthy blood"),
                                   tags$p("To explore the expression levels and summary statistics for your gene in our data, you’ll want to click the ", tags$b("Gene Symbol Search"), " tab in the navigation bar."),
                                   div(img(src = "tutorial9.png", height = "450px")),
                                   tags$em("This is the screen you will see upon clicking the tab."),
                                   tags$p("Let's learn more about our housekeeping gene,", tags$code("AP4B1"), ". We click the top left box, search the gene, and click it."),
                                   div(img(src = "tutorial10.png", height = "450px")),
                                   tags$p("The 'Select dataset' box allows us to choose from one of the seven studies described in the 'About' section. This will allow us to view our gene's expression level on different platforms (microarray or RNA-seq) and over differnt lengths of time (50 minutes to 16 weeks). Let's select the Rusch et al. dataset, and let's click the 'Average Expression' row on the table in the top left box. You should see this screen:"),
                                   div(img(src = "tutorial11.png", height = "450px")),
                                   tags$p("Let's explore what's on this screen:"),
                                   tags$ul(
                                     tags$li(tags$b("Top left:"), "the statistics we calculated for ", tags$code("AP4B1"), " in the study we selected. These statistics vary depending on the study you selected and whether the gene symbol of interest is represented by multiple probes (explained below). You can click on any statistic to visualize it in the bottom left box."),
                                     tags$li(tags$b("Top right:"), "a spaghetti plot of gene expression values for ", tags$code("AP4B1"), ", specifically probe ", tags$code("231714_s_at") ,". You can select from multiple probes on a microarray using the 'Select Probe' tool on the left. Each colored line represents a different individual from the study selected in the sidebar."),
                                     tags$li(tags$b("Bottom left:"), "a panel with 5 tabs to visualize your results. The stat density tab displays a density plot for the statistic of interest; the stat boxplot disaplys a boxplot for the statistic. The gene or probe you select will be represented by a blue line. The last three tabs are the scores from earlier and show where ", tags$code("AP4B1"), " lies among all other genes.", tutorial("If you hover over the question marks on these last three tabs, you will get an explanation of what each score signifies as a brief reminder.")),
                                     tags$li(tags$b("Bottom right:"), "a summary of the gene information and transcripts from the Ensembl Genome Browser. Some transcripts / genes may not have information in this box.")
                                   ),
                                   tags$p("We see on the bottom right that our gene is ", tags$code("adaptor related protein complex 4 subunit beta 1"), " and, as we expected, it is robustly expressed in blood, remains very stable over 12 weeks, and has modest variability between individuals.", tutorial("While we don't expect housekeeping genes to vary between individuals, the slight variation in AP4B1 could be explained by its numerous eQTLs.")),
                                   tags$p("The last option on this screen is adjusting the coloring of the spaghetti plot. By typing directly into the box, you can change how the lines are colored according to variables recorded by the study. Each dataset has different coloring options- in order to see all these options, go to 'Home', 'About' and explore each of the datasets. The table headers are the available options.", tutorial("To switch back to individual coloring, type 'subject' into the box. All studies have a 'subject' option.")),
                                   div(img(src = "tutorial12.png", height = "300px"),img(src = "tutorial13.png", height = "300px")),
                                   tags$p("Some of the interesting color options are 'treat' in the Obermoser dataset, which shows responses to vaccine or saline injections, and 'group' in the LaRocca dataset, which shows whether or not an increase in VO2 max was observed after 4 months of training."),
                                   tags$h4("Step 3: Collecting your biomarkers"),
                                   tags$p("This tool helps facilitate strong biomarker candidates for your trait of interest. You can search by condition to find new genes, explore genes you are interested in, or just observe patterns in healthy whole blood generally. We hope that these features allow you to make the best decisions for your whole-blood clincial biomarker research."),
                                   tags$p("For ", tags$code("AP4B1"), ", we can note the Average Expression as well as the Total Variation to make informed decisions about references ranges for this biomarker, and use that information to construct well-designed experiments for clinical validation of blood RNA biomarkers."),
                                   tags$p("Finally, all summary statistics for each dataset and scores", tutorial("Stable-polymorphic, flexibility, and housekeeping scores"), "are available in .csv format on the ", tags$b("Download Data"), " tab."),
                                   tags$h5("Conclusion and practical guidance"),
                                   tags$p("This concludes the tutorial. We recommend starting by investigating a disease or trait you are interested in using the eQTL variant search tab, and then filtering your results to find a suitable biomarker candidate for that disease. Then, search for that biomarker candidate in the Gene Symbol Search tab to find statistics for that gene and see how its expression levels vary from study to study."),
                                   tags$p("Thank you for using our application. We appreciate any feedback! Submit suggestions or comments using the 'Submit Feedback' tool on the left side of this page (you may need to scroll up to see it)."),
                                   tags$b("References can be found in the 'Acknowledgements' section")
                          ),
                          ### Acknowledgements----
                          tabPanel(title = "Acknowledgments",
                                   tags$h3("Credits"),
                                   tags$p("Designed by", tags$b("Will C Baltazar,"), "idea conceived by", tags$b("Laura B Ferguson,"), "supported and tested by", tags$b("The Messing Lab.")),
                                   tags$p("The Messing Lab"),
                                   tags$p("The Waggoner Center for Alcohol and Addiction Research"),
                                   tags$p("Department of Neurology"),
                                   tags$p("Department of Neuroscience"),
                                   tags$p("Dell Medical School"),
                                   tags$p("The University of Texas at Austin"),
                                   div(img(src = "wagg.jpeg", height = "40px"), img(src = "RGB_formal_Dell.jpg", height = "40px")),
                                   br(),
                                   tags$h3("References"),
                                   tags$p("Dusek, Jeffery A., Hasan H. Otu, Ann L. Wohlhueter, Manoj Dusek, Luiz F. Zerbini, Marie G. Joseph, Herbert Benson, and Towia A. Libermann. “Genomic Counter-Stress Changes Induced by the Relaxation Response.” PLoS ONE 3, no. 7 (July 2, 2008): e2576. https://doi.org/10.1371/journal.pone.0002576."),
                                   tags$p("Fergal J Martin, M Ridwan Amode, Alisha Aneja, Olanrewaju Austine-Orimoloye, Andrey G Azov, If Barnes, Arne Becker, Ruth Bennett, Andrew Berry, Jyothish Bhai, Simarpreet Kaur Bhurji, Alexandra Bignell, Sanjay Boddu, Paulo R Branco Lins, Lucy Brooks, Shashank Budhanuru Ramaraju, Mehrnaz Charkhchi, Alexander Cockburn, Luca Da Rin Fiorretto, Claire Davidson, Kamalkumar Dodiya, Sarah Donaldson, Bilal El Houdaigui, Tamara El Naboulsi, Reham Fatima, Carlos Garcia Giron, Thiago Genez, Gurpreet S Ghattaoraya, Jose Gonzalez Martinez, Cristi Guijarro, Matthew Hardy, Zoe Hollis, Thibaut Hourlier, Toby Hunt, Mike Kay, Vinay Kaykala, Tuan Le, Diana Lemos, Diego Marques-Coelho, José Carlos Marugán, Gabriela Alejandra Merino, Louisse Paola Mirabueno, Aleena Mushtaq, Syed Nakib Hossain, Denye N Ogeh, Manoj Pandian Sakthivel, Anne Parker, Malcolm Perry, Ivana Piližota, Irina Prosovetskaia, José G Pérez-Silva, Ahamed Imran Abdul Salam, Nuno Saraiva-Agostinho, Helen Schuilenburg, Dan Sheppard, Swati Sinha, Botond Sipos, William Stark, Emily Steed, Ranjit Sukumaran, Dulika Sumathipala, Marie-Marthe Suner, Likhitha Surapaneni, Kyösti Sutinen, Michal Szpak, Francesca Floriana Tricomi, David Urbina-Gómez, Andres Veidenberg, Thomas A Walsh, Brandon Walts, Elizabeth Wass, Natalie Willhoft, Jamie Allen, Jorge Alvarez-Jarreta, Marc Chakiachvili, Bethany Flint, Stefano Giorgetti, Leanne Haggerty, Garth R Ilsley, Jane E Loveland, Benjamin Moore, Jonathan M Mudge, John Tate, David Thybert, Stephen J Trevanion, Andrea Winterbottom, Adam Frankish, Sarah E Hunt, Magali Ruffier, Fiona Cunningham, Sarah Dyer, Robert D Finn, Kevin L Howe, Peter W Harrison, Andrew D Yates, and Paul Flicek
Ensembl 2023
Nucleic Acids Res. 2023, 51(D1):D933-D941
PMID: 36318249
doi:10.1093/nar/gkac958"),
                                   tags$p("Gómez-Carballa, A., Navarro, L., Pardo-Seco, J., Bello, X., Pischedda, S., Viz-Lasheras, S., Camino-Mera, A., et al. (2023). Music compensates for altered gene expression in age-related cognitive disorders. Scientific Reports, 13(1), 21259. https://doi.org/10.1038/s41598-023-48094-5"),
                                   tags$p("Gosch, A., Bhardwaj, A., & Courts, C. (2023). TrACES of time: Transcriptomic analyses for the contextualization of evidential stains – Identification of RNA markers for estimating time-of-day of bloodstain deposition. Forensic Science International: Genetics, 67, 102915. https://doi.org/10.1016/j.fsigen.2023.102915"),
                                   tags$p("LaRocca, T. J., Smith, M. E., Freeberg, K. A., Craighead, D. H., Helmuth, T., Robinson, M. M., Nair, K. S., Bryan, A. D., & Seals, D. R. (2023). Novel whole blood transcriptome signatures of changes in maximal aerobic capacity in response to endurance exercise training in healthy women. Physiological Genomics, 55(8), 338–344. https://doi.org/10.1152/physiolgenomics.00017.2023"),
                                   tags$p("Meaburn, E. L., Fernandes, C., Craig, I. W., Plomin, R., & Schalkwyk, L. C. (2009). Assessing individual differences in genome-wide gene expression in human whole blood: Reliability over four hours and stability over 10 months. Twin Research and Human Genetics: The Official Journal of the International Society for Twin Studies, 12(4), 372. https://doi.org/10.1375/twin.12.4.372"),
                                   tags$p("Obermoser, G., Presnell, S., Domico, K., Xu, H., Wang, Y., Anguiano, E., Thompson-Snipes, L., et al. (2013). Systems scale interactive exploration reveals quantitative and qualitative differences in response to influenza and pneumococcal vaccines. Immunity, 38(4), 831–844. https://doi.org/10.1016/j.immuni.2012.12.008"),
                                   tags$p("Rusch, H. L., Robinson, J., Yun, S., Osier, N. D., Martin, C., Brewin, C. R., & Gill, J. M. (2019). Gene expression differences in PTSD are uniquely related to the intrusion symptom cluster: A transcriptome-wide analysis in military service members. Brain, Behavior, and Immunity, 80, 904–908. https://doi.org/10.1016/j.bbi.2019.04.039"),
                                   tags$p("Sollis E, Mosaku A, Abid A, Buniello A, Cerezo M, Gil L, Groza T, Güneş O, Hall P, Hayhurst J, Ibrahim A, Ji Y, John S, Lewis E, MacArthur JAL, McMahon A, Osumi-Sutherland D, Panoutsopoulou K, Pendlington Z, Ramachandran S, Stefancsik R, Stewart J, Whetzel P, Wilson R, Hindorff L, Cunningham F, Lambert SA, Inouye M, Parkinson H, Harris LW.
The NHGRI-EBI GWAS Catalog: knowledgebase and deposition resource.
Nucleic Acids Res. 2022 Nov 9:gkac1010. doi: 10.1093/nar/gkac1010. Epub ahead of print. PMID: 36350656."),
                                   tags$p("The data used for the analyses described in this presentation were obtained from:", tags$link("https://gtexportal.org/home/downloads/adult-gtex/qtl") ,"the GTEx Portal on 05/04/24."),
                                   tags$p("FDA-NIH Biomarker Working Group. BEST (Biomarkers, EndpointS, and other Tools) Resource [Internet]. Silver Spring (MD): Food and Drug Administration (US); 2016-. Available from: https://www.ncbi.nlm.nih.gov/books/NBK326791/ Co-published by National Institutes of Health (US), Bethesda (MD).")
                          ),
                          ### Feedback----
                          tabPanel(title = "Submit Feedback",
                                   card(card_title("Submit your feedback!"),
                                        card_body(textInput(inputId = "feedbackName", "Email address"),
                                                  textInput(inputId = "feedbackTitle", "Subject"),
                                                  textInput(inputId = "feedback", "Please insert your feedback here.", width = "600px"),
                                                  actionButton(inputId = "submitFeedback", label = "Submit", width = "100px")
                                        )),
                                   textOutput(outputId = "feedbackResult"))
             )
    ),
    
    ## eQTL Search page----
    tabPanel(title = "eQTL Variant Search",
             ### Search options----
             sidebarLayout(
               sidebarPanel = sidebarPanel(width = 3,
                                           radioButtons(inputId = "gtexFilter", label = "Filter by", 
                                                        choices = c("Gene", "Disease/Trait","No filter"), selected = "No filter"),
                                           selectizeInput(inputId = "eqtlChoiceGene",
                                                          selected = NULL,
                                                          label = "Symbol to filter by",
                                                          choices = NULL),
                                           selectizeInput(inputId = "eqtlChoiceTrait",
                                                          selected = NULL,
                                                          label = "Disease or trait to filter by",
                                                          choices = NULL),
                                           sliderInput(inputId = "stableTableSlider",
                                                       label = "Stable-polymorphic score?",
                                                       min = -1, max = 22, value = -1),
                                           radioButtons(inputId = "stableInequality", label = "only show studies...",
                                                        choiceNames = c("=","<=",">="), 
                                                        choiceValues = c("equal","less","more")),
                                           sliderInput(inputId = "dynamicTableSlider",
                                                       label = "Number of DE studies?",
                                                       min = -1, max = 6, value = -1),
                                           radioButtons(inputId = "dynamicInequality", label = "only show studies...",
                                                        choiceNames = c("=","<=",">="), 
                                                        choiceValues = c("equal","less","more")),
                                           sliderInput(inputId = "houseSlider",
                                                       label = "Housekeeping score?",
                                                       min = -1, max = 39, value = -1),
                                           radioButtons(inputId = "houseInequality", label = "only show studies...",
                                                        choiceNames = c("=","<=",">="), 
                                                        choiceValues = c("equal","less","more")),
                                           tags$i("Note: setting sliders to -1 turns them off.", style = "font-size:12px"),
               ),
               
               ### Search output----
               mainPanel = mainPanel(width = 9,
                                     page_fillable(
                                       navset_card_underline(full_screen = T, height = "800px", title = "Search output",
                                                             nav_panel("Table of variants and associations",
                                                                       DT::dataTableOutput(outputId = "gtexTable")),
                                                             nav_panel("Plot number of variants in GTEx v8",
                                                                       plotOutput(outputId = "gtexGraph"),
                                                                       div(style = "display:inline-block",
                                                                           sliderInput(inputId = "gtexGraphSlider", label = "Number of elements on x-axis",
                                                                                       min = 4, max = 500, value = 100, step = 2, ticks = F),
                                                                           sliderInput(inputId = "gtexGraphSlider2", label = "Text size",
                                                                                       min = 3, max = 15, value = 10))
                                                             )
                                       )
                                       
                                     )
               )
             )
    ),
    
    ## Search for gene symbol----
    tabPanel("Gene Symbol Search",
             ### Gene symbol and dataset search options----
             sidebarLayout(
               sidebarPanel(
                 selectizeInput(inputId = "geneName", label = "Select Gene Symbol:", choices = NULL),
                 selectInput("pData", "Select dataset:",
                             choices = c("Gomez-Carballa et al. (50 minutes, finger, listened to music, RNA-seq)" = "pheno_gomez",
                                         "Meaburn et al. (4 hours/day 1, vein, 12 year olds, microarray)" = "pheno_meaburn1",
                                         "Meaburn et al. (4 hours/day 2, vein, 12 year olds, microarray)" = "pheno_meaburn2",
                                         "Gosch et al. (24 hours, finger, every 3 hrs w/o sleep, RNA-seq)" = "pheno_gosch",
                                         "Obermoser et al. (5 weeks, vein, vaccine cohort 1, microarray)" = "pheno_obermoser1",
                                         "Obermoser et al. (2 weeks, finger, vaccine cohort 2, microarray)" = "pheno_obermoser2",
                                         "Obermoser et al. (5 weeks, vein, vaccine cohort 2, microarray)" = "pheno_obermoser3",
                                         "Obermoser et al. (9 days, finger, vaccine cohort 2, microarray)" = "pheno_obermoser4",
                                         "Dusek et al. (8 weeks, relaxation response training, microarray)" = "pheno_dusek",
                                         "Rusch et al. (12 weeks, vein, males in military w/ sleeping disturbances, microarray)" = "pheno_rusch",
                                         "LaRocca et al. (16 weeks, women in endurance training, RNA-seq)" = "pheno_larocca"),
                             selected = "Gosch et al. (24 hours)"),
                 textInput("column", "Change Coloring (optional):", "subject"),
                 # uiOutput(outputId = "columnSelectUI"), # this will replace the textInput on the line above.
                 uiOutput(outputId = "probeSelectUI"),
                 card(card_title("Statistics legend"),
                      card_body(
                        tags$b(tags$p("Within Variation:")), tags$p("average standard deviation of a transcript's expression levels over time for a single individual"),
                        tags$b(tags$p("Total Variation:")), tags$p("standard deviation of a transcript's expression level across all individuals in the dataset"),
                        tags$b(tags$p("Rs:")), tags$p("within variation divided by total variation"),
                        tags$b(tags$p("Average Expression")), tags$p("the mean expression level of the transcript after normalization"),
                        tags$b(tags$p("Repeatability:")), tags$p("the upper bound of broad-sense heritability, H^2; or, the percent of variation in repeated measurements explained by individual and not environment determined using the repeatability function from the heritability R package"),
                        tags$b(tags$p("Genetic Variance:")), tags$p("quantifies the distance of expression measurements between individuals, or how much of the variation is explained by genetic influences. Determined using the repeatability function"),
                        tags$b(tags$p("Residual Variance:")), tags$p("quantifies the residual (unexplained) variance of expression measurements from an individual using the repeatability function"),
                        tags$b(tags$p("Subject VP:")), tags$p("proportion of variance explained by individual alone using the variancePartition package"),
                        tags$b(tags$p("Time VP:")), tags$p("proportion of variance explained by timepoint alone using variancePartition"),
                        tags$b(tags$p("(variable) VP:")), tags$p("proportion of variance in the gene's expression explained by (vaiable) using variancePartition"),
                        tags$b(tags$p("Residual VP:")), tags$p("proportion of residual (unexplained) variance using variancePartition")
                      )
                 )
               ),
               ### Gene symbol search output----
               mainPanel(
                 layout_columns(
                   card(card_title("Statistics"),
                        card_body(
                          DT::dataTableOutput(outputId = "vdata")
                        )),
                   card(
                     card_body(
                       plotOutput("genePlot")
                       
                     )
                   )
                 ),
                 layout_columns(
                   navset_card_tab(nav_panel(title = span("Stat Density", uiOutput("noStatSelected1", inline = TRUE)),
                                             plotOutput(outputId = "statDensityPlot")),
                                   nav_panel(title = span("Stat Boxplot", uiOutput("noStatSelected2", inline = TRUE)),
                                             plotOutput(outputId = "statBoxplot")),
                                   nav_panel(title = span(
                                     "Stable-polymorphic score",
                                     tooltip(
                                       icon("circle-question"),
                                       "The stable-polymorphic score is a measure of how stable a gene is within an individual over time and variable between individuals. It has a maximum value of 22 (ERAP2 is the only gene with this score)",
                                       placement = "right"
                                     )
                                   ),
                                   plotOutput("stablePolyScorePlot"),
                                   tags$b(tags$p("Number of stable-polymorphic thresholds passed:")),
                                   tags$p(textOutput("stableNum")), 
                                   tags$i(textOutput("stableStudies")),
                                   ),
                                   nav_panel(title = span(
                                     "Flexibility score",
                                     tooltip(
                                       icon("circle-question"),
                                       "The flexibility score is how many studies in which a gene was called differentially expressed with time. The maximum score is a 6.",
                                       placement = "right"
                                     )
                                   ),
                                   plotOutput("flexibilityScorePlot"),
                                   tags$b(tags$p("Number of studies where this gene was called differntially expressed over time:")),
                                   tags$p(textOutput("dynamicNum")), 
                                   tags$i(textOutput("dynamicStudies")),
                                   ),
                                   nav_panel(title = span(
                                     "Housekeeping score",
                                     tooltip(
                                       icon("circle-question"),
                                       "The housekeeping score is a measure of how stable a gene is within an individual, the same across individuals, and highly expressed. It has a maximum value of 39 (RANBP3 is the only gene with this score)",
                                       placement = "right"
                                     )
                                   ),
                                   plotOutput("housekeepingScorePlot"),
                                   tags$b(tags$p("Number of housekeeping thresholds passed:")),
                                   tags$p(textOutput("houseNum")), 
                                   tags$i(textOutput("houseStudies")),
                                   ),
                   ),
                   card(card_title("Gene/transcript information"),
                        card_body(
                          dataTableOutput(outputId = "geneInfo")
                        )
                   )
                 )
               )
             )
    ),
    
    ## Download page ----
    tabPanel("Download Data",
             sidebarLayout(
               sidebarPanel = sidebarPanel(
                 tags$h3("Download data"),
                 selectInput(inputId = "downloadSelectInput", label = "Choose CSV file", 
                             choices = c("Gomez 50 minute transcript statistics (3.4 MB)" = "gomez",
                                         "Meaburn 4 hour transcript statistics, day 1 (9.7 MB)" = "m1",
                                         "Meaburn 4 hour transcript statistics, day 2 (10.2 MB)" = "m2",
                                         "Gosch 24 hour transcript statistics (3.2 MB)" = "gosch",
                                         "Obermoser C1, vein 5 week transcript statistics (5.1 MB)" = "o1",
                                         "Obermoser C2, finger 2 week transcript statistics (5.9 MB)" = "o2",
                                         "Obermoser C2, vein 5 week transcript statistics (6.5 MB)" = "o3",
                                         "Obermoser C2, finger 9 day transcript statistics (5.2 MB)" = "o4",
                                         "Dusek 8 week transcript statistics (10.5 MB)" = "dusek",
                                         "Rusch 12 week transcript statistics (12.3 MB)" = "rusch",
                                         "LaRocca 16 week transcript statistics (2.8 MB)" = "larocca",
                                         "Stable-polymorphic scores (1.3 MB)" = "stable",
                                         "Flexibility scores (4.6 MB)" = "flex",
                                         "Housekeeping scores (3.4 MB)" = "house"), selected = ""),
                 downloadButton("downloadData", "Download")
               ),
               mainPanel = mainPanel(
                 dataTableOutput("downloadSelectPreview")
               )
             )
    ),
    tags$script(HTML("var header = $('.navbar > .container-fluid');
header.append('<div style=\"float:right\"><img src=\"Logo.png\" alt=\"alt\" style=\"float:right;width:50px;height:50px;padding-top:0px;\"> </a></div>');
    console.log(header)")
    )
  )
)


# Define the server logic ----
server <- function(input, output, session) {
  
  ## Search bar options ----
  updateSelectizeInput(session, inputId = "geneName", choices = symbolOptions, server = TRUE)
  updateSelectizeInput(session, inputId = "eqtlChoiceGene", choices = eqtlOptionsGene, server = TRUE)
  updateSelectizeInput(session, inputId = "eqtlChoiceTrait", choices = eqtlOptionsTrait, server = TRUE)
  
  ## Gene symbol search code----
  ### Conditional UI for if there are multiple probes for a gene ----
  #### Is it an Affy array? ----
  isMicroarray <- reactive({
    input$pData %in% c("pheno_dusek", "pheno_rusch", "pheno_meaburn1", "pheno_meaburn2")
  })
  #### Is it an Ilmn array? ----
  isBeadChip <- reactive({
    input$pData %in% c("pheno_obermoser1", "pheno_obermoser2", "pheno_obermoser3", "pheno_obermoser4")
  })
  #### Provide probe options for the UI ----
  probesForGene <- reactive({
    if (isMicroarray()) {
      return(e2$Probe[grep(input$geneName, rownames(e2))])
    } else if (isBeadChip()) {
      k <- switch(input$pData,
                  "pheno_obermoser1" = e5,
                  "pheno_obermoser2" = e6,
                  "pheno_obermoser3" = e7,
                  "pheno_obermoser4" = e8)
      return(k$Probe[grep(input$geneName, rownames(k))])
    } else {
      NULL
    }
  })
  #### Conditionally output the UI ----
  output$probeSelectUI <- renderUI({
    if (!is.null(probesForGene())) {
      selectInput("probeInput", "Select Probe:",
                  choices = probesForGene(),
                  selected = probesForGene()[1])
    } else {
      NULL
    }
  })
  
  ### Spaghetti plot output ----
  output$genePlot <- renderPlot({
    
    # Determine the expr object based on the user's selection for pData
    expr <- switch(input$pData,
                   "pheno_gomez" = e1,
                   "pheno_meaburn1" = e2,
                   "pheno_meaburn2" = e3,
                   "pheno_gosch" = e4,
                   "pheno_obermoser1" = e5,
                   "pheno_obermoser2" = e6,
                   "pheno_obermoser3" = e7,
                   "pheno_obermoser4" = e8,
                   "pheno_dusek" = e9,
                   "pheno_rusch" = e10,
                   "pheno_larocca" = e11
    )
    pdata <- switch(input$pData,
                    "pheno_gomez" = pheno_data[[1]],
                    "pheno_meaburn1" = pheno_data[[2]],
                    "pheno_meaburn2" = pheno_data[[3]],
                    "pheno_gosch" = pheno_data[[4]],
                    "pheno_obermoser1" = pheno_data[[5]],
                    "pheno_obermoser2" = pheno_data[[6]],
                    "pheno_obermoser3" = pheno_data[[7]],
                    "pheno_obermoser4" = pheno_data[[8]],
                    "pheno_dusek" = pheno_data[[9]],
                    "pheno_rusch" = pheno_data[[10]],
                    "pheno_larocca" = pheno_data[[11]]
    )
    checkprobe <- (isMicroarray() | isBeadChip()) & !is.null(input$probeInput)
    
    # Call the gene_graph function with user inputs
    if (checkprobe) {
      plot <- gene_graph_p(probeName = as.character(input$probeInput), pData = pdata,
                           column = input$column,
                           expr = expr) # Plot function for Illumina and Affymetrix microarrays
    } else {
      plot <- gene_graph(geneName = as.character(input$geneName), pData = pdata,
                         column = input$column,
                         expr = expr) # Plot function for RNA-seq experiments
    }
    
    # Return the generated plot
    return(plot)
  })
  
  ### Gene/Transcript Info output ----
  output$geneInfo <- renderDataTable({
    if (input$pData %in% c("pheno_rusch", "pheno_dusek", "pheno_meaburn1", "pheno_meaburn2")) {
      use <- anno_data$affy[anno_data$affy$affy_hg_u133_plus_2 == input$probeInput,]
    } else if (input$pData %in% c("pheno_obermoser1","pheno_obermoser2","pheno_obermoser3","pheno_obermoser4")) {
      use <- anno_data$ilmn[anno_data$ilmn$illumina_humanwg_6_v3 == input$probeInput,]
    } else if (input$pData == "pheno_gomez") {
      use <- anno_data$gomez[anno_data$gomez$ensembl_gene_id == attr(rownames(e1),"names")[which(rownames(e1) == input$geneName)],]
    } else if (input$pData == "pheno_gosch") {
      use <- anno_data$gosch[anno_data$gosch$ensembl_gene_id == attr(rownames(e4),"names")[which(rownames(e4) == input$geneName)],]
    } else if (input$pData == "pheno_larocca") {
      use <- anno_data$larocca[anno_data$larocca$ensembl_gene_id == attr(rownames(e11),"names")[which(rownames(e11) == input$geneName)],]
    }
    if(dim(use)[1] == 0 | is.na(use$ensembl_gene_id[1])) {
      return(data.frame(`Sorry` = "we couldn't find data on this gene, transcript, or probe :("))
    }
    
    ## Format the table properly
    ## Recognize if there is a "," in the gene description
    if (grepl(",", substr(use$description, start = 1, stop = str_locate(use$description, "\\[")))) {
      use$description <- gsub("\\],", "\\]/", use$description)
      skipdesc = T
    } else {
      skipdesc = F
    }
    ## Return result if there was a comma in the gene description
    if (skipdesc) {
      desc <- unlist(strsplit(use$description, "/"))
      use <- use[,-grep("description",colnames(use))]
      result <- apply(use, 2, function(x) {trimws(unlist(strsplit(x,",")))})
      ## Replace any empty columns
      result <- lapply(result, function(x) {
        if (length(x) > 0) {
          return(x)
        } else {
          return(rep("", times = length(result[["ensembl_gene_id"]])))
        }
      })
      result$description <- desc
      return(t(as.data.frame(result)))
    }
    ## Return the result if there wasn't a comma in the gene description
    result <- apply(use, 2, function(x) {trimws(unlist(strsplit(x,",")))})
    ## Replace any empty columns
    result <- lapply(result, function(x) {
      if (length(x) > 0) {
        return(x)
      } else {
        return(rep("", times = length(result[["ensembl_gene_id"]])))
      }
    })
    return(t(as.data.frame(result)))
  })
  
  ### Statistics output ----
  output$vdata <- DT::renderDataTable(expr = {
    vdata <- switch(input$pData,
                    "pheno_gomez" = variation_tables[[1]],
                    "pheno_meaburn1" = variation_tables[[2]],
                    "pheno_meaburn2" = variation_tables[[3]],
                    "pheno_gosch" = variation_tables[[4]],
                    "pheno_obermoser1" = variation_tables[[5]],
                    "pheno_obermoser2" = variation_tables[[6]],
                    "pheno_obermoser3" = variation_tables[[7]],
                    "pheno_obermoser4" = variation_tables[[8]],
                    "pheno_dusek" = variation_tables[[9]],
                    "pheno_rusch" = variation_tables[[10]],
                    "pheno_larocca" = variation_tables[[11]]
    )
    return(t(vdata[which(vdata$Symbol == input$geneName),]))
  }, selection = "single")
  
  ### Selected statistic and score outputs ----
  
  output$stableNum <- renderText({
    stable %>%
      dplyr::filter(Symbol == input$geneName) %>%
      dplyr::select(`Total`) %>%
      paste(.data, "thresholds", sep = " ")
  })
  output$stableStudies <- renderText({
    stable %>%
      dplyr::filter(Symbol == input$geneName) %>%
      dplyr::select(Thresholds) %>%
      paste()
  })
  
  output$dynamicNum <- renderText({
    dynamic %>%
      dplyr::filter(Symbol == input$geneName) %>%
      dplyr::select(Count) %>%
      paste(.data, "studies", sep = " ")
  })
  output$dynamicStudies <- renderText({
    dynamic %>%
      dplyr::filter(Symbol == input$geneName) %>%
      dplyr::select(Studies) %>%
      paste()
  })
  
  output$houseNum <- renderText({
    housekeeping %>%
      dplyr::filter(Symbol == input$geneName) %>%
      dplyr::select(Total) %>%
      paste(.data, "thresholds", sep = " ")
  })
  output$houseStudies <- renderText({
    housekeeping %>%
      dplyr::filter(Symbol == input$geneName) %>%
      dplyr::select(Thresholds) %>%
      paste()
  })
  
  #### Plot the selected statistic ----
  
  ## Displays warning if no input is detected
  output$noStatSelected1 <- renderUI({
    if (is.null(input$vdata_rows_selected)) {
      return(span(tooltip(icon("circle-exclamation"), 
                          "For Stat Density and Stat Boxplot to function, you must select a row in the table above.",
                          placement = "right")))
    } else {
      return(NULL)
    }
  })
  ## Repeat for good formatting
  output$noStatSelected2 <- renderUI({
    if (is.null(input$vdata_rows_selected)) {
      return(span(tooltip(icon("circle-exclamation"), 
                          "For Stat Density and Stat Boxplot to function, you must select a row in the table above.",
                          placement = "right")))
    } else {
      return(NULL)
    }
  })
  
  output$stablePolyScorePlot <- renderPlot({
    stnum <- stable %>%
      dplyr::filter(Symbol == input$geneName) %>%
      dplyr::select(Total) %>%
      as.numeric()
    ggplot(stable, aes(x = Total)) +
      geom_histogram(fill = "lightblue", binwidth = 1) +
      geom_vline(xintercept = stnum, color = "royalblue") +
      annotate(geom = "text", label = paste(input$geneName), x = stnum + 2, y = 5800) +
      theme_minimal_hgrid() +
      labs(x = "Stable-polymorphic score", y = "All genes")
  })
  
  output$flexibilityScorePlot <- renderPlot({
    dynum <- dynamic %>%
      dplyr::filter(Symbol == input$geneName) %>%
      dplyr::select(Count) %>%
      as.numeric()
    ggplot(dynamic, aes(x = Count)) +
      geom_histogram(fill = "lightpink1", binwidth = 1) +
      geom_vline(xintercept = dynum, color = "tomato") +
      annotate(geom = "text", label = paste(input$geneName), x = dynum + 0.75, y = 5000) +
      theme_minimal_hgrid() +
      labs(x = "Number of studies called DE", y = "All genes")
  })
  
  output$housekeepingScorePlot <- renderPlot({
    hsnum <- housekeeping %>%
      dplyr::filter(Symbol == input$geneName) %>%
      dplyr::select(Total) %>%
      as.numeric()
    ggplot(housekeeping, aes(x = Total)) +
      geom_histogram(fill = "palegreen1", binwidth = 1) +
      geom_vline(xintercept = hsnum, color = "green4") +
      annotate(geom = "text", label = paste(input$geneName), x = hsnum + 4, y = 1500) +
      theme_minimal_hgrid() +
      labs(x = "Houskeeping score", y = "All genes")
  })
  
  output$statDensityPlot <- renderPlot({
    vdata <- switch(input$pData,
                    "pheno_gomez" = variation_tables[[1]],
                    "pheno_meaburn1" = variation_tables[[2]],
                    "pheno_meaburn2" = variation_tables[[3]],
                    "pheno_gosch" = variation_tables[[4]],
                    "pheno_obermoser1" = variation_tables[[5]],
                    "pheno_obermoser2" = variation_tables[[6]],
                    "pheno_obermoser3" = variation_tables[[7]],
                    "pheno_obermoser4" = variation_tables[[8]],
                    "pheno_dusek" = variation_tables[[9]],
                    "pheno_rusch" = variation_tables[[10]],
                    "pheno_larocca" = variation_tables[[11]]
    )
    
    checkprobe <- (isMicroarray() | isBeadChip()) & !is.null(input$probeInput)
    if(!is.numeric(vdata[which(vdata$Symbol == input$geneName), colnames(vdata)[input$vdata_rows_selected]])) {
      return(NULL)
    }
    
    if (checkprobe) {
      ggplot(vdata) +
        geom_density(aes(x = .data[[colnames(vdata)[input$vdata_rows_selected]]])) +
        geom_vline(aes(xintercept = vdata[which(vdata$`Probe ID` == input$probeInput), colnames(vdata)[input$vdata_rows_selected]]), color = "royalblue") +
        theme_cowplot() +
        xlab(colnames(vdata)[input$vdata_rows_selected])
    } else {
      ggplot(vdata) +
        geom_density(aes(x = .data[[colnames(vdata)[input$vdata_rows_selected]]])) +
        geom_vline(aes(xintercept = vdata[which(vdata$Symbol == input$geneName), colnames(vdata)[input$vdata_rows_selected]]), color = "royalblue") +
        theme_cowplot() +
        xlab(colnames(vdata)[input$vdata_rows_selected])
    }
  })
  
  output$statBoxplot <- renderPlot({
    vdata <- switch(input$pData,
                    "pheno_gomez" = variation_tables[[1]],
                    "pheno_meaburn1" = variation_tables[[2]],
                    "pheno_meaburn2" = variation_tables[[3]],
                    "pheno_gosch" = variation_tables[[4]],
                    "pheno_obermoser1" = variation_tables[[5]],
                    "pheno_obermoser2" = variation_tables[[6]],
                    "pheno_obermoser3" = variation_tables[[7]],
                    "pheno_obermoser4" = variation_tables[[8]],
                    "pheno_dusek" = variation_tables[[9]],
                    "pheno_rusch" = variation_tables[[10]],
                    "pheno_larocca" = variation_tables[[11]]
    )
    
    checkprobe <- (isMicroarray() | isBeadChip()) & !is.null(input$probeInput)
    if(!is.numeric(vdata[which(vdata$Symbol == input$geneName), colnames(vdata)[input$vdata_rows_selected]])) {
      return(NULL)
    }
    
    if (checkprobe) {
      ggplot(vdata) +
        geom_boxplot(aes(y = .data[[colnames(vdata)[input$vdata_rows_selected]]])) +
        geom_hline(aes(yintercept = vdata[which(vdata$`Probe ID` == input$probeInput), colnames(vdata)[input$vdata_rows_selected]]), color = "royalblue") +
        theme_cowplot() +
        xlab(colnames(vdata)[input$vdata_rows_selected])
    } else {
      ggplot(vdata) +
        geom_boxplot(aes(y = .data[[colnames(vdata)[input$vdata_rows_selected]]])) +
        geom_hline(aes(yintercept = vdata[which(vdata$Symbol == input$geneName), colnames(vdata)[input$vdata_rows_selected]]), color = "royalblue") +
        theme_cowplot() +
        xlab(colnames(vdata)[input$vdata_rows_selected])
    }
  })
  
  ## eQTL search code ----
  ### Filter and display table ----
  output$gtexTable <- DT::renderDataTable({
    dynamicSliderOption <- if (input$dynamicTableSlider >= 0) {input$dynamicTableSlider} else {0:6}
    
    ## If gene is selected, return with no filtering by slider. Otherwise, filter by slider.
    if (input$gtexFilter == "Gene") {
      result <- gtex %>% 
        dplyr::filter(`Symbol of blood RNA` == input$eqtlChoiceGene)
      return(datatable(result, colnames = gtexColumnNames, escape = FALSE,
                       caption = tags$caption(style = 'caption-side: top; text-align: left;',
                                              tags$em('Hover over column names for more information.'))
      )
      )
    } else {
      if (input$stableInequality == "equal") {
        stableSliderOption <- if (input$stableTableSlider >= 0) {input$stableTableSlider} else {0:22}
      }
      if (input$stableInequality == "less") {
        stableSliderOption <- if (input$stableTableSlider >= 0) {0:input$stableTableSlider} else {0:22}
      }
      if (input$stableInequality == "more") {
        stableSliderOption <- if (input$stableTableSlider >= 0) {input$stableTableSlider:22} else {0:22}
      }
      if (input$dynamicInequality == "equal") {
        dynamicSliderOption <- if (input$dynamicTableSlider >= 0) {input$dynamicTableSlider} else {0:6}
      }
      if (input$dynamicInequality == "less") {
        stableSliderOption <- if (input$dynamicTableSlider >= 0) {0:input$dynamicTableSlider} else {0:6}
      }
      if (input$dynamicInequality == "more") {
        stableSliderOption <- if (input$dynamicTableSlider >= 0) {input$dynamicTableSlider:6} else {0:6}
      }
      if (input$houseInequality == "equal") {
        houseSliderOption <- if (input$houseSlider >= 0) {input$houseSlider} else {0:39}
      }
      if (input$houseInequality == "less") {
        houseSliderOption <- if (input$houseSlider >= 0) {0:input$houseSlider} else {0:39}
      }
      if (input$houseInequality == "more") {
        houseSliderOption <- if (input$houseSlider >= 0) {input$houseSlider:39} else {0:39}
      }
      
      ## If disease filter is selected, present all options that meet slider criteria and the trait of interest.
      if (input$gtexFilter == "Disease/Trait") {
        result <- gtex %>% 
          dplyr::filter(`GWAS Trait` == input$eqtlChoiceTrait, `Stable-polymorphic score` %in% stableSliderOption, `Flexibility score` %in% dynamicSliderOption, `Housekeeping score` %in% houseSliderOption)
      }
      
      ## If no filter is selected, present all options that meet slider criteria.
      if (input$gtexFilter == "No filter") {
        result <- gtex %>% 
          dplyr::filter(`Stable-polymorphic score` %in% stableSliderOption, `Flexibility score` %in% dynamicSliderOption, `Housekeeping score` %in% houseSliderOption)
      }
      ## Return gtex table with updated column names
      return(datatable(result, colnames = gtexColumnNames, escape = FALSE,
                       caption = tags$caption(style = 'caption-side: top; text-align: left;',
                                              tags$em('Hover over column names for more information.'))
      )
      )
    }
  })
  
  ### Count the number of variants for gene / trait ----
  output$gtexGraph <- renderPlot({
    if (input$gtexFilter %in% c("Gene", "No filter")) {
      index <- which(gtexplot$Var1 == input$eqtlChoiceGene)
      index <- if (index-input$gtexGraphSlider <= 0) {seq(from = 1, to = index + (input$gtexGraphSlider - index), by = 1)} else {seq(from = index - round(input$gtexGraphSlider/2), to = index + round(input$gtexGraphSlider/2), by = 1)}
      ggplot(data = gtexplot[index,], aes(x = Var1, y = Freq)) +
        geom_col(aes(fill = ifelse(Var1 == input$eqtlChoiceGene, "Target", "other"))) +
        scale_color_brewer(palette = "Paired", aesthetics = "fill") +
        labs(y = "eQTL variants", fill = "", x = "Blood RNA") +
        theme_minimal_hgrid() +
        theme(axis.text.x = element_text(angle = 75, vjust = 0.5, color = "grey30", size = input$gtexGraphSlider2))
    } else {
      index <- which(gtexplot2$Var1 == input$eqtlChoiceTrait)
      index <- if (index-input$gtexGraphSlider <= 0) {seq(from = 1, to = index + (input$gtexGraphSlider - index), by = 1)} else {seq(from = index - round(input$gtexGraphSlider/2), to = index + round(input$gtexGraphSlider/2), by = 1)}
      ggplot(data = gtexplot2[index,], aes(x = Var1, y = Freq)) +
        geom_col(aes(fill = ifelse(Var1 == input$eqtlChoiceTrait, "Target", "other"))) +
        scale_color_brewer(palette = "Paired", aesthetics = "fill") +
        labs(y = "eQTL variants", fill = "", x = "Trait") +
        theme_minimal_hgrid() +
        theme(axis.text.x = element_text(angle = 75, vjust = 0.5, color = "grey30", size = input$gtexGraphSlider2))
    }
  })
  
  
  ## Downloads ----
  
  datasetInput <- reactive({
    switch(input$downloadSelectInput,
           "gomez" = variation_tables[[1]],
           "m1" = variation_tables[[2]],
           "m2" = variation_tables[[3]],
           "gosch" = variation_tables[[4]],
           "o1" = variation_tables[[5]],
           "o2" = variation_tables[[6]],
           "o3" = variation_tables[[7]],
           "o4" = variation_tables[[8]],
           "dusek" = variation_tables[[9]],
           "rusch" = variation_tables[[10]],
           "larocca" = variation_tables[[11]],
           "stable" = stable,
           "flex" = dynamic,
           "house" = housekeeping)
  })
  
  fn <- reactive({
    switch(input$downloadSelectInput,
           "gomez" = "gomez_variation",
           "m1" = "meaburn_day1_variation",
           "m2" = "meaburn_day2_variation",
           "gosch" = "gosch_variation",
           "o1" = "obermoser_cohort1_vein_variation",
           "o2" = "obermoser_cohort2_finger_2wk_variation",
           "o3" = "obermoser_cohort2_vein_variation",
           "o4" = "obermoser_cohort2_finger_9d_variation",
           "dusek" = "dusek_variation",
           "rusch" = "rusch_variation",
           "larocca" = "larocca_variation",
           "stable" = "stablepolymorphic_genes",
           "flex" = "flexible_genes",
           "house" = "housekeeping_genes")
  })
  
  ### Preview ----
  output$downloadSelectPreview <- renderDataTable({
    preview <- datasetInput()
    return(datatable(preview, caption = "Data preview"))
  })
  
  ### Download CSV ----
  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = function(
    ) {
      paste(fn(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(datasetInput(), file, row.names = FALSE)
    }
  )
  
  
  ## Code for homepage----
  
  # Render info tables
  output$gomezTableHead <- DT::renderDataTable(head(pheno_data[[1]]))
  output$goschTableHead <- DT::renderDataTable(head(pheno_data[[4]]))
  output$obermoserTableHead1 <- DT::renderDataTable(datatable(head(pheno_data[[5]]), caption = "Cohort 1, Vein"))
  output$obermoserTableHead2 <- DT::renderDataTable(datatable(head(pheno_data[[7]]), caption = "Cohort 2, Vein"))
  output$DusekTableHead <- DT::renderDataTable(head(pheno_data[[9]]))
  output$ruschTableHead <- DT::renderDataTable(head(pheno_data[[10]]))
  output$laroccaTableHead <- DT::renderDataTable(head(pheno_data[[11]]))
  output$meaburnTableHead <- DT::renderDataTable(head(pheno_data[[2]]))
  
  # Feedback code
  observeEvent(input$submitFeedback, {
    if(input$feedbackName == "" | !grepl("^[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\\.[A-Za-z]{2,}$", input$feedbackName)) {return(showNotification("You must provide a valid email address and a title for your request.", type = "warning"))}
    if(input$feedbackTitle == "") {return(showNotification("You must provide a valid email address and a title for your request.", type = "warning"))}
    sender <- "willbaltazar01@gmail.com"
    recipients <- "willbaltazar01@gmail.com"
    subject <- paste("RShiny feedback from", input$feedbackName, ":", input$feedbackTitle)
    body <- paste("Sender:", input$feedbackName, "\nSubject:", input$feedbackTitle, "\n\nFeedback:\n", input$feedback)
    
    send.mail(from = sender,
              to = recipients,
              subject = subject,
              body = body,
              smtp = list(host.name = "smtp.gmail.com", port = 465, user.name = "willbaltazar01@gmail.com", passwd = "iefa zkxi pkpi fsmd
", ssl = TRUE),
              authenticate = TRUE,
              send = TRUE)
    output$feedbackResult <- renderText("Thank you for your feedback! If you'd like to submit more, please refresh the page (prevents spam)")
    updateActionButton(session = session, inputId = "submitFeedback", disabled = TRUE)
  })
}

# Run the app
shinyApp(ui, server)