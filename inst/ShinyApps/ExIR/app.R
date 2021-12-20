# Required packages
library(shiny)
library(shinythemes)
library(DT)
library(shinyWidgets)
library(shinyjs)
library(shinycssloaders)
library(colourpicker)
library(janitor)
library(ranger)
library(coop)
library(influential)
library(ggplot2)
library(igraph)

options(shiny.maxRequestSize = Inf)
options(warn=-1)

####**********************************************####

navbarPageWithText <- function(..., text) {

    if(as.integer(paste(unlist(packageVersion(pkg = "shiny")), collapse = "")) <= 160) {
        navbar <- navbarPage(...)
        textEl <- tags$p(class = "navbar-text", text)
        navbar[[3]][[1]]$children[[1]] <- htmltools::tagAppendChild(
            navbar[[3]][[1]]$children[[1]],
            textEl)
        navbar
    } else {
        navbar <- navbarPage(...)
        textEl <- tags$p(class = "navbar-text", text)
        navbar[[4]][[1]][[1]]$children[[1]] <- htmltools::tagAppendChild(
            navbar[[4]][[1]][[1]]$children[[1]],
            textEl)
        navbar
    }
}

####**********************************************####

fixUploadedFilesNames <- function(x) {
    if (is.null(x)) {
        return()
    }

    oldNames = x$datapath
    newNames = file.path(dirname(x$datapath),
                         x$name)
    file.rename(from = oldNames, to = newNames)
    x$datapath <- newNames
    x
}

##***********************##

# Define diff_data.assembly for shiny app

diff_data.assembly_for_App <-
    function (List)
    {
        datasets <- lapply(List, as.data.frame)
        feature.names <- unique(unlist(lapply(X = datasets, FUN = rownames)))
        Diff_data <- data.frame(Diff_value1 = rep(0, length(feature.names)),
                                row.names = feature.names)
        for (i in 1:length(datasets)) {
            feature.names.index <- match(rownames(datasets[[i]]),
                                         rownames(Diff_data))
            if (ncol(datasets[[i]]) == 2) {
                Diff_data[, paste("Diff_value", i, sep = "")] <- 0
                Diff_data[, paste("Diff_value", i, sep = "")][feature.names.index] <- datasets[[i]][,
                                                                                                    1]
                Diff_data[, paste("Sig_value", i, sep = "")] <- 1
                Diff_data[, paste("Sig_value", i, sep = "")][feature.names.index] <- datasets[[i]][,
                                                                                                   2]
            }
            else if (ncol(datasets[[i]]) == 1) {
                Diff_data[, paste("Diff_value", i, sep = "")] <- 0
                Diff_data[, paste("Diff_value", i, sep = "")][feature.names.index] <- datasets[[i]][,
                                                                                                    1]
            }
        }
        return(Diff_data)
    }

####***********************************************************####

# The App

ui <- navbarPageWithText(id = "inTabset",
                         text = "ExIR: An Elixir for Biologists",
                         theme = shinythemes::shinytheme(theme = "sandstone"),
                         titlePanel(title = NULL,
                                    windowTitle = "Experimental data-based integrative ranking (ExIR)"),
                         # Prerequisites
                         tabPanel(title = "Prerequisites", value = "Prerequisites", icon = icon("elementor"),
                                  shinyjs::useShinyjs(),

                                  # Define vertical tabs
                                  verticalTabsetPanel(
                                      selected = NULL,
                                      id = NULL,
                                      color = "#1b0251",
                                      contentWidth = 10,
                                      menuSide = "left",

                                      # Normalized Experimental Data tab
                                      verticalTabPanel(title = "Normalized Experimental Data",
                                                       icon = NULL,
                                                       box_height = "100px",
                                  fluidRow(

                                      # Experimental data
                                      column(12, style = "background-color:#F6FFFE;",
                                             tags$h3(tags$b(icon("vial"), "Normalized Experimental Data")),
                                             p(style="text-align: justify",
                                             "In case you don't have access to the experimental data, you may either ask the bioinformatician/data manager
                                               at your lab/institute who has analyzed the raw data to run the ExIR model or ask them to send you the normalized
                                               experimental data."),
                                             br(),
                                             tags$h5(tags$b("1. What types of experimental data could be used for running the ExIR?")),
                                             p(style="text-align: justify",
                                             "Generally speaking, any type of experimental data such as transcriptomic,
                                               proteomic, and metabolomic data could be used for running the ExIR model."),
                                             p("The only requirements are:"),
                                             tags$ul(
                                             tags$li(style="text-align: justify",
                                                 "The inclusion of at least two conditions (e.g. Treatment vs Control) or time-points in the dataset"),
                                             tags$li(style="text-align: justify",
                                             "Previous normalization and log2 transformation of the experimental data. This should be done according to the exact type of
                                                     experimental data that you are working on. For example, the normalization of single-cell and bulk
                                                     RNA-seq data is different. In case your data is not log2 transformed, there is an option in the app to
                                                     automatically log2 transform the data prior to running the model. However, please note that this would not
                                                     sobstitute the data normalization. You may ask your bioinformatician to give you access to the normalized data,
                                                     not the raw data!")
                                             ),
                                             br(),
                                             tags$h5(style="text-align: justify",
                                                 tags$b("2. The experimental data generated by which platforms could be used for running the ExIR?")),
                                                p(style="text-align: justify",
                                                  "Literally, there is no limitation to the platform used for data generation and you may use any
                                                  platform for the generation of your transcriptomic, proteomic, or metabolomic data."),
                                             br(),
                                             tags$h5(style="text-align: justify",
                                                     tags$b("3. Is there any limitation to the size of the dataset?")),
                                             p(style="text-align: justify",
                                               "As far as we have tested the ExIR model, there is no limit to the size of the dataset, neither the
                                               number of samples nor the number of features (e.g. genes)"),
                                             br(),
                                             tags$h5(style="text-align: justify",
                                                     tags$b("4. What file types are acceptible for inputting into the ExIR model?")),
                                             p(style="text-align: justify",
                                               "Two file types are acceptible for inputting into the ExIR model as follows:"),
                                             tags$ul(
                                             tags$li(style="text-align: justify",
                                                     p(tags$b("CSV:"), "A comma-separted value (CSV) file (.csv)")),
                                             tags$li(style="text-align: justify",
                                                     p(tags$b("TXT:"), "A tab-delimited file (.txt)"))
                                             ),
                                             br(),
                                             tags$h5(style="text-align: justify",
                                                     tags$b("5. What are the required features of the dataset of the experimental data?")),
                                             p(style="text-align: justify",
                                               "There are three things to consider when preparing an experimental data file as follows:"),
                                             tags$ul(
                                                 tags$li(style="text-align: justify",
                                                         p(tags$b("Samples on the columns:"), "The samples/cells should come on
                                                           the columns of the dataset")),
                                                 tags$li(style="text-align: justify",
                                                         p(tags$b("Features on the rows:"), "The features (i.e. genes, proteins,
                                                           metabolites, etc.) should come on the rows of the dataset)")),
                                                 tags$li(style="text-align: justify",
                                                         p(tags$b("Condition of samples:"), "A row should be added to the dataset
                                                           for defining the condition/state of samples/cells (e.g. Tumor or Healthy)"))
                                             ),
                                             p("You can see below a simple example of a properly structured dataset."),
                                             tableOutput("exprt_sample_table")
                                             )
                                      )
                                  ),

                                  # Differential Data tab
                                  verticalTabPanel(title = "Differential Data",
                                                   icon = NULL,
                                                   box_height = "100px",
                                      column(12, style = "background-color:#FAF6FF;",
                                             tags$h3(tags$b(icon("chart-line"), "Differential Data")),
                                             p(style="text-align: justify",
                                               "In case you don't have access to the differential data, you may either ask the bioinformatician/data manager
                                               at your lab/institute who has analyzed the raw data to run the ExIR model or ask them to send you a table of
                                               significant (filtered) differential data."),
                                             br(),
                                             tags$h5(style="text-align: justify",
                                                     tags$b("1. What file types are acceptible for inputting into the ExIR model?")),
                                             p(style="text-align: justify",
                                               "Two file types are acceptible for inputting into the ExIR model as follows:"),
                                             tags$ul(
                                                 tags$li(style="text-align: justify",
                                                         p(tags$b("CSV:"), "A comma-separted value (CSV) file (.csv)")),
                                                 tags$li(style="text-align: justify",
                                                         p(tags$b("TXT:"), "A tab-delimited file (.txt)"))
                                             ),
                                             br(),
                                             tags$h5(style="text-align: justify",
                                                     tags$b("2. What are the required features of the differential dataset")),
                                             p(style="text-align: justify",
                                               "There are three things to consider when preparing differential data files as \
                                               follows:"),
                                             tags$ul(
                                                 tags$li(style="text-align: justify",
                                                         p(tags$b("Features on the rows:"), "The features (i.e. genes, proteins,
                                                           metabolites, etc.) should come on the rows of the datasets)")),
                                                 tags$li(style="text-align: justify",
                                                         p(tags$b("Fold changes and  significance values on
                                                         columns:"), "Fold changes and their significance
                                                           values (e.g. P-value, adjusted P-value, etc.) should come on the columns of the dataset")),
                                                 tags$li(style="text-align: justify",
                                                         p(tags$b("Separate differential files:"), "In case your study includes more than two
                                                           conditions (e.g. different doses of drug treatments) or time-points (TPs), you
                                                           should upload separate differential datasets for each step. For instance,
                                                           if your study includes three TPs, you should provide a differential dataset
                                                           for TP2 vs TP1 and another dataset for TP3 vs TP2. Please also note that
                                                           the significance (adjusted P-value) column is mandatory for differential
                                                           datasets. Also, note that different differential datasets include only the
                                                           significantly differentially expressed features (genes) and, consequently,
                                                           different differential datasets could include both unique and common genes.")),
                                                 p("You can see below a simple example of differential datasets of a three TP study."),
                                                 fluidRow(
                                                     column(6,
                                                            tags$h5(tags$b("TP2 vs. TP1")),
                                                            tableOutput("differential_table1")
                                                            ),
                                                     column(6,
                                                            tags$h5(tags$b("TP3 vs. TP2")),
                                                            tableOutput("differential_table2")
                                                            )
                                                 )
                                                 )
                                             )
                                      ),

                                  # Regression Data tab
                                  verticalTabPanel(title = "Regression Data",
                                                   icon = NULL,
                                                   box_height = "100px",
                                                   column(12, style = "background-color:#FFFAF6;",
                                                   tags$h3(tags$b(icon("chart-line"), "Regression Data (Optional)")),
                                                   p(style="text-align: justify",
                                                     "In case your study includes more than two conditions or time-points (TPs), you may
                                                     optionally upload a separate dataset for the regression data.
                                                     If you don't have access to the regression data, you
                                                     may either ask the bioinformatician/data manager
                                               at your lab/institute who has analyzed the raw data to run the ExIR model or ask them to
                                               send you a table of regression data, if applicable."),
                                                   br(),
                                                   tags$h5(style="text-align: justify",
                                                           tags$b("1. What file types are acceptible for inputting into the ExIR model?")),
                                                   p(style="text-align: justify",
                                                     "Two file types are acceptible for inputting into the ExIR model as follows:"),
                                                   tags$ul(
                                                       tags$li(style="text-align: justify",
                                                               p(tags$b("CSV:"), "A comma-separted value (CSV) file (.csv)")),
                                                       tags$li(style="text-align: justify",
                                                               p(tags$b("TXT:"), "A tab-delimited file (.txt)"))
                                                   ),
                                                   br(),
                                                   tags$h5(style="text-align: justify",
                                                           tags$b("2. What are the required features of the regression dataset")),
                                                   p(style="text-align: justify",
                                                     "There are two things to consider when preparing a regression dataset as follows:"),
                                                   tags$ul(
                                                       tags$li(style="text-align: justify",
                                                               p(tags$b("Features on the rows:"), "The features (i.e. genes, proteins,
                                                           metabolites, etc.) should come on the rows of the datasets.")),
                                                       tags$li(style="text-align: justify",
                                                               p(tags$b("Regression values and significance values on
                                                         columns:"), "Regression values as well as their significance
                                                           values (e.g. P-value, adjusted P-value, etc.) should come on the columns of the dataset.
                                                                 Please also note that the significance (adjusted P-value) column is not mandatory
                                                     for regression dataset.")
                                                               ),
                                                       ),
                                                 p("You can see below a simple example of regression (trajectory) data with and without
                                                   significance values (adjusted P-value)."),
                                                 fluidRow(
                                                     column(6,
                                                                tags$h5(tags$b("Regression dataset with significance values")),
                                                            tableOutput("regression_table1")
                                                     ),
                                                     column(6,
                                                            tags$h5(tags$b("Regression dataset without significance values")),
                                                            tableOutput("regression_table2")
                                                     )
                                             )
                                             )
                                             ),

                                  # Desired List of Features tab
                                  verticalTabPanel(title = "Desired List of Features",
                                                   icon = NULL,
                                                   box_height = "100px",
                                      column(12, style = "background-color:#FFF6F7;",
                                             tags$h3(tags$b(icon("list-alt"), "Desired List of Features (Optional)")),
                                             p(style="text-align: justify",
                                               "The desired list of features (e.g. genes, proteins, metabolites, etc.) is optional and you could just
                                               ignore it so that the ExIR model will be run and trained on the entire dataset. However, you could upload
                                               a list of desired features to train the model asccording to that list. This is useful when you have a
                                               handful of candidate features (e.g. based on your previous assasys or according to the literature) and you
                                               would like to functionally classify them (into drivers, biomarkers, and mediators) and prioritize them for
                                               downstream experimental functional validations."),
                                             br(),
                                             tags$h5(style="text-align: justify",
                                                     tags$b("1. What file types are acceptible for inputting into the ExIR model?")),
                                             p(style="text-align: justify",
                                               "Two file types are acceptible for inputting into the ExIR model as follows:"),
                                             tags$ul(
                                                 tags$li(style="text-align: justify",
                                                         p(tags$b("CSV:"), "A comma-separted value (CSV) file (.csv)")),
                                                 tags$li(style="text-align: justify",
                                                         p(tags$b("TXT:"), "A tab-delimited file (.txt)"))
                                             ),
                                             br(),
                                             tags$h5(style="text-align: justify",
                                                     tags$b("2. Is there any limitation to the size (number of features)
                                                            of the desired list of features?")),
                                             p(style="text-align: justify",
                                               "As far as we have tested the ExIR model, there is no limit to the size of the desired
                                               list of features."),
                                             br(),
                                             tags$h5(style="text-align: justify",
                                                     tags$b("3. What are the required characteristics of the desired list of features?")),
                                             p(style="text-align: justify",
                                               "There are four things to consider when preparing the desired list of features as follows:"),
                                             tags$ul(
                                                 tags$li(style="text-align: justify",
                                                         p(tags$b("A single column file:"), "The desired list of features should be in a single
                                                           column .txt or .csv file)")),
                                                 tags$li(style="text-align: justify",
                                                         p(tags$b("No header:"), "The single column of desired features should not include
                                                           any header")),
                                                 tags$li(style="text-align: justify",
                                                         p(tags$b("Included in the normalized experimental data:"), "Apparently, all of the desired feature
                                                           should be included in the normalized experimental data and no external gene without any
                                                           information is acceptable")),
                                                 tags$li(style="text-align: justify",
                                                         p(tags$b("Identical feature names:"), "synonym names are not allowed, and you should use
                                                           exactly the same names of features that you have used in the normalized expermental dataset"))
                                             ),
                                             p("You can see below a simple example of a properly structured desired list of features."),
                                             fluidRow(
                                                 column(6,
                                                        tags$h5(tags$b("Sample desired list of features")),
                                                        tableOutput("desiredFeaturesSampleTable")
                                                 ),
                                                 column(6,
                                                        tags$h5(tags$b("A real example of desired list of features including five genes")),
                                                        tableOutput("desiredFeaturesRealTable")
                                                 )
                                             )
                                      )
                                  ),

                                  # Synonyms table
                                  verticalTabPanel(title = "Synonyms Table",
                                                   icon = NULL,
                                                   box_height = "100px",                                      # Desired list
                                                   column(12, style = "background-color:#FCFFF6;",
                                                          tags$h3(tags$b(icon("table"), "Synonyms Table (Optionally required
                                                                         for the visualization of ExIR results)")),
                                                          p(style="text-align: justify",
                                                            "The synonyms table is not required for running the ExIR model and will be
                                                            used for the visualization of the ExIR output. The purpose of the sysnonyms
                                                            table is simple! Consider a situation when you/your bioinformatician have
                                                            done data processing, normalization, and differential expression analysis on
                                                            the dataset in which gene names are in Ensembl gene name format (Ensembl ID) and,
                                                            consequently, the feature names in all of your files and datasets are in the
                                                            Ensembl ID format. At the end of the day, any figure you generate based on
                                                            the results of the ExIR model will show the Ensembl ID of genes. However, you
                                                            might want to show the symbol of genes (gene names) in the figure instead of
                                                            Ensembl IDs. So, instead of revising all of the files and datasets before inputting
                                                            into the ExIR model, you can input the files as they are and, instead, generate a
                                                            synonyms table for visualization of the results.
                                                            "),
                                                          br(),
                                                          tags$h5(style="text-align: justify",
                                                                  tags$b("1. Where can I find the synonyms of my features (genes/proteins)?")),
                                                          p(style="text-align: justify",
                                                            "There are several web servers, online databases as well as R packages for this purpose.
                                                            Some of these online tools and R packages are exampled below:"),
                                                          tags$ul(
                                                              tags$li(style="text-align: justify",
                                                                      p(tags$b("g:Profiler:"), tags$a("g:Profiler",
                                                                      href = "https://biit.cs.ut.ee/gprofiler/convert", style = "color:blue"),
                                                                      " is a public web server for characterising
                                                                        and manipulating gene lists.")),
                                                              tags$li(style="text-align: justify",
                                                                      p(tags$b("DAVID Gene ID Conversion Tool:"), tags$a("DAVID Gene ID Conversion Tool",
                                                                                                                         href = "https://david.ncifcrf.gov/home.jsp",
                                                                                                                         style = "color:blue"),
                                                                                                                         " is a gene ID conversion tool on the Database for
                                                                                                                         Annotation, Visualization and Integrated Discovery (DAVID) web server.")),
                                                              tags$li(style="text-align: justify",
                                                                      p(tags$b("bioDBnet: db2db:"), tags$a("bioDBnet: db2db",
                                                                                                                         href = "https://biodbnet.abcc.ncifcrf.gov/db/db2db.php",
                                                                                                                         style = "color:blue"), " allows for conversions of identifiers
                                                                        from one database to other database identifiers or annotations.")),
                                                              tags$li(style="text-align: justify",
                                                                      p(tags$b("R package mygene:"), tags$a("R package mygene",
                                                                                                           href = "https://bioconductor.org/packages/release/bioc/html/mygene.html",
                                                                                                           style = "color:blue"), " is an easy-to-use R wrapper to access
                                                                        MyGene.Info_ services. MyGene.Info_ provides simple-to-use REST web services to query/retrieve gene annotation data.")),
                                                              tags$li(style="text-align: justify",
                                                                      p(tags$b("R package gprofiler2:"), tags$a("R package gprofiler2",
                                                                                                           href = "https://cran.r-project.org/web/packages/gprofiler2/index.html",
                                                                                                           style = "color:blue"), " is a toolset for functional enrichment analysis
                                                                        and visualization, gene/protein/SNP identifier conversion and mapping orthologous genes across species."))

                                                          ),
                                                          br(),
                                                          tags$h5(style="text-align: justify",
                                                                  tags$b("2. What file types are acceptible for inputting into the visualization tool?")),
                                                          p(style="text-align: justify",
                                                            "Two file types are acceptible for inputting into the visualization tool as follows:"),
                                                          tags$ul(
                                                              tags$li(style="text-align: justify",
                                                                      p(tags$b("CSV:"), "A comma-separted value (CSV) file (.csv)")),
                                                              tags$li(style="text-align: justify",
                                                                      p(tags$b("TXT:"), "A tab-delimited file (.txt)"))
                                                          ),
                                                          br(),
                                                          tags$h5(style="text-align: justify",
                                                                  tags$b("3. What should be the number of rows in the synonyms table?")),
                                                          p(style="text-align: justify",
                                                            "The number of rows should be equal to the number of features on the columns of the normalized
                                                            experimental data."),
                                                          br(),
                                                          tags$h5(style="text-align: justify",
                                                                  tags$b("4. What are the required features of the synonyms table?")),
                                                          p(style="text-align: justify",
                                                            "There are four things to consider when preparing the synonyms table as follows:"),
                                                          tags$ul(
                                                              tags$li(style="text-align: justify",
                                                                      p(tags$b("Features on the rows:"), "The features (i.e. genes, proteins,
                                                           metabolites, etc.) should come on the rows of the table")),
                                                              tags$li(style="text-align: justify",
                                                                      p(tags$b("A two column file:"), "The synonyms table should be in a two
                                                           column .txt or .csv file)")),
                                                              tags$li(style="text-align: justify",
                                                                      p(tags$b("First column as the original names:"), "The first column of the table should
                                                                        include the original names of features (the same ones used in the normalized
                                                                        experimental data).")),
                                                              tags$li(style="text-align: justify",
                                                                      p(tags$b("Second column as the synonym names:"), "The second column of the table should
                                                                        include the the synonyms of the original names."))
                                                          ),
                                                          p("You can see below a simple example of a properly structured synonyms table."),
                                                          fluidRow(
                                                              column(6,
                                                                     tags$h5(tags$b("Sample synonyms table showing only the first five genes")),
                                                                     tableOutput("synonymsSampleTable")
                                                              ),
                                                              column(6,
                                                                     tags$h5(tags$b("A real example of a synonyms table showing only the first five genes")),
                                                                     tableOutput("synonymsRealTable")
                                                              )
                                                          )
                                                   )
                                  )
                                      )
                                  ),

                         ####*****************************************####

                         # Running ExIR
                         tabPanel(title = "Run ExIR", value = "runExIR", icon = icon("robot"),
                                  shinyjs::useShinyjs(),
                                  fluidRow(
                                      column(4,
                                             tags$h5("Running the Experimental-data-based Integrative Ranking (ExIR)"),
                                             sidebarPanel(width = 12,
                                                          style = "overflow-y:scroll; max-height: 90vh; position:relative;",
                                                          ## Input for running ExIR

                                                          ### Upload the files
                                                          fileInput("restore_exir_dataset", label = "Restore ExIR results from file (.rds) ",
                                                                    placeholder = ".rds file",
                                                                    accept = ".rds",
                                                                    buttonLabel = "Browse"),
                                                          fileInput("normExptlData", label = p(icon("upload"), "Upload the Normalized Experimental Data", style = "padding:0px; margin:0;"),
                                                                    accept = c(".csv", ".txt"),
                                                                    buttonLabel = "Browse"),

                                                          numericInput(inputId = "n_DifferentialDatasets",
                                                                       label = "Number of Differential Datasets:",
                                                                       value = 1, min = 1, step = 1, max = 10),

                                                          fileInput("DifferentialDataset1", label = p(icon("upload"), "Upload the Differential Dataset (#1)", style = "padding:0px; margin:0;"),
                                                                    accept = c(".csv", ".txt"),
                                                                    multiple = FALSE,
                                                                    buttonLabel = "Browse"),

                                                          conditionalPanel(condition = "input.n_DifferentialDatasets == 2 |
                                                                           input.n_DifferentialDatasets == 3 |
                                                                           input.n_DifferentialDatasets == 4 |
                                                                           input.n_DifferentialDatasets == 5 |
                                                                           input.n_DifferentialDatasets == 6 |
                                                                           input.n_DifferentialDatasets == 7 |
                                                                           input.n_DifferentialDatasets == 8 |
                                                                           input.n_DifferentialDatasets == 9 |
                                                                           input.n_DifferentialDatasets == 10",

                                                                           fileInput("DifferentialDataset2", label = p(icon("upload"), "Upload the Differential Dataset (#2)", style = "padding:0px; margin:0;"),
                                                                                     accept = c(".csv", ".txt"),
                                                                                     multiple = FALSE,
                                                                                     buttonLabel = "Browse")
                                                          ),

                                                          conditionalPanel(condition = "input.n_DifferentialDatasets == 3 |
                                                                           input.n_DifferentialDatasets == 4 |
                                                                           input.n_DifferentialDatasets == 5 |
                                                                           input.n_DifferentialDatasets == 6 |
                                                                           input.n_DifferentialDatasets == 7 |
                                                                           input.n_DifferentialDatasets == 8 |
                                                                           input.n_DifferentialDatasets == 9 |
                                                                           input.n_DifferentialDatasets == 10",

                                                                           fileInput("DifferentialDataset3", label = p(icon("upload"), "Upload the Differential Dataset (#3)", style = "padding:0px; margin:0;"),
                                                                                     accept = c(".csv", ".txt"),
                                                                                     multiple = FALSE,
                                                                                     buttonLabel = "Browse")
                                                          ),

                                                          conditionalPanel(condition = "input.n_DifferentialDatasets == 4 |
                                                                           input.n_DifferentialDatasets == 5 |
                                                                           input.n_DifferentialDatasets == 6 |
                                                                           input.n_DifferentialDatasets == 7 |
                                                                           input.n_DifferentialDatasets == 8 |
                                                                           input.n_DifferentialDatasets == 9 |
                                                                           input.n_DifferentialDatasets == 10",

                                                                           fileInput("DifferentialDataset4", label = p(icon("upload"), "Upload the Differential Dataset (#4)", style = "padding:0px; margin:0;"),
                                                                                     accept = c(".csv", ".txt"),
                                                                                     multiple = FALSE,
                                                                                     buttonLabel = "Browse")
                                                          ),

                                                          conditionalPanel(condition = "input.n_DifferentialDatasets == 5 |
                                                                           input.n_DifferentialDatasets == 6 |
                                                                           input.n_DifferentialDatasets == 7 |
                                                                           input.n_DifferentialDatasets == 8 |
                                                                           input.n_DifferentialDatasets == 9 |
                                                                           input.n_DifferentialDatasets == 10",

                                                                           fileInput("DifferentialDataset5", label = p(icon("upload"), "Upload the Differential Dataset (#5)", style = "padding:0px; margin:0;"),
                                                                                     accept = c(".csv", ".txt"),
                                                                                     multiple = FALSE,
                                                                                     buttonLabel = "Browse")
                                                          ),

                                                          conditionalPanel(condition = "input.n_DifferentialDatasets == 6 |
                                                                           input.n_DifferentialDatasets == 7 |
                                                                           input.n_DifferentialDatasets == 8 |
                                                                           input.n_DifferentialDatasets == 9 |
                                                                           input.n_DifferentialDatasets == 10",

                                                                           fileInput("DifferentialDataset6", label = p(icon("upload"), "Upload the Differential Dataset (#6)", style = "padding:0px; margin:0;"),
                                                                                     accept = c(".csv", ".txt"),
                                                                                     multiple = FALSE,
                                                                                     buttonLabel = "Browse")
                                                          ),

                                                          conditionalPanel(condition = "input.n_DifferentialDatasets == 7 |
                                                                           input.n_DifferentialDatasets == 8 |
                                                                           input.n_DifferentialDatasets == 9 |
                                                                           input.n_DifferentialDatasets == 10",

                                                                           fileInput("DifferentialDataset7", label = p(icon("upload"), "Upload the Differential Dataset (#7)", style = "padding:0px; margin:0;"),
                                                                                     accept = c(".csv", ".txt"),
                                                                                     multiple = FALSE,
                                                                                     buttonLabel = "Browse")
                                                          ),

                                                          conditionalPanel(condition = "input.n_DifferentialDatasets == 8 |
                                                                           input.n_DifferentialDatasets == 9 |
                                                                           input.n_DifferentialDatasets == 10",

                                                                           fileInput("DifferentialDataset8", label = p(icon("upload"), "Upload the Differential Dataset (#8)", style = "padding:0px; margin:0;"),
                                                                                     accept = c(".csv", ".txt"),
                                                                                     multiple = FALSE,
                                                                                     buttonLabel = "Browse")
                                                          ),

                                                          conditionalPanel(condition = "input.n_DifferentialDatasets == 9 |
                                                                           input.n_DifferentialDatasets == 10",

                                                                           fileInput("DifferentialDataset9", label = p(icon("upload"), "Upload the Differential Dataset (#9)", style = "padding:0px; margin:0;"),
                                                                                     accept = c(".csv", ".txt"),
                                                                                     multiple = FALSE,
                                                                                     buttonLabel = "Browse")
                                                          ),

                                                          conditionalPanel(condition = "input.n_DifferentialDatasets == 10",
                                                                           fileInput("DifferentialDataset10", label = p(icon("upload"), "Upload the Differential Dataset (#10)", style = "padding:0px; margin:0;"),
                                                                                     accept = c(".csv", ".txt"),
                                                                                     multiple = FALSE,
                                                                                     buttonLabel = "Browse")
                                                          ),

                                                          fileInput("RegressionDataset", label = p(icon("upload"), "Upload the Regression Dataset (optional)", style = "padding:0px; margin:0;"),
                                                                    accept = c(".csv", ".txt"),
                                                                    buttonLabel = "Browse"),

                                                          fileInput("desiredFeaturesList", label = p(icon("upload"), "Upload the Desired Features file (optional)", style = "padding:0px; margin:0;"),
                                                                    accept = c(".csv", ".txt"),
                                                                    buttonLabel = "Browse"),

                                                          ### Specify the condition row
                                                          textInput(inputId = "conditionRowname",
                                                                    label = "Condition row name:",
                                                                    value = NULL, placeholder = "Condition row name"),

                                                          ### Specify the r
                                                          sliderInput(inputId = "correlationCoeff",
                                                                      label = "Correlation coefficient threshold (for association network reconstruction):",
                                                                      min = 0, max = 1, value = 0.5, step = 0.1,
                                                                      ticks = TRUE, animate = FALSE),

                                                          # max.connections
                                                          numericInput(inputId = "max.connections",
                                                                       label = "Maximum number of connections (for association network reconstruction):",
                                                                       value = 20000, min = 10000, max = Inf, step = 5000),

                                                          ### alpha
                                                          sliderInput(inputId = "alpha",
                                                                      label = "Statistical significance threshold (alpha):",
                                                                      min = 0.01, max = 0.1, value = 0.05, step = 0.01,
                                                                      ticks = TRUE, animate = FALSE),

                                                          # Number of trees
                                                          numericInput(inputId = "num_trees",
                                                                       label = "Number of trees (for random forests classification):",
                                                                       value = 10000, min = 1000, max = 100000, step = 5000),

                                                          # mtry_option
                                                          tags$b(p("Number of features in each node (for random forests classification):")),
                                                          switchInput(
                                                              value = TRUE,
                                                              inputId = "mtry_option",
                                                              label = "Automatic",
                                                              labelWidth = "80px"
                                                          ),

                                                          # mtry
                                                          numericInput(inputId = "mtry",
                                                                       label = "Determine the numer of features:",
                                                                       value = NULL),

                                                          # num_permutations
                                                          numericInput(inputId = "num_permutations",
                                                                       label = "Number of permutations (for calculation of P-values):",
                                                                       value = 100, min = 10, max = Inf, step = 50),

                                                          # normalize
                                                          materialSwitch(
                                                              inputId = "Normalize",
                                                              label = "Log2 transform the Experimental Data?",
                                                              value = FALSE,
                                                              status = "primary"
                                                          ),

                                                          # seed
                                                          numericInput(inputId = "seed",
                                                                       label = "Seed number (for random processes):",
                                                                       value = 1234, min = -Inf, max = Inf, step = 10)
                                             )),

                                      ## Output ExIR model
                                      column(8,
                                             tags$h5("After retrieving the results proceed to the",
                                                     tags$em(tags$b("ExIR-based visualization")), "tab",
                                                     style = "text-align: center;"),
                                             fluidRow(
                                                 column(6,
                                                        uiOutput("runExIRbutton")
                                                        ),
                                                 column(6,
                                                        uiOutput("jumpToVisTab")
                                                 )

                                             ),
                                             br(),
                                             tags$h4(tags$b("ExIR Output Tables"), style = "text-align: center;"),

                                             # Define vertical tabs
                                             verticalTabsetPanel(
                                                 selected = "Drivers",
                                                 id = "exirOutputTabset",
                                                 color = "#2C0383",
                                                 contentWidth = 10,
                                                 menuSide = "right",

                                                 # Drivers
                                                 verticalTabPanel(title = tags$b("Drivers"),
                                                                  icon = icon("car"),
                                                                  box_height = "100px",
                                                                  fluidRow(
                                                                      htmlOutput("nullDriversTable"),
                                                                      column(11,

                                                                             DT::dataTableOutput("driversTable"),
                                                                             ),
                                                                      column(1)
                                                 )
                                                 ),

                                                 # Biomarkers
                                                 verticalTabPanel(title = tags$b("Biomarkers"),
                                                                  icon = icon("highlighter"),
                                                                  box_height = "100px",
                                                                  fluidRow(
                                                                      htmlOutput("nullBiomarkersTable"),
                                                                      column(11,

                                                                             DT::dataTableOutput("biomarkersTable"),
                                                                      ),
                                                                      column(1)
                                                                  )
                                                 ),

                                                 # Non-DE-Mediators
                                                 verticalTabPanel(title = tags$b("Non-DE-Mediators"),
                                                                  icon = icon("link"),
                                                                  box_height = "120px",
                                                                  fluidRow(
                                                                      htmlOutput("nullNonDEMediatorsTable"),
                                                                      column(11,

                                                                             DT::dataTableOutput("NonDEMediatorsTable"),
                                                                      ),
                                                                      column(1)
                                                                  )
                                                 ),

                                                 # DE-Mediators
                                                 verticalTabPanel(title = tags$b("DE-Mediators"),
                                                                  icon = icon("link"),
                                                                  box_height = "110px",
                                                                  fluidRow(
                                                                      htmlOutput("nullDEMediatorsTable"),
                                                                      column(11,

                                                                             DT::dataTableOutput("DEMediatorsTable"),
                                                                      ),
                                                                      column(1)
                                                                  )
                                                 )
                                             ),
                                             fluidRow(tags$h5(br(), "")),
                                             fluidRow(tags$h5(br(), "")),
                                             panel(id = "save_for_restore",
                                                   tags$h5(icon("bell"), "You may save the results dataset (as an RDS file) using the
                                                                          'SAVE DATASET' button for later use or sharing with a colleague.
                                                                          This would save you from redoing the analysis and
                                                                          would speed up your research!",
                                                           style = "background-color: white !important; text-align: center;"),
                                                   heading = NULL,
                                                   extra = actionButton(inputId = "closeRestoreNotif",
                                                                        label = " Close!",
                                                                        icon = icon("times-circle"),
                                                                        class = "btn-sm btn-primary",
                                                                        width = "100%"),
                                                   status = "default"
                                             ),
                                             downloadButton("save_exir_results", class = "btn-sm btn-block", "Save dataset", icon = icon("download"),
                                                            width = "100%"),
                                             fluidRow(tags$h5(br(), ""))
                                      )
                                      )
                                  ),

                         ####*****************************************####

                         # ExIR visualization

                         tabPanel(title = "ExIR-based visualization", value = "ExIRVis",
                                  icon = icon("project-diagram"),
                                  shinyjs::useShinyjs(),
                                  fluidRow(
                                      # Input for ExIR visualization
                                      column(4,
                                             tags$h5("Visualization options"),
                                             sidebarPanel(width = 12,
                                                          style = "overflow-y:scroll; max-height: 800px; position:relative;",
                                                          panel(
                                                                     fileInput("synonymsTable", label = p(icon("upload"), "Upload the Synonyms Table (optional)",
                                                                                                          style = "padding:0px; margin:0;"),
                                                                               accept = c(".csv", ".txt"),
                                                                               buttonLabel = "Browse"),
                                                                     actionButton(inputId = "resetNames",
                                                                                  label = " Reset names",
                                                                                  icon = icon("undo-alt"),
                                                                                  class = "btn-sm btn-primary",
                                                                                  width = "100%"),
                                                                     status = "default"
                                                                     ),
                                                          ### Specify the n
                                                          sliderInput(inputId = "numberOfFeatures",
                                                                      label = "Specify the number of top candidates for visualization of each class of features:",
                                                                      min = 1, max = 100, value = 10, step = 1,
                                                                      ticks = FALSE, animate = FALSE),
                                                          ### Specify driver type
                                                          switchInput(
                                                              value = TRUE,
                                                              inputId = "showDrivers",
                                                              label = "Show Drivers",
                                                              labelWidth = "150px"
                                                          ),
                                                          selectInput(inputId = "driverType",
                                                                      choices = list("Accelerator"="accelerator", "Decelerator" = "decelerator",
                                                                                     "Combined"="combined"),
                                                                      label = "Specify the type of drivers to be visualized:", selected = "combined"),
                                                          ### Specify biomarker type
                                                          switchInput(
                                                              value = TRUE,
                                                              inputId = "showBiomarkers",
                                                              label = "Show Biomarkers",
                                                              labelWidth = "150px"
                                                          ),
                                                          selectInput(inputId = "biomarkerType",
                                                                      choices = list("Up-regulated"="up-regulated", "Down-regulated" = "down-regulated",
                                                                                     "Combined"="combined"),
                                                                      label = "Specify the type of biomarkers to be visualized:", selected = "combined"),
                                                          switchInput(
                                                              value = TRUE,
                                                              inputId = "showDEMediators",
                                                              label = "Show DE-Mediators",
                                                              labelWidth = "150px"
                                                          ),
                                                          switchInput(
                                                              value = TRUE,
                                                              inputId = "showNonDEMediators",
                                                              label = "Show NonDE-Mediators",
                                                              labelWidth = "150px"
                                                          ),
                                                          selectInput(inputId = "selectionBasis",
                                                                      choices = list("Rank", "Adjusted p-value"),
                                                                      label = "Basis for the selection of top candidates for visualization:",
                                                                      selected = "Rank"),
                                                          prettyRadioButtons(inputId = "labelPosition", label = "Label position",
                                                                             choices = list("Left" = "left", "Right" = "right", "Bottom" = "bottom", "Top" = "top"),
                                                                             icon = icon("check"),
                                                                             bigger = TRUE,
                                                                             status = "info",
                                                                             animation = "jelly",
                                                                             selected = "top"),
                                                          # nrow
                                                          numericInput(inputId = "exirVis_nrow",
                                                                       label = "Determine the numer of rows of the plot:",
                                                                       value = 1),
                                                          sliderInput(inputId = "dotSizeMin",
                                                                      label = "Specify the size of dots with the lowest statistical significance:",
                                                                      min = 1, max = 10, value = 2, step = 1,
                                                                      ticks = FALSE, animate = FALSE),
                                                          sliderInput(inputId = "dotSizeMax",
                                                                      label = "Specify the size of dots with the highest statistical statistical significance:",
                                                                      min = 1, max = 10, value = 5, step = 1,
                                                                      ticks = FALSE, animate = FALSE),
                                                          selectInput(inputId = "typeColor",
                                                                      choices = list("Magma"="magma", "Inferno" = "inferno",
                                                                                     "Plasma"="plasma", "Viridis"="viridis",
                                                                                     "Cividis"="cividis"),
                                                                      label = "Specify the color palette to be used for the
                                                                      visualization of different types of candidates",
                                                                      selected = "viridis"),
                                                          numericInput(inputId = "strokeSize",
                                                                       label = "The size of stroke (border) around the dots:",
                                                                       value = 1.5, min = 0, max = 5, step = 0.5),
                                                          sliderInput(inputId = "strokeAlpha", label = "Stroke opacity (alpha):",
                                                                      min = 0, max = 1, value = 1, step = 0.1,
                                                                      ticks = FALSE, animate = FALSE),
                                                          colourInput(inputId = "dotColorLow",
                                                                      value =  "#2D29E6",
                                                                      label = "The color to be used for the visualization of dots
                                                                      (features) with the lowest Z-score values:",
                                                                      returnName = FALSE,
                                                                      closeOnClick = TRUE, allowTransparent = FALSE),
                                                          colourInput(inputId = "dotColorHigh",
                                                                      value =  "#B32C29",
                                                                      label = "The color to be used for the visualization of dots
                                                                      (features) with the highest Z-score values:",
                                                                      returnName = FALSE,
                                                                      closeOnClick = TRUE, allowTransparent = FALSE),
                                                          prettyRadioButtons(inputId = "legendPosition", label = "Legend position",
                                                                             choices = list("Left" = "left", "Right" = "right", "Bottom" = "bottom",
                                                                                            "Top" = "top", "Remove legend" = "none"),
                                                                             icon = icon("check"),
                                                                             bigger = TRUE,
                                                                             status = "info",
                                                                             animation = "jelly",
                                                                             selected = "bottom"),
                                                          prettyRadioButtons(inputId = "legendDirection", label = "Legend direction",
                                                                             choices = list("Horizontal" = "horizontal", "Vertical" = "vertical"),
                                                                             selected = "vertical",
                                                                             icon = icon("check"),
                                                                             bigger = TRUE,
                                                                             status = "info",
                                                                             animation = "jelly"),
                                                          prettyRadioButtons(inputId = "legendLayout", label = "Legend layout",
                                                                             choices = list("Horizontal" = "horizontal", "Vertical" = "vertical"),
                                                                             selected = "horizontal",
                                                                             icon = icon("check"),
                                                                             bigger = TRUE,
                                                                             status = "info",
                                                                             animation = "jelly"),
                                                          prettyCheckbox(inputId = "boxedLegend", label = "Boxed legend", value = TRUE,
                                                                         icon = icon("check"),
                                                                         status = "info"),
                                                          prettyCheckbox(inputId = "showPlotTitle", label = "Show plot title", value = TRUE,
                                                                         icon = icon("check"),
                                                                         status = "info"),
                                                          switchInput(
                                                              value = TRUE,
                                                              inputId = "autoPlotTitle",
                                                              label = "Auto-generate the plot title",
                                                              labelWidth = "200px"
                                                          ),
                                                          textInput(inputId = "plotTitle", label = "Plot title:",
                                                                    value = "ExIR-based prioritized feachers", placeholder = "Write your desired plot title"),
                                                          prettyRadioButtons(inputId = "titlePosition", label = "Title position",
                                                                             choices = list("Left" = "left", "Center" = "center", "Right" = "right"),
                                                                             selected = "left",
                                                                             icon = icon("check"),
                                                                             bigger = TRUE,
                                                                             status = "info",
                                                                             animation = "jelly"),
                                                          prettyCheckbox(inputId = "showPlotSubTitle", label = "Show plot subtitle", value = TRUE,
                                                                         icon = icon("check"),
                                                                         status = "info"),
                                                          switchInput(
                                                              value = TRUE,
                                                              inputId = "autoPlotSubTitle",
                                                              label = "Auto-generate the plot subtitle",
                                                              labelWidth = "200px"
                                                          ),
                                                          textInput(inputId = "plotSubTitle", label = "Plot subtitle:",
                                                                    value = "Top candidates", placeholder = "Write your desired plot subtitle"),
                                                          prettyRadioButtons(inputId = "subTitlePosition", label = "Subtitle position",
                                                                             choices = list("Left" = "left", "Center" = "center", "Right" = "right"),
                                                                             selected = "left",
                                                                             icon = icon("check"),
                                                                             bigger = TRUE,
                                                                             status = "info",
                                                                             animation = "jelly"),
                                                          numericInput(inputId = "plotTitleSize",
                                                                       label = "The scale of the plot title and subtitle:",
                                                                       value = 12, min = 1, max = Inf, step = 1),
                                                          textInput(inputId = "yAxisTitle", label = "Y axis title:",
                                                                    value = "Features", placeholder = "Write your desired Y axis title"),
                                                          prettyCheckbox(inputId = "show.y.axisGrid", label = "Draw Y axis grid lines", value = TRUE,
                                                                         icon = icon("check"),
                                                                         status = "info")
                                             )
                                      ),
                                      # Output of ExIR visualization
                                      column(8,
                                             column(6,
                                                    downloadButton("download_network_PDF", "Download PDF file", icon = icon("download"),
                                                                   class = "btn-sm btn-block")),
                                             column(6,
                                                    downloadButton("download_network_PNG", "Download PNG file", icon = icon("download"),
                                                                   class = "btn-sm btn-block")
                                             ),
                                             br(),
                                             br(),
                                             br(),
                                             column(12,
                                                    column(4,
                                                           numericInput(inputId = "figure.width", label = "Figure width (in)", value = 8,
                                                                        min = 1, max = Inf, step = 1, width = "100%")
                                                    ),
                                                    column(4,
                                                           numericInput(inputId = "figure.height", label = "Figure height (in)", value = 6,
                                                                        min = 1, max = Inf, step = 1, width = "100%")
                                                    ),
                                                    column(4,
                                                           numericInput(inputId = "PNG.resolution", label = "PNG file resolution (dpi)", value = 300,
                                                                        min = 72, max = Inf, step = 1, width = "100%"))
                                             ),
                                             br(),
                                             br(),
                                             br(),
                                             br(),
                                             plotOutput("ExIR_figure") %>% withSpinner(type = 4)
                                      )
                                  )
                         ),

                         ####*****************************************####

                         # Computational manipulation of cells
                         tabPanel(title = "Computational Cell Manipulation", value = "compManipulate", icon = icon("syringe"),

                                  shinyjs::useShinyjs(),
                                  fluidRow(
                                      column(4,
                                             tags$h5("Options"),
                                             sidebarPanel(width = 12,
                                                          style = "overflow-y:scroll; max-height: 800px; position:relative;",
                                                          ## Input for running compManipulate model

                                                          ### Specify the ko_vertices
                                                          pickerInput(
                                                              inputId = "ko_vertices",
                                                              label = "Select the features to assess their knockout:",
                                                              choices = NULL,
                                                              options = list(
                                                                  style = "btn-sm btn-primary",
                                                                  size = 10,
                                                                  `live-search` = TRUE,
                                                                  `selected-text-format` = "count > 3"),
                                                              multiple = TRUE
                                                          ),

                                                          ### Specify the upregulate_vertices
                                                          pickerInput(
                                                              inputId = "upregulate_vertices",
                                                              label = "Select the features to assess their up-regulation:",
                                                              choices = NULL,
                                                              options = list(
                                                                  style = "btn-sm btn-danger",
                                                                  size = 10,
                                                                  `live-search` = TRUE,
                                                                  `selected-text-format` = "count > 3"),
                                                              multiple = TRUE
                                                          ),
                                                          # no.sim_option
                                                          tags$b(p("Number of simulation runs (corresponding to the SIRIR model):")),
                                                          switchInput(
                                                              value = TRUE,
                                                              inputId = "no.sim_option",
                                                              label = "Automatic",
                                                              labelWidth = "80px"
                                                          ),

                                                          # no.sim
                                                          numericInput(inputId = "no.sim",
                                                                       label = "Determine the numer of simulation runs:",
                                                                       value = NULL),
                                                          # seed
                                                          numericInput(inputId = "seed_forCompMan",
                                                                       label = "Seed number (for random number generation):",
                                                                       value = 1234, min = -Inf, max = Inf, step = 10),
                                                          uiOutput("runCompManipulate")
                                             )),

                                      ## Output computational manipulation
                                      column(8,
                                             panel(tags$h3(tags$b("Computational Manipulation of Cells"),
                                                           style = "margin-top: 0px; font-family: 'sans-serif';"),
                                             tags$p("The computational manipulation of cells works based on the SIRIR
                                                           (SIR-based Influence Ranking) model and could be applied on the output
                                                           of the ExIR model. For feature (gene/protein/etc.) knockout the SIRIR
                                                           model is used to remove the feature from the network and assess its impact
                                                           on the flow of information (signaling) within the network. On the other
                                                           hand, in case of up-regulation a node similar to the desired node is
                                                           added to the network with exactly the same connections (edges) as of the
                                                           original node. Next, the SIRIR model is used to evaluate the difference
                                                           in the flow of information/signaling after adding (up-regulating) the
                                                           desired feature/node compared with the original network. Accordingly,
                                                           you may note that as the gene/protein knockout would impact on the
                                                           integrity of the under-investigation network as well as the networks
                                                           of other overlapping biological processes/pathways, it is recommended to
                                                           select those features that simultaneously have the highest
                                                           (most significant) ExIR-based rank and lowest knockout rank. In
                                                           contrast, as the up-regulation would not affect the integrity of the
                                                           network, you may select the features with highest (most significant)
                                                           ExIR-based and up-regulation-based ranks.",
                                                           style = "text-align: justify; margin-bottom: 0px; font-family:  'Source Sans Pro'; font-size: 16px;"),
                                             status = "success"),
                                             br(),
                                             tags$h4(tags$b("Ranking of the Impact of computational Manipulation of Features on the Network"),
                                                     style = "text-align: center;"),

                                             # Define vertical tabs
                                             verticalTabsetPanel(
                                                 selected = "Knockout Rankings",
                                                 id = "compManipRankings",
                                                 color = "#832C03",
                                                 contentWidth = 10,
                                                 menuSide = "right",

                                                 # Knockout
                                                 verticalTabPanel(title = tags$b("Knockout Rankings"),
                                                                  icon = NULL,
                                                                  box_height = "100px",
                                                                  fluidRow(
                                                                      htmlOutput("nullKnockoutRankings"),
                                                                      column(11,

                                                                             DT::dataTableOutput("knockoutRankingsTable"),
                                                                      ),
                                                                      column(1)
                                                                  )
                                                 ),

                                                 # Up-regulation
                                                 verticalTabPanel(title = tags$b("Up-regulation Rankings"),
                                                                  icon = NULL,
                                                                  box_height = "100px",
                                                                  fluidRow(
                                                                      htmlOutput("nullUpregulationRankings"),
                                                                      column(11,

                                                                             DT::dataTableOutput("upregulationRankingsTable"),
                                                                      ),
                                                                      column(1)
                                                                  )
                                                 ),

                                                 # Combined
                                                 verticalTabPanel(title = tags$b("Combined Rankings"),
                                                                  icon = NULL,
                                                                  box_height = "100px",
                                                                  fluidRow(
                                                                      htmlOutput("nullCombinedRankings"),
                                                                      column(11,

                                                                             DT::dataTableOutput("combinedRankingsTable"),
                                                                      ),
                                                                      column(1)
                                                                  )
                                                 )

                                             ),
                                             fluidRow(tags$h5(br(), ""))
                                      )
                                  )
                         ),

                         ####*****************************************####

                         # Citation
                         tabPanel("Citation", icon = icon("edit"), value = "citationTab",
                                  shinyjs::useShinyjs(),
                                  fluidRow(
                                      column(2),
                                      column(8,
                                             br(),
                                             h4(icon("scroll"), "Please cite the following two papers if you used this shiny app in your study."),
                                             br(),
                                             panel(footer = "",heading = "", status = "primary",
                                      h3(tags$b("ExIR: a versatile one-stop model for the extraction, classification, and prioritization of candidate genes from experimental data"),
                                         style = "color:grey88;"),
                                      br(),
                                      h5(tags$b("The ExIR manuscript is still under review.")),

                                      p("Salavaty A, Ramialison M, Currie PD.", tags$i("ExIR: a versatile one-stop model for the extraction, classification, and prioritization of candidate genes from experimental data"),
                                        style = "text-decoration:underline;")
                                      ),
                                      panel(footer = "",heading = "", status = "success",
                                             h3(tags$b("Integrated Value of Influence: An Integrative Method for the Identification of the Most Influential Nodes within Networks"),
                                                style = "color:grey88;"),
                                             br(),
                                             h5(tags$b("The IVI is published in Patterns, a gold standard data science journal published by the Cell Press.")),

                                             p("Salavaty A, Ramialison M, Currie PD.", tags$i("Integrated Value of Influence: An Integrative Method for the Identification of the Most Influential Nodes within Networks."), "Patterns (N Y). 2020 Jun 22;1(5):100052.",
                                               style = "text-decoration:underline;"),

                                             tags$li(a("DOI: 10.1016/j.patter.2020.100052", href = "https://doi.org/10.1016/j.patter.2020.100052",
                                                       style = "color:blue;")),
                                             tags$li(a("PMID: 33205118", href = "https://pubmed.ncbi.nlm.nih.gov/33205118/",
                                                       style = "color:blue;")),
                                             tags$li(a("PMCID: PMC7660386", href = "http://www.ncbi.nlm.nih.gov/pmc/articles/pmc7660386/",
                                                       style = "color:blue;")),

                                             h4(tags$b("The Bigger Picture of IVI")),
                                             p("Decoding the information buried within the interconnection of components could have
                                               several benefits for the smart control of a complex system. One of the major
                                               challenges in this regard is the identification of the most influential individuals
                                               that have the potential to cause the highest impact on the entire network. This
                                               knowledge could provide the ability to increase network efficiency and reduce costs.
                                               In this article, we present a novel algorithm termed the Integrated Value of Influence
                                               (IVI) that combines the most important topological characteristics of the network to
                                               identify the key individuals within it. The IVI is a versatile method that could
                                               benefit several fields such as sociology, economics, transportation, biology, and
                                               medicine. In biomedical research, for instance, identification of the true influential
                                               nodes within a disease-associated network could lead to the discovery of novel
                                               biomarkers and/or drug targets, a process that could have a considerable impact
                                               on society.", style="text-align: justify")
                                             ),
                                      panel(footer = "",heading = "", status = "danger",
                                             h3(tags$b("R package influential")),
                                             p("The IVI function as well as the centrality-based visualization function are part of the",
                                               tags$code("R package", em(tags$strong("influential"))),
                                               "Additionally, several other functions have been provided for the calculation of some commonly used centrality measures as well as the extraction,
                          classification and ranking of top candidate features from experimental data. You may install the R package influential via either",
                                               a("CRAN", href = "https://cran.r-project.org/package=influential", style = "color:blue;"), "or its",
                                               a("GitHub repo", href = "https://github.com/asalavaty/influential", style = "color:blue;"), ".", style="text-align: justify"),
                                             column(6,
                                                    tags$ul(tags$li(p(tags$b("CRAN:"), tags$pre("install.packages('influential')", style = "color:black;"))),
                                                            tags$li(p(tags$b("GitHub repo:"), tags$pre("## install.packages('devtools')
devtools::install_github('asalavaty/influential',
            build_vignettes = TRUE)", style = "color:black;")))
                                                    )
                                             ),
                                             column(6,
                                                    tags$img(src = 'influential_logo.jpg', align = "left")
                                             )
                                      )
                                  ),
                                  column(2)
                                  )
                         ),

                         ####*****************************************####

                         # About
                         tabPanel("About", icon = icon("address-card"), value = "aboutTab",
                                  shinyjs::useShinyjs(),
                                  fluidRow(
                                      column(2),
                                      column(8,
                                             panel(footer = "",heading = "", status = "default",
                                                   h3(tags$b("Credits"), style = "color:darkcyan"),
                                                   p("The ExIR project was done by", a("Adrian (Abbas) Salavaty", href = "https://www.abbassalavaty.com/", style = "color:blue"),
                                                     "and was supervised by",
                                                     a("Prof. Peter Currie", href = "https://www.armi.org.au/about/our-people/peter-currie/",
                                                       style = "color:blue"),
                                                     "and ",
                                                     a("Assoc. Prof. Mirana Ramialison", href = "https://www.mcri.edu.au/users/mirana-ramialison",
                                                       style = "color:blue"), ". The ExIR shiny app was designed and developed by Adrian according to the ExIR function of the", tags$code("R package influential"),
                                                     ". Also, the visualization of the ExIR results has been rooted from another function of the",
                                                     a("influential R package", href = "https://github.com/asalavaty/influential", style = "color:blue"),
                                                     ". You may have a look at the", em("CITATION"), "tab to get more information regarding the influential R package and how to install it.", style="text-align: justify"),
                                                   br(),
                                                   p("To get more information about the", tags$b("Influential Software Package team"), "refer to the", tags$em("About"), "page of the",  a("Influential Software Package portal", href = "https://influential.erc.monash.edu/", style = "color:blue"), "."),
                                                   br(),
                                                   p("Also, there is a", a("Youtube channel", href = "https://www.youtube.com/playlist?list=PL38ZLo00h-YHu2SbnQ-lfh4iaIsMQ99Qj", style = "color:blue"),
                                                     "dedicated to tutorial videos of different functions of the influential R package.", style="text-align: justify")
                                             ),
                                             panel(footer = "",heading = "", status = "default",
                                                   h3(tags$b("Acknowledgments"), style = "color:darkcyan"),
                                                   tags$ul(
                                                       tags$li(p("We would like to thanks Lan Nguyen, PhD, for his constructive feedback on the ExIR manuscript.", style="text-align: justify")),
                                                       tags$li(p("A part of the results of the ExIR project is based on data generated by the TCGA Research Network.", style="text-align: justify")),
                                                       tags$li(p("The ExIR project was supported by Monash University and The Australian Regenerative Medicine Institute (itself supported by grants from the State Government of Victoria and the Australian Government).", style="text-align: justify"))
                                                   )
                                             )
                                      ),
                                      column(2)
                                  )
                         )

)

####**********************************************####

server <- function(input, output, session,
                   restored_exir_results) {

    # Correct the style of fileInput buttons
    runjs("$('#normExptlData').parent().removeClass('btn-default').addClass('btn-danger');")
    runjs("$('#DifferentialDataset1').parent().removeClass('btn-default').addClass('btn-danger');")
    runjs("$('#DifferentialDataset2').parent().removeClass('btn-default').addClass('btn-danger');")
    runjs("$('#DifferentialDataset3').parent().removeClass('btn-default').addClass('btn-danger');")
    runjs("$('#DifferentialDataset4').parent().removeClass('btn-default').addClass('btn-danger');")
    runjs("$('#DifferentialDataset5').parent().removeClass('btn-default').addClass('btn-danger');")
    runjs("$('#DifferentialDataset6').parent().removeClass('btn-default').addClass('btn-danger');")
    runjs("$('#DifferentialDataset7').parent().removeClass('btn-default').addClass('btn-danger');")
    runjs("$('#DifferentialDataset8').parent().removeClass('btn-default').addClass('btn-danger');")
    runjs("$('#DifferentialDataset9').parent().removeClass('btn-default').addClass('btn-danger');")
    runjs("$('#DifferentialDataset10').parent().removeClass('btn-default').addClass('btn-danger');")

    # Take care of the mtry parameter
    shinyjs::hide("mtry")
    observeEvent(input$mtry_option, {
        shinyjs::reset("mtry")
        if(input$mtry_option == FALSE) {
            shinyjs::show("mtry")
        } else {
            shinyjs::hide("mtry")
        }

    })

    ##***********************##

    # Generate the sample experimental dataset
    output$exprt_sample_table <- renderTable({
        set.seed(1234)
        exprt.sample.table <- matrix(data = c(rnorm(n = 16, mean = 9.5, sd = 0.5)),
                                     nrow = 4, ncol = 4)
        exprt.sample.table <- as.data.frame(exprt.sample.table)
        colnames(exprt.sample.table) <- c(paste0("Gene", c(1:4)))
        rownames(exprt.sample.table) <- c(paste0("Tumor", c(1:2)), paste0("Normal", c(1:2)))
        exprt.sample.table$Condition <- c(rep("Tumor", 2), rep("Normal", 2))
        exprt.sample.table <- as.data.frame(t(as.matrix(exprt.sample.table)))
        exprt.sample.table

    }, rownames = TRUE, hover = TRUE, digits = 2)

    ##***********************##

    # Generate the sample differential dataset 1
    output$differential_table1 <- renderTable({
        set.seed(1234)
        differential_table1 <- matrix(data = c(runif(n = 5, min = -3, max = 4)),
                                     nrow = 5, ncol = 1)
        differential_table1 <- as.data.frame(differential_table1)
        set.seed(1234)
        differential_table1$Padj <- runif(n = 5, min = 0.001, max = 0.02)
        colnames(differential_table1)[1] <- "Log2FC"
        rownames(differential_table1) <- c(paste0("Gene", c(1:5)))
        differential_table1

    }, rownames = TRUE, hover = TRUE, digits = 3)

    ##***********************##

    # Generate the sample differential dataset 2
    output$differential_table2 <- renderTable({
        set.seed(4321)
        differential_table2 <- matrix(data = c(runif(n = 5, min = -5, max = 3)),
                                      nrow = 5, ncol = 1)
        differential_table2 <- as.data.frame(differential_table2)
        set.seed(4321)
        differential_table2$Padj <- runif(n = 5, min = 0.0002, max = 0.04)
        colnames(differential_table2)[1] <- "Log2FC"
        rownames(differential_table2) <- c(paste0("Gene", c(3:7)))
        differential_table2

    }, rownames = TRUE, hover = TRUE, digits = 3)

    ##***********************##

    # Generate the sample regression dataset 1
    output$regression_table1 <- renderTable({
        set.seed(1234)
        regression_table1 <- matrix(data = c(runif(n = 5, min = 0.5, max = 0.9)),
                                      nrow = 5, ncol = 1)
        regression_table1 <- as.data.frame(regression_table1)
        set.seed(1234)
        regression_table1$Padj <- runif(n = 5, min = 0.002, max = 0.03)
        colnames(regression_table1)[1] <- "Log2FC"
        rownames(regression_table1) <- c(paste0("Gene", c(1:5)))
        regression_table1

    }, rownames = TRUE, hover = TRUE, digits = 3)

    ##***********************##

    # Generate the sample regression dataset 2
    output$regression_table2 <- renderTable({
        set.seed(4321)
        regression_table2 <- matrix(data = c(runif(n = 5, min = 0.5, max = 0.8)),
                                    nrow = 5, ncol = 1)
        regression_table2 <- as.data.frame(regression_table2)
        colnames(regression_table2)[1] <- "Log2FC"
        rownames(regression_table2) <- c(paste0("Gene", c(1:5)))
        regression_table2

    }, rownames = TRUE, hover = TRUE, digits = 3)

    ##***********************##

    # Generate the sample desired list of features
    output$desiredFeaturesSampleTable <- renderTable({
        set.seed(4321)
        desiredFeaturesSampleTable <- matrix(data = c(paste0("Gene", sample(x = c(1:1000), size = 5, replace = FALSE))),
                                                      nrow = 5, ncol = 1)
        desiredFeaturesSampleTable

    }, rownames = FALSE, colnames = FALSE, hover = TRUE, digits = 3)

    ##***********************##

    # Generate the real desired list of features
    output$desiredFeaturesRealTable <- renderTable({
        desiredFeaturesRealTable <- matrix(data = c("SERHL2",
                                                    "IFI30",
                                                    "FBXO4",
                                                    "RBM11",
                                                    "ADGRF2"),
                                             nrow = 5, ncol = 1)
        desiredFeaturesRealTable

    }, rownames = FALSE, colnames = FALSE, hover = TRUE, digits = 3)

    ##***********************##

    # Generate the sample synonyms table
    output$synonymsSampleTable <- renderTable({
        synonymsSampleTable <- data.frame(matrix(data = c(paste0("Gene", c(1:5)),
                                                 paste0("Synonym-Gene", c(1:5))),
                                             nrow = 5, ncol = 2))

        colnames(synonymsSampleTable) <- c("Original feature name", "Synonym name")

        synonymsSampleTable

    }, rownames = FALSE, colnames = TRUE, hover = TRUE, digits = 3)

    ##***********************##

    # Generate the real synonyms table
    output$synonymsRealTable <- renderTable({
        synonymsRealTable <- matrix(data = c(
            "ENSG00000183569",
            "ENSG00000216490",
            "ENSG00000151876",
            "ENSG00000185272",
            "ENSG00000164393",

                                             "SERHL2",
                                             "IFI30",
                                             "FBXO4",
                                             "RBM11",
                                             "ADGRF2"
                                             ),
                                           nrow = 5, ncol = 2)
        colnames(synonymsRealTable) <- c("Original (Ensembl) gene name", "Synonym gene name")
        synonymsRealTable

    }, rownames = FALSE, colnames = TRUE, hover = TRUE, digits = 3)

    ####*******************************************************####

    # Get the uploaded data

    ## Define read in batch with progress function
    read_batch_with_progress <- function(file_path, file_name, nrows, no_batches){

        progressSweetAlert(
            size = NULL,
            status = "info",
            striped = TRUE,
            session = session, id = "batchFileUpload",
            title = "Reading in progress ...",
            display_pct = TRUE, value = 1
        )

        seq_length = ceiling(seq.int(from = 2, to = nrows-2, length.out = no_batches+1))
        seq_length = seq_length[-length(seq_length)]

        #read the first line
        ext <- tools::file_ext(file_name)
        df = switch(ext,
                    csv = read.csv(file_path, skip = 0, nrows = 1, row.names = 1, header = TRUE),
                    txt = read.delim(file_path, sep = "\t", skip = 0, nrows = 1, row.names = 1, header = TRUE),
                    validate("Invalid file; Please upload a .csv, or .txt file")
        )
        col_names = colnames(df)

        for(i in seq_along(seq_length)){

            updateProgressBar(
                session = session,
                id = "batchFileUpload",
                value = i*7
            )

            if(i == no_batches) chunk_size = -1 else chunk_size = seq_length[i+1] - seq_length[i]

            df_temp = switch(ext,
                             csv = read.csv(file_path, skip = seq_length[i], nrows = chunk_size,header = FALSE,stringsAsFactors = FALSE, row.names = 1),
                             txt = read.delim(file_path, sep = "\t", skip = seq_length[i], nrows = chunk_size,header = FALSE,stringsAsFactors = FALSE, row.names = 1),
                             validate("Invalid file; Please upload a .csv, or .txt file")
            )
            colnames(df_temp) = col_names
            df = rbind(df,df_temp)
        }

        updateProgressBar(
            session = session,
            id = "batchFileUpload",
            value = 80
        )

        df <- as.data.frame(t(as.matrix(df)))

        updateProgressBar(
            session = session,
            id = "batchFileUpload",
            value = 90
        )

        updateProgressBar(
            session = session,
            id = "batchFileUpload",
            title = "Data transformation in progress ...",
            value = 95
        )

        df[,-ncol(df)] <- as.data.frame(apply(df[,-ncol(df)],
                                              2, as.numeric))

        updateProgressBar(
            session = session,
            id = "batchFileUpload",
            value = 100
        )

        closeSweetAlert(session = session)
        return(df)
    }

    ##############################

    properExpData <- reactiveValues(value = 0)

    ## Check duplications in row names
    observeEvent(input$normExptlData, {

        ext <- tools::file_ext(input$normExptlData$name)
        df4CheckingDuplicates = switch(ext,
                                       csv = read.csv(input$normExptlData$datapath, skip = 0, nrows = 1, row.names = NULL, header = TRUE),
                                       txt = read.delim(input$normExptlData$datapath, sep = "\t", skip = 0, nrows = 1, row.names = NULL, header = TRUE),
                                       validate("Invalid file; Please upload a .csv, or .txt file")
        )

        dfColNumbers <- ncol(df4CheckingDuplicates) - 1

        df4CheckingDuplicates <- switch(ext,
                                        csv = read.csv(input$normExptlData$datapath, skip = 0, row.names = NULL, colClasses = c("character", rep("NULL", dfColNumbers)), header = TRUE),
                                        txt = read.delim(input$normExptlData$datapath, sep = "\t", skip = 0, row.names = NULL, colClasses = c("character", rep("NULL", dfColNumbers)), header = TRUE),
                                        validate("Invalid file; Please upload a .csv, or .txt file")
        )

        if(any(duplicated(df4CheckingDuplicates[,1]))) {
            sendSweetAlert(
                session = session,
                title = "The input Experimental data has duplicated row (feature) names!",
                text = tags$p("Please clean your input experimental data and re-upload that."),
                type = "warning"
            )
            properExpData$value <- 0
        } else {
            properExpData$value <- 1
        }

    })

    ## Experimental data
    experimentalData <- reactive({

        if(properExpData$value == 1) {

        req(input$normExptlData)

        n_rows = length(count.fields(input$normExptlData$datapath))

        inFile <- input$normExptlData

        if (is.null(inFile))
            return(NULL)
        ext <- tools::file_ext(input$normExptlData$name)

        df_out = read_batch_with_progress(file_path = input$normExptlData$datapath,
                                          file_name = input$normExptlData$name,
                                          nrows = n_rows,
                                          no_batches = 10)

        return(df_out)

        }
    })

    ###############

    ## Differential data 1
    differentialDataset1 <- reactive({
        req(input$DifferentialDataset1)

        inFile <- input$DifferentialDataset1

        if (is.null(inFile))
            return(NULL)
        ext <- tools::file_ext(input$DifferentialDataset1$name)
        switch(ext,
               csv = read.csv(input$DifferentialDataset1$datapath, row.names=1),
               txt = read.delim(input$DifferentialDataset1$datapath, sep = "\t", row.names=1),
               validate("Invalid file; Please upload a .csv, or .txt file")
        )
        })

    ###############

    ## Differential dataset 2
    differentialDataset2 <- reactive({
        req(input$DifferentialDataset2)

        inFile <- input$DifferentialDataset2

        if (is.null(inFile))
            return(NULL)
        ext <- tools::file_ext(input$DifferentialDataset2$name)
        switch(ext,
               csv = read.csv(input$DifferentialDataset2$datapath, row.names=1),
               txt = read.delim(input$DifferentialDataset2$datapath, sep = "\t", row.names=1),
               validate("Invalid file; Please upload a .csv, or .txt file")
        )
    })

    ###############

    ## Differential dataset 3
    differentialDataset3 <- reactive({
        req(input$DifferentialDataset3)

        inFile <- input$DifferentialDataset3

        if (is.null(inFile))
            return(NULL)
        ext <- tools::file_ext(input$DifferentialDataset3$name)
        switch(ext,
               csv = read.csv(input$DifferentialDataset3$datapath, row.names=1),
               txt = read.delim(input$DifferentialDataset3$datapath, sep = "\t", row.names=1),
               validate("Invalid file; Please upload a .csv, or .txt file")
        )
    })

    ###############

    ## Differential dataset 4
    differentialDataset4 <- reactive({
        req(input$DifferentialDataset4)

        inFile <- input$DifferentialDataset4

        if (is.null(inFile))
            return(NULL)
        ext <- tools::file_ext(input$DifferentialDataset4$name)
        switch(ext,
               csv = read.csv(input$DifferentialDataset4$datapath, row.names=1),
               txt = read.delim(input$DifferentialDataset4$datapath, sep = "\t", row.names=1),
               validate("Invalid file; Please upload a .csv, or .txt file")
        )
    })

    ###############

    ## Differential dataset 5
    differentialDataset5 <- reactive({
        req(input$DifferentialDataset5)

        inFile <- input$DifferentialDataset5

        if (is.null(inFile))
            return(NULL)
        ext <- tools::file_ext(input$DifferentialDataset5$name)
        switch(ext,
               csv = read.csv(input$DifferentialDataset5$datapath, row.names=1),
               txt = read.delim(input$DifferentialDataset5$datapath, sep = "\t", row.names=1),
               validate("Invalid file; Please upload a .csv, or .txt file")
        )
    })

    ###############

    ## Differential dataset 6
    differentialDataset6 <- reactive({
        req(input$DifferentialDataset6)

        inFile <- input$DifferentialDataset6

        if (is.null(inFile))
            return(NULL)
        ext <- tools::file_ext(input$DifferentialDataset6$name)
        switch(ext,
               csv = read.csv(input$DifferentialDataset6$datapath, row.names=1),
               txt = read.delim(input$DifferentialDataset6$datapath, sep = "\t", row.names=1),
               validate("Invalid file; Please upload a .csv, or .txt file")
        )
    })

    ###############

    ## Differential dataset 7
    differentialDataset7 <- reactive({
        req(input$DifferentialDataset7)

        inFile <- input$DifferentialDataset7

        if (is.null(inFile))
            return(NULL)
        ext <- tools::file_ext(input$DifferentialDataset7$name)
        switch(ext,
               csv = read.csv(input$DifferentialDataset7$datapath, row.names=1),
               txt = read.delim(input$DifferentialDataset7$datapath, sep = "\t", row.names=1),
               validate("Invalid file; Please upload a .csv, or .txt file")
        )
    })

    ###############

    ## Differential dataset 8
    differentialDataset8 <- reactive({
        req(input$DifferentialDataset8)

        inFile <- input$DifferentialDataset8

        if (is.null(inFile))
            return(NULL)
        ext <- tools::file_ext(input$DifferentialDataset8$name)
        switch(ext,
               csv = read.csv(input$DifferentialDataset8$datapath, row.names=1),
               txt = read.delim(input$DifferentialDataset8$datapath, sep = "\t", row.names=1),
               validate("Invalid file; Please upload a .csv, or .txt file")
        )
    })

    ###############

    ## Differential dataset 9
    differentialDataset9 <- reactive({
        req(input$DifferentialDataset9)

        inFile <- input$DifferentialDataset9

        if (is.null(inFile))
            return(NULL)
        ext <- tools::file_ext(input$DifferentialDataset9$name)
        switch(ext,
               csv = read.csv(input$DifferentialDataset9$datapath, row.names=1),
               txt = read.delim(input$DifferentialDataset9$datapath, sep = "\t", row.names=1),
               validate("Invalid file; Please upload a .csv, or .txt file")
        )
    })

    ###############

    ## Differential dataset 10
    differentialDataset10 <- reactive({
        req(input$DifferentialDataset10)

        inFile <- input$DifferentialDataset10

        if (is.null(inFile))
            return(NULL)
        ext <- tools::file_ext(input$DifferentialDataset10$name)
        switch(ext,
               csv = read.csv(input$DifferentialDataset10$datapath, row.names=1),
               txt = read.delim(input$DifferentialDataset10$datapath, sep = "\t", row.names=1),
               validate("Invalid file; Please upload a .csv, or .txt file")
        )
    })

    ##***********************##

    ## Regression dataset
    regressionDataset <- reactive({
        req(input$RegressionDataset)

        inFile <- input$RegressionDataset

        if (is.null(inFile))
            return(NULL)
        ext <- tools::file_ext(input$RegressionDataset$name)
        switch(ext,
               csv = read.csv(input$RegressionDataset$datapath, row.names=1),
               txt = read.delim(input$RegressionDataset$datapath, sep = "\t", row.names=1),
               validate("Invalid file; Please upload a .csv, or .txt file")
        )
    })

    ##***********************##

    ## Desired list of features
    desiredFeatures <- reactive({

        inFile <- input$desiredFeaturesList

        if (is.null(inFile))
            return(NULL)
        ext <- tools::file_ext(input$desiredFeaturesList$name)
        temp_output <-
        switch(ext,
               csv = read.csv(input$desiredFeaturesList$datapath, header = FALSE),
               txt = read.delim(input$desiredFeaturesList$datapath, sep = "\t", header = FALSE),
               validate("Invalid file; Please upload a .csv, or .txt file")
        )

        temp_output <- as.character(temp_output[,1])
    })

    ##***********************##

    # Define the model running button
    output$runExIRbutton <- renderUI({
        if(!is.null(experimentalData()) && !is.null(differentialDataset1())) {
            actionButton("go", "Run the ExIR Model", icon = icon("calculator"),
                         width = "100%", class = "btn-sm btn-primary")
        }
    })

    ##***********************##


    # Define the model running button
    output$jumpToVisTab <- renderUI({
        if(!is.null(experimentalData()) && !is.null(differentialDataset1())) {
            actionButton('goToVisTab', 'Visualize the ExIR results',
                         class = "btn-sm btn-block", icon = icon("pencil-ruler"),
                         width = "100%")
        }
    })

    ##***********************##

    # Take care of ExIR output tables

    ## Hide in the beginning
    # hide("exirOutputTabset")

    ##***********************##

    # Assemble the differential/regression dataset

    ## Gather differential/regression datasets in a list
    assembledDiffReg <- eventReactive(input$go, {

        if(is.null(input$RegressionDataset)) {

        if (input$n_DifferentialDatasets == 1) {
            return(influential::diff_data.assembly(
                differentialDataset1()
                ))
        } else if(input$n_DifferentialDatasets == 2) {
            return(influential::diff_data.assembly(
                differentialDataset1(),
                differentialDataset2()
            ))
        } else if(input$n_DifferentialDatasets == 3) {
            return(influential::diff_data.assembly(
                differentialDataset1(),
                differentialDataset2(),
                differentialDataset3()
            ))
        } else if(input$n_DifferentialDatasets == 4) {
            return(influential::diff_data.assembly(
                differentialDataset1(),
                differentialDataset2(),
                differentialDataset3(),
                differentialDataset4()
            ))
        } else if(input$n_DifferentialDatasets == 5) {
            return(influential::diff_data.assembly(
                differentialDataset1(),
                differentialDataset2(),
                differentialDataset3(),
                differentialDataset4(),
                differentialDataset5()
            ))
        } else if(input$n_DifferentialDatasets == 6) {
            return(influential::diff_data.assembly(
                differentialDataset1(),
                differentialDataset2(),
                differentialDataset3(),
                differentialDataset4(),
                differentialDataset5(),
                differentialDataset6()
            ))
        } else if(input$n_DifferentialDatasets == 7) {
            return(influential::diff_data.assembly(
                differentialDataset1(),
                differentialDataset2(),
                differentialDataset3(),
                differentialDataset4(),
                differentialDataset5(),
                differentialDataset6(),
                differentialDataset7()
            ))
        } else if(input$n_DifferentialDatasets == 8) {
            return(influential::diff_data.assembly(
                differentialDataset1(),
                differentialDataset2(),
                differentialDataset3(),
                differentialDataset4(),
                differentialDataset5(),
                differentialDataset6(),
                differentialDataset7(),
                differentialDataset8()
            ))
        } else if(input$n_DifferentialDatasets == 9) {
            return(influential::diff_data.assembly(
                differentialDataset1(),
                differentialDataset2(),
                differentialDataset3(),
                differentialDataset4(),
                differentialDataset5(),
                differentialDataset6(),
                differentialDataset7(),
                differentialDataset8(),
                differentialDataset9()
            ))
        } else if(input$n_DifferentialDatasets == 10) {
            return(influential::diff_data.assembly(
                differentialDataset1(),
                differentialDataset2(),
                differentialDataset3(),
                differentialDataset4(),
                differentialDataset5(),
                differentialDataset6(),
                differentialDataset7(),
                differentialDataset8(),
                differentialDataset9(),
                differentialDataset10()
            ))
        }
        } else {

            if (input$n_DifferentialDatasets == 1) {
                return(influential::diff_data.assembly(
                    differentialDataset1(),
                    regressionDataset()
                ))
            } else if(input$n_DifferentialDatasets == 2) {
                return(influential::diff_data.assembly(
                    differentialDataset1(),
                    differentialDataset2(),
                    regressionDataset()
                ))
            } else if(input$n_DifferentialDatasets == 3) {
                return(influential::diff_data.assembly(
                    differentialDataset1(),
                    differentialDataset2(),
                    differentialDataset3(),
                    regressionDataset()
                ))
            } else if(input$n_DifferentialDatasets == 4) {
                return(influential::diff_data.assembly(
                    differentialDataset1(),
                    differentialDataset2(),
                    differentialDataset3(),
                    differentialDataset4(),
                    regressionDataset()
                ))
            } else if(input$n_DifferentialDatasets == 5) {
                return(influential::diff_data.assembly(
                    differentialDataset1(),
                    differentialDataset2(),
                    differentialDataset3(),
                    differentialDataset4(),
                    differentialDataset5(),
                    regressionDataset()
                ))
            } else if(input$n_DifferentialDatasets == 6) {
                return(influential::diff_data.assembly(
                    differentialDataset1(),
                    differentialDataset2(),
                    differentialDataset3(),
                    differentialDataset4(),
                    differentialDataset5(),
                    differentialDataset6(),
                    regressionDataset()
                ))
            } else if(input$n_DifferentialDatasets == 7) {
                return(influential::diff_data.assembly(
                    differentialDataset1(),
                    differentialDataset2(),
                    differentialDataset3(),
                    differentialDataset4(),
                    differentialDataset5(),
                    differentialDataset6(),
                    differentialDataset7(),
                    regressionDataset()
                ))
            } else if(input$n_DifferentialDatasets == 8) {
                return(influential::diff_data.assembly(
                    differentialDataset1(),
                    differentialDataset2(),
                    differentialDataset3(),
                    differentialDataset4(),
                    differentialDataset5(),
                    differentialDataset6(),
                    differentialDataset7(),
                    differentialDataset8(),
                    regressionDataset()
                ))
            } else if(input$n_DifferentialDatasets == 9) {
                return(influential::diff_data.assembly(
                    differentialDataset1(),
                    differentialDataset2(),
                    differentialDataset3(),
                    differentialDataset4(),
                    differentialDataset5(),
                    differentialDataset6(),
                    differentialDataset7(),
                    differentialDataset8(),
                    differentialDataset9(),
                    regressionDataset()
                ))
            } else if(input$n_DifferentialDatasets == 10) {
                return(influential::diff_data.assembly(
                    differentialDataset1(),
                    differentialDataset2(),
                    differentialDataset3(),
                    differentialDataset4(),
                    differentialDataset5(),
                    differentialDataset6(),
                    differentialDataset7(),
                    differentialDataset8(),
                    differentialDataset9(),
                    differentialDataset10(),
                    regressionDataset()
                ))
            }
        }
    })

    ## Reset input un running the model
    observeEvent(input$go, {
        reset("normExptlData")
        reset("DifferentialDataset1")
        reset("DifferentialDataset2")
        reset("DifferentialDataset3")
        reset("DifferentialDataset4")
        reset("DifferentialDataset5")
        reset("DifferentialDataset6")
        reset("DifferentialDataset7")
        reset("DifferentialDataset8")
        reset("DifferentialDataset9")
        reset("DifferentialDataset10")
        reset("RegressionDataset")
        reset("desiredFeaturesList")
    })

    ##*********************************************************************##

    # Define the Shiny ExIR function

    shinyExIR <- function(Diff_data, Regr_nCol, Exptl_data, Desired_list, Condition_colname,
                          Normalize, r, max.connections, alpha, num_trees,
                          num_permutations, inf_const, seed) {

        progressSweetAlert(
            size = NULL,
            status = "info",
            striped = TRUE,
            session = session, id = "runShinyExIR",
            title = "Running the RxIR model ...",
            display_pct = TRUE, value = 1
        )

        updateProgressBar(
            session = session,
            title = "Preparing the input data ...",
            id = "runShinyExIR",
            value = 5
        )

        #Define the input arguments

        Diff_data_nCol <- ifelse(is.null(input$RegressionDataset),
                                 ncol(Diff_data),
                                 ncol(Diff_data) - Regr_nCol)

        Diff_value <- seq(1, Diff_data_nCol, 2)

        Regr_value <- if(is.null(input$RegressionDataset)) {
            NULL
        } else {
            Diff_data_nCol + 1
        }

        Sig_value <- if(is.null(input$RegressionDataset)) {
            seq(2, Diff_data_nCol, 2)
        } else {
            seq(2, Diff_data_nCol + Regr_nCol, 2)
        }

        mtry <- if(input$mtry_option == TRUE) {
            NULL
        } else {
            input$mtry
        }

        #change the colnames of Diff_data
        base::colnames(Diff_data) <- base::paste("source", base::colnames(Diff_data),
                                                 sep = ".")

        #change the Inf/-Inf diff values (applicable to sc-Data)
        for(i in 1:base::length(Diff_value)) {

            if(any(base::is.infinite(Diff_data[,Diff_value[i]]))) {

                temp.max.abs.diff.value <-
                    base::max(base::abs(Diff_data[,Diff_value[i]][!base::is.infinite(Diff_data[,Diff_value[i]])]))

                temp.inf.index <- base::which(base::is.infinite(Diff_data[,Diff_value[i]]))

                Diff_data[temp.inf.index, Diff_value[i]] <-
                    base::ifelse(base::unlist(Diff_data[temp.inf.index,Diff_value[i]]) > 0,
                                 temp.max.abs.diff.value*inf_const,
                                 -1*temp.max.abs.diff.value*inf_const)
            }
        }

        # Get the column number of condition column
        condition.index <- match(Condition_colname, colnames(Exptl_data))

        # Transform the condition column to a factor
        Exptl_data[,condition.index] <- base::as.factor(Exptl_data[,condition.index])

        # Normalize the experimental data (if required)
        if(Normalize) {
            Exptl_data[,-condition.index] <- log2(Exptl_data[,-condition.index]+1)
        }

        #***********#

        #1 Calculate differential score

        updateProgressBar(
            session = session,
            title = "Calculating the differential score ...",
            id = "runShinyExIR",
            value = 8
        )

        Diff_data$sum.Diff_value <- base::abs(base::apply(Diff_data[,Diff_value, drop = FALSE],1,sum))
        #range normalize the differential score
        Diff_data$sum.Diff_value <- 1+(((Diff_data$sum.Diff_value-min(Diff_data$sum.Diff_value))*(100-1))/
                                           (max(Diff_data$sum.Diff_value)-min(Diff_data$sum.Diff_value)))

        updateProgressBar(
            session = session,
            title = "Calculating the differential score ...",
            id = "runShinyExIR",
            value = 10
        )

        #***********#

        #2 Calculate regression/time-course R-squared score (if provided)

        updateProgressBar(
            session = session,
            title = "Calculating the regression/time-course R-squared score ...",
            id = "runShinyExIR",
            value = 10
        )

        if (!is.null(Regr_value)) {
            Diff_data$sum.Regr_value <- base::apply(Diff_data[,Regr_value, drop = FALSE],1,sum)
            #range normalize the R-squared score
            Diff_data$sum.Regr_value <- 1+(((Diff_data$sum.Regr_value-min(Diff_data$sum.Regr_value))*(100-1))/
                                               (max(Diff_data$sum.Regr_value)-min(Diff_data$sum.Regr_value)))
        }

        updateProgressBar(
            session = session,
            title = "Calculating the regression/time-course R-squared score ...",
            id = "runShinyExIR",
            value = 15
        )

        #3 Calculate statistical significance of differential/regression factors

        #***********#

        updateProgressBar(
            session = session,
            title = "Calculating the collective statistical significance of differential/regression factors ...",
            id = "runShinyExIR",
            value = 15
        )

        if (max(Diff_data[,Sig_value]) > 1 | min(Diff_data[,Sig_value]) < 0) {
            stop("input Sig-values (p-value/padj) must all be in the range 0 to 1!")
        }

        for(m in 1:length(Sig_value)) {

            if(min(Diff_data[,Sig_value[m]])==0) {

                #range normalize the primitive Sig_value
                temp.min_Sig_value <- base::sort(base::unique(Diff_data[,Sig_value[m]]))[2]

                Diff_data[,Sig_value[m]] <- temp.min_Sig_value+
                    (((Diff_data[,Sig_value[m]]-min(Diff_data[,Sig_value[m]]))*(max(Diff_data[,Sig_value[m]])-temp.min_Sig_value))/
                         (max(Diff_data[,Sig_value[m]])-min(Diff_data[,Sig_value[m]])))
            }
        }

        Diff_data$sum.Sig_value <- base::apply(-log10(Diff_data[,Sig_value, drop = FALSE]),1,sum)
        #range normalize the statistical significance
        Diff_data$sum.Sig_value <- 1+(((Diff_data$sum.Sig_value-min(Diff_data$sum.Sig_value))*(100-1))/
                                          (max(Diff_data$sum.Sig_value)-min(Diff_data$sum.Sig_value)))

        updateProgressBar(
            session = session,
            title = "Calculating the collective statistical significance of differential/regression factors ...",
            id = "runShinyExIR",
            value = 20
        )

        #4 Calculation of the Integrated Value of Influence (IVI)

        #***********#

        #a Separate a transcriptomic profile of diff features

        updateProgressBar(
            session = session,
            title = "Performing the random forests classification (supervised machine learning) ...",
            id = "runShinyExIR",
            value = 20
        )

        if(!is.null(Desired_list)) {
            sig.diff.index <- stats::na.omit(base::unique(base::match(Desired_list,
                                                                      colnames(Exptl_data))))
        } else {
            sig.diff.index <- stats::na.omit(base::unique(base::match(rownames(Diff_data),
                                                                      colnames(Exptl_data))))
        }

        exptl.for.super.learn <- Exptl_data[,sig.diff.index]
        exptl.for.super.learn$condition <- Exptl_data[,condition.index]

        #correct the names of features
        #first preserve a copy of original names
        features.exptl.for.super.learn <- colnames(exptl.for.super.learn)[-ncol(exptl.for.super.learn)]

        colnames(exptl.for.super.learn) <- janitor::make_clean_names(colnames(exptl.for.super.learn))

        #b Perform random forests classification
        base::set.seed(seed = seed)
        rf.diff.exptl <- ranger::ranger(formula = condition ~ .,
                                        data = exptl.for.super.learn,
                                        num.trees = num_trees,
                                        mtry = mtry,
                                        importance = "impurity_corrected",
                                        write.forest = FALSE)

        base::set.seed(seed = seed)
        rf.diff.exptl.pvalue <- as.data.frame(ranger::importance_pvalues(x = rf.diff.exptl,
                                                                         formula = condition ~ .,
                                                                         num.permutations = num_permutations,
                                                                         data = exptl.for.super.learn,
                                                                         method = "altmann"))

        #replace feature names (rownames) with their original names
        rownames(rf.diff.exptl.pvalue) <- features.exptl.for.super.learn

        if(any(is.na(rf.diff.exptl.pvalue[,"pvalue"])) |
           any(is.nan(rf.diff.exptl.pvalue[,"pvalue"]))) {
            rf.diff.exptl.pvalue[c(which(is.na(rf.diff.exptl.pvalue[,"pvalue"])),
                                   which(is.nan(rf.diff.exptl.pvalue[,"pvalue"]))),
                                 "pvalue"] <- 1
        }

        # filtering the RF output data
        select.number <- ifelse(!is.null(Desired_list),
                                round(length(Desired_list)/2),
                                100)

        if(length(which(rf.diff.exptl.pvalue[,"pvalue"] <alpha)) >= select.number) {
            rf.diff.exptl.pvalue <- base::subset(rf.diff.exptl.pvalue, rf.diff.exptl.pvalue$pvalue < alpha)

        } else {
            rf.pval.select <- which(rf.diff.exptl.pvalue[,"pvalue"] <alpha)
            rf.nonSig <- seq(nrow(rf.diff.exptl.pvalue))[-rf.pval.select]
            required.pos.importance <- select.number - length(rf.pval.select)

            temp.rf.diff.exptl.pvalue <- rf.diff.exptl.pvalue[rf.nonSig,]

            rf.importance.select <- utils::tail(order(temp.rf.diff.exptl.pvalue[,"importance"]),
                                                n = required.pos.importance)

            temp.rf.diff.exptl.pvalue <- temp.rf.diff.exptl.pvalue[rf.importance.select,]

            #combine pvalue-based and importance-based tables
            rf.diff.exptl.pvalue <- rbind(rf.diff.exptl.pvalue[rf.pval.select,],
                                          temp.rf.diff.exptl.pvalue)
        }

        # negative importance values could be considered as 0
        if(any(rf.diff.exptl.pvalue[,"importance"] < 0)) {
            rf.diff.exptl.pvalue[which(rf.diff.exptl.pvalue[,"importance"] < 0),
                                 "importance"] <- 0
        }

        # taking care of zero p-values
        if(min(rf.diff.exptl.pvalue[,"pvalue"])==0) {

            #range normalize the primitive pvalue
            temp.min_pvalue <- base::sort(base::unique(rf.diff.exptl.pvalue[,"pvalue"]))[2]

            rf.diff.exptl.pvalue[,"pvalue"] <- temp.min_pvalue+
                (((rf.diff.exptl.pvalue[,"pvalue"]-min(rf.diff.exptl.pvalue[,"pvalue"]))*(max(rf.diff.exptl.pvalue[,"pvalue"])-temp.min_pvalue))/
                     (max(rf.diff.exptl.pvalue[,"pvalue"])-min(rf.diff.exptl.pvalue[,"pvalue"])))
        }

        updateProgressBar(
            session = session,
            title = "Performing the random forests classification (supervised machine learning) ...",
            id = "runShinyExIR",
            value = 35
        )

        #***********#

        #5 Unsupervised machine learning (PCA)

        updateProgressBar(
            session = session,
            title = "Performing PCA (unsupervised machine learning) ...",
            id = "runShinyExIR",
            value = 35
        )

        Exptl_data.for.PCA.index <- stats::na.omit(base::match(base::rownames(rf.diff.exptl.pvalue),
                                                               base::colnames(Exptl_data)))
        temp.Exptl_data.for.PCA <- Exptl_data[,Exptl_data.for.PCA.index]

        temp.PCA <- stats::prcomp(temp.Exptl_data.for.PCA)
        temp.PCA.r <- base::abs(temp.PCA$rotation[,1])

        #range normalize the rotation values
        temp.PCA.r <- 1+(((temp.PCA.r-min(temp.PCA.r))*(100-1))/
                             (max(temp.PCA.r)-min(temp.PCA.r)))

        updateProgressBar(
            session = session,
            title = "Performing PCA (unsupervised machine learning) ...",
            id = "runShinyExIR",
            value = 40
        )

        #***********#

        #6 Performing the first round of association analysis

        updateProgressBar(
            session = session,
            title = "Performing the first round of association analysis ...",
            id = "runShinyExIR",
            value = 40
        )

        #c Performing correlation analysis

        #make ranked experimental data
        Exptl_data[,-condition.index] <- base::t(base::apply(X = Exptl_data[,-condition.index],
                                                             MARGIN = 1, base::rank))

        temp.corr <- coop::pcor(Exptl_data[,-condition.index])

        #reshape the cor matrix
        flt.Corr.Matrix <- function(cormat) {
            ut <- base::upper.tri(cormat)
            data.frame(
                row = base::rownames(cormat)[base::row(cormat)[ut]],
                column = base::rownames(cormat)[base::col(cormat)[ut]],
                cor  = (cormat)[ut]
            )
        }

        temp.corr <- flt.Corr.Matrix(cormat=temp.corr)

        #save a second copy of all cor data
        temp.corr.for.sec.round <- temp.corr

        #filter corr data for only those corr between diff features and themselves/others
        filter.corr.index <- stats::na.omit(base::unique(c(base::which(temp.corr$row %in% rownames(rf.diff.exptl.pvalue)),
                                                           base::which(temp.corr$column %in% rownames(rf.diff.exptl.pvalue)))))
        temp.corr <- temp.corr[filter.corr.index,]

        #filtering low level correlations
        cor.thresh <- r
        temp.corr <- base::subset(temp.corr, temp.corr[,3]>cor.thresh)

        if(nrow(temp.corr)> (max.connections*0.95)) {

            temp.corr.select.index <- utils::tail(order(temp.corr$cor),
                                                  n = round(max.connections*0.95))

            temp.corr <- temp.corr[temp.corr.select.index,]

        }

        diff.only.temp.corr <- temp.corr

        #getting the list of diff features and their correlated features
        diff.plus.corr.features <- base::unique(c(base::as.character(temp.corr[,1]),
                                                  base::as.character(temp.corr[,2])))

        #find the diff features amongst diff.plus.corr.features
        diff.only.features.index <- stats::na.omit(base::unique(base::match(rownames(rf.diff.exptl.pvalue),
                                                                            diff.plus.corr.features)))
        non.diff.only.features <- diff.plus.corr.features[-diff.only.features.index]

        updateProgressBar(
            session = session,
            title = "Performing the first round of association analysis ...",
            id = "runShinyExIR",
            value = 50
        )

        #***********#

        #6 Performing the first round of association analysis

        updateProgressBar(
            session = session,
            title = "Performing the second round of association analysis ...",
            id = "runShinyExIR",
            value = 50
        )

        if(base::length(non.diff.only.features) > 0) {

            #redo correlation analysis
            temp.corr <- temp.corr.for.sec.round
            rm(temp.corr.for.sec.round)

            #filter corr data for only those corr between non.diff.only.features and themselves/others
            filter.corr.index <- stats::na.omit(base::unique(c(base::which(temp.corr$row %in% non.diff.only.features),
                                                               base::which(temp.corr$column %in% non.diff.only.features))))
            temp.corr <- temp.corr[filter.corr.index,]

            #filtering low level correlations
            cor.thresh <- r
            temp.corr <- base::subset(temp.corr, temp.corr[,3]>cor.thresh)

            # separate non.diff.only features
            temp.corr.diff.only.index <- stats::na.omit(base::unique(c(base::which(temp.corr$row %in% rownames(rf.diff.exptl.pvalue)),
                                                                       base::which(temp.corr$column %in% rownames(rf.diff.exptl.pvalue)))))

            if(base::length(temp.corr.diff.only.index)>0) {
                temp.corr <- temp.corr[-temp.corr.diff.only.index,]
            }

            if(nrow(temp.corr)>(max.connections-nrow(diff.only.temp.corr))) {

                temp.corr.select.index <- utils::tail(order(temp.corr$cor),
                                                      n = (max.connections-nrow(diff.only.temp.corr)))

                temp.corr <- temp.corr[temp.corr.select.index,]
            }

            # recombine the diff.only.temp.corr data and temp.corr
            temp.corr <- base::rbind(temp.corr, diff.only.temp.corr)

        } else {
            temp.corr <- diff.only.temp.corr
            rm(temp.corr.for.sec.round, diff.only.temp.corr)
        }

        updateProgressBar(
            session = session,
            title = "Performing the second round of association analysis ...",
            id = "runShinyExIR",
            value = 60
        )

        #***********#

        # Network reconstruction

        updateProgressBar(
            session = session,
            title = "Network reconstruction ...",
            id = "runShinyExIR",
            value = 60
        )

        #d Graph reconstruction
        temp.corr.graph <- igraph::graph_from_data_frame(temp.corr[,c(1:2)])

        updateProgressBar(
            session = session,
            title = "Network reconstruction ...",
            id = "runShinyExIR",
            value = 65
        )

        #***********#

        # Calculation of the integrated value of influence (IVI)

        updateProgressBar(
            session = session,
            title = "Calculation of the integrated value of influence (IVI) ...",
            id = "runShinyExIR",
            value = 65
        )

        #e Calculation of IVI
        temp.corr.ivi <- influential::ivi(temp.corr.graph)

        updateProgressBar(
            session = session,
            title = "Calculation of the integrated value of influence (IVI) ...",
            id = "runShinyExIR",
            value = 70
        )

        #***********#

        # Calculation of the primitive driver score

        updateProgressBar(
            session = session,
            title = "Calculation of the primitive driver score ...",
            id = "runShinyExIR",
            value = 70
        )

        ## Driver score and ranking

        #a calculate first level driver score based on #3 and #4

        Diff_data$IVI <- 0
        Diff_data.IVI.index <- stats::na.omit(match(names(temp.corr.ivi),
                                                    rownames(Diff_data)))

        temp.corr.ivi.for.Diff_data.IVI.index <- stats::na.omit(match(rownames(Diff_data)[Diff_data.IVI.index],
                                                                      names(temp.corr.ivi)))

        Diff_data$IVI[Diff_data.IVI.index] <- temp.corr.ivi[temp.corr.ivi.for.Diff_data.IVI.index]

        #range normalize the IVI
        Diff_data$IVI <- 1+(((Diff_data$IVI-min(Diff_data$IVI))*(100-1))/
                                (max(Diff_data$IVI)-min(Diff_data$IVI)))

        Diff_data$first.Driver.Rank <- 1
        for (i in 1:nrow(Diff_data)) {
            if(c(any(Diff_data[i,Diff_value, drop = FALSE]<0) & any(Diff_data[i,Diff_value, drop = FALSE]>0))) {
                Diff_data$first.Driver.Rank[i] <- 0
            } else {
                Diff_data$first.Driver.Rank[i] <- Diff_data$sum.Sig_value[i]*Diff_data$IVI[i]
            }
        }
        #range normalize the first driver rank
        if(any(Diff_data$first.Driver.Rank == 0)) {
            Diff_data$first.Driver.Rank <- 0+(((Diff_data$first.Driver.Rank-min(Diff_data$first.Driver.Rank))*(100-0))/
                                                  (max(Diff_data$first.Driver.Rank)-min(Diff_data$first.Driver.Rank)))
        } else {
            Diff_data$first.Driver.Rank <- 1+(((Diff_data$first.Driver.Rank-min(Diff_data$first.Driver.Rank))*(100-1))/
                                                  (max(Diff_data$first.Driver.Rank)-min(Diff_data$first.Driver.Rank)))
        }

        updateProgressBar(
            session = session,
            title = "Calculation of the primitive driver score ...",
            id = "runShinyExIR",
            value = 75
        )

        #***********#

        updateProgressBar(
            session = session,
            title = "Calculation of the neighborhood driver score ...",
            id = "runShinyExIR",
            value = 75
        )

        #b (#6) calculate neighborhood score

        #get the list of network nodes
        network.nodes <- base::as.character(igraph::as_ids(V(temp.corr.graph)))

        neighborehood.score.table <- data.frame(node = network.nodes,
                                                N.score = 0)
        for (n in 1:nrow(neighborehood.score.table)) {
            first.neighbors <- igraph::as_ids(igraph::neighbors(graph = temp.corr.graph,
                                                                v = neighborehood.score.table$node[n],
                                                                mode = "all"))
            first.neighbors.index <- stats::na.omit(match(first.neighbors,
                                                          rownames(Diff_data)))

            neighborehood.score.table$N.score[n] <- sum(Diff_data$first.Driver.Rank[first.neighbors.index])
        }

        updateProgressBar(
            session = session,
            title = "Calculation of the neighborhood driver score ...",
            id = "runShinyExIR",
            value = 80
        )

        #***********#

        updateProgressBar(
            session = session,
            title = "Preparation of the driver table ...",
            id = "runShinyExIR",
            value = 80
        )

        Diff_data$N.score <- 0
        Diff_data.N.score.index <- stats::na.omit(match(neighborehood.score.table$node,
                                                        rownames(Diff_data)))

        neighborehood.score.table.for.Diff_data.N.score.index <- stats::na.omit(match(rownames(Diff_data)[Diff_data.N.score.index],
                                                                                      neighborehood.score.table$node))

        Diff_data$N.score[Diff_data.N.score.index] <- neighborehood.score.table$N.score[neighborehood.score.table.for.Diff_data.N.score.index]

        #range normalize (1,100) the neighborhood score
        Diff_data$N.score <- ifelse(sum(Diff_data$N.score) == 0, 1,
                                    1+(((Diff_data$N.score-min(Diff_data$N.score))*(100-1))/
                                           (max(Diff_data$N.score)-min(Diff_data$N.score))))

        #c calculate the final driver score

        Diff_data$final.Driver.score <- (Diff_data$first.Driver.Rank)*(Diff_data$N.score)
        Diff_data$final.Driver.score[Diff_data$final.Driver.score==0] <- NA

        # Create the Drivers table

        Driver.table <- Diff_data

        #remove the rows/features with NA in the final driver score
        Driver.table <- Driver.table[stats::complete.cases(Driver.table),]

        #filter the driver table by the desired list (if provided)
        if(!is.null(Desired_list)) {
            Driver.table.row.index <- stats::na.omit(match(Desired_list,
                                                           rownames(Driver.table)))
            Driver.table <- Driver.table[Driver.table.row.index,]
        }

        if(nrow(as.data.frame(Driver.table))==0) {Driver.table <- NULL} else {

            #range normalize final driver score
            ifelse(length(unique(Driver.table$final.Driver.score)) > 1,
                   Driver.table$final.Driver.score <- 1+(((Driver.table$final.Driver.score-min(Driver.table$final.Driver.score))*(100-1))/
                                                             (max(Driver.table$final.Driver.score)-min(Driver.table$final.Driver.score))),
                   Driver.table$final.Driver.score <- 1)

            #add Z.score
            Driver.table$Z.score <- base::scale(Driver.table$final.Driver.score)

            #add driver rank
            Driver.table$rank <- rank(-Driver.table$final.Driver.score,
                                      ties.method = "min")

            #add P-value
            Driver.table$p.value <- stats::pnorm(Driver.table$Z.score,
                                                 lower.tail = FALSE)

            #add adjusted pvalue
            Driver.table$padj <- stats::p.adjust(p = Driver.table$p.value,
                                                 method = "BH")

            #add driver type
            Driver.table$driver.type <- ""

            for (d in 1:nrow(Driver.table)) {

                if(sum(Driver.table[d,Diff_value])<0) {
                    Driver.table$driver.type[d] <- "Decelerator"

                } else if(sum(Driver.table[d,Diff_value])>0) {
                    Driver.table$driver.type[d] <- "Accelerator"
                } else {
                    Driver.table$driver.type[d] <- NA
                }
            }

            Driver.table <- Driver.table[stats::complete.cases(Driver.table),]

            #remove redundent columns
            Driver.table <- Driver.table[,c("final.Driver.score",
                                            "Z.score",
                                            "rank",
                                            "p.value",
                                            "padj",
                                            "driver.type")]

            #rename column names
            colnames(Driver.table) <- c("Score", "Z.score",
                                        "Rank", "P.value",
                                        "P.adj", "Type")

            #filtering redundant (NaN) results
            Driver.table <- Driver.table[stats::complete.cases(Driver.table),]

            if(nrow(as.data.frame(Driver.table))==0) {Driver.table <- NULL}

        }

        updateProgressBar(
            session = session,
            title = "Preparation of the driver table ...",
            id = "runShinyExIR",
            value = 85
        )

        #***********#

        updateProgressBar(
            session = session,
            title = "Preparation of the biomarker table ...",
            id = "runShinyExIR",
            value = 85
        )

        # Create the Biomarker table

        Biomarker.table <- Diff_data

        #remove the rows/features with NA in the final driver score
        Biomarker.table <- Biomarker.table[stats::complete.cases(Biomarker.table),]

        #filter the biomarker table by the desired list (if provided)
        if(!is.null(Desired_list)) {
            Biomarker.table.row.index <- stats::na.omit(match(Desired_list,
                                                              rownames(Biomarker.table)))
            Biomarker.table <- Biomarker.table[Biomarker.table.row.index,]
        }

        if(nrow(as.data.frame(Biomarker.table))==0) {Biomarker.table <- NULL} else {

            #add RF importance score and p-value
            Biomarker.table$rf.importance <- 0
            Biomarker.table$rf.pvalue <- 1

            Biomarker.table.rf.index <- stats::na.omit(match(rownames(rf.diff.exptl.pvalue),
                                                             rownames(Biomarker.table)))

            rf.for.Biomarker.table <- stats::na.omit(match(rownames(Biomarker.table)[Biomarker.table.rf.index],
                                                           rownames(rf.diff.exptl.pvalue)))

            Biomarker.table$rf.importance[Biomarker.table.rf.index] <- rf.diff.exptl.pvalue$importance[rf.for.Biomarker.table]
            Biomarker.table$rf.pvalue[Biomarker.table.rf.index] <- rf.diff.exptl.pvalue$pvalue[rf.for.Biomarker.table]

            #range normalize rf.importance and rf.pvalue
            Biomarker.table$rf.importance <- 1+(((Biomarker.table$rf.importance-min(Biomarker.table$rf.importance))*(100-1))/
                                                    (max(Biomarker.table$rf.importance)-min(Biomarker.table$rf.importance)))

            Biomarker.table$rf.pvalue <- -log10(Biomarker.table$rf.pvalue)
            Biomarker.table$rf.pvalue <- 1+(((Biomarker.table$rf.pvalue-min(Biomarker.table$rf.pvalue))*(100-1))/
                                                (max(Biomarker.table$rf.pvalue)-min(Biomarker.table$rf.pvalue)))

            #add rotation values
            Biomarker.table$rotation <- 0
            Biomarker.table.for.rotation <- stats::na.omit(match(names(temp.PCA.r),
                                                                 rownames(Biomarker.table)))

            Biomarker.table.rotation.index <- stats::na.omit(match(rownames(Biomarker.table)[Biomarker.table.for.rotation],
                                                                   names(temp.PCA.r)))

            Biomarker.table$rotation[Biomarker.table.for.rotation] <- temp.PCA.r[Biomarker.table.rotation.index]

            #range normalize rotation values
            Biomarker.table$rotation <- 1+(((Biomarker.table$rotation-min(Biomarker.table$rotation))*(100-1))/
                                               (max(Biomarker.table$rotation)-min(Biomarker.table$rotation)))

            #calculate biomarker score
            if(!is.null(Regr_value)) {
                Biomarker.table$final.biomarker.score <- (Biomarker.table$sum.Diff_value)*
                    (Biomarker.table$sum.Regr_value)*(Biomarker.table$sum.Sig_value)*
                    (Biomarker.table$rf.pvalue)*(Biomarker.table$rf.importance)*
                    (Biomarker.table$rotation)
            } else {
                Biomarker.table$final.biomarker.score <- (Biomarker.table$sum.Diff_value)*
                    (Biomarker.table$sum.Sig_value)*
                    (Biomarker.table$rf.pvalue)*(Biomarker.table$rf.importance)*
                    (Biomarker.table$rotation)
            }

            #range normalize biomarker score
            ifelse(length(unique(Biomarker.table$final.biomarker.score)) > 1,
                   Biomarker.table$final.biomarker.score <- 1+(((Biomarker.table$final.biomarker.score-min(Biomarker.table$final.biomarker.score))*(100-1))/
                                                                   (max(Biomarker.table$final.biomarker.score)-min(Biomarker.table$final.biomarker.score))),
                   Biomarker.table$final.biomarker.score <- 1)


            #add biomarker Z.score
            Biomarker.table$Z.score <- base::scale(Biomarker.table$final.biomarker.score)

            #add biomarker rank
            Biomarker.table$rank <- rank(-Biomarker.table$final.biomarker.score, ties.method = "min")

            #add biomarker P-value
            Biomarker.table$P.value <- stats::pnorm(Biomarker.table$Z.score,
                                                    lower.tail = FALSE)

            #add biomarker adjusted p-value
            Biomarker.table$padj <- stats::p.adjust(p = Biomarker.table$P.value,
                                                    method = "BH")

            #add biomarker type
            Biomarker.table$type <- ""

            for (b in 1:nrow(Biomarker.table)) {

                if(sum(Biomarker.table[b,Diff_value])<0) {
                    Biomarker.table$type[b] <- "Down-regulated"

                } else if(sum(Biomarker.table[b,Diff_value])>0) {
                    Biomarker.table$type[b] <- "Up-regulated"
                } else {
                    Biomarker.table$type[b] <- NA
                }
            }

            #remove redundent columns
            Biomarker.table <- Biomarker.table[,c("final.biomarker.score",
                                                  "Z.score", "rank",
                                                  "P.value", "padj", "type")]

            #rename column names
            colnames(Biomarker.table) <- c("Score", "Z.score",
                                           "Rank", "P.value",
                                           "P.adj", "Type")

            #filtering redundant (NaN) results
            Biomarker.table <- Biomarker.table[stats::complete.cases(Biomarker.table),]

            if(nrow(as.data.frame(Biomarker.table))==0) {Biomarker.table <- NULL}

        }

        updateProgressBar(
            session = session,
            title = "Preparation of the biomarker table ...",
            id = "runShinyExIR",
            value = 90
        )

        #***********#

        updateProgressBar(
            session = session,
            title = "Preparation of the DE-mediator table ...",
            id = "runShinyExIR",
            value = 90
        )

        # Create the DE mediators table

        DE.mediator.table <- Diff_data

        #include only rows/features with NA in the final driver score (which are mediators)
        DE.mediator.index <- which(is.na(DE.mediator.table$final.Driver.score))

        DE.mediator.table <- DE.mediator.table[DE.mediator.index,]

        DE.mediator.table$DE.mediator.score <- DE.mediator.table$sum.Sig_value*
            DE.mediator.table$IVI*
            DE.mediator.table$N.score

        #filter the DE mediators table by either the desired list or the list of network nodes

        DE.mediator.row.index <- stats::na.omit(match(network.nodes,
                                                      rownames(DE.mediator.table)))

        if(!is.null(Desired_list)) {
            desired.DE.mediator.row.index <- stats::na.omit(match(Desired_list,
                                                                  rownames(DE.mediator.table)[DE.mediator.row.index]))
            DE.mediator.row.index <- DE.mediator.row.index[desired.DE.mediator.row.index]
        }

        DE.mediator.table <- DE.mediator.table[DE.mediator.row.index,]
        if(nrow(as.data.frame(DE.mediator.table))==0) {DE.mediator.table <- NULL} else {

            #range normalize DE mediators score
            ifelse(length(unique(DE.mediator.table$DE.mediator.score)) > 1,
                   DE.mediator.table$DE.mediator.score <- 1+(((DE.mediator.table$DE.mediator.score-min(DE.mediator.table$DE.mediator.score))*(100-1))/
                                                                 (max(DE.mediator.table$DE.mediator.score)-min(DE.mediator.table$DE.mediator.score))),
                   DE.mediator.table$DE.mediator.score <- 1)

            #add DE mediators Z score
            DE.mediator.table$Z.score <- base::scale(DE.mediator.table$DE.mediator.score)

            #add DE mediators rank
            DE.mediator.table$rank <- rank(-DE.mediator.table$DE.mediator.score, ties.method = "min")

            #add DE mediators P-value
            DE.mediator.table$P.value <- stats::pnorm(DE.mediator.table$Z.score,
                                                      lower.tail = FALSE)

            #add DE mediators adjusted P-value
            DE.mediator.table$padj <- stats::p.adjust(p = DE.mediator.table$P.value,
                                                      method = "BH")

            #remove redundent columns
            DE.mediator.table <- DE.mediator.table[,c("DE.mediator.score", "Z.score",
                                                      "rank", "P.value", "padj")]

            #rename column names
            colnames(DE.mediator.table) <- c("Score", "Z.score",
                                             "Rank", "P.value",
                                             "P.adj")

            #filtering redundant (NaN) results
            DE.mediator.table <- DE.mediator.table[stats::complete.cases(DE.mediator.table),]

            if(nrow(as.data.frame(DE.mediator.table))==0) {DE.mediator.table <- NULL}

        }

        updateProgressBar(
            session = session,
            title = "Preparation of the DE-mediator table ...",
            id = "runShinyExIR",
            value = 95
        )

        #***********#

        updateProgressBar(
            session = session,
            title = "Preparation of the nonDE-mediator table ...",
            id = "runShinyExIR",
            value = 95
        )

        # Create the non-DE mediators table
        non.DE.mediators.index <- stats::na.omit(unique(match(rownames(Diff_data),
                                                              neighborehood.score.table$node)))

        non.DE.mediators.table <- neighborehood.score.table[-c(non.DE.mediators.index),]
        if(nrow(as.data.frame(non.DE.mediators.table))==0) {non.DE.mediators.table <- NULL} else {

            #filter the non-DE mediators table by either the desired list
            if(!is.null(Desired_list)) {
                non.DE.mediator.row.index <- stats::na.omit(match(Desired_list,
                                                                  non.DE.mediators.table$node))
                non.DE.mediators.table <- non.DE.mediators.table[non.DE.mediator.row.index,]
            }
        }

        if(nrow(as.data.frame(non.DE.mediators.table))==0) {non.DE.mediators.table <- NULL} else {

            rownames(non.DE.mediators.table) <- non.DE.mediators.table$node

            non.DE.mediators.ivi.index <- stats::na.omit(match(rownames(non.DE.mediators.table),
                                                               names(temp.corr.ivi)))

            non.DE.mediators.table$ivi <- temp.corr.ivi[non.DE.mediators.ivi.index]

            non.DE.mediators.table$non.DE.mediator.score <- non.DE.mediators.table$N.score*non.DE.mediators.table$ivi

            #range normalize nonDE mediators score
            ifelse(length(unique(non.DE.mediators.table$non.DE.mediator.score)) > 1,
                   non.DE.mediators.table$non.DE.mediator.score <- 1+(((non.DE.mediators.table$non.DE.mediator.score-min(non.DE.mediators.table$non.DE.mediator.score))*(100-1))/
                                                                          (max(non.DE.mediators.table$non.DE.mediator.score)-min(non.DE.mediators.table$non.DE.mediator.score))),
                   non.DE.mediators.table$non.DE.mediator.score <- 1)

            #add non-DE mediators Z.score
            non.DE.mediators.table$Z.score <- base::scale(non.DE.mediators.table$non.DE.mediator.score)

            #add non-DE mediators P-value
            non.DE.mediators.table$P.value <- stats::pnorm(non.DE.mediators.table$Z.score,
                                                           lower.tail = FALSE)

            #add non-DE mediators adjusted p-value
            non.DE.mediators.table$padj <- stats::p.adjust(p = non.DE.mediators.table$P.value,
                                                           method = "BH")

            #add non-DE mediators rank
            non.DE.mediators.table$rank <- rank(-non.DE.mediators.table$non.DE.mediator.score, ties.method = "min")

            #remove redundent columns
            non.DE.mediators.table <- non.DE.mediators.table[,c("non.DE.mediator.score",
                                                                "Z.score", "rank",
                                                                "P.value", "padj")]

            #rename column names
            colnames(non.DE.mediators.table) <- c("Score", "Z.score",
                                                  "Rank", "P.value",
                                                  "P.adj")

            #filtering redundant (NaN) results
            non.DE.mediators.table <- non.DE.mediators.table[stats::complete.cases(non.DE.mediators.table),]

            if(nrow(as.data.frame(non.DE.mediators.table))==0) {non.DE.mediators.table <- NULL}

        }

        updateProgressBar(
            session = session,
            title = "Preparation of the nonDE-mediator table ...",
            id = "runShinyExIR",
            value = 100
        )

        Results <- list("Driver table" = Driver.table,
                        "DE-mediator table" = DE.mediator.table,
                        "nonDE-mediator table" = non.DE.mediators.table,
                        "Biomarker table" = Biomarker.table,
                        "Graph" = temp.corr.graph)

        Results.Non.Null.vector <- vector()

        for (i in 1:length(Results)) {
            Results.Non.Null.vector[i] <- !is.null(Results[[i]])
        }

        Results <- Results[Results.Non.Null.vector]

        # set the class of Results
        base::class(Results) <- "ExIR_Result"

        shinyjs::show("save_for_restore")
        shinyjs::show("closeRestoreNotif")

        closeSweetAlert(session = session)

        return(Results)

    }

        ##***********************##

        # Run the ExIR model

        observeEvent(input$go, {

            if(input$conditionRowname == "") {
                sendSweetAlert(
                    session = session,
                    title = "The condition row name filed is empty!",
                    text = tags$p("Please specify the condition row name of your experimental data file. To get
                                  more information on that click on the PREREQUISITES tab."),
                    type = "warning"
                )
            } else if(input$mtry_option == FALSE && is.na(input$mtry)) {
                sendSweetAlert(
                    session = session,
                    title = "The filed defining the number of features in each node is empty!",
                    text = tags$p("Please either select Automatic or specify the number of features manually."),
                    type = "warning"
                )
            } else if(any(input$conditionRowname == colnames(experimentalData())) == FALSE) {
                sendSweetAlert(
                    session = session,
                    title = "The condition row name was not found!",
                    text = tags$p("The condition row name should match one of the row names of the experimental data.
                                  To get more information on that click on the PREREQUISITES tab."),
                    type = "warning"
                )
            }

        })

        ExIR_results <- eventReactive(input$go, ignoreInit = TRUE, {

            if(input$conditionRowname == "") {
                NULL
            } else if (input$mtry_option == FALSE && is.na(input$mtry)) {
                NULL
            } else if (any(input$conditionRowname == colnames(experimentalData())) == FALSE) {
                NULL
            } else
                shinyExIR(Diff_data = assembledDiffReg(),
                    Regr_nCol = ncol(regressionDataset()),
                    Exptl_data =  experimentalData(),
                    Desired_list = desiredFeatures(),
                    Condition_colname = input$conditionRowname,
                    Normalize = input$Normalize,
                    r = input$correlationCoeff,
                    max.connections = input$max.connections,
                    alpha = input$alpha,
                    num_trees = input$num_trees,
                    num_permutations = input$num_permutations,
                    inf_const = 10^10,
                    seed = input$seed)
        })

        input_ExIR_results <- reactive({
            if(!is.null(restored_exir_results()) && values$upload_state == 'restored_data') {
                restored_exir_results()
            } else {
                ExIR_results()
            }
        })

        reactive({
            class(input_ExIR_results()) <- "ExIR_Result"
        })

    ##***********************##

    # Take care of ExIR outputs

        ## Drivers
        output$driversTable <- DT::renderDataTable(server = FALSE,
            if(!is.null(input_ExIR_results()) && !is.null(input_ExIR_results()$`Driver table`)) {

                input$restore_exir_dataset # Update on file upload

                brks_drivers <- quantile(input_ExIR_results()$`Driver table`$Rank, probs = seq(.05, .95, .05), na.rm = TRUE)
                clrs_drivers <- rev(round(seq(200, 40, length.out = length(brks_drivers) + 1), 0) %>%
                    {paste0("rgb(205,", ., ",", ., ")")})

            DT::datatable(
                { input_ExIR_results()$`Driver table` },

                extensions = c('Buttons', 'SearchPanes', 'FixedHeader'),

                options = list(
                    columnDefs=list(list(targets=1:6, class="dt-right")),
                    paging = TRUE,
                    searching = TRUE,
                    pageLength = 15,
                    fixedColumns = TRUE,
                    autoWidth = TRUE,
                    ordering = TRUE,
                    dom = 'Bfrtip',
                    buttons = list(list(extend = 'csv', filename= paste0("ExIR model-based drivers ", "(", Sys.Date(), ")")),
                                   list(extend = 'excel', filename= paste0("ExIR model-based drivers ", "(", Sys.Date(), ")")),
                                   list(extend = 'print', filename= paste0("ExIR model-based drivers ", "(", Sys.Date(), ")"))
                                   )
                ),

                class = "display"
            ) %>%
                formatRound(columns=c("Score", "Z.score", "P.value", "P.adj"), digits=2) %>%
                    formatRound(columns=c("Rank"), digits=0) %>%
                    formatStyle('Rank', fontWeight = styleInterval(1, c('bold', 'bold'))) %>%
                    formatStyle(
                        'Type',
                        color = styleEqual(
                            c("Accelerator", "Decelerator"), c('#b22222', '#4520bb'))
                        ) %>% formatStyle('Rank', backgroundColor = styleInterval(brks_drivers, clrs_drivers))
            } else {
                NULL
            }
            )

        # Hide the nullDriversTable when the ExIR results are generated
        observe({
            if(!is.null(input_ExIR_results()) && !is.null(input_ExIR_results()$`Driver table`)) {
                shinyjs::hide("nullDriversTable")
            } else {
                shinyjs::show("nullDriversTable")
            }
        })

        # Define the nullDriversTable
        output$nullDriversTable <- renderUI({

            input$restore_exir_dataset # Update on file upload
                    tags$h5(br(),
                            "The ExIR model either has not run yet or didn't find any driver
                            according to your input data and settings!",
                            br(),
                            br(),
                            style = "background-color: lightgrey !important; text-align: center;
                                              padding:12px; margin:48px; border: 1px solid black; border-radius: 12px;")
            })

        #***********#

        ## Biomarkers
        output$biomarkersTable <- DT::renderDataTable(server = FALSE,
            if(!is.null(input_ExIR_results()) && !is.null(input_ExIR_results()$`Biomarker table`)) {

                input$restore_exir_dataset # Update on file upload

                brks_biomarkers <- quantile(input_ExIR_results()$`Biomarker table`$Rank, probs = seq(.05, .95, .05), na.rm = TRUE)
                clrs_biomarkers <- rev(round(seq(200, 40, length.out = length(brks_biomarkers) + 1), 0) %>%
                                {paste0("rgb(205,", ., ",", ., ")")})

                DT::datatable(
                    { input_ExIR_results()$`Biomarker table` },

                    extensions = c('Buttons', 'SearchPanes', 'FixedHeader'),

                    options = list(
                        columnDefs=list(list(targets=1:6, class="dt-right")),
                        paging = TRUE,
                        searching = TRUE,
                        pageLength = 15,
                        fixedColumns = TRUE,
                        autoWidth = TRUE,
                        ordering = TRUE,
                        dom = 'Bfrtip',
                        buttons = list(list(extend = 'csv', filename= paste0("ExIR model-based biomarkers ", "(", Sys.Date(), ")")),
                                       list(extend = 'excel', filename= paste0("ExIR model-based biomarkers ", "(", Sys.Date(), ")")),
                                       list(extend = 'print', filename= paste0("ExIR model-based biomarkers ", "(", Sys.Date(), ")"))
                        )
                    ),

                    class = "display"
                ) %>%
                    formatRound(columns=c("Score", "Z.score", "P.value", "P.adj"), digits=2) %>%
                    formatRound(columns=c("Rank"), digits=0) %>%
                    formatStyle('Rank', fontWeight = styleInterval(1, c('bold', 'bold'))) %>%
                    formatStyle(
                        'Type',
                        color = styleEqual(
                            c("Up-regulated", "Down-regulated"), c('#b22222', '#4520bb'))
                    ) %>% formatStyle('Rank', backgroundColor = styleInterval(brks_biomarkers, clrs_biomarkers))
            } else {
                NULL
            }
        )

        # Hide the nullBiomarkersTable when the ExIR results are generated
        observe({
            if(!is.null(input_ExIR_results()) && !is.null(input_ExIR_results()$`Biomarker table`)) {
                shinyjs::hide("nullBiomarkersTable")
            } else {
                shinyjs::show("nullBiomarkersTable")
            }
        })

        # Define the nullBiomarkersTable
        output$nullBiomarkersTable <- renderUI({

            input$restore_exir_dataset # Update on file upload

            tags$h5(br(),
                    "The ExIR model either has not run yet or didn't find any biomarker
                            according to your input data and settings!",
                    br(),
                    br(),
                    style = "background-color: lightgrey !important; text-align: center;
                                              padding:12px; margin:48px; border: 1px solid black; border-radius: 12px;")
        })

        #***********#

        ## NonDE Mediators
        output$NonDEMediatorsTable <- DT::renderDataTable(server = FALSE,
            if(!is.null(input_ExIR_results()) && !is.null(input_ExIR_results()$`nonDE-mediator table`)) {

                input$restore_exir_dataset # Update on file upload

                brks_nonDE_Mediators <- quantile(input_ExIR_results()$`nonDE-mediator table`$Rank, probs = seq(.05, .95, .05), na.rm = TRUE)
                clrs_nonDE_Mediators <- rev(round(seq(200, 40, length.out = length(brks_nonDE_Mediators) + 1), 0) %>%
                                                {paste0("rgb(205,", ., ",", ., ")")})

                DT::datatable(
                    { input_ExIR_results()$`nonDE-mediator table` },

                    extensions = c('Buttons', 'SearchPanes', 'FixedHeader'),

                    options = list(
                        columnDefs=list(list(targets=1:5, class="dt-right")),
                        paging = TRUE,
                        searching = TRUE,
                        pageLength = 15,
                        fixedColumns = TRUE,
                        autoWidth = TRUE,
                        ordering = TRUE,
                        dom = 'Bfrtip',
                        buttons = list(list(extend = 'csv', filename= paste0("ExIR model-based nonDE-Mediators ", "(", Sys.Date(), ")")),
                                       list(extend = 'excel', filename= paste0("ExIR model-based nonDE-Mediators ", "(", Sys.Date(), ")")),
                                       list(extend = 'print', filename= paste0("ExIR model-based nonDE-Mediators ", "(", Sys.Date(), ")"))
                        )
                    ),

                    class = "display"
                ) %>%
                    formatRound(columns=c("Score", "Z.score", "P.value", "P.adj"), digits=2) %>%
                    formatRound(columns=c("Rank"), digits=0) %>%
                    formatStyle('Rank', fontWeight = styleInterval(1, c('bold', 'bold'))) %>%
                    formatStyle('Rank', backgroundColor = styleInterval(brks_nonDE_Mediators, clrs_nonDE_Mediators))
            } else {
                NULL
            }
        )

        # Hide the nullNonDEMediatorsTable when the ExIR results are generated
        observe({
            if(!is.null(input_ExIR_results()) && !is.null(input_ExIR_results()$`nonDE-mediator table`)) {
                shinyjs::hide("nullNonDEMediatorsTable")
            } else {
                shinyjs::show("nullNonDEMediatorsTable")
            }
        })

        # Define the nullNonDEMediatorsTable
        output$nullNonDEMediatorsTable <- renderUI({

            input$restore_exir_dataset # Update on file upload

            tags$h5(br(),
                    "The ExIR model either has not run yet or didn't find any nonDE-Mediator
                            according to your input data and settings!",
                    br(),
                    br(),
                    style = "background-color: lightgrey !important; text-align: center;
                                              padding:12px; margin:48px; border: 1px solid black; border-radius: 12px;")
        })

        #***********#

        ## DE Mediators
        output$DEMediatorsTable <- DT::renderDataTable(server = FALSE,
            if(!is.null(input_ExIR_results()) && !is.null(input_ExIR_results()$`DE-mediator table`)) {

                input$restore_exir_dataset # Update on file upload

                brks_DE_Mediators <- quantile(input_ExIR_results()$`DE-mediator table`$Rank, probs = seq(.05, .95, .05), na.rm = TRUE)
                clrs_DE_Mediators <- rev(round(seq(200, 40, length.out = length(brks_DE_Mediators) + 1), 0) %>%
                                             {paste0("rgb(205,", ., ",", ., ")")})

                DT::datatable(
                    { input_ExIR_results()$`DE-mediator table` },

                    extensions = c('Buttons', 'SearchPanes', 'FixedHeader'),

                    options = list(
                        columnDefs=list(list(targets=1:5, class="dt-right")),
                        paging = TRUE,
                        searching = TRUE,
                        pageLength = 15,
                        fixedColumns = TRUE,
                        autoWidth = TRUE,
                        ordering = TRUE,
                        dom = 'Bfrtip',
                        buttons = list(list(extend = 'csv', filename= paste0("ExIR model-based DE-Mediators ", "(", Sys.Date(), ")")),
                                       list(extend = 'excel', filename= paste0("ExIR model-based DE-Mediators ", "(", Sys.Date(), ")")),
                                       list(extend = 'print', filename= paste0("ExIR model-based DE-Mediators ", "(", Sys.Date(), ")"))
                        )
                    ),

                    class = "display"
                ) %>%
                    formatRound(columns=c("Score", "Z.score", "P.value", "P.adj"), digits=2) %>%
                    formatRound(columns=c("Rank"), digits=0) %>%
                    formatStyle('Rank', fontWeight = styleInterval(1, c('bold', 'bold'))) %>%
                    formatStyle('Rank', backgroundColor = styleInterval(brks_DE_Mediators, clrs_DE_Mediators))
            } else {
                NULL
            }
        )

        # Hide the nullDEMediatorsTable when the ExIR results are generated
        observe({
            if(!is.null(input_ExIR_results()) && !is.null(input_ExIR_results()$`DE-mediator table`)) {
                shinyjs::hide("nullDEMediatorsTable")
            } else {
                shinyjs::show("nullDEMediatorsTable")
            }
        })

        # Define the nullDEMediatorsTable
        output$nullDEMediatorsTable <- renderUI({

            input$restore_exir_dataset # Update on file upload

            tags$h5(br(),
                    "The ExIR model either has not run yet or didn't find any DE-Mediator
                            according to your input data and settings!",
                    br(),
                    br(),
                    style = "background-color: lightgrey !important; text-align: center;
                                              padding:12px; margin:48px; border: 1px solid black; border-radius: 12px;")
        })


        ####**********************************************####

        # Take care of visualization

        observeEvent(input$goToVisTab, {
            updateTabsetPanel(session, inputId = "inTabset", selected = "ExIRVis")
        })

        ## Synonyms table

            synonymValuesState <- reactiveValues(
            upload_state = 'reseted_names'
            )

        observeEvent(input$synonymsTable, {
            synonymValuesState$upload_state <- 'syn_names'
        })

        observeEvent(input$resetNames, {
            synonymValuesState$upload_state <- 'reseted_names'
        })

        # Reset input file on clicking reset
        observeEvent(input$synonymsTable, {
            reset("resetNames")
        })

        observeEvent(input$resetNames, {
            reset("synonymsTable")

        })

        synonyms_Table <- reactive({

            if(synonymValuesState$upload_state == 'reseted_names') {
                NULL
            } else {

            inFile <- input$synonymsTable

            if (is.null(inFile))
                return(NULL)
            ext <- tools::file_ext(input$synonymsTable$name)
            switch(ext,
                   csv = read.csv(input$synonymsTable$datapath),
                   txt = read.delim(input$synonymsTable$datapath, sep = "\t"),
                   validate("Invalid file; Please upload a .csv, or .txt file")
            )
            }
        })

        ## rendering the output plot
        finalplot <- reactive({

            req(input_ExIR_results())

            influential::exir.vis(exir.results = input_ExIR_results(),
                                  synonyms.table = synonyms_Table(),
                                  n = input$numberOfFeatures,
                                  driver.type = input$driverType,
                                  biomarker.type = input$biomarkerType,
                                  show.drivers = input$showDrivers,
                                  show.biomarkers = input$showBiomarkers,
                                  show.de.mediators = input$showDEMediators,
                                  show.nonDE.mediators = input$showNonDEMediators,
                                  basis = input$selectionBasis,
                                  label.position = input$labelPosition,
                                  nrow = input$exirVis_nrow,
                                  dot.size.min = input$dotSizeMin,
                                  dot.size.max = input$dotSizeMax,
                                  type.color = input$typeColor,
                                  stroke.size = input$strokeSize,
                                  stroke.alpha = input$strokeAlpha,
                                  dot.color.low = input$dotColorLow,
                                  dot.color.high = input$dotColorHigh,
                                  legend.position = input$legendPosition,
                                  legend.direction = input$legendDirection,
                                  legends.layout = input$legendLayout,
                                  boxed.legend = input$boxedLegend,
                                  show.plot.title = input$showPlotTitle,
                                  plot.title = ifelse(input$autoPlotTitle == TRUE, "auto", input$plotTitle),
                                  title.position = input$titlePosition,
                                  plot.title.size = input$plotTitleSize,
                                  show.plot.subtitle = input$showPlotSubTitle,
                                  plot.subtitle = ifelse(input$autoPlotSubTitle == TRUE, "auto", input$plotSubTitle),
                                  subtitle.position = input$subTitlePosition,
                                  y.axis.title = input$yAxisTitle,
                                  show.y.axis.grid = input$show.y.axisGrid
            )
        })


        output$ExIR_figure <- renderPlot({
            finalplot()
        }, res = 96, height = 600)

        observe({
            if (!is.null(finalplot()) == TRUE) {
                # enable the download buttons
                shinyjs::enable("download_network_PDF")
                shinyjs::enable("download_network_PNG")
                shinyjs::enable("figure.width")
                shinyjs::enable("figure.height")
                shinyjs::enable("PNG.resolution")
            }
        })

        # Show/hide visualization options
        observe({
            if(input$showDrivers == FALSE) {
                shinyjs::hide(id = "driverType", anim = TRUE, animType = "slide")
            } else {
                shinyjs::show(id = "driverType", anim = TRUE, animType = "slide")
            }
        })

        observe({
            if(input$showBiomarkers == FALSE) {
                shinyjs::hide(id = "biomarkerType", anim = TRUE, animType = "slide")
            } else {
                shinyjs::show(id = "biomarkerType", anim = TRUE, animType = "slide")
            }
        })

        observe({
            if(input$showPlotTitle == FALSE) {
                shinyjs::hide(id = "autoPlotTitle", anim = TRUE, animType = "slide")
                shinyjs::hide(id = "plotTitle", anim = TRUE, animType = "slide")
                shinyjs::hide(id = "titlePosition", anim = TRUE, animType = "slide")
            } else {
                shinyjs::show(id = "autoPlotTitle", anim = TRUE, animType = "slide")
                shinyjs::show(id = "titlePosition", anim = TRUE, animType = "slide")
            }
        })

        observe({
            if(input$autoPlotTitle == TRUE) {
                shinyjs::hide(id = "plotTitle", anim = TRUE, animType = "slide")
            } else {
                shinyjs::show(id = "plotTitle", anim = TRUE, animType = "slide")
            }
        })

        observe({
            if(input$showPlotSubTitle == FALSE) {
                shinyjs::hide(id = "autoPlotSubTitle", anim = TRUE, animType = "slide")
                shinyjs::hide(id = "plotSubTitle", anim = TRUE, animType = "slide")
                shinyjs::hide(id = "subTitlePosition", anim = TRUE, animType = "slide")
            } else {
                shinyjs::show(id = "autoPlotSubTitle", anim = TRUE, animType = "slide")
                shinyjs::show(id = "subTitlePosition", anim = TRUE, animType = "slide")
            }
        })

        observe({
            if(input$autoPlotSubTitle == TRUE) {
                shinyjs::hide(id = "plotSubTitle", anim = TRUE, animType = "slide")
            } else {
                shinyjs::show(id = "plotSubTitle", anim = TRUE, animType = "slide")
            }
        })

        observe({
            if(input$showPlotTitle == FALSE && input$showPlotSubTitle == FALSE) {
                shinyjs::hide(id = "plotTitleSize", anim = TRUE, animType = "slide")
            } else {
                shinyjs::show(id = "plotTitleSize", anim = TRUE, animType = "slide")
            }
        })

        # Take care of download buttons
        output$download_network_PDF <- downloadHandler(
            filename = paste0("ExIR model-based classified and prioritized features ", "(", Sys.Date(), ").pdf"),
            content = function(file) {
                ggsave(plot = finalplot(), filename = file, device = "pdf",
                       width = input$figure.width, height = input$figure.height, units = "in")
            }
        )

        output$download_network_PNG <- downloadHandler(
            filename = paste0("ExIR model-based classified and prioritized features ", "(", Sys.Date(), ").png"),
            content = function(file) {
                ggsave(plot = finalplot(), filename = file, device = "png",
                       width = input$figure.width, height = input$figure.height, units = "in", dpi = input$PNG.resolution)
            }
        )

        ##********************##
        ##********************##

        # Taking care of different file restoring and uploads

        # Hide the save dataset recommendation and its button until progress is done (progressbar)
        shinyjs::hide("save_for_restore")
        shinyjs::hide("closeRestoreNotif")

        observeEvent(input$closeRestoreNotif, {
            shinyjs::hide("save_for_restore")
            shinyjs::hide("closeRestoreNotif")
        })

        values <- reactiveValues(
            upload_state = NULL
        )

        inputToListen <- reactive({
            list(input$normExptlData ,
                     input$DifferentialDataset1 ,
                     input$DifferentialDataset2 ,
                     input$DifferentialDataset3 ,
                     input$DifferentialDataset4 ,
                     input$DifferentialDataset5 ,
                     input$DifferentialDataset6 ,
                     input$DifferentialDataset7 ,
                     input$DifferentialDataset8 ,
                     input$DifferentialDataset9 ,
                     input$DifferentialDataset10 ,
                     input$RegressionDataset ,
                     input$desiredFeaturesList)
        })

        observeEvent(inputToListen(), {
            values$upload_state <- 'raw_data'
        })

        observeEvent(input$restore_exir_dataset, {
            values$upload_state <- 'restored_data'
        })

        # Reset input file on uploading of another one
        observeEvent(inputToListen(), {
            reset("restore_exir_dataset")
        })

        observeEvent(input$restore_exir_dataset, {
            reset("normExptlData")
            reset("DifferentialDataset1")
            reset("DifferentialDataset2")
            reset("DifferentialDataset3")
            reset("DifferentialDataset4")
            reset("DifferentialDataset5")
            reset("DifferentialDataset6")
            reset("DifferentialDataset7")
            reset("DifferentialDataset8")
            reset("DifferentialDataset9")
            reset("DifferentialDataset10")
            reset("RegressionDataset")
            reset("desiredFeaturesList")

        })

        # Saving and restoring the ExIR results

        # Add restored ExIR reactive values
        restored_exir_results <- reactiveVal()

        # Save ExIR results
        output$save_exir_results <- downloadHandler(
            filename = function() {
                paste0("Session_dataset-", Sys.Date(), ".rds")
            },
            content = function(file) {
                exir_ds <- isolate(ExIR_results())

                save_list <- list(
                    exir_ds = exir_ds
                )
                saveRDS(save_list, file)
            }
        )

        # Reactive restore file
        restore_exir <- reactive({
            validate(
                need(input$restore_exir_dataset, message = FALSE),
                need(tools::file_ext(input$restore_exir_dataset$name) == "rds" ,
                     message = FALSE)
            )
            input$restore_exir_dataset
        })

        observeEvent(input$restore_exir_dataset, {
            if(tools::file_ext(input$restore_exir_dataset$name) != "rds") {
                showNotification("Invalid file; Please upload a .rds file", type = "warning")
            }
        })

        # Reactive to store restored information
        restoreD_exir <- reactive({
            rs_exir <- readRDS(restore_exir()$datapath)
            rs_exir
        })

        # Restore state
        observeEvent(restoreD_exir(), {
            rs_exir <- restoreD_exir()
            restored_exir_results(rs_exir$exir_ds)
        })

        observe({
            if (!is.null(ExIR_results()) == 1) {

                # enable the Save dataset button
                shinyjs::enable("save_exir_results")

                # change the html of the download button
                # shinyjs::html("save_exir_results",
                #               sprintf("<i class='fa fa-download'></i>
                #                   Save dataset (file size: %s)",
                #                       format(object.size(ivi_results()), units = "Kb", digits = 0)
                #               )
                # )
            }
        })

        ##********************##
        ##********************##

        # disable or hide the buttons and inputs on page load
        shinyjs::disable("save_exir_results")
        shinyjs::disable("download_network_PDF")
        shinyjs::disable("download_network_PNG")
        shinyjs::disable("figure.width")
        shinyjs::disable("figure.height")
        shinyjs::disable("PNG.resolution")

        ####**********************************************####

        # Take care of computational manipulation

        # Take care of the no.sim parameter
        shinyjs::hide("no.sim")
        observeEvent(input$no.sim_option, {
            shinyjs::reset("no.sim")
            if(input$no.sim_option == FALSE) {
                shinyjs::show("no.sim")
            } else {
                shinyjs::hide("no.sim")
            }

        })

        # Define the runCompManipulate button
        output$runCompManipulate <- renderUI({
            if(!is.null(input_ExIR_results()) && as.logical(sum(!is.null(input$ko_vertices),
                                                                !is.null(input$upregulate_vertices)))) {
                actionButton("compManipulate", "Run the Manipulation Model", icon = icon("calculator"),
                             width = "100%", class = "btn-sm btn-primary")
            }
        })

        # Take care of inputs

        ## Define the choices of vertices

        observe({
            if(!is.null(input_ExIR_results())) {

                if(is.null(isolate(desiredFeatures()))) {
                    vertexChoices <- igraph::as_ids(V(isolate(input_ExIR_results()$Graph)))
                } else {
                    vertexChoices <- igraph::as_ids(V(isolate(input_ExIR_results()$Graph)))[which(igraph::as_ids(V(isolate(input_ExIR_results()$Graph))) %in% isolate(desiredFeatures()))]
                }

                updatePickerInput(
                    session = session,
                    inputId = "ko_vertices",
                    choices = vertexChoices
                )

                updatePickerInput(
                    session = session,
                    inputId = "upregulate_vertices",
                    choices = vertexChoices
                )
            }
        })

        ##*********************************##

        # Define the shinySirir function

        shinySirir <- function(graph, vertices = V(graph),
                               beta = 0.5, gamma = 1,
                               no.sim = igraph::vcount(graph)*10,  seed = 1234) {

            #Define a data frame
            temp.loocr.table <- data.frame(difference.value = vector("numeric", length = length(vertices)),
                                           rank = vector("integer", length = length(vertices)))

            if(class(vertices) == "character") {
                rownames(temp.loocr.table) <- vertices
            } else if(class(vertices) == "igraph.vs") {
                rownames(temp.loocr.table) <- igraph::as_ids(vertices)
            }


            #Model the spreading based on all nodes
            set.seed(seed)
            all.included.spread <- igraph::sir(graph = graph, beta = beta,
                                               gamma = gamma, no.sim = no.sim)

            #Getting the mean of spread in each independent experiment
            all.mean.spread <- vector("numeric", length = length(all.included.spread))

            for (i in 1:length(all.included.spread)) {
                all.mean.spread[i] <- max(all.included.spread[[i]]$NR)
            }
            all.mean.spread <- mean(all.mean.spread)

            #Model the spread based on leave one out cross ranking (LOOCR)

            for(s in 1:length(vertices)) {

                incProgress(
                    amount = 90/(length(input$ko_vertices) + length(input$upregulate_vertices)),
                    message = NULL,
                    detail = NULL
                )

                temp.graph <- igraph::delete_vertices(graph, unlist(vertices[s]))

                set.seed(seed)

                loocr.spread <- igraph::sir(graph = temp.graph, beta = beta,
                                            gamma = gamma, no.sim = no.sim)

                loocr.mean.spread <- vector("numeric", length = length(loocr.spread))

                for (h in 1:length(loocr.spread)) {
                    loocr.mean.spread[h] <- max(loocr.spread[[h]]$NR)
                }
                loocr.mean.spread <- mean(loocr.mean.spread)
                temp.loocr.table$difference.value[s] <- all.mean.spread - loocr.mean.spread
            }

            temp.loocr.table$rank <- rank(-temp.loocr.table$difference.value, ties.method = "min")

            return(temp.loocr.table)

        }

        ##*********************************##

        # Define the shinyComp_manipulate function

        shinyComp_manipulate <- function(exir_output = NULL,
                                         ko_vertices = NULL,
                                         upregulate_vertices = NULL,
                                         beta = 0.5,
                                         gamma = 1,
                                         no.sim = igraph::vcount(graph) * 10,
                                         seed = 1234) {

            ##**************************##
            # Take care of input graph
            graph <- exir_output$Graph

            ##**************************##
            # Define the progress bar
            withProgress(
                min = 0,
                max = 100,
                value = 10,
                message = "Computational Manipulation of Cells is in progress ...",
                detail = "This may take a while depending on the size of the biological network under-investigation and
                the number of knockout and up-regulation features selected!",
                quoted = FALSE, {

            ##**************************##
            # Over-expression function

            overexpr <- function (graph, vertices = V(graph), beta = 0.5, gamma = 1,
                                  no.sim = igraph::vcount(graph) * 10, seed = 1234)
            {
                temp.loocr.table <- data.frame(difference.value = vector("numeric",
                                                                         length = length(vertices)), rank = vector("integer",
                                                                                                                   length = length(vertices)))
                if (class(vertices) == "character") {
                    rownames(temp.loocr.table) <- vertices
                } else if (class(vertices) == "igraph.vs") {
                    rownames(temp.loocr.table) <- igraph::as_ids(vertices)
                }
                set.seed(seed)
                all.included.spread <- igraph::sir(graph = graph, beta = beta,
                                                   gamma = gamma, no.sim = no.sim)
                all.mean.spread <- vector("numeric", length = length(all.included.spread))
                for (i in 1:length(all.included.spread)) {
                    all.mean.spread[i] <- max(all.included.spread[[i]]$NR)
                }
                all.mean.spread <- mean(all.mean.spread)
                for (s in 1:length(vertices)) {

                    incProgress(
                        amount = 90/(length(input$ko_vertices) + length(input$upregulate_vertices)),
                        message = NULL,
                        detail = NULL
                    )

                    all.edges <- as.data.frame(igraph::as_edgelist(graph))
                    all.desired.edges.index <- unlist(apply(X = all.edges, MARGIN = 2,
                                                            FUN = function(i) which(i %in% rownames(temp.loocr.table)[s])))

                    all.edges <- all.edges[c(all.desired.edges.index),]

                    all.edges.from.indices <- unlist(lapply(X = all.edges[,1],
                                                            FUN = function(i) which(igraph::as_ids(V(graph)) %in% i)))

                    all.edges.to.indices <- unlist(lapply(X = all.edges[,2],
                                                          FUN = function(i) which(igraph::as_ids(V(graph)) %in% i)))

                    all.edges.desired.indices <- data.frame(from = all.edges.from.indices, to = all.edges.to.indices)

                    temp.graph <- igraph::add_vertices(graph = graph, nv = 1,
                                                       name = paste0(rownames(temp.loocr.table)[s], "Fold1"))

                    desired.vertex.index <- match(rownames(temp.loocr.table)[s], igraph::as_ids(V(temp.graph)))

                    desired.vertexFold1.index <- match(paste0(rownames(temp.loocr.table)[s], "Fold1"),
                                                       igraph::as_ids(V(temp.graph)))

                    all.edges.desired.indices[which(all.edges.desired.indices[,1] == desired.vertex.index),1] <- desired.vertexFold1.index
                    all.edges.desired.indices[which(all.edges.desired.indices[,2] == desired.vertex.index),2] <- desired.vertexFold1.index
                    all.edges.desired.indices <- apply(X = all.edges.desired.indices, MARGIN = 1,
                                                       FUN = function(i) paste(i[1], i[2], sep = ","))

                    all.edges.desired.indices <- paste(all.edges.desired.indices, collapse = ",")

                    all.edges.desired.indices <- as.numeric(unlist(strsplit(x = all.edges.desired.indices, split = ",")))

                    temp.graph <- igraph::add_edges(graph = temp.graph, edges = all.edges.desired.indices)

                    set.seed(seed)
                    loocr.spread <- igraph::sir(graph = temp.graph, beta = beta,
                                                gamma = gamma, no.sim = no.sim)
                    loocr.mean.spread <- vector("numeric", length = length(loocr.spread))
                    for (h in 1:length(loocr.spread)) {
                        loocr.mean.spread[h] <- max(loocr.spread[[h]]$NR)
                    }
                    loocr.mean.spread <- mean(loocr.mean.spread)
                    temp.loocr.table$difference.value[s] <- loocr.mean.spread - all.mean.spread
                    loocr.mean.spread
                }
                temp.loocr.table$rank <- rank(-temp.loocr.table$difference.value,
                                              ties.method = "min")
                return(temp.loocr.table)
            }

            ##**************************##
            # Knockout results
            if(!is.null(ko_vertices)) {
                base::suppressWarnings(
                    ko_results <- influential::sirir(
                        graph = graph,
                        vertices = ko_vertices,
                        beta = beta,
                        gamma = gamma,
                        no.sim = no.sim,
                        seed = seed
                    )
                )

                ko_results <- cbind(Feature_name = rownames(ko_results),
                                    ko_results,
                                    Manipulation_type = "Knockout")

                rownames(ko_results) <- NULL
                colnames(ko_results)[3] <- "Rank"

                # correct the orders
                ko_results <- ko_results[order(ko_results[,3]),]

            } else {ko_results <- NULL}

            ##**************************##

            # Over-expression results
            if(!is.null(upregulate_vertices)) {
                base::suppressWarnings(
                    overexpr_results <- overexpr(
                        graph = graph,
                        vertices = upregulate_vertices,
                        beta = beta,
                        gamma = gamma,
                        no.sim = no.sim,
                        seed = seed
                    )
                )

                overexpr_results <- cbind(Feature_name = rownames(overexpr_results),
                                          overexpr_results,
                                          Manipulation_type = "Up-regulation")

                rownames(overexpr_results) <- NULL
                colnames(overexpr_results)[3] <- "Rank"

                # correct the orders
                overexpr_results <- overexpr_results[order(overexpr_results[,3]),]

            } else {overexpr_results <- NULL}

            ##**************************##

            # Combined results
            if(!is.null(ko_results) && !is.null(overexpr_results)) {
                combined_results <- rbind(ko_results,
                                          overexpr_results)

                # Correct the combined ranks
                combined_results[,3] <- rank(-combined_results[,2],
                                             ties.method = "min")

                # correct the orders
                combined_results <- combined_results[order(combined_results[,3]),]

            } else {
                combined_results <- NULL
            }

            ##**************************##

            # Correct results data frames

            ko_results <- ko_results[,-2]
            overexpr_results <- overexpr_results[,-2]
            combined_results <- combined_results[,-2]

            ##**************************##

            # Final results
            final.results <- list(Knockout = ko_results,
                                  Up_regulation = overexpr_results,
                                  Combined = combined_results
            )
            if(sum(sapply(final.results, is.null)) == 3) {
                cat("You should input the name of at least a vertix (feature/gene/etc)
        in the 'ko_vertices' or 'upregulate_vertices' argument.")
            } else {
                non_null_index <- which(sapply(final.results, function(i) (!is.null(i))))
                final.results <- final.results[non_null_index]
                names(final.results) <- names(non_null_index)
                return(final.results)
            }
                }
            )
        }


        ##*********************************##

        ## run the computational manipulation model
        compManiResults <- eventReactive(input$compManipulate, ignoreInit = TRUE, {

            req(input_ExIR_results())

            ExIRGraph <- input_ExIR_results()$Graph

            shinyComp_manipulate(
                exir_output = input_ExIR_results(),
                ko_vertices = input$ko_vertices,
                upregulate_vertices = input$upregulate_vertices,
                beta = 0.5,
                gamma = 1,
                no.sim = ifelse(input$no.sim_option == TRUE,
                                igraph::vcount(ExIRGraph) * 10,
                                input$no.sim),
                seed = input$seed_forCompMan
            )
        })

        ##***********************##

        # Take care of output of computational manipulation

        ## Knockouts
        output$knockoutRankingsTable <- DT::renderDataTable(server = FALSE,
            if(!is.null(compManiResults()) && !is.null(compManiResults()$Knockout)) {

                brks_ko <- quantile(compManiResults()$Knockout$Rank, probs = seq(.05, .95, .05), na.rm = TRUE)
                clrs_ko <- rev(round(seq(200, 40, length.out = length(brks_ko) + 1), 0) %>%
                                        {paste0("rgb(205,", ., ",", ., ")")})

                DT::datatable(
                    { compManiResults()$Knockout },

                    extensions = c('Buttons', 'SearchPanes', 'FixedHeader'),

                    options = list(
                        columnDefs=list(list(targets=1:3, class="dt-right")),
                        paging = TRUE,
                        searching = TRUE,
                        pageLength = 10,
                        fixedColumns = TRUE,
                        autoWidth = TRUE,
                        ordering = TRUE,
                        dom = 'Bfrtip',
                        buttons = list(list(extend = 'csv', filename= paste0("Computational manipulation based knockout rankings ", "(", Sys.Date(), ")")),
                                       list(extend = 'excel', filename= paste0("Computational manipulation based knockout rankings ", "(", Sys.Date(), ")")),
                                       list(extend = 'print', filename= paste0("Computational manipulation based knockout rankings ", "(", Sys.Date(), ")"))
                        )
                    ),

                    class = "display"
                ) %>%
                    formatRound(columns=c("Rank"), digits=0) %>%
                    formatStyle('Rank', fontWeight = styleInterval(1, c('bold', 'bold'))) %>%
                    formatStyle(
                        'Manipulation_type',
                        color = styleEqual(
                            c("Knockout"), c('#4520bb'))
                    ) %>% formatStyle('Rank', backgroundColor = styleInterval(brks_ko, clrs_ko))
            } else {
                NULL
            }
        )

        # Hide the nullKnockoutRankings when the manipulation results are generated
        observe({
            if(!is.null(compManiResults()) && !is.null(compManiResults()$Knockout)) {
                shinyjs::hide("nullKnockoutRankings")
            } else {
                shinyjs::show("nullKnockoutRankings")
            }
        })

        # Define the nullKnockoutRankings
        output$nullKnockoutRankings <- renderUI({

            tags$h5(br(),
                    "The computational manipulation model either has not run yet or no feature has been selected
                    to test its knockout effect!",
                    br(),
                    br(),
                    style = "background-color: lightgrey !important; text-align: center;
                                              padding:12px; margin:48px; border: 1px solid black; border-radius: 12px;")
        })

        #***********#

        ## Up-regulation
        output$upregulationRankingsTable <- DT::renderDataTable(server = FALSE,
            if(!is.null(compManiResults()) && !is.null(compManiResults()$Up_regulation)) {

                brks_upreg <- quantile(compManiResults()$Up_regulation$Rank, probs = seq(.05, .95, .05), na.rm = TRUE)
                clrs_upreg <- rev(round(seq(200, 40, length.out = length(brks_upreg) + 1), 0) %>%
                                   {paste0("rgb(205,", ., ",", ., ")")})

                DT::datatable(
                    { compManiResults()$Up_regulation },

                    extensions = c('Buttons', 'SearchPanes', 'FixedHeader'),

                    options = list(
                        columnDefs=list(list(targets=1:3, class="dt-right")),
                        paging = TRUE,
                        searching = TRUE,
                        pageLength = 10,
                        fixedColumns = TRUE,
                        autoWidth = TRUE,
                        ordering = TRUE,
                        dom = 'Bfrtip',
                        buttons = list(list(extend = 'csv', filename= paste0("Computational manipulation based up-regulation rankings ", "(", Sys.Date(), ")")),
                                       list(extend = 'excel', filename= paste0("Computational manipulation based up-regulation rankings ", "(", Sys.Date(), ")")),
                                       list(extend = 'print', filename= paste0("Computational manipulation based up-regulation rankings ", "(", Sys.Date(), ")"))
                        )
                    ),

                    class = "display"
                ) %>%
                    formatRound(columns=c("Rank"), digits=0) %>%
                    formatStyle('Rank', fontWeight = styleInterval(1, c('bold', 'bold'))) %>%
                    formatStyle(
                        'Manipulation_type',
                        color = styleEqual(
                            c("Up-regulation"), c('#b22222'))
                    ) %>% formatStyle('Rank', backgroundColor = styleInterval(brks_upreg, clrs_upreg))
            } else {
                NULL
            }
        )

        # Hide the nullUpregulationRankings when the manipulation results are generated
        observe({
            if(!is.null(compManiResults()) && !is.null(compManiResults()$Up_regulation)) {
                shinyjs::hide("nullUpregulationRankings")
            } else {
                shinyjs::show("nullUpregulationRankings")
            }
        })

        # Define the nullUpregulationRankings
        output$nullUpregulationRankings <- renderUI({

            tags$h5(br(),
                    "The computational manipulation model either has not run yet or no feature has been selected
                    to test its up-regulation effect!",
                    br(),
                    br(),
                    style = "background-color: lightgrey !important; text-align: center;
                                              padding:12px; margin:48px; border: 1px solid black; border-radius: 12px;")
        })

        #***********#

        ## Combined
        output$combinedRankingsTable <- DT::renderDataTable(server = FALSE,
            if(!is.null(compManiResults()) && !is.null(compManiResults()$Combined)) {

                brks_combined <- quantile(compManiResults()$Combined$Rank, probs = seq(.05, .95, .05), na.rm = TRUE)
                clrs_combined <- rev(round(seq(200, 40, length.out = length(brks_combined) + 1), 0) %>%
                                      {paste0("rgb(205,", ., ",", ., ")")})

                DT::datatable(
                    { compManiResults()$Combined },

                    extensions = c('Buttons', 'SearchPanes', 'FixedHeader'),

                    options = list(
                        columnDefs=list(list(targets=1:3, class="dt-right")),
                        paging = TRUE,
                        searching = TRUE,
                        pageLength = 10,
                        fixedColumns = TRUE,
                        autoWidth = TRUE,
                        ordering = TRUE,
                        dom = 'Bfrtip',
                        buttons = list(list(extend = 'csv', filename= paste0("Computational manipulation based combined rankings ", "(", Sys.Date(), ")")),
                                       list(extend = 'excel', filename= paste0("Computational manipulation based combined rankings ", "(", Sys.Date(), ")")),
                                       list(extend = 'print', filename= paste0("Computational manipulation based combined rankings ", "(", Sys.Date(), ")"))
                        )
                    ),

                    class = "display"
                ) %>%
                    formatRound(columns=c("Rank"), digits=0) %>%
                    formatStyle('Rank', fontWeight = styleInterval(1, c('bold', 'bold'))) %>%
                    formatStyle(
                        'Manipulation_type',
                        color = styleEqual(
                            c("Up-regulation", "Knockout"), c('#b22222', '#4520bb'))
                    ) %>% formatStyle('Rank', backgroundColor = styleInterval(brks_combined, clrs_combined))
            } else {
                NULL
            }
        )

        # Hide the nullCombinedRankings when the manipulation results are generated
        observe({
            if(!is.null(compManiResults()) && !is.null(compManiResults()$Combined)) {
                shinyjs::hide("nullCombinedRankings")
            } else {
                shinyjs::show("nullCombinedRankings")
            }
        })

        # Define the nullCombinedRankings
        output$nullCombinedRankings <- renderUI({

            tags$h5(br(),
                    "The computational manipulation model either has not run yet or only a single class of features
                    has been selected!",
                    br(),
                    br(),
                    style = "background-color: lightgrey !important; text-align: center;
                                              padding:12px; margin:48px; border: 1px solid black; border-radius: 12px;")
        })

        ####*****************************************************####

        if (!interactive()) {
            session$onSessionEnded(function() {
                stopApp()
                q("no")
            })
        }

}

# Run the application
shinyApp(ui = ui, server = server)



