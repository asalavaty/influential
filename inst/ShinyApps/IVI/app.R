# Required packages
library(shiny)
library(shinythemes)
library(shinyWidgets)
library(shinyjs)
library(shinycssloaders)
library(colourpicker)
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

####***********************************************************####

# The App

ui <- navbarPageWithText(id = "inTabset",
                         text = "IVI: An Integrative Method for the Identification of the Most Influential Nodes within Networks",
    theme = shinythemes::shinytheme(theme = "sandstone"),
    titlePanel(title = NULL,
               windowTitle = "Integrated Value of Infuence (IVI)"),
    # IVI calculation
    tabPanel(title = "Calculation of IVI", value = "CalIVI", icon = icon("calculator"),
             tags$head(
               tags$link(
                 rel = "icon",
                 type = "image/x-icon",
                 href = "influential_logo.ico")
             ),
             shinyjs::useShinyjs(),
             fluidRow(
        column(4,
               tags$h5("Calculation of the Integrated Value of Influence (IVI)"),
        sidebarPanel(width = 12,
                     style = "overflow-y:scroll; max-height: 800px; position:relative;",
        ## Input for IVI calculation
        ### Upload the file
        fileInput("file", label = "Upload your data",
                  accept = c(".csv", ".txt", ".tsv", ".sif"),
                  buttonLabel = "Browse"),
               ### Select the type of input file
        prettyRadioButtons("fileType", "Select your file type:",
                            choices = list("Table of node pairs" = "pairs",
                                           "Adjacency matrix" = "adjacency",
                                           "Incidence matrix" = "incidence",
                                           "Cytoscape SIF file" = "sif"),
                           icon = icon("check"),
                           bigger = TRUE,
                           status = "info",
                           animation = "jelly"),
        ## Input for IVI calculation
               ### Specify weightedness
        prettyRadioButtons("weighted", "Weightedness of the network:",
                            choices = list(Unweighted = FALSE, Weighted = TRUE),
                           icon = icon("check"),
                           bigger = TRUE,
                           status = "info",
                           animation = "jelly"),
               ### If require a weighted graph, specify the column number
               numericInput("weightColumn", "Column number of edge weights (if applicable):",
                            value = 0),
               ### Specify directedness
        prettyRadioButtons("directed", "Directedness of the network:",
                            choices = list(Undirected = FALSE, Directed = TRUE),
                            icon = icon("check"),
                            bigger = TRUE,
                            status = "info",
                            animation = "jelly"),
        ### Specify mode
        prettyRadioButtons("mode", "Mode (applicable to directed networks):",
                     choices = list("All edges" = "all",
                                    "Incoming edges" = "in",
                                    "Outgoing edges" = "out"),
                     icon = icon("check"),
                     bigger = TRUE,
                     status = "info",
                     animation = "jelly"),
               ### Loops
        materialSwitch(inputId = "loops", label = "Include the loop edges for IVI calculation",
                      value = TRUE,
                      status = "info"),
        ### If require a weighted graph, specify the column number
        numericInput("d", "Distance value corresponding to the calculation of collective influence (optimal: 3,4):",
                     value = 3, step = 1, min = 1),
               ### Scale
        materialSwitch(inputId = "scaled", label = "Scale the IVI values",
                      value = TRUE,
                      status = "info")
               )),

        ## Output IVI calculation
        column(8,
               tags$h5("After retrieving the results proceed to the",
                       tags$em(tags$b("IVI-based network visualization")), "tab"),
               dataTableOutput("ivi_table"),
               column(4,
                      uiOutput("runIVIbutton"),br()
               ),
               column(4,
                      downloadButton("download_IVI_table", "Download", icon = icon("download"),
                                     width = "100%", class = "btn-info")
                      ),
               column(4,
                      actionButton('jumpToVisTab', 'Visualize the network',
                                   class = "btn-success", icon = icon("pencil-ruler"),
                                   width = "100%")
                      )
               )
        )),
    # IVI visualization
    tabPanel(title = "IVI-based network visualization", value = "NetVis",
             icon = icon("project-diagram"),
             shinyjs::useShinyjs(),
    fluidRow(
        # Input for IVI visualization
        column(4,
               tags$h5("Visualization options"),
               sidebarPanel(width = 12,
                            style = "overflow-y:scroll; max-height: 800px; position:relative;",
                            selectInput(inputId = "layout",
                                        choices = list("KK"="kk", "Star" = "star", "Tree"="tree",
                                                       "Components"="components", "Circle"="circle",
                                                       "Automatic"="automatic", "Grid"="grid",
                                                       "Sphere"="sphere", "Random"="random",
                                                       "DRL"="drl", "FR"="fr", "Gem"="gem",
                                                       "Graphopt"="graphopt", "LGL"="lgl",
                                                       "MDS"="mds"),
                                        label = "Layout", selected = "grid"),
                            selectInput(inputId = "node.color",
                                        choices = list("Magma"="magma", "Inferno" = "inferno",
                                                       "Plasma"="plasma", "Viridis"="viridis",
                                                       "Cividis"="cividis"),
                                        label = "Node fill colors", selected = "viridis"),
                            numericInput(inputId = "node.min.size",
                                         label = "Node minimum size (corresponding to nodes with lowest IVI values)",
                                         value = 3, min = 1, max = Inf),
                            numericInput(inputId = "node.max.size",
                                         label = "Node maximum size (corresponding to nodes with highest IVI values)",
                                         value = 15, min = 2, max = Inf),
                            numericInput(inputId = "dist.power",
                                         label = "Distinction power",
                                         value = 1, min = 1, max = 5),
                            selectInput(inputId = "node.shape", label = "Node shape",
                                        choices = list("Circle" = "circle", "Square" = "square",
                                                       "Diamond"="diamond", "Triangle"="triangle",
                                                       "Inverted triangle"="inverted triangle"),
                                        selected = "circle"),
                            numericInput(inputId = "stroke.size",
                                         label = "Stroke (node border) size",
                                         value = 1.5, min = 0, max = 5, step = 0.5),
                            prettyRadioButtons(inputId = "stroke.color", label = "Stroke (node border) color:",
                                         choices = list("Identical to node color" = "identical",
                                                        "Select desired color" = NULL),
                                         icon = icon("check"),
                                         bigger = TRUE,
                                         status = "info",
                                         animation = "jelly"),
                            colourInput(inputId = "pick.stroke.color",
                                        value =  "black", label = NULL, returnName = TRUE,
                                        closeOnClick = TRUE, allowTransparent = TRUE),
                            sliderInput(inputId = "stroke.alpha", label = "Stroke opacity (alpha)",
                                        min = 0, max = 1, value = 0.6, step = 0.1,
                                        ticks = TRUE, animate = FALSE),
                            prettyCheckbox(inputId = "show.labels", label = "Show node labels",
                                          value = TRUE,
                                          icon = icon("check"),
                                          status = "info"),
                            sliderInput(inputId = "label.cex", label = "Relative label size",
                                        min = 0.1, max = 2, value = 0.4, step = 0.2,
                                        ticks = TRUE, animate = FALSE),
                            colourInput(inputId = "label.color",
                                        value =  "black", label = "Label color", returnName = TRUE,
                                        closeOnClick = TRUE, allowTransparent = TRUE),
                            prettyCheckbox(inputId = "directed.edges", label = "Visualize the network as directed (if applicable)",
                                          value = FALSE,
                                          icon = icon("check"),
                                          status = "info"),
                            sliderInput(inputId = "arrow.width", label = "Arrow width (applicable to directed networks)",
                                        min = 1, max = 100, value = 25, step = 5,
                                        ticks = TRUE, animate = FALSE),
                            sliderInput(inputId = "arrow.length", label = "Arrow length (applicable to directed networks)",
                                        min = 0.01, max = 1, value = 0.07, step = 0.05,
                                        ticks = TRUE, animate = FALSE),
                            sliderInput(inputId = "edge.width", label = "Edge width",
                                        min = 0.1, max = 5, value = 0.5, step = 0.5,
                                        ticks = TRUE, animate = FALSE),
                            prettyCheckbox(inputId = "weighted.edges", label = "Visualize the network as weighted (if applicable)",
                                          value = FALSE,
                                          icon = icon("check"),
                                          status = "info"),
                            sliderInput(inputId = "edge.width.min.max", label = "Relative min and max of edges",
                                        min = 0.1, max = 5, value = c(0.2, 1), step = 0.5,
                                        dragRange = TRUE, ticks = TRUE, animate = FALSE),
                            colourInput(inputId = "edge.color",
                                        value =  "grey75", label = "Edge color", returnName = TRUE,
                                        closeOnClick = TRUE, allowTransparent = TRUE),
                            selectInput(inputId = "edge.linetype", label = "Edge linetype",
                                        choices = list("Two-dash" = "twodash", "Long-dash" = "longdash",
                                                       "Dot-dash" = "dotdash", "Dotted" = "dotted",
                                                       "Dashed" = "dashed", "Solid" = "solid"),
                                        selected = "solid"),
                            prettyRadioButtons(inputId = "legend.position", label = "Legend position",
                                         choices = list("Left" = "left", "Right" = "right", "Bottom" = "bottom",
                                                        "Top" = "top", "Remove legend" = "none"),
                                         icon = icon("check"),
                                         bigger = TRUE,
                                         status = "info",
                                         animation = "jelly",
                                         selected = "right"),
                            prettyRadioButtons(inputId = "legend.direction", label = "Legend direction",
                                         choices = list("Horizontal" = "horizontal", "Vertical" = "vertical"),
                                         selected = "vertical",
                                         icon = icon("check"),
                                         bigger = TRUE,
                                         status = "info",
                                         animation = "jelly"),
                            textInput(inputId = "legend.title", label = "Legend title",
                                      value = "IVI value", placeholder = "Write your desired legend title"),
                            prettyCheckbox(inputId = "boxed.legend", label = "Boxed legend", value = TRUE,
                                           icon = icon("check"),
                                           status = "info"),
                            prettyCheckbox(inputId = "show.plot.title", label = "Show plot title", value = TRUE,
                                           icon = icon("check"),
                                           status = "info"),
                            textInput(inputId = "plot.title", label = "Plot title",
                                      value = "IVI-based Network", placeholder = "Write your desired plot title"),
                            prettyRadioButtons(inputId = "title.position", label = "Title position",
                                         choices = list("Left" = "left", "Center" = "center", "Right" = "right"),
                                         selected = "center",
                                         icon = icon("check"),
                                         bigger = TRUE,
                                         status = "info",
                                         animation = "jelly"),
                            prettyCheckbox(inputId = "show.bottom.border", label = "Show bottom border line", value = TRUE,
                                          icon = icon("check"),
                                          status = "info"),
                            prettyCheckbox(inputId = "show.left.border", label = "Show left border line", value = TRUE,
                                          icon = icon("check"),
                                          status = "info")
               )
               ),
        # Output of IVI visualization
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
                   plotOutput("network_figure") %>% withSpinner(type = 4)
        )
        )
    ),

    ####*****************************************####

    # Citation
    tabPanel("Citation", icon = icon("edit"), value = "citationTab",
             shinyjs::useShinyjs(),
             fluidRow(
               column(7,
                      h2(tags$b("Integrated Value of Influence: An Integrative Method for the Identification of the Most Influential Nodes within Networks"),
                         style = "color:darkcyan;"),
                      br(),
                      h5(tags$b("Please use the corresponding paper of the IVI if you used this shiny app in your study. The IVI is published in Patterns, a gold standard data science journal published by the Cell Press.")),

                      p("Salavaty A, Ramialison M, Currie PD.", tags$i("Integrated Value of Influence: An Integrative Method for the Identification of the Most Influential Nodes within Networks."), "Patterns (N Y). 2020 Jun 22;1(5):100052.",
                        style = "text-decoration:underline;"),

                      tags$li(a("DOI: 10.1016/j.patter.2020.100052", href = "https://doi.org/10.1016/j.patter.2020.100052",
                                style = "color:blue;")),
                      tags$li(a("PMID: 33205118", href = "https://pubmed.ncbi.nlm.nih.gov/33205118/",
                                style = "color:blue;")),
                      tags$li(a("PMCID: PMC7660386", href = "http://www.ncbi.nlm.nih.gov/pmc/articles/pmc7660386/",
                                style = "color:blue;")),

                      h3(tags$b("The Bigger Picture")),
                      p("Decoding the information buried within the interconnection of components could have several
                                             benefits for the smart control of a complex system. One of the major challenges in this regard
                                             is the identification of the most influential individuals that have the potential to cause the
                                             highest impact on the entire network. This knowledge could provide the ability to increase network
                                             efficiency and reduce costs. In this article, we present a novel algorithm termed the Integrated
                                             Value of Influence (IVI) that combines the most important topological characteristics of the
                                             network to identify the key individuals within it. The IVI is a versatile method that could
                                             benefit several fields such as sociology, economics, transportation, biology, and medicine.
                                             In biomedical research, for instance, identification of the true influential nodes within a
                                             disease-associated network could lead to the discovery of novel biomarkers and/or drug targets,
                                             a process that could have a considerable impact on society.", style="text-align: justify"),
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
               ),
               column(5,
                      img(src='IVI.Patterns.cover.jpg', align = "right")
               )
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
                            p("The IVI project was done by", a("Adrian (Abbas) Salavaty", href = "https://www.abbassalavaty.com/", style = "color:blue"),
                              "and was supervised by",
                              a("Prof. Peter Currie", href = "https://www.armi.org.au/about/our-people/peter-currie/",
                                style = "color:blue"),
                              "and ",
                              a("Assoc. Prof. Mirana Ramialison", href = "https://www.mcri.edu.au/users/mirana-ramialison",
                                style = "color:blue"), ". The IVI shiny app was designed and developed by Adrian according to the IVI function of the", tags$code("R package influential"),
                              ". Also, the visualization of the networks based on IVI values has been rooted from another function of the",
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
                              tags$li(p("We would like to thanks Ehsan Rezaei-Darzi (MSc Biostatistics) for his consultations regarding proper use of statistical methods and evaluations in the IVI project.
                          Also, thanks to the constructive feedback on the IVI manuscript by Hieu Tri Nim.", style="text-align: justify")),
                              tags$li(p("A part of the results of the IVI project is based on data generated by the TCGA Research Network.", style="text-align: justify")),
                              tags$li(p("The IVI project was supported by Monash University and The Australian Regenerative Medicine Institute (itself supported by grants from the State Government of Victoria and the Australian Government).", style="text-align: justify"))
                            )
                      )
               ),
               column(2)
             )
    )

    )

####**********************************************####

server <- function(input, output, session) {

  # Correct the style of fileInput button
  runjs("$('#file').parent().removeClass('btn-default').addClass('btn-danger');")

    # Define the calculation button
    output$runIVIbutton <- renderUI({
        if(!is.null(data())) {
            actionButton("go", "Calculate IVI", icon = icon("calculator"),
                         width = "100%", class = "btn-primary")
        }
    })

    # Get the uploaded data

    data <- reactive({
        inFile <- input$file

        if (is.null(inFile))
            return(NULL)
        ext <- tools::file_ext(input$file$name)
        switch(ext,
               csv = read.csv(input$file$datapath, sep = ","),
               tsv = read.delim(input$file$datapath, sep = "\t"),
               txt = read.delim(input$file$datapath, sep = "\t"),
               sif = read.delim(input$file$datapath, header = FALSE, sep = "\t")[,c(1,3)],
               validate("Invalid file; Please upload a .csv, .tsv, .txt, or .sif file")
        )
    })

    ####**********************************************####

        # Reconstruct the graph

    final.graph <- eventReactive(input$go, {
        req(input$file)

        fileType <- input$fileType
        temp.graph <-
            switch(fileType,
                   pairs = igraph::graph_from_data_frame(d = data(), directed = input$directed),
                   adjacency = igraph::graph_from_adjacency_matrix(adjmatrix = data(), mode = input$mode, weighted = input$weighted),
                   incidence = igraph::graph_from_incidence_matrix(incidence = data(), mode = input$mode, directed = input$directed, weighted = input$weighted),
                   sif = igraph::graph_from_data_frame(d = data(), directed = input$directed))

        temp.graph <-
            set.edge.attribute(graph = temp.graph,
                               name = "weight",
                               index= E(temp.graph),
                               value = data()[,input$weightColumn])
        return(temp.graph)
    })

    observe({
        if(input$weighted == FALSE || input$weightColumn == 0) {

            final.graph <- eventReactive(input$go, {
                req(input$file)

                fileType <- input$fileType
                switch(fileType,
                       pairs = igraph::graph_from_data_frame(d = data(), directed = input$directed),
                       adjacency = igraph::graph_from_adjacency_matrix(adjmatrix = data(), mode = input$mode, weighted = input$weighted),
                       incidence = igraph::graph_from_incidence_matrix(incidence = data(), mode = input$mode, directed = input$directed, weighted = input$weighted),
                       sif = igraph::graph_from_data_frame(d = data(), directed = input$directed))
            })
        }
    })

    observeEvent(input$weighted, {
        reset("weightColumn")
    })

        observeEvent(input$fileType, {
            if(input$fileType != "pairs") {
                hide("weightColumn", anim = TRUE, animType = "slide")
                reset("weightColumn")
            } else {
              shinyjs::show("weightColumn", anim = TRUE, animType = "slide")
            }
        })

        observeEvent(input$weighted, {
            if(input$weighted == FALSE) {
                disable("weightColumn")
            } else {
                enable("weightColumn")
            }
        })

        ####**********************************************####

        # Calculation of IVI

        observeEvent(input$go, {
            sendSweetAlert(
                session = session,
                title = "Network reconstruction and IVI calculation in progress!",
                text = tags$p("This may take a while depending on the size of the dataset..."),
                type = "info"
            )
        })

        ivi_results <- eventReactive(input$go, {
            req(final.graph())

            withProgress(message = "Progress ...",
                         detail = NULL,
                         value = 65, min = 1, max = 100, {

                ivi(graph = final.graph(), weights = NULL, directed = input$directed,
                    mode = input$mode, loops = input$loops, d = input$d, scaled = input$scaled)
            })
        })

        sorted_ivi_results <- eventReactive(input$go, {
            req(input$file)
            rev(sort(ivi_results()))
        })

        final.IVI.table <- eventReactive(input$go, {
            data.frame(Node_name = names(sorted_ivi_results()),
                       IVI_value = round(sorted_ivi_results(), digits = 2))
        })

    output$ivi_table <- renderDataTable({
        input$go # Update on file upload
        final.IVI.table()
    },
    options = list(pageLength = 15))

    observe({
        if (!is.null(sorted_ivi_results()) == 1) {
            # enable the download button
            shinyjs::enable("download_IVI_table")
            # change the html of the download button
            shinyjs::html("download_IVI_table",
                          sprintf("<i class='fa fa-download'></i>
                              Download (file size: %s)",
                                  format(object.size(sorted_ivi_results()), units = "Kb", digits = 0)
                          )
            )
            # enable the visualization button
            shinyjs::enable("jumpToVisTab")
        }
    })


    output$download_IVI_table <- downloadHandler(
        filename = function() {
            paste0(gsub(pattern = paste0(".", tools::file_ext(fixUploadedFilesNames(input$file)$name)),
                        replacement = "", x = fixUploadedFilesNames(input$file)$name, ignore.case = TRUE),
                   "_IVI_values", ".tsv")
        },
        content = function(file) {
            write.table(x = data.frame(Node_name = names(sorted_ivi_results()),
                                       IVI_value = sorted_ivi_results()),
                        file = file, row.names = FALSE,
                        quote = FALSE, sep = "\t")
        }
    )

    # disable the download button on page load
    shinyjs::disable("download_IVI_table")

    ####*****************************************************####

    # rendering the output plot

        finalplot <- reactive({

        observeEvent(input$stroke.color, {
            if(input$stroke.color == "identical"){
                hide("pick.stroke.color")
                disable("pick.stroke.color")
            }else{
                enable("pick.stroke.color")
              shinyjs::show("pick.stroke.color", anim = TRUE, animType = "fade")
            }
        })

        observe({
            if(!is.weighted(final.graph()) | input$weighted == FALSE) {
                disable("weighted.edges")
            } else if(is.weighted(final.graph()) | input$weighted == TRUE) {
                enable("weighted.edges")
            }
        })


                    cent_network.vis(graph = final.graph(),
                                     cent.metric = ivi_results(),
                                     layout = input$layout,
                                     node.color = input$node.color,
                                     node.size.min = input$node.min.size,
                                     node.size.max = input$node.max.size,
                                     dist.power = input$dist.power,
                                     node.shape = input$node.shape,
                                     stroke.size = input$stroke.size,
                                     stroke.color = ifelse(input$stroke.color == "identical",
                                                           input$stroke.color, input$pick.stroke.color),
                                     stroke.alpha = input$stroke.alpha,
                                     show.labels = input$show.labels,
                                     label.cex = input$label.cex,
                                     label.color = input$label.color,
                                     directed = input$directed.edges,
                                     arrow.width = input$arrow.width,
                                     arrow.length = input$arrow.length,
                                     edge.width = input$edge.width,
                                     weighted = input$weighted.edges,
                                     edge.width.min = input$edge.width.min.max[1],
                                     edge.width.max = input$edge.width.min.max[2],
                                     edge.color = input$edge.color,
                                     edge.linetype = input$edge.linetype,
                                     legend.position = input$legend.position,
                                     legend.direction = input$legend.direction,
                                     legend.title = input$legend.title,
                                     boxed.legend = input$boxed.legend,
                                     show.plot.title = input$show.plot.title,
                                     plot.title = input$plot.title,
                                     title.position = input$title.position,
                                     show.bottom.border = input$show.bottom.border,
                                     show.left.border = input$show.left.border,
                                     seed = 1234)
        })

    observe({
        if(input$weighted.edges == FALSE) {
            hide("edge.width.min.max")
        } else {
          shinyjs::show("edge.width.min.max", anim = TRUE, animType = "slide")
        }
    })

    output$network_figure <- renderPlot({
        req(final.graph())
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


    observeEvent(input$jumpToVisTab, {
        updateTabsetPanel(session, inputId = "inTabset", selected = "NetVis")
    })

    observe({
        if(input$directed.edges == TRUE) {
          shinyjs::show(id = "arrow.width", anim = TRUE, animType = "slide")
          shinyjs::show(id = "arrow.length", anim = TRUE, animType = "slide")
        } else {
          shinyjs::hide(id = "arrow.width", anim = TRUE, animType = "slide")
          shinyjs::hide(id = "arrow.length", anim = TRUE, animType = "slide")
        }
    })

    observe({
        if(input$show.plot.title == FALSE) {
          shinyjs::hide(id = "plot.title", anim = TRUE, animType = "slide")
            hshinyjs::hide(id = "title.position", anim = TRUE, animType = "slide")
        } else {
          shinyjs::show(id = "plot.title", anim = TRUE, animType = "slide")
          shinyjs::show(id = "title.position", anim = TRUE, animType = "slide")
        }
    })

    output$download_network_PDF <- downloadHandler(
        filename = function() {
            paste0(gsub(pattern = paste0(".", tools::file_ext(fixUploadedFilesNames(input$file)$name)),
                        replacement = "", x = fixUploadedFilesNames(input$file)$name, ignore.case = TRUE),
                   "_IVI_based_network", ".pdf")
        },
        content = function(file) {
            ggsave(plot = finalplot(), filename = file, device = "pdf",
                   width = input$figure.width, height = input$figure.height, units = "in")
        }
    )

    output$download_network_PNG <- downloadHandler(
        filename = function() {
            paste0(gsub(pattern = paste0(".", tools::file_ext(fixUploadedFilesNames(input$file)$name)),
                        replacement = "", x = fixUploadedFilesNames(input$file)$name, ignore.case = TRUE),
                   "_IVI_based_network", ".png")
        },
        content = function(file) {
            ggsave(plot = finalplot(), filename = file, device = "png",
                   width = input$figure.width, height = input$figure.height, units = "in", dpi = input$PNG.resolution)
        }
    )

    # disable or hide the buttons and inputs on page load
    shinyjs::disable("download_network_PDF")
    shinyjs::disable("download_network_PNG")
    shinyjs::disable("jumpToVisTab")
    shinyjs::disable("figure.width")
    shinyjs::disable("figure.height")
    shinyjs::disable("PNG.resolution")
    shinyjs::hide("arrow.width")
    shinyjs::hide("arrow.length")

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



