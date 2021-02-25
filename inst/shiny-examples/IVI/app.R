# Required packages
library(shiny)
library(shinythemes)
library(shinyWidgets)
library(shinyjs)
library(shinycssloaders)
library(colourpicker)
library(ggplot2)
library(igraph)

options(shiny.maxRequestSize = 100 * 1024^2)
options(warn=-1)

####**********************************************####

navbarPageWithText <- function(..., text) {
    navbar <- navbarPage(...)
    textEl <- tags$p(class = "navbar-text", text)
    navbar[[3]][[1]]$children[[1]] <- htmltools::tagAppendChild(
        navbar[[3]][[1]]$children[[1]], textEl)
    navbar
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

####**********************************************####

# Neighborhood Connectivity function

neighborhood.connectivity <- function(graph, vertices = V(graph), mode = "all") {

    # Getting the names of vertices
    if(class(vertices) == "igraph.vs") {
        node.names <- as.character(igraph::as_ids(vertices))
    } else {
        node.names <- as.character(vertices)
    }


    # Getting the first neighbors of each node
    node.neighbors <- sapply(as.list(node.names),
                             FUN = function(i) as.character(igraph::as_ids(igraph::neighbors(graph = graph,
                                                                                             v = i,
                                                                                             mode = mode))))

    # Getting the neighborhood size of each node
    first.neighbors.size <- sapply(node.neighbors,
                                   function(s) igraph::neighborhood.size(graph = graph,
                                                                         nodes = s,
                                                                         mode = mode,
                                                                         order = 1) - 1)

    # Calculation of neighborhood connectivity

    if(length(vertices) == 1) {
        first.neighbors.size.sum <- sum(first.neighbors.size)
        temp.nc <- first.neighbors.size.sum/nrow(node.neighbors)
    } else {

        first.neighbors.size.sum <- sapply(first.neighbors.size, sum)
        temp.nc <- vector(mode = "numeric", length = length(vertices))

        for (i in 1:length(vertices)) {
            temp.nc[i] <- first.neighbors.size.sum[i]/length(node.neighbors[[i]])
        }
    }

    temp.nc[c(which(is.nan(temp.nc)), which(is.na(temp.nc)))] <- 0

    names(temp.nc) <- node.names

    return(temp.nc)

}

####**********************************************####

# H-index function

h_index <- function(graph, vertices = V(graph), mode = "all") {

    # Getting the first neighbors of each node
    first.neighbors <- igraph::neighborhood(graph, nodes = vertices, mode = mode)

    # Getting the neighbors of each vertex
    node.neighbors <- sapply(first.neighbors, function(n) rownames(as.matrix(n[[]][-1])), simplify = F)

    # Getting the neighborhood size of each node
    first.neighbors.size <- lapply(node.neighbors, function(s) igraph::neighborhood.size(graph, s,
                                                                                         mode = mode, order = 1) - 1)

    # Calculation of H-index
    hindex <- vector(mode = "integer", length = length(vertices))

    for (i in 1:length(vertices)) {

        temp.neighbors.size <- unlist(first.neighbors.size[i])

        temp.neighbors.size <- temp.neighbors.size[order(temp.neighbors.size,
                                                         decreasing = TRUE)]
        if(length(temp.neighbors.size) == 0) { hindex[i] <- 0

        } else if(max(temp.neighbors.size) == 0) {
            hindex[i] <- 0
        } else {hindex[i] <- utils::tail(which(temp.neighbors.size >=
                                                   seq_along(temp.neighbors.size)), 1)
        }

        rm(temp.neighbors.size)
    }

    # Getting the names of vertices
    if(class(vertices) == "igraph.vs") {
        node.names <- as.character(igraph::as_ids(vertices))
    } else {
        node.names <- as.character(vertices)
    }
    names(hindex) <- node.names
    return(hindex)
}

####**********************************************####

# Local H-index function

lh_index <- function(graph, vertices = V(graph), mode = "all") {

    # Getting the first neighbors of each node
    first.neighbors <- igraph::neighborhood(graph, nodes = vertices, mode = mode)

    # Calculation of local H-index (LH-index)
    lhindex <- vector(mode = "integer", length = length(vertices))

    for (i in 1:length(vertices)) {

        lhindex[i] <- sum(h_index(graph = graph,
                                  vertices = unlist(first.neighbors[i]),
                                  mode = mode))
    }
    # Getting the names of vertices
    if(class(vertices) == "igraph.vs") {
        node.names <- as.character(igraph::as_ids(vertices))
    } else {
        node.names <- as.character(vertices)
    }
    names(lhindex) <- node.names
    return(lhindex)
}

####**********************************************####

# Collective Influence (CI) function

collective.influence <- function(graph, vertices = V(graph), mode = "all", d=3) {

    ci <- vector(mode="numeric", length=length(vertices))  # collective influence output

    reduced.degrees <- degree(graph = graph,
                              v = vertices,
                              mode = mode) - 1

    # Only identify nodes at distance d
    nodes.at.distance <- igraph::neighborhood(graph = graph, nodes = vertices,
                                              mode = mode, order=d, mindist=d)

    for (i in 1:length(nodes.at.distance)) {
        rd <- reduced.degrees[i]  # i is the index of the node
        rd.list <- reduced.degrees[igraph::as_ids(nodes.at.distance[[i]])]
        # Setting 0 as default in case the graph doesn't have 2nd-order neighbours
        rd.neighbours <- ifelse(length(rd.list) > 0, sum(rd.list), 0)
        ci[i] <- rd * rd.neighbours
    }

    # Getting the names of vertices
    if(class(vertices) == "igraph.vs") {
        node.names <- as.character(igraph::as_ids(vertices))
    } else {
        node.names <- as.character(vertices)
    }
    names(ci) <- node.names
    return(ci)
}

####**********************************************####

# ClusterRank function

clusterRank <- function(graph, vids = V(graph),
                        directed = FALSE, loops = TRUE) {

    vertex.transitivity <- vector(mode = "numeric")

    if(directed) {

        cl.Rank.mode <- "out"

        for(i in V(graph)) {
            vertex.neighborhood <- igraph::neighborhood(graph = graph,
                                                        order = 1, nodes=i,
                                                        mode = cl.Rank.mode)[[1]][-1]
            if(length(vertex.neighborhood) < 2){
                vertex.transitivity <- base::append(vertex.transitivity, NaN)
            } else {
                indc.subgraph <- igraph::induced.subgraph(graph = graph, vertex.neighborhood)
                vertex.transitivity <- base::append(vertex.transitivity,
                                                    igraph::ecount(indc.subgraph)/(igraph::vcount(indc.subgraph)*(igraph::vcount(indc.subgraph)-1)))
            }
        }

    } else {

        cl.Rank.mode <- "all"

        vertex.transitivity <- igraph::transitivity(graph = graph, type = "local")

    }

    if(class(vids) == "igraph.vs") {
        vertices.index <- stats::na.omit(match(vids, V(graph)))
    } else {
        vertices.index <- stats::na.omit(match(vids, igraph::as_ids(V(graph))))
    }

    cl.Rank <- vector(mode = "numeric")

    for (i in V(graph)[vertices.index]) {
        if (is.nan(vertex.transitivity[i])) {
            cl.Rank <- append(cl.Rank, NaN)
        }
        else {

            selected.v.neighborhood <- igraph::neighborhood(graph = graph,
                                                            order = 1, nodes = i,
                                                            mode = cl.Rank.mode)[[1]][-1]
            temp.cl.Rank <- 0
            for (j in selected.v.neighborhood) {
                temp.cl.Rank <- temp.cl.Rank + igraph::degree(graph = graph,
                                                              v = j, mode = cl.Rank.mode,
                                                              loops = loops) + 1
            }
            cl.Rank <- append(cl.Rank, temp.cl.Rank * vertex.transitivity[i])
        }
    }

    if (igraph::is.named(graph)) {
        names(cl.Rank) <- igraph::V(graph)$name[vertices.index]
    }

    return(cl.Rank)
}

####**********************************************####

# IVI function
ivi <- function(graph, vertices = V(graph), weights = NULL, directed = FALSE,
                mode = "all", loops = TRUE, d = 3, scaled = TRUE) {

    #Calculation of required centrality measures

    DC <- igraph::degree(graph = graph, v = vertices, mode = mode, loops = loops)
    CR <- clusterRank(graph = graph, vids = vertices, directed = directed, loops = loops)
    LH_index <- lh_index(graph = graph, vertices = vertices, mode = mode)
    NC <- neighborhood.connectivity(graph = graph, vertices = vertices, mode = mode)
    BC <- igraph::betweenness(graph = graph, v = vertices, directed = directed, weights = weights)
    CI <- collective.influence(graph = graph, vertices = vertices, mode = mode, d = d)

    #Generating temporary measures

    temp.DC <- DC
    temp.CR <- CR
    temp.LH_index <- LH_index
    temp.NC <- unlist(NC)
    temp.BC <- BC
    temp.CI <- CI

    #Removing the NAN and NA values

    temp.DC[c(which(is.nan(temp.DC)), which(is.na(temp.DC)))] <- 0
    temp.CR[c(which(is.nan(temp.CR)), which(is.na(temp.CR)))] <- 0
    temp.LH_index[c(which(is.nan(temp.LH_index)), which(is.na(temp.LH_index)))] <- 0
    temp.NC[c(which(is.nan(temp.NC)), which(is.na(temp.NC)))] <- 0
    temp.BC[c(which(is.nan(temp.BC)), which(is.na(temp.BC)))] <- 0
    temp.CI[c(which(is.nan(temp.CI)), which(is.na(temp.CI)))] <- 0

    #1-100 normalization of centrality measures

    if(any(temp.DC > 0)) {
        temp.DC <- 1+(((temp.DC-min(temp.DC))*(100-1))/(max(temp.DC)-min(temp.DC)))
    }

    if(any(temp.CR > 0)) {
        temp.CR <- 1+(((temp.CR-min(temp.CR))*(100-1))/(max(temp.CR)-min(temp.CR)))
    }

    if(any(temp.LH_index > 0)) {
        temp.LH_index <- 1+(((temp.LH_index-min(temp.LH_index))*(100-1))/(max(temp.LH_index)-min(temp.LH_index)))
    }

    if(any(temp.NC > 0)) {
        temp.NC <- 1+(((temp.NC-min(temp.NC))*(100-1))/(max(temp.NC)-min(temp.NC)))
    }

    if(any(temp.BC > 0)) {
        temp.BC <- 1+(((temp.BC-min(temp.BC))*(100-1))/(max(temp.BC)-min(temp.BC)))
    }

    if(any(temp.CI > 0)) {
        temp.CI <- 1+(((temp.CI-min(temp.CI))*(100-1))/(max(temp.CI)-min(temp.CI)))
    }

    #Calculation of IVI

    spreading.rank <- ((temp.NC+temp.CR)*(temp.BC+temp.CI))

    if(sum(spreading.rank) == 0) {
        spreading.rank[] <- 1
    }

    hubness.rank <- (temp.DC+temp.LH_index)

    if(sum(hubness.rank) == 0) {
        hubness.rank[] <- 1
    }
    temp.ivi <- (hubness.rank)*(spreading.rank)

    #1-100 normalization of IVI

    if(scaled == TRUE) {

        if(length(unique(temp.ivi)) > 1) {
            temp.ivi <- 1+(((temp.ivi-min(temp.ivi))*(100-1))/(max(temp.ivi)-min(temp.ivi)))
        }
    }

    return(temp.ivi)
}

####**********************************************####

# Define the visualization function
cent_network.vis <- function(graph,
                             cent.metric,
                             layout = "kk",
                             node.color = "viridis",
                             node.size.min = 3,
                             node.size.max = 15,
                             dist.power = 1,
                             node.shape = "circle",
                             stroke.size = 1.5,
                             stroke.color = "identical",
                             stroke.alpha = 0.6,
                             show.labels = TRUE,
                             label.cex = 0.4,
                             label.color = "black",
                             directed = FALSE,
                             arrow.width = 25,
                             arrow.length = 0.07,
                             edge.width = 0.5,
                             weighted = FALSE,
                             edge.width.min = 0.2,
                             edge.width.max = 1,
                             edge.color = "grey75",
                             edge.linetype = "solid",
                             legend.position = "right",
                             legend.direction = "vertical",
                             legend.title = "IVI value",
                             boxed.legend = TRUE,
                             show.plot.title = TRUE,
                             plot.title = "IVI-based Network",
                             title.position = "center",
                             show.bottom.border = TRUE,
                             show.left.border = TRUE,
                             seed = 1234) {

    # preparing the layout
    if(layout == "kk") {
        base::set.seed(seed = seed)
        plotcord <- base::data.frame(igraph::layout_with_kk(graph))
    } else if(layout == "star") {
        base::set.seed(seed = seed)
        plotcord <- base::data.frame(igraph::layout_as_star(graph))
    } else if(layout == "tree") {
        base::set.seed(seed = seed)
        plotcord <- base::data.frame(igraph::layout_as_tree(graph))
    } else if(layout == "components") {
        base::set.seed(seed = seed)
        plotcord <- base::data.frame(igraph::layout_components(graph))
    } else if(layout == "circle") {
        base::set.seed(seed = seed)
        plotcord <- base::data.frame(igraph::layout_in_circle(graph))
    } else if(layout == "automatic") {
        base::set.seed(seed = seed)
        plotcord <- base::data.frame(igraph::layout_nicely(graph))
    } else if(layout == "grid") {
        base::set.seed(seed = seed)
        plotcord <- base::data.frame(igraph::layout_on_grid(graph))
    } else if(layout == "dh") {
        base::set.seed(seed = seed)
        plotcord <- base::data.frame(igraph::layout_with_dh(graph))
    } else if(layout == "sphere") {
        base::set.seed(seed = seed)
        plotcord <- base::data.frame(igraph::layout_on_sphere(graph))
    } else if(layout == "random") {
        base::set.seed(seed = seed)
        plotcord <- base::data.frame(igraph::layout_randomly(graph))
    } else if(layout == "drl") {
        base::set.seed(seed = seed)
        plotcord <- base::data.frame(igraph::layout_with_drl(graph))
    } else if(layout == "fr") {
        base::set.seed(seed = seed)
        plotcord <- base::data.frame(igraph::layout_with_fr(graph))
    } else if(layout == "gem") {
        base::set.seed(seed = seed)
        plotcord <- base::data.frame(igraph::layout_with_gem(graph))
    } else if(layout == "graphopt") {
        base::set.seed(seed = seed)
        plotcord <- base::data.frame(igraph::layout_with_graphopt(graph))
    } else if(layout == "lgl") {
        base::set.seed(seed = seed)
        plotcord <- base::data.frame(igraph::layout_with_lgl(graph))
    } else if(layout == "mds") {
        base::set.seed(seed = seed)
        plotcord <- base::data.frame(igraph::layout_with_mds(graph))
    } else if(layout == "sugiyama") {
        base::set.seed(seed = seed)
        plotcord <- base::data.frame(igraph::layout_with_sugiyama(graph))
    }

    ####*******************************####

    # preparing the plotcord table
    base::rownames(plotcord) <- igraph::as_ids(V(graph))
    base::colnames(plotcord) = c("X","Y")

    # add the centrality measure
    plotcord$cent.metric <- cent.metric

    # range normalize the Node size based on the centrality measure
    plotcord$Node.size <- node.size.min+((((cent.metric)^dist.power-min((cent.metric)^dist.power))*(node.size.max-node.size.min))/
                                             (max((cent.metric)^dist.power)-min((cent.metric)^dist.power)))

    # add the Node name
    plotcord$Node.name <- base::as.character(igraph::as_ids(V(graph)))

    ####*******************************####

    # get edges (pairs of node IDs)
    edgelist <- base::data.frame(igraph::get.edgelist(graph))

    # prepare a four column edge data frame with source and destination coordinates
    edges.cord <- base::data.frame(base::matrix(nrow = nrow(edgelist),
                                                ncol = 4), stringsAsFactors = FALSE)
    base::colnames(edges.cord) <- c("X1", "Y1", "X2", "Y2")
    for(i in 1:nrow(edges.cord)) {
        edges.cord$X1[i] <- plotcord[which(igraph::as_ids(V(graph)) %in% edgelist[i,1]),1]
        edges.cord$Y1[i] <- plotcord[which(igraph::as_ids(V(graph)) %in% edgelist[i,1]),2]
        edges.cord$X2[i] <- plotcord[which(igraph::as_ids(V(graph)) %in% edgelist[i,2]),1]
        edges.cord$Y2[i] <- plotcord[which(igraph::as_ids(V(graph)) %in% edgelist[i,2]),2]
    }

    # refine end positions of edges for directerd networks
    if(directed) {
        for (i in 1:nrow(edges.cord)) {

            # correct the x coordinate of arrow
            if(edges.cord$X1[i]>edges.cord$X2[i]) {
                edges.cord$X2[i] <- edges.cord$X2[i]+0.115
            } else if(edges.cord$X1[i]<edges.cord$X2[i]) {
                edges.cord$X2[i] <- edges.cord$X2[i]-0.115
            }

            # correct the y coordinate of arrow
            if(edges.cord$Y1[i]>edges.cord$Y2[i]) {
                edges.cord$Y2[i] <- edges.cord$Y2[i]+0.115
            } else if(edges.cord$Y1[i]<edges.cord$Y2[i]) {
                edges.cord$Y2[i] <- edges.cord$Y2[i]-0.115
            }
        }
    }

    # set the edge width
    if(weighted) {
        edges.cord$Weight <- igraph::E(graph)$weight
        #range normalize the weight
        edges.cord$Weight <- edge.width.min+(((edges.cord$Weight-min(edges.cord$Weight))*(edge.width.max-edge.width.min))/
                                                 (max(edges.cord$Weight)-min(edges.cord$Weight)))
    } else {
        edges.cord$Weight <- edge.width
    }

    ####*******************************####

    # draw the plot
    temp.plot <- ggplot2::ggplot(data = plotcord, ggplot2::aes(x = X, y = Y))

    ##***********##

    # add the edges
    if(directed) {
        temp.plot <- temp.plot +
            ggplot2::geom_segment(data=edges.cord, ggplot2::aes(x=X1, y=Y1, xend = X2, yend = Y2),
                                  size = edges.cord$Weight,
                                  arrow = ggplot2::arrow(angle = arrow.width,
                                                         length = ggplot2::unit(arrow.length, "in"),
                                                         type = "closed"),
                                  colour = edge.color,
                                  linetype = edge.linetype)
    } else {
        temp.plot <- temp.plot +
            ggplot2::geom_segment(data=edges.cord, ggplot2::aes(x=X1, y=Y1, xend = X2, yend = Y2),
                                  size = edges.cord$Weight,
                                  colour = edge.color,
                                  linetype = edge.linetype)
    }

    ##***********##

    # add nodes

    # define node shapes
    node.shape <- base::as.data.frame(node.shape, stringsAsFactors = FALSE)
    for (i in 1:nrow(node.shape)) {
        if(node.shape[i,1] == "circle") {
            node.shape[i,1] <- 21
        } else if(node.shape[i,1] == "square") {
            node.shape[i,1] <- 22
        } else if(node.shape[i,1] == "diamond") {
            node.shape[i,1] <- 23
        } else if(node.shape[i,1] == "triangle") {
            node.shape[i,1] <- 24
        } else if(node.shape[i,1] == "inverted triangle") {
            node.shape[i,1] <- 25
        }
    }

    node.shape <- base::as.numeric(node.shape[,1])

    # add stroke color
    base::suppressWarnings(
        if(stroke.color == "identical") {
            temp.plot <- temp.plot +
                ggplot2::geom_point(data = plotcord, ggplot2::aes(x = X, y = Y, colour = cent.metric),
                                    shape = node.shape,
                                    size = plotcord$Node.size,
                                    stroke = stroke.size,
                                    alpha = stroke.alpha,
                                    show.legend = FALSE) +
                ggplot2::scale_color_viridis_c(option = node.color,
                                               begin = 0.15)
        } else {
            temp.plot <- temp.plot +
                ggplot2::geom_point(data = plotcord, ggplot2::aes(x = X, y = Y),
                                    shape = node.shape,
                                    colour = stroke.color,
                                    size = plotcord$Node.size,
                                    stroke = stroke.size,
                                    alpha = stroke.alpha,
                                    show.legend = FALSE)
        }
    )

    # add node objects
    temp.plot <- temp.plot +
        ggplot2::geom_point(data = plotcord, ggplot2::aes(x = X, y = Y, fill = cent.metric),
                            shape = node.shape,
                            stroke = 0,
                            size = plotcord$Node.size) +

        ##***********##

        # add node color
        ggplot2::scale_fill_viridis_c(option = node.color,
                                      begin = 0.15)

    ##***********##

    # add node labels
    if(show.labels) {
        temp.plot <- temp.plot +
            ggplot2::geom_text(data = plotcord,
                               ggplot2::aes(x = X, y = Y, label=Node.name),
                               size = plotcord$Node.size*label.cex,
                               color = label.color)
    }

    ##***********##
    # expand the x and y limits
    temp.plot <- temp.plot +
        ggplot2::scale_x_continuous(expand=c(0,1)) +
        ggplot2::scale_y_continuous(expand=c(0,1))

    ##***********##

    # add title
    if(show.plot.title) {
        temp.plot <- temp.plot +
            ggplot2::ggtitle(label = plot.title)
    }

    # define title position
    if(title.position == "left") {
        title.position <- 0
    } else if(title.position == "center") {
        title.position <- 0.5
    } else if(title.position == "right") {
        title.position <- 1
    }

    title.position <- base::as.numeric(title.position)

    ##***********##

    # add theme elements

    # add main plot theme
    temp.plot <- temp.plot +
        ggplot2::theme_void() +

        # add legend specifications
        ggplot2::theme(legend.position = legend.position,
                       legend.direction = legend.direction,
                       legend.spacing.y = ggplot2::unit(0.12, "in"),
                       plot.margin = ggplot2::unit(c(0.2,0.2,0.2,0.2),
                                                   units = "cm"),
                       panel.border = ggplot2::element_blank()) +
        ggplot2::labs(fill = legend.title)

    if(boxed.legend) {
        temp.plot <- temp.plot +
            ggplot2::theme(legend.box.background = ggplot2::element_rect(color="black", size=0.5),
                           legend.margin = ggplot2::margin(c(3,3,3,3)),
                           legend.box.margin = ggplot2::margin(c(3,1,3,3)),
                           legend.box.spacing = ggplot2::unit(0, "cm"))
    }

    # add border lines
    if(show.bottom.border) {
        temp.plot <- temp.plot +
            ggplot2::theme(axis.line.x.bottom = ggplot2::element_line(color = 'black'))
    }

    if(show.left.border) {
        temp.plot <- temp.plot +
            ggplot2::theme(axis.line.y.left = ggplot2::element_line(color = 'black'))
    }

    # set title position
    if(show.plot.title) {
        temp.plot <- temp.plot +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = title.position))
    }
    return(temp.plot)
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
                        p("Decoding the information buried within the interconnection of components could have several benefits for the smart control of a complex system. One of the major challenges in this regard is the identification of the most influential individuals that have the potential to cause the highest impact on the entire network. This knowledge could provide the ability to increase network efficiency and reduce costs. In this article, we present a novel algorithm termed the Integrated Value of Influence (IVI) that combines the most important topological characteristics of the network to identify the key individuals within it. The IVI is a versatile method that could benefit several fields such as sociology, economics, transportation, biology, and medicine. In biomedical research, for instance, identification of the true influential nodes within a disease-associated network could lead to the discovery of novel biomarkers and/or drug targets, a process that could have a considerable impact on society."),
                        h3(tags$b("R package influential")),
                        p("The IVI function as well as the centrality-based visualization function are part of the",
                          tags$code("R package", em(tags$strong("influential"))),
                          "Additionally, several other functions have been provided for the calculation of some commonly used centrality measures as well as the extraction,
                          classification and ranking of top candidate features from experimental data. You may install the R package influential via either",
                          a("CRAN", href = "https://cran.r-project.org/package=influential", style = "color:blue;"), "or its",
                          a("GitHub repo", href = "https://github.com/asalavaty/influential", style = "color:blue;"), "."),
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

    # About
    tabPanel("About", icon = icon("address-card"), value = "aboutTab",
             shinyjs::useShinyjs(),
             fluidRow(
                 column(7,
                        h3(tags$b("Credits"), style = "color:darkcyan"),
                        p("The IVI project was done by", a("Adrian (Abbas) Salavaty", href = "https://www.abbassalavaty.com/", style = "color:blue"),
                          "and was supervised by",
                          a("Prof. Peter Currie", href = "https://www.armi.org.au/about/our-people/peter-currie/",
                            style = "color:blue"),
                          "and ",
                          a("Assoc. Prof. Mirana Ramialison", href = "https://www.armi.org.au/about/our-people/mirana-ramialison/",
                          style = "color:blue"), ". The IVI shiny app was designed and developed by Adrian according to the IVI function of the", tags$code("R package influential"),
                          ". Also, the visualization of the networks based on IVI values has been rooted from another function of the",
                          a("influential R package", href = "https://github.com/asalavaty/influential", style = "color:blue"),
                          ". You may have a look at the", em("CITATION"), "tab to get more information regarding the influential R package and how to install it."),
                        br(),
                        p("Also, there is a", a("Youtube channel", href = "https://www.youtube.com/playlist?list=PL38ZLo00h-YHu2SbnQ-lfh4iaIsMQ99Qj", style = "color:blue"),
                          "dedicated to tutorial videos of different functions of the influential R package.")
                        ),
                 column(5,
                        h3(tags$b("Acknowledgments"), style = "color:darkcyan"),
                        tags$ul(
                        tags$li(p("We would like to thanks Ehsan Rezaei-Darzi (MSc Biostatistics) for his consultations regarding proper use of statistical methods and evaluations in the IVI project.
                          Also, thanks to the constructive feedback on the IVI manuscript by Hieu Tri Nim.")),
                        tags$li(p("A part of the results of the IVI project is based on data generated by the TCGA Research Network.")),
                        tags$li(p("The IVI project was supported by Monash University and The Australian Regenerative Medicine Institute (itself supported by grants from the State Government of Victoria and the Australian Government)."))
                        )
                        )
             )
    )

    )

####**********************************************####

server <- function(input, output, session) {

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
                show("weightColumn", anim = TRUE, animType = "slide")
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
                show("pick.stroke.color", anim = TRUE, animType = "fade")
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
            show("edge.width.min.max", anim = TRUE, animType = "slide")
        }
    })

    output$network_figure <- renderPlot({
        req(final.graph())
        Sys.sleep(2)
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
            show(id = "arrow.width", anim = TRUE, animType = "slide")
            show(id = "arrow.length", anim = TRUE, animType = "slide")
        } else {
            hide(id = "arrow.width", anim = TRUE, animType = "slide")
            hide(id = "arrow.length", anim = TRUE, animType = "slide")
        }
    })

    observe({
        if(input$show.plot.title == FALSE) {
            hide(id = "plot.title", anim = TRUE, animType = "slide")
            hide(id = "title.position", anim = TRUE, animType = "slide")
        } else {
            show(id = "plot.title", anim = TRUE, animType = "slide")
            show(id = "title.position", anim = TRUE, animType = "slide")
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
                   width = 8, height = 6, units = "in")
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
                   width = 8, height = 6, units = "in", dpi = 300)
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

    # citation tab

}

# Run the application
shinyApp(ui = ui, server = server)



