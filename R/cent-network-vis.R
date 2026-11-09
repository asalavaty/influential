#=============================================================================
#
#    Visualization of a graph based on centrality measures
#
#=============================================================================

#' Centrality-based network visualization
#'
#' This function has been developed for the visualization of a network based on
#' applying a centrality measure to the size and color of network nodes. You are
#' also able to adjust the directedness and weight of connections. Some of the documentations
#' of the arguments of this function have been adapted from ggplot2 and igraph packages.
#' A shiny app has also been developed for the calculation of IVI as well as IVI-based network
#' visualization, which is accessible using the `influential::runShinyApp("IVI")` command.
#' You can also access the shiny app online at https://influential.erc.monash.edu/.
#' @param graph A graph (network) of the igraph class.
#' @param cent.metric A numeric vector of the desired centrality measure previously
#' calculated by any means. For example, you may use the function \code{\link[influential]{ivi}}
#' for the calculation of the Integrated Value of Influence (IVI) of network nodes. Please note that
#' if the centrality measure has been calculated by any means other than the \code{influential} package, make
#' sure that the order of the values in the \code{cent.metric} vector is consistent with the order of vertices
#' in the network \code{(V(graph))}.
#' @param layout The layout to be used for organizing network nodes. Current available layouts include
#' \code{"kk", "star", "tree", "components", "circle", "automatic", "grid",
#' "sphere", "random", "dh", "drl", "fr", "gem", "graphopt", "lgl", "mds", and "sugiyama"}
#' (default is set to "kk"). For a complete description of different layouts and their
#' underlying algorithms please refer to the function \code{\link[igraph]{layout_}}.
#' @param node.group A vector of the same length as the number of network nodes defining the group each node of the network belongs to.
#' @param node.color A character string indicating the colormap option to use.
#' Five options are available: "magma" (or "A"), "inferno" (or "B"), "plasma"
#' (or "C"), "viridis" (or "D", the default option) and "cividis" (or "E").
#' @param node.size.min The size of nodes with the lowest value of the centrality measure (default is set to 3).
#' @param node.size.max The size of nodes with the highest value of the centrality measure (default is set to 15).
#' @param dist.power The power to be used to visualize more distinction between nodes with high and low
#' centrality measure values. The higher the power, the smaller the nodes with lower values of the centrality
#' measure will become. Default is set to 1, meaning the relative sizes of nodes are reflective of their
#' actual centrality measure values.
#' @param node.shape The shape of nodes. Current available shapes include \code{"circle",
#' "square", "diamond", "triangle", and "inverted triangle"} (default is set to "circle"). You can also
#' set different shapes to different groups of nodes by providing a character vector of shapes of nodes with
#' the same length and order of network vertices. This is useful when plotting a network that include different
#' type of node (for example, up- and down-regulated features).
#' @param stroke.size The size of stroke (border) around the nodes (default is set to 1.5).
#' @param stroke.color The color of stroke (border) around the nodes (default is set to "identical" meaning that the
#' stroke color of a node will be identical to its corresponding node color). You can also
#' set different colors to different groups of nodes by providing a character vector of colors of nodes with
#' the same length and order of network vertices. This is useful when plotting a network that include different
#' type of node (for example, up- and down-regulated features).
#' @param stroke.alpha The transparency of the stroke (border) around the nodes which should
#' be a number between 0 and 1 (default is set to 0.6).
#' @param show.labels Logical scalar, whether to show node labels or not (default is set to TRUE).
#' @param label.cex The amount by which node labels should be scaled relative to the node sizes (default is set to 0.4).
#' @param label.color The color of node labels (default is set to "black").
#' @param directed Logical scalar, whether to draw the network as directed or not (default is set to FALSE).
#' @param arrow.width The width of arrows in the case the network is directed (default is set to 25).
#' @param arrow.length The length of arrows in inch in the case the network is directed (default is set to 0.07).
#' @param edge.width The constant width of edges if the network is unweighted (default is set to 0.5).
#' @param weighted Logical scalar, whether the network is a weighted network or not (default is set to FALSE).
#' @param edge.width.min The width of edges with the lowest weight (default is set to 0.2).
#' This parameter is ignored for unweighted networks.
#' @param edge.width.max The width of edges with the highest weight (default is set to 1).
#' This parameter is ignored for unweighted networks.
#' @param edge.color The color of edges (default is set to "grey75").
#' @param edge.linetype The line type of edges. Current available linetypes include
#' \code{"twodash", "longdash", "dotdash", "dotted", "dashed", and "solid"} (default is set to "solid").
#' @param legend.position The position of legends ("none", "left", "right",
#' "bottom", "top", or two-element numeric vector). The default is set to "right".
#' @param legend.direction layout of items in legends ("horizontal" or "vertical").
#' The default is set to "vertical".
#' @param legend.title The legend title in the string format (default is set to "Centrality measure").
#' @param boxed.legend Logical scalar, whether to draw a box around the legend or not (default is set to TRUE).
#' @param show.plot.title Logical scalar, whether to show the plot title or not (default is set to TRUE).
#' @param plot.title The plot title in the string format (default is set to "Centrality Measure-based Network").
#' @param title.position The position of title ("left", "center", or "right"). The default is set to "center".
#' @param show.bottom.border Logical scalar, whether to draw the bottom border line (default is set to TRUE).
#' @param show.left.border Logical scalar, whether to draw the left border line (default is set to TRUE).
#' @param seed A single value, interpreted as an integer to be used for random number generation for preparing
#' the network layout (default is set to 1234).
#' @return A plot with the class ggplot.
#' @keywords cent_network.vis
#' @family visualization functions
#' @seealso \code{\link[influential]{ivi}}
#' @export cent_network.vis
#' @examples
#' \dontrun{
#' MyData <- coexpression.data
#' My_graph <- graph_from_data_frame(MyData)
#' Graph_IVI <- ivi(graph = My_graph, mode = "all")
#' Graph_IVI_plot <- cent_network.vis(graph = My_graph, cent.metric = Graph_IVI,
#'                                    legend.title = "IVI",
#'                                    plot.title = "IVI-based Network")
#' }
cent_network.vis <- function(graph,
                             cent.metric,
                             layout = "kk",
                             node.group = NULL,
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
                             legend.title = "Centrality\nmeasure",
                             boxed.legend = TRUE,
                             show.plot.title = TRUE,
                             plot.title = "Centrality Measure-based Network",
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
  
  # add node.group
  if(!is.null(node.group)) {
    plotcord$Group <- node.group
  } else {
    plotcord$Group <- plotcord$Node.name
  }
  
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
    if(length(stroke.color) == 1 && stroke.color == "identical") {
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
        ggplot2::geom_point(data = plotcord, ggplot2::aes(x = X, y = Y, color = Group),
                            shape = node.shape,
                            size = plotcord$Node.size,
                            stroke = stroke.size,
                            alpha = stroke.alpha,
                            show.legend = ifelse(length(stroke.color) == 1, FALSE, TRUE)) +
        ggplot2::scale_color_manual(values = stroke.color)
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
