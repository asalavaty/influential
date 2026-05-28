#=============================================================================
#
#    Visualization of ExIR results
#
#=============================================================================

#' Visualization of ExIR results
#'
#' This function has been developed for the visualization of ExIR results. Some of the documentations
#' of the arguments of this function have been adapted from ggplot2 package.
#' A shiny app has also been developed for Running the ExIR model, visualization of its results as well as computational
#' simulation of knockout and/or up-regulation of its top candidate outputs, which is accessible using
#' the `influential::runShinyApp("ExIR")` command.
#' You can also access the shiny app online at https://influential.erc.monash.edu/.
#' @param exir.results An object of class \code{"ExIR_Result"} which is the output of the function \code{"exir"}.
#' @param synonyms.table (Optional) A data frame or matrix with two columns including a column for the used feature
#' names in the input data of the \code{"exir"} model and the other column their synonyms. Note, the original feature names should
#' always come as the first column and the synonyms as the second one. For example, if
#' the original feature names used for running the \code{"exir"} model are Ensembl gene
#' symbols, you can use their HGNC synonyms in the second column to be used for the visualization of the ExIR results
#' @param n An integer specifying the number of top candidates to be selected from each category of ExIR results (default is set to 10).
#' @param driver.type A string specifying the type of drivers to be used for the selection of top N candidates. The possible types
#' include \code{"combined"} (meaning both driver types), \code{"accelerator"} and \code{"decelerator"} (default is set to "combined").
#' @param biomarker.type A string specifying the type of biomarkers to be used for the selection of top N candidates. Possible types
#' include \code{"combined"} (meaning both biomarker types), \code{"up-regulated"} and \code{"down-regulated"} (default is set to "combined").
#' @param show.drivers Logical scalar, whether to show Drivers or not (default is set to TRUE).
#' @param show.biomarkers Logical scalar, whether to show Biomarkers or not (default is set to TRUE).
#' @param show.de.mediators Logical scalar, whether to show DE-mediators or not (default is set to TRUE).
#' @param show.nonDE.mediators Logical scalar, whether to show nonDE-mediators or not (default is set to TRUE).
#' @param basis A string specifying the basis for the selection of top N candidates from each category of the results. Possible options include
#' \code{"Rank"} and \code{"Adjusted p-value"} (default is set to "Rank").
#' @param label.position By default, the labels are displayed on the top of the plot. Using label.position it is possible
#' to place the labels on either of the four sides by setting label.position = c("top", "bottom", "left", "right").
#' @param nrow Number of rows of the plot (default is set to 1).
#' @param dot.size.min The size of dots with the lowest statistical significance (default is set to 2).
#' @param dot.size.max The size of dots with the highest statistical significance (default is set to 5).
#' @param type.color A character string or function indicating the color palette to be used for the visualization of
#' different types of candidates. You may choose one of the Viridis palettes including "magma" (or "A"),
#' "inferno" (or "B"), "plasma" (or "C"), "viridis" (or "D", the default option) and "cividis" (or "E"), use a function specifying
#' your desired palette, or manually specify the vector of colors for different types.
#' @param stroke.size The size of stroke (border) around the dots (default is set to 1.5).
#' @param stroke.alpha The transparency of the stroke (border) around the dots which should
#' be a number between 0 and 1 (default is set to 1).
#' @param dot.color.low The color to be used for the visualization of dots (features) with the lowest Z-score values (default is set to "blue").
#' @param dot.color.high The color to be used for the visualization of dots (features) with the highest Z-score values (default is set to "red").
#' @param legend.position The position of legends ("none", "left", "right",
#' "bottom", "top", or two-element numeric vector). The default is set to "bottom".
#' @param legend.direction Layout of items in legends ("horizontal" or "vertical").
#' The default is set to "vertical".
#' @param legends.layout Layout of different legends of the plot ("horizontal" or "vertical").
#' The default is set to "horizontal".
#' @param boxed.legend Logical scalar, whether to draw a box around the legend or not (default is set to TRUE).
#' @param show.plot.title Logical scalar, whether to show the plot title or not (default is set to TRUE).
#' @param plot.title The plot title in the string format (default is set to "auto" which automatically generates a title for the plot).
#' @param title.position The position of title ("left", "center", or "right"). The default is set to "left".
#' @param plot.title.size The font size of the plot title (default is set to 12).
#' @param show.plot.subtitle Logical scalar, whether to show the plot subtitle or not (default is set to TRUE).
#' @param plot.subtitle The plot subtitle in the string format (default is set to "auto" which automatically generates a subtitle for the plot).
#' @param subtitle.position The position of subtitle ("left", "center", or "right"). The default is set to "left".
#' @param y.axis.title The title of the y axis (features title). Default is set to "Features".
#' @param show.y.axis.grid Logical scalar, whether to draw y axis grid lines (default is set to TRUE).
#' @return A plot with the class ggplot.
#' @keywords exir.vis
#' @family visualization functions
#' @seealso \code{\link[influential]{exir}}
#' @export exir.vis
#' @examples
#' \dontrun{
#' MyResults <- exir.results
#' ExIR.plot <- exir.vis(exir.results = MyResults, n = 5)
#' }
exir.vis <- function(exir.results,
                     synonyms.table = NULL,
                     n = 10,
                     driver.type = "combined",
                     biomarker.type = "combined",
                     show.drivers = TRUE,
                     show.biomarkers = TRUE,
                     show.de.mediators = TRUE,
                     show.nonDE.mediators = TRUE,
                     basis = "Rank",
                     label.position = "top",
                     nrow = 1,
                     dot.size.min = 2,
                     dot.size.max = 5,
                     type.color = "viridis",
                     stroke.size = 1.5,
                     stroke.alpha = 1,
                     dot.color.low = "blue",
                     dot.color.high = "red",
                     legend.position = "bottom",
                     legend.direction = "vertical",
                     legends.layout = "horizontal",
                     boxed.legend = TRUE,
                     show.plot.title = TRUE,
                     plot.title = "auto",
                     title.position = "left",
                     plot.title.size = 12,
                     show.plot.subtitle = TRUE,
                     plot.subtitle = "auto",
                     subtitle.position = "left",
                     y.axis.title = "Feature",
                     show.y.axis.grid = TRUE) {
  
  if(base::identical(base::inherits(exir.results, "ExIR_Result"), FALSE)) {
    stop("The provided ExIR model-result is wrong. The exir.results should be
           the output of the exir model and have a 'ExIR_Result' class.",
         call. = FALSE)
  }
  
  # prepare a list for storing the results
  exir.for.plot <- base::list()
  
  # select top N features
  
  # for driver table
  if(any(base::names(exir.results) == "Driver table")) {
    top.N.driver.table <- exir.results[[which(names(exir.results) %in% "Driver table")]]
    top.N.driver.table$Feature <- base::rownames(top.N.driver.table)
    top.N.driver.table$Class <- "Driver"
    
    # top N combined
    if(driver.type == "combined" & nrow(top.N.driver.table)> n) {
      
      if(basis == "Rank") {
        top.drivers.index <- utils::head(order(top.N.driver.table$Rank),
                                         n = n)
        
      } else if(basis == "Adjusted p-value") {
        top.drivers.index <- utils::head(order(top.N.driver.table$P.adj),
                                         n = n)
      }
      top.N.driver.table <- top.N.driver.table[top.drivers.index,]
      top.N.driver.table$Rank <- base::rank(top.N.driver.table$Rank,
                                            ties.method = "min")
      
      # top N accelerator
    } else if(driver.type == "accelerator") {
      top.N.driver.table <- base::subset(top.N.driver.table,
                                         Type == "Accelerator")
      
      if(basis == "Rank") {
        top.drivers.index <- utils::head(order(top.N.driver.table$Rank),
                                         n = n)
        
      } else if(basis == "Adjusted p-value") {
        top.drivers.index <- utils::head(order(top.N.driver.table$P.adj),
                                         n = n)
      }
      top.N.driver.table <- top.N.driver.table[top.drivers.index,]
      top.N.driver.table$Rank <- base::rank(top.N.driver.table$Rank,
                                            ties.method = "min")
      
      # top N decelerator
    } else if(driver.type == "decelerator") {
      top.N.driver.table <- base::subset(top.N.driver.table,
                                         Type == "Decelerator")
      if(basis == "Rank") {
        top.drivers.index <- utils::head(order(top.N.driver.table$Rank),
                                         n = n)
        
      } else if(basis == "Adjusted p-value") {
        top.drivers.index <- utils::head(order(top.N.driver.table$P.adj),
                                         n = n)
      }
      top.N.driver.table <- top.N.driver.table[top.drivers.index,]
      top.N.driver.table$Rank <- base::rank(top.N.driver.table$Rank,
                                            ties.method = "min")
    }
    base::rownames(top.N.driver.table) <- NULL
    top.N.driver.table$Type[top.N.driver.table$Type == "Accelerator"] <- "Accelerator\ndriver"
    top.N.driver.table$Type[top.N.driver.table$Type == "Decelerator"] <- "Decelerator\ndriver"
    exir.for.plot <- base::append(x = exir.for.plot,
                                  values = base::list(top.N.driver.table))
  }
  
  ####********************####
  
  # for biomarker table
  if(any(base::names(exir.results) == "Biomarker table")) {
    top.N.biomarker.table <- exir.results[[which(names(exir.results) %in% "Biomarker table")]]
    top.N.biomarker.table$Feature <- base::rownames(top.N.biomarker.table)
    top.N.biomarker.table$Class <- "Biomarker"
    
    # top N combined
    if(biomarker.type == "combined" & nrow(top.N.biomarker.table)> n) {
      
      if(basis == "Rank") {
        top.biomarkers.index <- utils::head(order(top.N.biomarker.table$Rank),
                                            n = n)
        
      } else if(basis == "Adjusted p-value") {
        top.biomarkers.index <- utils::head(order(top.N.biomarker.table$P.adj),
                                            n = n)
      }
      top.N.biomarker.table <- top.N.biomarker.table[top.biomarkers.index,]
      top.N.biomarker.table$Rank <- base::rank(top.N.biomarker.table$Rank,
                                               ties.method = "min")
      
      # top N up-regulated
    } else if(biomarker.type == "up-regulated") {
      top.N.biomarker.table <- base::subset(top.N.biomarker.table,
                                            Type == "Up-regulated")
      if(basis == "Rank") {
        top.biomarkers.index <- utils::head(order(top.N.biomarker.table$Rank),
                                            n = n)
        
      } else if(basis == "Adjusted p-value") {
        top.biomarkers.index <- utils::head(order(top.N.biomarker.table$P.adj),
                                            n = n)
      }
      top.N.biomarker.table <- top.N.biomarker.table[top.biomarkers.index,]
      top.N.biomarker.table$Rank <- base::rank(top.N.biomarker.table$Rank,
                                               ties.method = "min")
      
    } else if(biomarker.type == "down-regulated") {
      top.N.biomarker.table <- base::subset(top.N.biomarker.table,
                                            Type == "Down-regulated")
      if(basis == "Rank") {
        top.biomarkers.index <- utils::head(order(top.N.biomarker.table$Rank),
                                            n = n)
        
      } else if(basis == "Adjusted p-value") {
        top.biomarkers.index <- utils::head(order(top.N.biomarker.table$P.adj),
                                            n = n)
      }
      top.N.biomarker.table <- top.N.biomarker.table[top.biomarkers.index,]
      top.N.biomarker.table$Rank <- base::rank(top.N.biomarker.table$Rank,
                                               ties.method = "min")
      
    }
    base::rownames(top.N.biomarker.table) <- NULL
    top.N.biomarker.table$Type[top.N.biomarker.table$Type == "Up-regulated"] <- "Up-regulated\nbiomarker"
    top.N.biomarker.table$Type[top.N.biomarker.table$Type == "Down-regulated"] <- "Down-regulated\nbiomarker"
    
    exir.for.plot <- base::append(x = exir.for.plot,
                                  values = base::list(top.N.biomarker.table))
  }
  
  ####********************####
  
  # for nonDE-mediator table
  if(any(base::names(exir.results) == "nonDE-mediator table")) {
    top.N.nonDE.mediator.table <- exir.results[[which(names(exir.results) %in% "nonDE-mediator table")]]
    top.N.nonDE.mediator.table$Type <- "nonDE-mediator"
    top.N.nonDE.mediator.table$Feature <- base::rownames(top.N.nonDE.mediator.table)
    top.N.nonDE.mediator.table$Class <- "nonDE-mediator"
    
    # top N combined
    
    if(basis == "Rank") {
      top.nonDE.mediators.index <- utils::head(order(top.N.nonDE.mediator.table$Rank),
                                               n = n)
      
    } else if(basis == "Adjusted p-value") {
      top.nonDE.mediators.index <- utils::head(order(top.N.nonDE.mediator.table$P.adj),
                                               n = n)
    }
    top.N.nonDE.mediator.table <- top.N.nonDE.mediator.table[top.nonDE.mediators.index,]
    top.N.nonDE.mediator.table$Rank <- base::rank(top.N.nonDE.mediator.table$Rank,
                                                  ties.method = "min")
    base::rownames(top.N.nonDE.mediator.table) <- NULL
    exir.for.plot <- base::append(x = exir.for.plot,
                                  values = base::list(top.N.nonDE.mediator.table))
  }
  
  ####********************####
  
  # for DE-mediator table
  if(any(base::names(exir.results) == "DE-mediator table")) {
    top.N.DE.mediator.table <- exir.results[[which(names(exir.results) %in% "DE-mediator table")]]
    top.N.DE.mediator.table$Type <- "DE-mediator"
    top.N.DE.mediator.table$Feature <- base::rownames(top.N.DE.mediator.table)
    top.N.DE.mediator.table$Class <- "DE-mediator"
    
    # top N combined
    
    if(basis == "Rank") {
      top.DE.mediators.index <- utils::head(order(top.N.DE.mediator.table$Rank),
                                            n = n)
      
    } else if(basis == "Adjusted p-value") {
      top.DE.mediators.index <- utils::head(order(top.N.DE.mediator.table$P.adj),
                                            n = n)
    }
    top.N.DE.mediator.table <- top.N.DE.mediator.table[top.DE.mediators.index,]
    top.N.DE.mediator.table$Rank <- base::rank(top.N.DE.mediator.table$Rank,
                                               ties.method = "min")
    base::rownames(top.N.DE.mediator.table) <- NULL
    exir.for.plot <- base::append(x = exir.for.plot,
                                  values = base::list(top.N.DE.mediator.table))
  }
  
  # combine the results for plotting
  exir.for.plot <- base::Reduce(function(dtf1, dtf2) base::merge(dtf1, dtf2,
                                                                 all = TRUE,
                                                                 all.x = TRUE),
                                exir.for.plot)
  
  # correct the features names
  if(!is.null(synonyms.table)) {
    synonyms.table <- base::as.data.frame(synonyms.table, stringsAsFactors = FALSE)
    synonyms.index <- base::match(exir.for.plot$Feature,
                                  synonyms.table[,1])
    exir.for.plot$Feature <- synonyms.table[synonyms.index,2]
  }
  
  # correct the Type levels
  driver.levels <- base::unique(exir.for.plot$Type[base::grep("driver", exir.for.plot$Type)])
  biomarker.levels <- base::unique(exir.for.plot$Type[base::grep("biomarker", exir.for.plot$Type)])
  mediators.levels <- base::unique(exir.for.plot$Type[base::seq_along(rownames(exir.for.plot))[-c(base::grep("driver", exir.for.plot$Type),
                                                                                                  base::grep("biomarker", exir.for.plot$Type))]])
  exir.for.plot$Type <- base::factor(exir.for.plot$Type,
                                     levels = c(driver.levels,
                                                biomarker.levels,
                                                mediators.levels))
  
  # correct the Class levels
  mediators.class.levels <- base::unique(exir.for.plot$Class[base::seq_along(rownames(exir.for.plot))[-c(base::grep("Driver", exir.for.plot$Class),
                                                                                                         base::grep("Biomarker", exir.for.plot$Class))]])
  exir.for.plot$Class <- base::factor(exir.for.plot$Class,
                                      levels = c("Driver",
                                                 "Biomarker",
                                                 mediators.class.levels))
  
  # remove undesired classes
  if(isFALSE(show.drivers) & any(exir.for.plot$Class == "Driver")) {
    exir.for.plot <- exir.for.plot[-c(which(exir.for.plot$Class == "Driver")),]
  }
  
  if(isFALSE(show.biomarkers) & any(exir.for.plot$Class == "Biomarker")) {
    exir.for.plot <- exir.for.plot[-c(which(exir.for.plot$Class == "Biomarker")),]
  }
  
  if(isFALSE(show.de.mediators) & any(exir.for.plot$Class == "DE-mediator")) {
    exir.for.plot <- exir.for.plot[-c(which(exir.for.plot$Class == "DE-mediator")),]
  }
  
  if(isFALSE(show.nonDE.mediators) & any(exir.for.plot$Class == "nonDE-mediator")) {
    exir.for.plot <- exir.for.plot[-c(which(exir.for.plot$Class == "nonDE-mediator")),]
  }
  
  # correct the P.adj to be used as the dot size
  
  exir.for.plot$P.value[is.nan(exir.for.plot$P.value)] <- 1
  exir.for.plot$P.adj[is.nan(exir.for.plot$P.adj)] <- 1
  exir.for.plot$P.value[is.na(exir.for.plot$P.value)] <- 1
  exir.for.plot$P.adj[is.na(exir.for.plot$P.adj)] <- 1
  
  if(min(exir.for.plot$P.adj)==0) {
    
    #range normalize the primitive P.adj
    temp.min_P.adj <- base::sort(base::unique(exir.for.plot$P.adj))[2]
    
    exir.for.plot$P.adj <- temp.min_P.adj+
      (((exir.for.plot$P.adj-min(exir.for.plot$P.adj))*(max(exir.for.plot$P.adj)-temp.min_P.adj))/
         (max(exir.for.plot$P.adj)-min(exir.for.plot$P.adj)))
  }
  
  # Set the P.adj based on min and max arguments
  if(length(unique(exir.for.plot$P.adj)) == 1) {
    exir.for.plot$P.adj <- mean(c(dot.size.min, dot.size.max))
    
  } else {
    exir.for.plot$P.adj <- dot.size.min+(((-log10(exir.for.plot$P.adj)-min(-log10(exir.for.plot$P.adj)))*(dot.size.max-dot.size.min))/
                                           (max(-log10(exir.for.plot$P.adj))-min(-log10(exir.for.plot$P.adj))))
  }
  
  # Correct the levels of Features based on each class
  visClassLength <- base::length(base::unique(exir.for.plot$Class))
  
  visFeatureLevels <- base::unique(base::as.character(base::unlist(base::sapply(X = 1:visClassLength,
                                                                                FUN = function(i) {
                                                                                  exir.for.plot$Feature[which(exir.for.plot$Class %in% base::unique(exir.for.plot$Class)[i])]
                                                                                }))))
  
  exir.for.plot$Feature <- base::factor(exir.for.plot$Feature,
                                        levels = visFeatureLevels)
  
  ####*******************************####
  
  # draw the plot
  temp.exir.plot <- ggplot2::ggplot(data = exir.for.plot,
                                    ggplot2::aes(x = Rank, y = Feature)) +
    
    # add node objects
    ggplot2::geom_point(ggplot2::aes(fill = Z.score,
                                     colour = Type,
                                     size = P.adj),
                        shape = 21,
                        stroke = 0,
                        alpha = 1) +
    
    ##***********##
    
    # add color of Type
    ggplot2::geom_point(ggplot2::aes(colour = Type,
                                     size = P.adj),
                        shape = 21,
                        stroke = stroke.size,
                        alpha = stroke.alpha)
  
  
  ##***********##
  
  # add node and stroke colors
  if(base::inherits(type.color, "character") & base::length(type.color) == 1) {
    if(base::length(base::grep(type.color,
                               c("magma", "inferno",
                                 "plasma", "viridis", "cividis",
                                 "A", "B", "C", "D", "E"))) == 1) {
      temp.exir.plot <- temp.exir.plot +
        ggplot2::scale_colour_viridis_d(option = type.color)
    } else {
      temp.exir.plot <- temp.exir.plot +
        ggplot2::scale_colour_manual(values = type.color)
    }
  } else {
    temp.exir.plot <- temp.exir.plot +
      ggplot2::scale_colour_manual(values = type.color)
  }
  
  temp.exir.plot <- temp.exir.plot +
    ggplot2::scale_fill_gradient(name = "Z-score",
                                 low = dot.color.low,
                                 high = dot.color.high) +
    
    ##***********##
    
    # correct size identity inside of aes
    ggplot2::scale_size_identity(guide = "legend") +
    
    
    ##***********##
    
    # add y axis title
    ggplot2::ylab(y.axis.title) +
    
    
    ##***********##
    
    # add facets
    ggplot2::facet_wrap(. ~ Class,
                        scales = "free_x",
                        strip.position = label.position,
                        nrow = nrow) +
    
    
    ##***********##
    
    # add theme elements
    
    ggplot2::theme_bw()
  
  # set the x axis breaks
  rank.uniqueness <- base::vector(mode = "integer")
  for(i in base::unique(exir.for.plot$Class)) {
    rank.uniqueness <- base::append(x = rank.uniqueness, values =
                                      length(unique(exir.for.plot$Rank[exir.for.plot$Class == i]))
    )
  }
  by.x_continuous <- base::ifelse(any(rank.uniqueness == 1), 1, base::round(n/5))
  
  # set the x axis numbers and start
  x.axis.numbers <- function(x) {
    if(by.x_continuous %% 2 == 0 & n > 7) {
      base::seq(2, max(x), by = by.x_continuous)
    } else if(n > 7) {
      base::seq(1, max(x), by = by.x_continuous)
    } else {
      base::seq(1, max(x), by = 1)
    }
  }
  
  # set the x axis limits
  x.axis.limits <- function(x) {
    c((min(x)-(n/50)), (max(x)+(n/50)))
  }
  
  temp.exir.plot <- temp.exir.plot +
    ggplot2::scale_x_continuous(breaks = x.axis.numbers,
                                limits = x.axis.limits) +
    
    ggplot2::theme(legend.title = ggplot2::element_text(size = 10),
                   legend.box = legends.layout,
                   legend.position = legend.position,
                   legend.direction = legend.direction) +
    
    ggplot2::guides(colour = ggplot2::guide_legend(keyheight = 1.5),
                    fill = ggplot2::guide_colorbar(frame.colour = "black",
                                                   barwidth = 1.5),
                    size = ggplot2::guide_legend(title = "Statistical\nsignificance"))
  
  if(boxed.legend) {
    temp.exir.plot <- temp.exir.plot +
      ggplot2::theme(legend.box.background = ggplot2::element_rect(color="black", size=0.5),
                     legend.margin = ggplot2::margin(c(3,3,3,3)),
                     legend.box.margin = ggplot2::margin(c(3,1,3,3)))
  }
  
  ##***********##
  
  # add title
  
  if(plot.title == "auto") {
    plot.title <- "ExIR model-based prioritized features"
  }
  
  if(show.plot.title) {
    temp.exir.plot <- temp.exir.plot +
      ggplot2::labs(title = plot.title)
  }
  
  # add subtitle
  
  if(plot.subtitle == "auto") {
    plot.subtitle <- base::paste(basis,
                                 "-based selection of top ",
                                 n,
                                 " candidates", sep = "")
  }
  
  if(show.plot.subtitle) {
    temp.exir.plot <- temp.exir.plot +
      ggplot2::labs(subtitle = plot.subtitle)
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
  
  # define subtitle position
  if(subtitle.position == "left") {
    subtitle.position <- 0
  } else if(subtitle.position == "center") {
    subtitle.position <- 0.5
  } else if(subtitle.position == "right") {
    subtitle.position <- 1
  }
  
  subtitle.position <- base::as.numeric(subtitle.position)
  
  title.size <- plot.title.size - 2
  
  # set title position
  if(show.plot.title) {
    temp.exir.plot <- temp.exir.plot +
      ggplot2::theme(plot.title = ggplot2::element_text(size = title.size,
                                                        hjust = title.position))
  }
  
  # set subtitle position
  if(show.plot.subtitle) {
    subtitle.size <- title.size - 2
    temp.exir.plot <- temp.exir.plot +
      ggplot2::theme(plot.subtitle = ggplot2::element_text(size = subtitle.size,
                                                           hjust = subtitle.position))
  }
  
  ##***********##
  
  # Set the order of legends
  
  temp.exir.plot <- temp.exir.plot + ggplot2::guides(color = ggplot2::guide_legend(order = 1),
                                                     size = ggplot2::guide_legend(order = 2))
  
  ##***********##
  
  if(base::identical(show.y.axis.grid, FALSE)) {
    temp.exir.plot <- temp.exir.plot +
      ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
  }
  
  return(temp.exir.plot)
}
