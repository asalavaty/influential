# Visualization of ExIR results

This function has been developed for the visualization of ExIR results.
Some of the documentations of the arguments of this function have been
adapted from ggplot2 package. A shiny app has also been developed for
Running the ExIR model, visualization of its results as well as
computational simulation of knockout and/or up-regulation of its top
candidate outputs, which is accessible using the
\`influential::runShinyApp("ExIR")\` command. You can also access the
shiny app online at https://influential.erc.monash.edu/.

## Usage

``` r
exir.vis(
  exir.results,
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
  show.y.axis.grid = TRUE
)
```

## Arguments

- exir.results:

  An object of class `"ExIR_Result"` which is the output of the function
  `"exir"`.

- synonyms.table:

  (Optional) A data frame or matrix with two columns including a column
  for the used feature names in the input data of the `"exir"` model and
  the other column their synonyms. Note, the original feature names
  should always come as the first column and the synonyms as the second
  one. For example, if the original feature names used for running the
  `"exir"` model are Ensembl gene symbols, you can use their HGNC
  synonyms in the second column to be used for the visualization of the
  ExIR results

- n:

  An integer specifying the number of top candidates to be selected from
  each category of ExIR results (default is set to 10).

- driver.type:

  A string specifying the type of drivers to be used for the selection
  of top N candidates. The possible types include `"combined"` (meaning
  both driver types), `"accelerator"` and `"decelerator"` (default is
  set to "combined").

- biomarker.type:

  A string specifying the type of biomarkers to be used for the
  selection of top N candidates. Possible types include `"combined"`
  (meaning both biomarker types), `"up-regulated"` and
  `"down-regulated"` (default is set to "combined").

- show.drivers:

  Logical scalar, whether to show Drivers or not (default is set to
  TRUE).

- show.biomarkers:

  Logical scalar, whether to show Biomarkers or not (default is set to
  TRUE).

- show.de.mediators:

  Logical scalar, whether to show DE-mediators or not (default is set to
  TRUE).

- show.nonDE.mediators:

  Logical scalar, whether to show nonDE-mediators or not (default is set
  to TRUE).

- basis:

  A string specifying the basis for the selection of top N candidates
  from each category of the results. Possible options include `"Rank"`
  and `"Adjusted p-value"` (default is set to "Rank").

- label.position:

  By default, the labels are displayed on the top of the plot. Using
  label.position it is possible to place the labels on either of the
  four sides by setting label.position = c("top", "bottom", "left",
  "right").

- nrow:

  Number of rows of the plot (default is set to 1).

- dot.size.min:

  The size of dots with the lowest statistical significance (default is
  set to 2).

- dot.size.max:

  The size of dots with the highest statistical significance (default is
  set to 5).

- type.color:

  A character string or function indicating the color palette to be used
  for the visualization of different types of candidates. You may choose
  one of the Viridis palettes including "magma" (or "A"), "inferno" (or
  "B"), "plasma" (or "C"), "viridis" (or "D", the default option) and
  "cividis" (or "E"), use a function specifying your desired palette, or
  manually specify the vector of colors for different types.

- stroke.size:

  The size of stroke (border) around the dots (default is set to 1.5).

- stroke.alpha:

  The transparency of the stroke (border) around the dots which should
  be a number between 0 and 1 (default is set to 1).

- dot.color.low:

  The color to be used for the visualization of dots (features) with the
  lowest Z-score values (default is set to "blue").

- dot.color.high:

  The color to be used for the visualization of dots (features) with the
  highest Z-score values (default is set to "red").

- legend.position:

  The position of legends ("none", "left", "right", "bottom", "top", or
  two-element numeric vector). The default is set to "bottom".

- legend.direction:

  Layout of items in legends ("horizontal" or "vertical"). The default
  is set to "vertical".

- legends.layout:

  Layout of different legends of the plot ("horizontal" or "vertical").
  The default is set to "horizontal".

- boxed.legend:

  Logical scalar, whether to draw a box around the legend or not
  (default is set to TRUE).

- show.plot.title:

  Logical scalar, whether to show the plot title or not (default is set
  to TRUE).

- plot.title:

  The plot title in the string format (default is set to "auto" which
  automatically generates a title for the plot).

- title.position:

  The position of title ("left", "center", or "right"). The default is
  set to "left".

- plot.title.size:

  The font size of the plot title (default is set to 12).

- show.plot.subtitle:

  Logical scalar, whether to show the plot subtitle or not (default is
  set to TRUE).

- plot.subtitle:

  The plot subtitle in the string format (default is set to "auto" which
  automatically generates a subtitle for the plot).

- subtitle.position:

  The position of subtitle ("left", "center", or "right"). The default
  is set to "left".

- y.axis.title:

  The title of the y axis (features title). Default is set to
  "Features".

- show.y.axis.grid:

  Logical scalar, whether to draw y axis grid lines (default is set to
  TRUE).

## Value

A plot with the class ggplot.

## See also

[`exir`](https://asalavaty.github.io/influential/reference/exir.md)

Other visualization functions:
[`cent_network.vis()`](https://asalavaty.github.io/influential/reference/cent_network.vis.md)

## Examples

``` r
if (FALSE) { # \dontrun{
MyResults <- exir.results
ExIR.plot <- exir.vis(exir.results = MyResults, n = 5)
} # }
```
