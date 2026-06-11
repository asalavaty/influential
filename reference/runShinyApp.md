# Run shiny app

Run shiny apps included in the influential R package. Also, a web-based
[Influential Software Package](https://influential.erc.monash.edu/) with
a convenient user-interface (UI) has been developed for the comfort of
all users including those without a coding background.

## Usage

``` r
runShinyApp(shinyApp)
```

## Arguments

- shinyApp:

  The name of the shiny app you want to run. You can get the exact name
  of the available shiny apps via the following command.
  *list.files(system.file("ShinyApps", package = "influential"))*.
  Please also note this function is case-sensitive.

## Value

A shiny app.

## Examples

``` r
if (FALSE) { # \dontrun{
runShinyApp(shinyApp = "IVI")
} # }
```
