# Assessment of innate features and associations of two network centrality measures (dependent and independent)

This function assesses innate features and the association of two
centrality measures (or any two other continuous variables) from the
aspect of distribution mode, dependence, linearity, monotonicity,
partial-moments based correlation, and conditional probability of
deviating from corresponding means in opposite direction. This function
assumes one variable as dependent and the other as independent for
regression analyses. The non-linear nature of the association of two
centrality measures is evaluated based on generalized additive models
(GAM). The monotonicity of the association is evaluated based on
comparing the squared coefficient of Spearman correlation and R-squared
of rank regression analysis. Also, the correlation between two variables
is assessed via non-linear non-parametric statistics (NNS). For the
conditional probability assessment, the independent variable is
considered as the condition variable.

## Usage

``` r
double.cent.assess(
  data,
  nodes.colname,
  dependent.colname,
  independent.colname,
  plot = FALSE
)
```

## Arguments

- data:

  A data frame containing the values of two continuous variables and the
  name of observations (nodes).

- nodes.colname:

  The character format (quoted) name of the column containing the name
  of observations (nodes).

- dependent.colname:

  The character format (quoted) name of the column containing the values
  of the dependent variable.

- independent.colname:

  The character format (quoted) name of the column containing the values
  of the independent variable.

- plot:

  logical; FALSE (default) Plots quadrant means of NNS correlation
  analysis.

## Value

A list of 11 objects including:

\- Summary of the basic statistics of two centrality measures (or any
two other continuous variables).

\- The results of normality assessment of two variable (p-value \> 0.05
imply that the variable is normally distributed).

\- Description of the normality assessment of the dependent variable.

\- Description of the normality assessment of the independent variable.

\- Results of the generalized additive modeling (GAM) of the data.

\- The association type based on simultaneous consideration of normality
assessment, GAM Computation with smoothness estimation, Spearman
correlation, and ranked regression analysis of splines.

\- The Hoeffding's D Statistic of dependence (ranging from -0.5 to 1).

\- Description of the dependence significance.

\- Correlation between variables based on the NNS method.

\- The last two objects are the conditional probability of deviation of
two centrality measures from their corresponding means in opposite
directions based on both the entire network and the split-half random
sample of network nodes.

## See also

[`ad.test`](https://rdrr.io/pkg/nortest/man/ad.test.html) for
Anderson-Darling test for normality,
[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) for Generalized additive
models with integrated smoothness estimation,
[`lm`](https://rdrr.io/r/stats/lm.html) for Fitting Linear Models,
[`hoeffd`](https://rdrr.io/pkg/Hmisc/man/hoeffd.html) for Matrix of
Hoeffding's D Statistics, and
[`NNS.dep`](https://rdrr.io/pkg/NNS/man/NNS.dep.html) for NNS Dependence

Other centrality association assessment functions:
[`cond.prob.analysis()`](https://asalavaty.github.io/influential/reference/cond.prob.analysis.md),
[`double.cent.assess.noRegression()`](https://asalavaty.github.io/influential/reference/double.cent.assess.noRegression.md)

## Examples

``` r
if (FALSE) { # \dontrun{
MyData <- centrality.measures
My.metrics.assessment <- double.cent.assess(data = MyData,
                                            nodes.colname = rownames(MyData),
                                            dependent.colname = "BC",
                                            independent.colname = "NC")
} # }
```
