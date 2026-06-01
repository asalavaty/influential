# Assessment of innate features and associations of two network centrality measures

This function assesses innate features and the association of two
centrality measures (or any two other continuous variables) from the
aspect of distribution mode, dependence, linearity, partial-moments
based correlation, and conditional probability of deviating from
corresponding means in opposite direction (centrality2 is used as the
condition variable). This function doesn't consider which variable is
dependent and which one is independent and no regression analysis is
done. Also, the correlation between two variables is assessed via
non-linear non-parametric statistics (NNS). For the conditional
probability assessment, the centrality2 variable is considered as the
condition variable.

## Usage

``` r
double.cent.assess.noRegression(
  data,
  nodes.colname,
  centrality1.colname,
  centrality2.colname
)
```

## Arguments

- data:

  A data frame containing the values of two continuous variables and the
  name of observations (nodes).

- nodes.colname:

  The character format (quoted) name of the column containing the name
  of observations (nodes).

- centrality1.colname:

  The character format (quoted) name of the column containing the values
  of the Centrality_1 variable.

- centrality2.colname:

  The character format (quoted) name of the column containing the values
  of the Centrality_2 variable.

## Value

A list of nine objects including:

\- Summary of the basic statistics of two centrality measures (or any
two other continuous variables).

\- The results of normality assessment of two variable (p-value \> 0.05
imply that the variable is normally distributed).

\- Description of the normality assessment of the centrality1 (first
variable).

\- Description of the normality assessment of the centrality2 (second
variable).

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
[`hoeffd`](https://rdrr.io/pkg/Hmisc/man/hoeffd.html) for Matrix of
Hoeffding's D Statistics, and
[`NNS.dep`](https://rdrr.io/pkg/NNS/man/NNS.dep.html) for NNS Dependence

Other centrality association assessment functions:
[`cond.prob.analysis()`](https://asalavaty.github.io/influential/reference/cond.prob.analysis.md),
[`double.cent.assess()`](https://asalavaty.github.io/influential/reference/double.cent.assess.md)

## Examples

``` r
if (FALSE) { # \dontrun{
MyData <- centrality.measures
My.metrics.assessment <- double.cent.assess.noRegression(data = MyData,
                                            nodes.colname = rownames(MyData),
                                            centrality1.colname = "BC",
                                            centrality2.colname = "NC")
} # }
```
