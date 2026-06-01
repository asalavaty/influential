# Conditional probability of deviation from means

This function calculates the conditional probability of deviation of two
centrality measures (or any two other continuous variables) from their
corresponding means in opposite directions.

## Usage

``` r
cond.prob.analysis(data, nodes.colname, Desired.colname, Condition.colname)
```

## Arguments

- data:

  A data frame containing the values of two continuous variables and the
  name of observations (nodes).

- nodes.colname:

  The character format (quoted) name of the column containing the name
  of observations (nodes).

- Desired.colname:

  The character format (quoted) name of the column containing the values
  of the desired variable.

- Condition.colname:

  The character format (quoted) name of the column containing the values
  of the condition variable.

## Value

A list of two objects including the conditional probability of deviation
of two centrality measures (or any two other continuous variables) from
their corresponding means in opposite directions based on both the
entire network and the split-half random sample of network nodes.

## See also

Other centrality association assessment functions:
[`double.cent.assess()`](https://asalavaty.github.io/influential/reference/double.cent.assess.md),
[`double.cent.assess.noRegression()`](https://asalavaty.github.io/influential/reference/double.cent.assess.noRegression.md)

## Examples

``` r
if (FALSE) { # \dontrun{
MyData <- centrality.measures
My.conditional.prob <- cond.prob.analysis(data = MyData,
                                          nodes.colname = rownames(MyData),
                                          Desired.colname = "BC",
                                          Condition.colname = "NC")
                                          } # }
```
