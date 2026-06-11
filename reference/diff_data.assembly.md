# Assembling the differential/regression data

This function assembles a dataframe required for running the
***`ExIR`*** model. You may provide as many differential/regression data
as you wish. Also, the datasets should be filtered beforehand according
to your desired thresholds and, consequently, should only include the
significant data. Each dataset provided should be a dataframe with one
or two columns. The first column should always include
differential/regression values and the second one (if provided) the
significance values. Please also note that the significance (adjusted
P-value) column is mandatory for differential datasets.

## Usage

``` r
diff_data.assembly(...)
```

## Arguments

- ...:

  Desired datasets/dataframes.

## Value

A dataframe including the collective list of features in rows and all of
the differential/regression data and their statistical significance in
columns with the same order provided by the user.

## See also

[`exir`](https://asalavaty.github.io/influential/reference/exir.md)

## Examples

``` r
if (FALSE) { # \dontrun{
my.Diff_data <- diff_data.assembly(Differential_data1,
                                   Differential_data2,
                                   Regression_data1)
} # }
```
