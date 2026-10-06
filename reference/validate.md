# Validates the accuracy of the calibrated model on out of sample data.

This function uses the quantity, allocation, and configuration
disagreement to validate the model across the landscape using the
parameters from the calibrate function. Ideally the model is calibrated
with 2 or more years of data and validated for the last year or if you
have 6 or more years of data then the model can be validated for the
final 2 years.

## Usage

``` r
validate(config_rds_file)
```

## Arguments

- config_rds_file:

  Path to config file produced when calling \`configuration\`. The
  config file includes all data necessary used to set up c++ PoPS model

## Value

a data frame of statistical measures of model performance.
