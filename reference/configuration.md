# PoPS (configuration

Function with a single input and output list for parsing, transforming,
and performing all checks for all functions to run the pops c++ model

## Usage

``` r
configuration(config_file, testing = FALSE)
```

## Arguments

- config_file:

  yaml or csv file with paths and data necessary for formatting data to
  be used in the c++ model

- testing:

  only used during tests otherwise ignored

## Value

config list with all data ready for pops C++ or error message
