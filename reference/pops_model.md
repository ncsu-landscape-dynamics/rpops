# PoPS (Pest or Pathogen Spread) model

Wrapper for pops_model_cpp only used internally in other functions.

## Usage

``` r
pops_model(config)
```

## Arguments

- config:

  a list of all of the necessary data, parameters, files, and paths
  formatted for the C++ model created from configuration and modified
  slightly by the other functions for each simulation.

## Value

list of vector matrices of infected and susceptible hosts per simulated
year and associated statistics (e.g. spread rate)
