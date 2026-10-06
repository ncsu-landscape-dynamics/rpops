# PoPS (Pest or Pathogen Spread) model

A dynamic species distribution model for pest or pathogen spread in
forest or agricultural ecosystems. The model is process based meaning
that it uses understanding of the effect of weather and other
environmental factors on reproduction and survival of the pest/pathogen
in order to forecast spread of the pest/pathogen into the future. This
function performs a single stochastic realization of the model and is
predominantly used for automated tests of model features.

## Usage

``` r
pops(config)
```

## Arguments

- config:

  Path to config file produced when calling \`configuration\`. The
  config file includes all data necessary used to set up c++ PoPS model

## Value

list of infected and susceptible per year
