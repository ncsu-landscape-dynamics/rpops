# Compares quantity and allocation disagreement of two raster data sets

Uses quantity and allocation disagreement metrics by Pontius and
Millones (2014) and omission and commission errors of comparing a
modeled raster data set to a reference raster data set.

## Usage

``` r
quantity_allocation_disagreement(
  reference,
  comparison,
  use_configuration = FALSE,
  mask = NULL,
  use_distance = FALSE
)
```

## Arguments

- reference:

  the raster with ground truth data. For all metrics expect RMSE these
  values are reclassified with values \> 1 becoming 1, values \< 1 going
  to 0, and NA values staying NA.

- comparison:

  the raster with simulated data. For all metrics expect RMSE these
  values are reclassified with values \> 1 becoming 1, values \< 1 going
  to 0.

- use_configuration:

  Boolean if you want to use configuration disagreement for comparing
  model runs. Default is FALSE.

- mask:

  Used to provide a mask to remove 0's that are not true negatives from
  comparisons.

- use_distance:

  Boolean if you want to compare distance between simulations and
  observations. Default is FALSE.

## Value

A data frame with spatial configuration metrics. Particularly quantity,
allocation, and total disagreement, omission and commission, and
directional disagreement where directional disagreement.
