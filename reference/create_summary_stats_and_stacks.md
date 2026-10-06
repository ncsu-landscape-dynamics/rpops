# Create summary stats and raster summaries

This function takes the outputs from PoPS lite and calculates the mean,
max, min, and standard deviation of each cell and calculates the summary
statistics for the study area.

## Usage

``` r
create_summary_stats_and_stacks(config_rds_file)
```

## Arguments

- config_rds_file:

  Path to config rds file produced when calling \`configuration\`. The
  config file includes all data necessary used to set up c++ PoPS model.
  This is the same file called when calling \`simulate\`.

## Value

creates and writes raster stacks from all pops_lite runs in the outputs
folder. Also reates and writes mean, standard deviation, median, min,
and max runs, and summary statistics.
