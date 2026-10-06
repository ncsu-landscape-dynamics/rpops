# Package index

## Calibration and validation of the model to observed data

Functions for calibrating and validating a model built with PoPS.

- [`calibrate()`](calibrate.md) : Calibrates the reproductive rate and
  dispersal scales of the pops model.
- [`validate()`](validate.md) : Validates the accuracy of the calibrated
  model on out of sample data.

## Running the model

Functions for running a PoPS model.

- [`pops()`](pops.md) : PoPS (Pest or Pathogen Spread) model
- [`pops_simulate()`](pops_simulate.md) : PoPS Simulate

## Internal data handlers

Functions for error checks and data handling for pops setup

- [`configuration()`](configuration.md) : PoPS (configuration
- [`pops_model()`](pops_model.md) : PoPS (Pest or Pathogen Spread) model
- [`quantity_allocation_disagreement()`](quantity_allocation_disagreement.md)
  : Compares quantity and allocation disagreement of two raster data
  sets
- [`create_summary_stats_and_stacks()`](create_summary_stats_and_stacks.md)
  : Create summary stats and raster summaries
