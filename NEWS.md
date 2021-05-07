# rpops (development version)

* Exposed and Infected populations can both be present at the start of a simulation

* Internal functions for data handling have switched from the raster package to the terra package (greatly speeds up data handling and preparation).

* Overpopulation function (individuals in areas of high population leave the area and disperse longer distances on average).

* New dispersal kernels added: Uniform, Power-law, Deterministic neighbor, Hyperbolic-Secant, Gamma, Weibull, Normal, and Logistic.

# rpops (1.0.2)

* Multirun and validate draw a parameter set from the distribution for each run as intended

* Exposed populations now exported each time-step

* Fixed error when calibrating for less than 4 parameters

# rpops (1.0.0)

Version 1.0.0 of the _PoPS Core_ C++ library and its interfaces: _rpops_ R package and _r.pops.spread_ GRASS GIS module. The release of rpops includes:

* Susceptible-infected (`SI`) and susceptible-exposed-infected (`SEI`) host phases (`model_type`, `latency_period`).

* Host mortality tracking (`mortality_rate`, `mortality`).

* Host removal and pesticide application treatments (`treatments`, `treatment_date`, `pesticide_duration`).

* Host resistance based on pesticide application treatments (`pesticide_duration` > 0).

* Treatments applied only to a ratio of hosts (`treatment_application`).

* Yearly pest removal based on lethal temperature (`lethal_temperature`, `lethal_month`).

* Two different dispersal kernels (`natural_dispersal_kernel`, `anthropogenic_dispersal_kernel`).

* Cauchy and exponential radial dispersal kernels use Von Mises distribution.

* Seasonal spread (`seasonality` in months).

* Multiple stochastic runs (`pops_multirun`).

* Parallel execution of multiple runs (`number_of_cores` in `pops_multirun`).

* Output of average and standard deviation of infected hosts across multiple runs and for a single stochastic run (`simulation_mean`,  `single_run`,  `simulation_sd`).

* Average and standard deviations for output averages (`number_infecteds`, `infected_areas`).

* Infection probability output in percent (`probability`).

* Spread rate measurement in 4 cardinal directions (`west_rate`, `east_rate`, `south_rate`, `north_rate` ).

* Distance to quarantine in 4 cardinal directions (`north_distance_to_quarantine`, `south_distance_to_quarantine`, `east_distance_to_quarantine`, `west_distance_to_quarantine` ).

* Probability of quarantine escape (`escape_probability`).