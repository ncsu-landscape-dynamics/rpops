# Lists for checks and tests
# calibration success metric option
success_metric_options <- c("quantity", "allocation", "configuration", "quantity and allocation",
                            "quantity and configuration", "allocation and configuration",
                            "quantity, allocation, and configuration", "accuracy", "precision",
                            "recall", "specificity", "accuracy and precision",
                            "accuracy and specificity", "accuracy and recall",
                            "precision and recall", "precision and specificity",
                            "recall and specificity", "accuracy, precision, and recall",
                            "accuracy, precision, and specificity",
                            "accuracy, recall, and specificity",
                            "precision, recall, and specificity",
                            "accuracy, precision, recall, and specificity",
                            "rmse", "distance", "mcc", "mcc and quantity", "mcc and distance",
                            "rmse and distance", "mcc and configuration", "mcc and RMSE",
                            "mcc, quantity, and configuration")

quantity_list <- c("quantity", "quantity and allocation", "quantity and configuration",
                   "quantity, allocation, and configuration", "mcc and quantity",
                   "mcc, quantity, and configuration")

allocation_list <- c("allocation", "quantity and allocation", "allocation and configuration",
                     "quantity, allocation, and configuration")

configuration_list <- c("configuration",  "quantity and configuration",
                        "allocation and configuration", "quantity, allocation, and configuration",
                        "mcc and configuration", "mcc, quantity, and configuration")

accurracy_list <- c("accuracy", "accuracy and precision", "accuracy and specificity",
                    "accuracy and recall", "accuracy, precision, and recall",
                    "accuracy, precision, and specificity",
                    "accuracy, recall, and specificity",
                    "accuracy, precision, recall, and specificity")

precision_list <- c("precision", "accuracy and precision", "precision and recall",
                    "precision and specificity", "accuracy, precision, and recall",
                    "accuracy, precision, and specificity",
                    "precision, recall, and specificity",
                    "accuracy, precision, recall, and specificity")

recall_list <- c("recall", "accuracy and recall", "precision and recall", "recall and specificity",
                 "accuracy, precision, and recall", "accuracy, recall, and specificity",
                 "precision, recall, and specificity",
                 "accuracy, precision, recall, and specificity")

specificity_list <- c("specificity", "accuracy and specificity", "precision and specificity",
                      "recall and specificity", "accuracy, precision, and specificity",
                      "accuracy, recall, and specificity",  "precision, recall, and specificity",
                      "accuracy, precision, recall, and specificity")

rmse_list <- c("rmse", "rmse and distance", "mcc and RMSE")

distance_list <- c("distance", "mcc and distance", "rmse and distance")

mcc_list <- c("mcc", "mcc and quantity", "mcc and distance", "mcc and configuration",
              "mcc and RMSE", "mcc, quantity, and configuration")

weather_type_list <- c("deterministic", "probabilistic", "none")

# Check for correct kernel options
kernel_list <- c(
  "cauchy",
  "Cauchy",
  "exponential",
  "Exponential",
  "uniform",
  "Uniform",
  "deterministic neighbor",
  "deterministic-neighbor",
  "power law",
  "power-law",
  "Power-law",
  "Power-Law",
  "Power Law",
  "Power law",
  "hyperbolic secant",
  "hyperbolic-secant",
  "Hyperbolic-secant",
  "Hyperbolic-Secant",
  "Hyperbolic secant",
  "Hyperbolic Secant",
  "gamma",
  "Gamma",
  # "exponential power",
  # "exponential-power",
  # "Exponential-power",
  # "Exponential-Power",
  # "Exponential power",
  "weibull",
  "Weibull",
  "normal",
  # "log normal",
  # "log-normal",
  # "Log-normal",
  # "Log-Normal",
  # "Log normal",
  # "Log Normal",
  "logistic",
  "Logistic",
  "network",
  "Network"
)

output_list <- c("all_simulations", "summary_outputs", "None")
output_write_list <- c("all_simulations", "summary_outputs")

si_list <- c("SEI", "susceptible-exposed-infected", "susceptible_exposed_infected",
             "Susceptible-Exposed-Infected", "Susceptible_Exposed_Infected")

sei_list <- c("SI", "susceptible-infected", "susceptible_infected", "Susceptible-Infected",
              "Susceptible_Infected")

treatment_list <- c("ratio", "all infected")
network_movement_options <- c("walk", "jump", "teleport")
aws_bucket_list <- c("casestudy_creation", "model_api")
parallel_function_list <- c("validate", "multirun", "sensitivity")
parameter_draw_list <- c("validate", "pops", "multirun", "sensitivity", "casestudy_creation")
val_cal_list <- c("validate", "calibrate")
raster_list <- c("grd", "tif", "img", "vrt")
failed_check_list <- c("checks_passed", "failed_check")
output_frequency_list <-
  c("week", "month", "day", "year", "time_step", "every_n_steps", "final_step")
csv_list <- c("csv", "txt")
pest_host_table_colnames <- c("host", "susceptibility_mean", "susceptibility_sd",
                          "mortality_rate_mean", "morality_rate_sd", "mortality_time_lag")
competency_table_colnames <- c("competency_mean", "competency_sd")
