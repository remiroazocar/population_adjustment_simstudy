

# Methods for population adjustment with limited access to individual patient data - a review and simulation study: code

### Antonio Remiro-Azócar, Anna Heath, Gianluca Baio
### *remiroantonio@gmail.com*
### *2020*

This repository contains the `R` code used for my paper [Methods for population adjustment with limited access to individual patient data: A review and simulation study][1], co-authored with [Prof. Gianluca Baio][2] and [Prof. Anna Heath][3]. 

## Utilizing the Scripts

In order to use this repository, the user must first download a copy to their local machine. The user must set the working directory to the location where the download was made. To run the pipeline, the user should then open and run the scripts in the following order:

|          Script           | Explanation                                                  |
| :-----------------------: | ------------------------------------------------------------ |
|   `survival_settings.R`   | Specifies the settings of the simulation settings and saves them in `"./survival_settings.RData"` |
|       `gen_data.R`        | Loads the simulation settings and generates the data for the study (saving the data to the `"./Data/"` subdirectory) |
| `population_adjustment.R` | Performs the indirect comparison methods on the simulated data (saving the resulting means and variances to the `"./Results/"` subdirectory) |
|       `analysis.R`        | Processes the results of the simulation study and computes and graphs the relevant performance metrics (the analyses are saved to the `"./Analysis/"` subdirectory) |

In addition, the `functions.R` script contains user-defined MAIC functions, functions to evaluate the performance measures of interest and functions to present simulation results in a nested loop plot by Rücker and Schwarzer 2014.<sup>1</sup> The file `./Analysis/scenarios.csv`  presents the parameter values or settings for each scenario and summarizes the key performance measures/results associated with each, as presented in the paper. 

The `./Example` subdirectory features example `R` code implementing MAIC, STC and the Bucher method on a simulated example (as per Appendix F of the Supplementary Material). 

The data generation process takes about 5 hours and the indirect comparison methods take about 2 hours, using an Intel Core i7-8650 CPU (1.90 GHz) processor. The `doSNOW` package is used to parallelize the performance of the indirect comparison methods, distributing the tasks to different cores of the computer. 

The code presented here was prepared in `RStudio` using `R` version `3.6.3` in a Windows architecture, with 64-bit operating system. The following packages and version were used:

* `doSNOW 1.0.18` used in combination with `foreach()` to start up local clusters that distribute parallel tasks to different cores
* `dplyr 0.8.5` for data manipulation
* `parallel 3.6.3` to detect the number of CPU cores
* `simstudy 0.1.15` draws correlated (or uncorrelated) covariates from a Gaussian copula in the data-generating process
* `survival 3.1.8` to fit the standard (STC) and weighted (MAIC) Cox proportional hazards regressions
* `tidyr 1.0.2` for data manipulation

## References
1. [**Presenting simulation results in a nested loop plot**][4]  Rücker, G., Schwarzer, G.; BMC Med Res Methodol 14, 129 (2014); 
doi:10.1186/1471-2288-14-129

[1]: https://doi.org/10.1002/jrsm.1511
[2]: http://www.statistica.it/gianluca/
[3]: https://sites.google.com/site/annaheathstats/
[4]: https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-14-12
