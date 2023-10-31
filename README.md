# PS_survival
Simulation codes for Principal Stratification Analysis of Noncompliance with Time-to-Event Outcomes

## File Structure

- `Simulation util.R` contains utility functions used in the simulation studies, including drawing samples, calculating survival probability and RACE and creating plots.
- `gen_datasets.R` generates the datasets.
- `simulation_x_codes.R` runs the corresponding simulation and makes plots as they appear in the paper.

## How to run

Before running the codes, make sure that you have `rstan` configured and package `PStrata` installed.

Guide to get started with RStan: [RStan Website](mc-stan.org/users/interfaces/rstan)

To install PStrata, please run the following script in `R`.
```{r}
install.packages("PStrata")
```

In cases when the package cannot be installed in the above way, please try installing it from the GitHub repository:
```{r}
install.packages("devtools")
devtools::install_github("LauBok/PStrata")
```

To replicate simulation 1, run the following script (change 1 to 2 or 3 for the other two simulations):
```{r}
source("simulation_1_code.R")
```
