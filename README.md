# HydroState
_HydroState_ is an R-package that idenifies regime changes in streamflow time-series that is not explained by variations in precipitation. The package allows a flexible set of Hidden Markov Models of annual, seasonal or monthly time-step to be built and which includes precipitation as a predictor of streamflow. Suites of models can be build for a single catchment, ranging from from one to three states and each with differing combinations of error models and auto-correlation terms, allowing the most parsimineous model to easily be idenified (by AIC). The entire package is written in R S4 object oriented code. Future versions will include documentation of the funcations, classes and methods.

The package was developed for the analysis of long-term streamflow data during and after meteorlogical drought. The package development was funded by the Victorian Government Climate and Water Initiate (https://www.water.vic.gov.au/climate-change/research/vicwaci). The package was developed and applied to 162 unregulated catchments across Victoria (Australia). The annual and monthly observed data is included in the package and is available using the following R commands:
```R
data(streamflow_annual)
data(streamflow_monthly)
```

# Example 1. A two-state annual model

This _HydroState_ example builds and calibrates one two-state annual model and plots the results. The plot below shows the input data, the estimated Viterbi states over time and the states with the cumulative rainfall residual (D) and the conditional state probabilities (E).

![example1](https://user-images.githubusercontent.com/8623994/90325344-c24f2080-dfbd-11ea-956e-c64b7820342f.png)


```R
#----------------
rm(list=ls())
library(hydroState)
library(DEoptim)
library(truncnorm)

# Load flow data
data(streamflow_annual)

# Extract one catchment
gaugeID = 221201;
streamflow_annual = streamflow_annual[streamflow_annual$gauge==gaugeID,]

# Convert to format for hydroState
streamflow_annual = data.frame(year = streamflow_annual$hy_year, flow=streamflow_annual$q, precipitation=streamflow_annual$p)
filt = streamflow_annual$year>=1900
streamflow_annual = streamflow_annual[filt,]

# Build input objects.
# Note a 2-state linear model with first-order serial correlation is built. A truncated normal distribution is used for
# each Markov state.
transition.graph=matrix(TRUE,2,2)
Qhat = new('Qhat.log', input.data=streamflow_annual)
QhatModel = new('QhatModel.homo.normal.linear.AR1', input.data=streamflow_annual, transition.graph=transition.graph)
markov = new('markov.annualHomogeneous', transition.graph=transition.graph)

# Build HydroState object
model = new('hydroState',input.data=streamflow_annual, Qhat.object=Qhat, QhatModel.object=QhatModel, markov.model.object=markov)

# Fit the model using DEoptim
model <- hydroState::fit(model,pop.size.perParameter = 10, max.generations=500)

# Name the states names with 1990 being defined as a 'norma' runoff year.
model <- setStateNames(model, 1990)

# Plot Viterbi states
viterbi(model)

# Plot pseduo residuals
check.PseudoResiduals(model)
```
