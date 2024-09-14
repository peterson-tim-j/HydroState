# hydroState
_hydroState_ is an R-package that identifies regime changes in streamflow time-series that is not explained by variations in precipitation. For details of the mathematics, and its use in understanding catchment drought non-recover, see:

Peterson TJ, Saft M, Peel MC & John A (2021), Watersheds may not recover from drought, Science, DOI: 10.1126/science.abd5085

The package allows a flexible set of Hidden Markov Models of annual, seasonal or monthly time-step to be built and which includes precipitation as a predictor of streamflow. Suites of models can be build for a single catchment, ranging from from one to three states and each with differing combinations of error models and auto-correlation terms, allowing the most parsimonious model to easily be identified (by AIC). The entire package is written in R S4 object oriented code. Future versions will include documentation of the functions, classes and methods.

The package development was funded by the Victorian Government The Department of Environment, Land, Water and Planning _Climate and Water Initiate_ (https://www.water.vic.gov.au/climate-change/research/vicwaci). 

The package was developed for the analysis of long-term streamflow data during and after meteorological drought. The package has been applied to 161 unregulated catchments across Victoria (Australia). Stream gauge control information and annual and monthly observed data is included in the package and is available using the following R commands:
```R
data(streamflow_annual)
data(streamflow_monthly)
data(streamflow_gaugeInfo)
data(fitted.model)
```

Below are examples for (i) building a single annual model with a pre-defined model structure and (ii) building 64 different annuals models and selecting the most parsimonious model (by AIC) - as was undertaken for the paper. Additional examples are provided in the _vignette folder_(https://github.com/peterson-tim-j/HydroState/tree/master/vignette).

# Example 1. A two-state annual model

This _hydroState_ example builds and calibrates one two-state annual model and plots the results. The plot below shows the input data, the estimated Viterbi states over time and the states with the cumulative rainfall residual (D) and the conditional state probabilities (E).

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

# Example 2. Building 64 models and selecting the most parsimonious model

This _hydroState_ example creates 64 models of varying complexity and then calibrates each one using multiple trials (see paper supplemental materials). Importantly, the calibration of all of the models takes many hours. The plot below shows the input data, the estimated Viterbi states over time and the states with the cumulative rainfall residual (D) and the conditional state probabilities (E). Given that the Millennium Drought ended in 2010, the plots shows that the runoff had not recovered by water years starting in March 2016.

![example2](https://user-images.githubusercontent.com/8623994/118099548-bb542200-b418-11eb-821e-a9817f170e1e.png)


```R
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

# Build all combinations of annual models for this gauge.
all.Models <- new('hydroState.allModels',as.character(gaugeID), streamflow_annual, allow.flickering=F)

# Calibrate (using MLW) each of the models.
all.Models <- fit(all.Models, pop.size.perParameter=10, max.generations=500, doParallel=F)

# Select the best model (byt AIC)
best.model = getAIC.bestModel(all.Models)

# Name the states names with 1990 being defined as a 'normal' runoff year.
best.model <- setStateNames(best.model, 1990)

# Plot Viterbi states
viterbi(best.model)

# Plot pseduo residuals
check.PseudoResiduals(best.model)
```
