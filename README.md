# hydroState
_hydroState_ is an R-package that identifies regime changes in streamflow time-series that is not explained by variations in precipitation. For details of the mathematics, and its use in understanding catchment drought non-recover, see:

Peterson TJ, Saft M, Peel MC & John A (2021), Watersheds may not recover from drought, Science, DOI: 10.1126/science.abd5085

The package allows a flexible set of Hidden Markov Models of annual, seasonal or monthly time-step to be built and which includes precipitation as a predictor of streamflow. Suites of models can be build for a single catchment, ranging from from one to three states and each with differing combinations of error models and auto-correlation terms, allowing the most parsimonious model to easily be identified (by AIC). The entire package is written in R S4 object oriented code with documentation of the user facing functions.

The package development was funded by the Victorian Government The Department of Energy, Environment, and Climate Action (https://www.water.vic.gov.au/). 

The package was developed for the analysis of long-term streamflow data during and after meteorological drought. The package has been applied to 161 unregulated catchments across Victoria (Australia). A subset of this dataset and fitted models is provide for several catchments to demonstrate how to use hydroState. Below are examples for (i) building a single annual model with a pre-defined default model structure and (ii) building 36 different annual models and selecting the most parsimonious model (by AIC) - as was undertaken for the paper. Additional examples including a seasonal analysis are provided in the _vignette folder_(https://github.com/peterson-tim-j/HydroState/tree/master/vignette).

# Example 1. A two-state annual model

This _hydroState_ example builds and calibrates one two-state annual model and plots the results. The plot below shows the input data, the estimated Viterbi states over time and the conditional state probabilities.

![example1](https://github.com/peterson-tim-j/HydroState/tree/master/vignette/annual.two.state/annual.two.state.timeseries.png)

```R
#----------------
rm(list=ls())
library(hydroState)

# Load flow data
data(streamflow_annual_221201)

# Build the default model, 2-state linear model
model = build(input.data = streamflow_annual_221201)

# Fit the model
model = fit.hydroState(model, pop.size.perParameter = 10, max.generations = 500)

# Review residual plots
plot(model, pse.residuals = TRUE, siteID = '221201', do.pdf = FALSE)

# Set the reference year to name the states
model = setInitialYear(model, 1990)

# Plot all state plots (Example 1 figure)
plot(model) 

```

# Example 2. Building 36 models and selecting the most parsimonious model

This _hydroState_ example creates 36 models of varying complexity and then calibrates each one using multiple trials (see paper supplemental materials). This is also further explained with the calibration procedure in the _build.all_ function. Importantly, the calibration of all of the models takes many hours. The plot below shows the input data, the estimated Viterbi states over time and the conditional state probabilities for the best model. Given that the Millennium Drought ended in 2010, the plots shows that the runoff had not recovered by water years starting in March 2016.

![example2](https://github.com/peterson-tim-j/HydroState/tree/master/vignette/adjusted.state.model/gamma.annual.model.timeseries.407211.png)


```R
rm(list=ls())
library(hydroState)

# Load flow data
data(streamflow_annual_407211)

# Build all model combinations with 'a0' and 'std' as state-dependent parameters 
all.models = build.all(input.data = streamflow_annual_407211,
                           state.shift.parameters =c('a0','std'),
                           siteID = "407211")

# Fit all models (takes many hours)
<!-- all.models = fit.hydroState(all.models, pop.size.perParameter = 10, max.generations = 5000, doParallel = T) -->

# get table of AIC values for each model (from existing fitted models)
get.AIC(all.models.annual.fitted.407211)

# get best model with lowest AIC
best.model = all.models.annual.fitted.407211@models[[which(names(all.models.annual.fitted.407211@models) == rownames(get.AIC(all.models.annual.fitted.407211))[1])]]

# Review residual plots
plot(best.model, pse.residuals = TRUE, siteID = '407211', do.pdf = FALSE)

# Set the reference year to name the states
best.model = setInitialYear(best.model, 1991)

# Plot all result plots (Example 2 figure)
plot(best.model) 
```
