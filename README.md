# hydroState

![CRAN/METACRAN Version](https://img.shields.io/cran/v/hydroState?style=flat)


### Overview
hydroState is an R-package that identifies regime changes in streamflow time-series that is not explained by variations in precipitation. For details of the mathematics, and its use in understanding catchment drought non-recovery, see:

Peterson TJ, Saft M, Peel MC & John A (2021), Watersheds may not recover from drought, Science, DOI: 10.1126/science.abd5085

The package allows a flexible set of Hidden Markov Models of annual, seasonal or monthly time-step to be built and which includes precipitation as a predictor of streamflow. Suites of models can be build for a single catchment, ranging from from one to three states and each with differing combinations of error models and auto-correlation terms, allowing the most parsimonious model to easily be identified (by AIC). The entire package is written in R S4 object oriented code with documentation of the user facing functions. 

### Installation

You can install the development version of hydroState from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("peterson-tim-j/HydroState")
```

You can also get the official release version from CRAN _soon_:

``` r
install.packages("hydroState")
```
### Usage

Be sure to see the [Getting Started](articles/hydroState.html) article for an example of this workflow. The package contains a default hydroState model object that explains streamflow as a function of precipitation using a linear model. Once the model object is built ([build()](reference/build.html)), the model is fitted ([fit.hydroState()](reference/fit.hydroState.html)) to determine the most likely rainfall-runoff state at each time-step. To assess the adequacy of the fit, the residuals are plotted ([plot()](reference/plot.hydroState.html)), and an adequate fit requires the residuals to be normally distributed, uniform, with minimal correlation and minimal trends. The resulting runoff states from the fitted model can then be evaluated over time ([plot()](reference/plot.hydroState.html)) and even exported ([get.states()](reference/get.states.html)) with the state values, confidence intervals, and conditional probabilities at each time-step. Input data requires a dataframe with catchment average runoff and precipitation at annual, seasonal, or monthly timesteps, and gaps with missing data are permitted. 

#### Adjusting the default model

To better explain the rainfall-runoff relationship, the default model can be adjusted by selecting various items within the [build()](reference/build.html) function. These include:

* `data.transform`: transform streamflow observations in order to reduce skew: `boxcox`, `log`, `burbidge`, or `none`
* `parameters`: account for auto-correlation through including the degree of auto-correlation: `AR1`, `AR2`, or `AR3`
* `state.shift.parameters`: assign either the intercept, slope, standard deviation, or auto-correlation parameter as state dependent parameter
* `error.distribution`: adjust the error distribution with `normal`, `gamma`, or `truc.normal`
* `seasonal.parameters`: account for intra-annual varation within the rainfall-runoff relationship
* `transition.graph`: set the number of possible states in the model (1, 2, or 3).

An example of how to adjust the default model is demonstrated within the [Adjust the default state model](articles/adjust.state.model.html) article and [Seasonal and monthly models](articles/subAnnual.models.html) article. 

#### Finding the best model

There is an additional option to construct all possible types of models using the [build.all()](reference/build.all.html) function, and compare them using the same [fit.hydroState()](reference/fit.hydroState.html) function. The most likely model can be selected based on the AIC where the best model will have the lowest AIC. An example of this is demonstrated at the end of the [Adjust the default state model](articles/adjust.state.model.html) article and [Seasonal and monthly models](articles/subAnnual.models.html) article. To get stated, it is recommended to evaluate the default model at first with one state and again with two states as shown in the [Getting Started](articles/hydroState.html) article.

### How to cite hydroState:

Peterson TJ, Saft M, Peel MC & John A (2021), Watersheds may not recover from drought, Science, DOI: 10.1126/science.abd5085. 

### Acknoledgement

The package development was funded by the Victorian Government The Department of Energy, Environment, and Climate Action (https://www.water.vic.gov.au/). 

