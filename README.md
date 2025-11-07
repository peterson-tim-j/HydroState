# hydroState

[![CRAN/METACRAN Version](https://img.shields.io/cran/v/hydroState?logo=R&logoColor=blue&label=CRAN%20version)](https://cran.r-project.org/package=hydroState)
[![](http://cranlogs.r-pkg.org/badges/grand-total/hydroState)](https://cran.r-project.org/package=hydroState)
[![](http://cranlogs.r-pkg.org/badges/last-month/hydroState)](https://cran.r-project.org/package=hydroState)
![Static Badge](https://img.shields.io/badge/Scence_Magazine-red?logo=Clarivate&label=Top%201%25%20Most%20Cited%20in&link=https%3A%2F%2Fdoi.org%2F10.1126%2Fscience.abd5085)
![](https://img.shields.io/github/stars/peterson-tim-j/HydroState?style=social&labelColor=yellow&color=yellow)
[![GitHub forks](https://img.shields.io/github/forks/peterson-tim-j/HydroState)](https://github.com/peterson-tim-j/HydroState/network)

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

Be sure to see the [Getting Started](https://peterson-tim-j.github.io/HydroState/articles/hydroState.html) article for an example of this workflow. The package contains a default hydroState model object that explains streamflow as a function of precipitation using a linear model. Once the model object is built ([build()](https://peterson-tim-j.github.io/HydroState/reference/build.html)), the model is fitted ([fit.hydroState()](https://peterson-tim-j.github.io/HydroState/reference/fit.hydroState.html)) to determine the most likely rainfall-runoff state at each time-step. To assess the adequacy of the fit, the residuals are plotted ([plot()](https://peterson-tim-j.github.io/HydroState/reference/plot.hydroState.html)), and an adequate fit requires the residuals to be normally distributed, uniform, with minimal correlation and minimal trends. The resulting runoff states from the fitted model can then be evaluated over time ([plot()](https://peterson-tim-j.github.io/HydroState/reference/plot.hydroState.html)) and even exported ([get.states()](https://peterson-tim-j.github.io/HydroState/reference/get.states.html)) with the state values, confidence intervals, and conditional probabilities at each time-step. Input data requires a dataframe with catchment average runoff and precipitation at annual, seasonal, or monthly timesteps, and gaps with missing data are permitted. 

#### Adjusting the default model

To better explain the rainfall-runoff relationship, the default model can be adjusted by selecting various items within the [build()](https://peterson-tim-j.github.io/HydroState/reference/build.html) function. These include:

* `data.transform`: transform streamflow observations in order to reduce skew: `boxcox`, `log`, `burbidge`, or `none`
* `parameters`: account for auto-correlation through including the degree of auto-correlation: `AR1`, `AR2`, or `AR3`
* `state.shift.parameters`: assign either the intercept, slope, standard deviation, or auto-correlation parameter as state dependent parameter
* `error.distribution`: adjust the error distribution with `normal`, `gamma`, or `truc.normal`
* `seasonal.parameters`: account for intra-annual varation within the rainfall-runoff relationship
* `transition.graph`: set the number of possible states in the model (1, 2, or 3).

An example of how to adjust the default model is demonstrated within the [Adjust the default state model](https://peterson-tim-j.github.io/HydroState/articles/adjust.state.model.html) article and [Seasonal and monthly models](https://peterson-tim-j.github.io/HydroState/articles/subAnnual.models.html) article. 

#### Finding the best model

There is an additional option to construct all possible types of models using the [build.all()](https://peterson-tim-j.github.io/HydroState/reference/build.all.html) function, and compare them using the same [fit.hydroState()](https://peterson-tim-j.github.io/HydroState/reference/fit.hydroState.html) function. The most likely model can be selected based on the AIC where the best model will have the lowest AIC. An example of this is demonstrated at the end of the [Adjust the default state model](https://peterson-tim-j.github.io/HydroState/articles/adjust.state.model.html) article and [Seasonal and monthly models](https://peterson-tim-j.github.io/HydroState/articles/subAnnual.models.html) article. To get stated, it is recommended to evaluate the default model at first with one state and again with two states as shown in the [Getting Started](https://peterson-tim-j.github.io/HydroState/articles/hydroState.html) article.

### How to cite hydroState:

Peterson TJ, Saft M, Peel MC & John A (2021), Watersheds may not recover from drought, Science, DOI: 10.1126/science.abd5085. 

### Acknowledgments:

The package development was funded by the Victorian Government The Department of Energy, Environment, and Climate Action (https://www.water.vic.gov.au/). 

