#' hydroState: Hidden Markov modeling for hydrological state change
#'
#' \code{hydroState} is an R-package that identifies regime changes in streamflow time-series that is not explained by variations in precipitation.
#'
#' For details of the mathematics, and its use in understanding catchment drought non-recovery, see:
#'
#' Peterson TJ, Saft M, Peel MC & John A (2021), Watersheds may not recover from drought, Science, DOI: 10.1126/science.abd5085
#'
#' The package allows a flexible set of Hidden Markov Models (HMMs) of annual, seasonal or monthly time-step to be built and which includes precipitation as a predictor of streamflow. Suites of models can be build for a single catchment, ranging from from one to three states and each with differing combinations of error models and auto-correlation terms, allowing the most parsimonious model to easily be identified (by AIC). The entire package is written in R S4 object oriented code with user facing functions and vignettes.
#'
#' @section Functions:
#' \itemize{
#'   \item \code{{build()}} — Builds hydroState model
#'   \item \code{{build.all()}} — Builds all hydroState models
#'   \item \code{{fit.hydroSata()}} — Fit hydroState model(s)
#'   \item \code{{setInitialYear()}} — Sets state names given initial year
#'   \item \code{{plot()}} — Plot states or pseudo residuals over time
#'   \item \code{{get.residuals}} — Get pseudo residuals
#'   \item \code{{get.states}} — Get states
#'   \item \code{{check()}} — Check reliability of state predictions
#'   \item \code{{get.AIC()}} — Get AIC
#'   \item \code{{get.seasons}} — Get pseudo residuals
#' }
#'
#' @section Vignettes:
#' \itemize{
#'   \item {vignette("hydroState", package = "hydroState")}
#'   \item {vignette("Adjust the default state model", package = "hydroState")}
#'   \item {vignette("Seasonal and monthly models", package = "hydroState")}
#' }
#'
#' @section Usage:
#' The package contains a default \code{hydroState} model object that explains streamflow as a function of precipitation using a linear model. Once the model object is built (\code{{build()}}), the model is fitted (\code{{fit.hydroState()}}) to determine the most likely rainfall-runoff state at each time-step. To assess the adequacy of the fit, the residuals are plotted (\code{{plot()}}), and an adequate fit requires the residuals to be normally distributed, uniform, with minimal correlation and minimal trends. The resulting runoff states from the fitted model can then be evaluated over time (\code{{plot()}}) and even exported ({\code{get.states()}}) with the state values, confidence intervals, and conditional probabilities at each time-step. Input data requires a dataframe with catchment average runoff and precipitation at annual, seasonal, or monthly timesteps, and gaps with missing data are permitted. An example of this workflow with the default model is demonstrated within the {vignette("hydroState", package = "hydroState")} vignette.
#'
#' To better explain the rainfall-runoff relationship, the default model can be adjusted by selecting various items within the \code{{build()}} function. These include:
#'  \itemize{
#'   \item{\code{data.transform}:} transform streamflow observations in order to reduce skew: 'boxcox', 'log', 'burbidge', or 'none'
#'   \item{\code{parameters}:} account for auto-correlation through including the degree of auto-correlation: 'AR1', 'AR2', or 'AR3'
#'    \item{\code{state.shift.parameters}:} assign either the intercept, slope, or auto-correlation parameter as state dependent parameter
#'    \item{\code{error.distribution}:} adjust the error distribution with 'normal', 'gamma', or 'truc.normal'
#'    \item{\code{seasonal.parameters}:} account for intra-annual varation within the rainfall-runoff relationship
#'    \item{\code{transition.graph}:} set the number of possible states in the model (1, 2, or 3).
#'  }
#'
#' An example of how to adjust the default model is demonstrated within the {vignette("Adjust the default state model", package = "hydroState")}
#'   and {vignette("Seasonal and monthly models", package = "hydroState")}.
#'
#' There is an additional option to construct all possible types of models using the \code{{build.all()}}, and compare them using the same \code{{fit.hydroState()}} function. The most likely model can be selected based on the AIC where the best model will have the lowest AIC. An example of this is demonstrated at the end of the {vignette("Adjust the default state model", package = "hydroState")}
#'   and {vignette("Seasonal and monthly models", package = "hydroState")}. To get stated, it is recommended to evaluate the default model at first with one state and again with two states.
#'
#'
#' @section Acknowledgments:
#'  The package development was funded by the Victorian Government The Department of Energy, Environment, and Climate Action (https://www.water.vic.gov.au/).
#'
#'
#'
#' @docType package
#' @name hydroState-package
"_PACKAGE"
