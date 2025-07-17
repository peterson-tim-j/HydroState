\name{hydroState-package}

\alias{hydroState-package}

\alias{hydroState}

\docType{package}

\title{Overview of methods and procedures}

\description{

\code{hydroState} is an R-package that identifies regime changes in streamflow time-series that is not explained by variations in precipitation. For details of the mathematics, and its use in understanding catchment drought non-recovery, see:

Peterson TJ, Saft M, Peel MC & John A (2021), Watersheds may not recover from drought, Science, DOI: 10.1126/science.abd5085

The package allows a flexible set of Hidden Markov Models (HMMs) of annual, seasonal or monthly time-step to be built and which includes precipitation as a predictor of streamflow. Suites of models can be build for a single catchment, ranging from from one to three states and each with differing combinations of error models and auto-correlation terms, allowing the most parsimonious model to easily be identified (by AIC). The entire package is written in R S4 object oriented code and is accesible with the following documentation of specific functions and vignettes. This documentation is organized with the general workflow: Build - Fit - Review - Evaluate.

The package contains a default \code{hydroState} model object that explains streamflow as a function of precipitation using a linear model. Once the model object is built (\code{build}), the model is fitted (\code{fit.hydroState}) to determine the most likly rainfall-runoff state at each time-step. To assess the adequacy of the fit, the residuals are plotted (\code{plot}), and an adequate fit requires the residuals to be normaly distributed, uniform, with minimal correlation and minimal trends. The resulting runoff states from the fitted model can then be evaluated over time (\code{plot}) and even exported (\code{get.states}) with the state values, confidence intervals, and conditional probabilities at each time-step. Input data requires a dataframe with catchment average runoff and precipitation at annual, seasonal, or monthly timesteps, and gaps with missing data are permitted. An example of this workflow with the default model is demonstrated within the \code{vignette("annual.two.state", package = "hydroState")}.

To better explain the rainfall-runoff relationship, the default model can be adjusted by selecting various items within the \code{build} function. These include:
\itemize{
\item{\code{data.transform}:} transform streamflow observations in order to reduce skew: 'boxcox', 'log', 'burbidge', or 'none'
\item{\code{parameters}:} account for auto-correlation through including the degree of auto-correlation: 'AR1', 'AR2', or 'AR3'
\item{\code{state.shift.parameters}:} assign either the intercept, slope, or auto-correlation parameter as state dependent parameter
\item{\code{error.distribution}:} adjust the error distribution with 'normal', 'gamma', or 'truc.normal'
\item{\code{seasonal.parameters}:} account for intra-annual varation within the rainfall-runoff relationship
\item{\code{transition.graph}:} set the number of possible states in the model (1, 2, or 3).
}
An example of how to adjust the default model is demonstrated within the \code{vignette("adjust.state.model", package = "hydroState")} and \code{vignette("subAnnual.models", package = "hydroState")}

There is an additional option to construct all possible types of models using the \code{build.all}, and compare them using the same \code{fit.hydroState} function. The most likely model can be selected based on the AIC where the best model will have the lowest AIC. An example of this is demonstrated at the end of the \code{vignette("adjust.state.model", package = "hydroState")} and \code{vignette("subAnnual.models", package = "hydroState")}. To begin, it is recommended to evaluate the default model at first with one state and again with two states.
}

\section{Authors}{
Except where indicated otherwise, the methods and functions in this package were written by Tim Peterson and Thomas Westfall.
}

\section{Acknowledgments}{
The package development was funded by the Victorian Government The Department of Environment, Land, Water and Planning Climate and Water Initiate (https://www.water.vic.gov.au/climate-change/research/vicwaci).
}


\keyword{ package }
\keyword{ rainfall-runoff }
\keyword{ states }
