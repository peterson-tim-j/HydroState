\name{hydroState-package}

\alias{hydroState-package}

\alias{hydroState}

\docType{package}

\title{Overview of methods and procedures}

\description{

\code{hydroState} provides methods to construct and evaluate hidden Markov models (HMM) of annual, seasonal, or monthly streamflow runoff with precipitation as a predictor. The state of the relationship between these observations is evaluated overtime. The package contains a default \code{hydroState} model to evaluate probable shifts in the intercept of the rainfall-runoff relationship or shifts in other terms can be evaluated. The default model can also be expanded with various terms to better define the relationship as discussed below, but the general workflow is as follows. Once the default or constructed model is built (\code{buildModel}), the model is fitted (\code{fitModel}) to determine the most likly rainfall-runoff state at each time-step. To assess the adequacy of the fit, the residuals are plotted (\code{plot.residuals}), and an adequate fit requires the residuals to be normaly distributed and not have any correlation or trends. The resulting runoff states from the fitted model can then be evaluated overtime (\code{plot.states}) and even exported (\code{get.states}) with the state values, confidence intervals, and conditional probabilities at each time-step. Input data requires a dataframe with flow and precipitation at annual, seasonal, or monthly timesteps, and gaps with missing data are permitted. An example of this workflow with the default model is demonstrated within the \code{annual.two.state} vignette ([annual.two.state](https://github.com/peterson-tim-j/HydroState/tree/master/vignettes/annual.two.state/annual.two.state.html))

To better explain the rainfall-runoff relationship, the default model can be adjusted by selecting various items within the \code{buildModel} function. These include: \code{data.transform}, \code{parameters}, \code{state.shift.parameters}, \code{error.distribution}, and \code{transition.graph}. The \code{data.transform} offers various transformaitons of the observations in order to reduce skew. These transformaitons include: 'boxcox', 'log', 'burbidge', or 'none'. Auto-correlation can be acconted for through including the degree of auto-correlation ('AR1', 'AR2', or 'AR3') in \code{parameters}. The default model only evaluates a state shift in the intercept of the relationship, but the slope and/or auto-correlation terms can also be investigated with assigning them as \code{state.shift.parameters}. In order to better explain the error distribution, the \code{error.distribution} allows selection of any of the following distributions: 'normal', 'gamma', or 'truc.normal'. When monthly or seasonal data is included, intra-annual varation within the rainfall-runoff relationship can be accounted for with assigning parameters as \code{seasonal.parameters}. Lastly, adjusting the \code{transition.graph} sets the number of possible states in the model. An example of how to adjust the default model is demonstrated within the \code{adjust.state.model} vignette ([adjust.state.model](https://github.com/peterson-tim-j/HydroState/tree/master/vignettes/adjust.state.model/adjust.state.model.html)) and

There is an additional option to construct all possible types of models using the \code{buildModelAll}, and compare them using the same \code{fitModel} function. The most likely model is selected based on the AIC where the best model will have the lowest AIC. An example of this is demonstrated at the end of the \code{adjust.state.model} vignette ([adjust.state.model](https://github.com/peterson-tim-j/HydroState/tree/master/vignettes/adjust.state.model/adjust.state.model.html)) and the \code{subAnnual.models} vignette ([subAnnual.models](https://github.com/peterson-tim-j/HydroState/tree/master/vignettes/subAnnual.models/subAnnual.models.html)). To begin, it is recommended to evaluate the default model at first with one state and again with two states. This documentation is organized with four sections that define the general workflow for using \code{hydroState}: Build - Fit - Review - Evaluate.


-----------------------------------------------------------------------------------
}

  \section{I. Build}{

    \tabular{ll}{
      \code{\link{buildModel}}\tab build hydroState model\cr
      \code{\link{buildModelAll}}\tab build all possible models\cr
      \code{\link{showModelAll}}\tab show reference table of all modelss\cr
      ------------------------\tab------------------------------------------------
  }
  }
   \section{II. Fit}{

    \tabular{ll}{
      \code{\link{fitModel}}\tab fit built hydroState model(s)\cr
      ------------------------\tab------------------------------------------------
  }
  }
  \section{III. Review}{

    \tabular{ll}{
      \code{\link{plot.residuals}}\tab plot residuals\cr
      \code{\link{get.residuals}}\tab get residuals\cr
      ------------------------\tab------------------------------------------------
  }
  }
  \section{IV. Evaluate}{

    \tabular{ll}{
      \code{\link{setInitialYear}}\tab set initial year for assigning state names\cr
      \code{\link{plot.states}}\tab plot states\cr
      \code{\link{get.states}}\tab get states\cr
      ------------------------\tab------------------------------------------------
  }
  }

\section{Authors}{
Except where indicated otherwise, the methods and functions in this package were written by Tim Peterson and Thomas Westfall.
}

\section{Acknowledgments}{
...

}


\keyword{ package }
\keyword{ rainfall-runoff }
\keyword{ states }
