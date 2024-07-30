\name{hydroState-package}

\alias{hydroState-package}

\alias{hydroState}

\docType{package}

\title{Description of the methods in the hydroState package}

\description{

\code{hydroState} provides methods to construct and evaluate Hidden Markov Models of annual, seasonal or monthly streamflow with precipitation as a predictor. The state of the relationship between these observations is evaluated overtime. The user has options for building the HMM models such as: tansforming streamflow (boxcox, burbidge, or log), selecting the number of states (one to three), residual distribution (normal, gamma, skewed), and degree of auto-correlation (one to three). Details for constructing models with these options are described below. The package allows for the evaluation of multiple models through finding the most likely model by AIC, and the results, the most likely states, can be visualized overtime in a timeseries plot with the observations.


This package is written in R S4 object oriented code. \code{hydroState} is the main class or object. This object requires input from three additional classes in order to build a 'hydroState' model. These classes are known as \code{Qhat}, \code{QhatModel}, and \code{markov}. This is where the user has options for constructing the model. \code{Qhat} allows for transformaiton of streamflow observations. \code{QhatModel} offers a selection of distributions and serial-correlation at the annual and seasonal time-step. \code{markov} provides the type of markov model to run. The combination of these classes creates a 'hydroState' object or model that is evaluated. Vignettes \(\vignette{example1}, \vignette{example2}) show how to construct this \code{hydroState} object.

hydroState evaluates this model overtime at every time-step \eqn{t} where runoff (\eqn{_iQ_t}) is simulated as being in

For evaluating multiple models, a \code{hydroState.allModels} class exist. This class is a child of \code{hydroState} and depends on methods within this main parent class. Similarly, \code{hydroState.subAnnual.allModels} is a child of \code{hydroState.allModels} and provides the evaluation of multiple models at a seasonal timescale. Vignettes \(\vignette{example3}, \vignette{example4}) show how to use the \code{hydroState.allModels} and \code{hydroState.subAnnual.allModels} objects.


Below is a list of these classes and most important methods to build and evaluate hydroState models.

-----------------------------------------------------------------------------------
}



  \section{I. Transforming Observations}{

    \tabular{ll}{
      \code{\link{Qhat.none}}\tab No transformation\cr
      \code{\link{Qhat.log}}\tab Transformation in natural log-space\cr
      \code{\link{Qhat.burbidge}}\tab Burbige transformation\cr
      \code{\link{Qhat.boxcox}}\tab Boxcox transformation\cr
      ------------------------\tab------------------------------------------------
  }
  }

\section{Authors}{
Except where indicated otherwise, the methods and functions in this package were written by Tim Peterson.
}

\section{Acknowledgments}{
...

}


\keyword{ package }
\keyword{ rainfall-runoff }
\keyword{ states }
