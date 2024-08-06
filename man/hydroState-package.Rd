\name{hydroState-package}

\alias{hydroState-package}

\alias{hydroState}

\docType{package}

\title{Overview of methods and procedures}

\description{

\code{hydroState} provides methods to construct and evaluate hidden Markov models (HMM) of annual or subannual streamflow with precipitation as a predictor. The state of the relationship between these observations is evaluated overtime. The package contains a default \code{hydroState} model to evaluate probable shifts in this rainfall-runoff relationship or the model can be constructed with various items to better define the relationship as discussed below. Once the default or constructed model is built (\code{buildModel}), the model is fitted (\code{fitModel}) to determine the most likly rainfall-runoff state at each time-step. To assess the adequacy of the fit, the residuals are plotted (\code{plot.residuals}), and an adequate fit requires the residuals to be normaly distributed and not have any correlation or trends. The resulting runoff states from the fitted model can then be evaluated overtime (\code{plot.states}) and even exported (\code{get.states}) with the state values, confidence intervals, and conditional probabilities of each state.

A \code{hydroState} model operates in S4, object oriented programming, and constructing a model requires the selection of three objects. The first object is the \code{data.transform}, or 'Qhat' object, from \code{select.transform} that offers various transformaitons of the observations including: 'boxcox', 'log', 'burbidge', or 'none'. The second object is the \code{state.model}, or 'QhatModel' object, from \code{select.stateModel} that offers a variety of items to better define the rainfall-runoff relationship including: the state number (1, 2, or 3), residual distribution ('normal', 'gamma', or 'truc.normal'), and degree of auto-correlation ('AR1', 'AR2', or 'AR3'). Furthermore \code{select.stateModel} provides an option to select the parameters in the model expected to shift (\code{state.shift.parameters}) and apply seasonal variation within model parameters (\code{seasonal.paraneters}) if monthly data is provided. The third object is the form of the \code{Markov}, or 'markov' object, from \code{select.Markov} which currently only provides a homogeneous Markov model. Further details on the selection of options within these objects for constructing a \code{hydroState} model are explained within each topic.

There is an additional option to construct all possible types of models using the \code{buildModelAll}, and compare them using the same \code{fitModel} function. The most likely model is selected based on the AIC where the best model will have the lowest AIC. To begin, it is recommended to evaluate the default model at first with one state and again with two states. This documentation is organized with four sections that define the general workflow for using \code{hydroState}: Build - Fit - Review - Evaluate.


-----------------------------------------------------------------------------------
}



  \section{I. Build}{

    \tabular{ll}{
      \code{\link{buildModel}}\tab build hydroState model\cr
      \code{\link{select.transform}}\tab select transformation of observations\cr
      \code{\link{select.stateModel}}\tab select type of state model\cr
      \code{\link{select.Markov}}\tab select form of Markov\cr
      \code{\link{buildModelAll}}\tab build all possible models\cr
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
Except where indicated otherwise, the methods and functions in this package were written by Tim Peterson.
}

\section{Acknowledgments}{
...

}


\keyword{ package }
\keyword{ rainfall-runoff }
\keyword{ states }
