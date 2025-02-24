##' @include abstracts.R QhatModel.subAnnual.homo.gamma.linear.R
## @export
QhatModel.subAnnual.homo.gamma.linear.AR1 <- setClass(
  # Set the name for the class
  "QhatModel.subAnnual.homo.gamma.linear.AR1",

  package='hydroState',

  contains=c('QhatModel.subAnnual.homo.gamma.linear'),

  # Set the default values for the slots. (optional)
  prototype=list(
    parameters = new('parameters',c('mean.a0.amp', 'mean.a0.phase', 'mean.a0.disp',
                                    'mean.a1','mean.AR1','std.a0'),c(1,1,1,1,1,1))

  )
)

# Valid object?
validObject <- function(object) {
  TRUE

}
setValidity("QhatModel.subAnnual.homo.gamma.linear.AR1", validObject)

# Initialise object
setMethod("initialize","QhatModel.subAnnual.homo.gamma.linear.AR1", function(.Object, input.data, transition.graph=matrix(T,2,2),
                                                                         state.dependent.mean.a0=T,state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                                         state.dependent.mean.AR1=F,
                                                                         subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F,subAnnual.dependent.std.a0=F) {
  .Object@input.data <- input.data
  .Object@use.truncated.dist <- F
  .Object@nStates = ncol(transition.graph)
  .Object@precip.delta = getStartEndIndex(input.data)


  # Check and set definition of seasons.
  .Object <- setSeasons(.Object, input.data)

  # Setup parameter definitions.
  .Object <- setupParameters(.Object, state.dependent.mean.a0, state.dependent.mean.a1, state.dependent.std.a0,
                             state.dependent.mean.AR1, state.dependent.mean.AR2=NA, state.dependent.mean.AR3=NA,
                             subAnnual.dependent.mean.a0, subAnnual.dependent.mean.a1, subAnnual.dependent.std.a0)

  validObject(.Object)
  .Object
}
)
#' @export getMean
setMethod(f="getMean",signature=c("QhatModel.subAnnual.homo.gamma.linear.AR1","data.frame"),definition=function(.Object, data)
{
  # Get object parameter list
  parameters = getParameters(.Object@parameters)

  # Get number of data rows
  nrows = length(data$Qhat.precipitation);

  # Set up AR parameters
  ncols.AR1 = length(parameters$mean.AR1)
  if (ncols.AR1==1 || ncols.AR1 ==.Object@nStates) {
    AR1.est = matrix(rep(parameters$mean.AR1,each=nrows),nrows,.Object@nStates);
  } else if (ncols.AR1<.Object@nStates) {
    stop(paste('The number of parameters for the AR1 term of the mean model must must equal 1 or the number of states of ',.Object@nStates))
  }

  # print('subAR1')
  # Get non-AR estimates.
  Qhat.model <- callNextMethod()

  # Now run AR for each continuous period
  for(j in 1:NROW(.Object@precip.delta)){ # could make for each in the future

    # Add AR components
    for (i in ((.Object@precip.delta[j,1]+1):.Object@precip.delta[j,2])) {

      Qhat.model[i,] <- Qhat.model[i,] + Qhat.model[i-1,] * AR1.est[i,]

    }


  }

  return(Qhat.model)
}
)
