##' @include abstracts.R QhatModel.subAnnual.homo.gamma.linear.R
##' @export
QhatModel.subAnnual.homo.gamma.linear.AR2 <- setClass(
  # Set the name for the class
  "QhatModel.subAnnual.homo.gamma.linear.AR2",

  package='hydroState',

  contains=c('QhatModel.subAnnual.homo.gamma.linear'),

  # Set the default values for the slots. (optional)
  prototype=list(
    parameters = new('parameters',c('mean.a0.amp', 'mean.a0.phase', 'mean.a0.disp',
                                    'mean.a1','mean.AR1','mean.AR2','std.a0'),c(1,1,1,1,1,1,1))

  )
)

# Valid object?
validObject <- function(object) {
  TRUE

}
setValidity("QhatModel.subAnnual.homo.gamma.linear.AR2", validObject)

# Initialise object
setMethod("initialize","QhatModel.subAnnual.homo.gamma.linear.AR2", function(.Object, input.data, transition.graph=matrix(T,2,2),
                                                                         state.dependent.mean.a0=T,state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                                         state.dependent.mean.AR1=F, state.dependent.mean.AR2=F,
                                                                         subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F,subAnnual.dependent.std.a0=F) {
  .Object@use.truncated.dist <- F
  .Object@nStates = ncol(transition.graph)

  # Check and set definition of seasons.
  .Object <- setSeasons(.Object, input.data)

  # Setup parameter definitions.
  .Object <- setupParameters(.Object, state.dependent.mean.a0, state.dependent.mean.a1, state.dependent.std.a0,
                             state.dependent.mean.AR1, state.dependent.mean.AR2, state.dependent.mean.AR3=NA,
                             subAnnual.dependent.mean.a0, subAnnual.dependent.mean.a1, subAnnual.dependent.std.a0)

  validObject(.Object)
  .Object
}
)

setMethod(f="getMean",signature=c("QhatModel.subAnnual.homo.gamma.linear.AR2","data.frame"),definition=function(.Object, data)
{
  # Get object parameter list
  parameters = getParameters(.Object@parameters)

  # Get number of data rows
  nrows = length(data$Qhat.precipitation);

  # Set up AR parameters
  ncols.AR1 = length(parameters$mean.AR1)
  ncols.AR2 = length(parameters$mean.AR2)
  if (ncols.AR1==1 || ncols.AR1 ==.Object@nStates) {
    AR1.est = matrix(rep(parameters$mean.AR1,each=nrows),nrows,.Object@nStates);
  } else if (mean.AR1<.Object@nStates) {
    stop(paste('The number of parameters for the a0 term of the mean model must must equal 1 or the number of states of ',.Object@nStates))
  }
  if (ncols.AR2==1 || ncols.AR2 ==.Object@nStates) {
    AR2.est = matrix(rep(parameters$mean.AR2,each=nrows),nrows,.Object@nStates);
  } else if (mean.AR2<.Object@nStates) {
    stop(paste('The number of parameters for the AR2 term of the mean model must must equal 1 or the number of states of ',.Object@nStates))
  }

  # Get non-AR estimates.
  Qhat.model <- callNextMethod()

  # Add AR components
  Qhat.model[2,] <-Qhat.model[2,] + Qhat.model[1,] * AR1.est[2,]
  for (i in 3:nrows) {
    Qhat.model[i,] <-Qhat.model[i,] + Qhat.model[i-1,] * AR1.est[i,] + Qhat.model[i-2,] * AR2.est[i,]
  }

  return(Qhat.model)
}
)
