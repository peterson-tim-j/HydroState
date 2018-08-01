##' @include abstracts.R QhatModel.homo.gamma.linear.R QhatModel.homo.normal.linear.AR2.R
##' @export
QhatModel.homo.gamma.linear.AR2 <- setClass(
  # Set the name for the class
  "QhatModel.homo.gamma.linear.AR2",

  package='hydroState',

  contains=c('QhatModel.homo.gamma.linear', 'QhatModel.homo.normal.linear.AR2'),

  # Set the default values for the slots. (optional)
  prototype=list(
    parameters =  new('parameters',c('mean.a0', 'mean.a1','mean.AR1','mean.AR2','std.a0'),c(1,1,1,1,1))
  )
)

# Valid object?
validObject <- function(object) {
  TRUE

}
setValidity("QhatModel.homo.gamma.linear.AR2", validObject)

# Initialise object
#setGeneric(name="initialize",def=function(.Object,input.data){standardGeneric("initialize")})
setMethod("initialize","QhatModel.homo.gamma.linear.AR2", function(.Object,input.data, transition.graph=matrix(T,2,2),state.dependent.mean.a0=T,
                                                                  state.dependent.mean.a1=F,state.dependent.mean.AR1=F,state.dependent.mean.AR2=F,
                                                                  state.dependent.std.a0=T) {

  .Object@use.truncated.dist <- F
  .Object@nStates = ncol(transition.graph)

  # Set the number of parameter values per parameter name.
  parameter.length <- as.numeric(c(state.dependent.mean.a0, state.dependent.mean.a1, state.dependent.mean.AR1, state.dependent.mean.AR2,
                                   state.dependent.std.a0)) * (.Object@nStates-1) + 1

  # Set up model terms for mean and standard deviation.
  .Object@parameters = new('parameters', c('mean.a0', 'mean.a1', 'mean.AR1', 'mean.AR2', 'std.a0'), parameter.length)

  validObject(.Object)
  .Object
}
)

setMethod(f="getMean",signature=c("QhatModel.homo.gamma.linear.AR2","data.frame"),definition=function(.Object, data)
{
  return(getMean.AR2(.Object, data))
}
)
