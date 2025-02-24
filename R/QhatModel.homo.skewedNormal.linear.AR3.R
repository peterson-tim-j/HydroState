##' @include abstracts.R QhatModel.homo.skewedNormal.linear.R QhatModel.homo.normal.linear.AR3.R
## @export
QhatModel.homo.skewedNormal.linear.AR3 <- setClass(
  # Set the name for the class
  "QhatModel.homo.skewedNormal.linear.AR3",

  package='hydroState',

  contains=c('QhatModel.homo.skewedNormal.linear','QhatModel.homo.normal.linear.AR3'),

  # Set the default values for the slots. (optional)
  prototype=list(
    parameters =  new('parameters',c('mean.a0', 'mean.a1','mean.AR1','mean.AR2','mean.AR3','std.a0', 'shape.a0'),c(1,1,1,1,1,1,1))
  )
)

# Valid object?
validObject <- function(object) {
  TRUE

}
setValidity("QhatModel.homo.skewedNormal.linear.AR3", validObject)

# Initialise object
#setGeneric(name="initialize",def=function(.Object,input.data){standardGeneric("initialize")})
setMethod("initialize","QhatModel.homo.skewedNormal.linear.AR3", function(.Object, input.data,transition.graph=matrix(T,2,2),state.dependent.mean.a0=T,
                                                                          state.dependent.mean.a1=F,state.dependent.mean.AR1=F,
                                                                          state.dependent.mean.AR2=F, state.dependent.mean.AR3=F,
                                                                          state.dependent.std.a0=T, state.dependent.shape.a0=T) {
  .Object@input.data <- input.data
  .Object@use.truncated.dist <- F
  .Object@nStates = ncol(transition.graph)

  # Set the number of parameter values per parameter name.
  parameter.length <- as.numeric(c(state.dependent.mean.a0, state.dependent.mean.a1, state.dependent.mean.AR1,
                                   state.dependent.mean.AR2, state.dependent.mean.AR3, state.dependent.std.a0,
                                   state.dependent.shape.a0)) * (.Object@nStates-1) + 1

  # Set up model terms for mean and standard deviation.
  .Object@parameters = new('parameters', c('mean.a0', 'mean.a1', 'mean.AR1', 'mean.AR2', 'mean.AR3','std.a0','shape.a0'), parameter.length)


  validObject(.Object)
  .Object
}
)

setMethod(f="getMean",signature=c("QhatModel.homo.skewedNormal.linear.AR3","data.frame"),definition=function(.Object, data)
{
  Qhat.model.NAs = matrix(NA,NROW(data),.Object@nStates)

  for(i in 1:NROW(.Object@precip.delta)){
    Qhat.model.NAs[.Object@precip.delta[i,1]:.Object@precip.delta[i,2],] = getMean.AR3(.Object, data[.Object@precip.delta[i,1]:.Object@precip.delta[i,2],])
  }

  return(Qhat.model.NAs)
}

)
