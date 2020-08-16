##' @include abstracts.R
##' @export
Qhat.burbidge <- setClass(
  # Set the name for the class
  "Qhat.burbidge",

  package='hydroState',

  contains=c('Qhat'),

  # Define the slots
  slots = c(
    input.data = "data.frame",
    parameters = "parameters"
  ),

  # Set the default values for the slots. (optional)
  prototype=list(
    input.data = data.frame(year=c(0),month=c(0),precipitation=c(0)),
    parameters= new('parameters',c('lambda.burbidge'),c(10))
  )
)

# Valid object?
validObject <- function(object) {
  if(object@parameters$lambda.burbidge >=0) TRUE
  else warning("parameters$lambda.burbidge must be >=0")
}
setValidity("Qhat.burbidge", validObject)

setMethod("initialize","Qhat.burbidge", function(.Object, input.data) {
  .Object@input.data <- input.data
  validObject(.Object)
  .Object
}
)

# Calculate the transformed flow
#setGeneric(name="getQhat",def=function(.Object, data){standardGeneric("getQhat")})
setMethod(f="getQhat",signature=c("Qhat.burbidge",'data.frame'),definition=function(.Object, data)
{
  if (!is.data.frame(data))
    stop('"Data" must be a data.frame.')


  # Get object parameter list
  parameters = getParameters(.Object@parameters)

  data$Qhat.flow <- asinh(data$flow/abs(parameters$lambda.burbidge))/asinh(1.0/abs(parameters$lambda.burbidge))
  data$Qhat.precipitation <- data$precipitation

  return(data)

}
)

# Calculate the transformed flow using the object data
#setGeneric(name="getQhat",def=function(.Object){standardGeneric("getQhat")})
setMethod(f="getQhat",signature="Qhat.burbidge",definition=function(.Object)
{
  data = .Object@input.data$flow
  return(getQhat(.Object, data))
}
)

setMethod(f="getQ.backTransformed",signature=c("Qhat.burbidge",'data.frame'),definition=function(.Object, data)
{
  if (!is.data.frame(data))
    stop('"Data" must be a data.frame.')

  # Get object parameter list
  parameters = getParameters(.Object@parameters)

  data$flow.modelled <- sinh(data$Qhat.flow * asinh(1.0/abs(parameters$lambda.burbidge)))*abs(parameters$lambda.burbidge)
  return(data)

}
)
