##' @include abstracts.R
## @export
Qhat.burbidge <- setClass(
  # Set the name for the class
  "Qhat.burbidge",

  package='hydroState',

  contains=c('Qhat','Qhat.boxcox'),

  # Define the slots
  slots = c(
    input.data = "data.frame",
    parameters = "model_parameters",
    constant = 'numeric'
  ),

  # Set the default values for the slots. (optional)
  prototype=list(
    input.data = data.frame(year=c(0),month=c(0),precipitation=c(0)),
    parameters= new('model_parameters',c('lambda.burbidge'),c(10)),
    constant = 1
  )
)

# Valid object?
validObject <- function(object) {
  if(object@parameters$lambda.burbidge >=0) TRUE
  else warning("parameters$lambda.burbidge must be >=0")
}
setValidity("Qhat.burbidge", validObject)

setMethod("initialize","Qhat.burbidge", function(.Object, input.data, constant = 1) {
  .Object@input.data <- input.data
  .Object@constant <- constant
  validObject(.Object)
  .Object
}
)

# Calculate the transformed flow
setMethod(f="getQhat",signature=c("Qhat.burbidge",'data.frame'),definition=function(.Object, data)
{
  if (!is.data.frame(data))
    stop('"Data" must be a data.frame.')

  # Get object parameter list
  parameters <- parameters(.Object@parameters)
  data$q <- asinh((data$q + .Object@constant)/abs(parameters$lambda.burbidge))/asinh(1.0/abs(parameters$lambda.burbidge))
  data$p <- data$p

  return(data)
}
)

# Calculate the transformed flow using the object data
setMethod(f="getQhat",signature="Qhat.burbidge",definition=function(.Object)
{
  data = .Object@input.data$q
  return(getQhat(.Object, data))
}
)

setMethod(f="getQ.backTransformed",signature=c("Qhat.burbidge",'data.frame'),definition=function(.Object, data)
{
  if (!is.data.frame(data))
    stop('"Data" must be a data.frame.')

  parameters <- parameters(.Object@parameters)

  data$flow.modelled <- sinh(data$q * asinh(1.0/abs(parameters$lambda.burbidge)))*abs(parameters$lambda.burbidge) - .Object@constant
  return(data)

}
)

setMethod(f="get.zeroFlow",signature=c("Qhat.burbidge"),definition=function(.Object)
{
  data = data.frame(p = 0, q= 0)
  data = getQhat(.Object, data)
  return( as.numeric(data$q[1]))
}
)
