##' @include abstracts.R parameters.R

## @export
Qhat.boxcox <- setClass(
  # Set the name for the class
  "Qhat.boxcox",

  package='hydroState',

  contains=c('Qhat'),

  # Define the slots
  slots = c(
    input.data = "data.frame",
    parameters = "model_parameters",
    constant = 'numeric'
  ),

  # Set the default values for the slots. (optional)
  prototype=list(
    input.data = data.frame(year=c(0),month=c(0),precipitation=c(0)),
    parameters = new('model_parameters',c('lambda'),c(1)),
    constant = 1
  )
)

# Valid object?
validObject <- function(object) {
  if(object@parameters$lambda >=0) TRUE
  else warning("parameters$lambda must be >=0")
}
setValidity("Qhat.boxcox", validObject)

# Initialise object
setMethod("initialize","Qhat.boxcox", function(.Object, input.data, constant = 1) {
  .Object@input.data <- input.data
  .Object@constant <- constant
  validObject(.Object)
  .Object
}
)
# Calculate the transformed flow
setMethod(f="getQhat",signature=c("Qhat.boxcox",'data.frame'),definition=function(.Object, data)
          {
            if (!is.data.frame(data))
              stop('"Data" must be a data.frame.')

            # Get object parameter list
            parameters <- parameters(.Object@parameters)

            if (parameters$lambda>1e-8) {
              data$q <- ((data$q + .Object@constant)^parameters$lambda - 1)/parameters$lambda
            } else {
              data$q <- log(data$q + .Object@constant)
            }
            data$p <- data$p

            return(data)
          }
)

# Calculate the transformed flow using the object data
setMethod(f="getQhat",signature="Qhat.boxcox",definition=function(.Object)
          {
             data = .Object@input.data
             return(getQhat(.Object, data))
          }
)

setMethod(f="getQ.backTransformed",signature=c("Qhat.boxcox",'data.frame'),definition=function(.Object, data)
{
  if (!is.data.frame(data))
    stop('"Data" must be a data.frame.')

  parameters <- parameters(.Object@parameters)

  if (parameters$lambda>1e-8) {
    data$flow.modelled <- ( data$q * parameters$lambda + 1) ^ (1/parameters$lambda) - .Object@constant
  } else {
    data$flow.modelled <- exp(data$q)-1
  }
  return(data)
}
)

setMethod(f="get.zeroFlow",signature=c("Qhat.boxcox"),definition=function(.Object)
{
  data = data.frame(p = 0, q= 0)
  data = getQhat(.Object, data)
  return( as.numeric(data$q[1]))
}
)
