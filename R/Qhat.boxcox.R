##' @include abstracts.R parameters.R
##' @export
Qhat.boxcox <- setClass(
  # Set the name for the class
  "Qhat.boxcox",

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
    parameters= new('parameters',c('lambda'),c(1))
  )
)

# Valid object?
validObject <- function(object) {
  if(object@parameters$lambda >=0) TRUE
  else warning("parameters$lambda must be >=0")
}
setValidity("Qhat.boxcox", validObject)

# Initialise object
#setGeneric(name="initialize",def=function(.Object,input.data){standardGeneric("initialize")})
setMethod("initialize","Qhat.boxcox", function(.Object, input.data) {
  .Object@input.data <- input.data
  validObject(.Object)
  .Object
}
)
# Calculate the transformed flow
#setGeneric(name="getQhat",def=function(.Object, data){standardGeneric("getQhat")})
setMethod(f="getQhat",signature=c("Qhat.boxcox",'data.frame'),definition=function(.Object, data)
          {
            if (!is.data.frame(data))
              stop('"Data" must be a data.frame.')


            # Get object parameter list
            parameters = getParameters(.Object@parameters)

            if (parameters$lambda>1e-8) {
              data$Qhat.flow <- ((data$flow+1)^parameters$lambda-1)/parameters$lambda
            } else {
              data$Qhat.flow <- log(data$flow +1)
            }
            data$Qhat.precipitation <- data$precipitation

            return(data)
          }
)

# Calculate the transformed flow using the object data
#setGeneric(name="getQhat",def=function(.Object){standardGeneric("getQhat")})
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

  # Get object parameter list
  parameters = getParameters(.Object@parameters)

  if (parameters$lambda>1e-8) {
    data$flow.modelled <- ( data$Qhat.flow * parameters$lambda + 1) ^ (1/parameters$lambda)-1
  } else {
    data$flow.modelled <- exp(data$Qhat.flow)-1
  }
  return(data)
}
)
