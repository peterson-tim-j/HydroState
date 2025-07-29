##' @include abstracts.R parameters.R
## @export
Qhat.log <- setClass(
  # Set the name for the class
  "Qhat.log",

  package='hydroState',

  contains=c('Qhat'),

  # Define the slots
  slots = c(
    input.data = "data.frame",
    parameters = 'parameters'
  ),

  # Set the default values for the slots. (optional)
  prototype=list(
    input.data = data.frame(year=c(0),month=c(0),precipitation=c(0)),
    parameters= new('parameters',c(),c())
  )


)

# Initialise object
#setGeneric(name="initialize",def=function(.Object,input.data){standardGeneric("initialize")})
setMethod("initialize","Qhat.log", function(.Object, input.data) {
  .Object@input.data <- input.data
  validObject(.Object)
  .Object
}
)

# Calculate the transformed flow
#setGeneric(name="getQhat",def=function(.Object, data){standardGeneric("getQhat")})
setMethod(f="getQhat",signature=c("Qhat.log",'data.frame'),definition=function(.Object, data)
{

  if (!is.data.frame(data))
    stop('"data" must be a data.frame.')

  data$Qhat.flow <- log(data$flow +1)
  data$Qhat.precipitation <- data$precipitation
  return(data)
}
)

# Calculate the transformed flow using the object data
#setGeneric(name="getQhat",def=function(.Object){standardGeneric("getQhat")})
setMethod(f="getQhat",signature="Qhat.log",definition=function(.Object)
{
   data = .Object@input.data
   return(getQhat(.Object, data))
}
)


setMethod(f="getQ.backTransformed",signature=c("Qhat.log",'data.frame'),definition=function(.Object, data)
{

  if (!is.data.frame(data))
    stop('"data" must be a data.frame.')

  data$flow.modelled <- exp(data$Qhat.flow)-1
  return(data)
}
)
