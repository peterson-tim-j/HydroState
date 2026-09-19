##' @include abstracts.R parameters.R
## @export
Qhat.log <- setClass(
  # Set the name for the class
  "Qhat.log",

  package='hydroState',

  contains=c('Qhat', 'Qhat.boxcox'),

  # Define the slots
  slots = c(
    input.data = "data.frame",
    parameters = 'model_parameters',
    constant = 'numeric'
  ),

  # Set the default values for the slots. (optional)
  prototype=list(
    input.data = data.frame(year=c(0),month=c(0),precipitation=c(0)),
    parameters= new('model_parameters',c(),c()),
    constant = 1
  )
)

# Initialise object
#setGeneric(name="initialize",def=function(.Object,input.data){standardGeneric("initialize")})
setMethod("initialize","Qhat.log", function(.Object, input.data, constant = 1) {
  .Object@input.data <- input.data
  .Object@constant <- constant
  validObject(.Object)
  .Object
}
)

# Calculate the transformed flow
setMethod(f="getQhat",signature=c("Qhat.log",'data.frame'),definition=function(.Object, data)
{

  if (!is.data.frame(data))
    stop('"data" must be a data.frame.')

  data$Qhat.flow <- log(data$q + .Object@constant)
  data$Qhat.precipitation <- data$precipitation
  return(data)
}
)

# Calculate the transformed flow using the object data
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

  data$flow.modelled <- exp(data$q) - .Object@constant
  return(data)
}
)

setMethod(f="get.zeroFlow",signature=c("Qhat.log"),definition=function(.Object)
{
  data = data.frame(p = 0, q= 0)
  data = getQhat(.Object, data)
  return( as.numeric(data$q[1]))
}
)
