##' @include abstracts.R
##' @export
Qhat.log <- setClass(
  # Set the name for the class
  "Qhat.log",

  package='hydroState',

  contains=c('Qhat'),

  # Define the slots
  slots = c(
    input.data = "data.frame",
    parameters = "list"
  ),

  # Set the default values for the slots. (optional)
  prototype=list(
    input.data = data.frame(year=c(0),month=c(0),precipitation=c(0)),
    parameters= list()
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

# create a method to assign the parameter
#setGeneric(name="setParameters",def=function(.Object,parameters){standardGeneric("setParameters")})
setMethod(f="setParameters",
          signature="Qhat.log",
          definition=function(.Object,parameters)
          {
            return(.Object)
          }
)
setMethod(f="setParameters.fromTransformed",
          signature="Qhat.log",
          definition=function(.Object,parameters)
          {
            return(.Object)          }
)

# create a method to assign the parameter
#setGeneric(name="getParameters",def=function(.Object){standardGeneric("getParameters")})
setMethod(f="getParameters",signature="Qhat.log",definition=function(.Object)
{
  return(list())
}
)


setMethod(f="getTransformedParameterBounds",
          signature="Qhat.log",
          definition=function(.Object)
          {
            lowerBound = list()
            upperBound = list()
            return(list(lower = lowerBound, upper = upperBound))
          }
)

# Calculate the transformed flow
#setGeneric(name="getQhat",def=function(.Object, data){standardGeneric("getQhat")})
setMethod(f="getQhat",signature=c("Qhat.log",'data.frame'),definition=function(.Object, data)
          {

            if (is.data.frame(data))
              data = data$flow

              return( log(data + 1))
          }
)

# Calculate the transformed flow using the object data
#setGeneric(name="getQhat",def=function(.Object){standardGeneric("getQhat")})
setMethod(f="getQhat",signature="Qhat.log",definition=function(.Object)
          {
             data = .Object@input.data$flow
             return(getQhat(.Object, data))
          }
)
