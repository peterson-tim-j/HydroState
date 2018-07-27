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
    parameters = "list"
  ),

  # Set the default values for the slots. (optional)
  prototype=list(
    input.data = data.frame(year=c(0),month=c(0),precipitation=c(0)),
    parameters= list(lambda = 0.5)
  )
)

# Valid object?
validObject <- function(object) {
  if(object@parameters$lambda >=0) TRUE
  else warning("parameters$lambda must be >=0")
}
setValidity("Qhat.burbidge", validObject)

# # Initialise object
# #setGeneric(name="initialize",def=function(.Object,input.data){standardGeneric("initialize")})
# setMethod("initialize","Qhat.burbidge", function(.Object, input.data) {
#   .Object@input.data <- input.data
#   validObject(.Object)
#   .Object
# }
# )
#
# # create a method to assign the parameter
# #setGeneric(name="setParameters",def=function(.Object,parameters){standardGeneric("setParameters")})
# setMethod(f="setParameters",
#           signature="Qhat.burbidge",
#           definition=function(.Object,parameters)
#           {
#             .Object@parameters$lambda <- parameters$lambda
#             validObject(.Object)
#             return(.Object)
#           }
# )
setMethod(f="setParameters.fromTransformed",
          signature="Qhat.burbidge",
          definition=function(.Object,parameters)
          {
            .Object@parameters$lambda <- 10^parameters$lambda
            validObject(.Object)
            return(.Object)          }
)

# # create a method to assign the parameter
# #setGeneric(name="getParameters",def=function(.Object){standardGeneric("getParameters")})
# setMethod(f="getParameters",signature="Qhat.burbidge",definition=function(.Object)
#           {
#             return(list(lambda = .Object@parameters$lambda))
#           }
# )
setMethod(f="getTransformedParameterBounds",
          signature="Qhat.burbidge",
          definition=function(.Object)
          {
            parameters = getParameters(.Object)
            lowerBound = list(lambda = rep(-6, length(.Object@parameters$lambda)))
            upperBound = list(lambda = rep(3, length(.Object@parameters$lambda)))
            return(list(lower = lowerBound, upper = upperBound))
          }
)

# Calculate the transformed flow
#setGeneric(name="getQhat",def=function(.Object, data){standardGeneric("getQhat")})
setMethod(f="getQhat",signature=c("Qhat.burbidge",'data.frame'),definition=function(.Object, data)
          {
            if (is.data.frame(data))
              data = data$flow

            return(asinh(data/abs(Object@parameters$lambda))/abs(Object@parameters$lambda))

          }
)
#
# # Calculate the transformed flow using the object data
# #setGeneric(name="getQhat",def=function(.Object){standardGeneric("getQhat")})
# setMethod(f="getQhat",signature="Qhat.burbidge",definition=function(.Object)
#           {
#              data = .Object@input.data$flow
#              return(getQhat(.Object, data))
#           }
# )
