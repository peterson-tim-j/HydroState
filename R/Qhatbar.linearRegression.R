##' @include abstracts.R
##' @export
Qhatbar.linearRegression <- setClass(
  # Set the name for the class
  "Qhatbar.linearRegression",

  package='hydroState',

  contains=c('Qhatbar'),

  # Define the slots
  slots = c(
    input.data = "data.frame",
    Pbar = 'numeric',
    parameters = "list"
  ),

  # Set the default values for the slots. (optional)
  prototype=list(
    input.data = data.frame(year=c(0),month=c(0),precipitation=c(0)),
    Pbar=Inf,
    parameters = list(
      slope = 1.0,
      intercept = 1.0)
  )
)

# Valid object?
validObject <- function(object) {
  if(object@Pbar >=0) TRUE
  else warning("Pbar must be >=0")
}
setValidity("Qhatbar.linearRegression", validObject)

# Initialise object
#setGeneric(name="initialize",def=function(.Object,input.data){standardGeneric("initialize")})
setMethod("initialize","Qhatbar.linearRegression", function(.Object, input.data) {
  .Object@input.data <- input.data
  .Object@Pbar <- mean(input.data$precipitation)
  validObject(.Object)
  .Object
}
)


# create a method to assign the parameter
#setGeneric(name="setParameters",def=function(.Object,parameters) {standardGeneric("setParameters")})
setMethod(f="setParameters",
          signature="Qhatbar.linearRegression",
          definition=function(.Object,parameters)
          {
            .Object@parameters$slope <- parameters$slope
            .Object@parameters$intercept <- parameters$intercept
            validObject(.Object)
            return(.Object)
          }
)
setMethod(f="setParameters.fromTransformed",
          signature="Qhatbar.linearRegression",
          definition=function(.Object,parameters)
          {
            .Object@parameters$slope <- parameters$slope
            .Object@parameters$intercept <- parameters$intercept
            return(.Object)
          }
)
# create a method to assign the parameter
#setGeneric(name="getParameters",def=function(.Object,position) {standardGeneric("getParameters")})
setMethod(f="getParameters",
          signature="Qhatbar.linearRegression",
          definition=function(.Object)
          {
            return(list(slope=.Object@parameters$slope,intercept=.Object@parameters$intercept))
          }
)
setMethod(f="getTransformedParameterBounds",
          signature="Qhatbar.linearRegression",
          definition=function(.Object)
          {
            parameters = getParameters(.Object)
            lowerBound = list(slope = rep(0.05, length(.Object@parameters$slope)),
                              intercept = rep(2, length(.Object@parameters$intercept))
            )
            upperBound = list(slope = rep(0.05, length(.Object@parameters$slope)),
                              intercept = rep(2, length(.Object@parameters$intercept))
            )
            return(list(lower = lowerBound, upper = upperBound))
          }
)


# Calculate the transformed flow at the mean annual precip
#setGeneric(name="getQhatbar",def=function(.Object, Qhat) {standardGeneric("getQhatbar")})
setMethod(f="getQhatbar",signature=c("Qhatbar.linearRegression","numeric"),definition=function(.Object, data)
          {

            # Calculate the residuals between the input transformed flow and the regression estimate.
            resid = data - getQhatbarModelled(.Object, data)

            # Get Qhat at the man ennual precip.
            Qbar = getQhatbarModelled.atPbar(.Object, data)
            # Calculate the streamflow at the mean annual precip.
            return(data.frame(Qhatbar=Qbar+resid))
          }
)

# Calculate the regerssion estimate of the transformed flowfor the object annual precip
#' @exportMethod getQhatbarModelled
setGeneric(name="getQhatbarModelled",def=function(.Object,data) {standardGeneric("getQhatbarModelled")})
setMethod(f="getQhatbarModelled",signature="Qhatbar.linearRegression",definition=function(.Object, data)
          {
            df = cbind.data.frame(.Object@input.data, Qhat=data)
            lm.fit = lm(Qhat ~ precipitation, df);
            lm.intercept = lm.fit$coefficients[[1]]
            lm.slope = lm.fit$coefficients[[2]]

            # Calculate the residuals between the input transformed flow and the regression estimate.
            return(.Object@input.data$precipitation * lm.slope * .Object@parameters$slope + lm.intercept * .Object@parameters$intercept)
          }
)
setGeneric(name="getQhatbarModelled.atPbar",def=function(.Object,data) {standardGeneric("getQhatbarModelled.atPbar")})
setMethod(f="getQhatbarModelled.atPbar",signature=c("Qhatbar.linearRegression","numeric"),definition=function(.Object, data)
        {
          df = cbind.data.frame(.Object@input.data, Qhat=data)
          lm.fit = lm(Qhat ~ precipitation, df);
          lm.intercept = lm.fit$coefficients[[1]]
          lm.slope = lm.fit$coefficients[[2]]

          # Calculate the residuals between the input transformed flow and the regression estimate.
          return(.Object@Pbar * lm.slope * .Object@parameters$slope + lm.intercept * .Object@parameters$intercept)
        }
)

