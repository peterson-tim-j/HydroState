##' @include abstracts.R QhatModel.homo.normal.linear.R
##' @export
QhatModel.homo.normal.linear.AR1 <- setClass(
  # Set the name for the class
  "QhatModel.homo.normal.linear.AR1",

  package='hydroState',

  contains=c('QhatModel.homo.normal.linear'),

  # Set the default values for the slots. (optional)
  prototype=list(
    parameters = list(
      mean.a1 = 1.0,
      mean.a0 = 1.0,
      mean.correlation = 0,
      std.a0 = 1.0
      )
  )
)

# Valid object?
validObject <- function(object) {
  TRUE

}
setValidity("QhatModel.homo.normal.linear.AR1", validObject)

# Initialise object
#setGeneric(name="initialize",def=function(.Object,input.data){standardGeneric("initialize")})
setMethod("initialize","QhatModel.homo.normal.linear.AR1", function(.Object, input.data, transition.graph=matrix(T,2,2),state.dependent.mean.a0=T, state.dependent.mean.a1=F,state.dependent.mean.correlation=F, state.dependent.std.a0=T) {
  .Object@input.data <- input.data

  .Object@nStates = ncol(transition.graph)

  # Set up model terms for mean and standard deviation.
  .Object@parameters$mean.a0 = 0;
  .Object@parameters$mean.a1 = 0.1;
  .Object@parameters$std.a0 = 100;
  .Object@parameters$mean.correlation = 0.1;
  if (state.dependent.mean.a0)
    .Object@parameters$mean.a0 = rep(.Object@parameters$mean.a0, .Object@nStates);
  if (state.dependent.mean.a1)
    .Object@parameters$mean.a1 = rep(.Object@parameters$mean.a1, .Object@nStates);
  if (state.dependent.mean.correlation)
    .Object@parameters$mean.correlation = rep(.Object@parameters$mean.correlation, .Object@nStates);
  if (state.dependent.std.a0)
    .Object@parameters$std.a0 = rep(.Object@parameters$std.a0, .Object@nStates);

  validObject(.Object)
  .Object
}
)


# create a method to assign the parameter
#setGeneric(name="setParameters",def=function(.Object,parameters) {standardGeneric("setParameters")})
setMethod(f="setParameters",
          signature="QhatModel.homo.normal.linear.AR1",
          definition=function(.Object,parameters)
          {
            .Object@parameters$mean.a1 <- parameters$mean.a1
            .Object@parameters$mean.a0 <- parameters$mean.a0
            .Object@parameters$std.a0 <- parameters$std.a0
            .Object@parameters$mean.correlation <- parameters$mean.correlation;
            validObject(.Object)
            return(.Object)
          }
)
setMethod(f="setParameters.fromTransformed",
          signature="QhatModel.homo.normal.linear.AR1",
          definition=function(.Object,parameters)
          {
            .Object@parameters$mean.a1 <- 10^parameters$mean.a1
            .Object@parameters$mean.a0 <- parameters$mean.a0*100
            .Object@parameters$std.a0 <- 10^parameters$std.a0
            .Object@parameters$mean.correlation <- parameters$mean.correlation;
            return(.Object)
          }
)
# create a method to assign the parameter
#setGeneric(name="getParameters",def=function(.Object,position) {standardGeneric("getParameters")})
setMethod(f="getParameters",
          signature="QhatModel.homo.normal.linear.AR1",
          definition=function(.Object)
          {
            return(list(mean.a1=.Object@parameters$mean.a1,mean.a0=.Object@parameters$mean.a0, mean.correlation=.Object@parameters$mean.correlation, std.a0=.Object@parameters$std.a0))
          }
)

setMethod(f="getTransformedParameterBounds",
          signature="QhatModel.homo.normal.linear.AR1",
          definition=function(.Object)
          {
            parameters = getParameters(.Object)
            lowerBound = list(mean.a1 = rep(-6, length(.Object@parameters$mean.a1)),
                              mean.a0 = rep(-3, length(.Object@parameters$mean.a0)),
                              mean.correlation = rep(-1, length(.Object@parameters$mean.correlation)),
                              std.a0 = rep(log10(0.01), length(.Object@parameters$std.a0))
            )
            upperBound = list(mean.a1 = rep(log10(0.1), length(.Object@parameters$mean.a1)),
                              mean.a0 = rep(3, length(.Object@parameters$mean.a0)),
                              mean.correlation = rep(1, length(.Object@parameters$mean.correlation)),
                              std.a0 = rep(log10(1000), length(.Object@parameters$std.a0))
            )
            return(list(lower = lowerBound, upper = upperBound))
          }
)


# Calculate the transformed flow at the mean annual precip
#setGeneric(name="getMean",def=function(.Object, data) {standardGeneric("getMean")})
setMethod(f="getMean",signature=c("QhatModel.homo.normal.linear.AR1","data.frame"),definition=function(.Object, data)
          {

            ncols.a1 = length(.Object@parameters$mean.a1)
            ncols.a0 = length(.Object@parameters$mean.a0)
            ncols.correlation = length(.Object@parameters$mean.correlation)
            nrows = length(data$precipitation);
            ncols.max = max(c(ncols.a0 ,ncols.a1,ncols.correlation))
            if (ncols.max > .Object@nStates)
                stop(paste('The number of parameters for each term of the mean model must must equal 1 or the number of states of ',.Object@nStates))

            # Check which terms are uniform for all states and whic terms are unique
            # to each state.
            if (ncols.a0==1 || ncols.a0==.Object@nStates) {
                a0.est = matrix(rep(.Object@parameters$mean.a0,each=nrows),nrows,.Object@nStates);
            } else if (ncols.a0<.Object@nStates) {
              stop(paste('The number of parameters for the a0 term of the mean model must must equal 1 or the number of states of ',.Object@nStates))
            }
            if (ncols.a1==1 || ncols.a1==.Object@nStates) {
              a1.est = matrix(rep(.Object@parameters$mean.a1,each=nrows),nrows,.Object@nStates);
            } else if (ncols.a1<.Object@nStates) {
              stop(paste('The number of parameters for the a1 term of the mean model must must equal 1 or the number of states of ',.Object@nStates))
            }
            if (ncols.correlation==1 || ncols.correlation ==.Object@nStates) {
              correlation.est = matrix(rep(.Object@parameters$mean.correlation,each=nrows),nrows,.Object@nStates);
            } else if (mean.correlation<.Object@nStates) {
              stop(paste('The number of parameters for the a0 term of the mean model must must equal 1 or the number of states of ',.Object@nStates))
            }

            precip.data = matrix(data$precipitation,nrows,.Object@nStates);

            # Calculate the mean with the AR1 componants
            Qhat.model <- matrix(0, nrows, .Object@nStates)
            Qhat.model[1,] <- precip.data[1,] * a1.est[1,] + a0.est[1,]
            for (i in 2:nrows) {
              Qhat.model[i,] <- precip.data[i,] * a1.est[i,] + a0.est[i,] + Qhat.model[i-1,] * correlation.est[i,]
            }

            return(Qhat.model)
          }
)

# setMethod(f="getPseudoResiduals",signature="getPseudoResiduals",definition=function(.Object, data)
# {
#
#   corr = .Object@parameters$mean.correlation
#   if (length(corr)>1)
#     stop("Pseudo residuals cannot yet be calaculated for state dependent correlation structures.")
#
#   # This model DOES HAVE serial correlation, so  apply method from 6.3.2 of
#   # Zucchini, McDonald and Langrock, 2016, Hidden Markov Models for Time Series, 2nd Edision
#   #---------------------------------------------------------------------
#   nrows = length(data$precipitation);
#   uniform.pseudo.residuals = rep(0, nrows-2)
#   normal.pseudo.residuals = rep(0, nrows-2)
#   for (i in 2:(nrows-1)) {
#     normal.pseudo.residuals[i-1] <- (data$Qhat[i] - corr/(1+corr^2) * (data$Qhat[i-1] + data$Qhat[i+1]))*sqrt(1+corr^2)
#
#     uniform.pseudo.residuals[i-1] <- dnorm(normal.pseudo.residuals[i-1],mean = corr*(data$Qhat[i-1] + data$Qhat[i+1])/(1+corr^2), sd = sqrt(1/(1+corr^2)))
#   }
#   return(list(uniform=uniform.pseudo.residuals, normal=normal.pseudo.residuals))
# }
# )
