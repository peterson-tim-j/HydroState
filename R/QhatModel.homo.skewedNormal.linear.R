##' @include abstracts.R QhatModel.homo.skewedNormal.linear.R
##' @export
QhatModel.homo.skewedNormal.linear <- setClass(
  # Set the name for the class
  "QhatModel.homo.skewedNormal.linear",

  package='hydroState',

  contains=c('QhatModel.homo.normal.linear'),

    # Set the default values for the slots. (optional)
  prototype=list(
    input.data = data.frame(year=c(0),month=c(0),precipitation=c(0)),
    nStates = Inf,
    parameters = list(
      location.a0 = 1.0,
      location.a1 = 1.0,
      scale.a0 = 1.0,
      shape.a0 = 0
      )
  )
)

# Valid object?
validObject <- function(object) {
  TRUE

}
setValidity("QhatModel.homo.skewedNormal.linear", validObject)

# Initialise object
#setGeneric(name="initialize",def=function(.Object,input.data){standardGeneric("initialize")})
setMethod("initialize","QhatModel.homo.skewedNormal.linear", function(.Object, input.data, transition.graph=matrix(T,2,2),state.dependent.location.a0=T, state.dependent.location.a1=F, state.dependent.scale.a0=T, state.dependent.shape.a0=F) {
  .Object@input.data <- input.data

  .Object@nStates = ncol(transition.graph)

  # Set up model terms for mean and standard deviation.
  .Object@parameters$location.a0 = 0;
  .Object@parameters$location.a1 = 0.1;
  .Object@parameters$scale.a0 = 1;
  .Object@parameters$shape.a0 = 0;

  if (state.dependent.location.a0)
    .Object@parameters$location.a0 = rep(.Object@parameters$location.a0, .Object@nStates);
  if (state.dependent.location.a1)
    .Object@parameters$location.a1= rep(.Object@parameters$location.a1, .Object@nStates);
  if (state.dependent.scale.a0)
    .Object@parameters$scale.a0 = rep(.Object@parameters$scale.a0, .Object@nStates);
  if (state.dependent.shape.a0)
    .Object@parameters$shape.a0 = rep(.Object@parameters$shape.a0, .Object@nStates);

  validObject(.Object)
  .Object
}
)


# create a method to assign the parameter
#setGeneric(name="setParameters",def=function(.Object,parameters) {standardGeneric("setParameters")})
setMethod(f="setParameters",
          signature="QhatModel.homo.skewedNormal.linear",
          definition=function(.Object,parameters)
          {
            .Object@parameters$location.a1 <- parameters$location.a1
            .Object@parameters$location.a0 <- parameters$location.a0
            .Object@parameters$scale.a0 <- parameters$scale.a0
            .Object@parameters$shape.a0 <- parameters$shape.a0
            validObject(.Object)
            return(.Object)
          }
)
setMethod(f="setParameters.fromTransformed",
          signature="QhatModel.homo.skewedNormal.linear",
          definition=function(.Object,parameters)
          {
            .Object@parameters$location.a1 <- 10^parameters$location.a1
            .Object@parameters$location.a0 <- parameters$location.a0
            .Object@parameters$scale.a0 <- 10^parameters$scale.a0
            .Object@parameters$shape.a0 <- parameters$shape.a0
            #print(paste('Input loc. a0:',parameters$location.a0))

            return(.Object)
          }
)
# create a method to assign the parameter
setMethod(f="getParameters",
          signature="QhatModel.homo.skewedNormal.linear",
          definition=function(.Object)
          {
            return(list(location.a1=.Object@parameters$location.a1,
                        location.a0=.Object@parameters$location.a0,
                        scale.a0=.Object@parameters$scale.a0,
                        shape.a0 = .Object@parameters$shape.a0))
          }
)

setMethod(f="getDistributionPercentiles",
          signature=c("QhatModel.homo.skewedNormal.linear","data.frame","numeric"),
          definition=function(.Object, data, precentiles=c(0.5, 0.95))
          {

            # Check data is a data frame
            if (!is.data.frame(data))
              stop('Input "data" must be a a data frame.')

            # Check Qhat is in data
            if (!any(names(data)=="Qhat"))
              stop('Input "data" must be a a data frame with a variable named "Qhat".')

            # Get the skew dist location, shape, scale
            markov.mean = getMean(.Object, data)
            markov.scale = getScale(.Object, data)
            markov.shape = getShape(.Object, data)

            # Derive location parameters.
            markov.location <- markov.mean - markov.scale * markov.shape/sqrt(markov.shape^2+1) * sqrt(2/pi)

            # Calculate probabilities
            est <- vector('list',length(precentiles))
            for (k in 1:length(precentiles)) {
              est[[k]] = matrix(NA,nrow(data),.Object@nStates)
            }
            for (i in 1:nrow(data)) {
              for (j in 1:.Object@nStates) {
                for (k in 1:length(precentiles)) {
                  est[[k]][i,j] = sn::qsn(precentiles[k], xi=markov.location[i,j], omega=markov.scale[i,j], alpha=markov.shape[i,j])
                }
              }
            }
            #est = lapply(precentiles, sn::qsn, xi=markov.location, omega=markov.scale, alpha=markov.shape)
            #stop("DBG')")
            names(est) = precentiles
            return(est)
          }
)

setMethod(f="getTransformedParameterBounds",
          signature="QhatModel.homo.skewedNormal.linear",
          definition=function(.Object)
          {
            parameters = getParameters(.Object)
            lowerBound = list(location.a1 = rep(-6, length(.Object@parameters$location.a1)),
                              location.a0 = rep(-10, length(.Object@parameters$location.a0)),
                              scale.a0 = rep(log10(0.01), length(.Object@parameters$scale.a0)),
                              shape.a0 = rep(-5, length(.Object@parameters$shape.a0))
                              #shape.a0 = rep(-1e-6, length(.Object@parameters$shape.a0))
            )
            upperBound = list(location.a1 = rep(log10(0.1), length(.Object@parameters$location.a1)),
                              location.a0 = rep(10, length(.Object@parameters$location.a0)),
                              scale.a0 = rep(log10(1000), length(.Object@parameters$scale.a0)),
                              shape.a0 = rep(5, length(.Object@parameters$shape.a0))
                              #shape.a0 = rep(1e-6, length(.Object@parameters$shape.a0))
            )
            return(list(lower = lowerBound, upper = upperBound))
          }
)

# Get transition matrix with no input data.
setMethod(f="getEmissionDensity",
          signature=c("QhatModel.homo.skewedNormal.linear","data.frame"),
          definition=function(.Object, data)
          {

            # Check Qhat is in data
            if (!any(names(data)=="Qhat"))
              stop('Input "data" must be a a data frame with a variable named "Qhat".')

            # Get the skew dist location, shape, scale
            markov.mean = getMean(.Object, data)
            markov.scale = getScale(.Object, data)
            markov.shape = getShape(.Object, data)

            # Derive location parameters.
            markov.location <- markov.mean - markov.scale * markov.shape/sqrt(markov.shape^2+1) * sqrt(2/pi)
            #print(paste('Set   loc. a0:',.Object@parameters$location.a0))
            # Calculate probabilities.
            #P = dnorm(data$Qhat,markov.means,markov.stds)
            P <- matrix(NA, nrow(data),.Object@nStates)
            for (i in 1:.Object@nStates) {

              P[,i] <- sn::dsn(data$Qhat, xi=markov.location[,i], omega=markov.scale[,i], alpha=markov.shape[,i], tau=0, dp=NULL, log=FALSE)
            }
            #stop("DBG')")
            return(P)
          }
)
# Calculate the transformed flow at the mean annual precip
#setGeneric(name="getMean",def=function(.Object, data) {standardGeneric("getMean")})
setMethod(f="getMean",signature=c("QhatModel.homo.skewedNormal.linear","data.frame"),definition=function(.Object, data)
{
            ncols.a1 = length(.Object@parameters$location.a1)
            ncols.a0 = length(.Object@parameters$location.a0)
            nrows = length(data$precipitation);
            ncols.max = max(c(ncols.a0 ,ncols.a1))
            if (ncols.max > .Object@nStates)
              stop(paste('The number of parameters for each term of the mean model must must equal 1 or the number of states of ',.Object@nStates))

            # Check which terms are uniform for all states and whic terms are unique
            # to each state.
            if (ncols.a0==1 || ncols.a0==.Object@nStates) {
              a0.est = matrix(rep(.Object@parameters$location.a0,each=nrows),nrows,.Object@nStates);
            } else if (ncols.a0<.Object@nStates) {
              stop(paste('The number of parameters for the a0 term of the mean model must must equal 1 or the number of states of ',.Object@nStates))
            }
            if (ncols.a1==1 || ncols.a1==.Object@nStates) {
              a1.est = matrix(rep(.Object@parameters$location.a1,each=nrows),nrows,.Object@nStates);
            } else if (ncols.a1<.Object@nStates) {
              stop(paste('The number of parameters for the a1 term of the mean model must must equal 1 or the number of states of ',.Object@nStates))
            }

            precip.data = matrix(data$precipitation,nrows,.Object@nStates);

            # Calculate the non-AR1 componants
            Qhat.model <- precip.data * a1.est + a0.est

            return(Qhat.model)
          }
)

setGeneric(name="getScale",def=function(.Object, data) {standardGeneric("getScale")})
setMethod(f="getScale",signature=c("QhatModel.homo.skewedNormal.linear","data.frame"),definition=function(.Object, data)
{

  ncols.a0 = length(.Object@parameters$scale.a0)
  nrows = length(data$precipitation);

  a0.est = matrix(rep(.Object@parameters$scale.a0,each=nrows),nrows,.Object@nStates);

  return(a0.est)
}
)

setGeneric(name="getShape",def=function(.Object, data) {standardGeneric("getShape")})
setMethod(f="getShape",signature=c("QhatModel.homo.skewedNormal.linear","data.frame"),definition=function(.Object, data)
{
  ncols.a0 = length(.Object@parameters$shape.a0)
  nrows = length(data$precipitation);

  a0.est = matrix(rep(.Object@parameters$shape.a0,each=nrows),nrows,.Object@nStates);

  return(a0.est)
}
)

setMethod(f="generate.sample.Qhat",signature="QhatModel.homo.skewedNormal.linear",definition=function(.Object, data, sample.states)
{

  # Generate a synthtic series of Qhat usng the input random series of states
  nSamples = length(sample.states)
  sample.Qhat = rep(0,nSamples)
  moments = getMoments(.Object, data)
  if (nrow(moments$mean)!=nSamples)
    stop('The input vector od sample states must equal the numer of rows of the input data.')

  for (i in 1:nSamples)
    sample.Qhat[i] <-  sn::rsn(data$Qhat, xi=markov.location[i,sample.states[i]], omega=markov.scale[i,sample.states[i]], alpha=markov.shape[i,sample.states[i]], tau=0, dp=NULL)

  return(sample.Qhat)
}
)
