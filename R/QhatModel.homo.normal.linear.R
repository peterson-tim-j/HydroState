##' @include abstracts.R
##' @export
QhatModel.homo.normal.linear <- setClass(
  # Set the name for the class
  "QhatModel.homo.normal.linear",

  package='hydroState',

  contains=c('QhatModel'),

  # Define the slots
  slots = c(
    input.data = "data.frame",
    nStates = 'numeric',
    parameters = "list"
  ),

  # Set the default values for the slots. (optional)
  prototype=list(
    input.data = data.frame(year=c(0),month=c(0),precipitation=c(0)),
    nStates = Inf,
    parameters = list(
      mean.a1 = 1.0,
      mean.a0 = 1.0,
      std.a0 = 1.0
      )
  )
)

# Valid object?
validObject <- function(object) {
  TRUE

}
setValidity("QhatModel.homo.normal.linear", validObject)

# Initialise object
#setGeneric(name="initialize",def=function(.Object,input.data){standardGeneric("initialize")})
setMethod("initialize","QhatModel.homo.normal.linear", function(.Object, input.data, transition.graph=matrix(T,2,2),state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T) {
  .Object@input.data <- input.data

  .Object@nStates = ncol(transition.graph)

  # Set up model terms for mean and standard deviation.
  .Object@parameters$mean.a0 = 0;
  .Object@parameters$mean.a1 = 0.1;
  .Object@parameters$std.a0 = 100;
  if (state.dependent.mean.a0)
    .Object@parameters$mean.a0 = rep(.Object@parameters$mean.a0, .Object@nStates);
  if (state.dependent.mean.a1)
    .Object@parameters$mean.a1 = rep(.Object@parameters$mean.a1, .Object@nStates);
  if (state.dependent.std.a0)
    .Object@parameters$std.a0 = rep(.Object@parameters$std.a0, .Object@nStates);

  validObject(.Object)
  .Object
}
)


# create a method to assign the parameter
#setGeneric(name="setParameters",def=function(.Object,parameters) {standardGeneric("setParameters")})
setMethod(f="setParameters",
          signature="QhatModel.homo.normal.linear",
          definition=function(.Object,parameters)
          {
            .Object@parameters$mean.a1 <- parameters$mean.a1
            .Object@parameters$mean.a0 <- parameters$mean.a0
            .Object@parameters$std.a0 <- parameters$std.a0
            validObject(.Object)
            return(.Object)
          }
)
setMethod(f="setParameters.fromTransformed",
          signature="QhatModel.homo.normal.linear",
          definition=function(.Object,parameters)
          {
            .Object@parameters$mean.a1 <- 10^parameters$mean.a1
            .Object@parameters$mean.a0 <- parameters$mean.a0*100
            .Object@parameters$std.a0 <- 10^parameters$std.a0
            return(.Object)
          }
)
# create a method to assign the parameter
#setGeneric(name="getParameters",def=function(.Object,position) {standardGeneric("getParameters")})
setMethod(f="getParameters",
          signature="QhatModel.homo.normal.linear",
          definition=function(.Object)
          {
            return(list(mean.a1=.Object@parameters$mean.a1,mean.a0=.Object@parameters$mean.a0, std.a0=.Object@parameters$std.a0))
          }
)

# Get transition matrix with no input data.
setMethod(f="getEmissionDensity",
          signature=c("QhatModel.homo.normal.linear","data.frame"),
          definition=function(.Object, data)
          {

            # Check Qhat is in data
            if (!any(names(data)=="Qhat"))
              stop('Input "data" must be a a data frame with a variable named "Qhat".')

            # Get the moments
            markov.mean = getMean(.Object, data)
            markov.variance = getVariance(.Object, data)
            markov.stds <- sqrt(markov.variance)

            # Calculate probabilities.
            #P = dnorm(data$Qhat,markov.means,markov.stds)
            P <- matrix(NA, nrow(data),.Object@nStates)
            for (i in 1:.Object@nStates) {
              P[,i] = dnorm(data$Qhat,markov.mean[,i],markov.stds[,i])
            }

            return(P)
          }
)

setMethod(f="getDistributionPercentiles",
          signature=c("QhatModel.homo.normal.linear","data.frame","numeric"),
          definition=function(.Object, data, precentiles=c(0.5, 0.95))
          {

            # Check data is a data frame
            if (!is.data.frame(data))
              stop('Input "data" must be a a data frame.')

            # Check Qhat is in data
            if (!any(names(data)=="Qhat"))
              stop('Input "data" must be a a data frame with a variable named "Qhat".')

            # Get the moments
            markov.mean = getMean(.Object, data)
            markov.variance = getVariance(.Object, data)
            markov.stds <- sqrt(markov.variance)

            # Calculate probabilities.
            est = lapply(precentiles, qnorm, mean= markov.mean, sd=markov.stds)
            names(est) = precentiles
            return(est)
          }
)

setMethod(f="getTransformedParameterBounds",
          signature="QhatModel.homo.normal.linear",
          definition=function(.Object)
          {
            parameters = getParameters(.Object)
            lowerBound = list(mean.a1 = rep(-6, length(.Object@parameters$mean.a1)),
                              mean.a0 = rep(-3, length(.Object@parameters$mean.a0)),
                              std.a0 = rep(log10(0.01), length(.Object@parameters$std.a0))
            )
            upperBound = list(mean.a1 = rep(log10(0.1), length(.Object@parameters$mean.a1)),
                              mean.a0 = rep(3, length(.Object@parameters$mean.a0)),
                              std.a0 = rep(log10(1000), length(.Object@parameters$std.a0))
            )
            return(list(lower = lowerBound, upper = upperBound))
          }
)


# Calculate the transformed flow at the mean annual precip
setGeneric(name="getMean",def=function(.Object, data) {standardGeneric("getMean")})
setMethod(f="getMean",signature=c("QhatModel.homo.normal.linear","data.frame"),definition=function(.Object, data)
{
            ncols.a1 = length(.Object@parameters$mean.a1)
            ncols.a0 = length(.Object@parameters$mean.a0)
            nrows = length(data$precipitation);
            ncols.max = max(c(ncols.a0 ,ncols.a1))
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

            precip.data = matrix(data$precipitation,nrows,.Object@nStates);

            # Calculate the non-AR1 componants
            Qhat.model <- precip.data * a1.est + a0.est

            return(Qhat.model)
          }
)

setGeneric(name="getVariance",def=function(.Object, data) {standardGeneric("getVariance")})
setMethod(f="getVariance",signature=c("QhatModel.homo.normal.linear","data.frame"),definition=function(.Object, data)
{

  ncols.a0 = length(.Object@parameters$std.a0)
  nrows = length(data$precipitation);

  a0.est = matrix(rep(.Object@parameters$std.a0,each=nrows),nrows,.Object@nStates);
  #precip.data = rep(.Object@input.data$precipitation,1,ncols.a0);

  return(a0.est^2)

}
)

setMethod(f="generate.sample.Qhat",signature="QhatModel.homo.normal.linear",definition=function(.Object, data, sample.states)
{

  # Generate a synthtic series of Qhat usng the input random series of states
  nSamples = length(sample.states)
  sample.Qhat = rep(0,nSamples)

  markov.mean = getMean(.Object, data)
  markov.variance = getVariance(.Object, data)

  if (nrow(markov.mean)!=nSamples)
    stop('The input vector of sample states must equal the numer of rows of the input data.')

  for (i in 1:nSamples)
    sample.Qhat[i] = rnorm(1, mean = markov.mean[i,sample.states[i]], sd = sqrt(markov.variance[i,sample.states[i]]))

  return(sample.Qhat)
}
)

# setMethod(f="getZScore.atObs",signature="QhatModel.homo.normal.linear",definition=function(.Object, data)
# {
#   # Get the emision probs from the input data
#   emissionProbs <- getEmissionProbabilities(.Object, data)
#
#   # Get the standard deviation for each state
#   sigma2 <- getVariance(.Object, data)
#   sigma2 <- unique(sigma2)
#   if (nrow(sigma2)>1)
#     stop("getStateDistributionProb: This function requires that the variance for each state is constant over time.")
#   sigma = as.vector(sqrt(sigma2))
#
#   # With the asumption of a zero mean, derive the Z-score value from the emission probs.
#   emission.zscores = matrix(0,nrow(data),.Object@nStates)
#   for (i in 1:.Object@nStates)
#     emission.zscores[,i] <- qnorm(emissionProbs[,i], mean=0, sd = sigma[i])
#
#   return(emission.zscores)
# }
# )

# setMethod(f="getZScoreProbabilities.atIncrements",signature="QhatModel.homo.normal.linear",definition=function(.Object, data, nIncrements)
# {
#   # PROBLEMS: This function should give the prob. of increments of Qhat (from obs data) for each state.
#
#   # Get range in precipitation and derive vector of length nIncrements
#
#   #
#
#   #----------------------------------------
#   # Get the emision probs from the input data
#   emissionProbs <- getEmissionProbabilities(.Object, data)
#
#   # Get the standard deviation for each state
#   sigma2 <- getVariance(.Object, data)
#   sigma2 <- unique(sigma2)
#   if (nrow(sigma2)>1)
#     stop("getStateDistributionProb: This function requires that the variance for each state is constant over time.")
#   sigma = as.vector(sqrt(sigma2))
#
#   # With the asumption of a zero mean, derive the Z-score value from the emission probs.
#   emission.zscores = matrix(0,nrow(data),.Object@nStates)
#   for (i in 1:.Object@nStates)
#     emission.zscores[,i] <- qnorm(emissionProbs[,i], mean=0, sd = sigma[i])
#
#   # Derive a vector of z-score values that covers the range from min to max with length of nIncrements
#   zscore.max = max(emission.zscores,na.rm = TRUE)
#   zscore.min = min(emission.zscores,na.rm = TRUE)
#   zscore.increments = seq(zscore.min, zscore.max, length.out = nIncrements)
#
#   # Derive the probabilities of the z-score increments for each state (assuming zero mean)
#   zscore.probs = matrix(0,length(zscore.increments),.Object@nStates)
#   for (i in 1:.Object@nStates)
#     zscore.probs[,i] <- pnorm(zscore.increments, mean=0, sd = sigma[i])
#
#   if (any(is.na(zscore.probs)))
#     stop('NAs in emission.zscores')
#
#   # Return zscore increments and probs
#   zscore.probs = cbind(zscore.increments, zscore.probs)
#   return(zscore.probs)
# }
# )
