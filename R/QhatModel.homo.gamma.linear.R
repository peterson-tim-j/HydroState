##' @include abstracts.R QhatModel.homo.normal.linear.R
## @export
QhatModel.homo.gamma.linear <- setClass(
  # Set the name for the class
  "QhatModel.homo.gamma.linear",

  package='hydroState',

  contains=c('QhatModel.homo.normal.linear'),

    # Set the default values for the slots. (optional)
  prototype=list(
    input.data = data.frame(year=c(0),month=c(0),precipitation=c(0)),
    nStates = Inf,
    use.truncated.dist=F,
    parameters = new('parameters',c('mean.a0', 'mean.a1','std.a0'),c(1,1,1))
  )
)

# Valid object?
validObject <- function(object) {
  TRUE

}
setValidity("QhatModel.homo.gamma.linear", validObject)

# Initialise object
#setGeneric(name="initialize",def=function(.Object,input.data){standardGeneric("initialize")})
setMethod("initialize","QhatModel.homo.gamma.linear", function(.Object, input.data, transition.graph=matrix(T,2,2),state.dependent.mean.a0=T,
                                                                      state.dependent.mean.a1=F, state.dependent.mean.trend=NA,state.dependent.std.a0=T) {

  .Object@input.data <- input.data
  .Object@use.truncated.dist <- F
  .Object@nStates = ncol(transition.graph)

  # Set the number of parameter values per parameter name.
  parameter.length <- as.numeric(c(state.dependent.mean.a0, state.dependent.mean.a1, state.dependent.std.a0)) * (.Object@nStates-1) + 1

  # Set up model terms for mean and standard deviation.
  .Object@parameters = new('parameters', c('mean.a0', 'mean.a1', 'std.a0'), parameter.length)
  validObject(.Object)
  .Object
}
)

setMethod(f="getDistributionPercentiles",
          signature=c("QhatModel.homo.gamma.linear","data.frame","numeric"),
          definition=function(.Object, data, precentiles=c(0.5, 0.95))
          {

            # Check data is a data frame
            if (!is.data.frame(data))
              stop('Input "data" must be a a data frame.')

            # Check Qhat is in data
            if (!any(names(data)=="Qhat.flow"))
              stop('Input "data" must be a a data frame with a variable named "Qhat".')

            # Get distribution mean, dispersion and familty parameters.
            markov.mean = getMean(.Object, data)
            markov.variance = getVariance(.Object, data)

            # Limit mean to >0
            filt = is.finite(markov.mean) & markov.mean <= sqrt(.Machine$double.eps)*1000
            markov.mean[filt] = sqrt(.Machine$double.eps)*1000

            # Derive the Gamme parametesr from the modelled mean and variance.
            markov.shape <- markov.mean^2/markov.variance
            markov.scale <- markov.variance/markov.mean

            # Calculate probabilities
            est <- vector('list',length(precentiles))
            for (k in 1:length(precentiles)) {
              est[[k]] = matrix(NA,nrow(data),.Object@nStates)
            }
            for (i in 1:nrow(data)) {
              for (j in 1:.Object@nStates) {
                for (k in 1:length(precentiles)) {
                  if (is.na(markov.mean[i,j])) {
                    est[[k]][i,j] = NA
                  } else {
                    est[[k]][i,j] = qgamma(precentiles[k], shape=markov.shape[i,j], scale=markov.scale[i,j])
                  }

                }
              }
            }

            names(est) = precentiles
            return(est)
          }
)

# Get transition matrix with no input data.
setMethod(f="getEmissionDensity",
          signature=c("QhatModel.homo.gamma.linear","data.frame"),
          definition=function(.Object, data, cumProb.threshold.Qhat)
          {

            # Check Qhat is in data
            if (!any(names(data)=="Qhat.flow"))
              stop('Input "data" must be a a data frame with a variable named "Qhat.flow".')

            if (!any(names(data)=="Qhat.precipitation"))
              stop('Input "data" must be a a data frame with a variable named "Qhat.precipitation".')

            # Get the Burr mean, dispersion and familty parameters.
            markov.mean = getMean(.Object, data)
            markov.variance = getVariance(.Object, data)

            # Limit mean to >0
            filt = is.finite(markov.mean) & markov.mean <= sqrt(.Machine$double.eps)*1000
            markov.mean[filt] = sqrt(.Machine$double.eps)*1000

            # Derive the Gamme parametesr from the modelled mean and variance.
            markov.shape <- markov.mean^2/markov.variance
            markov.scale <- markov.variance/markov.mean

            # Find non NAs
            filt <- !is.na(data$Qhat.flow)

             # Initialise returned probs., P
            P <- matrix(NA, nrow(data),.Object@nStates)

            # Increase zero values to machine precision.
            data$Qhat.flow[data$Qhat.flow==0] = sqrt(.Machine$double.eps)

            # Calculate probabilities.
            if (!all(is.na(cumProb.threshold.Qhat))) {

              if (length(cumProb.threshold.Qhat)!=nrow(data))
                stop('The length of cumProb.threshold.Qhat must equal the number of rows of input data')

              # Increase zero values to machine precision.
              cumProb.threshold.Qhat[cumProb.threshold.Qhat==0] = sqrt(.Machine$double.eps)

              for (i in 1:.Object@nStates) {
                for (j in 1:nrow(data)) {
                  if (is.na(cumProb.threshold.Qhat[j])) {
                    P[j,i] <- NA
                  } else {
                    P[j,i] <- pgamma(cumProb.threshold.Qhat[j], shape=markov.shape[j,i], scale=markov.scale[j,i])
                  }
                }
              }
            } else {
                for (i in 1:.Object@nStates) {
                  P[,i] <- dgamma(data$Qhat.flow, shape=markov.shape[,i], scale=markov.scale[,i])
                }
            }

            return(P)
          }
)

setMethod(f="generate.sample.Qhat.fromViterbi",signature=c("QhatModel.homo.gamma.linear",'data.frame','numeric'),definition=function(.Object, data, viterbi.states)
{

  if (length(viterbi.states) != nrow(data))
    stop('The length of the viterbi.states must equal the number of rows in data.')

  # Generate a synthtic series of Qhat usng the input random series of states
  nSamples = length(viterbi.states)
  sample.Qhat = rep(0,nSamples)

  # Get the Burr mean, dispersion and familty parameters.
  markov.mean = getMean(.Object, data)
  markov.variance = getVariance(.Object, data)

  # Limit mean to >0
  filt = is.finite(markov.mean) & markov.mean <= sqrt(.Machine$double.eps)*1000
  markov.mean[filt] = sqrt(.Machine$double.eps)*1000

  # Derive the Gamme parametesr from the modelled mean and variance.
  markov.shape <- markov.mean^2/markov.variance
  markov.scale <- markov.variance/markov.mean

  if (nrow(markov.shape)!=nSamples)
    stop('The input vector of sample states must equal the numer of rows of the input data.')

  for (i in 1:nSamples)
    sample.Qhat[i] <-  rgamma(1, shape=markov.shape[i,viterbi.states[i]], scale=markov.scale[i,viterbi.states[i]])

  return(sample.Qhat)
}
)


setMethod(f="generate.sample.Qhat",signature="QhatModel.homo.gamma.linear",definition=function(.Object, data, nSamples)
{

  if (length(nSamples)!=1 || nSamples<=0)
    stop('The input nSamples must be a single number >=1.')

  # Get the Burr mean, dispersion and familty parameters.
  markov.mean = getMean(.Object, data)
  markov.variance = getVariance(.Object, data)

  # Limit mean to >0
  filt = is.finite(markov.mean) & markov.mean <= sqrt(.Machine$double.eps)*1000
  markov.mean[filt] = sqrt(.Machine$double.eps)*1000

  # Derive the Gamme parametesr from the modelled mean and variance.
  markov.shape <- markov.mean^2/markov.variance
  markov.scale <- markov.variance/markov.mean

  sample.Qhat = vector('list',.Object@nStates)
  for (i in 1:.Object@nStates) {
    sample.Qhat[[i]] = matrix(NA,nrow(data),nSamples)
    for (j in 1:nrow(data)) {
      sample.Qhat[[i]][j,] = rgamma(nSamples, shape=markov.shape[j,i], scale=markov.scale[j,i])
    }
  }

  return(sample.Qhat)
}
)
