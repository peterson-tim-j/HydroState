##' @include abstracts.R QhatModel.homo.normal.linear.R
## @export
QhatModel.homo.skewedNormal.linear <- setClass(
  # Set the name for the class
  "QhatModel.homo.skewedNormal.linear",

  package='hydroState',

  contains=c('QhatModel.homo.normal.linear'),

    # Set the default values for the slots. (optional)
  prototype=list(
    input.data = data.frame(year=c(0),month=c(0),precipitation=c(0)),
    nStates = Inf,
    use.truncated.dist=F,
    parameters = new('parameters',c('mean.a0', 'mean.a1','std.a0', 'shape.a0'),c(1,1,1,1))
  )
)

# Valid object?
validObject <- function(object) {
  TRUE

}
setValidity("QhatModel.homo.skewedNormal.linear", validObject)

# Initialise object
#setGeneric(name="initialize",def=function(.Object,input.data){standardGeneric("initialize")})
setMethod("initialize","QhatModel.homo.skewedNormal.linear", function(.Object, input.data, transition.graph=matrix(T,2,2),state.dependent.mean.a0=T,
                                                                      state.dependent.mean.a1=F, state.dependent.std.a0=T, state.dependent.shape.a0=T) {
  .Object@input.data <- input.data
  .Object@use.truncated.dist <- F
  .Object@nStates = ncol(transition.graph)

  # Set the number of parameter values per parameter name.
  parameter.length <- as.numeric(c(state.dependent.mean.a0, state.dependent.mean.a1, state.dependent.std.a0, state.dependent.shape.a0)) * (.Object@nStates-1) + 1

  # Set up model terms for mean and standard deviation.
  .Object@parameters = new('parameters', c('mean.a0', 'mean.a1', 'std.a0','shape.a0'), parameter.length)
  validObject(.Object)
  .Object
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
            if (!any(names(data)=="Qhat.flow"))
              stop('Input "data" must be a a data frame with a variable named "Qhat".')

            # Get the skew dist location, shape, scale
            markov.mean = getMean(.Object, data)
            markov.variance = getVariance(.Object, data)
            markov.shape = getShape(.Object, data)

            # Derive location parameters.
            markov.location <- markov.mean - markov.variance * markov.shape/sqrt(markov.shape^2+1) * sqrt(2/pi)

            # Calculate probabilities
            est <- vector('list',length(precentiles))
            for (k in 1:length(precentiles)) {
              est[[k]] = matrix(NA,nrow(data),.Object@nStates)
            }
            for (i in 1:nrow(data)) {
              for (j in 1:.Object@nStates) {
                for (k in 1:length(precentiles)) {
                  est[[k]][i,j] = sn::qsn(precentiles[k], xi=markov.location[i,j], omega=markov.variance[i,j], alpha=markov.shape[i,j])
                }
              }
            }
            #est = lapply(precentiles, sn::qsn, xi=markov.location, omega=markov.scale, alpha=markov.shape)
            #stop("DBG')")
            names(est) = precentiles
            return(est)
          }
)

# Get transition matrix with no input data.
setMethod(f="getEmissionDensity",
          signature=c("QhatModel.homo.skewedNormal.linear","data.frame"),
          definition=function(.Object, data, cumProb.threshold.Qhat)
          {

            # Check Qhat is in data
            if (!any(names(data)=="Qhat.flow"))
              stop('Input "data" must be a a data frame with a variable named "Qhat.flow".')

            if (!any(names(data)=="Qhat.precipitation"))
              stop('Input "data" must be a a data frame with a variable named "Qhat.precipitation".')

            # Get the skew dist location, shape, scale
            markov.mean = getMean(.Object, data)
            markov.variance = getVariance(.Object, data)
            markov.shape = getShape(.Object, data)

            # Derive location parameters.
            markov.location <- markov.mean - markov.variance * markov.shape/sqrt(markov.shape^2+1) * sqrt(2/pi)

            # Calculate probabilities.
            P <- matrix(NA, nrow(data),.Object@nStates)
            if (!all(is.na(cumProb.threshold.Qhat))) {

              if (length(cumProb.threshold.Qhat)!=nrow(data))
                stop('The length of cumProb.threshold.Qhat must equal the number of rows of input data')

              for (i in 1:.Object@nStates) {
                for (j in 1:nrow(data)) {
                  if (is.na(data$Qhat.flow[j])) {
                    P[j,i] <- NA
                  } else {
                    P[j,i] <- sn::psn(cumProb.threshold.Qhat[j], xi=markov.location[j,i], omega=markov.variance[j,i], alpha=markov.shape[j,i], tau=0, dp=NULL, log=FALSE)
                  }
                }
              }
            } else {
                for (i in 1:.Object@nStates) {
                  P[,i] <- sn::dsn(data$Qhat.flow, xi=markov.location[,i], omega=markov.variance[,i], alpha=markov.shape[,i], tau=0, dp=NULL, log=FALSE)
                }
            }

            return(P)
          }
)

setGeneric(name="getShape",def=function(.Object, data) {standardGeneric("getShape")})
setMethod(f="getShape",signature=c("QhatModel.homo.skewedNormal.linear","data.frame"),definition=function(.Object, data)
{
  parameters = getParameters(.Object@parameters)
  ncols.a0 = length(parameters$shape.a0)
  nrows = length(data$Qhat.precipitation);

  a0.est = matrix(rep(parameters$shape.a0,each=nrows),nrows,.Object@nStates);

  return(a0.est)
}
)

setMethod(f="generate.sample.Qhat.fromViterbi",signature=c("QhatModel.homo.skewedNormal.linear",'data.frame','numeric'),definition=function(.Object, data, viterbi.states)
{

  if (length(viterbi.states) != nrow(data))
    stop('The length of the viterbi.states must equal the number of rows in data.')

  # Generate a synthtic series of Qhat usng the input random series of states
  nSamples = length(viterbi.states)
  sample.Qhat = rep(0,nSamples)

  # Get the skew dist location, shape, scale
  markov.mean = getMean(.Object, data)
  markov.variance = getVariance(.Object, data)
  markov.shape = getShape(.Object, data)

  # Derive location parameters.
  markov.location <- markov.mean - markov.variance * markov.shape/sqrt(markov.shape^2+1) * sqrt(2/pi)

  for (i in 1:nSamples)
    sample.Qhat[i] <-  sn::rsn(1, xi=markov.location[i,viterbi.states[i]],
                               omega=markov.variance[i,viterbi.states[i]], alpha=markov.shape[i,viterbi.states[i]], tau=0, dp=NULL)

  return(sample.Qhat)
}
)


setMethod(f="generate.sample.Qhat",signature="QhatModel.homo.skewedNormal.linear",definition=function(.Object, data, nSamples)
{

  if (length(nSamples)!=1 || nSamples<=0)
    stop('The input nSamples must be a single number >=1.')

  # Get the skew dist location, shape, scale
  markov.mean = getMean(.Object, data)
  markov.variance = getVariance(.Object, data)
  markov.shape = getShape(.Object, data)

  # Derive location parameters.
  markov.location <- markov.mean - markov.variance * markov.shape/sqrt(markov.shape^2+1) * sqrt(2/pi)

  sample.Qhat = vector('list',.Object@nStates)
  for (i in 1:.Object@nStates) {
    sample.Qhat[[i]] = matrix(NA,nrow(data),nSamples)
    for (j in 1:nrow(data)) {
      sample.Qhat[[i]][j,] = sn::rsn(nSamples, xi=markov.location[j, i],
                                     omega=markov.variance[j, i], alpha=markov.shape[j, i], tau=0, dp=NULL)
    }
  }

  return(sample.Qhat)
}
)
