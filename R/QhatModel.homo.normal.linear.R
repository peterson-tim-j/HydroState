##' @include abstracts.R parameters.R
## @export
QhatModel.homo.normal.linear <- setClass(
  # Set the name for the class
  "QhatModel.homo.normal.linear",

  package='hydroState',

  contains=c('QhatModel'),

  # Define the slots
  slots = c(
    input.data = "data.frame",
    nStates = 'numeric',
    use.truncated.dist = 'logical',
    parameters = "parameters"
  ),

  # Set the default values for the slots. (optional)
  prototype=list(
    input.data = data.frame(year=c(0),month=c(0),precipitation=c(0)),
    nStates = Inf,
    use.truncated.dist=T,
    parameters = new('parameters',c('mean.a0', 'mean.a1','std.a0'),c(1,1,1))

  )
)

# Valid object?
validObject <- function(object) {
  TRUE

}
setValidity("QhatModel.homo.normal.linear", validObject)

# Initialise object
#setGeneric(name="initialize",def=function(.Object,input.data){standardGeneric("initialize")})
setMethod("initialize","QhatModel.homo.normal.linear", function(.Object, input.data, use.truncated.dist=T, transition.graph=matrix(T,2,2),
                                                                state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.mean.trend=NA, state.dependent.std.a0=T) {

  .Object@input.data <- input.data

  .Object@use.truncated.dist = use.truncated.dist

  .Object@nStates = ncol(transition.graph)

  # Set the number of parameter values per parameter name and set up model terms for mean and standard deviation and trend.
  if (is.na(state.dependent.mean.trend)) {
    parameter.length <- as.numeric(c(state.dependent.mean.a0, state.dependent.mean.a1, state.dependent.std.a0)) * (.Object@nStates-1) + 1
    .Object@parameters = new('parameters', c('mean.a0', 'mean.a1', 'std.a0'), parameter.length)
  } else {
    parameter.length <- as.numeric(c(state.dependent.mean.a0, state.dependent.mean.a1, state.dependent.mean.trend, state.dependent.std.a0)) * (.Object@nStates-1) + 1
    .Object@parameters = new('parameters', c('mean.a0', 'mean.a1', 'mean.trend', 'std.a0'), parameter.length)
  }

  validObject(.Object)
  .Object
}
)

setGeneric(name="get.SeaonalityPeriod",def=function(.Object,input.data){standardGeneric("get.SeaonalityPeriod")})
setMethod("get.SeaonalityPeriod","QhatModel.homo.normal.linear", function(.Object, input.data) {

  if (any(names(input.data)=="month")) {
  # Find seasonal period.
  # Adapted from https://stackoverflow.com/questions/12824931/measure-the-periodicity-of-a-sequence-of-numbers-r
  ii <- 0
  while (TRUE) {
    ii <- ii + 1
    LAG <- sum((diff(input.data$month, lag = ii) == 0) - 1)
    if (LAG == 0) { break }
  }
  seasonal.step.size = ii

  # Check seaonality is number of AR terms.
  if (seasonal.step.size<=1)
    stop(paste('The period of seasonality is', seasonal.step.size, ' input data rows. It must be greater than the maximum AR trm for the model.'))

  # Check seaonality is sensible.
  if (seasonal.step.size>12)
    warning(paste('The period of seasonality is', seasonal.step.size, ' input data rows. The module is designed for <=12 (i.e. monthly or less frequent).'))

  } else {
    seasonal.step.size = 0
  }

  return(seasonal.step.size)
}
)

# Get transition matrix with no input data.
setMethod(f="getEmissionDensity",
          signature=c("QhatModel.homo.normal.linear","data.frame"),
          definition=function(.Object, data, cumProb.threshold.Qhat)
          {

            # Check Qhat is in data
            if (!any(names(data)=="Qhat.flow"))
              stop('Input "data" must be a a data frame with a variable named "Qhat.flow".')

            if (!any(names(data)=="Qhat.precipitation"))
              stop('Input "data" must be a a data frame with a variable named "Qhat.precipitation".')

            # Get the moments
            markov.mean = getMean(.Object, data)
            markov.variance = getVariance(.Object, data)
            markov.stds <- sqrt(markov.variance)

            # For truncated flow, the mean  must >=0 else very negative means can arise
            if (.Object@use.truncated.dist && any(markov.mean<0,na.rm = T) ) {
                P = matrix(Inf, nrow(markov.mean), .Object@nStates)
                #markov.mean = pmax(0, markov.mean,na.rm = T)
            }

            # Calculate probabilities.
            P <- matrix(NA, nrow(data),.Object@nStates)

            # Set the lower limit of a truncated normal dist.
            lowerX = -Inf;
            if (.Object@use.truncated.dist)
              lowerX = 0;

            # get the prob.
            if (!all(is.na(cumProb.threshold.Qhat))) {
              if (length(cumProb.threshold.Qhat)!=nrow(data))
                stop('The length of cumProb.threshold.Qhat must equal the number of rows of input data')

              for (i in 1:.Object@nStates) {
                P[,i] = ptruncnorm(cumProb.threshold.Qhat[i], a=lowerX, mean=markov.mean[,i], sd=markov.stds[,i])
              }
            } else {
              for (i in 1:.Object@nStates) {
                P[,i] = dtruncnorm(data$Qhat.flow, a=lowerX, mean=markov.mean[,i], sd=markov.stds[,i])
              }
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
            if (!any(names(data)=="Qhat.flow"))
              stop('Input "data" must be a a data frame with a variable named "Qhat.flow".')

            # Get the moments
            markov.mean = getMean(.Object, data)
            markov.variance = getVariance(.Object, data)
            markov.stds <- sqrt(markov.variance)

            # Initialise returned variable
            est = vector('list',length(precentiles))

            # Set the lower limit of a truncated normal dist.
            lowerX = -Inf;
            if (.Object@use.truncated.dist)
              lowerX = 0;

            # Calculate probabilities.
            for (i in 1:length(precentiles)) {
              est.tmp = matrix(0,nrow(markov.mean),ncol(markov.mean))
              for (j in 1:.Object@nStates) {
                est.tmp[,j] <- qtruncnorm(precentiles[i], a=lowerX, mean= markov.mean[,j], sd=markov.stds[,j])
              }
              est[[i]] = est.tmp
            }


            names(est) = precentiles
            return(est)
          }
)

# Calculate the transformed flow at the mean annual precip
setGeneric(name="getMean",def=function(.Object, data) {standardGeneric("getMean")})
setMethod(f="getMean",signature=c("QhatModel.homo.normal.linear","data.frame"),definition=function(.Object, data)
{

            # Get object parameter list
            parameters = getParameters(.Object@parameters)

            ncols.a1 = length(parameters$mean.a1)
            ncols.a0 = length(parameters$mean.a0)
            ncols.trend = 0
            if ('mean.trend' %in% names(parameters)) {
              ncols.trend = length(parameters$mean.trend)
            }
            nrows = length(data$Qhat.precipitation);
            ncols.max = max(c(ncols.a0 ,ncols.a1, ncols.trend))

            if (ncols.max > .Object@nStates)
              stop(paste('The number of parameters for each term of the mean model must must equal 1 or the number of states of ',.Object@nStates))

            # Check which terms are uniform for all states and whic terms are unique
            # to each state.
            if (ncols.a0==1 || ncols.a0==.Object@nStates) {
              a0.est = matrix(rep(parameters$mean.a0,each=nrows),nrows,.Object@nStates);
            } else if (ncols.a0<.Object@nStates) {
              stop(paste('The number of parameters for the a0 term of the mean model must must equal 1 or the number of states of ',.Object@nStates))
            }
            if (ncols.a1==1 || ncols.a1==.Object@nStates) {
              a1.est = matrix(rep(parameters$mean.a1,each=nrows),nrows,.Object@nStates);
            } else if (ncols.a1<.Object@nStates) {
              stop(paste('The number of parameters for the a1 term of the mean model must must equal 1 or the number of states of ',.Object@nStates))
            }
            if (ncols.trend==1 || ncols.trend==.Object@nStates) {
              trend.est = matrix(rep(parameters$mean.trend,each=nrows),nrows,.Object@nStates);
            } else {
              trend.est = 0
            }

            time.vals = matrix(data$year - data$year[1],nrows,.Object@nStates)
            precip.data = matrix(data$Qhat.precipitation,nrows,.Object@nStates);

            # Calculate the non-AR1 componants
            a0.est <- 100 * a0.est
            Qhat.model <- precip.data * a1.est + a0.est + time.vals * trend.est

            # print(paste('...DBG getMean.AR0 nrows Qhat.model.NAs:',nrow(Qhat.model)))

            return(Qhat.model)
          }
)

# Get variance. It is calculates as a parameter times the Qhat variance.
setGeneric(name="getVariance",def=function(.Object, data) {standardGeneric("getVariance")})
setMethod(f="getVariance",signature=c("QhatModel.homo.normal.linear","data.frame"),definition=function(.Object, data)
{
  # Get object parameter list
  parameters = getParameters(.Object@parameters)

  ncols.a0 = length(parameters$std.a0)
  nrows = length(data$Qhat.precipitation);

  # Get variance of the Qhat
  Qhat.var = var(data$Qhat.flow, na.rm=T)

  a0.est = Qhat.var * matrix(rep(parameters$std.a0,each=nrows),nrows,.Object@nStates);
  #precip.data = rep(.Object@input.data$Qhat.precipitation,1,ncols.a0);

  return(a0.est)

}
)

setMethod(f="generate.sample.Qhat.fromViterbi",signature=c("QhatModel.homo.normal.linear",'data.frame','numeric'),definition=function(.Object, data, viterbi.states)
{

  if (length(viterbi.states) != nrow(data))
    stop('The length of the viterbi.states must equal the number of rows in data.')

  # Generate a synthtic series of Qhat usng the input random series of states
  nSamples = length(viterbi.states)
  sample.Qhat = rep(0,nSamples)

  markov.mean = getMean(.Object, data)
  markov.variance = getVariance(.Object, data)

  if (nrow(markov.mean)!=nSamples)
    stop('The input vector of sample states must equal the numer of rows of the input data.')

  # Set the lower limit of a truncated normal dist.
  lowerX = -Inf;
  if (.Object@use.truncated.dist)
    lowerX = 0;

  for (i in 1:nSamples)
    sample.Qhat[i] = rtruncnorm(1, a=lowerX, mean = markov.mean[i,viterbi.states[i]], sd = sqrt(markov.variance[i,viterbi.states[i]]))

  return(sample.Qhat)
}
)

setMethod(f="generate.sample.Qhat",signature="QhatModel.homo.normal.linear",definition=function(.Object, data, nSamples)
{

  if (length(nSamples)!=1 || nSamples<=0)
    stop('The input nSamples must be a single number >=1.')

  # Get distribution for each state
  markov.mean = getMean(.Object, data)
  markov.variance = getVariance(.Object, data)

  # Set the lower limit of a truncated normal dist.
  lowerX = -Inf;
  if (.Object@use.truncated.dist)
    lowerX = 0;

  sample.Qhat = vector('list',.Object@nStates)
  for (i in 1:.Object@nStates) {
    sample.Qhat[[i]] = matrix(NA,nrow(data),nSamples)
    for (j in 1:nrow(data)) {
      sample.Qhat[[i]][j,] = rtruncnorm(nSamples, a=lowerX, mean = markov.mean[j,i], sd = sqrt(markov.variance[j,i]))
    }
  }

  return(sample.Qhat)
}
)

