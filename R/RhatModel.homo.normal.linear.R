##' @include abstracts.R parameters.R
## @export
RhatModel.homo.normal.linear <- setClass(
  # Set the name for the class
  "RhatModel.homo.normal.linear",

  package='hydroState',

  contains=c('QhatModel'),

  # Define the slots
  slots = c(
    input.data = "data.frame",
    nStates = 'numeric',
    use.truncated.dist = 'logical',
    parameters = "parameters",
    subAnnual.Monthly.Steps = 'numeric'
  ),

  # Set the default values for the slots. (optional)
  prototype=list(
    input.data = data.frame(year=c(0),month=c(0),flow=c(0)),
    nStates = Inf,
    use.truncated.dist=T,
    parameters = new('parameters',c('mean.a0','mean.trend', 'std.a0', 'std.a1'),c(1,1,1,1)),
    subAnnual.Monthly.Steps = c(2, 5, 8, 11)
  )
)

# Valid object?
validObject <- function(object) {
  TRUE

}
setValidity("RhatModel.homo.normal.linear", validObject)

# Initialise object
#setGeneric(name="initialize",def=function(.Object,input.data){standardGeneric("initialize")})
setMethod("initialize","RhatModel.homo.normal.linear", function(.Object, use.truncated.dist=T, input.data, transition.graph=matrix(T,2,2),
                                                                state.dependent.mean.a0=T,
                                                                state.dependent.mean.trend = NA,
                                                                state.dependent.std.a0=T,
                                                                state.dependent.std.a1=0) {
  .Object@input.data <- input.data
  .Object@use.truncated.dist = use.truncated.dist
  .Object@nStates = ncol(transition.graph)

  # Set the number of parameter values per parameter name and set up model terms for mean and standard deviation and trend.
  if (state.dependent.std.a1==0) {
    if (is.na(state.dependent.mean.trend)) {
      parameter.length <- as.numeric(c(state.dependent.mean.a0, state.dependent.std.a0)) * (.Object@nStates-1) + 1
      .Object@parameters = new('parameters', c('mean.a0', 'std.a0'), parameter.length)
    } else {
      parameter.length <- as.numeric(c(state.dependent.mean.a0, state.dependent.mean.trend, state.dependent.std.a0)) * (.Object@nStates-1) + 1
      .Object@parameters = new('parameters', c('mean.a0', 'mean.trend', 'std.a0'), parameter.length)
    }
  } else {
    if (is.na(state.dependent.mean.trend)) {
      parameter.length <- as.numeric(c(state.dependent.mean.a0, state.dependent.std.a0, state.dependent.std.a1)) * (.Object@nStates-1) + 1
      .Object@parameters = new('parameters', c('mean.a0', 'std.a0', 'std.a1'), parameter.length)
    } else {
      parameter.length <- as.numeric(c(state.dependent.mean.a0, state.dependent.mean.trend, state.dependent.std.a0, state.dependent.std.a1)) * (.Object@nStates-1) + 1
      .Object@parameters = new('parameters', c('mean.a0', 'mean.trend', 'std.a0', 'std.a1'), parameter.length)
    }
  }
    validObject(.Object)
  .Object
}
)

# setGeneric(name="is.stationary",def=function(.Object) {standardGeneric("is.stationary")})
# setMethod(f="is.stationary",signature=c("RhatModel.homo.normal.linear"),definition=function(.Object)
# {
#   # Get object parameter list
#   parameters = getParameters(.Object@parameters)
#
#   if (all(parameters$mean.AR1 >-1 & parameters$mean.AR1 <1)) {
#     return(TRUE)
#   } else {
#     return(FALSE)
#   }
# }
# )

# Calculate the transformed flow at the mean annual precip
# setGeneric(name="getMean",def=function(.Object, data) {standardGeneric("getMean")})
setMethod(f="getMean",signature=c("RhatModel.homo.normal.linear","data.frame"),definition=function(.Object, data)
{

  # Get object parameter list
  parameters = getParameters(.Object@parameters)

  ncols.a0 = length(parameters$mean.a0)
  ncols.trend = 0

  if ('mean.trend' %in% names(parameters)) {
    ncols.trend = length(parameters$mean.trend)
  }
  nrows = length(data$Qhat.flow);
  ncols.max = max(c(ncols.a0, ncols.trend))

  if (ncols.max > .Object@nStates)
    stop(paste('The number of parameters for each term of the mean model must must equal 1 or the number of states of ',.Object@nStates))

  # Check which terms are uniform for all states and whic terms are unique
  # to each state.
  if (ncols.a0==1 || ncols.a0==.Object@nStates) {
    a0.est = matrix(rep(parameters$mean.a0,each=nrows),nrows,.Object@nStates);
  } else if (ncols.a0<.Object@nStates) {
    stop(paste('The number of parameters for the a0 term of the mean model must must equal 1 or the number of states of ',.Object@nStates))
  }
  if (ncols.trend==1 || ncols.trend==.Object@nStates) {
    trend.est = matrix(rep(parameters$mean.trend,each=nrows),nrows,.Object@nStates);
  } else {
    trend.est = 0
  }

  # time.vals = matrix(data$year - data$year[1],nrows,.Object@nStates)
  # Q.data = matrix(data$Qhat.flow,nrows,.Object@nStates);

  # Calculate the mean log(concentration) over time and with each state
  # Calculate the non-AR1 componants
  # assume a and b vary seasonally...
  # a1.est = a1.amp * sin(a1.phase+time.vals/12) + a1.drop
  # a0.est = a0.amp * sin(a0.phase+time.vals/12) + a0.drop

  # Calculate the estimate of mean c, bounds are -5 to 5 for all parameters
  # L0.est = 10^L0.est
  # V0.est = 10^V0.est
  # Rq.est = 10^Rq.est
  # yq.est = 10^yq.est
  a0.est = a0.est

  # just residuals
  Qhat.model <- a0.est

  # print(paste('...DBG getMean.AR0 nrows Qhat.model.NAs:',nrow(Qhat.model)))

  return(Qhat.model)
}
)

# Get variance. It is calculates as a parameter times the Qhat variance.
# setGeneric(name="getVariance",def=function(.Object, data) {standardGeneric("getVariance")})
setMethod(f="getVariance",signature=c("RhatModel.homo.normal.linear","data.frame"),definition=function(.Object, data)
{
  # Get object parameter list
  parameters = getParameters(.Object@parameters)

  ncols.a0 = length(parameters$std.a0)
  nrows = length(data$Qhat.flow);

  if ('std.a1' %in% names(parameters)){
    ncols.a1 = length(parameters$std.a1)
    nrows = length(data$Qhat.precipitation);
  } else {
    ncols.a1=0
  }

  # Get variance of the Qhat # flow is residuals and precipitation is flow,
  Qhat.var = var(data$Qhat.flow, na.rm=T)

  if (ncols.a1 >0){
    a0.est = Qhat.var * matrix(rep(parameters$std.a0,each=nrows),nrows,.Object@nStates) +
      matrix(rep(parameters$std.a1,each=nrows),nrows,.Object@nStates) * data$Qhat.precipitation;
  }else{
    a0.est = Qhat.var * matrix(rep(parameters$std.a0,each=nrows),nrows,.Object@nStates)
  }

  #precip.data = rep(.Object@input.data$Qhat.precipitation,1,ncols.a0);
  # print(head(a0.est))
  return(a0.est)

}
)


# # Get transition matrix with no input data.
setMethod(f="getEmissionDensity",
          signature=c("RhatModel.homo.normal.linear","data.frame","logical"),
          definition=function(.Object, data, cumProb.threshold.Qhat =NA)
          {
            # print('getEmissionDensity')
            # browse()
            # Check Qhat is in data
            if (!any(names(data)=="Qhat.flow"))
              stop('Input "data" must be a a data frame with a variable named "Qhat.flow".')

            # if (!any(names(data)=="Qhat.precipitation"))
            #   stop('Input "data" must be a a data frame with a variable named "Qhat.precipitation".')

            # Get the moments
            markov.mean = getMean(.Object, data)
            markov.variance = getVariance(.Object, data)
            markov.stds <- sqrt(markov.variance)

            # message(paste('markov.mean:',head(markov.mean),'markov.variance', head(markov.variance),'markov.stds', head(markov.stds)))

            # print('markov.mean')
            # print((markov.mean[9746:9750]))
            # print('markov.variance')
            # print(head(markov.variance))
            # print('markov.stds')
            # print((markov.stds[9746:9750]))

            # message(paste('... number of Markov mean < 0. Number =',sum(markov.mean<0)))

            # comment out for residual analysis because, we can have a negative mean..
            # For truncated flow, the mean  must >=0 else very negative means can arise
            # if (.Object@use.truncated.dist && any(markov.mean<0,na.rm = T) ) {
            #   P = matrix(Inf, nrow(markov.mean), .Object@nStates)
            #   # print('P1')
            #   # print(head(P))
            #
            #   #markov.mean = pmax(0, markov.mean,na.rm = T)
            # }

            # Calculate probabilities.
            P <- matrix(NA, nrow(data),.Object@nStates)
            # print('P2')
            # print(head(P))
            # Set the lower limit of a truncated normal dist.
            lowerX = -Inf;
            if (.Object@use.truncated.dist)
              lowerX = -Inf;

            # index to only run on observed C and Q, ignore gaps
            # indC <- !is.na(data$C) & !is.na(data$flow)

            # get the prob.
            if (!all(is.na(cumProb.threshold.Qhat))) { #are all of the values not equal to na?
              if (length(cumProb.threshold.Qhat)!=nrow(data))
                stop('The length of cumProb.threshold.Qhat must equal the number of rows of input data')

              for (i in 1:.Object@nStates) {
                P[,i] = ptruncnorm(cumProb.threshold.Qhat[i], a=-Inf, mean=markov.mean[,i], sd=markov.stds[,i])
                # print('ptruncnorm P[,i]')
                # print(P[,i])
              }
            } else {
              for (i in 1:.Object@nStates) { # 'a' is the lower limit of the truncated normal distribution, making -Inf for residual analysis
                P[,i] = dtruncnorm(data$Qhat.flow, a=-Inf, mean=markov.mean[,i], sd=markov.stds[,i])
                # print("data$Qhat.C")
                # print(head(data$Qhat.C))
                # print('dtruncnorm P[,i]')
                # print(P[,i][1:111])
                # message(paste("data$Qhat.C[111]= ",data$Qhat.C[111]," lowerX= ",lowerX, "markov.mean[111]= ",markov.mean[111]," markov.stds[111]= ",markov.stds[111]))

              }
            }

            # print("P")
            # print(sum(P))

            return(P)
          }
)

setMethod(f="getDistributionPercentiles",
          signature=c("RhatModel.homo.normal.linear","data.frame","numeric"),
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
              lowerX = -Inf;

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

setMethod(f="generate.sample.Qhat.fromViterbi",signature=c("RhatModel.homo.normal.linear",'data.frame','numeric'),definition=function(.Object, data, viterbi.states)
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
    lowerX = -Inf;

  for (i in 1:nSamples)
    sample.Qhat[i] = rtruncnorm(1, a=lowerX, mean = markov.mean[i,viterbi.states[i]], sd = sqrt(markov.variance[i,viterbi.states[i]]))

  return(sample.Qhat)
}
)

setMethod(f="generate.sample.Qhat",signature="RhatModel.homo.normal.linear",definition=function(.Object, data, nSamples)
{

  if (length(nSamples)!=1 || nSamples<=0)
    stop('The input nSamples must be a single number >=1.')

  # Get distribution for each state
  markov.mean = getMean(.Object, data)
  markov.variance = getVariance(.Object, data)

  # Set the lower limit of a truncated normal dist.
  lowerX = -Inf;
  if (.Object@use.truncated.dist)
    lowerX = -Inf;

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
