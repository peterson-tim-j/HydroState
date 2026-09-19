##' @include abstracts.R parameters.R
## @export
QhatModel <- setClass(
  # Set the name for the class
  "QhatModel",

  package='hydroState',

  # Define the slots
  slots = c(
    data = "data.frame",
    formula.mean = 'formula',
    formula.var = 'formula',
    MAX_AR_TERMS = 'numeric',
    formula.mean.AR = 'data.frame',
    formula.var.AR = 'data.frame',
    error.distribution = 'character',
    use.truncated.dist = 'logical',
    parameters = 'model_parameters',
    nStates = 'numeric'
  ),
  # Set the default values for the slots. (optional)
  prototype=list(
    data = data.frame(year=numeric(), p=numeric(), q=numeric()),
    formula.mean = as.formula('q ~ a0 + a1*p'),
    formula.var = as.formula('e ~ b0'),
    MAX_AR_TERMS = 3,
    formula.mean.AR = data.frame(AR = rep('', 3),
                                 is_used = rep(F, 3),
                                 param = rep('', 3)),
    formula.var.AR = data.frame(AR = rep('', 3),
                                is_used = rep(F, 3),
                                param = rep('', 3)),
    error.distribution = 'trunc.norm',
    use.truncated.dist = T
  )
)


# Initialise object
setMethod("initialize","QhatModel", function(.Object,
                                             data = data.frame(year=numeric(), p=numeric(), q=numeric()),
                                             formula.mean = as.formula('q ~ a0 + a1*p'),
                                             formula.var = as.formula('e ~ b0'),
                                             shift.pars = c('a0', 'b0'),
                                             error.distribution = 'trunc.norm',
                                             transition.graph=matrix(T,2,2)) {

  # Add inputs to object
  .Object@data <- data
  .Object@nStates = ncol(transition.graph)
  .Object@formula.mean <- formula.mean
  .Object@formula.var <- formula.var
  .Object@error.distribution = error.distribution
  .Object@use.truncated.dist = F
  if (error.distribution == 'truncnorm' || error.distribution == 'gamma')
    .Object@use.truncated.dist = T

  # Get the parameter names from all formula and exclude predictors
  pars.mean = all.vars(formula.mean[[3]])     # RHS of formula
  ind = !(pars.mean %in% colnames(data))      # Exclude predictors
  pars.mean = pars.mean[ind]
  ind = !(grepl('.obs', pars.mean))            # Exclude estimates (e.g. estimated q)
  pars.mean = pars.mean[ind]

  pars.err = all.vars(formula.var[[3]])     # RHS of formula
  ind = !(pars.err %in% colnames(data))      # Exclude predictors
  pars.err = pars.err[ind]
  ind = !(grepl('.obs', pars.err))            # Exclude estimates (e.g. estimated q)
  pars.err = pars.err[ind]

  # Get any auto correlation terms in the formula and record the parameter(s) used
  # in the AR terms. This is used by the formula eval to calculate the AR terms.
  f.terms = attr(terms(.Object@formula.mean), "term.labels")
  formula.mean.AR = data.frame(AR = rep('', .Object@MAX_AR_TERMS),
                               is_used = rep(F, .Object@MAX_AR_TERMS),
                               param = rep('', .Object@MAX_AR_TERMS))
  formula.var.AR = .Object@formula.mean.AR
  for (i in 1:.Object@MAX_AR_TERMS) {
    ind = grepl( paste0(':AR',i), f.terms)
    if (any(ind)) {
      formula.mean.AR[i,1] = paste0('AR',i)

      formula.mean.AR[i,2] = T

      # Record parameter multiplying AR term
      formula.mean.AR[i,3] = sub(":.*", "", f.terms[ind])

      # Remove AR term from parameters list
      ind = !(pars.mean %in% paste0('AR',i) )
      if (length(ind)>0)
        pars.mean = pars.mean[ind]
    }
  }
  .Object@formula.mean.AR <- subset(formula.mean.AR, formula.mean.AR$is_used==T)
  f.terms = attr(terms(.Object@formula.var), "term.labels")
  for (i in 1:.Object@MAX_AR_TERMS) {
    ind = grepl( paste0(':AR',i), f.terms)
    if (any(ind)) {
      formula.var.AR[i,1] = paste0('AR',i)

      .Object@formula.var.AR[i,2] = T

      # Record parameter multiplying AR term
      .Object@formula.var.AR[i,3] = sub(":.*", "", f.terms[ind])

      # Remove AR term from parameters list
      ind = !(pars.err %in% paste0('AR',i) )
      if (length(ind)>0)
        pars.err = pars.err[ind]
    }
  }
  .Object@formula.var.AR <- subset(formula.var.AR, formula.var.AR$is_used==T)

  # Get unique params and remove pi.
  pars = unique(c(pars.mean, pars.err))
  ind = !('pi' %in% pars)
  pars = pars[ind]

  # Build vector of number of states in each parameter. Used below to create the parameters obj
  pars.nstates =  rep(1, length(pars))
  names(pars.nstates) = pars
  if (length(shift.pars)>0)
    pars.nstates[shift.pars] = .Object@nStates

  # Build obj for parameters
  .Object@parameters = new('model_parameters', pars, pars.nstates)

  .Object
}
)

# Get start and end indices for independent variable
# @exportMethod getStartEndIndex
setGeneric(name = "getStartEndIndex", def = function(data) {standardGeneric("getStartEndIndex")})
setMethod(f = "getStartEndIndex", signature = c("data.frame"), function(data){


  end.index = which(diff(is.finite(data$precipitation)) == -1)

  #if no gaps, just return, first obs. precip to last obs. precip
  if(length(end.index) == 0){
    return(data.frame(start.index = checkmate::wf(is.finite(data$precipitation),TRUE),
                      end.index = checkmate::wl(is.finite(data$precipitation), TRUE)))
  }

  end.index[NROW(end.index) + 1] = checkmate::wl(is.finite(data$precipitation), TRUE)

  start.index = rep(NA,length(end.index))

  if(checkmate::wf(is.finite(data$precipitation),TRUE) ==1){ #if starts at first time-step
    start.index[1] = checkmate::wf(is.finite(data$precipitation),TRUE)
    start.index[2:NROW(end.index)] = which(diff(is.finite(data$precipitation)) == 1) +1
  }else{ # if not,
    start.index[1:NROW(end.index)] = which(diff(is.finite(data$precipitation)) == 1) +1
  }

  delta <- data.frame(start.index,end.index)

  return(delta)

}
)

# Get transition matrix with no input data.
setGeneric(name="getEmissionDensity",def=function(.Object, data, zero.Flow, cumProb.threshold.Qhat) {standardGeneric("getEmissionDensity")})
setMethod(f="getEmissionDensity",signature=c("QhatModel","data.frame", 'numeric', 'missing'),
definition=function(.Object, data, zero.Flow, cumProb.threshold.Qhat)
{
  return( getEmissionDensity(.Object, data, zero.Flow, -999) )
}
)
setMethod(f="getEmissionDensity",signature=c("QhatModel","data.frame", 'numeric', 'numeric'),
definition=function(.Object, data, zero.Flow, cumProb.threshold.Qhat)
{
  # Get obs dependent data
  dep.name <- all.vars(.Object@formula.mean[[2]])
  dep.data <- data[, dep.name]

  # Get the moments
  mu <- state.mean(.Object, data)
  sigma2 <- state.variance(.Object, data, mu)
  sigma <- sqrt(sigma2)

  # Calculate probabilities.
  P <- matrix(NA, nrow(data),.Object@nStates)

  # Set the lower limit of a truncated normal dist.
  lowerX = -Inf;
  if (.Object@use.truncated.dist)
    lowerX = zero.Flow;

  # Check if getting cumulative probs or prob density.
  get_probDensity = T
  if (!missing(cumProb.threshold.Qhat) && sum(cumProb.threshold.Qhat) != -999) {
    if (length(cumProb.threshold.Qhat)!=nrow(data))
      stop('The length of cumProb.threshold.Qhat must equal the number of rows of input data')

    get_probDensity = F
  }

  # Get the prob. for user defined parametric distribution.
  if (.Object@error.distribution == 'norm') {
    if (get_probDensity)
      for (i in 1:.Object@nStates)
        P[,i] = dnorm(dep.data, mean=mu[,i], sd=sigma[,i])
    else
      for (i in 1:.Object@nStates)
        P[,i] = pnorm(cumProb.threshold.Qhat[i], mean=mu[,i], sd=sigma[,i])

  } else if (.Object@error.distribution == 'trunc.norm') {
    if (get_probDensity)
      for (i in 1:.Object@nStates)
        P[,i] = dtruncnorm(dep.data, a=lowerX, mean=mu[,i], sd=sigma[,i])
    else
      for (i in 1:.Object@nStates)
        P[,i] = ptruncnorm(cumProb.threshold.Qhat[i], a=lowerX, mean=mu[,i], sd=sigma[,i])

  } else if (.Object@error.distribution == 'gamma') {
      # Limit mean to >0
      dep.data[dep.data==0] = sqrt(.Machine$double.eps)
      filt = is.finite(mu) & mu <= sqrt(.Machine$double.eps)*1000
      mu[filt] = sqrt(.Machine$double.eps)*1000

      # Derive the gamma parameters from the modelled mean and variance.
      markov.shape <- mu^2/sigma2
      markov.scale <- sigma2/mu

      if (get_probDensity)
        for (i in 1:.Object@nStates)
          P[,i] <- suppressWarnings(dgamma(dep.data, shape=markov.shape[,i], scale=markov.scale[,i]))#},warning = function(s){message('AR parameters not stationary')})
      else
        for (i in 1:.Object@nStates)
          for (j in 1:nrow(data)) {
            if (is.na(cumProb.threshold.Qhat[j]))
              P[j,i] <- NA
            else
              P[j,i] <- pgamma(cumProb.threshold.Qhat[j], shape=markov.shape[j,i], scale=markov.scale[j,i])
        }
  }

  return(P)
}
)

# @exportMethod getDistributionPercentiles
setGeneric(name="getDistributionPercentiles",def=function(.Object, data, zero.Flow, precentiles) {standardGeneric("getDistributionPercentiles")})
setMethod(f="getDistributionPercentiles", signature=c("QhatModel","data.frame","numeric","numeric"),
definition=function(.Object, data, zero.Flow, precentiles=c(0.5, 0.95))
{
  # Get moments over time
  mu = state.mean(.Object, data)
  sigma2 = state.variance(.Object, data, mu)
  sigma <- sqrt(sigma2)

  # Initialise returned variable
  est = vector('list',length(precentiles))

  # Set the lower limit of a truncated normal dist.
  lowerX = -Inf;
  if (.Object@use.truncated.dist)
    lowerX = zero.Flow

  # Calculate probabilities.
  if (.Object@error.distribution == 'norm') {
    for (i in 1:length(precentiles)) {
      est.tmp = matrix(0,nrow(mu),ncol(mu))
      for (j in 1:.Object@nStates) {
        est.tmp[,j] <- qnorm(precentiles[i], mean= mu[,j], sd=sigma[,j])
      }
      est[[i]] = est.tmp
    }
  } else if (.Object@error.distribution == 'trunc.norm') {
    for (i in 1:length(precentiles)) {
      est.tmp = matrix(0,nrow(mu),ncol(mu))
      for (j in 1:.Object@nStates)
        est.tmp[,j] <- qtruncnorm(precentiles[i], a=lowerX, mean= mu[,j], sd=sigma[,j])

      est[[i]] = est.tmp
    }
  } else if (.Object@error.distribution == 'gamma') {
    # Limit mean to >0
    filt = is.finite(mu) & mu <= sqrt(.Machine$double.eps)*1000
    mu[filt] = sqrt(.Machine$double.eps)*1000

    # Derive the gamma parameters from the modelled mean and variance.
    markov.shape <- mu^2/sigma2
    markov.scale <- sigma2/mu

    est.tmp = matrix(0,nrow(mu),ncol(mu))

    for (i in 1:length(precentiles)) {
      est.tmp = matrix(0,nrow(mu),ncol(mu))
      for (j in 1:.Object@nStates) {
        for (k in 1:nrow(data)) {
          if (is.na(mu[k,j]))
            est.tmp[k,j] = NA
          else
            est.tmp[k,j] = qgamma(precentiles[i], shape=markov.shape[k,j], scale=markov.scale[k,j])
        }
      }
      est[[i]] = est.tmp
    }
  }

  names(est) = precentiles
  return(est)
}
)


# Calculate the mean predictor for each state over time
setGeneric(name="state.mean",def=function(.Object, data) {standardGeneric("state.mean")})
setMethod(f="state.mean",signature=c("QhatModel", "data.frame"), definition=function(.Object, data)
  {
    # Check predictants and predictors columns in mean formula are in the columns of data.
    f.terms = all.vars(.Object@formula.mean)
    param.names = names(parameters(.Object@parameters))
    f.terms = f.terms[!(f.terms %in% param.names)]
    ind = !(grepl('.obs', f.terms))            # Exclude estimates (e.g. estimated q)
    f.terms = f.terms[ind]
    if (length(f.terms)>0 && !any(f.terms %in% names(data)))
      stop('Input "data" must be a a data frame with columns for each input and output of the mean formula.')

    # Get obs dependent data
    dep.name = all.vars(.Object@formula.mean[[2]])
    dep.data = data[, dep.name]

    # Remove obs dependent data from data passed to formula. This is done
    # to ensure it is not used in the RHS of the formulas.
    data[, dep.name] <- NULL

    # Evaluate formula
    mu = eval.formula(.Object,
                      .Object@formula.mean,
                      .Object@formula.mean.AR,
                      data,
                      dep.name,
                      dep.data)
    return(mu)
  }
)

# Calculate the variance predictor for each state over time
setGeneric(name="state.variance",def=function(.Object, data, mu) {standardGeneric("state.variance")})
setMethod(f="state.variance",signature=c("QhatModel", "data.frame", 'matrix'), definition=function(.Object, data, mu)
  {
    # Check predictants and predictors columns in variance formula are in the columns of data
    f.terms = all.vars(.Object@formula.var[[3]])
    param.names = names(parameters(.Object@parameters))
    f.terms = f.terms[!(f.terms %in% param.names)]
    ind = !(grepl('.obs', f.terms))            # Exclude estimates (e.g. estimated q)
    f.terms = f.terms[ind]
    if (length(f.terms)>0 && !any(f.terms %in% names(data)))
      stop('Input "data" must be a a data frame with columns for each input and output of the variance formula.')

    # Get obs dependent data
    dep.name = all.vars(.Object@formula.mean[[2]])
    dep.data = data[, dep.name]

    # Remove obs dependent data from data passed to formula. This is done
    # to ensure it is not used in the RHS of the formulas.
    data[, dep.name] <- NULL

    # Get the moments
    sigma2 = eval.formula(.Object, .Object@formula.var,
                          .Object@formula.var.AR,
                          data,
                          dep.name,
                          dep.data,
                          mu )

    return(sigma2)
  }
)

# Calculate the transformed flow at the mean annual precip
setGeneric(name="eval.formula",def=function(.Object, f, f.AR_terms, data.indep, pred.name, pred.obs, pred.est) {standardGeneric("eval.formula")})
setMethod(f="eval.formula",signature=c("QhatModel", 'formula', "data.frame", "data.frame"),
          definition=function(.Object, f, f.AR_terms, data.indep, pred.name, pred.obs, pred.est)
{
  # Define base statistical functions that remove NAs
  var <- function(x, ...) stats::var(x, na.rm = TRUE, ...)
  IQR <- function(x, ...) stats::IQR(x, na.rm = TRUE, ...)
  mean <- function(x, ...) stats::mean(x, na.rm = TRUE, ...)
  median <- function(x, ...) stats::median(x, na.rm = TRUE, ...)

  # Get matrix of object parameter vals with each state in a different column.
  params.df = parameters(.Object@parameters,  transformed = F, as_df =T)

  # Evaluate formula for each state, but omit the auto correlation terms.
  # The first if() adds the predicted mean dependent variable for state i to
  # the inputs enviro. The second if() added the observed dependent variable.
  nsteps = nrow(data.indep)
  vals = matrix(NA, nsteps, .Object@nStates)
  for (i in 1:.Object@nStates) {
    ind = params.df[ ,'state'] == 0 | params.df[ , 'state'] == i
    inputs.bystate = as.list( t(params.df[ind ,3, drop=F]) )
    names(inputs.bystate) = params.df[ind, 1]

    # Add each independent variable
    inputs.bystate = c( inputs.bystate, as.vector(data.indep) )
    if (!missing(pred.est) && ncol(pred.est)==.Object@nStates && nrow(pred.est)==nrow(data.indep)) {
      inputs.bystate = c( inputs.bystate, list(mu = pred.est[,i]))
      names(inputs.bystate)[length(inputs.bystate)] <- pred.name
    }
    # Add observed data - if used in RHS of formula. This allows variance
    # of independent variable to be used to scale the variance formula
    if (!missing(pred.obs) && length(pred.obs)==nrow(data.indep)) {
      inputs.bystate = c( inputs.bystate, list(predictant = pred.obs))
      names(inputs.bystate)[length(inputs.bystate)] <- paste0(pred.name,'.obs')
    }

    # Set AR inputs to zero. The AR terms are calculated below.
    if (nrow(f.AR_terms)>0) {
      inputs.AR = as.list( rep(0, nrow(f.AR_terms)))
      names(inputs.AR) <- f.AR_terms[,'AR']
      inputs.bystate = c( inputs.bystate, inputs.AR)
    }
    vals[, i] = eval(f[[3]], envir = inputs.bystate)
  }

  # Loop through all time steps and add auto-regressive terms given in the formula
  if (nrow(f.AR_terms)>0) {
    has.AR1 = F
    has.AR2 = F
    has.AR3 = F
    # Find the index to each AR parameter
    if ('AR1' %in% f.AR_terms$AR) {
      ind = f.AR_terms$AR %in% 'AR1'
      AR.param.name = f.AR_terms[ind,'param']
      ind.AR1 = params.df$parameter %in% AR.param.name
      has.AR1 = T
    }
    if ('AR2' %in% f.AR_terms$AR) {
      ind = f.AR_terms$AR %in% 'AR2'
      AR.param.name = f.AR_terms[ind,'param']
      ind.AR2 = params.df$parameter %in% AR.param.name
      has.AR2 = T
    }
    if ('AR3' %in% f.AR_terms$AR) {
      ind = f.AR_terms$AR %in% 'AR3'
      AR.param.name = f.AR_terms[ind,'param']
      ind.AR3 = params.df$parameter %in% AR.param.name
      has.AR3 = T
    }

    # Add AR 1 to 3 terms to formula evaluation.
    vals.mean = colMeans(vals, na.rm = FALSE)
    for (i in 2:nrow(vals)) {
      if (has.AR1 && (i<3 || (!has.AR2 && !has.AR3)))
        vals[i, ] = vals[i, ] + params.df[ind.AR1,'value'] * ( vals[i-1, ] - vals.mean)
      else if (has.AR1 && has.AR2 && (i<4 || !has.AR3))
        vals[i, ] = vals[i, ] + params.df[ind.AR1,'value'] * ( vals[i-1, ] - vals.mean) +
                                params.df[ind.AR2,'value'] * ( vals[i-2, ] - vals.mean)
      else if (has.AR1 && has.AR2 && has.AR3 && i<5)
        vals[i, ] = vals[i, ] + params.df[ind.AR1,'value'] * ( vals[i-1, ] - vals.mean) +
                                params.df[ind.AR2,'value'] * ( vals[i-2, ] - vals.mean) +
                                params.df[ind.AR3,'value'] * ( vals[i-3, ] - vals.mean)
    }
  }
  return(vals)
  }
)

setGeneric(name="generate.sample.Qhat.fromViterbi",def=function(.Object, data, viterbi.states) {standardGeneric("generate.sample.Qhat.fromViterbi")})
setMethod(f="generate.sample.Qhat.fromViterbi",signature=c("QhatModel",'data.frame','numeric'),definition=function(.Object, data, viterbi.states)
{

  if (length(viterbi.states) != nrow(data))
    stop('The length of the viterbi.states must equal the number of rows in data.')

  # Generate a synthtic series of Qhat usng the input random series of states
  nSamples = length(viterbi.states)
  sample.Qhat = rep(0,nSamples)

  # Get moments over time
  mu = state.mean(.Object, data)
  sigma2 = state.variance(.Object, data, mu)

  if (nrow(mu)!=nSamples)
    stop('The input vector of sample states must equal the numer of rows of the input data.')

  # Set the lower limit of a truncated normal dist.
  lowerX = -Inf;
  if (.Object@use.truncated.dist)
    lowerX = 0;

  for (i in 1:nSamples)
    sample.Qhat[i] = rtruncnorm(1, a=lowerX, mean = mu[i,viterbi.states[i]], sd = sqrt(sigma2[i,viterbi.states[i]]))

  return(sample.Qhat)
}
)

setGeneric(name="generate.sample.Qhat",def=function(.Object, data, nSamples) {standardGeneric("generate.sample.Qhat")})
setMethod(f="generate.sample.Qhat",signature="QhatModel",definition=function(.Object, data, nSamples)
{

  if (length(nSamples)!=1 || nSamples<=0)
    stop('The input nSamples must be a single number >=1.')

  # Get moments over time
  mu = state.mean(.Object, data)
  sigma2 = state.variance(.Object, data, mu)

  # Set the lower limit of a truncated normal dist.
  lowerX = -Inf;
  if (.Object@use.truncated.dist)
    lowerX = 0;

  sample.Qhat = vector('list',.Object@nStates)
  for (i in 1:.Object@nStates) {
    sample.Qhat[[i]] = matrix(NA,nrow(data),nSamples)
    for (j in 1:nrow(data)) {
      sample.Qhat[[i]][j,] = rtruncnorm(nSamples, a=lowerX, mean = mu[j,i], sd = sqrt(sigma2[j,i]))
    }
  }

  return(sample.Qhat)
}
)

