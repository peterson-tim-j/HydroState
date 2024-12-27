##' @include Qhat.boxcox.R
##' @include QhatModel.homo.normal.linear.R
##' ##' @include QhatModel.homo.normal.linear.AR1..R
##' @include markov.annualHomogeneous.R
## @export
hydroState <- setClass(
  # Set the name for the class
  "hydroState",

  package='hydroState',

  # Define the slots
  slots = c(
    input.data = "data.frame",
    Qhat.object   = "Qhat",
    QhatModel.object   = "QhatModel",
    markov.model.object = "markov",
    calibration.results = "list",
    state.labels = "character"
  ),

  # Set the default values for the slots. (optional)
  prototype=list(
    input.data = data.frame(year=c(0),month=c(0),precipitation=c(0),flow=c(0)),
    Qhat.object  = new('Qhat.boxcox',input.data = data.frame(year=c(0),month=c(0),precipitation=c(0),flow=c(0))),
    QhatModel.object = new('QhatModel.homo.normal.linear',input.data = data.frame(year=c(0),month=c(0),precipitation=c(0),flow=c(0))),
    markov.model.object = new('markov.annualHomogeneous',transition.graph = matrix(TRUE,2,2)),
    calibration.results = list(optim=0,member=0),
    state.labels = c('')
  )
)

# Valid object?
setGeneric(name="validObject",def=function(.Object) {standardGeneric("validObject")})
setMethod(f="validObject",
          signature="hydroState",
          definition=function(.Object)
          {

          Tprob <- getTransitionProbabilities(.Object)
          P0 <- getInitialStateProbabilities(.Object)


          is.valid <- T
          if (any(Tprob>1) || any(Tprob<0)) {
            is.valid <- F
            warning('The tranistion probabilities ARE NOT between zero and one.')
          }

          if (getNumStates(.Object@markov.model.object)>1) {
            if (any( abs(rowSums(Tprob)-1)>sqrt(.Machine$double.eps))) {
              is.valid <- F
              warning('The tranistion probabilities for at least one state DO NOT sum to 1.0.')
            }
          }

          if (abs(sum(P0)-1)>sqrt(.Machine$double.eps)) {
            is.valid <- F
            warning('The initial probabilities DO NOT sum to 1.0.')
          }

          #Check the number of states in the QHat model equal that in the markov object
          if (getNumStates(.Object@markov.model.object) != .Object@QhatModel.object@nStates)
            stop('The number of states in the markov.model input differs from that within the QhatModel.object input.')
  }
)

# Initialise the object.
#setGeneric(name="initialize",def=function(.Object,input.data, Qhat.object, QhatModel.object, markov.model.object, ...){standardGeneric("initialize")})
setMethod(f="initialize",signature="hydroState",definition=function(.Object, input.data, Qhat.object, QhatModel.object, markov.model.object)
          {
          # Check the input data has columns 'precipitation' and 'flow' and the 'recipitation has no gaps.
          if (!is.data.frame(input.data))
            stop('The input input.data must be a data frame.')
          if (!any(names(input.data)=='precipitation'))
            stop('The input input.data must contain the column "precipitation".')
          if (!any(names(input.data)=='flow'))
            stop('The input input.data must contain the column "flow".')
          if (!any(names(input.data)=='year'))
            stop('The input input.data must contain the column "year".')

          # if (length(diff(which(!is.finite(input.data$precipitation)))) >= 1)
          #   message(paste('The independent varaible (precip.) contains gaps at ', sum(!is.finite(input.data$precipitation)), ' timesteps.',' Model built ignoring gaps.',sep=""))

          if (all(!is.finite(input.data$flow)))
            stop('The input input.data$flow does not contain any finite values".')

          if (max(diff(unique(input.data$year)), na.rm = TRUE) !=1){
            message(paste('There are ',sum(ifelse(diff(na.omit(unique(input.data$year))) !=1, diff(na.omit(unique(input.data$year))), 0)), ' years missing.',' Model built ignoring missing years.', sep=""))
          }

          if (any(!is.numeric(input.data$flow)))
            stop('The input input.data$flow must contain only numeric data.')
          if (any(!is.numeric(input.data$precipitation)))
            stop('The input input.data$precipitation must contain only numeric data.')

          .Object@input.data = input.data
          .Object@Qhat.object = Qhat.object
          .Object@QhatModel.object = QhatModel.object
          .Object@markov.model.object = markov.model.object
          .Object
          }
)

# Get the number of states
setGeneric(name="getNumMarkovStates",def=function(.Object) {standardGeneric("getNumMarkovStates")})
setMethod(f="getNumMarkovStates",
          signature="hydroState",
          definition=function(.Object)
          {
            return(nrow(.Object@Markov.transition.graph))
          }
)

# Get the full set of model parameters as a list.
#setGeneric(name="getParameters",def=function(.Object) {standardGeneric("getParameters")})
setMethod(f="getParameters",
          signature="hydroState",
          definition=function(.Object)
          {
            return(list(
              Qhat = getParameters(.Object@Qhat.object@parameters),
              QhatModel = getParameters(.Object@QhatModel.object@parameters),
              markov = getParameters(.Object@markov.model.object@parameters)
            ))
          }
)

# Get the full set of model parameters as a list.
# @exportMethod getParameters.asVector
setGeneric(name="getParameters.asVector",def=function(.Object) {standardGeneric("getParameters.asVector")})
setMethod(f="getParameters.asVector",
          signature="hydroState",
          definition=function(.Object)
          {
            parameters.all = getParameters(.Object);

            return( c(
              unlist(parameters.all$Qhat.flow),
              unlist(parameters.all$QhatModel),
              unlist(parameters.all$markov)
              ))
          }
)

# #' @exportMethod getParameters.transformed
#setGeneric(name="getParameters.transformed",def=function(.Object) {standardGeneric("getParameters.transformed")})
setMethod(f="getParameters.transformed",
          signature="hydroState",
          definition=function(.Object)
          {
            Qhat = getParameters.transformed(.Object@Qhat.object@parameters)
            QhatModel = getParameters.transformed(.Object@QhatModel.object@parameters)
            markov = getParameters.transformed(.Object@markov.model.object@parameters)

            return(list(
              Qhat = Qhat,
              QhatModel = QhatModel,
              markov = markov
            ))
          }
)

# @exportMethod getParameters.transformed.asVector
setGeneric(name="getParameters.transformed.asVector",def=function(.Object) {standardGeneric("getParameters.transformed.asVector")})
setMethod(f="getParameters.transformed.asVector",
          signature="hydroState",
          definition=function(.Object)
          {
            parameters.all = getParameters.transformed(.Object);

            return( c(
              unlist(parameters.all$Qhat),
              unlist(parameters.all$QhatModel),
              unlist(parameters.all$markov)
            ))
          }
)

# @exportMethod getBounds.transformed.asVector
setGeneric(name="getBounds.transformed.asVector",def=function(.Object) {standardGeneric("getBounds.transformed.asVector")})
setMethod(f="getBounds.transformed.asVector",
          signature="hydroState",
          definition=function(.Object)
          {
            # Get the parameter bounds for each object
            Qhat = getBounds.transformed(.Object@Qhat.object@parameters)
            QhatModel = getBounds.transformed(.Object@QhatModel.object@parameters)
            markov = getBounds.transformed(.Object@markov.model.object@parameters)

            lowerBound = c( unlist(Qhat$lower),
                            unlist(QhatModel$lower),
                            unlist(markov$lower))

            upperBound = c( unlist(Qhat$upper),
                            unlist(QhatModel$upper),
                            unlist(markov$upper))
            return( cbind(lowerBound, upperBound))
          }
)

# Get the full set of model parameters as a list.
#setGeneric(name="setParameters",def=function(.Object,parameters) {standardGeneric("setParameters")})
setMethod(f="setParameters",
          signature=c("hydroState","list"),
          definition=function(.Object,parameters)
          {

            if (is.list(parameters)) {
              .Object@Qhat.object@parameters <- setParameters(.Object@Qhat.object@parameters, parameters$Qhat)
              .Object@QhatModel.object@parameters <- setParameters(.Object@QhatModel.object@parameters, parameters$QhatModel)
              .Object@markov.model.object@parameters <- setParameters(.Object@markov.model.object@parameters, parameters$markov)
            # } else if(is.numeric(parameters)) {
            #   .Object@Qhat.object@parameters <- setParameters.fromVector(.Object@Qhat.object@parameters, parameters$Qhat)
            #   .Object@QhatModel.object@parameters <- setParameters.fromVector(.Object@QhatModel.object@parameters, parameters$QhatModel)
            #   .Object@markov.model.object@parameters <- setParameters.fromVector(.Object@markov.model.object@parameters, parameters$markov)

            } else{
              message("setParameters error")
            }
            return(.Object)
          }
)
# Get the full set of model parameters as a list.
#setGeneric(name="setParameters",def=function(.Object,parameters) {standardGeneric("setParameters")})
setMethod(f="setParameters",
          signature=c("hydroState","numeric"),
          definition=function(.Object,parameters)
          {
            # Get the parameter structure for each object
            Qhat = getParameters(.Object@Qhat.object@parameters)
            QhatModel = getParameters(.Object@QhatModel.object@parameters)
            markov = getParameters(.Object@markov.model.object@parameters)

            # Calc. number parameters per object
            nQhat <- sum(lengths(Qhat))
            nQhatModel <- sum(lengths(QhatModel))
            nMarkov <- sum(lengths(markov))

            if (nQhat>0) {
              from =1;
              to = nQhat;
              parameters.aslist = relist(parameters[from:to],Qhat)
              .Object@Qhat.object@parameters < -setParameters(.Object@Qhat.object@parameters, parameters.aslist)
            }


            if (nQhatModel>0)
            {
              from =1+nQhat;
              to = nQhat+nQhatModel;
              parameters.aslist = relist(parameters[from:to],QhatModel)
              .Object@QhatModel.object@parameters <- setParameters(.Object@QhatModel.object@parameters, parameters.aslist)
            }

            if (nMarkov>0) {
              from =1+nQhat+nQhatModel;
              to = nQhat+nQhatModel+nMarkov;
              parameters.aslist = relist(parameters[from:to],markov)
              .Object@markov.model.object@parameters <- setParameters(.Object@markov.model.object@parameters, parameters.aslist)
            }

            return(.Object)

          }
)

# @exportMethod setParameters.fromTransformed.asVector
setGeneric(name="setParameters.fromTransformed.asVector",def=function(.Object,parameters.asVector) {standardGeneric("setParameters.fromTransformed.asVector")})
setMethod(f="setParameters.fromTransformed.asVector",
          signature=c("hydroState","numeric"),
          definition=function(.Object,parameters.asVector)
          {
            # Get the parameter structure for each object
            Qhat = getParameters(.Object@Qhat.object@parameters)
            QhatModel = getParameters(.Object@QhatModel.object@parameters)
            markov = getParameters(.Object@markov.model.object@parameters)

            # Calc. number parameters per object
            nQhat <- sum(lengths(Qhat))
            nQhatModel <- sum(lengths(QhatModel))
            nMarkov <- sum(lengths(markov))

            if (nQhat>0) {
              from =1;
              to = nQhat;
              parameters.aslist = relist(parameters.asVector[from:to],Qhat)
              .Object@Qhat.object@parameters <- setParameters.fromTransformed(.Object@Qhat.object@parameters, parameters.aslist)
            }

            if (nQhatModel>0)
            {
              from =1+nQhat;
              to = nQhat+nQhatModel;
              parameters.aslist = relist(parameters.asVector[from:to],QhatModel)
              .Object@QhatModel.object@parameters <- setParameters.fromTransformed(.Object@QhatModel.object@parameters, parameters.aslist)
            }

            if (nMarkov>0) {
              from =1+nQhat+nQhatModel;
              to = nQhat+nQhatModel+nMarkov;
              parameters.aslist = relist(parameters.asVector[from:to],markov)
              .Object@markov.model.object@parameters <- setParameters.fromTransformed(.Object@markov.model.object@parameters, parameters.aslist)
            }

            return(.Object)
          }
)

# Get a vector of initial state probabilies
setMethod(f="getInitialStateProbabilities",signature="hydroState",  definition=function(.Object){})
setMethod(f="getInitialStateProbabilities",
          signature="hydroState",
          definition=function(.Object)
          {
            return(getInitialStateProbabilities(.Object@markov.model.object))
          }
)

# Get a vector of transiton probs
setMethod(f="getTransitionProbabilities",signature="hydroState",  definition=function(.Object){})
setMethod(f="getTransitionProbabilities",
          signature="hydroState",
          definition=function(.Object)
          {
            return(getTransitionProbabilities(.Object@markov.model.object))
          }
)

# Get start and end indices for independent variable

# @exportMethod getStartEndIndex
setGeneric(name = "getStartEndIndex", def = function(data,...) {standardGeneric("getStartEndIndex")})
setMethod(f="getStartEndIndex", signature = c(data = "data.frame"), definition = function(data){

          # # Get the transformed flow, Qhat
          # data = .Object@input.data
          if('day' %in% colnames(data)){

            end.index = which(diff(as.Date(lubridate::make_datetime(data$year, data$month, data$day, tz = "Australia/Melbourne"))) != 1)

            if(length(end.index) > 1){ # if more than one continous period
              end.index[NROW(end.index) + 1] = sum(is.finite(data$precipitation))
            } else {
              end.index = sum(is.finite(data$precipitation))
            }

            start.index = rep(NA,length(end.index))
            start.index[1] = 1

            if(length(start.index) > 1){
              start.index[2:NROW(end.index)] = end.index[1:(NROW(end.index)-1)] + 1
            }

          }else if('month' %in% colnames(data)){

            # check if data is monthly or seasonal by finding the longest sequence of years (12 = monthly, 4 = seasonal)
            if(with(rle(data$year), max(lengths)) > 4){
              end.index = which(diff(data$month) != 1 & diff(data$month) != -11)
            }else{
              end.index = which(diff(data$month) != 3 & diff(data$month) != -9)
            }

            if(length(end.index) > 1){ # if more than one continuous period
              end.index[NROW(end.index) + 1] = sum(is.finite(data$precipitation))
            } else {
              end.index = sum(is.finite(data$precipitation))
            }


            start.index = rep(NA,length(end.index))
            start.index[1] = 1

            if(length(start.index) > 1){
              start.index[2:NROW(end.index)] = end.index[1:(NROW(end.index)-1)] + 1
            }

          }else{#('year' %in% colnames(data)){ assume annual

            end.index = which(diff(data$year != 1) & diff(data$year) != 0)

            if(length(end.index) > 1){ # if more than one continous period
              end.index[NROW(end.index) + 1] = sum(is.finite(data$precipitation))
            } else {
              end.index = sum(is.finite(data$precipitation))
            }

            start.index = rep(NA,length(end.index))
            start.index[1] = 1

            if(length(start.index) > 1){
              start.index[2:NROW(end.index)] = end.index[1:(NROW(end.index)-1)] + 1
            }
          }

          # print(delta)
          delta <- cbind(start.index,end.index)

          #re-compute index, if observations are less than minimum_period
          # if(any(delta[,2] - delta[,1] < minimum_period)){
          #   delta <- delta[-which(delta[,2] - delta[,1] < minimum_period),]
          # }

          return(delta)

}
)


# Get the model negative log liklihood.
# @exportMethod getNegLogLikelihood
setGeneric(name="getNegLogLikelihood",def=function(.Object, parameters, ...) {standardGeneric("getNegLogLikelihood")})
setMethod(f="getNegLogLikelihood",signature=c(.Object="hydroState",parameters="list"),definition=function(.Object,parameters)
          {
            # Set the parameters.
           .Object <- setParameters(.Object, parameters)

            return(getNegLogLikelihood(.Object))
          }
          )
# @exportMethod getNegLogLikelihood
#setGeneric(name="getNegLogLikelihood",def=function(.Object) {standardGeneric("getNegLogLikelihood")})
setMethod(f="getNegLogLikelihood",signature=c(.Object="hydroState",parameters='missing'),definition=function(.Object)
          {

            # Get the transformed flow, Qhat
            data = getQhat(.Object@Qhat.object, .Object@input.data)

            # Add Qhat to the input data
            #data = cbind.data.frame(.Object@input.data,Qhat=Qhat)

            # Set-up to run for all periods with continuous observations of independent variable (precipitation)
            delta = getStartEndIndex(data)

            # Get the markov likelihood and return. Importantly, the object QhatBar is passed so that
            # the markov object can get the model estimates of the transformed flow mean, standard deviation etc.
            if(NROW(delta)>1){

              # Get the probabiity of the observed Qhat for each state at each time point.
              emission.probs = lapply(1:NROW(delta), function(i) getEmissionDensity(.Object@QhatModel.object, data[delta[i,1]:delta[i,2],], NA))

                if (all(is.na(unlist(emission.probs))) || max(unlist(emission.probs), na.rm=T)==0) {
                    return(Inf)
                  }

                nll <- lapply(1:NROW(delta), function(i) getLogLikelihood(.Object@markov.model.object, data[delta[i,1]:delta[i,2],], emission.probs[delta[i,1]:delta[i,2],]))
                nll <- sum(unlist(nll))
              }else{

                emission.probs = getEmissionDensity(.Object@QhatModel.object, data, NA)

                if (all(is.na((emission.probs))) || max((emission.probs), na.rm=T)==0) {
                  return(Inf)
                }

                nll <- getLogLikelihood(.Object@markov.model.object, data, emission.probs)
              }

            if (!is.finite(nll)) {
              return(Inf)
            } else {
              return( -nll)
            }


          }
)

# @exportMethod getNegLogLikelihood.fromTransformedVector
setGeneric(name="getNegLogLikelihood.fromTransformedVector",def=function(parameters,.Object) {standardGeneric("getNegLogLikelihood.fromTransformedVector")})
setMethod(f="getNegLogLikelihood.fromTransformedVector",signature=c(parameters="numeric",.Object="hydroState"),definition=function(parameters,.Object)
          {
            # Set the parameters.
            .Object <- setParameters.fromTransformed.asVector(.Object, parameters)
            nll = getNegLogLikelihood(.Object)
            return(nll)
          }
)

# @exportMethod getAIC
setGeneric(name="getAIC",def=function(.Object, ...) {standardGeneric("getAIC")})
setMethod(f="getAIC",signature="hydroState",definition=function(.Object)
{
  # Set neg. log liklihood
  nll <- getNegLogLikelihood(.Object)

  # Get the number of parameters
  np <- length(getParameters.asVector(.Object))

  return(2*(nll+np))
}
)

# @exportMethod getAICc
setGeneric(name="getAICc",def=function(.Object, ...) {standardGeneric("getAICc")})
setMethod(f="getAICc",signature="hydroState",definition=function(.Object)
{
  # Set neg. log liklihood
  nll <- getNegLogLikelihood(.Object)

  # Get the number of parameters
  np <- length(getParameters.asVector(.Object))

  # Calculate small sample size correction.
  n <- nrow(getQhat(.Object@Qhat.object, .Object@input.data))
  corr <- 2*np*(np + 1)/(n - np - 1)

  return(2*(nll+np+corr))
}
)


# @exportMethod fit
setGeneric(name="fit",def=function(.Object,
                                   DEstrategy=NA,
                                   pop.size.perParameter=NA,
                                   max.generations=NA,
                                   Domains=NA,
                                   reltol=NA,
                                   steptol=NA,
                                   print.iterations=NA,
                                   use.initial.parameters=NA,
                                   doParallel=NA,
                                   ...) {standardGeneric("fit")})
setMethod(f = "fit",signature="hydroState",definition=function(.Object,
                                                               DEstrategy=3,
                                                               pop.size.perParameter=25,
                                                               max.generations=10000,
                                                               Domains,
                                                               reltol=1e-8,
                                                               steptol=50,
                                                               print.iterations = 25,
                                                               use.initial.parameters=F,
                                                               doParallel = F,
                                                               ...)
{


  # Get initial parameters and obj valu
  par.initial <- getParameters.transformed.asVector(.Object)
  obj.initial <- getNegLogLikelihood.fromTransformedVector(par.initial, .Object)

  # Get parameter boounds in transformed space.
  if (is.na(Domains))
    Domains <- getBounds.transformed.asVector(.Object)

  # No. parameters
  nvars <- nrow(Domains)

  # Set population size
  NP = nvars*pop.size.perParameter

  # Create and initial population WITH the existing model parameters.
  # DEstrategy <- 3
  # DEstrategy <- 2
  Fweight=1.5
  if (use.initial.parameters){
    par.initial = as.matrix(t(par.initial),1,nvars)
    for (i in 2:NP) {
      par.initial = rbind(par.initial, runif(nvars,t(Domains[,1]), t(Domains[,2])))
    }

    # Set optimizer options
    #if parallel
    if (doParallel==T){
      controls = list(initialpop=par.initial,reltol=reltol, steptol=steptol, itermax=max.generations, trace=print.iterations, NP=NP,
                      c=0.01, strategy=DEstrategy, parallelType = "auto",...)
    }else{
      controls = list(initialpop=par.initial,reltol=reltol, steptol=steptol, itermax=max.generations, trace=print.iterations, NP=NP,
                      c=0.01, strategy=DEstrategy)
    }
  } else {
    if (doParallel==T){
      controls = list(reltol=reltol, steptol=steptol, itermax=max.generations, trace=print.iterations, NP=NP, c=0.01,
                      strategy=DEstrategy, parallelType = "auto", ...)
    }else{
      controls = list(reltol=reltol, steptol=steptol, itermax=max.generations, trace=print.iterations, NP=NP, c=0.01,
                      strategy=DEstrategy)
      # message("controls here")
    }
  }


  message('... Starting calibration using the following settings:')
  message(paste('    - Initial parameter set neg. log liklihood:',obj.initial))
  message(paste('    - total population size:',nvars*pop.size.perParameter))
  message(paste('    - relative tolerance:',reltol))
  message(paste('    - iterations required that meet the tolerance:',steptol))
  message(paste('    - maximum iterations allowed:',max.generations))
  message(paste('    - DEoptim strategy type:',DEstrategy))


  # Run optimiser and assign solution to the object.
  calib.results <- DEoptim(getNegLogLikelihood.fromTransformedVector,
                           lower = as.vector(Domains[,1]),
                           upper = as.vector(Domains[,2]),
                           control=  controls, .Object=.Object)


  # Add calibration outputs to the object
  .Object@calibration.results <- calib.results$optim

  # Set calibrated parameters
  .Object <- setParameters.fromTransformed.asVector(.Object, .Object@calibration.results$bestmem)

  message('... Finished Calibration.')
  message(paste('    Best solution:',.Object@calibration.results$bestval))

  # Check the modle is valid
  if (!validObject(.Object))
    warning('The model parameters produced an INVALID MODEL.')

  return(.Object)
}
)


# @exportMethod setStateNames
setGeneric(name="setStateNames",def=function(.Object, year.normalFlow) {standardGeneric("setStateNames")})
setMethod(f="setStateNames",signature=c("hydroState","numeric"),definition=function(.Object, year.normalFlow)
{

  if (!validObject(.Object))
    stop('The model parameters produced an INVALID MODEL.')

  # Get the viterbi state at the normal year.
  states.viterbi <- viterbi(.Object, do.plot=F)

  # Find the year defining normal flow and precipitation.
  for (i in 1:length(year.normalFlow)) {
    filt = states.viterbi[,1] == year.normalFlow[i];

    if (!any(is.na(states.viterbi[filt,'Viterbi State Number']))){
      states.viterbi <- states.viterbi[filt,]
      year.normalFlow = year.normalFlow[i]
      break
    }
  }


  # Assess if the model has a monthly or annal time step and get the
  # 'normal' precipitation. If the monthly, then
  # the monthly precip will be taken rather than the annual average.
  year.filt = .Object@input.data$year ==year.normalFlow
  if (any(names(.Object@input.data)=="month")) {
    is.monthly = T
    Pavg = cbind(.Object@input.data$year[year.filt], .Object@input.data$month[year.filt], .Object@input.data$precipitation[year.filt])
  } else {
    is.monthly = F
    Pavg = c(year.normalFlow, .Object@input.data$precipitation[year.filt])
  }

  # If monthly, get the most frequent state for the months of the year.
  if (is.monthly) {
    states.viterbi = c(year.normalFlow, as.numeric(names(sort(table(states.viterbi[,3]),decreasing=TRUE)[1])))
  }

  # Derive the Qhat value at each state at the precip defined as a normal year.
  # To achieve this the input data precipitation is replaced by the above mean. This is
  # undertaken to allow for estimates from the QhatModel object that require a timeseries
  # rather than a time-series, eg those wih serial correlation.
  #-----------

  # Assign normal precip conditions to all obs data
  data  = .Object@input.data
  data = getQhat(.Object@Qhat.object, .Object@input.data)
  Qhat = data$Qhat.flow
  if (is.monthly) {
    for (i in 1:nrow(Pavg)) {
     filt = data$month == Pavg[i,2]
     data$precipitation[filt] = Pavg[i,3]
    }
  } else {
    data$precipitation = Pavg[2]
  }

  # get number of states
  nStates = getNumStates(.Object@markov.model.object)

  # Get the MEDIAN QhatModel value of Qhat at normal precip and at the normal year
  state.est <- getDistributionPercentiles(.Object@QhatModel.object, data, 0.5)
  state.est <- state.est[[1]]
  state.est = state.est[year.filt,]
  if (is.monthly) {
    if (nStates==1) {
      state.est = sum(state.est)
    } else {
      state.est = colSums(state.est, na.rm = TRUE)
    }
  }
  states.normal <- states.viterbi[[2]]
  state.est.normal <- state.est[states.normal]

  # Assign names to the states!
  .Object@state.labels = vector("character",nStates)
  .Object@state.labels[states.normal] = 'Normal'
  for (i in 1:nStates) {

    if (i==states.normal)
      next

    if (state.est[i] == state.est.normal) {
      # This cae is for 2+ states having the same state as mean as occurs for the normal year.
      .Object@state.labels[i] = 'Normal (duplicate)';
    } else  if (state.est[i] < state.est.normal && state.est[i] == min(state.est) && sum(state.est<state.est.normal)==2) {
      .Object@state.labels[i] = 'Very low';
    } else if ( state.est[i] < state.est.normal) {
      .Object@state.labels[i] = 'Low';
    } else if (state.est[i] > state.est.normal && state.est[i] == max(state.est) && sum(state.est>state.est.normal)==2) {
      .Object@state.labels[i] = 'Very high';
    } else if (state.est[i] > state.est.normal) {
      .Object@state.labels[i] = 'High';
    } else {
      #.Object@state.labels[i] = '(State label error)';
      stop('State label error')
    }
  }

  return(.Object)
}
)



# @export plot_graph
setGeneric(name="plot_graph",def=function(.Object, main=NA, relsize=NA) {standardGeneric("plot_graph")})
setMethod(f="plot_graph",signature="hydroState",definition=function(.Object, main='Transtion Probability Graph', relsize=0.8)
{

  # Get Transition probs. and round to 4 digits
  Tprob<- getTransitionProbabilities(.Object)
  Tprob <- round(Tprob,4)

  # Set Tprob to -Inf where graph
  T.graph <- .Object@markov.model.object@transition.graph;
  Tprob[!T.graph]=-Inf

  # Set the posotion of the self arrow
  nStates = getNumStates(.Object@markov.model.object)
  if (nStates==1) {
    self.shiftx =-0.1;
    self.shifty =0.1;
  } else if (nStates==2) {
    self.shiftx =c(-0.1,0.1);
    self.shifty =c(0.1,0.1)
  } else if (nStates==3) {
    self.shiftx =c(0.1,-0.1,-0.1);
    self.shifty =c(0.1,0.1,0.1);
  } else {
    self.shiftx = NULL;
  }

   if (length(.Object@state.labels)==0){
     plotmat(Tprob,absent=-Inf, , main=main,self.shiftx=self.shiftx, self.shifty=self.shifty,cex.txt=0.8, relsize=relsize)
   } else {
     plotmat(Tprob,name=.Object@state.labels,absent=-Inf, main=main,self.shiftx=self.shiftx,self.shifty=self.shifty, , cex.txt=0.8, relsize=relsize)
   }

}
)

# @exportMethod viterbiF\
setGeneric(name="viterbi",def=function(.Object, data, do.plot=NA, plot.percentiles=NA, plot.yearRange=NA, plot.options=NA) {standardGeneric("viterbi")})
setMethod(f="viterbi",signature=c("hydroState","missing","missing","missing","missing","missing"),
          definition=function(.Object, data, do.plot=T, plot.percentiles = c(0.05, 0.5, 0.95), plot.yearRange=numeric(),plot.options = c("A","B","C","D"))
{

  if (!validObject(.Object))
    stop('The model parameters produced an INVALID MODEL.')

  # Get the transformed flow, Qhat, and add Qhat to the input data
  data = getQhat(.Object@Qhat.object, .Object@input.data)
  # data = cbind.data.frame(.Object@input.data,Qhat)

  # Run the viterbi algorithm
  states <- viterbi(.Object, data, do.plot=T, plot.percentiles = c(0.05, 0.5, 0.95), plot.yearRange, plot.options)
  return(states)
}
)
setMethod(f="viterbi",signature=c("hydroState","missing","logical","missing","missing","character"),
          definition=function(.Object, data, do.plot, plot.percentiles, plot.yearRange, plot.options)
          {

            if (!validObject(.Object))
              stop('The model parameters produced an INVALID MODEL.')

            # Get the transformed flow, Qhat, and add Qhat to the input data
            data = getQhat(.Object@Qhat.object, .Object@input.data)

            # Run the viterbi algorithm
            states <- viterbi(.Object, data, do.plot, c(0.05, 0.5, 0.95), numeric(),plot.options)
            return(states)
          }
)
setMethod(f="viterbi",signature=c("hydroState","missing","logical","numeric","numeric","character"),
          definition=function(.Object, data, do.plot=T, plot.percentiles = c(0.05, 0.5, 0.95), plot.yearRange=numeric(),plot.options = c("A","B","C","D"))
{
  if (!validObject(.Object))
    stop('The model parameters produced an INVALID MODEL.')

  # Get the transformed flow, Qhat, and add Qhat to the input data
  data = getQhat(.Object@Qhat.object, .Object@input.data)
  # data = cbind.data.frame(.Object@input.data,Qhat)

  # Run the viterbi algorithm
  states <- viterbi(.Object, data, do.plot, plot.percentiles, plot.yearRange,plot.options)
  return(states)
}
)

setMethod(f="viterbi",signature=c("hydroState","missing","logical","missing","missing","missing"),
          definition=function(.Object, data, do.plot, plot.percentiles, plot.yearRange, plot.options)
          {

            if (!validObject(.Object))
              stop('The model parameters produced an INVALID MODEL.')

            # Get the transformed flow, Qhat, and add Qhat to the input data
            data = getQhat(.Object@Qhat.object, .Object@input.data)

            # Run the viterbi algorithm
            states <- viterbi(.Object, data, do.plot, c(0.05, 0.5, 0.95), numeric())
            return(states)
          }
)
setMethod(f="viterbi",signature=c("hydroState","missing","logical","numeric","numeric","missing"),
          definition=function(.Object, data, do.plot=T, plot.percentiles = c(0.05, 0.5, 0.95), plot.yearRange=numeric(),plot.options = c("A","B","C","D"))
          {
            if (!validObject(.Object))
              stop('The model parameters produced an INVALID MODEL.')

            # Get the transformed flow, Qhat, and add Qhat to the input data
            data = getQhat(.Object@Qhat.object, .Object@input.data)
            # data = cbind.data.frame(.Object@input.data,Qhat)

            # Run the viterbi algorithm
            states <- viterbi(.Object, data, do.plot, plot.percentiles, plot.yearRange,plot.options)
            return(states)
          }
)

setMethod(f="viterbi",signature=c("hydroState","data.frame","logical","numeric","numeric","missing"),
          definition=function(.Object, data, do.plot=T, plot.percentiles = c(0.05, 0.5, 0.95), plot.yearRange=numeric(),plot.options)
  {

  if (!validObject(.Object))
    stop('The model parameters produced an INVALID MODEL.')

  # Handle monthly and yearly time steps
  if (any(names(data)=="month")) {
    obsDates.asISO = as.Date(ISOdate(data$year,data$month,1))
    obsDates = cbind(data$year,data$month)
    obsDates.Precip.asISO = as.Date(ISOdate(.Object@input.data$year,.Object@input.data$month,1))
    plot.units = 'month'
  } else {
    obsDates.asISO = as.Date(ISOdate(data$year,1,1))
    obsDates = data$year
    obsDates.Precip.asISO = as.Date(ISOdate(.Object@input.data$year,1,1))
    plot.units = 'year'
  }

  # Get the transformed flow, Qhat
  if (!any(names(data)=='Qhat.flow'))
    stop('Input "data" must contain a column named "Qhat.flow".')

  data.withNAs <- data
  Qhat = data$Qhat.flow
  nQhat = length(Qhat)

  # get number of states
  nStates = getNumStates(.Object@markov.model.object)

  # Built filter for non NAs.
  filt <- !is.na(data$Qhat.flow)

  # Exit if modle has one state
  if (nStates==1) {
    viterbiPath <- rep(1,length(Qhat))

    if (is.vector(obsDates)) {
      results <- matrix(NA,length(filt), 3)
      results[,1] <- obsDates
      results[filt,2] <-viterbiPath[filt]
      results[,3] <-data$flow
      colnames(results) <- c('Year','Viterbi State Number', 'Obs. flow')
    } else {
      results <- matrix(NA,length(filt), 4)
      results[,1:2] <- obsDates
      results[filt,3] <-viterbiPath[filt]
      results[,4] <-data$flow
      colnames(results) <- c('Year','Month','Viterbi State Number', 'Obs. flow')
    }
  } else {
    # Get transition probs.
    transProbs = getTransitionProbabilities(.Object@markov.model.object)

    # get emiision probs.
    emissionProbs = getEmissionDensity(.Object@QhatModel.object, data, NA)

    # Get initial states
    startProbs = getInitialStateProbabilities(.Object@markov.model.object)
    States = 1:nStates

    # Set probs to zero if the is no obs. data
    transProbs[is.na(transProbs)]       = 0
    #emissionProbs[is.na(emissionProbs)] = 0
    emissionProbs[!filt,] = 1
    nQhat <- nrow(data)

    # Set zero emmision probs to machine percision
    for(state in States) {
      filt.zeros = emissionProbs[,state]==0
      emissionProbs[filt.zeros,state] = .Machine$double.eps
    }

    #Run Viterbi algorithm. Adapted from https://cran.r-project.org/web/packages/HMM
    v <- array(NA,c(nStates,nQhat))
    dimnames(v) <- list(states=States,index=1:nQhat)
    # Init
    for(state in States)
    {
      v[state,1] = log(startProbs[state] * emissionProbs[1,state])
    }
    # Iteration
    for(k in 2:nQhat)
    {
      for(state in States)
      {
        maxi = NULL
        for(previousState in States)
        {
          temp = v[previousState,k-1] + log(transProbs[previousState,state])
          maxi = max(maxi, temp)
        }
        v[state,k] = log(emissionProbs[k,state]) + maxi
      }
    }
    # Traceback
    viterbiPath = rep(NA,nQhat)
    for(state in States)
    {
      if(max(v[,nQhat])==v[state,nQhat])
      {
        viterbiPath[nQhat] = state
        break
      }
    }
    for(k in (nQhat-1):1)
    {
      for(state in States)
      {
        if(max(v[,k]+log(transProbs[,viterbiPath[k+1]]))
           ==v[state,k]+log(transProbs[state,viterbiPath[k+1]]))
        {
          viterbiPath[k] = state
          break
        }
      }
    }
  }

  # Caculate the flow at the percentile AND if the catchment was at the normal state.
  # This requires the state names to be set.
  #--------------
  if (.Object@state.labels[1]!='') {

    if (length(plot.percentiles)!=3)
      stop('The input "percentiles" must contain only three values.')
    plot.percentiles =  sort(plot.percentiles)

    # Get the precentiles each state at each time point and remove rows with no obs Qhat
    state.est = getDistributionPercentiles(.Object@QhatModel.object, data, plot.percentiles)

    # Remove NAs from the input data and Qhat
    data <- data[filt,]
    Qhat <- Qhat[filt]
    obsDates.withNAs <- obsDates
    if (is.vector(obsDates)) {
      obsDates <- obsDates[filt]
    } else {
      obsDates <- obsDates[filt,]
    }
    obsDates.asISO.withNAs <- obsDates.asISO
    obsDates.asISO <- obsDates.asISO[filt]
    viterbiPath <- viterbiPath[filt]
    nQhat = sum(filt)
    for (i in 1:3) {
      state.est[[i]] <- as.matrix(state.est[[i]][filt,], ncol=nStates)
    }

    # Keep the state.ext for the 'Normal' flow state.
    ind.stateNames.normal <- which(.Object@state.labels=='Normal')
    state.est.normal <- matrix(0, nQhat,3)
    for (i in 1:3) {
      for (j in 1:nQhat) {
        state.est.normal[j,i] = state.est[[i]][j,ind.stateNames.normal]
      }
    }

    # Create an index for if the catchment is in a normal flow state.
    ind.flow.normal <- viterbiPath==ind.stateNames.normal

    # Extract the percentile Qhat values for each percentile and each time step
    viterbi.est =matrix(0, nQhat,3)

    for (i in 1:3) {
      for (j in 1:nQhat) {
        viterbi.est[j,i] = state.est[[i]][j,viterbiPath[j]]
      }
    }

    # Back Transform Qhat models estimates to flow
    flow.viterbi.est =matrix(0, nQhat,3)
    flow.normal.est =matrix(0, nQhat,3)
    for (i in 1:length(plot.percentiles)) {
      data.tmp <- data.withNAs;
      data.tmp$Qhat.flow <- NA
      data.tmp$Qhat.flow[filt] = viterbi.est[,i]
      flow.viterbi.est[,i] <- getQ.backTransformed(.Object@Qhat.object,data.tmp)$flow.modelled[filt]

      data.tmp <- data.withNAs;
      data.tmp$Qhat.flow <- NA
      data.tmp$Qhat.flow[filt] = state.est.normal[,i]
      flow.normal.est[,i] <- getQ.backTransformed(.Object@Qhat.object,data.tmp)$flow.modelled[filt]
    }

    # Get the conditional probabilities.
    emissionProbs = getEmissionDensity(.Object@QhatModel.object, data, NA)
    state.probs = getConditionalStateProbabilities(.Object@markov.model.object, data, emissionProbs)

    # Collate returned data.
    if (is.vector(obsDates)) {
      results <- matrix(NA,length(filt), 9+2*nStates)
      results[,1] <- obsDates.withNAs
      results[filt,2] <-viterbiPath
      results[,3] <-data.withNAs$flow
      results[filt,4:9] <- cbind(flow.viterbi.est, flow.normal.est)
      results[filt, 10:(10+nStates-1)] = t(state.probs)
      results[filt, (10+nStates):(10+2*nStates-1)] = t(emissionProbs)
      colnames(results) <- c('Year','Viterbi State Number', 'Obs. flow',
                             paste('Viterbi Flow -',plot.percentiles*100,'%ile',sep=''),
                             paste('Normal State Flow-',plot.percentiles*100,'%ile',sep=''),
                             paste('Conditional Prob.-',.Object@state.labels,sep=''),
                             paste('Emission Density-',.Object@state.labels,sep=''))
    } else {
      results <- matrix(NA,length(filt), 10+2*nStates)
      results[,1:2] <- obsDates.withNAs
      results[filt,3] <-viterbiPath
      results[,4] <-data.withNAs$flow
      results[filt,5:10] <- cbind(flow.viterbi.est, flow.normal.est)
      results[filt, 11:(11+nStates-1)] = t(state.probs)
      results[filt, (11+nStates):(11+2*nStates-1)] = t(emissionProbs)
      colnames(results) <- c('Year','Month','Viterbi State Number', 'Obs. flow',
                             paste('Viterbi Flow -',plot.percentiles*100,'%ile',sep=''),
                             paste('Normal State Flow -',plot.percentiles*100,'%ile',sep=''),
                             paste('Conditional Prob.-',.Object@state.labels,sep=''),
                             paste('Emission Density-',.Object@state.labels,sep=''))
    }

  } else {
    # Collate returned data.
    if (is.vector(obsDates)) {
      results <- matrix(NA,length(filt), 3)
      results[,1] <- obsDates
      results[filt,2] <-viterbiPath[filt]
      results[,3] <-data$flow
      colnames(results) <- c('Year','Viterbi State Number', 'Obs. flow')
    } else {
      results <- matrix(NA,length(filt), 4)
      results[,1:2] <- obsDates
      results[filt,3] <-viterbiPath[filt]
      results[,4] <-data$flow
      colnames(results) <- c('Year','Month','Viterbi State Number', 'Obs. flow')
    }
  }

  # Do plotting
  if (do.plot) {
    if (length(.Object@state.labels)==1 && .Object@state.labels=='') {
      warning('State names must be set for plotting.')
      return(results)
    }

    # Define colours for the states (derived from command: brewer.pal(5,"Spectral"))
    state.colours = rep("grey",nStates)
    if (length(.Object@state.labels)==nStates) {
      state.colours.all = c("#D7191C","#FDAE61", "#ABDDA4",'#7c9fb6',"#496c83")
      for ( i in 1:nStates) {
        if (.Object@state.labels[i]=='Very low')
          state.colours[i] = state.colours.all[1]

        if (.Object@state.labels[i]=='Low')
          state.colours[i] = state.colours.all[2]

        if (.Object@state.labels[i]=='Normal')
          state.colours[i] = state.colours.all[3]

        if (.Object@state.labels[i]=='High')
          state.colours[i] = state.colours.all[4]

        if (.Object@state.labels[i]=='Very high')
          state.colours[i] = state.colours.all[5]
      }

    }

    # Derive a matrix for the start and end of each bar for the plotting of the range in the percentiles
    line.matrix <- cbind(viterbi.est[,1], viterbi.est[,3])
    flow.line.matrix <- cbind(flow.viterbi.est[,1], flow.viterbi.est[,3])

    # Get input grapics settings
    op <- par(no.readonly = T);

    # Change graphics settings
    nrow.plots = 4
    layout(matrix(c(1,2,2,2,3,3,4,4), 8, 1, byrow = TRUE))
    par(mar = c(0.2,5,0.2,5))


    # Calc axis limits.
    ylim.flow <- c(0, ceiling(max( c(max(data$flow,na.rm=T),max(flow.viterbi.est,na.rm=T)))))
    ylim.precip.max = max(ylim.flow)*3
    ylim.precip <- c(-ylim.precip.max,0)

    # Setup year range for plotitng
    if (length(plot.yearRange)==2 && all(is.numeric(plot.yearRange)) && all(plot.yearRange>0) && all(plot.yearRange<=as.numeric(format(Sys.Date(), "%Y")))) {
      xlim = as.Date(c(ISOdate(plot.yearRange[1],1,1), ISOdate(plot.yearRange[2],12,31)))
    } else {
      xlim = c(min(obsDates.asISO)-7, max(obsDates.asISO)+7)
    }

    # Plot obs precip
    # pframe = padr::thicken(padr::pad(data.frame(obsDates.asISO.withNAs, .Object@input.data$precipitation), interval = plot.units), "3 month")
    pframe = padr::pad(data.frame(obsDates.asISO.withNAs, .Object@input.data$precipitation), interval = plot.units)
    plot(pframe, type='s',col='grey', lwd=1,
         xlim=xlim, xlab='', ylab='', main='', xaxt='n')
    mtext("Precip.",side=2,line=3)
    mtext(paste("[mm/",plot.units,"]",sep=''),side=2,line=2, cex=0.85)

    xaxis.ticks = as.Date(ISOdate(seq(1900,2020,by=10),1,1))
    abline(v=xaxis.ticks, col = "lightgray", lty = "dotted",lwd = par("lwd"))
    grid(NA,NULL)
    plot.range=par("usr")
    text(plot.range[1]+diff(plot.range[1:2])*0.025, plot.range[3]+diff(par("usr")[3:4])*0.95, labels="", font=1, cex=2,pos=1)

    # Plot obs flow
    pframe = padr::pad(data.frame(obsDates.asISO.withNAs, data.withNAs$flow), interval = plot.units)
    plot(pframe, type='l',col='grey', lwd=1, ylim=ylim.flow, xlim=xlim,
         xlab='', ylab='',xaxt='n')
    # ryticks.min = signif(1.2*max(data$precipitation,na.rm=T),1)
    # ryticks = seq(-ryticks.min, 0, by=signif(ryticks.min/2,1))
    # print(-1*(data$precipitation))
    # twoord.plot(obsDates.asISO, data$flow,
    #             obsDates.asISO, -1*(data$precipitation),
    #             lylim=ylim.flow,rylim= ylim.precip,type=c("l","s"),
    #             xlab="Year", ylab=paste("Flow [mm/",plot.units,"]",sep=''),rylab=paste("Precip [mm/",plot.units,"]",sep=''),
    #             lytickpos=seq(0,max(ylim.flow), by=signif(max(ylim.flow)/5,1)), rytickpos=ryticks,
    #             lcol='black',rcol='blue')



    # Plot Markov states as boxes
    for (i in 1:nQhat) {
      points(obsDates.asISO[i], flow.viterbi.est[i,2],col=state.colours[viterbiPath[i]], bg=state.colours[viterbiPath[i]], pch=21)
      lines(rep(obsDates.asISO[i],2), flow.line.matrix[i,],col=state.colours[viterbiPath[i]], lwd=1)

      # Plot the normal flow in years when the Viterbi state is not normal.
      if (!ind.flow.normal[i]) {
        points(obsDates.asISO[i], flow.normal.est[i,2],col=state.colours.all[3], bg=state.colours.all[3], pch=1)
        lines(rep(obsDates.asISO[i],2), c(flow.normal.est[i,1],flow.normal.est[i,3]),col=state.colours.all[3], lwd=1, lty=1)
      }
    }


    # Add axis labels and legend
    mtext("Flow",side=2,line=3)
    mtext(paste("[mm/",plot.units,"]"),side=2,line=2, cex=0.85)
    #mtext('Year',side=1,line=2)
    #legend('topleft', legend=.Object@state.labels, pch=21, col=state.colours, pt.bg=state.colours, xjust=0)
    abline(v=xaxis.ticks, col = "lightgray", lty = "dotted",lwd = par("lwd"))
    grid(NA,NULL)
    legend('topright', legend=c('Obs. flow',paste(.Object@state.labels,' (5th - 50th - 95th)',sep=''),'Est. normal (5th - 50th - 95th)'),
           lty=c(1,1,1,1),lwd=1,pch=c(NA,21,21,1), col=c('grey',state.colours,state.colours.all[3]),
           pt.bg=c(NA,state.colours,NA), xjust=0, cex=1, bg='white')
    plot.range=par("usr")
    text(plot.range[1]+diff(plot.range[1:2])*0.025, plot.range[3]+diff(par("usr")[3:4])*0.95, labels="", font=1, cex=2,pos=1)

    # Calc axis limits.
    ylim.qhat <- c(floor(min( c(min(data$Qhat.flow,na.rm=T),min(viterbi.est,na.rm=T)))) ,
                   ceiling(max( c(max(data$Qhat.flow,na.rm=T),max(viterbi.est,na.rm=T)))))

    # Plot obs Qhat
    pframe = padr::pad(data.frame(obsDates.asISO.withNAs, data.tmp$Qhat.flow), interval = plot.units)
    plot(pframe, type='l',col='grey', lwd=1,
         ylim=ylim.qhat, xlim=xlim, xlab='', ylab='',xaxt='n')

    # Plot Markov states as boxes
    for (i in 1:nQhat) {
      points(obsDates.asISO[i], viterbi.est[i,2],col=state.colours[viterbiPath[i]], bg=state.colours[viterbiPath[i]], pch=21)
      lines(rep(obsDates.asISO[i],2), line.matrix[i,],col=state.colours[viterbiPath[i]], lwd=1)
    }

    # Add axis labels
    mtext("Transformed flow",side=2,line=3)
    mtext(paste("[f(mm) /",plot.units,"]"),side=2,line=2, cex=0.85)
    # mtext('Year',side=1,line=2)
    abline(v=xaxis.ticks, col = "lightgray", lty = "dotted",lwd = par("lwd"))
    grid(NA,NULL)
    plot.range=par("usr")
    text(plot.range[1]+diff(plot.range[1:2])*0.025, plot.range[3]+diff(par("usr")[3:4])*0.95, labels="", font=1, cex=2,pos=1)


    # Plot the cummulative rainfall residual.
    #--------------

    # # Calculate the means and residuals
    # if (plot.units == 'yr') {
    #   P.mean = mean(data$precipitation)
    #   P.resid = data$precipitation - P.mean;
    # } else {
    #   P.mean = rep(NA,length(.Object@QhatModel.object@subAnnual.Monthly.Steps))
    #   P.resid = rep(NA, length(data$precipitation))
    #   for (i in 1:length(.Object@QhatModel.object@subAnnual.Monthly.Steps)) {
    #    filt = data$month == .Object@QhatModel.object@subAnnual.Monthly.Steps[i]
    #    P.mean[i] =  mean(data$precipitation[filt])
    #    P.resid[filt] = data$precipitation[filt] - P.mean[i]
    #   }
    # }
    #
    # # # Plot the residuals
    # # plot(obsDates.asISO, P.resid,type='p',xlim=xlim, col='white', bg='white', pch=21, xlab='', ylab='')
    # # for (i in 1:nQhat) {
    # #   points(obsDates.asISO[i], P.resid[i],col=state.colours[viterbiPath[i]], bg=state.colours[viterbiPath[i]], pch=21)
    # # }
    # # lines(obsDates.asISO, rep(0,length(obsDates.asISO)),col='grey')
    # # mtext("Rainfall residual [mm]",side=2,line=3)
    # # mtext('Year',side=1,line=2)
    #
    # # Calculate the cumulative residuals.
    # P.cumResid = cumsum(P.resid)
    #
    # # Plot the cumm residuals
    # plot(obsDates.asISO, P.cumResid, type='l',col='grey', lwd=1, xlim=xlim, xlab='', ylab='',xaxt='n')
    # abline(v=xaxis.ticks, col = "lightgray", lty = "dotted",lwd = par("lwd"))
    # grid(NA,NULL)
    # legend('topright', legend=c('Cum. residual ',.Object@state.labels),
    #        lty=c(1,NA,NA),pch=c(NA,21,21), col=c('grey',state.colours),
    #        pt.bg=c(NA,state.colours), xjust=0, cex=1.5, bg='white')
    #
    # # Colour the points by the Viterbi state.
    # for (i in 1:nQhat) {
    #   points(obsDates.asISO[i], P.cumResid[i],col=state.colours[viterbiPath[i]], bg=state.colours[viterbiPath[i]], pch=21)
    # }
    # mtext("Cum. rainfall resid.",side=2,line=3)
    # mtext(paste("[mm]"),side=2,line=2, cex=0.85)
    # plot.range=par("usr")
    # text(plot.range[1]+diff(plot.range[1:2])*0.025, plot.range[3]+diff(par("usr")[3:4])*0.95, labels="D", font=1, cex=2,pos=1)

    # Plot the conditional state probability for each state
    #if (nStates>1) {
      emissionProbs = getEmissionDensity(.Object@QhatModel.object, data, NA)
      state.probs = getConditionalStateProbabilities(.Object@markov.model.object, data, emissionProbs)

      # Plot bar graph
      par(mar = c(4,5,0.2,5))
      plot(obsDates.asISO, state.probs[1,], type = 'l', xlim=xlim, col=state.colours[1],
           ylim=c(0,1),xlab='', ylab='', lwd=1, xaxt='n')
      if (nStates>1) {
        for ( i in 2:nStates) {
          lines(obsDates.asISO, state.probs[i,], col=state.colours[i])
        }
      }
      mtext("State Prob.",side=2,line=3)
      mtext(paste("[-]"),side=2,line=2, cex=0.85)
      mtext('Year',side=1,line=3)
      axis(1,at= xaxis.ticks,labels=seq(1900,2020,by=10))
      abline(v=xaxis.ticks, col = "lightgray", lty = "dotted",lwd = par("lwd"))
      grid(NA,NULL)
      legend('bottomleft', legend=.Object@state.labels,
             lty=c(1,1),pch=c(NA,NA), col=state.colours,lwd=1,
             xjust=0, cex=1, bg='white')
      plot.range=par("usr")
      text(plot.range[1]+diff(plot.range[1:2])*0.025, plot.range[3]+diff(par("usr")[3:4])*0.95, labels="", font=1, cex=2,pos=1)


    #}

    # Reset graphics options
    #--------------
    par(op)

  }

  # Return data
  return(as.data.frame(results))
  }
)

setMethod(f="viterbi",signature=c("hydroState","data.frame","logical","numeric","numeric","character"),
          definition=function(.Object, data, do.plot=T, plot.percentiles = c(0.05, 0.5, 0.95), plot.yearRange=numeric(),plot.options = c())
          {

            if (!validObject(.Object))
              stop('The model parameters produced an INVALID MODEL.')

            # Handle monthly and yearly time steps
            if (any(names(data)=="month")) {
              obsDates.asISO = as.Date(ISOdate(data$year,data$month,1))
              obsDates = cbind(data$year,data$month)
              obsDates.Precip.asISO = as.Date(ISOdate(.Object@input.data$year,.Object@input.data$month,1))
              plot.units = 'month'
            } else {
              obsDates.asISO = as.Date(ISOdate(data$year,1,1))
              obsDates = data$year
              obsDates.Precip.asISO = as.Date(ISOdate(.Object@input.data$year,1,1))
              plot.units = 'year'
            }

            # message(obsDates.Precip.asISO[1:27])
            # Get the transformed flow, Qhat
            if (!any(names(data)=='Qhat.flow'))
              stop('Input "data" must contain a column named "Qhat.flow".')

            data.withNAs <- data
            Qhat = data$Qhat.flow
            nQhat = length(Qhat)

            # get number of states
            nStates = getNumStates(.Object@markov.model.object)

            # Built filter for non NAs.
            filt <- !is.na(data$Qhat.flow)

            # Exit if modle has one state
            if (nStates==1) {
              viterbiPath <- rep(1,length(Qhat))

              if (is.vector(obsDates)) {
                results <- matrix(NA,length(filt), 3)
                results[,1] <- obsDates
                results[filt,2] <-viterbiPath[filt]
                results[,3] <-data$flow
                colnames(results) <- c('Year','Viterbi State Number', 'Obs. flow')
              } else {
                results <- matrix(NA,length(filt), 4)
                results[,1:2] <- obsDates
                results[filt,3] <-viterbiPath[filt]
                results[,4] <-data$flow
                colnames(results) <- c('Year','Month','Viterbi State Number', 'Obs. flow')
              }
            } else {
              # Get transition probs.
              transProbs = getTransitionProbabilities(.Object@markov.model.object)

              # get emiision probs.
              emissionProbs = getEmissionDensity(.Object@QhatModel.object, data, NA)

              # Get initial states
              startProbs = getInitialStateProbabilities(.Object@markov.model.object)
              States = 1:nStates

              # Set probs to zero if the is no obs. data
              transProbs[is.na(transProbs)]       = 0
              #emissionProbs[is.na(emissionProbs)] = 0
              emissionProbs[!filt,] = 1
              nQhat <- nrow(data)

              # Set zero emmision probs to machine percision
              for(state in States) {
                filt.zeros = emissionProbs[,state]==0
                emissionProbs[filt.zeros,state] = .Machine$double.eps
              }

              #Run Viterbi algorithm. Adapted from https://cran.r-project.org/web/packages/HMM
              v <- array(NA,c(nStates,nQhat))
              dimnames(v) <- list(states=States,index=1:nQhat)
              # Init
              for(state in States)
              {
                v[state,1] = log(startProbs[state] * emissionProbs[1,state])
              }
              # Iteration
              for(k in 2:nQhat)
              {
                for(state in States)
                {
                  maxi = NULL
                  for(previousState in States)
                  {
                    temp = v[previousState,k-1] + log(transProbs[previousState,state])
                    maxi = max(maxi, temp)
                  }
                  v[state,k] = log(emissionProbs[k,state]) + maxi
                }
              }
              # Traceback
              viterbiPath = rep(NA,nQhat)
              for(state in States)
              {
                if(max(v[,nQhat])==v[state,nQhat])
                {
                  viterbiPath[nQhat] = state
                  break
                }
              }
              for(k in (nQhat-1):1)
              {
                for(state in States)
                {
                  if(max(v[,k]+log(transProbs[,viterbiPath[k+1]]))
                     ==v[state,k]+log(transProbs[state,viterbiPath[k+1]]))
                  {
                    viterbiPath[k] = state
                    break
                  }
                }
              }
            }

            # Caculate the flow at the percentile AND if the catchment was at the normal state.
            # This requires the state names to be set.
            #--------------
            if (.Object@state.labels[1]!='') {

              if (length(plot.percentiles)!=3)
                stop('The input "percentiles" must contain only three values.')
              plot.percentiles =  sort(plot.percentiles)

              # Get the precentiles each state at each time point and remove rows with no obs Qhat
              state.est = getDistributionPercentiles(.Object@QhatModel.object, data, plot.percentiles)

              # Remove NAs from the input data and Qhat
              data <- data[filt,]
              Qhat <- Qhat[filt]
              obsDates.withNAs <- obsDates
              if (is.vector(obsDates)) {
                obsDates <- obsDates[filt]
              } else {
                obsDates <- obsDates[filt,]
              }
              obsDates.asISO.withNAs <- obsDates.asISO
              obsDates.asISO <- obsDates.asISO[filt]
              viterbiPath <- viterbiPath[filt]
              nQhat = sum(filt)
              for (i in 1:3) {
                state.est[[i]] <- as.matrix(state.est[[i]][filt,], ncol=nStates)
              }

              # Keep the state.ext for the 'Normal' flow state.
              ind.stateNames.normal <- which(.Object@state.labels=='Normal')
              state.est.normal <- matrix(0, nQhat,3)
              for (i in 1:3) {
                for (j in 1:nQhat) {
                  state.est.normal[j,i] = state.est[[i]][j,ind.stateNames.normal]
                }
              }

              # Create an index for if the catchment is in a normal flow state.
              ind.flow.normal <- viterbiPath==ind.stateNames.normal

              # Extract the percentile Qhat values for each percentile and each time step
              viterbi.est =matrix(0, nQhat,3)

              for (i in 1:3) { # code for baseflow
                for (j in 1:nQhat) {
                  viterbi.est[j,i] = state.est[[i]][j,viterbiPath[j]]
                }
              }

              # Back Transform Qhat models estimates to flow
              flow.viterbi.est =matrix(0, nQhat,3)
              flow.normal.est =matrix(0, nQhat,3)
              for (i in 1:length(plot.percentiles)) {
                data.tmp <- data.withNAs;
                data.tmp$Qhat.flow <- NA
                data.tmp$Qhat.flow[filt] = viterbi.est[,i]
                flow.viterbi.est[,i] <- getQ.backTransformed(.Object@Qhat.object,data.tmp)$flow.modelled[filt]

                data.tmp <- data.withNAs;
                data.tmp$Qhat.flow <- NA
                data.tmp$Qhat.flow[filt] = state.est.normal[,i]
                flow.normal.est[,i] <- getQ.backTransformed(.Object@Qhat.object,data.tmp)$flow.modelled[filt]
              }

              # Get the conditional probabilities.
              emissionProbs = getEmissionDensity(.Object@QhatModel.object, data, NA)
              state.probs = getConditionalStateProbabilities(.Object@markov.model.object, data, emissionProbs)

              # Collate returned data.
              if (is.vector(obsDates)) {
                results <- matrix(NA,length(filt), 9+2*nStates)
                results[,1] <- obsDates.withNAs
                results[filt,2] <-viterbiPath
                results[,3] <-data.withNAs$flow
                results[filt,4:9] <- cbind(flow.viterbi.est, flow.normal.est)
                results[filt, 10:(10+nStates-1)] = t(state.probs)
                results[filt, (10+nStates):(10+2*nStates-1)] = t(emissionProbs)
                colnames(results) <- c('Year','Viterbi State Number', 'Obs. flow',
                                       paste('Viterbi Flow -',plot.percentiles*100,'%ile',sep=''),
                                       paste('Normal State Flow-',plot.percentiles*100,'%ile',sep=''),
                                       paste('Conditional Prob.-',.Object@state.labels,sep=''),
                                       paste('Emission Density-',.Object@state.labels,sep=''))
              } else {
                results <- matrix(NA,length(filt), 10+2*nStates)
                results[,1:2] <- obsDates.withNAs
                results[filt,3] <-viterbiPath
                results[,4] <-data.withNAs$flow
                results[filt,5:10] <- cbind(flow.viterbi.est, flow.normal.est)
                results[filt, 11:(11+nStates-1)] = t(state.probs)
                results[filt, (11+nStates):(11+2*nStates-1)] = t(emissionProbs)
                colnames(results) <- c('Year','Month','Viterbi State Number', 'Obs. flow',
                                       paste('Viterbi Flow -',plot.percentiles*100,'%ile',sep=''),
                                       paste('Normal State Flow -',plot.percentiles*100,'%ile',sep=''),
                                       paste('Conditional Prob.-',.Object@state.labels,sep=''),
                                       paste('Emission Density-',.Object@state.labels,sep=''))
              }

            } else {
              # Collate returned data.
              if (is.vector(obsDates)) {
                results <- matrix(NA,length(filt), 3)
                results[,1] <- obsDates
                results[filt,2] <-viterbiPath[filt]
                results[,3] <-data$flow
                colnames(results) <- c('Year','Viterbi State Number', 'Obs. flow')
              } else {
                results <- matrix(NA,length(filt), 4)
                results[,1:2] <- obsDates
                results[filt,3] <-viterbiPath[filt]
                results[,4] <-data$flow
                colnames(results) <- c('Year','Month','Viterbi State Number', 'Obs. flow')
              }
            }

            # Do plotting  ### filter to adjust index if water year so plotting is not messed up?
            if (do.plot) {
              if (length(.Object@state.labels)==1 && .Object@state.labels=='') {
                warning('State names must be set for plotting.')
                return(results)
              }

              # Define colours for the states (derived from command: brewer.pal(5,"Spectral"))
              state.colours = rep("grey",nStates)
              if (length(.Object@state.labels)==nStates) {
                state.colours.all = c("#D7191C","#FDAE61", "#ABDDA4",'#7c9fb6',"#496c83")
                for ( i in 1:nStates) {
                  if (.Object@state.labels[i]=='Very low')
                    state.colours[i] = state.colours.all[1]

                  if (.Object@state.labels[i]=='Low')
                    state.colours[i] = state.colours.all[2]

                  if (.Object@state.labels[i]=='Normal')
                    state.colours[i] = state.colours.all[3]

                  if (.Object@state.labels[i]=='High')
                    state.colours[i] = state.colours.all[4]

                  if (.Object@state.labels[i]=='Very high')
                    state.colours[i] = state.colours.all[5]
                }

              }

              # Derive a matrix for the start and end of each bar for the plotting of the range in the percentiles
              line.matrix <- cbind(viterbi.est[,1], viterbi.est[,3])
              flow.line.matrix <- cbind(flow.viterbi.est[,1], flow.viterbi.est[,3])

              # Get input grapics settings
              op <- par(no.readonly = T);

              # Change graphics settings
              nrow.plots = length(plot.options)

              if(nrow.plots == 4){
                layout(matrix(c(1,1,2,2,3,3,4,4), 8, 1, byrow = TRUE))
              }
              if(nrow.plots == 3){
                layout(matrix(c(1,1,2,2,2,3,3,3), 8, 1, byrow = TRUE))
              }
              if(nrow.plots == 2){
                layout(matrix(c(1,1,1,1,2,2,2,2), 8, 1, byrow = TRUE))
              }
              if(nrow.plots == 1){
                layout(matrix(c(1,1,1,1,1,1,1,1), 8, 1, byrow = TRUE))
              }
              par(mar = c(0.2,5,0.2,5))


              # Calc axis limits.
              ylim.flow <- c(0, ceiling(max( c(max(data$flow,na.rm=T),max(flow.viterbi.est,na.rm=T)))))
              ylim.precip.max = max(ylim.flow)*3
              ylim.precip <- c(-ylim.precip.max,0)

              # Setup year range for plotitng
              if (length(plot.yearRange)==2 && all(is.numeric(plot.yearRange)) && all(plot.yearRange>0) && all(plot.yearRange<=as.numeric(format(Sys.Date(), "%Y")))) {
                xlim = as.Date(c(ISOdate(plot.yearRange[1],1,1), ISOdate(plot.yearRange[2],12,31)))
              } else {
                xlim = c(min(obsDates.asISO), max(obsDates.asISO))
              }
              # message(paste("xlim = ", xlim,sep=""))

              # if seasonal observations, adjust plot type...
              if(plot.units == "month"){
                if(with(rle(data$year), max(lengths)) > 4){
                  plot.type = "s"
                }else{
                  plot.type = "p"
                }
              }else{
                plot.type = "s"
              }

              if("A" %in% plot.options){
                if(tail(plot.options, n=1) =="A"){
                  par(mar = c(4,5,0.2,5))
                }
                # Plot obs precip

                pframe = padr::pad(data.frame(obsDates.asISO.withNAs, .Object@input.data$precipitation), interval = plot.units)
                if(plot.units == "month"){ #if seasonal... adjust to get plot with connecting lines..
                  if(with(rle(data$year), max(lengths)) <= 4){
                    pframe$.Object.input.data.precipitation = zoo::na.approx(object = replace(pframe$.Object.input.data.precipitation, is.na(pframe$.Object.input.data.precipitation), NA), maxgap = 2)
                  }
                }

                plot(pframe, type="s",col='grey', lwd=1,
                     xlim=xlim, xlab='', ylab='', main='', xaxt='n')
                mtext("Precip.",side=2,line=3, cex = 0.7)
                mtext(paste("[mm/",plot.units,"]",sep=''),side=2,line=2, cex=0.6)

                xaxis.ticks = as.Date(ISOdate(seq(1900,2020,by=10),1,1))
                abline(v=xaxis.ticks, col = "lightgray", lty = "dotted",lwd = par("lwd"))
                if(tail(plot.options, n=1) =="A"){
                  axis(1,at= xaxis.ticks,labels=seq(1900,2020,by=10))
                  mtext('Year',side=1,line=3)
                }
                grid(NA,NULL)
                plot.range=par("usr")
                text(plot.range[1]+diff(plot.range[1:2])*0.025, plot.range[3]+diff(par("usr")[3:4])*0.95, labels="", font=1, cex=2,pos=1)
              }

              if("B" %in% plot.options){
                if(tail(plot.options, n=1) =="B"){
                  par(mar = c(4,5,0.2,5))
                }

                # Plot obs flow
                pframe = padr::pad(data.frame(obsDates.asISO.withNAs, data.withNAs$flow), interval = plot.units)
                if(plot.units == "month"){ #if seasonal... adjust to get plot with connecting lines..
                  if(with(rle(data$year), max(lengths)) <= 4){
                    pframe$data.withNAs.flow = zoo::na.approx(object = replace(pframe$data.withNAs.flow, is.na(pframe$data.withNAs.flow), NA), maxgap = 2)
                  }
                }

                plot(pframe, type="l", col='grey', lwd=1, ylim=ylim.flow, xlim=xlim,
                     xlab='', ylab='',xaxt='n')
                # ryticks.min = signif(1.2*max(data$precipitation,na.rm=T),1)
                # ryticks = seq(-ryticks.min, 0, by=signif(ryticks.min/2,1))
                # print(-1*(data$precipitation))
                # twoord.plot(obsDates.asISO, data$flow,
                #             obsDates.asISO, -1*(data$precipitation),
                #             lylim=ylim.flow,rylim= ylim.precip,type=c("l","s"),
                #             xlab="Year", ylab=paste("Flow [mm/",plot.units,"]",sep=''),rylab=paste("Precip [mm/",plot.units,"]",sep=''),
                #             lytickpos=seq(0,max(ylim.flow), by=signif(max(ylim.flow)/5,1)), rytickpos=ryticks,
                #             lcol='black',rcol='blue')



                # Plot Markov states as boxes
                for (i in 1:nQhat) {
                  points(obsDates.asISO[i], flow.viterbi.est[i,2],col=state.colours[viterbiPath[i]], bg=state.colours[viterbiPath[i]], pch=21)
                  lines(rep(obsDates.asISO[i],2), flow.line.matrix[i,],col=state.colours[viterbiPath[i]], lwd=1)

                  # Plot the normal flow in years when the Viterbi state is not normal.
                  if (!ind.flow.normal[i]) {
                    points(obsDates.asISO[i], flow.normal.est[i,2],col=state.colours.all[3], bg=state.colours.all[3], pch=1)
                    lines(rep(obsDates.asISO[i],2), c(flow.normal.est[i,1],flow.normal.est[i,3]),col=state.colours.all[3], lwd=1, lty=1)
                  }
                }


                # Add axis labels and legend
                mtext("Flow",side=2,line=3, cex = 0.7)
                mtext(paste("[mm/",plot.units,"]"),side=2,line=2, cex=0.6)
                #mtext('Year',side=1,line=2)
                #legend('topleft', legend=.Object@state.labels, pch=21, col=state.colours, pt.bg=state.colours, xjust=0)
                xaxis.ticks = as.Date(ISOdate(seq(1900,2020,by=10),1,1))
                abline(v=xaxis.ticks, col = "lightgray", lty = "dotted",lwd = par("lwd"))
                if(tail(plot.options, n=1) =="B"){
                  axis(1,at= xaxis.ticks,labels=seq(1900,2020,by=10))
                  mtext('Year',side=1,line=3)
                }
                grid(NA,NULL)

                if(tail(plot.options, n=1) =="B"){
                  legend('topleft', legend=c('Obs. flow',paste(.Object@state.labels,' (5th - 50th - 95th)',sep=''),'Est. normal (5th - 50th - 95th)'),
                         lty=c(1,1,1,1),lwd=1,pch=c(NA,21,21,1), col=c('grey',state.colours,state.colours.all[3]),
                         pt.bg=c(NA,state.colours,NA), xjust=0, cex=1.3, bg='transparent')
                }else{
                  legend('topleft', legend=c('Obs. flow',paste(.Object@state.labels,' (5th - 50th - 95th)',sep=''),'Est. normal (5th - 50th - 95th)'),
                         lty=c(1,1,1,1),lwd=1,pch=c(NA,21,21,1), col=c('grey',state.colours,state.colours.all[3]),
                         pt.bg=c(NA,state.colours,NA), xjust=0, cex=1.1, bg='transparent')
                }

                plot.range=par("usr")
                text(plot.range[1]+diff(plot.range[1:2])*0.025, plot.range[3]+diff(par("usr")[3:4])*0.95, labels="", font=1, cex=2,pos=1)
              }

              if("C" %in% plot.options){
                if(tail(plot.options, n=1) =="C"){
                  par(mar = c(4,5,0.2,5))
                }
                # Calc axis limits.
                ylim.qhat <- c(floor(min( c(min(data$Qhat.flow,na.rm=T),min(viterbi.est,na.rm=T)))) ,
                               ceiling(max( c(max(data$Qhat.flow,na.rm=T),max(viterbi.est,na.rm=T)))))

                # Plot obs Qhat
                pframe = padr::pad(data.frame(obsDates.asISO.withNAs, data.tmp$Qhat.flow), interval = plot.units)
                if(plot.units == "month"){ #if seasonal... adjust to get plot with connecting lines..
                  if(with(rle(data$year), max(lengths)) <= 4){
                    pframe$data.tmp.Qhat.flow = zoo::na.approx(object = replace(pframe$data.tmp.Qhat.flow, is.na(pframe$data.tmp.Qhat.flow), NA), maxgap = 2)
                  }
                }

                plot(pframe, type="l",col='grey', lwd=1,
                     ylim=ylim.qhat, xlim=xlim, xlab='', ylab='',xaxt='n')

                # Plot Markov states as boxes
                for (i in 1:nQhat) {
                  points(obsDates.asISO[i], viterbi.est[i,2],col=state.colours[viterbiPath[i]], bg=state.colours[viterbiPath[i]], pch=21)
                  lines(rep(obsDates.asISO[i],2), line.matrix[i,],col=state.colours[viterbiPath[i]], lwd=1)
                }

                # Add axis labels
                mtext("Transformed flow",side=2,line=3, cex = 0.7)
                mtext(paste("[f(mm) /",plot.units,"]"),side=2,line=2, cex=0.6)
                # mtext('Year',side=1,line=2)
                xaxis.ticks = as.Date(ISOdate(seq(1900,2020,by=10),1,1))
                abline(v=xaxis.ticks, col = "lightgray", lty = "dotted",lwd = par("lwd"))
                if(tail(plot.options, n=1) =="C"){
                  axis(1,at= xaxis.ticks,labels=seq(1900,2020,by=10))
                  mtext('Year',side=1,line=3)
                }
                grid(NA,NULL)

                if(tail(plot.options, n=1) =="C"){
                  legend('bottomleft', legend=c('Obs. flow',paste(.Object@state.labels,' (5th - 50th - 95th)',sep=''),'Est. normal (5th - 50th - 95th)'),
                         lty=c(1,1,1,1),lwd=1,pch=c(NA,21,21,1), col=c('grey',state.colours,state.colours.all[3]),
                         pt.bg=c(NA,state.colours,NA), xjust=0, cex=1.3, bg='transparent')
                }

                plot.range=par("usr")
                text(plot.range[1]+diff(plot.range[1:2])*0.025, plot.range[3]+diff(par("usr")[3:4])*0.95, labels="", font=1, cex=2,pos=1)
              }

              # Plot the cummulative rainfall residual.
              #--------------

              # # Calculate the means and residuals
              # if (plot.units == 'yr') {
              #   P.mean = mean(data$precipitation)
              #   P.resid = data$precipitation - P.mean;
              # } else {
              #   P.mean = rep(NA,length(.Object@QhatModel.object@subAnnual.Monthly.Steps))
              #   P.resid = rep(NA, length(data$precipitation))
              #   for (i in 1:length(.Object@QhatModel.object@subAnnual.Monthly.Steps)) {
              #     filt = data$month == .Object@QhatModel.object@subAnnual.Monthly.Steps[i]
              #     P.mean[i] =  mean(data$precipitation[filt])
              #     P.resid[filt] = data$precipitation[filt] - P.mean[i]
              #   }
              # }
              #
              # # # Plot the residuals
              # # plot(obsDates.asISO, P.resid,type='p',xlim=xlim, col='white', bg='white', pch=21, xlab='', ylab='')
              # # for (i in 1:nQhat) {
              # #   points(obsDates.asISO[i], P.resid[i],col=state.colours[viterbiPath[i]], bg=state.colours[viterbiPath[i]], pch=21)
              # # }
              # # lines(obsDates.asISO, rep(0,length(obsDates.asISO)),col='grey')
              # # mtext("Rainfall residual [mm]",side=2,line=3)
              # # mtext('Year',side=1,line=2)
              #
              # # Calculate the cumulative residuals.
              # P.cumResid = cumsum(P.resid)
              #
              # # Plot the cumm residuals
              # plot(obsDates.asISO, P.cumResid, type='l',col='grey', lwd=1, xlim=xlim, xlab='', ylab='',xaxt='n')
              # abline(v=xaxis.ticks, col = "lightgray", lty = "dotted",lwd = par("lwd"))
              # grid(NA,NULL)
              # legend('topright', legend=c('Cum. residual ',.Object@state.labels),
              #        lty=c(1,NA,NA),pch=c(NA,21,21), col=c('grey',state.colours),
              #        pt.bg=c(NA,state.colours), xjust=0, cex=1.5, bg='white')
              #
              # # Colour the points by the Viterbi state.
              # for (i in 1:nQhat) {
              #   points(obsDates.asISO[i], P.cumResid[i],col=state.colours[viterbiPath[i]], bg=state.colours[viterbiPath[i]], pch=21)
              # }
              # mtext("Cum. rainfall resid.",side=2,line=3)
              # mtext(paste("[mm]"),side=2,line=2, cex=0.85)
              # plot.range=par("usr")
              # text(plot.range[1]+diff(plot.range[1:2])*0.025, plot.range[3]+diff(par("usr")[3:4])*0.95, labels="D", font=1, cex=2,pos=1)

              if("D" %in% plot.options){
                if(tail(plot.options, n=1) =="D"){
                  par(mar = c(4,5,0.2,5))
                }
                # Plot the conditional state probability for each state
                # if (nStates>1) {
                emissionProbs = getEmissionDensity(.Object@QhatModel.object, data, NA)
                state.probs = getConditionalStateProbabilities(.Object@markov.model.object, data, emissionProbs)

                # Plot bar graph
                # par(mar = c(4,5,0.2,5))
                plot(obsDates.asISO, state.probs[1,], type = 'l', xlim=xlim, col=state.colours[1],
                     ylim=c(0,1),xlab='', ylab='', lwd=1, xaxt='n')
                if (nStates>1) {
                  for ( i in 2:nStates) {
                    lines(obsDates.asISO, state.probs[i,], col=state.colours[i])
                  }
                }
                mtext("State Prob.",side=2,line=3, cex = 0.7)
                mtext(paste("[-]"),side=2,line=2, cex=0.6)

                xaxis.ticks = as.Date(ISOdate(seq(1900,2020,by=10),1,1))
                if(tail(plot.options, n=1) =="D"){
                  axis(1,at= xaxis.ticks,labels=seq(1900,2020,by=10))
                  mtext('Year',side=1,line=3)
                }
                abline(v=xaxis.ticks, col = "lightgray", lty = "dotted",lwd = par("lwd"))
                grid(NA,NULL)
                legend('bottomleft', legend=.Object@state.labels,
                       lty=c(1,1),pch=c(NA,NA), col=state.colours,lwd=1,
                       xjust=0, cex=1.1, bg='white')
                plot.range=par("usr")
                text(plot.range[1]+diff(plot.range[1:2])*0.025, plot.range[3]+diff(par("usr")[3:4])*0.95, labels="", font=1, cex=2,pos=1)
              }

              #}

              # Reset graphics options
              #--------------
              par(op)

            }

            # Return data
            if(do.plot){
              return(invisible())
            }else{
              return(as.data.frame(results))
            }

          }
)

# @exportMethod check.viterbi
setGeneric(name="check.viterbi",def=function(.Object, nSamples=NA) {standardGeneric("check.viterbi")})
setMethod(f="check.viterbi",signature="hydroState",definition=function(.Object, nSamples=100000)
{
  if (!validObject(.Object))
    stop('The model parameters produced an INVALID MODEL.')

  # Duplicate the inout data 100 times.
  data = .Object@input.data
  while (nrow(data)<nSamples) {
    date.tmp <- .Object@input.data
    date.tmp$year <- date.tmp$year + (max(data$year)-min(data$year))+1
    data <- rbind.data.frame(data,date.tmp)
  }

  # Get the transformed precipitation.
  data = getQhat(.Object@Qhat.object, data)

  # Generate a synthtic series from the HMM.
  states.sample = generate.sample.states(.Object@markov.model.object, data)

  # Get a time series of Qhat using the sample states
  Qhat.sample = generate.sample.Qhat.fromViterbi(.Object@QhatModel.object,data, states.sample)

  data$Qhat.flow <- NULL
  data = cbind.data.frame(data, Qhat.flow = Qhat.sample)

  # Find the Viterbi states for the sample Qhat data
  states.viterbi = viterbi(.Object, data, do.plot = F, plot.percentiles=c(0.05, 0.5, 0.95),plot.yearRange=numeric())

  # Filer out NAs
  states.sample = states.sample[is.finite(Qhat.sample)]
  data = data[is.finite(Qhat.sample),]

  # Assess probability that the inferred state equals the 'known' state.
  nStates = getNumStates(.Object@markov.model.object)
  nsamples = nrow(data)
  Pr = matrix(Inf, nStates,nStates)
  for (i in 1:nStates) {
    filt = states.sample == i
    states.viterbi.filt = states.viterbi[filt,'Viterbi State Number']
    for (j in 1:nStates) {
      Pr[i,j] = sum(states.viterbi.filt==j,na.rm =T)/sum(filt)
    }
  }

  Pr.df = data.frame(Pr)
  state.names = .Object@state.labels
  if (length(.Object@state.labels)==0 || nchar(state.names)==0) {
    names(Pr.df) = paste("Inferred State",1:nStates)
    row.names(Pr.df) = paste("Known State",1:nStates)
  } else {
    names(Pr.df) = paste("Inferred State",state.names)
    row.names(Pr.df) = paste("Known State",state.names)

  }

  return(Pr.df)
}
)


# @exportMethod check.PseudoResiduals
setGeneric(name="check.PseudoResiduals",def=function(.Object, nIncrements=20, do.plot=T) {standardGeneric("check.PseudoResiduals")})
setMethod(f="check.PseudoResiduals",signature="hydroState",definition=function(.Object, nIncrements, do.plot)
{

  if (!validObject(.Object))
    stop('The model parameters produced an INVALID MODEL.')

  # This function derives the uniform and normal psedu residuals.
  # The approach is based upon Zucchini, McDonald and Langrock, 2016 then section 6.3.2 P 108, 2nd edition
  # but modified to allow fow a time varying markov state mean. This required the
  # use of the Viterbi state to determine which value to take in the conditional
  # distribtion.

    # Gte Qhat and add to input.data
    data  = .Object@input.data
    n = nrow(data)
    data = getQhat(.Object@Qhat.object, .Object@input.data)
    Qhat <- data$Qhat.flow
    # data = cbind.data.frame(.Object@input.data,Qhat)
    # get number of states
    nStates = getNumStates(.Object@markov.model.object)

    # Built filter for non NAs.
    filt <- !is.na(data$Qhat.flow)

    # get emission densities
    emissionDensity <- getEmissionDensity(.Object@QhatModel.object, data, NA)
    emissionDensity[!filt,] <- NA

    # Set the range in Qhat values at which to derive the conditional probs.
    Qhat.increments = seq(floor(min(Qhat[filt])),ceiling(max(Qhat[filt])),length.out=100)

    # Get the emmision cumulative prob (not density) at each Qhat.increments value at each tiem step
    cumProb.increments = array(NA, dim=c(n,nStates, length(Qhat.increments)))
    for (j in 1:length(Qhat.increments)) {
        cumProb.increments[,,j] <- getEmissionDensity(.Object@QhatModel.object, data, cumProb.threshold.Qhat= rep(Qhat.increments[j], nrow(data)))
      }


    # # Get the emission cumulative probs and sort for each state.
    # emissionCumProbs <- getEmissionDensity(.Object@QhatModel.object, data, getCumProb=T)
    # cumProb.increments = matrix(0,sum(filt),nStates)
    # for (i  in 1:nStates) {
    #   cumProb.increments[,i] = sort(emissionCumProbs[filt,i])
    # }
    #
    # # Identify the most likely emission cumulative probs over time using the viterbi states.
    # # This is used in the latter selection of a row within the conditional distribution.
    # viterbi.states <- viterbi(.Object, do.plot=F, plot.percentiles = c(0.05, 0.5, 0.95),plot.yearRange=numeric())
    # viterbi.emissionCumProbs = rep(0,nrow(data))
    # for (i in 1:nrow(data)) {
    #   if (all(is.finite(emissionCumProbs[i,])))
    #     viterbi.emissionCumProbs[i] <- emissionCumProbs[i,viterbi.states[i,'Viterbi State Number']]
    # }

    # Adapted from Zucchini, McDonald and Langrock, 2016, Hidden Markov Models for Time Series, 2nd Edision, Appendix A.
    #---------------------------------------------------------------------
    # NOTE, unlilke Zucchini, here the cumulative IS NOT USED. This is because cumProb.increments is
    # already the cumulative prob.

    cum.prob <- getConditionalProbabilities(.Object@markov.model.object, data, emissionDensity, cumProb.increments)
    #cumdists <- rbind(rep(0,n), cdists)
    #cumProb.increments <- rbind(rep(0,nStates), cumProb.increments)

    P.est <- rep(NA,n)
    for (i in 1:n)
    {
      if (!filt[i])
        next

      # Interpolate the conditional cumulative distribution from the discrete values to Qhat at time point i
      # NOTE, to avoid P.est[i]==0 (and hence -inf pseudo normal residuals), data$Qhat.flow[i] is limit to >0
      # machine preision.
      x = max(data$Qhat.flow[i], sqrt(.Machine$double.eps))
      ind <- which(Qhat.increments>=x)[1]

      if (ind==1) {
        P.est[i] <- cum.prob[1,i]
      } else if (ind<length(Qhat.increments)) {
        P.lo  <- cum.prob[ind-1,i]
        P.hi  <- cum.prob[ind,i]
        P.est[i] <- (P.lo + (P.hi - P.lo)/(Qhat.increments[ind] - Qhat.increments[ind-1])*(x - Qhat.increments[ind-1]))
      } else{
        P.est[i] <- cum.prob[nrow(cum.prob),i]
      }
      # Limit to <0.999 to avoid problems withj acf().
      P.est[i] <- min(c(0.999,P.est[i]))

    }
    npsr     <- qnorm(P.est)
    data <- cbind.data.frame(data, unif.pseuou.resid=P.est, norm.pseudo.resid=npsr)
    #
    #---------------------------------------------------------------------

    # Plot
    if (do.plot) {

      # Handle monthly and yearly time steps
      if (any(names(data)=="month")) {
        obsDates = ISOdate(data$year,data$month,1)
      } else {
        obsDates = data$year
      }

      # Get input grapics settings
      op <- par(no.readonly = T);

      plot.titles=F

      # Change graphics settings
      #par(mar = rep(4, 4))
      #par(mfrow=c(3,2),cex=1.2)
      #par(mfrow=c(3,2), cex.axis=1.2, cex.lab=1.2)
      par(mfrow=c(5,1), cex.axis=1.2, cex.lab=1.2)

      # Plot time series of pseduo normal resdiuals
      filt.isfinite = is.finite(npsr)
      ylim = c( min(c(npsr[filt.isfinite],-3),na.rm = T), max(c(npsr[filt.isfinite],3),na.rm = T))
      if (plot.titles) {
        plot(obsDates[filt.isfinite],npsr[filt.isfinite], type='p', xlab='Date',ylab="Normal-pseudo resid.", ylim=ylim,main='(A) Time-series of normal-pseudo residuals')
      } else {
        plot(obsDates[filt.isfinite],npsr[filt.isfinite], type='p', xlab='Date',ylab="Normal-pseudo resid.", ylim=ylim)
      }
      abline(0, 0, col='blue')
      abline(qnorm(0.975), 0, col='blue', lwd=0.5, lty=2)
      abline(qnorm(0.995), 0, col='blue', lwd=0.5, lty=2)
      abline(-qnorm(0.975), 0, col='blue', lwd=0.5, lty=2)
      abline(-qnorm(0.995), 0, col='blue', lwd=0.5, lty=2)
      if (!plot.titles) {
        plot.range=par("usr")
        text(plot.range[1]+diff(plot.range[1:2])*0.05, plot.range[3]+diff(par("usr")[3:4])*0.95, labels="A", font=1, cex=2,pos=1)
      }

      # plot auto correlation funtion
      if (plot.titles) {
        acf(npsr, na.action=na.pass, xlab='Normal-pseudo resid.', main='(B) Auto-correlation of normal-pseudo residuals')
      } else {
        acf(npsr[filt.isfinite], na.action=na.pass, xlab='Normal-pseudo resid.', main='')
        plot.range=par("usr")
        text(plot.range[1]+diff(plot.range[1:2])*0.1, plot.range[3]+diff(par("usr")[3:4])*0.95, labels="B", font=1, cex=2,pos=1)
      }

      # Plot histogram of uniform pseduo residuals
      if (plot.titles) {
        hist(P.est[filt.isfinite], breaks=nIncrements,freq=F, plot=T, xlab="Uniform-pseudo resid.", ylab='Freq.', xlim=c(0,1),main='(C) Histogram of uniform-pseudo residuals')
      } else {
        hist(P.est[filt.isfinite], breaks=nIncrements,freq=F, plot=T, xlab="Uniform-pseudo resid.", ylab='Freq.', xlim=c(0,1),main='')
        plot.range=par("usr")
        text(plot.range[1]+diff(plot.range[1:2])*0.05, plot.range[3]+diff(par("usr")[3:4])*0.95, labels="C", font=1, cex=2,pos=1)
      }
      abline(1, 0, col='blue', lwd=0.75, lty=2)

      # Plot histogram of normal pseudo resid.
      x <- seq(-5, 5, length=10000)
      y <- dnorm(x, mean=0, sd=1)
      if (plot.titles) {
        hist(npsr[filt.isfinite], breaks=nIncrements,freq=F, plot=T, xlim=ylim,xlab="Normal-pseudo resid.", ylab='Freq.',main='(D) Histogram of normal-pseudo residuals')
      } else {
        hist(npsr[filt.isfinite], breaks=nIncrements,freq=F, plot=T, xlim=ylim,xlab="Normal-pseudo resid.", ylab='Freq.',main='')
        plot.range=par("usr")
        text(plot.range[1]+diff(plot.range[1:2])*0.05, plot.range[3]+diff(par("usr")[3:4])*0.95, labels="D", font=1, cex=2,pos=1)
      }
      lines(x, y,  ylim=c(0, max(c(y,npsr))*1.05), col='blue', lwd=0.75, lty=2)

      # # Plot qq plot of uniform pseduo residuals
      # qqplot(x=runif(length(P.est[filt.isfinite])),y=P.est[filt.isfinite],xlim=c(0,1),ylim=c(0,1),xlab = "Theoretical Quantiles",ylab="Uniform-pseudo resid.", main='(E) Q-Q plot of uniform-pseudo residuals')
      # #qqplot(x=runif(length(P.est[filt.isfinite])),y=P.est[filt.isfinite],xlim=c(0,1),ylim=c(0,1),xlab = "Theoretical Quantiles",ylab="Uniform-pseudo resid.", main='')
      # qqline(P.est[filt.isfinite], distribution = qunif)
      # plot.range=par("usr")
      # #text(plot.range[1]+diff(plot.range[1:2])*0.05, plot.range[3]+diff(par("usr")[3:4])*0.95, labels="J", font=1, cex=2,pos=1)

      # Plot histogram of normal pseudo resid.
      if (plot.titles) {
        qqnorm(npsr[filt.isfinite],ylim=ylim, xlab = "Theoretical Quantiles",ylab="Normal-pseudo resid.", main='(E) Q-Q plot of normal-pseudo residuals')
      } else {
        qqnorm(npsr[filt.isfinite],ylim=ylim, xlab = "Theoretical Quantiles",ylab="Normal-pseudo resid.", main='')
        plot.range=par("usr")
        text(plot.range[1]+diff(plot.range[1:2])*0.05, plot.range[3]+diff(par("usr")[3:4])*0.95, labels="E", font=1, cex=2,pos=1)
      }
      qqline(npsr, col='blue', lwd=0.75, lty=2)

      # Plot the Shapiro–Wilk p-value (normality test) in lower RHS corner.
      ShapiroWilk.pvalue = shapiro.test(npsr[filt.isfinite])$p.value
      plot.range=par("usr")
      text(plot.range[1]+diff(plot.range[1:2])*0.99, plot.range[3]+diff(par("usr")[3:4])*0.1, labels=paste("Shapiro-Wilk p-value=",round(ShapiroWilk.pvalue,3)), font=1, cex=1.5,pos=2)

      # Plot AIC in lower RHS corner.
      AIC = getAIC(.Object)
      plot.range=par("usr")
      text(plot.range[1]+diff(plot.range[1:2])*0.99, plot.range[3]+diff(par("usr")[3:4])*0.225, labels=paste("AIC=",round(AIC,2)), font=1, cex=1.5,pos=2)

      # Reset graphics options
      par(op)
    }

    if(do.plot){
      return(invisible())
    }else{
      return(data)
    }
}
)

# @exportMethod drought.resilience.index
setGeneric(name="drought.resilience.index",def=function(.Object, year.drought.start, year.drought.end,year.postdrought.end) {standardGeneric("drought.resilience.index")})
setMethod(f="drought.resilience.index",signature="hydroState",definition=function(.Object, year.drought.start, year.drought.end,year.postdrought.end)
{

  if (getNumStates(.Object@markov.model.object)==1)
    return(NA)

  # Get the viterbi states.
  states.viterbi <- viterbi(.Object, do.plot=F)

  # Get the transformed flow, Qhat
  data  = .Object@input.data
  data = getQhat(.Object@Qhat.object, data)

  # Generate 10,000 samples of Qhat at each time point and state.
  nSamples=10000
  Qhat.sample = generate.sample.Qhat(.Object@QhatModel.object,data, nSamples=nSamples)

  # Get viterbi years years during and after the user-defined drought
  filt.years.drought = states.viterbi[,1]>=year.drought.start & states.viterbi[,1]<=year.drought.end
  filt.years.postDrought = states.viterbi[,1]>year.drought.end & states.viterbi[,1]<=year.postdrought.end

  # Build indexes to each named state
  ind.verylow = which(.Object@state.labels=='Very low')
  ind.low = which(.Object@state.labels=='Low')
  ind.normal = which(.Object@state.labels=='Normal')
  ind.normal_dup = which(.Object@state.labels=='Normal (duplicate)')
  ind.high = which(.Object@state.labels=='High')
  ind.veryhigh = which(.Object@state.labels=='Very high')

  # Get the fraction of time the catchment is within a flow state lower than normal during the drought.
  if (length(ind.verylow)>0){
    filt.lowStates.drought = (states.viterbi[,'Viterbi State Number'] == ind.low | states.viterbi[,'Viterbi State Number'] == ind.verylow) &
      filt.years.drought
  } else {
    filt.lowStates.drought = states.viterbi[,'Viterbi State Number'] == ind.low & filt.years.drought
  }
  f.drought = sum(filt.lowStates.drought,na.rm =T)/sum(filt.years.drought,na.rm =T)

  # Get the fraction of time the catchment is within a flow state lower than normal after the drought.
  if (length(ind.verylow)>0){
    filt.lowStates.postDrought = (states.viterbi[,'Viterbi State Number'] == ind.low | states.viterbi[,'Viterbi State Number'] == ind.verylow) &
      filt.years.postDrought
  } else {
    filt.lowStates.postDrought = states.viterbi[,'Viterbi State Number'] == ind.low & filt.years.postDrought
  }
  f.postdrought = sum(filt.lowStates.postDrought,na.rm =T)/sum(filt.years.postDrought,na.rm =T)

  # For all low & very low flow drought year, collate the sampled normal and subnormal flow and the calculate the probability that the flow was below normal.
  Qhat.sample.normalState.drought = numeric()
  Qhat.sample.lowStates.drought = numeric()
  for (i in which(filt.lowStates.drought)) {
   if (states.viterbi[i,'Viterbi State Number']==ind.verylow || states.viterbi[i,'Viterbi State Number']==ind.low) {
     Qhat.sample.lowStates.drought = c(Qhat.sample.lowStates.drought, Qhat.sample[[states.viterbi[i,'Viterbi State Number']]][i,])
     Qhat.sample.normalState.drought = c(Qhat.sample.normalState.drought, Qhat.sample[[ind.normal]][i,])
   }
  }
  if (length(Qhat.sample.normalState.drought)==0) {
    P.droughtFlow.le.normal=0
  } else {
    P.droughtFlow.le.normal = sum(Qhat.sample.lowStates.drought<Qhat.sample.normalState.drought,na.rm =T)/length(Qhat.sample.normalState.drought)
  }

  # For all low & very low flow post-drought year, collate the sampled normal and subnormal flow and the calculate the probability that the flow was below normal.
  Qhat.sample.normalState.postDrought = numeric()
  Qhat.sample.lowStates.postDrought = numeric()
  for (i in which(filt.lowStates.postDrought)) {
    if (states.viterbi[i,'Viterbi State Number']==ind.verylow || states.viterbi[i,'Viterbi State Number']==ind.low) {
      Qhat.sample.lowStates.postDrought = c(Qhat.sample.lowStates.postDrought, Qhat.sample[[states.viterbi[i,'Viterbi State Number']]][i,])
      Qhat.sample.normalState.postDrought = c(Qhat.sample.normalState.postDrought, Qhat.sample[[ind.normal]][i,])
    }
  }
  if (length(Qhat.sample.normalState.postDrought)==0) {
    P.postDroughtFlow.le.normal=0
  } else {
    P.postDroughtFlow.le.normal = sum(Qhat.sample.lowStates.postDrought<Qhat.sample.normalState.postDrought,na.rm =T)/length(Qhat.sample.normalState.postDrought)
  }


  # Finally calculae resilience index.
  R <- sqrt((1-f.drought*P.droughtFlow.le.normal)^2 + (1-f.postdrought*P.postDroughtFlow.le.normal)^2)/sqrt(2)

  # Build matrix of componants
  R = data.frame(frac.drought = f.drought, Prob_Qdrought_LE_Qnormal = P.droughtFlow.le.normal, frac.postdrought = f.postdrought,  Prob_Qpostdrought_LE_Qnormal = P.postDroughtFlow.le.normal, resilience.index = R)

  return(R)
}
)

# @exportMethod forecast
setGeneric(name="forecast",def=function(.Object, t) {standardGeneric("forecast")})
setMethod(f="forecast",signature="hydroState",definition=function(.Object, t)
  {


  }
)
