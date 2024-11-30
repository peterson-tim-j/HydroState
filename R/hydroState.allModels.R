##' @include hydroState.R
## @export
hydroState.allModels <- setClass(
  # Set the name for the class
  "hydroState.allModels",

  package='hydroState',

  # contains=c('hydroState'),

  # Define the slots
  slots = c(
    siteID= 'character',
    calib.reference.model.name = 'list',
    calib.reference.criteria.met = 'logical',
    models = 'list'

  ),

  # Set the default values for the slots. (optional)
  prototype=list(
    siteID = '(not set)',
    calib.reference.model.name= vector('list',1),
    calib.reference.criteria.met= vector('logical',1),
    models  = vector('list',1)
  )
)

# Valid object?
validObject <- function(object) {
  TRUE
}
setValidity("hydroState.allModels", validObject)

# Initialise the object.
#setGeneric(name="initialize",def=function(.Object,input.data, Qhat.object, QhatModel.object, markov.model.object, ...){standardGeneric("initialize")})
setMethod(f="initialize",signature="hydroState.allModels",definition=function(.Object, models, siteID)
{

  .Object@siteID <- siteID

  # Build hydrostate models
  .Object@models <- unlist(models)

  # Define Refernce models. That is, the calibration obecjive function that this models needs to meet or exceed.
  .Object@calib.reference.model.name <- names(models)

  model.names = names(.Object@models)
  .Object@calib.reference.criteria.met = rep(F,length(model.names))
  names(.Object@calib.reference.criteria.met) <- model.names

  return(.Object)

}

)

setMethod(f = "fit",signature="hydroState.allModels",definition=function(.Object,
                                                                         pop.size.perParameter=25,
                                                                         max.generations=10000,
                                                                         Domains,
                                                                         reltol=1e-8,
                                                                         steptol=50,
                                                                         print.iterations = 25,
                                                                         use.initial.parameters=F,
                                                                         doParallel=F,
                                                                         ...)
{

  # Set the min and max numbe rof calibrations per model
  minTrialsPerModel=5
  maxTrialsPerModel=20

  # if(is.NUll(doParallel))
  #   doParallel == F


  model.names <- names(.Object@models)
  nModels = length(.Object@models)

  message(paste('Calibrating ',nModels,'Models.'))
  message(paste('... Minimum number of calibration per model: ',minTrialsPerModel))
  message(paste('... Maximum number of calibration per model: ',maxTrialsPerModel))

  for (i in 1:nModels) {

    # Get reference model calibration result.
    current.model.name <- model.names[i]
    ref.model.name <- .Object@calib.reference.model.name[[current.model.name]]
    if (nchar(ref.model.name)==0){
      max.calib.result=-Inf
    } else {
        ref.model.obj <- .Object@models[[ref.model.name]]
        max.calib.result <- ref.model.obj@calibration.results$bestval
    }

    # Print user summary
    message('************************************')
    message(paste('Calibrating Model:',current.model.name))
    message(paste('Reference model name for maximum acceptable calib. solution:',ref.model.name))
    message(paste('Reference model maximum acceptable calib. solution:',max.calib.result))
    message('')

    current.model.err = Inf
    j=1
    while (j<minTrialsPerModel || (j<maxTrialsPerModel && current.model.err>max.calib.result)){
      DEstrategy <- sample(c(2,3),1,replace=T)
      message('')
      message(paste('... Doing Calibration Trial ',j, ' using DEstrategy:',DEstrategy))

      # Get number of params
      nParams = length(getParameters.asVector(.Object@models[[i]]))

      # hydroState b.c. ran parallel based on number of params... now just do parallel if user assigns it...
      # if (nParams>=doParallel || doParallel==T) {
      # print(doParallel)

      if (doParallel==T) {
        model <- .Object@models[[i]]
        assign("model", model, envir=globalenv())
                 #new.env(parent = baseenv()))

        # error occurs after calling this...
        model <- fit(model, DEstrategy=DEstrategy, pop.size.perParameter=pop.size.perParameter, max.generations=max.generations, Domains=Domains, reltol=reltol, steptol=steptol, print.iterations=print.iterations,
                     use.initial.parameters = use.initial.parameters, parallelType= "auto", packages = list('hydroState','truncnorm'),parVar=list('model'))
      } else {
        model <- fit(.Object@models[[i]], DEstrategy=DEstrategy, pop.size.perParameter=pop.size.perParameter, max.generations=max.generations, reltol=reltol, steptol=steptol, print.iterations=print.iterations,
                     use.initial.parameters = use.initial.parameters,  doParallel = F, ...)
      }

      # Store the best model to date int the object in case the reference maximum obj cannot be met.
      if (model@calibration.results$bestval < current.model.err) {
        .Object@models[[i]] <- model
        current.model.err <- model@calibration.results$bestval
      }

      j <- j+1
    }

    if (current.model.err<=max.calib.result) {
      .Object@calib.reference.criteria.met[[current.model.name]] = T
      message('Best solution was <= reference model calibration solution!')
    } else {
      warning('Best solution was > reference model calibration solution!')
    }
    message('************************************')

  }
  return(.Object)
}
)

setMethod(f="getNegLogLikelihood",signature="hydroState.allModels",definition=function(.Object)
{
  return(getNegLogLikelihood(.Object,plot.data=T))
}
)
setMethod(f="getNegLogLikelihood",signature="hydroState.allModels",definition=function(.Object, plot.data)
{
  # Initialise outputs
  model.names <- names(.Object@models)
  nModels <- length(model.names)
  nLL <- rep(Inf, nModels)
  names(nLL) <- model.names

  # Get AIC for each model
  for (i in 1:nModels) {
    nLL[i] <- getNegLogLikelihood(.Object@models[[i]])
  }

  # Plot AIC vs model name
  if (missing(plot.data) || plot.data) {
    par(mar=c(15,2,2,2))
    plot(1:nModels, nLL, type='p',ylab = 'Neg. Log Liklihood')
    axis(1, at=1:nModels, labels=model.names,las=2)
  }

  return(nLL)

}
)

# #
# # Get vector of AICc.
# setMethod(f="getAICc",signature="hydroState.allModels",definition=function(.Object)
# {
#   return(getAICc(.Object,plot.data=T))
# }
# )
# # Get vector of AIC.
# setMethod(f="getAIC",signature="hydroState.allModels",definition=function(.Object)
# {
#   return(getAIC(.Object,plot.data=T))
# }
# )
setMethod(f="getAIC",signature="hydroState.allModels",definition=function(.Object, plot.data, use.sampleSize.penality)
{
  # Initialise outputs
  model.names <- names(.Object@models)
  nModels <- length(model.names)
  AIC <- rep(Inf, nModels)
  names(AIC) <- model.names

  # Get AIC for each model
  if (use.sampleSize.penality) {
    for (i in 1:nModels) {
      AIC[i] <- getAICc(.Object@models[[i]])
    }
  } else {
    for (i in 1:nModels) {
      AIC[i] <- getAIC(.Object@models[[i]])
    }

  }

  # Build vector of model parameters.
  nParameters <- rep(Inf, nModels)
  names(nParameters) <- model.names
  for (i in 1:nModels)
    nParameters[i] <- length(getParameters.asVector(.Object@models[[i]]))

  # Building a matrix of AIC and nParams.
  AIC <- cbind(AIC, nParameters)
  rownames(AIC) <- model.names
  colnames(AIC) <- c('AIC','nParameters')

  # Plot AIC vs model name
  if (plot.data) {
    par(mar=c(15,2,2,2))
    plot(1:nModels, AIC[,1], type='p',ylab = 'AIC')
    axis(1, at=1:nModels, labels=model.names,las=2)
  }

  return(AIC)
}
)


# Get best model (by mimimum AIC)
# @exportMethod getAIC.bestModel
setGeneric(name="getAIC.bestModel",def=function(.Object, use.calib.reference.criteria=NA, use.sampleSize.penality=NA,min.obs.per.state=NA) {standardGeneric("getAIC.bestModel")})
setMethod(f="getAIC.bestModel",signature="hydroState.allModels",definition=function(.Object, use.calib.reference.criteria=T, use.sampleSize.penality=T, min.obs.per.state=3)
{
  # Get AIC
  model.names <- names(.Object@models)
  nModels <- length(model.names)
  AIC <- getAIC(.Object, plot.data=F, use.sampleSize.penality)

  # Build filt to models that met calib. criteria
  filt = rep(T, nModels)
  if (use.calib.reference.criteria) {
    for (i in 1:nModels) {
      current.model.name <- model.names[i]
      filt[i] <- F
      if (.Object@calib.reference.criteria.met[[current.model.name]])
        filt[i] <- T
    }
  }
  AIC <- AIC[filt,]

  # Filter out models that have a non-finite AIC.
  nModels <- nrow(AIC)
  filt = rep(F, nModels)
  for (i in 1:nModels) {
    if(is.finite(AIC[i,1]))
       filt[i]=T
  }
  AIC <- AIC[filt,]

  # Find those models where (by Viterbi decoding) it is on each state >= min.obs.per.state times.
  nModels <- nrow(AIC)
  model.names <- rownames(AIC)
  filt = rep(T, nModels)
  for (i in 1:nModels) {

    current.model.name <- model.names[i]

    # get Viterbi states.
    viterbi.data <- viterbi(.Object@models[[current.model.name]], do.plot=F)

    # Count frequency of each state.
    if (any(table(viterbi.data[,'Viterbi State Number'])<min.obs.per.state)) {
      filt[i] <- F
      message(paste('... The following model has a state with only',min(table(viterbi.data[,'Viterbi State Number'])), ' time points: ',current.model.name ))
    }

  }
  AIC <- AIC[filt,]

  # Find best model by the lowest AIC.
  # If there are >1 models with the best AIC, then take that with the lowest number of parameters.
  AIC <- AIC[order(AIC[,1],AIC[,2]),]
  bestModel.name <- rownames(AIC)[1]

  # Get best 1 state model
  ind = grep("model.1State",rownames(AIC))
  oneStateModel.names = rownames(AIC)[ind]
  ind = which(AIC[,1] == min(AIC[ind,1]))[1]
  bestModel.1State.name <- rownames(AIC)[ind]

  # Caclulate the delta AIC weights. From, Kenneth P. Burnham,  David R. Anderson
  # Model Selection and Multimodel Inference:  A Practical Information-Theoretic Approach
  # Second Edition, section 2.9.
  deltaAIC = AIC[,1] - AIC[1,1]
  AIC.weights = exp(-deltaAIC/2)/sum(exp(-deltaAIC/2))
  AIC = cbind(AIC,AIC.weights)
  colnames(AIC) <- c('AIC','nParameters','AIC.weights')

  # Calculate the Evidence Ratio between the best model and the
  # best one-state model. The approach is from Burnham et al p78.
  # To calculate this, first the best one-state model is found.
  model.names <- rownames(AIC);
  ind = grep('model.1State', model.names, value=F)
  ind = ind[1]
  if (ind>1) {
    EF.oneState = AIC[1,3]/AIC[ind,3]
  } else {
    EF.oneState = 1;
  }

  # Get the best model
  bestModel <- .Object@models[[bestModel.name]]

  # Get the best 1 state model
  bestModel.1State <- .Object@models[[bestModel.1State.name]]

  return(list(model=bestModel, model.best1State = bestModel.1State, AIC=AIC, evidenceRatio.oneState=EF.oneState))
}
)
