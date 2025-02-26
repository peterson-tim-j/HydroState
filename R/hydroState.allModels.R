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
    models = 'list',
    models.summary = 'data.frame'

  ),

  # Set the default values for the slots. (optional)
  prototype=list(
    siteID = '(not set)',
    calib.reference.model.name= vector('list',1),
    calib.reference.criteria.met= vector('logical',1),
    models  = vector('list',1),
    models.summary = data.frame()
  )
)

# Valid object?
validObject <- function(object) {
  TRUE
}
setValidity("hydroState.allModels", validObject)

# Initialise the object.
# setGeneric(name="initialize",def=function(.Object, ...){standardGeneric("initialize")})
setMethod(f="initialize",signature="hydroState.allModels",definition=function(.Object, models, siteID)
  {


  .Object@siteID <- siteID

  # assign hydroState all models
  .Object@models <- models

  return(.Object)

}
)


# get summary of all models.
# @exportMethod get.summary.table
setGeneric(name="get.summary.table",def=function(.Object, models.summary = data.frame()){standardGeneric("get.summary.table")})
setMethod(f="get.summary.table",signature="hydroState.allModels",definition=function(.Object, models.summary)
{

  if(length(models.summary)>0){

    .Object@models.summary = models.summary

  }else{
    # Calc details for summary
    nparams = as.numeric(sapply(.Object@models, function(x) length(getParameters.asVector(x))))

    nstates = as.numeric(sapply(.Object@models, function(x) ncol(x@markov.model.object@transition.graph)))

    state.structure = sapply(.Object@models, function(x) ifelse(all(x@markov.model.object@transition.graph == TRUE),"","S"))


    data.transform = sapply(.Object@models, function(x) ifelse(length(x@Qhat.object@parameters@values)>0, 'boxcox','log'))

    #adjust number of model parameters depending on data transform
    nparams = as.numeric(sapply(1:length(nparams), function(x) ifelse(data.transform[x] == 'boxcox', nparams[x] +1, nparams[x])))


    error.distribution = sapply(.Object@models, function(x) ifelse(x@QhatModel.object@use.truncated.dist == TRUE,'truc.normal',
                                                                   ifelse(grepl('normal',class(x@QhatModel.object)),'normal',
                                                                          ifelse(grepl('gamma',class(x@QhatModel.object)),'gamma',''))))

    auto.correlation = as.numeric(sapply(.Object@models, function(x) ifelse('mean.AR3' %in% names(x@QhatModel.object@parameters@values), '3',
                                                                            ifelse('mean.AR2' %in% names(x@QhatModel.object@parameters@values), '2',
                                                                                   ifelse('mean.AR1' %in% names(x@QhatModel.object@parameters@values), '1',
                                                                                          "0")))))

    state.shift = sapply(.Object@models, function(x) names(x@QhatModel.object@parameters@values[which(c(length(x@QhatModel.object@parameters@values$mean.a0),
                                                                                                        length(x@QhatModel.object@parameters@values$mean.a1),
                                                                                                        length(x@QhatModel.object@parameters@values$mean.a0.amp),
                                                                                                        length(x@QhatModel.object@parameters@values$mean.a0.phase),
                                                                                                        length(x@QhatModel.object@parameters@values$mean.a0.disp),
                                                                                                        length(x@QhatModel.object@parameters@values$mean.a1.amp),
                                                                                                        length(x@QhatModel.object@parameters@values$mean.a1.phase),
                                                                                                        length(x@QhatModel.object@parameters@values$mean.a1.disp)
                                                                                                        ) > 1)]))
    # simplify to mean.a0 or mean.a1 if seasonal..
    state.shift = unlist(lapply(1:length(state.shift), function(x) ifelse(length(unlist(state.shift[[x]])) < 1, "none",substr(unlist(state.shift[[x]]),1,7))))


    # build matrix of all
    all.models.matrix = data.frame(model.index = 1:length(.Object@models), model.names = names(.Object@models),nparams = nparams,nstates = nstates, state.structure = state.structure, data.trans = data.transform, error.dist = error.distribution, auto.corr = auto.correlation, state.shift = state.shift, ref.model = NA, row.names = 1)

    all.models.set = 0
    # seperate, order, then rebine if multiple state shifts being investigated
    if('mean.a0' %in% unique(all.models.matrix$state.shift)){

      all.models.matrix.a0 = all.models.matrix[which(all.models.matrix$state.shift == 'mean.a0'),]

      all.models.matrix.a0 = all.models.matrix.a0[order(all.models.matrix.a0$nparams,
                                                        all.models.matrix.a0$nstates),]
      all.models.set = all.models.set +1
    }

    if('mean.a1' %in% unique(all.models.matrix$state.shift)){

      all.models.matrix.a1 = all.models.matrix[which(all.models.matrix$state.shift == 'mean.a1'),]

      all.models.matrix.a1 = all.models.matrix.a1[order(all.models.matrix.a1$nparams,
                                                        all.models.matrix.a1$nstates),]
      all.models.set = all.models.set +1
    }
    if('none' %in% unique(all.models.matrix$state.shift)){

      all.models.matrix.none = all.models.matrix[which(all.models.matrix$state.shift == 'none'),]

      all.models.matrix.none = all.models.matrix.none[order(all.models.matrix.none$nparams,
                                                            all.models.matrix.none$nstates),]
      all.models.set = all.models.set +1
    }

    # remove duplicate one.state models when multiple state.shifts are investigated..
    if('mean.a0' %in% unique(all.models.matrix$state.shift) &&
       'mean.a1' %in% unique(all.models.matrix$state.shift)){

      if('none' %in% unique(all.models.matrix$state.shift)){

        if('log' %in% all.models.matrix$data.trans){

          if('normal' %in% all.models.matrix$error.dist){
            dup.models = all.models.matrix.none[which(all.models.matrix.none$state.shift == 'none' &
                                                        all.models.matrix.none$data.trans== 'log' &
                                                        all.models.matrix.none$error.dist== 'normal'),]
            remove.models =  dup.models[unique(dup.models$auto.corr)*2+1,]

            all.models.matrix.none = all.models.matrix.none[is.na(match(all.models.matrix.none$model.names,remove.models$model.names)),]
          }

          if('truc.normal' %in% all.models.matrix$error.dist){
            dup.models = all.models.matrix.none[which(all.models.matrix.none$state.shift == 'none' &
                                                        all.models.matrix.none$data.trans== 'log' &
                                                        all.models.matrix.none$error.dist== 'truc.normal'),]
            remove.models =  dup.models[unique(dup.models$auto.corr)*2+1,]

            all.models.matrix.none = all.models.matrix.none[is.na(match(all.models.matrix.none$model.names,remove.models$model.names)),]
          }

          if('gamma' %in% all.models.matrix$error.dist){
            dup.models = all.models.matrix.none[which(all.models.matrix.none$state.shift == 'none' &
                                                        all.models.matrix.none$data.trans== 'log' &
                                                        all.models.matrix.none$error.dist== 'gamma'),]
            remove.models =  dup.models[unique(dup.models$auto.corr)*2+1,]

            all.models.matrix.none = all.models.matrix.none[is.na(match(all.models.matrix.none$model.names,remove.models$model.names)),]
          }
        }

        if('boxcox' %in% all.models.matrix$data.trans){

          if('normal' %in% all.models.matrix$error.dist){
            dup.models = all.models.matrix.none[which(all.models.matrix.none$state.shift == 'none' &
                                                        all.models.matrix.none$data.trans== 'boxcox' &
                                                        all.models.matrix.none$error.dist== 'normal'),]
            remove.models =  dup.models[unique(dup.models$auto.corr)*2+1,]

            all.models.matrix.none = all.models.matrix.none[is.na(match(all.models.matrix.none$model.names,remove.models$model.names)),]
          }

          if('truc.normal' %in% all.models.matrix$error.dist){
            dup.models = all.models.matrix.none[which(all.models.matrix.none$state.shift == 'none' &
                                                        all.models.matrix.none$data.trans== 'boxcox' &
                                                        all.models.matrix.none$error.dist== 'truc.normal'),]
            remove.models =  dup.models[unique(dup.models$auto.corr)*2+1,]

            all.models.matrix.none = all.models.matrix.none[is.na(match(all.models.matrix.none$model.names,remove.models$model.names)),]
          }

          if('gamma' %in% all.models.matrix$error.dist){
            dup.models = all.models.matrix.none[which(all.models.matrix.none$state.shift == 'none' &
                                                        all.models.matrix.none$data.trans== 'boxcox' &
                                                        all.models.matrix.none$error.dist== 'gamma'),]
            remove.models =  dup.models[unique(dup.models$auto.corr)*2+1,]

            all.models.matrix.none = all.models.matrix.none[is.na(match(all.models.matrix.none$model.names,remove.models$model.names)),]
          }
        }
      }
    }

    # make list to iterate through model sets
    if(all(c('none','mean.a0','mean.a1') %in% unique(all.models.matrix$state.shift))){
      all.models.matrix = rbind(all.models.matrix.none,all.models.matrix.a0,all.models.matrix.a1)
    }else if(all(c('none','mean.a0') %in% unique(all.models.matrix$state.shift))){
      all.models.matrix = rbind(all.models.matrix.none,all.models.matrix.a0)
    }else if(all(c('none','mean.a1') %in% unique(all.models.matrix$state.shift))){
      all.models.matrix = rbind(all.models.matrix.none,all.models.matrix.a1)
    }else if(all(c('mean.a0','mean.a1') %in% unique(all.models.matrix$state.shift))){
      all.models.matrix = rbind(all.models.matrix.a0,all.models.matrix.a1)
    }else if('none' %in% unique(all.models.matrix$state.shift)){
      all.models.matrix = all.models.matrix.none
    }else if('mean.a0' %in% unique(all.models.matrix$state.shift)){
      all.models.matrix = all.models.matrix.a0
    }else if('mean.a1' %in% unique(all.models.matrix$state.shift)){
      all.models.matrix = all.models.matrix.a1
    }

    # isolate each model set to set reference models
    for(i in 1:length(unique(all.models.matrix$state.shift))){
      all.models.matrix.set= all.models.matrix[which(all.models.matrix$state.shift == unique(all.models.matrix$state.shift)[i]),]

      if('log' %in% all.models.matrix.set$data.trans){
        log.models = all.models.matrix.set[which(all.models.matrix.set$data.trans == 'log'),]

        if('normal' %in% log.models$error.dist){
          temp.models = log.models[which(log.models$error.dist == 'normal'),]

          if('none' %in% temp.models$state.shift){
            temp.models$ref.model[1] = ""
            temp.models$ref.model[2:NROW(temp.models)] = temp.models$model.names[1:(NROW(temp.models)-1)]
          }else if('mean.a0' %in% temp.models$state.shift){
            temp.models$ref.model[1] =  ifelse(i == 1, "",max.temp.normal.log.model)
            temp.models$ref.model[2:NROW(temp.models)] = temp.models$model.names[1:(NROW(temp.models)-1)]
          }else if('mean.a1' %in% temp.models$state.shift){
            temp.models$ref.model[1] =  ifelse(i == 1, "",max.temp.normal.log.model)
            temp.models$ref.model[2:NROW(temp.models)] = temp.models$model.names[1:(NROW(temp.models)-1)]
          }

          # reestablish in all. models.
          all.models.matrix$ref.model[!is.na(match(all.models.matrix$model.names,temp.models$model.names))] <- temp.models$ref.model

          # set max model for state of log models
          max.temp.normal.log.model = temp.models$model.names[NROW(temp.models)]
        }

        if('truc.normal' %in% log.models$error.dist){
          temp.models = log.models[which(log.models$error.dist == 'truc.normal'),]

          if('none' %in% temp.models$state.shift){
            temp.models$ref.model[1] = ""
            temp.models$ref.model[2:NROW(temp.models)] = temp.models$model.names[1:(NROW(temp.models)-1)]
          }else if('mean.a0' %in% temp.models$state.shift){
            temp.models$ref.model[1] =  ifelse(i == 1, "",max.temp.truc.normal.log.model)
            temp.models$ref.model[2:NROW(temp.models)] = temp.models$model.names[1:(NROW(temp.models)-1)]
          }else if('mean.a1' %in% temp.models$state.shift){
            temp.models$ref.model[1] =  ifelse(i == 1, "",max.temp.truc.normal.log.model)
            temp.models$ref.model[2:NROW(temp.models)] = temp.models$model.names[1:(NROW(temp.models)-1)]
          }

          # reestablish in all. models.
          all.models.matrix$ref.model[!is.na(match(all.models.matrix$model.names,temp.models$model.names))] <- temp.models$ref.model

          # set max model for state of log models
          max.temp.truc.normal.log.model = temp.models$model.names[NROW(temp.models)]
        }


        if('gamma' %in% log.models$error.dist){
          temp.models <- log.models[which(log.models$error.dist == 'gamma'),]

          if('none' %in% temp.models$state.shift){
            temp.models$ref.model[1] = ""
            temp.models$ref.model[2:NROW(temp.models)] = temp.models$model.names[1:(NROW(temp.models)-1)]
          }else if('mean.a0' %in% temp.models$state.shift){
            temp.models$ref.model[1] =  ifelse(i == 1, "",max.temp.gamma.log.model)
            temp.models$ref.model[2:NROW(temp.models)] = temp.models$model.names[1:(NROW(temp.models)-1)]
          }else if('mean.a1' %in% temp.models$state.shift){
            temp.models$ref.model[1] =  ifelse(i == 1, "",max.temp.gamma.log.model)
            temp.models$ref.model[2:NROW(temp.models)] = temp.models$model.names[1:(NROW(temp.models)-1)]
          }

          # reestablish in all. models.
          all.models.matrix$ref.model[!is.na(match(all.models.matrix$model.names,temp.models$model.names))] <- temp.models$ref.model

          # set max model for state of log models
          max.temp.gamma.log.model = temp.models$model.names[NROW(temp.models)]

        }
      }

      if('boxcox' %in% all.models.matrix.set$data.trans){
        boxcox.models = all.models.matrix.set[which(all.models.matrix.set$data.trans == 'boxcox'),]

        if('normal' %in% boxcox.models$error.dist){
          temp.models = boxcox.models[which(boxcox.models$error.dist == 'normal'),]

          if('none' %in% temp.models$state.shift){
            temp.models$ref.model[1] = ifelse('log' %in% all.models.matrix.set$data.trans,log.models$model.names[which(log.models$error.dist == 'normal')][1],"")
            temp.models$ref.model[2:NROW(temp.models)] = temp.models$model.names[1:(NROW(temp.models)-1)]
          }else if('mean.a0' %in% temp.models$state.shift){
            temp.models$ref.model[1] =  ifelse(i == 1, "", max.temp.normal.boxcox.model)
            temp.models$ref.model[2:NROW(temp.models)] = temp.models$model.names[1:(NROW(temp.models)-1)]
          }else if('mean.a1' %in% temp.models$state.shift){
            temp.models$ref.model[1] =  ifelse(i == 1, "", max.temp.normal.boxcox.model)
            temp.models$ref.model[2:NROW(temp.models)] = temp.models$model.names[1:(NROW(temp.models)-1)]
          }

          # reestablish in all. models.
          all.models.matrix$ref.model[!is.na(match(all.models.matrix$model.names,temp.models$model.names))] <- temp.models$ref.model

          # set max model for state of boxcox models
          max.temp.normal.boxcox.model = temp.models$model.names[NROW(temp.models)]
        }

        if('truc.normal' %in% boxcox.models$error.dist){
          temp.models = boxcox.models[which(boxcox.models$error.dist == 'truc.normal'),]

          if('none' %in% temp.models$state.shift){
            temp.models$ref.model[1] = ifelse('log' %in% all.models.matrix.set$data.trans,log.models$model.names[which(log.models$error.dist == 'truc.normal')][1],"")
            temp.models$ref.model[2:NROW(temp.models)] = temp.models$model.names[1:(NROW(temp.models)-1)]
          }else if('mean.a0' %in% temp.models$state.shift){
            temp.models$ref.model[1] =  ifelse(i == 1, "", max.temp.truc.normal.boxcox.model)
            temp.models$ref.model[2:NROW(temp.models)] = temp.models$model.names[1:(NROW(temp.models)-1)]
          }else if('mean.a1' %in% temp.models$state.shift){
            temp.models$ref.model[1] =  ifelse(i == 1, "", max.temp.truc.normal.boxcox.model)
            temp.models$ref.model[2:NROW(temp.models)] = temp.models$model.names[1:(NROW(temp.models)-1)]
          }

          # reestablish in all. models.
          all.models.matrix$ref.model[!is.na(match(all.models.matrix$model.names,temp.models$model.names))] <- temp.models$ref.model

          # set max model for state of boxcox models
          max.temp.truc.normal.boxcox.model = temp.models$model.names[NROW(temp.models)]
        }


        if('gamma' %in% boxcox.models$error.dist){
          temp.models = boxcox.models[which(boxcox.models$error.dist == 'gamma'),]

          if('none' %in% temp.models$state.shift){
            temp.models$ref.model[1] = ifelse('log' %in% all.models.matrix.set$data.trans,log.models$model.names[which(log.models$error.dist == 'gamma')][1],"")
            temp.models$ref.model[2:NROW(temp.models)] = temp.models$model.names[1:(NROW(temp.models)-1)]
          }else if('mean.a0' %in% temp.models$state.shift){
            temp.models$ref.model[1] =  ifelse(i == 1, "",max.temp.gamma.boxcox.model)
            temp.models$ref.model[2:NROW(temp.models)] = temp.models$model.names[1:(NROW(temp.models)-1)]
          }else if('mean.a1' %in% temp.models$state.shift){
            temp.models$ref.model[1] =  ifelse(i == 1, "",max.temp.gamma.boxcox.model)
            temp.models$ref.model[2:NROW(temp.models)] = temp.models$model.names[1:(NROW(temp.models)-1)]
          }

          # reestablish in all. models.
          all.models.matrix$ref.model[!is.na(match(all.models.matrix$model.names,temp.models$model.names))] <- temp.models$ref.model

          # set max model for state of boxcox models
          max.temp.gamma.boxcox.model = temp.models$model.names[NROW(temp.models)]

        }
      }

    }


    #export summary
    .Object@models.summary <- all.models.matrix
  }

  # reorder to ensure empty calib ref are first...
  empty.ref = which(.Object@models.summary$ref.model =="")
  full.ref = which(.Object@models.summary$ref.model !="")
  .Object@models.summary = rbind(.Object@models.summary[empty.ref,],.Object@models.summary[full.ref,])

  # clean all.models (remove the duplicate one-state if so) this happens if investigating slope and intercept seperately
  .Object@models <- .Object@models[match(.Object@models.summary$model.names,names(.Object@models))]


  ### Then adjust calib.reference details based on models.summary from user or machine

  # Define Refernce models. That is, the calibration obecjive function that this models needs to meet or exceed.
  .Object@calib.reference.model.name <- as.list(.Object@models.summary$ref.model)
  names(.Object@calib.reference.model.name) <- .Object@models.summary$model.names


  model.names = names(.Object@models.summary$model.names)
  .Object@calib.reference.criteria.met = rep(F,length(model.names))
  names(.Object@calib.reference.criteria.met) <- model.names

  return(.Object)

}
)

# @exportMethod fit
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

        assign("model", model, envir=new.env())
                 #myEnv  (parent = baseenv()))

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
      if (!(current.model.name %in% .Object@calib.reference.criteria.met)){
        filt[i] <- T
      }else if(.Object@calib.reference.criteria.met[[current.model.name]]){
        filt[i] <- T
      }
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

