##' @include abstracts.R QhatModel.homo.gamma.linear.R
## @export
QhatModel.subAnnual.homo.gamma.linear <- setClass(
  # Set the name for the class
  "QhatModel.subAnnual.homo.gamma.linear",

  package='hydroState',

  contains=c('QhatModel.homo.gamma.linear'),

  slots = c(
    subAnnual.Monthly.Steps = 'numeric',
    subAnnual.dependent.mean.a0 = 'logical',
    subAnnual.dependent.mean.a1 = 'logical',
    subAnnual.dependent.std.a0 = 'logical'
  ),

  # Set the default values for the slots. (optional)
  prototype=list(
    input.data = data.frame(year=c(0),month=c(0),precipitation=c(0)),
    nStates = Inf,
    use.truncated.dist=F,
    subAnnual.Monthly.Steps = c(2, 5, 8, 11),
    subAnnual.dependent.mean.a0 = T,
    subAnnual.dependent.mean.a1 = F,
    subAnnual.dependent.std.a0 = F,
    parameters = new('parameters',c('mean.a0.amp', 'mean.a0.phase', 'mean.a0.disp',
                                    'mean.a1','std.a0'),c(1,1,1,1,1))
  )
)

# Valid object?
validObject <- function(object) {
  TRUE

}
setValidity("QhatModel.subAnnual.homo.gamma.linear", validObject)

# Initialise object
setMethod("initialize","QhatModel.subAnnual.homo.gamma.linear", function(.Object, input.data, transition.graph=matrix(T,2,2),
                                                                      state.dependent.mean.a0=T,state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                                      subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F,subAnnual.dependent.std.a0=F) {
  .Object@input.data <- input.data
  .Object@use.truncated.dist <- F
  .Object@nStates = ncol(transition.graph)

  # Check and set definition of seasons.
  .Object <- setSeasons(.Object, input.data)

  # Setup parameter definitions.
  .Object <- setupParameters(.Object, state.dependent.mean.a0, state.dependent.mean.a1, state.dependent.std.a0,
    state.dependent.mean.AR1=NA, state.dependent.mean.AR2=NA, state.dependent.mean.AR3=NA,
    subAnnual.dependent.mean.a0, subAnnual.dependent.mean.a1, subAnnual.dependent.std.a0)

  validObject(.Object)
  .Object
}
)
#' @exportMethod setSeasons
setGeneric(name="setSeasons",def=function(.Object, input.data) {standardGeneric("setSeasons")})
setMethod(f="setSeasons",signature=c("QhatModel.subAnnual.homo.gamma.linear",'data.frame'),
          definition=function(.Object, input.data)
{

  # Set-up to run for all periods with continuous observations of independent variable (precipitation)
  delta = getStartEndIndex(input.data)

  subAnnual.Monthly.Steps = list()

  for(j in 1:NROW(delta)){

    data = input.data[delta[j,1]:delta[j,2],]

    if(delta[j,2] - delta[j,1] > 12){
      # Get the seasons from the input data.
      subAnnual.Monthly.Steps[[j]] = sort(unique(data$month))
      subAnnual.Monthly.StepSize.min = rep(NA, length(subAnnual.Monthly.Steps[[j]]))
      subAnnual.Monthly.StepSize.max = rep(NA, length(subAnnual.Monthly.Steps[[j]]))


      # Loop though each unique month and calc. the time step size.
      for (i in 1:length(subAnnual.Monthly.Steps[[j]])) {
        step.sizes = diff(which(data$month==subAnnual.Monthly.Steps[[j]][i]))
        subAnnual.Monthly.StepSize.min[i] = min(step.sizes)
        subAnnual.Monthly.StepSize.max[i] = max(step.sizes)
      }

      # Check that there is only one time-step size for each month.
      for (i in 1:length(subAnnual.Monthly.Steps[[j]])) {
        if (subAnnual.Monthly.StepSize.min[i] != subAnnual.Monthly.StepSize.max[i])
          stop(paste('The time steps for month', subAnnual.Monthly.Steps[[j]][i],'are not constant. Check input data.'))
      }

      # Check that each month has the same time step size.
      if (length(unique(subAnnual.Monthly.StepSize.min))!=1 || length(unique(subAnnual.Monthly.StepSize.max))!=1)
        warning(paste('The time step sizes differ across the year. Consider checking the input data.'))

   }else{
     subAnnual.Monthly.Steps[[j]] = sort(unique(data$month))
     step.sizes = diff(subAnnual.Monthly.Steps[[j]])
     subAnnual.Monthly.StepSize.min =  min(step.sizes)
     subAnnual.Monthly.StepSize.max =  max(step.sizes)

     if(!(subAnnual.Monthly.StepSize.min %in% c(1,-11)) || !(subAnnual.Monthly.StepSize.max %in% c(1,-11))){
       warning(paste('The monthly time step size differ across the year. Consider checking the input data.'))
     }
   }


  }

  #Assign the definitions of seasons to the object
  .Object@subAnnual.Monthly.Steps = unlist(subAnnual.Monthly.Steps)

  return(.Object)
}
)

setGeneric(name="setupParameters",def=function(.Object, state.dependent.mean.a0, state.dependent.mean.a1, state.dependent.std.a0,
                                               state.dependent.mean.AR1, state.dependent.mean.AR2, state.dependent.mean.AR3,
                                               subAnnual.dependent.mean.a0, subAnnual.dependent.mean.a1, subAnnual.dependent.std.a0) {standardGeneric("setupParameters")})
setMethod(f="setupParameters",signature=c("QhatModel.subAnnual.homo.gamma.linear",rep('logical',9)),
definition=function(.Object, state.dependent.mean.a0, state.dependent.mean.a1, state.dependent.std.a0,
                    state.dependent.mean.AR1, state.dependent.mean.AR2, state.dependent.mean.AR3,
                    subAnnual.dependent.mean.a0, subAnnual.dependent.mean.a1, subAnnual.dependent.std.a0)
{

            # Set seasonal slots
            .Object@subAnnual.dependent.mean.a0 <- subAnnual.dependent.mean.a0
            .Object@subAnnual.dependent.mean.a1 <- subAnnual.dependent.mean.a1
            .Object@subAnnual.dependent.std.a0 <- subAnnual.dependent.std.a0

            # Setup AR terms.
            AR.param.lengths = numeric()
            AR.param.names = c();
            if (!is.na(state.dependent.mean.AR1)) {
              if (state.dependent.mean.AR1) {
                AR.param.lengths = c(AR.param.lengths, 1)
              } else {
                AR.param.lengths = c(AR.param.lengths, 0)
              }
              AR.param.names = c(AR.param.names, 'mean.AR1')
            }
            if (!is.na(state.dependent.mean.AR2)) {
              if (state.dependent.mean.AR2) {
                AR.param.lengths = c(AR.param.lengths, 1)
              } else {
                AR.param.lengths = c(AR.param.lengths, 0)
              }
              AR.param.names = c(AR.param.names, 'mean.AR2')
            }
            if (!is.na(state.dependent.mean.AR3)) {
              if (state.dependent.mean.AR3) {
                AR.param.lengths = c(AR.param.lengths, 1)
              } else {
                AR.param.lengths = c(AR.param.lengths, 0)
              }
              AR.param.names = c(AR.param.names, 'mean.AR3')
            }

            # Set the number of parameter values per parameter name. At this point the seasons are not accounted for.
            parameter.length <- as.numeric(c(rep(state.dependent.mean.a0, max(c(1, subAnnual.dependent.mean.a0*3))),
                                             rep(state.dependent.mean.a1, max(c(1, subAnnual.dependent.mean.a1*3))),
                                             AR.param.lengths,
                                             rep(state.dependent.std.a0, max(c(1, subAnnual.dependent.std.a0*3))))) * (.Object@nStates-1) + 1

            # Set up model terms for mean and standard deviation.
            if (subAnnual.dependent.std.a0) {
              if (subAnnual.dependent.mean.a0 && !subAnnual.dependent.mean.a1) {
                .Object@parameters = new('parameters', c('mean.a0.amp', 'mean.a0.phase', 'mean.a0.disp',
                                                         'mean.a1', AR.param.names,
                                                         'std.a0.amp', 'std.a0.phase', 'std.a0.disp'), parameter.length)
              } else if (!subAnnual.dependent.mean.a0 && subAnnual.dependent.mean.a1) {
                .Object@parameters = new('parameters', c('mean.a0', 'mean.a1.amp', 'mean.a1.phase', 'mean.a1.disp',
                                                         AR.param.names,
                                                         'std.a0.amp', 'std.a0.phase', 'std.a0.disp'), parameter.length)
              } else if (subAnnual.dependent.mean.a0 && subAnnual.dependent.mean.a1) {
                .Object@parameters = new('parameters', c('mean.a0.amp', 'mean.a0.phase', 'mean.a0.disp',
                                                         'mean.a1.amp', 'mean.a1.phase', 'mean.a1.disp',
                                                         AR.param.names,
                                                         'std.a0.amp', 'std.a0.phase', 'std.a0.disp'), parameter.length)
              } else {
                .Object@parameters = new('parameters', c('mean.a0', 'mean.a1', AR.param.names,
                                                         'std.a0.amp', 'std.a0.phase', 'std.a0.disp'), parameter.length)
              }

            } else {
              if (subAnnual.dependent.mean.a0 && !subAnnual.dependent.mean.a1) {
                .Object@parameters = new('parameters', c('mean.a0.amp', 'mean.a0.phase', 'mean.a0.disp',
                                                         'mean.a1', AR.param.names, 'std.a0'), parameter.length)
              } else if (!subAnnual.dependent.mean.a0 && subAnnual.dependent.mean.a1) {
                .Object@parameters = new('parameters', c('mean.a0',
                                                         'mean.a1.amp', 'mean.a1.phase', 'mean.a1.disp',
                                                         AR.param.names, 'std.a0'), parameter.length)
              } else if (subAnnual.dependent.mean.a0 && subAnnual.dependent.mean.a1) {
                .Object@parameters = new('parameters', c('mean.a0.amp', 'mean.a0.phase', 'mean.a0.disp',
                                                         'mean.a1.amp', 'mean.a1.phase', 'mean.a1.disp',
                                                         AR.param.names, 'std.a0'), parameter.length)
              } else {
                .Object@parameters = new('parameters', c('mean.a0', 'mean.a1',AR.param.names,'std.a0'), parameter.length)
              }
            }
            return(.Object)
}
)

setMethod(f="getMean",signature=c("QhatModel.subAnnual.homo.gamma.linear","data.frame"),definition=function(.Object, data)
{
  # Get object parameter list
  parameters = getParameters(.Object@parameters)

  # Get number of data rows
  nrows = length(data$Qhat.precipitation);

  # Setup mean.a0 parameter.
  if (.Object@subAnnual.dependent.mean.a0) {
    ncols.mean.a0.amp = length(parameters$mean.a0.amp)
    ncols.mean.a0.phase = length(parameters$mean.a0.phase)
    ncols.mean.a0.disp = length(parameters$mean.a0.disp)

    if (ncols.mean.a0.amp==1) {
      mean.a0.amp = rep(parameters$mean.a0.amp, length.out=.Object@nStates);
    } else {
      mean.a0.amp = parameters$mean.a0.amp
    }

    if (ncols.mean.a0.phase==1){
      mean.a0.phase = rep(parameters$mean.a0.phase, length.out=.Object@nStates);
    } else {
      mean.a0.phase = parameters$mean.a0.phase
    }

    if (ncols.mean.a0.disp==1){
      mean.a0.disp = rep(parameters$mean.a0.disp, length.out=.Object@nStates);
    } else {
      mean.a0.disp = parameters$mean.a0.disp
    }

    # Apply sinusoidal model.
    mean.a0 = matrix(rep(mean.a0.disp,each=nrows),nrows,.Object@nStates) + matrix(rep(mean.a0.amp,each=nrows),nrows,.Object@nStates) *
      sin( 2*pi*(matrix(data$month/12,nrows,.Object@nStates) + matrix(rep(mean.a0.phase,each=nrows),nrows,.Object@nStates)))

  } else {
    mean.a0 = matrix(rep(parameters$mean.a0,each=nrows),nrows,.Object@nStates);
  }

  # Setup mean.a1 parameter.
  if (.Object@subAnnual.dependent.mean.a1) {
    ncols.mean.a1.amp = length(parameters$mean.a1.amp)
    ncols.mean.a1.phase = length(parameters$mean.a1.phase)
    ncols.mean.a1.disp = length(parameters$mean.a1.disp)

    if (ncols.mean.a1.amp==1) {
      mean.a1.amp = rep(parameters$mean.a1.amp, length.out=.Object@nStates);
    } else {
      mean.a1.amp = parameters$mean.a1.amp
    }

    if (ncols.mean.a1.phase==1){
      mean.a1.phase = rep(parameters$mean.a1.phase, length.out=.Object@nStates);
    } else {
      mean.a1.phase = parameters$mean.a1.phase
    }

    if (ncols.mean.a1.disp==1){
      mean.a1.disp = rep(parameters$mean.a1.disp, length.out=.Object@nStates);
    } else {
      mean.a1.disp = parameters$mean.a1.disp
    }

    # Apply sinusoidal model.
    mean.a1 = matrix(rep(mean.a1.disp,each=nrows),nrows,.Object@nStates) + matrix(rep(mean.a1.amp,each=nrows),nrows,.Object@nStates) *
      sin( 2*pi*(matrix(data$month/12,nrows,.Object@nStates) + matrix(rep(mean.a1.phase,each=nrows),nrows,.Object@nStates)))

  } else {
    mean.a1 = matrix(rep(parameters$mean.a1,each=nrows),nrows,.Object@nStates);
  }

  # Expand precip to the number of states.
  precip.data = matrix(data$Qhat.precipitation,nrows,.Object@nStates);

  # Calculate the sesonal estimate of mean flow.
  Qhat.model <- precip.data * mean.a1 + 100 * mean.a0

  return(Qhat.model)
}
)

 setMethod(f="getVariance",signature=c("QhatModel.subAnnual.homo.gamma.linear","data.frame"),definition=function(.Object, data)
 {
    # Get object parameter list
   parameters = getParameters(.Object@parameters)

   # Get number of data rows
   nrows = length(data$Qhat.precipitation);

   # Setup mean.a0 parameter.
   if (.Object@subAnnual.dependent.std.a0) {
     ncols.std.a0.amp = length(parameters$std.a0.amp)
     ncols.std.a0.phase = length(parameters$std.a0.phase)
     ncols.std.a0.disp = length(parameters$std.a0.disp)

     if (ncols.std.a0.amp==1) {
       std.a0.amp = rep(parameters$std.a0.amp, length.out=.Object@nStates);
     } else {
       std.a0.amp = parameters$std.a0.amp
     }

     if (ncols.std.a0.phase==1){
       std.a0.phase = rep(parameters$std.a0.phase, length.out=.Object@nStates);
     } else {
       std.a0.phase = parameters$std.a0.phase
     }

     if (ncols.std.a0.disp==1){
       std.a0.disp = rep(parameters$std.a0.disp, length.out=.Object@nStates);
     } else {
       std.a0.disp = parameters$std.a0.disp
     }

     # Apply sinusoidal model.
     std.a0 = matrix(rep(std.a0.disp,each=nrows),nrows,.Object@nStates) + matrix(rep(std.a0.amp,each=nrows),nrows,.Object@nStates) *
       sin( 2*pi*(matrix(data$month/12,nrows,.Object@nStates) + matrix(rep(std.a0.phase,each=nrows),nrows,.Object@nStates)))

   } else {
     std.a0 = matrix(rep(parameters$std.a0,each=nrows),nrows,.Object@nStates);
   }

   # Scale by variance of transformed observed flow and return.
   return(std.a0 * var(data$Qhat.flow, na.rm=T))

 }
 )
