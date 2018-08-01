##' @include abstracts.R QhatModel.homo.gamma.linear.R
##' @export
QhatModel.seasonal.homo.gamma.linear <- setClass(
  # Set the name for the class
  "QhatModel.seasonal.homo.gamma.linear",

  package='hydroState',

  contains=c('QhatModel.homo.gamma.linear'),

  slots = c(
    summer.endMonth = 'numeric',
    autumn.endMonth = 'numeric',
    winter.endMonth = 'numeric',
    spring.endMonth = 'numeric',
    seasonal.dependent.mean.a0 = 'logical',
    seasonal.dependent.mean.a1 = 'logical',
    seasonal.dependent.std.a0 = 'logical'
  ),

  # Set the default values for the slots. (optional)
  prototype=list(
    input.data = data.frame(year=c(0),month=c(0),precipitation=c(0)),
    nStates = Inf,
    use.truncated.dist=F,
    summer.endMonth=2,
    autumn.endMonth=5,
    winter.endMonth=8,
    spring.endMonth=11,
    seasonal.dependent.mean.a0 = T,
    seasonal.dependent.mean.a1 = F,
    seasonal.dependent.std.a0 = F,
    parameters = new('parameters',c('mean.summer.a0', 'mean.autumn.a0', 'mean.winter.a0', 'mean.spring.a0',
                                    'mean.a1',
                                    'std.a0'),c(1,1,1,1,1,1))
  )
)

# Valid object?
validObject <- function(object) {
  TRUE

}
setValidity("QhatModel.seasonal.homo.gamma.linear", validObject)

# Initialise object
setMethod("initialize","QhatModel.seasonal.homo.gamma.linear", function(.Object, input.data, transition.graph=matrix(T,2,2),
                                                                      state.dependent.mean.a0=T,state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                                      seasonal.dependent.mean.a0=T, seasonal.dependent.mean.a1=F,seasonal.dependent.std.a0=F,
                                                                      summer.endMonth=2,autumn.endMonth=5, winter.endMonth=8, spring.endMonth=11) {
  .Object@use.truncated.dist <- F
  .Object@nStates = ncol(transition.graph)

  # Check and set definition of seasons.
  .Object <- setSeasons(.Object, input.data, summer.endMonth,autumn.endMonth, winter.endMonth, spring.endMonth)

  # Setup parameter definitions.
  .Object <- setupParameters(.Object, state.dependent.mean.a0, state.dependent.mean.a1, state.dependent.std.a0,
    state.dependent.mean.AR1=NA, state.dependent.mean.AR2=NA, state.dependent.mean.AR3=NA,
    seasonal.dependent.mean.a0, seasonal.dependent.mean.a1, seasonal.dependent.std.a0)

  validObject(.Object)
  .Object
}
)

setGeneric(name="setSeasons",def=function(.Object, input.data, summer.endMonth,autumn.endMonth, winter.endMonth, spring.endMonth) {standardGeneric("setSeasons")})
setMethod(f="setSeasons",signature=c("QhatModel.seasonal.homo.gamma.linear",'data.frame',rep("numeric",4)),
          definition=function(.Object, input.data, summer.endMonth,autumn.endMonth, winter.endMonth, spring.endMonth)
{
  # Check the definitions of the seasons is unique.
  all.seasons = c(summer.endMonth,autumn.endMonth, winter.endMonth, spring.endMonth)
  if ( sum(summer.endMonth==all.seasons)>1 ||
       sum(autumn.endMonth==all.seasons)>1 ||
       sum(winter.endMonth==all.seasons)>1 ||
       sum(spring.endMonth==all.seasons)>1)
    stop('The end months for each season are not unique.')

  #Assign the definitions of seasons to the object
  .Object@summer.endMonth = summer.endMonth
  .Object@autumn.endMonth = autumn.endMonth
  .Object@winter.endMonth = winter.endMonth
  .Object@spring.endMonth = spring.endMonth

  # Check that all input data can be mapped to a season.
  # Build a matrix of ao and a1 where each row is a different season.
  data.season = rep(F, nrow(input.data))
  data.season[input.data$month == .Object@summer.endMonth] = T
  data.season[input.data$month == .Object@autumn.endMonth] = T
  data.season[input.data$month == .Object@winter.endMonth] = T
  data.season[input.data$month == .Object@spring.endMonth] = T
  if (any(!data.season))
    stop('A season cannot be assigned to some of the input data. Check the input data is at a seasonal step that corresponds to the defined end months of each season.')

  return(.Object)
}
)

setGeneric(name="setupParameters",def=function(.Object, state.dependent.mean.a0, state.dependent.mean.a1, state.dependent.std.a0,
                                               state.dependent.mean.AR1, state.dependent.mean.AR2, state.dependent.mean.AR3,
                                               seasonal.dependent.mean.a0, seasonal.dependent.mean.a1, seasonal.dependent.std.a0) {standardGeneric("setupParameters")})
setMethod(f="setupParameters",signature=c("QhatModel.seasonal.homo.gamma.linear",rep('logical',9)),
definition=function(.Object, state.dependent.mean.a0, state.dependent.mean.a1, state.dependent.std.a0,
                    state.dependent.mean.AR1, state.dependent.mean.AR2, state.dependent.mean.AR3,
                    seasonal.dependent.mean.a0, seasonal.dependent.mean.a1, seasonal.dependent.std.a0)
{

            # Set seasonal slots
            .Object@seasonal.dependent.mean.a0 <- seasonal.dependent.mean.a0
            .Object@seasonal.dependent.mean.a1 <- seasonal.dependent.mean.a1
            .Object@seasonal.dependent.std.a0 <- seasonal.dependent.std.a0

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
            parameter.length <- as.numeric(c(rep(state.dependent.mean.a0, max(c(1, seasonal.dependent.mean.a0*4))),
                                             rep(state.dependent.mean.a1, max(c(1, seasonal.dependent.mean.a1*4))),
                                             AR.param.lengths,
                                             rep(state.dependent.std.a0, max(c(1, seasonal.dependent.std.a0*4))))) * (.Object@nStates-1) + 1

            # Set up model terms for mean and standard deviation.
            if (seasonal.dependent.std.a0) {
              if (seasonal.dependent.mean.a0 && !seasonal.dependent.mean.a1) {
                .Object@parameters = new('parameters', c('mean.summer.a0', 'mean.autumn.a0', 'mean.winter.a0', 'mean.spring.a0',
                                                         'mean.a1', AR.param.names,
                                                         'std.summer.a0', 'std.autumn.a0', 'std.winter.a0', 'std.spring.a0'), parameter.length)
              } else if (!seasonal.dependent.mean.a0 && seasonal.dependent.mean.a1) {
                .Object@parameters = new('parameters', c('mean.a0', 'mean.summer.a1', 'mean.autumn.a1', 'mean.winter.a1', 'mean.spring.a1',
                                                         AR.param.names,
                                                         'std.summer.a0', 'std.autumn.a0', 'std.winter.a0', 'std.spring.a0'), parameter.length)
              } else if (seasonal.dependent.mean.a0 && seasonal.dependent.mean.a1) {
                .Object@parameters = new('parameters', c('mean.summer.a0', 'mean.autumn.a0', 'mean.winter.a0', 'mean.spring.a0',
                                                         'mean.summer.a1', 'mean.autumn.a1', 'mean.winter.a1', 'mean.spring.a1',
                                                         AR.param.names,
                                                         'std.summer.a0', 'std.autumn.a0', 'std.winter.a0', 'std.spring.a0'), parameter.length)
              } else {
                .Object@parameters = new('parameters', c('mean.a0', 'mean.a1', AR.param.names,
                                                         'std.summer.a0', 'std.autumn.a0', 'std.winter.a0', 'std.spring.a0'), parameter.length)
              }

            } else {
              if (seasonal.dependent.mean.a0 && !seasonal.dependent.mean.a1) {
                .Object@parameters = new('parameters', c('mean.summer.a0', 'mean.autumn.a0', 'mean.winter.a0', 'mean.spring.a0',
                                                         'mean.a1', AR.param.names, 'std.a0'), parameter.length)
              } else if (!seasonal.dependent.mean.a0 && seasonal.dependent.mean.a1) {
                .Object@parameters = new('parameters', c('mean.a0', 'mean.summer.a1', 'mean.autumn.a1', 'mean.winter.a1', 'mean.spring.a1',
                                                         AR.param.names, 'std.a0'), parameter.length)
              } else if (seasonal.dependent.mean.a0 && seasonal.dependent.mean.a1) {
                .Object@parameters = new('parameters', c('mean.summer.a0', 'mean.autumn.a0', 'mean.winter.a0', 'mean.spring.a0',
                                                         'mean.summer.a1', 'mean.autumn.a1', 'mean.winter.a1', 'mean.spring.a1',
                                                         AR.param.names, 'std.a0'), parameter.length)
              } else {
                .Object@parameters = new('parameters', c('mean.a0', 'mean.a1',AR.param.names,'std.a0'), parameter.length)
              }
            }
            return(.Object)
}
)

setMethod(f="getMean",signature=c("QhatModel.seasonal.homo.gamma.linear","data.frame"),definition=function(.Object, data)
{
  # Get object parameter list
  parameters = getParameters(.Object@parameters)

  # Get number of data rows
  nrows = length(data$Qhat.precipitation);

  # Setup mean.a0 parameter.
  if (.Object@seasonal.dependent.mean.a0) {
    ncols.summer.a0 = length(parameters$mean.summer.a0)
    ncols.autumn.a0 = length(parameters$mean.autumn.a0)
    ncols.winter.a0 = length(parameters$mean.winter.a0)
    ncols.spring.a0 = length(parameters$mean.spring.a0)

    if (ncols.summer.a0==1) {
      mean.summer.a0 = rep(parameters$mean.summer.a0, length.out=.Object@nStates);
    } else {
      mean.summer.a0 = parameters$mean.summer.a0
    }

    if (ncols.autumn.a0==1){
      mean.autumn.a0 = rep(parameters$mean.autumn.a0, length.out=.Object@nStates);
    } else {
      mean.autumn.a0 = parameters$mean.autumn.a0
    }

    if (ncols.winter.a0==1){
      mean.winter.a0 = rep(parameters$mean.winter.a0, length.out=.Object@nStates);
    } else {
      mean.winter.a0 = parameters$mean.winter.a0
    }

    if (ncols.spring.a0==1){
      mean.spring.a0 = rep(parameters$mean.spring.a0, length.out=.Object@nStates);
    } else {
      mean.spring.a0 = parameters$mean.spring.a0
    }

    # Arrange seaonal parameters to a matrix
    mean.a0 = matrix(NA, nrows, .Object@nStates)
    filt = data$month == .Object@summer.endMonth
    mean.a0[filt, ]  = rep(mean.summer.a0, each=sum(filt))
    filt = data$month == .Object@autumn.endMonth
    mean.a0[filt, ]  = rep(mean.autumn.a0, each=sum(filt))
    filt = data$month == .Object@winter.endMonth
    mean.a0[filt, ]  = rep(mean.winter.a0, each=sum(filt))
    filt = data$month == .Object@spring.endMonth
    mean.a0[filt, ]  = rep(mean.spring.a0, each=sum(filt))

  } else {
    mean.a0 = matrix(rep(parameters$mean.a0,each=nrows),nrows,.Object@nStates);
  }

  # Setup mean.a1 parameter.
  if (.Object@seasonal.dependent.mean.a1) {
    ncols.summer.a1 = length(parameters$mean.summer.a1)
    ncols.autumn.a1 = length(parameters$mean.autumn.a1)
    ncols.winter.a1 = length(parameters$mean.winter.a1)
    ncols.spring.a1 = length(parameters$mean.spring.a1)

    if (ncols.summer.a1==1) {
      mean.summer.a1 = rep(parameters$mean.summer.a1, length.out=.Object@nStates);
    } else {
      mean.summer.a1 = parameters$mean.summer.a1
    }

    if (ncols.autumn.a1==1){
      mean.autumn.a1 = rep(parameters$mean.autumn.a1, length.out=.Object@nStates);
    } else {
      mean.autumn.a1 = parameters$mean.autumn.a1
    }

    if (ncols.winter.a1==1){
      mean.winter.a1 = rep(parameters$mean.winter.a1, length.out=.Object@nStates);
    } else {
      mean.winter.a1 = parameters$mean.winter.a1
    }

    if (ncols.spring.a1==1){
      mean.spring.a1 = rep(parameters$mean.spring.a1, length.out=.Object@nStates);
    } else {
      mean.spring.a1 = parameters$mean.spring.a1
    }

    # Arrange seaonal parameters to a matrix
    mean.a1 = matrix(NA, nrows, .Object@nStates)
    filt = data$month == .Object@summer.endMonth
    mean.a1[filt, ]  = rep(mean.summer.a1, each=sum(filt))
    filt = data$month == .Object@autumn.endMonth
    mean.a1[filt, ]  = rep(mean.autumn.a1, each=sum(filt))
    filt = data$month == .Object@winter.endMonth
    mean.a1[filt, ]  = rep(mean.winter.a1, each=sum(filt))
    filt = data$month == .Object@spring.endMonth
    mean.a1[filt, ]  = rep(mean.spring.a1, each=sum(filt))

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

 setMethod(f="getVariance",signature=c("QhatModel.seasonal.homo.gamma.linear","data.frame"),definition=function(.Object, data)
 {
    # Get object parameter list
   parameters = getParameters(.Object@parameters)

   # Get number of data rows
   nrows = length(data$Qhat.precipitation);

   # Setup std.a0 parameter.
   if (.Object@seasonal.dependent.std.a0) {
     ncols.summer.a0 = length(parameters$std.summer.a0)
     ncols.autumn.a0 = length(parameters$std.autumn.a0)
     ncols.winter.a0 = length(parameters$std.winter.a0)
     ncols.spring.a0 = length(parameters$std.spring.a0)

     if (ncols.summer.a0==1) {
       std.summer.a0 = rep(parameters$std.summer.a0, length.out=.Object@nStates);
     } else {
       std.summer.a0 = parameters$std.summer.a0
     }

     if (ncols.autumn.a0==1){
       std.autumn.a0 = rep(parameters$std.autumn.a0, length.out=.Object@nStates);
     } else {
       std.autumn.a0 = parameters$std.autumn.a0
     }

     if (ncols.winter.a0==1){
       std.winter.a0 = rep(parameters$std.winter.a0, length.out=.Object@nStates);
     } else {
       std.winter.a0 = parameters$std.winter.a0
     }

     if (ncols.spring.a0==1){
       std.spring.a0 = rep(parameters$std.spring.a0, length.out=.Object@nStates);
     } else {
       std.spring.a0 = parameters$std.spring.a0
     }

     # Arrange seaonal parameters to a matrix
     std.a0 = matrix(NA, nrows, .Object@nStates)
     filt = data$month == .Object@summer.endMonth
     std.a0[filt, ]  = rep(std.summer.a0, each=sum(filt))
     filt = data$month == .Object@autumn.endMonth
     std.a0[filt, ]  = rep(std.autumn.a0, each=sum(filt))
     filt = data$month == .Object@winter.endMonth
     std.a0[filt, ]  = rep(std.winter.a0, each=sum(filt))
     filt = data$month == .Object@spring.endMonth
     std.a0[filt, ]  = rep(std.spring.a0, each=sum(filt))

   } else {
     std.a0 = matrix(rep(parameters$std.a0,each=nrows),nrows,.Object@nStates);
   }

   #
   #
   # ncols.a0 = length(parameters$std.a0)
   # nrows = length(data$Qhat.precipitation);
   #
   #  Calculate the variance by multiplying the state scaling parameter by the variance in transformed flow FOR EACH SEASON.
   # seaonal.var = matrix(NA, nrow(data), .Object@nStates)
   # filt = data$month == .Object@summer.endMonth
   # seaonal.var[filt,] = var(data$Qhat.flow[filt], na.rm=T)
   # filt = data$month == .Object@autumn.endMonth
   # seaonal.var[filt,] = var(data$Qhat.flow[filt], na.rm=T)
   # filt = data$month == .Object@winter.endMonth
   # seaonal.var[filt,] = var(data$Qhat.flow[filt], na.rm=T)
   # filt = data$month == .Object@spring.endMonth
   # seaonal.var[filt,] = var(data$Qhat.flow[filt], na.rm=T)
   # std.a0 = matrix(rep(parameters$std.a0,each=nrows),nrows,.Object@nStates);

   # Scale by variance of transformed observed flow and return.
   return(std.a0 * var(data$Qhat.flow, na.rm=T))

 }
 )
