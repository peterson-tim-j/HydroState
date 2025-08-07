#' Get seasons
#'
#' \code{get.seasons}
#'
#' @description
#' Aggregates monthly data to 4 seasons in a year.
#'
#' @details
#' This function takes sums monthly runoff and precipitation observations into 4 seasons of a year.
#'
#' @param input.data dataframe of monthly runoff and precipitation observations. Gaps with missing data in either streamflow or precipitation are permitted, and the handling of them is further discussed in \code{build}. Monthly data is required when using \code{seasonal.parameters} that assumes selected model parameters are better defined with a sinusoidal function.
#'
#' @return
#' A dataframe of seasonal observations with an additional column counting the number of months in each season.
#'
#' @keywords seasons
#'
#'
#' @export get.seasons
#'
#' @examples
#'
#' # Load data
#' data(streamflow_monthly_221201)
#'
#' # aggregate monthly data to seasonal
#' streamflow_seasonal_221201 = get.seasons(streamflow_monthly_221201)
#'
#'

get.seasons <- function(input.data=data.frame(year=c(), month = c(), flow=c(), precip=c())){

#Validate
  # WILL NEED TO CHANGE SUM TO AVERAGE, and maybe middle month for concentration...

  if(!('year' %in% colnames(input.data))){
    stop("'input.data' must contain a 'year' column with an integer of years")
  }


  if(!('month' %in% colnames(input.data))){
    stop("'input.data' must contain a 'month' column with an integer of months")
  }

  # If monthly data, sort in ascending order by year and month
  if('month' %in% colnames(input.data)){
    input.data = input.data[order(input.data[,'year'],input.data[,'month']),]
  }else if('day' %in% colnames(input.data)){
    input.data = input.data[order(input.data[,'year'],input.data[,'month'],input.data[,'day']),]
  }

  delta = getStartEndIndex(input.data)


  # Aggregate monthly data to seasons
  season.months = matrix(0,4,2)
  season.months[1,] = c(12,2)
  season.months[2,] = c(3,5)
  season.months[3,] = c(6,8)
  season.months[4,] = c(9,11)
  nrow=0
  season.end.month.current = 0
  #streamflow_monthly.seasonal = data.frame(year=c(), month=c(), flow=c(), precipitation=c(), nmonths=c())
  streamflow_monthly.seasonal = matrix(NA,nrow(input.data),5)

  for(j in 1:nrow(delta)){

    data = input.data[delta[j,1]:delta[j,2],]

    for (i in 1:nrow(data)) {

      # Get season index.
      if (data$month[i]>=season.months[1,1] || data$month[i]<=season.months[1,2]) {
        season.ind = 1
      } else {
        season.ind = which(data$month[i]>=season.months[,1] & data$month[i]<=season.months[,2])
      }
      season.end.month = season.months[season.ind,2]

      if (season.end.month == season.end.month.current) {
        streamflow_monthly.seasonal[nrow,3] = streamflow_monthly.seasonal[nrow,3] + data$flow[i]
        streamflow_monthly.seasonal[nrow,4] = streamflow_monthly.seasonal[nrow,4] + data$precipitation[i]
        streamflow_monthly.seasonal[nrow,5] = streamflow_monthly.seasonal[nrow,5] + 1
      } else {
        nrow = nrow + 1
        if (season.end.month==2) {
          streamflow_monthly.seasonal[nrow,1] = data$year[i]+1
        } else {
          streamflow_monthly.seasonal[nrow,1] = data$year[i]
        }
        streamflow_monthly.seasonal[nrow,2] = season.end.month
        streamflow_monthly.seasonal[nrow,3] = data$flow[i]
        streamflow_monthly.seasonal[nrow,4] = data$precipitation[i]
        streamflow_monthly.seasonal[nrow,5] = 1

        season.end.month.current = season.end.month
      }
    }
  }

  filt = !is.na(streamflow_monthly.seasonal[,1]) & streamflow_monthly.seasonal[,5]==3
  streamflow_monthly.seasonal = streamflow_monthly.seasonal[filt,]
  streamflow_monthly.seasonal = data.frame(year=streamflow_monthly.seasonal[,1], month=streamflow_monthly.seasonal[,2], flow=streamflow_monthly.seasonal[,3], precipitation=streamflow_monthly.seasonal[,4], nmonths=streamflow_monthly.seasonal[,5])

  # Remove the first 5 years if precip to minimise AR1 spin up effects. THIS IS A TEMP FIX!!!
  # filt = streamflow_monthly.seasonal$year <= (min(streamflow_monthly.seasonal$year)+4)
  # streamflow_monthly.seasonal$flow[filt] = NA

  return(streamflow_monthly.seasonal)
}

#' @keywords internal
select.transform <- function(func = 'boxcox', input.data=data.frame(year=c(), flow=c(), precip=c())){

  #Validate
  func = paste('Qhat.',func,sep='')
  if(func %in% c('Qhat.none','Qhat.log','Qhat.burbidge','Qhat.boxcox')){

    # # If monthly data, sort in ascending order by year and month
    # if('month' %in% colnames(input.data)){
    #   input.data = input.data[order(input.data[,'year'],input.data[,'month']),]
    # }else if('day' %in% colnames(input.data)){
    #   input.data = input.data[order(input.data[,'year'],input.data[,'month'],input.data[,'day']),]
    # }


    return(new(func, input.data))

  }else{

    stop('Please choose an appropriate function: none, log, burbidge, or boxcox')
  }


}

#' @keywords internal
select.stateModel <- function(input.data = data.frame(year=c(), flow=c(), precip=c()),
                              parameters = c('a0','a1','std'),
                              seasonal.parameters,
                              state.shift.parameters = c('a0','std'),
                              error.distribution,
                              transition.graph = matrix(TRUE,2,2)){

  #Validate input.data
  if(!is.data.frame(input.data)){
    stop("'input.data' is not a dataframe")
  }
  if(!('year' %in% colnames(input.data))){
    stop("'input.data' must contain a 'year' column with an integer of years")
  }

  #Validate parameters
  if(!is.character(parameters)){
    stop("'parameters' are not in a character vector")
  }

  if(!('a0' %in% parameters && 'a1' %in% parameters && 'std' %in% parameters)){
    stop("'parameters' must contain 'a0', 'a1', 'std'")
  }
  if(('AR1' %in% parameters && 'AR2' %in% parameters && 'AR3' %in% parameters)||
     ('AR1' %in% parameters && 'AR2' %in% parameters)  ||
     ('AR1' %in% parameters && 'AR3' %in% parameters) ||
     ('AR2' %in% parameters && 'AR3' %in% parameters)){
    stop("'parameters' can contain only one auto-correlation term, either 'AR1', 'AR2', or 'AR3'")
  }

  #Validate seasonal.parameters
  if(length(seasonal.parameters) > 0 && !is.character(seasonal.parameters)){
    stop("'seasonal.parameters' are not characters")
  }

  if('AR1' %in% seasonal.parameters || 'AR2' %in% seasonal.parameters || 'AR3' %in% seasonal.parameters){
    stop("'seasonal.parameters' cannot contain auto-correlation terms")
  }
  if(!('month' %in% colnames(input.data)) && length(seasonal.parameters) > 0){
    stop("Please include monthly data with a 'month' column in the 'input.data' for subannual analysis.
         hydroState defaults to subannual analysis when 'seasonal.parameters' are selected, but 'input.data' does not contain subannual data.
         Columns with data for 'month' requried, else consider annual analysis with no 'seasonal.parameters'")
  }

  #Validate state.shift.parameters
  if(!is.character(state.shift.parameters)){
    stop("'state.shift.parameters' are not character vector")
  }

  # If state shift parameters contain auto-correlation terms, but no auto-correlation in model parameters, STOP
  if('AR1' %in% state.shift.parameters || 'AR2' %in% state.shift.parameters || 'AR3' %in% state.shift.parameters){
    if(!(any(parameters == 'AR1') || any(parameters == 'AR2') ||  any(parameters == 'AR3'))){
      stop("'state.shift.parameters' cannot contain auto-correlation terms because 'parameters' in state.model does not contain auto-correlation")
    }
  }

  #Validate error.distribution
  if(!is.character(error.distribution)){
    stop("'error.distribution' is not a character string")
  }
  if(!(error.distribution %in% c('normal','truc.normal','gamma'))){
    stop("'error.distribution' must contain either 'normal', 'truc.normal' or 'gamma'")
  }

  if(length(seasonal.parameters) > 0 && error.distribution == 'normal'){
    stop("Please select 'gamma' for error.distribution. At the moment, hydroState does not support seasonal parameters with 'error.distribution' = 'normal' or 'truc.normal'")

  }else if(length(seasonal.parameters) > 0 && error.distribution == 'truc.normal'){

    stop("Please select 'gamma' for error.distribution. At the moment, hydroState does not support seasonal parameters with 'error.distribution' = 'normal' or 'truc.normal'")
  }

  #Validate transition.graph
  if(!is.matrix(transition.graph)){
    stop("'transition.graph' must be a matrix")
  }

  # If monthly data, sort in ascending order by year and month
  if('month' %in% colnames(input.data)){
    input.data = input.data[order(input.data[,'year'],input.data[,'month']),]
  }else if('day' %in% colnames(input.data)){
    input.data = input.data[order(input.data[,'year'],input.data[,'month'],input.data[,'day']),]
  }

  # Curate state model
  if('AR1' %in% parameters || 'AR2' %in% parameters || 'AR3' %in% parameters){
    auto.array = c(length(which(parameters == 'AR1')) > 0,
                   length(which(parameters == 'AR2')) > 0,
                   length(which(parameters == 'AR3')) > 0)
    auto.term = paste('AR',which(auto.array == TRUE),sep='')

    # Assume seasonal model if seasona paramaters
    if(length(seasonal.parameters) > 0){

      # Add if statement here if error.distribution for seasonal accepts normal or truc.normal; gamma for now
      func = paste('QhatModel.subAnnual.homo.',error.distribution,'.linear.',auto.term,sep='')

    }else{ # assume annual

      # if truncated normal distribution, set item to true
      if(error.distribution == 'truc.normal'){
        use.truncated.dist = TRUE
        func = paste('QhatModel.homo.','normal','.linear.',auto.term,sep='')
      }else{
        use.truncated.dist = FALSE
        func = paste('QhatModel.homo.',error.distribution,'.linear.',auto.term,sep='')
      }

    }
    # Else, no auto-correlation
  }else{

    # Assume seasonal model if seasonal paramaters
    if(length(seasonal.parameters) > 0){

      # Add if statement here if error.distribution for seasonal accepts normal or truc.normal; gamma for now
      func = paste('QhatModel.subAnnual.homo.',error.distribution,'.linear',sep='')

    }else{

      # if truncated normal distribution, set item to true
      if(error.distribution == 'truc.normal'){
        use.truncated.dist = TRUE
        func = paste('QhatModel.homo.','normal','.linear',sep='')
      }else{
        use.truncated.dist = FALSE
        func = paste('QhatModel.homo.',error.distribution,'.linear',sep='')
      }
    }

  }

    # Check parameters that make model, to assign state shifts, no auto-correlation
    if(length(parameters) == 3 && all(parameters %in% c('a0','a1','std'))){

      state.array = c(length(which(state.shift.parameters == 'a0')) > 0,
                      length(which(state.shift.parameters == 'a1')) > 0,
                      length(which(state.shift.parameters == 'std')) > 0)

      # If seasonal parameters, assume seasonal
      if(length(seasonal.parameters) > 0){

        #(one-day, have option to choose state shift for particular sinusoidal param (amp, phase, disp))
        seasonal.array = c(length(which(seasonal.parameters == 'a0')) > 0,
                           length(which(seasonal.parameters == 'a1')) > 0,
                           length(which(seasonal.parameters == 'std')) > 0)

        return(new(func, input.data, transition.graph,
                     state.dependent.mean.a0 = state.array[1],
                     state.dependent.mean.a1 = state.array[2],
                     state.dependent.std.a0 = state.array[3],
                     subAnnual.dependent.mean.a0 = seasonal.array[1],
                     subAnnual.dependent.mean.a1 = seasonal.array[2],
                     subAnnual.dependent.std.a0 = seasonal.array[3]))

        # If not seasonal model, proceed with only option to select state shifts in model parameters
      }else{

        if(error.distribution %in% c('truc.normal','normal')){

          return(new(func, input.data, use.truncated.dist, transition.graph,
                     state.dependent.mean.a0 = state.array[1],
                     state.dependent.mean.a1 = state.array[2],
                     state.dependent.mean.trend=NA,
                     state.dependent.std.a0 = state.array[3]))

        }else{ #gamma does not need use.truncated.dist input

          return(new(func, input.data, transition.graph,
                     state.dependent.mean.a0 = state.array[1],
                     state.dependent.mean.a1 = state.array[2],
                     state.dependent.mean.trend=NA,
                     state.dependent.std.a0 = state.array[3]))
        }
      }
    }

   # AR1 - now the same for auto-correlation, seasonal? seasonal state shifts for amp, phase, or disp? Or seasonal state shifts for all? ELSE non-seasonal state shifts
  if(length(parameters) == 4 && all(parameters %in% c('a0','a1','std','AR1'))){

    state.array = c(length(which(state.shift.parameters == 'a0')) > 0,
                    length(which(state.shift.parameters == 'a1')) > 0,
                    length(which(state.shift.parameters == 'AR1')) > 0,
                    length(which(state.shift.parameters == 'std')) > 0)

    # If seasonal parameters, assume seasonal
    if(length(seasonal.parameters) > 0){

      #(one-day, have option to choose state shift for particular sinusoidal param (amp, phase, disp))
      seasonal.array = c(length(which(seasonal.parameters == 'a0')) > 0,
                         length(which(seasonal.parameters == 'a1')) > 0,
                         length(which(seasonal.parameters == 'std')) > 0)

        return(new(func, input.data, transition.graph,
                   state.dependent.mean.a0 = state.array[1],
                   state.dependent.mean.a1 = state.array[2],
                   state.dependent.mean.AR1 = state.array[3],
                   state.dependent.std.a0 = state.array[4],
                   subAnnual.dependent.mean.a0 = seasonal.array[1],
                   subAnnual.dependent.mean.a1 = seasonal.array[2],
                   subAnnual.dependent.std.a0 = seasonal.array[3]))

      # If not seasonal model, proceed with only option to select state shifts in model parameters
    }else{

      if(error.distribution %in% c('truc.normal','normal')){

        return(new(func, use.truncated.dist, input.data, transition.graph,
                   state.dependent.mean.a0 = state.array[1],
                   state.dependent.mean.a1 = state.array[2],
                   state.dependent.mean.trend=NA,
                   state.dependent.mean.AR1 = state.array[3],
                   state.dependent.std.a0 = state.array[4]))
      }else{ # gamma

        return(new(func, input.data, transition.graph,
                   state.dependent.mean.a0 = state.array[1],
                   state.dependent.mean.a1 = state.array[2],
                   state.dependent.mean.trend=NA,
                   state.dependent.mean.AR1 = state.array[3],
                   state.dependent.std.a0 = state.array[4]))
      }
    }
  }
  # AR2
  if(length(parameters) == 4 && all(parameters %in% c('a0','a1','std','AR2'))){

    state.array = c(length(which(state.shift.parameters == 'a0')) > 0,
                    length(which(state.shift.parameters == 'a1')) > 0,
                    length(which(state.shift.parameters == 'AR1')) > 0,
                    length(which(state.shift.parameters == 'AR2')) > 0,
                    length(which(state.shift.parameters == 'std')) > 0)

    # If seasonal parameters, assume seasonal
    if(length(seasonal.parameters) > 0){

      #(one-day, have option to choose state shift for particular sinusoidal param (amp, phase, disp))
      seasonal.array = c(length(which(seasonal.parameters == 'a0')) > 0,
                         length(which(seasonal.parameters == 'a1')) > 0,
                         length(which(seasonal.parameters == 'std')) > 0)

      return(new(func, input.data, transition.graph,
                 state.dependent.mean.a0 = state.array[1],
                 state.dependent.mean.a1 = state.array[2],
                 state.dependent.mean.AR1 = state.array[3],
                 state.dependent.mean.AR2 = state.array[4],
                 state.dependent.std.a0 = state.array[5],
                 subAnnual.dependent.mean.a0 = seasonal.array[1],
                 subAnnual.dependent.mean.a1 = seasonal.array[2],
                 subAnnual.dependent.std.a0 = seasonal.array[3]))

      # If not seasonal model, proceed with only option to select state shifts in model parameters
    }else{

      if(error.distribution %in% c('truc.normal','normal')){

        return(new(func, input.data,use.truncated.dist, transition.graph,
                   state.dependent.mean.a0 = state.array[1],
                   state.dependent.mean.a1 = state.array[2],
                   state.dependent.mean.trend=NA,
                   state.dependent.mean.AR1 = state.array[3],
                   state.dependent.mean.AR2 = state.array[4],
                   state.dependent.std.a0 = state.array[5]))
      }else{ # gamma

        return(new(func, input.data, transition.graph,
                   state.dependent.mean.a0 = state.array[1],
                   state.dependent.mean.a1 = state.array[2],
                   state.dependent.mean.trend=NA,
                   state.dependent.mean.AR1 = state.array[3],
                   state.dependent.mean.AR2 = state.array[4],
                   state.dependent.std.a0 = state.array[5]))
      }
    }
  }

  # AR3
  if(length(parameters) == 4 && all(parameters %in% c('a0','a1','std','AR3'))){

    state.array = c(length(which(state.shift.parameters == 'a0')) > 0,
                    length(which(state.shift.parameters == 'a1')) > 0,
                    length(which(state.shift.parameters == 'AR1')) > 0,
                    length(which(state.shift.parameters == 'AR2')) > 0,
                    length(which(state.shift.parameters == 'AR3')) > 0,
                    length(which(state.shift.parameters == 'std')) > 0)

    # If seasonal parameters, assume seasonal
    if(length(seasonal.parameters) > 0){

      #(one-day, have option to choose state shift for particular sinusoidal param (amp, phase, disp))
      seasonal.array = c(length(which(seasonal.parameters == 'a0')) > 0,
                         length(which(seasonal.parameters == 'a1')) > 0,
                         length(which(seasonal.parameters == 'std')) > 0)

      return(new(func, input.data, transition.graph,
                 state.dependent.mean.a0 = state.array[1],
                 state.dependent.mean.a1 = state.array[2],
                 state.dependent.mean.AR1 = state.array[3],
                 state.dependent.mean.AR2 = state.array[4],
                 state.dependent.mean.AR3 = state.array[5],
                 state.dependent.std.a0 = state.array[6],
                 subAnnual.dependent.mean.a0 = seasonal.array[1],
                 subAnnual.dependent.mean.a1 = seasonal.array[2],
                 subAnnual.dependent.std.a0 = seasonal.array[3]))

      # If not seasonal model, proceed with only option to select state shifts in model parameters
    }else{

      if(error.distribution %in% c('truc.normal','normal')){

        return(new(func, input.data,  use.truncated.dist, transition.graph,
                   state.dependent.mean.a0 = state.array[1],
                   state.dependent.mean.a1 = state.array[2],
                   state.dependent.mean.trend=NA,
                   state.dependent.mean.AR1 = state.array[3],
                   state.dependent.mean.AR2 = state.array[4],
                   state.dependent.mean.AR3 = state.array[5],
                   state.dependent.std.a0 = state.array[6]))
      }else{ # gamma

        return(new(func, input.data, transition.graph,
                   state.dependent.mean.a0 = state.array[1],
                   state.dependent.mean.a1 = state.array[2],
                   state.dependent.mean.trend=NA,
                   state.dependent.mean.AR1 = state.array[3],
                   state.dependent.mean.AR2 = state.array[4],
                   state.dependent.mean.AR3 = state.array[5],
                   state.dependent.std.a0 = state.array[6]))
      }
    }
  }

}

#' @keywords internal
select.Markov <- function(flickering = FALSE,
                          transition.graph){

  #Validate
  if(flickering == TRUE){

    Markov = new('markov.annualHomogeneous.flickering', transition.graph)

  }else{

    Markov = new('markov.annualHomogeneous',transition.graph)

  }

}

#' Builds hydroState model
#'
#' \code{build}
#'
#' @description
#' \code{build} builds a hydrostate model object with either a default model or the model can be specified with options from below. Every model depends on a linear base model where streamflow, \eqn{Q}, is a function of precipitation, \eqn{P}: \eqn{Q = Pa_1 + a_0}. The default model is a variant of this base linear model with state shifts expected in the intercept, \eqn{a_0}, and standard deviation, \eqn{std}, of the rainfall-runoff relationship. There are additional options to adjust this model with auto-correlation terms and seasonal parameters that can be both independent of state of change with state. The number of states and assumed error distribution can also be selected. After the model is built, the hydroState model is ready to be fitted with \code{fit.hydroState()}
#'
#' @details
#' There are a selection of items to consider when defining the rainfall-runoff relationship and investigating state shifts in this relationship. hydroState provides various options for modelling the rainfall-runoff relationship.
#' \itemize{
#' \item{Data gaps  with \code{input.data}}: When there is missing \code{input.data} in either the dependent variable, streamflow, or independent variable, precipitation, the emissions probability of the missing time-step is set equal to one. This essentially ignores the missing periods. The time step after the missing period has a state probability dependent on the length of the gap. The larger the gap, the closer the state probability gets to approaching a finite probability near zero (as the transition probabilities are recursively multiplied). When the model has auto-correlation terms and there are gaps in the dependent variable, the auto-correlation function restarts at the beginning of each continuous period after the gap. This ignores auto-correlation at the first time steps after the gap. For instance, an 'AR1' model would ignore the contribution of the prior time step for the first (1) observation after the gap.
#' \item{Transform Observations with \code{data.transform}}: Transforms streamflow observations to remove heteroscedasticity. Often there is skew within hydrologic data. When defining relationships between rainfall-runoff, this skew results in an unequal variance in the residuals, heteroscedasticity. Transforming streamflow observations is often required. There are several options to transform observations. Since the degree of transformation is not typically known, \code{boxcox} is the default. Other options include: \code{log}, \code{burbidge}, and of course, \code{none} when no transformation is performed.
#' \item{Model Structure with \code{parameters} and \code{seasonal.parameters}}: The structure of the model depends on the \code{parameters}. hydroState simulates runoff, \eqn{Q}, as being in one of a finite states, \eqn{i}, at every time-step, \eqn{t}, depending on the distribution of states at prior time steps. This results in a runoff distribution for each state that can vary over time (\eqn{\widehat{_tQ_i}}). The model defines the relationship that is susceptible to state shifts with precipitation, \eqn{P_t}, as a predictor. This takes the form as a simple linear model \eqn{\widehat{_tQ_i} = f(P_t)}:
#'
#'
#' \eqn{\widehat{_tQ_i} = P_ta_1 + a_0}
#'
#' where \eqn{a_0} and \eqn{a_1} are constant parameters. These parameters and the model error, \eqn{std}, establish the rainfall-runoff relationship and are required parameters for every built model object. These are the default parameters: c(``a0",``a1",``std").
#' \item{Auto-correlation \code{parameters}}: The relationship may contain serial correlation and would be better defined with an auto-regressive term:
#'
#' \eqn{\widehat{_tQ_i} = P_ta_1 + a_0 + AR1\widehat{_{t-1}Q}}
#'
#' where \eqn{AR1} is the lag-1 auto-correlation term. Either, lag-1: \eqn{AR1}, lag-2: \eqn{AR2}, and lag-3: \eqn{AR3} auto-correlation coefficients are an option as additional parameters to better define the rainfall-runoff relationship.
#' \item{Sub-annual analysis with \code{seasonal.parameters}}: Additional options include explaining the seasonal rainfall-runoff relationship with a sinusoidal function that better defines either of the constant parameters or error (\eqn{a_0, a_1, std}) throughout the year, i.e:
#'
#' \eqn{a_0 = a_{0.disp} + a_{0.amp} * sin(2\pi(\frac{M_t}{12} + a_{0.phase}))}
#'
#' where \eqn{M_t} is an integer month at \eqn{t}. Monthly streamflow and precipitation are required as \code{input.data} for the sub-annual analysis.
#'
#' \item{State Dependent Parameters with \code{state.shift.parameters}}: These are state dependent parameters where they are subject to shift in order to better explain the state of streamflow over time. Any or all of the previously chosen parameters can be selected (\eqn{a_0, a_1, std, AR1, AR2, AR3}). The default model evaluates shifts in the rainfall-runoff relationship with \eqn{a_0} \eqn{std} as state dependent parameters.
#' \item{Distribution of the Residuals with \code{error.distribution}}: The distribution of the residuals (error) within a state of the model can be chosen to reduce skew and assist with making models statistically adequate (see \code{plot(pse.residuals = TRUE)}). Either normal: \code{normal}, truncated normal: \code{truc.normal}, or gamma: \code{gamma} distributions are acceptable. These error distribution ensures streamflow is greater than zero \eqn{Q > 0}, and specifically for \code{truc.normal} greater than or equal to zero \eqn{Q >= 0}. The default is \code{truc.normal}. Sub-annual models are restricted to only a \code{gamma} distribution.
#' \item{Markov flickering with \code{flickering}}: When flickering is \code{FALSE}, the markov avoids state shifts for very short duration, and hence for a state shift to occur it should last for an extended period. The default is FALSE. If TRUE, flickering between more states is more sensitive. For further explanation on this method, see: \href{https://hess.copernicus.org/articles/7/652/2003/}{Lambert et al., 2003}. The current form of the Markov model is homogeneous where the transition probabilities are time-invariant.
#' \item{Number of States with \code{transition.graph}}: The number of possible states in the rainfall-runoff relationship and transition between the states is selected with the transition.graph. The default is a 2-state model in a 2 by 2 unstructured matrix with a TRUE transition to and from each state (i.e. \code{matrix(TRUE,2,2)}). hydroState accepts 1-state up to 3-states (i.e. for 3-state unstructured transition graph: \code{matrix(TRUE,3,3)}). The unstructured transition graph allows either state to remain in the current state or transition between any state. For the 3-state transition graph, one may want to assume the transitions can only occur in a particular order, as in (Very low -> Low -> Normal->) rather than (Very low <-> Low <-> Normal <-> Very low). Thus, a structured graph is also acceptable (3-state structured: \code{matrix(c(TRUE,TRUE,FALSE,FALSE,TRUE,TRUE,TRUE,FALSE,TRUE),3,3)}). More details in Supplementary Materials of (Peterson TJ, Saft M, Peel MC & John A (2021), Watersheds may not recover from drought, Science, DOI: \doi{10.1126/science.abd5085}).
#'}
#'
#' @param input.data dataframe of annual, seasonal, or monthly runoff and precipitation observations. Gaps with missing data in either streamflow or precipitation are permitted. Monthly data is required when using \code{seasonal.parameters}.
#' @param data.transform character sting with the method of transformation. The default is 'boxcox'. Other options: 'log', 'burbidge', 'none'
#' @param parameters character vector of parameters to construct model. Required and default: \code{a0}, \code{a1}, \code{std}. Auto-correlation terms optional: \code{AR1}, \code{AR2}, or \code{AR3}.
#' @param seasonal.parameters character vector of one or all parameters (\code{a0}, \code{a1}, \code{std}) defined as a sinusoidal function to represent seasonal variation. Requires monthly or seasonal data. Default is empty, no seasonal parameters.
#' @param state.shift.parameters character vector of one or all parameters (\code{a0}, \code{a1}, \code{std}, \code{AR1}, \code{AR2}, \code{AR3}) able to shift as dependent on state. Default is \code{a0} and \code{std}.
#' @param error.distribution character string of the distribution in the HMM error. Default is 'truc.normal'. Others include: 'normal' or 'gamma'
#' @param flickering logical \code{TRUE}/\code{FALSE}. \code{TRUE} = allows more sensitive markov flickering between states over time. When \code{FALSE} (default), state needs to persist for at least three time steps before state shift can occur.
#' @param transition.graph matrix given the number of states. Default is a 2-state matrix (2 by 2): \code{matrix(TRUE,2,2)}.
#'
#' @return
#' A built hydroState model object ready to be fitted with \code{fit.hydroState()}
#'
#' @keywords hydroState build
#'
#'
#' @export
#' @importFrom methods new
#'
#' @examples
#' # Load data
#' data(streamflow_annual_221201)
#'
#' ## Build default annual hydroState model
#' model = build(input.data = streamflow_annual_221201)
#'
#' # OR
#'
#' ## Build annual hydroState model with specified objects

#' # Build hydroState model with: 2-state, normal error distribution,
#' # 1-lag of auto-correlation, and state dependent parameters ('a1', 'std')
#' model = build(input.data = streamflow_annual_221201,
#'                    data.transform = 'boxcox',
#'                    parameters = c('a0','a1','std','AR1'),
#'                    state.shift.parameters = c('a1','std'),
#'                    error.distribution = 'normal',
#'                    flickering = FALSE,
#'                    transition.graph = matrix(TRUE,2,2))
#'



build <- function(input.data = data.frame(year=c(), flow=c(), precip=c()),
                       data.transform = 'boxcox',
                       parameters = c('a0','a1','std'),
                       seasonal.parameters = NULL,
                       state.shift.parameters = c('a0','std'),
                       error.distribution = 'truc.normal',
                       flickering = FALSE,
                       transition.graph = matrix(TRUE,2,2)){


  #Validate input.data
  if(!is.data.frame(input.data)){
    stop("'input.data' is not a dataframe")
  }
  if(!('year' %in% colnames(input.data))){
    stop("'input.data' must contain a 'year' column with an integer of years")
  }

  # If monthly data, sort in ascending order by year and month, assume calender year given
  if('month' %in% colnames(input.data)){
    input.data = input.data[order(input.data[,'year'],input.data[,'month']),]
  }else if('day' %in% colnames(input.data)){
    input.data = input.data[order(input.data[,'year'],input.data[,'month'],input.data[,'day']),]
  }

  # default flickering to FALSE if NULL
  if(is.null(flickering))
    flickering = FALSE

  ######################################################

  #######################################################

  # default transition graph if NULL
  if(is.null(transition.graph))
    transition.graph = matrix(TRUE,2,2)

  # default data.transform if NULL
  if(missing(data.transform)){
    data.transform = select.transform(func = 'boxcox', input.data)
  }else{
    data.transform = select.transform(data.transform, input.data)
  }

  # default state model error.distribution if NULL
  if(missing(error.distribution) && length(seasonal.parameters) < 1){
    error.distribution = 'truc.normal'
  }else if(missing(error.distribution) && length(seasonal.parameters) > 0){
    error.distribution = 'gamma'
  }

# If no inputs except input.data, just run default analysis
    if(is.null(parameters) && is.null(seasonal.parameters) && is.null(state.shift.parameters)){

      # create state model
      stateModel = select.stateModel(input.data,
                                     parameters = c('a0','a1','std'),
                                     seasonal.parameters,
                                     state.shift.parameters = c('a0','std'),
                                     error.distribution = error.distribution,
                                     transition.graph = transition.graph)

      # create Markov model
      Markov = select.Markov(flickering, transition.graph)

      # build default model
      return(new('hydroState',input.data, data.transform, stateModel, Markov))

    }else if(!is.null(parameters) && is.null(seasonal.parameters) && is.null(state.shift.parameters)){

      # create state model
      stateModel = select.stateModel(input.data,
                                     parameters,
                                     seasonal.parameters,
                                     state.shift.parameters = c('a0','std'),
                                     error.distribution = error.distribution,
                                     transition.graph = transition.graph)

      # create Markov model
      Markov = select.Markov(flickering, transition.graph)

      # build default model
      return(new('hydroState',input.data, data.transform, stateModel, Markov))

    }else if(!is.null(parameters) && !is.null(seasonal.parameters) && is.null(state.shift.parameters)){


      # create state model
      stateModel = select.stateModel(input.data,
                                     parameters,
                                     seasonal.parameters,
                                     state.shift.parameters = c('a0','std'),
                                     error.distribution = error.distribution,
                                     transition.graph = transition.graph)

      # create Markov model
      Markov = select.Markov(flickering, transition.graph)

      # build default model
      return(new('hydroState',input.data, data.transform, stateModel, Markov))

    }else if(!is.null(parameters) && !is.null(seasonal.parameters) && !is.null(state.shift.parameters)){

      # create state model
      stateModel = select.stateModel(input.data,
                                     parameters,
                                     seasonal.parameters,
                                     state.shift.parameters,
                                     error.distribution = error.distribution,
                                     transition.graph = transition.graph)

      # create Markov model
      Markov = select.Markov(flickering, transition.graph)

      # build default model
      return(new('hydroState',input.data, data.transform, stateModel, Markov))

    }else if(is.null(parameters) && is.null(seasonal.parameters) && !is.null(state.shift.parameters)){

      # create state model
      stateModel = select.stateModel(input.data,
                                     parameters = c('a0','a1','std'),
                                     seasonal.parameters,
                                     state.shift.parameters,
                                     error.distribution = error.distribution,
                                     transition.graph = transition.graph)

      # create Markov model
      Markov = select.Markov(flickering, transition.graph)

      # build default model
      return(new('hydroState',input.data, data.transform, stateModel, Markov))

    }else if(is.null(parameters) && !is.null(seasonal.parameters) && !is.null(state.shift.parameters)){

      # create state model
      stateModel = select.stateModel(input.data,
                                     parameters = c('a0','a1','std'),
                                     seasonal.parameters,
                                     state.shift.parameters,
                                     error.distribution = error.distribution,
                                     transition.graph = transition.graph)

      # create Markov model
      Markov = select.Markov(flickering, transition.graph)

      # build default model
      return(new('hydroState',input.data, data.transform, stateModel, Markov))

    }else if(!is.null(parameters) && is.null(seasonal.parameters) && !is.null(state.shift.parameters)){

      # create state model
      stateModel = select.stateModel(input.data,
                                     parameters,
                                     seasonal.parameters,
                                     state.shift.parameters,
                                     error.distribution = error.distribution,
                                     transition.graph = transition.graph)

      # create Markov model
      Markov = select.Markov(flickering, transition.graph)

      # build default model
      return(new('hydroState',input.data, data.transform, stateModel, Markov))

    }else if(is.null(parameters) && !is.null(seasonal.parameters) && is.null(state.shift.parameters)){

      # create state model
      stateModel = select.stateModel(input.data,
                                     parameters = c('a0','a1','std'),
                                     seasonal.parameters,
                                     state.shift.parameters = c('a0','std'),
                                     error.distribution = error.distribution,
                                     transition.graph = transition.graph)

      # create Markov model
      Markov = select.Markov(flickering, transition.graph)

      # build default model
      return(new('hydroState',input.data, data.transform, stateModel, Markov))

    }

}


#' Builds all hydroState models
#'
#' \code{build.all}
#'
#' @description
#' \code{build.all} builds all possible combinations of hydroState models. The same fields are available as in \code{build()} in order to specify the type of models to be built. After all models are built, they are fitted using the same \code{fit.hydroState()} function.
#'
#' @details
#' All possible combinations of hydroState models are built for each auto-correlation lag and residual distribution from 1 to 3 states for a specified data transformation. This allows for investigation of state changes in the \code{state.shift.parameters}: the intercept \code{c('a0', 'std')} or slope \code{c('a1', 'std')}. To reduce the number of models in the search, specify which field(s) to remain constant. For example, to investigate the best model with the number of auto-correlation terms and number of states with a \code{boxcox} data transform and \code{gamma} distribution of the residuals, set \code{data.transform} to \code{boxcox} and \code{error.distribution} to \code{gamma}. If no fields are specified, all possible model combinations are built. If investigating state shifts in the intercept \code{a0} and slope \code{a1}, it is recommended to build and fit the model combinations separately.
#'
#' @param input.data dataframe of annual, seasonal, or monthly runoff and precipitation observations. Gaps with missing data in either streamflow or precipitation are permitted, and the handling of them is further discussed in \code{build}. Monthly data is required when using \code{seasonal.parameters} that assumes selected model parameters are better defined with a sinusoidal function.
#' @param data.transform character string of method of transformation. If empty, the default builds all possible combinations of models with \code{boxcox} data transformation.
#' @param parameters character vector of parameters to determine model form. If empty, the default builds all possible combinations of model forms.
#' @param seasonal.parameters character vector of parameters with sinusoidal function to represent seasonal variation. Requires monthly or seasonal data. If empty and monthly or seasonal data is given, the default builds all possible combinations of models with a seasonal parameter for each and all parameters.
#' @param state.shift.parameters character vector of one or all parameters to identify state dependent parameters. Only one set of parameters permitted. If empty, the default builds all possible model combinations with \code{c('a0','std')} as state shift parameters.
#' @param error.distribution character string of the distribution in the HMM error. If empty, the default builds models with all possible combinations of error distribution: \code{c('truc.normal', 'normal','gamma')}
#' @param flickering logical \code{TRUE}/\code{FALSE}. \code{TRUE} = allows more sensitive markov flickering between states over time. When \code{FALSE} (default), state needs to persist for at least three time steps before state shift can occur.
#' @param \item{transition.graph}{
#'      matrix given the number of states. If empty, the default builds models with all possible combinations of states:
#'      1-state matrix (1 by 1): \code{matrix(TRUE,1,1)},
#'      2-state matrix (2 by 2): \code{matrix(TRUE,2,2)},
#'      3-state matrix (3 by 3): \code{matrix(TRUE,3,3)}.
#'      }
#' @param siteID character string of site identifier.
#' @param summary.table data frame with a table summarizing all built models and corresponding reference model. From function \code{summary()}. If empty, summary table will be built automatically.
#'
#'
#' @return
#' A list of built hydroState models with every combination of objects ready to be fitted
#'
#' @keywords hydroState build all
#'
#'
#' @export
#' @importFrom methods new
#'
#' @examples
#' # Load data
#' data(streamflow_annual_221201)
#'
#' # Build all annual models with state shift in intercept 'a0'
#' all.annual.models = build.all(input.data = streamflow_annual_221201,
#'                                     state.shift.parameters = c('a0','std'),
#'                                      siteID = '221201')
#'
#' # OR
#'
#' # Build all annual models with state shift in slope 'a1'
#' all.annual.models = build.all(input.data = streamflow_annual_221201,
#'                                     state.shift.parameters = c('a1','std'),
#'                                      siteID = '221201')
#'


build.all <-function(input.data = data.frame(year=c(), flow=c(), precip=c()),
                         data.transform = NULL,
                         parameters = NULL,
                         seasonal.parameters = NULL,
                         state.shift.parameters = NULL,
                         error.distribution = NULL,
                         flickering = FALSE,
                         transition.graph = NULL,
                         summary.table = NULL,
                         siteID = NULL){

  #######################################################

  #called from build()


  #site id
  if(is.null(siteID)){
    siteID = ""
  }

  # get data.transform (all possible options)

  if(is.null(data.transform)){
    data.transform = list('boxcox')
    names(data.transform) = c("boxcox")
  }else if(length(data.transform)>1){
    stop("Only one data.transform is permitted.")
  }else{
    data.transform = list(data.transform)
  }

  if(is.null(parameters)){
    parameters = list(list('a0','a1','std'), list('a0','a1','std','AR1'), list('a0','a1','std','AR2'), list('a0','a1','std','AR3'))
    names(parameters) = c("AR0","AR1","AR2","AR3")
  }else if(length(parameters)>length(unique(parameters))){
    stop("parameters requiers only one set of unique paremeters as character vector (i.e. c('param1', 'param2')) per build.all")
  }else{
    parameters = list(parameters)
    names(parameters) = ifelse("AR1" %in% unlist(parameters),"AR1",ifelse("AR2" %in% unlist(parameters),"AR2", ifelse("AR3" %in% unlist(parameters), "AR3","AR0")))
  }

  if(is.null(transition.graph)){
    transition.graph = list(matrix(TRUE,1,1),matrix(TRUE,2,2),matrix(TRUE,3,3))#,matrix(c(TRUE,TRUE,FALSE,FALSE,TRUE,TRUE,TRUE,FALSE,TRUE),3,3))
    names(transition.graph) = c("1State","2State","3State")
  }else{
    transition.graph = list(transition.graph)
    names(transition.graph) = ifelse(transition.graph == matrix(TRUE,1,1),"1State", ifelse(transition.graph == matrix(TRUE,2,2), "2State",ifelse(transition.graph == matrix(TRUE,3,3),"3State",paste(NCOL(transition.graph),"UserState",sep=""))))
  }

  if(is.null(state.shift.parameters)){
    state.shift.parameters = list(c('a0','std'))
    names(state.shift.parameters) = paste(unlist(state.shift.parameters), collapse = '')
  }else if(length(state.shift.parameters)>length(unique(state.shift.parameters))){
    stop("state.shift.parameters requiers only one set of unique parameters as character vector (i.e. c('param1', 'param2')) per build.all")
  }else{
    state.shift.parameters = list(state.shift.parameters)
    names(state.shift.parameters) =  paste(unlist(state.shift.parameters), collapse = '')  #ifelse("a0" %in% unlist(state.shift.parameters),"a0",ifelse("a1" %in% unlist(state.shift.parameters),"a1", ""))

  }

  if('month' %in% colnames(input.data)){

    if(length(seasonal.parameters) < 1){
      seasonal.parameters = list(list('a0'), list('a0','std'),list('a0','a1','std'),
                                 list('a1'), list('a1','std'), list('a0','a1'),list('std'))
      names(seasonal.parameters) =list('a0.seasonal', 'a0.std.seasonal', 'a0.a1.std.seasonal',
                                       'a1.seasonal', 'a1.std.seasonal', 'a0.a1.seasonal', 'std.seasonal')
    }else{
      seasonal.parameters = list(seasonal.parameters)
      names(seasonal.parameters) = ifelse("a0" %in% unlist(seasonal.parameters),"a0.seasonal",
                                          ifelse("a1" %in% unlist(seasonal.parameters),"a1.seasonal",
                                                 ifelse("a0" %in% unlist(seasonal.parameters) && "a1" %in% unlist(seasonal.parameters) ,"a0.a1.seasonal",
                                                        ifelse("a0" %in% unlist(seasonal.parameters) && "std" %in% unlist(seasonal.parameters) ,"a0.std.seasonal",
                                                               ifelse("a1" %in% unlist(seasonal.parameters) && "std" %in% unlist(seasonal.parameters) ,"a1.std.seasonal",
                                                                      ifelse("a0" %in% unlist(seasonal.parameters) && "a1" %in% unlist(seasonal.parameters) && "std" %in% unlist(seasonal.parameters) ,"a0.a1.std.seasonal",
                                                                             ifelse("std" %in% unlist(seasonal.parameters),"std.seasonal","")))))))

    }
  }

  if('month' %in% colnames(input.data)){

    if(length(seasonal.parameters) > 0){

      error.distribution = list('gamma')
      message('Subannual models built with gamma error.distribution')

    }else if(missing(error.distribution)){
      error.distribution = list('truc.normal','normal','gamma')
    }else{
      error.distribution = list(error.distribution)
    }

  }else{

    if(missing(error.distribution)){
      error.distribution = list('truc.normal','normal','gamma')
    }else{
      error.distribution = list(error.distribution)
    }

  }

  # build all models

  if(length(seasonal.parameters) > 0){
    build.all.model = vector(mode = "list", length = length(transition.graph) * length(error.distribution) * length(state.shift.parameters) * length(seasonal.parameters) * length(data.transform))
    build.all.model.names = vector(mode = "list", length = length(transition.graph) * length(error.distribution) * length(state.shift.parameters) * length(seasonal.parameters) * length(data.transform))

    build.all.model.count = 1

    # seasonal first... if seasonal parameters
    for(i in 1:length(parameters)){

      temp.parameters = parameters[[i]]

      temp.parameters.name = names(parameters[i])

      for(j in 1:length(data.transform)){

        temp.data.transform = data.transform[[j]]

        for(k in 1:length(state.shift.parameters)){

          temp.state.shift.parameters = state.shift.parameters[[k]]

          temp.state.shift.parameters.name = names(state.shift.parameters[k])

          for(m in 1:length(seasonal.parameters)){

            temp.seasonal.parameters = seasonal.parameters[[m]]

            temp.seasonal.parameters.name = names(seasonal.parameters[m])

            for(z in 1:length(error.distribution)){

              temp.error.distribution = error.distribution[[z]]

              for(y in 1:length(transition.graph)){

                temp.transition.graph = transition.graph[[y]]

                temp.transition.graph.name = names(transition.graph[y])


                build.all.model[[build.all.model.count]] = hydroState::build(input.data,
                                                                      data.transform = temp.data.transform,
                                                                      parameters = unlist(temp.parameters),
                                                                      seasonal.parameters = unlist(temp.seasonal.parameters),
                                                                      state.shift.parameters = unlist(temp.state.shift.parameters),
                                                                      error.distribution = unlist(temp.error.distribution),
                                                                      flickering = flickering,
                                                                      transition.graph = temp.transition.graph)

                build.all.model.names[[build.all.model.count]] <- paste("model.",temp.transition.graph.name,".",temp.error.distribution,".",temp.data.transform,".",temp.parameters.name,".",temp.seasonal.parameters.name,".",temp.state.shift.parameters.name,sep="")

                build.all.model.count = build.all.model.count + 1

              }
            }
          }
        }
      }
    }

    # if not seasonal, just do annual..
  }else{
    build.all.model = vector(mode = "list", length = length(transition.graph) * length(error.distribution) * length(state.shift.parameters) * length(data.transform))
    build.all.model.names = vector(mode = "list", length = length(transition.graph) * length(error.distribution) * length(state.shift.parameters) * length(data.transform))

    build.all.model.count = 1

    for(i in 1:length(parameters)){

      temp.parameters = parameters[[i]]

      temp.parameters.name = names(parameters[i])

      for(j in 1:length(data.transform)){

        temp.data.transform = data.transform[[j]]

        for(k in 1:length(state.shift.parameters)){

          temp.state.shift.parameters = state.shift.parameters[[k]]

          temp.state.shift.parameters.name = names(state.shift.parameters[k])

          for(z in 1:length(error.distribution)){

            temp.error.distribution = error.distribution[[z]]

            for(y in 1:length(transition.graph)){

              temp.transition.graph = transition.graph[[y]]

              temp.transition.graph.name = names(transition.graph[y])


              build.all.model[[build.all.model.count]] = hydroState::build(input.data,
                                                                    data.transform = temp.data.transform,
                                                                    parameters = unlist(temp.parameters),
                                                                    state.shift.parameters = unlist(temp.state.shift.parameters),
                                                                    error.distribution = unlist(temp.error.distribution),
                                                                    flickering = flickering,
                                                                    transition.graph = temp.transition.graph)

              build.all.model.names[[build.all.model.count]] <- paste("model.",temp.transition.graph.name,".",temp.error.distribution,".",temp.data.transform,".",temp.parameters.name,".",temp.state.shift.parameters.name,sep="")

              build.all.model.count = build.all.model.count + 1

            }
          }
        }
      }
    }
  }

  build.all.model <- setNames(build.all.model, unlist(build.all.model.names))

  # make all.models a all.Models object
  build.all.model <- new('hydroState.allModels',models = build.all.model, siteID = siteID)

  if(is.null(summary.table)){

    build.all.model <- get.summary.table(build.all.model)

  }else{

    build.all.model <- get.summary.table(build.all.model, summary.table)
  }
  # now create summary table and assign calibration items


  # build.all.model.table <- showModelAll(build.all.model)
  # order and make a hydroState.allModels object
  # debugonce("get.summary.table","hydroState.allModels")

  return(build.all.model)

}

#' Summarize all models
#'
#'\code{summary}
#'
#' @description
#' \code{summary} outputs a summary table of all built hydroState models and allows users to edit the reference models for calibration.
#'
#' @details
#' For every model object in \code{build.all}, there is a reference model for calibration. The reference model is a slightly simpler model with one less parameter. During calibration with \code{fit}, the model performance must exceed the the performance of the reference model else the model is rejected. This function is used to output the \code{summary.table} and adjust the reference models. Afterwards, all models can be re-build with including this \code{summary.table} in the \code{build.all} function.
#'
#'
#' @param object hydroState.allModels object with a list of models from \code{build.all}
#' @param ... No additional input required
#'
#' @return
#' A data.frame with a summary table of all models and reference models
#'
#' @keywords hydroState summary all models
#'
#' @export
#'
#' @examples
#' # Show summary table of all fitted model details and reference models
#'
#' \dontrun{
#' all.models.ref.table = summary(all.models)
#' }

summary.hydroState.allModels <- function(object, ...){

  all.models = object


  if(class(all.models)[1] == 'hydroState.allModels'){

    all.models.summary <- get.summary.table(all.models, all.models@models.summary)

    return(as.data.frame(all.models.summary@models.summary))

  }else{

    stop("'all.models' is not a 'hydroState.allModels' object")
  }
}


#' Fit hydroState model
#'
#' \code{fit.hydroState}
#'
#' @description
#' \code{fit.hydroState} fits a single hydroState model (\code{build}) or multiple models (\code{build.all}) using global optimization by differential evolution \href{https://cran.r-project.org/package=DEoptim}{DEoptim} library. If fitting all models be sure to install and load the \href{https://cran.r-project.org/package=parallelly}{parallelly} library. The fitting of all models may take hours or days, but the calibration can occur in parallel if the \href{https://cran.r-project.org/package=parallelly}{parallelly} library is installed and loaded.
#'
#' @details
#' After a hydroState model object is built, the model is ready to fit to the observed streamflow through minimizing the negative log-likelihood function too calibrate model parameters. The only required input is the given built hydroState model object or hydroState.allModels object (all models from \code{build.all = TRUE}). When fitting all models, the models are fitted from least to most complex (least to max amount of parameters). Each model has a minimum of 5 and maximum of 20 calibration attempts to outperform the prior reference model else the model is rejected. For instance 'model.1State.normal.log.AR0' contains 3-parameters and is the reference model for 'model.1State.normal.log.AR1' which contains 4-parameters. The objective function of 'model.1State.normal.log.AR1' must calibrate the model with a lower negative log-likelihood than 'model.1State.normal.log.AR0'. These reference models are pre-defined, but this function allows the user to edit the reference models in the data.frame if needed using \code{summary}. Details on the likelihood function is as follows:
#'
#' The likelihood function is estimated as:
#'
#' \eqn{L_{T} = \delta P(x_{1}) + \Gamma \delta P(x_{2})...\Gamma \delta P(x_{T})1'}
#'
#' where:
#' \itemize{
#'  \item{\eqn{\delta}}{ is the initial state distribution, the initial probability of being in each state: \eqn{\delta = \begin{pmatrix} \delta_{1} \\ 1- \delta_{1} \end{pmatrix}}}
#'  \item{\eqn{P(x)}}{ is the \eqn{m} x \eqn{m} diagonal emissions matrix of the probability density for each state using a lower tail truncated Gaussian distribution or a two-parameter Gamma distribution}
#'    \itemize{
#'    \item{}{ \eqn{f_{Gau}(x=\widehat{_{obs}q_{t}}; \mu = \widehat{_{t}q_{i}}, \sigma = \sigma_{i}, a = 0) = \frac{\phi(\frac{x-\mu}{\sigma})}{\sigma(1-\Phi(\frac{a-\mu}{\sigma}))}}}
#'    \item{}{ \eqn{f_{Gam}(x = \widehat{_{obs}q_{t}};k = \frac{\widehat{_{t}q_{i}}^2}{\sigma_{i}^2}, \theta = \frac{\sigma_{i}^2}{\widehat{_{t}q_{i}}}) = \frac{x^{k-1}e^{\frac{x}{\theta}}}{\theta^{k}\Gamma(k)}}}
#'    \item{}{ where \eqn{\phi} is the probability density function for the standard normal distribution, \eqn{\Phi} is the cumulative distribution function for the standard normal distribution, \eqn{k} is the shape parameter, \eqn{\theta} is the scale parameter, and \eqn{\Gamma(k)} is the gamma function.
#'    For more details, refer to pg. 8-17 in Supplementary Materials of (Peterson TJ, Saft M, Peel MC & John A (2021), Watersheds may not recover from drought, Science, DOI: \doi{10.1126/science.abd5085}).}
#'    }
#'  \item{\eqn{\Gamma} is the transition matrix}
#'  \item{\eqn{T} is the number of time-steps.}
#'  }
#'
#' @param model built hydroState model object, hydroState.allModels object, or hydroState.subAnnual.allModels object
#' @param pop.size.perParameter integer that should be greater than or equal to the number of parameters in the model. The default is '10' and is sufficient for all models.
#' @param max.generations integer that will stop the optimizer when set number of generations are reached. The default is '500'.
#' @param doParallel TRUE/FALSE to perform fitting in parallel on all computer cores. Default is FALSE
#'
#' @return
#' A fitted hydroState model
#'
#' @keywords hydroState fit
#'
#'
#' @export
#'
#' @import DEoptim
#' @import methods
#' @import stats
#' @import truncnorm
#' @importFrom utils relist
#'
#'
#'
#' @examples
#'
#' # Load data
#' data(streamflow_annual_221201)
#'
#' ## Build default annual hydroState model
#' model = build(input.data = streamflow_annual_221201)
#'
#' ## Fit built model (runtime ~ 14 sec)
#' \dontrun{
#' model = fit.hydroState(model)
#' }
#'
#' ## Fit all built models (runtime > several hours)
#' # Load data
#' data(streamflow_annual_221201)
#'
#' ## Build all annual models
#' all.annual.models = build.all(input.data = streamflow_annual_221201, siteID = '221201')
#'
#' \dontrun{
#' # Fit all (runtime > several hours)
#' all.annual.models = fit.hydroState(all.annual.models)
#'
#' }
#'

fit.hydroState <- function(model,
                pop.size.perParameter = 10,
                max.generations = 500,
                doParallel = FALSE
                ) {

  if (methods::is(model, "hydroState")) {
    return(fit(model,
                        DEstrategy = 3,
                        pop.size.perParameter = pop.size.perParameter,
                        max.generations = max.generations,
                        Domains = NA,
                        reltol = 1e-8,
                        steptol = 50,
                        print.iterations = 25,
                        use.initial.parameters = FALSE,
                        doParallel = doParallel
    ))
  }

  if (methods::is(model, "hydroState.allModels") ||
      methods::is(model, "hydroState.subAnnual.allModels")) {
    return(fit(model,
                        DEstrategy = 3,
                        pop.size.perParameter = pop.size.perParameter,
                        max.generations = 10000,
                        Domains = NA,
                        reltol = 1e-8,
                        steptol = 50,
                        print.iterations = 25,
                        use.initial.parameters = FALSE,
                        doParallel = doParallel
    ))
  }

  stop("Invalid model class. Expecting 'hydroState', 'hydroState.allModels', or 'hydroState.subAnnual.allModels'.")
}




#' Get AIC
#'
#' \code{get.AIC}
#'
#' @description
#' \code{get.AIC} retrieves Akaike information criteria from a fitted hydroState model object or all models.
#'
#' @details
#' The AIC is the negative log-likelihood of the model plus a penalty for model parameters. This function can be performed on a single model or a selection of models to find the lowest AIC of the set.
#'
#' @param model fitted \code{hydroState} model object.
#'
#' @return AIC value of a single model or a list variable of AIC values for al models
#'
#' @keywords AIC
#'
#'
#' @export get.AIC
#'
#' @examples
#' # Load fitted model
#' data(model.annual.fitted.221201)
#'
#' ## AIC of a single model
#' get.AIC(model.annual.fitted.221201)
#'
#' ## Lowest AIC of a model set
#' get.AIC(all.models.annual.fitted.407211)
#'


get.AIC <- function(model){


  if(!(class(model)[1] %in% c("hydroState","hydroState.allModels", "hydroState.subAnnual.allModels")))
    stop('model is not an appropriate class. Please ensure input model is a built hydroState model with the class as either of the following: "hydroState", "hydroState.allModels", or "hydroState.subAnnual.allModels"')


  # Validate
  if((class(model)[1] %in% c("hydroState.allModels", "hydroState.subAnnual.allModels"))){

      return(getAIC.bestModel(model)$AIC)

  }else{


      return(getAIC(model))

  }

}



#' Get pseudo residuals
#'
#' \code{get.residuals}
#'
#' @description
#' The pseudo residuals were derived from the conditional probabilities of the observations. At each time-step, the pseudo residual is the probability of an observation occurring given the prior observations and latter observations.
#'
#' @details
#' \code{get.residuals} retrieves residuals from the fitted model and exports them as a data frame.
#'
#'
#' @param model fitted hydroState model object.
#'
#' @return
#' Data frame of residuals for each time-step
#'
#' @keywords residuals
#'
#'
#' @export get.residuals
#'
#' @examples
#' # Load fitted model
#' data(model.annual.fitted.221201)
#'
#' ## Get residuals in a dataframe
#' get.residuals(model = model.annual.fitted.221201)
#'
#'


get.residuals <- function(model){

    return(check.PseudoResiduals(model, do.plot = F))

}



#'Sets state names given initial year
#'
#' \code{setInitialYear}
#'
#' @description
#' sets the state names for each time-step relative to the initial year given
#'
#' @details
#' hydroState assigns names to the computed states. This requires choosing an initial year where the state value from that year will be named 'Normal'. Other state values will be given names relative to the state value in the initial year. The choice of the initial year does not affect results. It is a means to more easily interpret the difference in state values relative to each other. It is best to choose a year based on the question being asked. For example, in testing the impact of drought, a year before the beginning of the drought, 1990, was selected as an initial year when conditions were considered 'Normal' (Peterson TJ, Saft M, Peel MC & John A (2021), Watersheds may not recover from drought, Science, DOI: \doi{10.1126/science.abd5085})
#'
#'
#' @param model fitted hydroState model object.
#' @param initial.year integer with year (YYYY). Default is first year in input.data.
#'
#' @return
#' A fitted hydroState model object with state names for each time-step ready for \code{plot}
#'
#' @keywords state names
#'
#' @export
#'
#' @examples
#' # Load fitted model
#' data(model.annual.fitted.221201)
#'
#' ## Set initial year to set state names
#' model.annual.fitted.221201 =
#'                 setInitialYear(model = model.annual.fitted.221201,
#'                                initial.year = 1990)
#'
#'


setInitialYear <- function(model, initial.year){ #make go to first year of dataframe

  #set state names
  if(is.null(initial.year)){
    stop("Please provide an initial.year to set the state names")
  }
  return(setStateNames(model, initial.year))


}


#' Plot states or pseudo residuals over time
#'
#' \code{plot}
#'
#' @description
#' \code{plot} produces several figures to visualize pseudo residuals or results of the markov states over time. \code{setInitialYear} is required before \code{plot}. It is recommend to evaluate the pseudo residuals before the markov states. The pseudo residuals are the probability of an observation occurring at each time-step given the prior observations and latter observations, and these are derived from the conditional probabilities of the observations. The markov states are from the Viterbi algorithm globally decoding the model to estimate the most probable sequence of states.
#'
#' @details
#' \code{plot} produces five figures of psuedo residuals OR up to four figures of the results from the fitted hydroState model. When the \code{pse.residuals} is FALSE, the default \code{plot} produces all four result figures. Figures are more easily viewed as a pdf exported to the current working directory (\code{do.pdf = TRUE}).
#' \itemize{
#'  \item{psuedo residual figures}{}
#'    \itemize{
#'     \item{A)}{ Time-series of normal-pseudo residuals to ensure the residuals each year are within the confidence intervals.}
#'     \item{B)}{ Auto-correlation function (ACF) of normal-pseudo residuals to ensure there is minimal serial correlation in residuals. Lag spikes should be below confidence interval at each lag (except 0).}
#'     \item{C)}{ Histogram of uniform-pseudo residuals should show uniform distribution (equal frequency for each residual value)}
#'     \item{D)}{ Histogram of normal-pseudo residuals should show normal distribution centered on zero and with no skew}
#'     \item{E)}{ Quantile-Quantile (Q-Q) plot where normal-pseudo residuals vs. theoretical quantities should align on the diagonal line. The last plot contains the Akaike information criterion (AIC) and Shapiro-Wilk p-value. The AIC is an estimator to determine the most parsimonious, best performing model given the number of parameters. When comparing models, the lowest AIC is the best performing model. Shapiro-Wilks test for normality in the residuals and a p-value greater than 0.05 (chosen alpha level) indicates the residuals are normally distributed; the null hypothesis that the residuals are normally distributed is not rejected.}
#'  }
#'  \item{markov state figures}{}
#'    \itemize{
#'      \item{A)}{ independent variable: precipitation}
#'      \item{B)}{ dependent variable and states: streamflow observations, most likely state, and relative normal state estimate}
#'      \item{C)}{ transformed dependent variable and states: transformed streamflow observations and most likely state}
#'      \item{D)}{ conditional state probabilities for each state: probability of hydroState model remaining in given state}
#'    }
#'   }
#'
#'   These figures are often large, and below are a few common errors when the plotting window is too small. Exporting the plots as a pdf is recommend for the pseudo residual figure (\code{do.pdf = TRUE}).
#'  \itemize{
#'  \item{"Error in plot.new() : figure margins too large":} reset plot window with "dev.off()", enlarge plot area and re-run \code{plot.residuals}.
#'  \item{"Error in par(op) : invalid value specified for graphical parameter "pin"} if the R plot window is not reset with "dev.off", an additional \code{plot.residuals} attempt will result in this error.
#'  }
#'
#' @param x is the fitted hydroState model object.
#' @param pse.residuals option to plot pseudo residuals. Default is FALSE.
#' @param ind.variable option to plot independent variable over time. Default is TRUE.
#' @param dep.variable option to plot dependent variable and states over time. Default is TRUE.
#' @param dep.variable.transformed option to plot transformed dependent variable and states over time. Default is TRUE.
#' @param cond.state.prob option to plot the conditional state probabilities over time for each state. Default is TRUE.
#' @param siteID character string of catchment identifier (i.e. gauge ID). Default is NULL. Only recommended when do.pdf = TRUE.
#' @param do.pdf option to export figures as a pdf. Default is FALSE.
#' @param ... additional arguments passed for plotting, none available at this time.
#'
#' @return
#' plots to evaluate rainfall-runoff states over time along with observations and the conditional probabilities of each state.
#'
#' @keywords plot states results
#'
#'
#' @export
#' @import diagram
#' @import graphics
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom utils tail
#' @importFrom zoo na.approx
#'
#' @examples
#' # Load fitted model
#' data(model.annual.fitted.221201)
#'
#' ## Set initial year to set state names
#' model.annual.fitted.221201 =
#'                   setInitialYear(model = model.annual.fitted.221201,
#'                   initial.year = 1990)
#'
#' ## Plot only residuals
#' plot(model.annual.fitted.221201, pse.residuals = TRUE)
#'
#' ## Plot all markov state figures
#' plot(model.annual.fitted.221201)
#'
#' ## Plot only dependent variable transformed with markov states
#' plot(model.annual.fitted.221201,
#'              ind.variable = FALSE,
#'              dep.variable = FALSE,
#'              dep.variable.transformed = TRUE,
#'              cond.state.prob = FALSE)
#'


plot.hydroState <- function(x, ...,
                       pse.residuals = FALSE,
                       ind.variable = TRUE,
                       dep.variable = TRUE,
                       dep.variable.transformed = TRUE,
                       cond.state.prob = TRUE,
                       siteID = NULL,
                       do.pdf = FALSE
                       ){

  model = x

    if(pse.residuals == TRUE & do.pdf == TRUE){

      pdf(paste(siteID,"_Residuals_",class(model@QhatModel.object)[1],"_states-",model@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
      check.PseudoResiduals(model, do.plot = T)
      title(siteID)
      dev.off()
      return(check.PseudoResiduals(model, do.plot = T))

    }else if(pse.residuals == TRUE & do.pdf == FALSE){

      par(mar=c(4,4,1,1))
      return(check.PseudoResiduals(model, do.plot = T))
    }


    # plot everything, default, export as pdf is option too..
    if(ind.variable == TRUE && dep.variable == TRUE && dep.variable.transformed == TRUE && cond.state.prob == TRUE){

      if(do.pdf == TRUE){

        pdf(paste(siteID,"_Results_ABCD_",class(model@QhatModel.object)[1],"_states-",model@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
        viterbi(model, do.plot = T, plot.options = c("A","B","C","D"))
        title(siteID)
        dev.off()
        return(viterbi(model, do.plot = T, plot.options = c("A","B","C","D")))

      } else {

          return(viterbi(model, do.plot = T, plot.options = c("A","B","C","D")))

      }

    }

    # plot only A
    if(ind.variable == TRUE && dep.variable != TRUE && dep.variable.transformed != TRUE && cond.state.prob != TRUE){

      if(do.pdf == TRUE){

        pdf(paste(siteID,"_Results_A_",class(model@QhatModel.object)[1],"_states-",model@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
        viterbi(model, do.plot = T, plot.options = c("A"))
        title(siteID)
        dev.off()
        return(viterbi(model, do.plot = T, plot.options = c("A")))

      } else {

            return(viterbi(model, do.plot = T, plot.options = c("A")))

        }

      }



    # plot only A and B
    if(ind.variable == TRUE && dep.variable == TRUE && dep.variable.transformed != TRUE && cond.state.prob != TRUE){
#
      if(do.pdf == TRUE){

        pdf(paste(siteID,"_Results_AB_",class(model@QhatModel.object)[1],"_states-",model@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
        viterbi(model, do.plot = T, plot.options = c("A","B"))
        title(siteID)
        dev.off()
        return(viterbi(model, do.plot = T, plot.options = c("A","B")))

      } else {
        return(viterbi(model, do.plot = T, plot.options = c("A","B")))
      }
    }

      # plot only A and B and C
      if(ind.variable == TRUE && dep.variable == TRUE && dep.variable.transformed == TRUE && cond.state.prob != TRUE){

        if(do.pdf == TRUE){

          pdf(paste(siteID,"_Results_ABC_",class(model@QhatModel.object)[1],"_states-",model@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
          viterbi(model, do.plot = T, plot.options = c("A","B","C"))
          title(siteID)
          dev.off()
          return(viterbi(model, do.plot = T, plot.options = c("A","B","C")))

        } else {
          return(viterbi(model, do.plot = T, plot.options = c("A","B","C")))
        }
      }

    # plot only A and D
    if(ind.variable == TRUE && dep.variable != TRUE && dep.variable.transformed != TRUE && cond.state.prob == TRUE){

      if(do.pdf == TRUE){

        pdf(paste(siteID,"_Results_AD_",class(model@QhatModel.object)[1],"_states-",model@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
        viterbi(model, do.plot = T, plot.options = c("A","D"))
        title(siteID)
        dev.off()
        return(viterbi(model, do.plot = T, plot.options = c("A","D")))

      } else {
        return(viterbi(model, do.plot = T, plot.options = c("A","D")))
      }
    }

    # plot only A and C
    if(ind.variable == TRUE && dep.variable != TRUE && dep.variable.transformed == TRUE && cond.state.prob != TRUE){

      if(do.pdf == TRUE){

        pdf(paste(siteID,"_Results_AC_",class(model@QhatModel.object)[1],"_states-",model@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
        viterbi(model, do.plot = T, plot.options = c("A","C"))
        title(siteID)
        dev.off()
        return(viterbi(model, do.plot = T, plot.options = c("A","C")))

      } else {
        return(viterbi(model, do.plot = T, plot.options = c("A","C")))
      }
    }

    # plot only A , C, D
    if(ind.variable == TRUE && dep.variable != TRUE && dep.variable.transformed == TRUE && cond.state.prob == TRUE){

      if(do.pdf == TRUE){

        pdf(paste(siteID,"_Results_ACD_",class(model@QhatModel.object)[1],"_states-",model@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
        viterbi(model, do.plot = T, plot.options = c("A","C","D"))
        title(siteID)
        dev.off()
        return(viterbi(model, do.plot = T, plot.options = c("A","C","D")))

      } else {
        return(viterbi(model, do.plot = T, plot.options = c("A","C","D")))
      }
    }

    # plot only B
    if(ind.variable != TRUE && dep.variable == TRUE && dep.variable.transformed != TRUE && cond.state.prob != TRUE){

      if(do.pdf == TRUE){

        pdf(paste(siteID,"_Results_B_",class(model@QhatModel.object)[1],"_states-",model@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
        viterbi(model, do.plot = T, plot.options = c("B"))
        title(siteID)
        dev.off()
        return(viterbi(model, do.plot = T, plot.options = c("B")))

      } else {
        return(viterbi(model, do.plot = T, plot.options = c("B")))
      }
    }

    # plot only B and C
    if(ind.variable != TRUE && dep.variable == TRUE && dep.variable.transformed == TRUE && cond.state.prob != TRUE){

      if(do.pdf == TRUE){

        pdf(paste(siteID,"_Results_BC_",class(model@QhatModel.object)[1],"_states-",model@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
        viterbi(model, do.plot = T, plot.options = c("B","C"))
        title(siteID)
        dev.off()
        return(viterbi(model, do.plot = T, plot.options = c("B","C")))

      } else {
        return(viterbi(model, do.plot = T, plot.options = c("B","C")))
      }
    }

    # plot only B and D
    if(ind.variable != TRUE && dep.variable == TRUE && dep.variable.transformed != TRUE && cond.state.prob == TRUE){
#
      if(do.pdf == TRUE){

        pdf(paste(siteID,"_Results_BD_",class(model@QhatModel.object)[1],"_states-",model@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
        viterbi(model, do.plot = T, plot.options = c("B","D"))
        title(siteID)
        dev.off()
        return(viterbi(model, do.plot = T, plot.options = c("B","D")))

      } else {
        return(viterbi(model, do.plot = T, plot.options = c("B","D")))
      }
    }

    # plot only B,  C and D
    if(ind.variable != TRUE && dep.variable == TRUE && dep.variable.transformed == TRUE && cond.state.prob == TRUE){

      if(do.pdf == TRUE){

        pdf(paste(siteID,"_Results_BCD_",class(model@QhatModel.object)[1],"_states-",model@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
        viterbi(model, do.plot = T, plot.options = c("B","C","D"))
        title(siteID)
        dev.off()
        return(viterbi(model, do.plot = T, plot.options = c("B","C","D")))

      } else {
        return(viterbi(model, do.plot = T, plot.options = c("B","C","D")))
      }
    }

    # plot only  C
    if(ind.variable != TRUE && dep.variable != TRUE && dep.variable.transformed == TRUE && cond.state.prob != TRUE){

      if(do.pdf == TRUE){

          pdf(paste(siteID,"_Results_C_",class(model@QhatModel.object)[1],"_states-",model@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
        viterbi(model, do.plot = T, plot.options = c("C"))
        title(siteID)
        dev.off()
        return(viterbi(model, do.plot = T, plot.options = c("C")))

      } else {
        return(viterbi(model, do.plot = T, plot.options = c("C")))
      }
    }

    # plot only  C & D
    if(ind.variable != TRUE && dep.variable != TRUE && dep.variable.transformed == TRUE && cond.state.prob == TRUE){

      if(do.pdf == TRUE){

        pdf(paste(siteID,"_Results_CD_",class(model@QhatModel.object)[1],"_states-",model@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
        viterbi(model, do.plot = T, plot.options = c("C","D"))
        title(siteID)
        dev.off()
        return(viterbi(model, do.plot = T, plot.options = c("C","D")))

      } else {
        return(viterbi(model, do.plot = T, plot.options = c("C","D")))
      }
    }

    # plot only   D
    if(ind.variable != TRUE && dep.variable != TRUE && dep.variable.transformed != TRUE && cond.state.prob == TRUE){

      if(do.pdf == TRUE){

        pdf(paste(siteID,"_Results_D_",class(model@QhatModel.object)[1],"_states-",model@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
        viterbi(model, do.plot = T, plot.options = c("D"))
        title(siteID)
        dev.off()
        return(viterbi(model, do.plot = T, plot.options = c("D")))

      } else {
        return(viterbi(model, do.plot = T, plot.options = c("D")))
      }
    }


}

#' Get states
#'
#' \code{get.states}
#'
#' @description
#' \code{get.states} uses the Viterbi algorithm to globally decode the model and estimate the most probable sequence of states.
#'
#' @details These dataframe of results include:
#' \itemize{
#' \item{time-step:} year and possibly either season or month for subannual analysis
#' \item{Viterbi State Number:} state number (i.e. 1, 2, or 3) to differentiate states
#' \item{Obs. flow:} streamflow observations
#' \item{Viterbi Flow:} flow values of the Viterbi state including the 5\% and 95\% confidence intervals. These are the most likely flow state values at each time-step of the given states.
#' \item{Normal State Flow:} flow values of the normal state including the 5\% and 95\% confidence intervals. These Normal state flow values are the values from the normal state at each time-step. When the most likely state is the Normal state for a time-step, the Viterbi flow state value equals the Normal flow state value. This Normal state can be visualized relative to the most likely Viterbi state in the "dep.variable" plot from \code{plot.states}.
#' \item{Conditional Prob:} conditional probabilities for each state show the probability of remaining in the given state. When the conditional probability is closer to 1, there is a higher probability that hydroState model remains in that state for the next time-step.
#' \item{Emission Density:} emission density for each state is the result of multiplying the conditional probabilities by the transition probabilities at each timestep.
#' }
#'
#'
#' @param model fitted \code{hydroState} model object.
#'
#' @return
#' data frame of results to evaluate the rainfall-runoff states over time
#'
#' @keywords get states results
#'
#'
#' @export get.states
#'
#' @examples
#' # Load fitted model
#' data(model.annual.fitted.221201)
#'
#' ## Set initial year to set state names
#' model.annual.fitted.221201 =
#'                 setInitialYear(model = model.annual.fitted.221201,
#'                                initial.year = 1990)
#'
#' ## Get states
#' model.annual.fitted.221201.states =
#'                 get.states(model = model.annual.fitted.221201)
#'


get.states <- function(model){

      return(viterbi(model, do.plot = F))

}

#' Check reliability of state predictions
#'
#' \code{check}
#'
#' @description
#' \code{check} After fitting the model, the reliability of the estimated states can be assessed by generating synthetic state sequences and then assessing how well the model identified them.
#'
#' @details
#' This validates the model's states at each time-step through re-sampling the input data and re-running the Viterbi algorithm. The input data is duplicated 100 times, and a synthetic series is generated from the model with sample states. This provides a time series of the transformed streamflow observations that can be compared with observations of the `known' state. The Viterbi states of the re-sampled transformed observations are inferred, and the probability of the inferred state equaling the 'known' state is calculated.
#'
#' @param model fitted \code{hydroState} model.
#' @param n.samples integer of samples to re-sample. Default is 100000.
#'
#' @return
#' A data frame is returned with a matrix depending on the number of states. For a 2 state model, a 2x2 matrix is returned. The diagonal cell estimates the probability of correctly identifying that state. The off diagonals estimate the probability of incorrectly identifying a state that in state 2.
#'
#' @keywords check model viterbi
#'
#'
#' @export check
#'
#' @examples
#' ## Check reliability of state predictions (>5s to run)
#' \dontrun{
#' check(model = model.annual.fitted.221201)
#'}


check <- function(model, n.samples = 100000){

  if(class(model)[1] == "hydroState"){

  return(check.viterbi(model, nSamples = n.samples))
  }else{
    stop('model is not a hydroState object')
  }

}



