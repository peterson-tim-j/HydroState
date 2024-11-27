#' Set Seasons
#'
#' \code{set.seasons}
#'
#' @description
#' Aggregates monthly data to 4 seasons in a year.
#'
#' @details
#' Sets 4 seasons
#'
#' @param input.data dataframe of monthly runoff and precipitation observations. Gaps with missing data in either streamflow or precipitation are permitted, and the handling of them is further discussed in \code{select.Markov}. Monthly data is required when using \code{seasonal.parameters} that assumes selected model parameters are better defined with a sinusoidal function.
#'
#' @return
#' A dataframe of seasonal observations with an additional column counting the number of months in each season.
#'
#' @keywords seasons
#'
#'
#' @export set.seasons
#'
#' @examples
#'
#' # Load data
#' data(streamflow_monthly_221201)
#'
#' # aggregate monthly data to seasonal
#' streamflow_seasonal_221201 = set.seasons(streamflow_monthly_221201)
#'
#'

set.seasons <- function(input.data=data.frame(year=c(), month = c(), flow=c(), precip=c())){

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


select.transform <- function(func = 'boxcox', input.data=data.frame(year=c(), flow=c(), precip=c())){

  #Validate
  func = paste('Qhat.',func,sep='')
  if(func %in% c('Qhat.none','Qhat.log','Qhat.burbidge','Qhat.boxcox')){

    # If monthly data, sort in ascending order by year and month
    if('month' %in% colnames(input.data)){
      input.data = input.data[order(input.data[,'year'],input.data[,'month']),]
    }else if('day' %in% colnames(input.data)){
      input.data = input.data[order(input.data[,'year'],input.data[,'month'],input.data[,'day']),]
    }


    return(new(func, input.data))

  }else{

    stop('Please choose an appropriate function: none, log, burbidge, or boxcox')
  }


}

select.stateModel <- function(input.data = data.frame(year=c(), flow=c(), precip=c()),
                              parameters = list('a0','a1','std'),
                              seasonal.parameters = list(),
                              state.shift.parameters = list('a0','std'),
                              error.distribution = 'truc.normal',
                              transition.graph = matrix(TRUE,2,2)){

  #Validate input.data
  if(!is.data.frame(input.data)){
    stop("'input.data' is not a dataframe")
  }
  if(!('year' %in% colnames(input.data))){
    stop("'input.data' must contain a 'year' column with an integer of years")
  }

  #Validate parameters
  if(!is.list(parameters)){
    stop("'parameters' are not in a list")
  }
  if(!is.character(unlist(parameters))){
    stop("'parameters' are not characters")
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
  if(length(seasonal.parameters) > 0 && !is.character(unlist(seasonal.parameters))){
    stop("'seasonal.parameters' are not characters")
  }
  if(!is.list(seasonal.parameters)){
    stop("'seasonal.parameters' are not in a list")
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
  if(!is.character(unlist(state.shift.parameters))){
    stop("'state.shift.parameters' are not characters")
  }
  if(!is.list(state.shift.parameters)){
    stop("'state.shift.parameters' are not in a list")
  }
  # If state shift parameters contain auto-correlation terms, but no auto-correlation in model parameters, STOP
  if('AR1' %in% state.shift.parameters || 'AR2' %in% state.shift.parameters || 'AR3' %in% state.shift.parameters){
    if(!(any(parameters == 'AR1') || any(parameters == 'AR2') ||  any(parameters == 'AR3'))){
      stop("'state.shift.parameters' cannot contain auto-correlation terms because 'parameters' in state.model does not contain auto-correlation")
    }
  }

  #Validate error.distribution
  if(!is.character(error.distribution)){
    stop("'error.distribution' is not a string of characters")
  }
  if(!(error.distribution %in% c('normal','truc.normal','gamma'))){
    stop("'error.distribution' must contain either 'normal', 'truc.normal' or 'gamma'")
  }
  if(length(seasonal.parameters) > 0 && (error.distribution == 'normal' || error.distribution == 'truc.normal')){
    stop("Please select 'gamma' for error.distribution. At the moment, hydroState does not support subannual analysis with 'error.distribution' = 'normal' or 'truc.normal'")
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


select.Markov <- function(flickering = FALSE,
                          transition.graph){

  #Validate
  if(flickering == TRUE){

    Markov = new('markov.annualHomogeneous.flickering', transition.graph)

  }else{

    Markov = new('markov.annualHomogeneous',transition.graph)

  }

}

#'Builds hydroState model
#'
#' \code{buildModel}
#'
#' @description
#' \code{buildModel} builds a hydrostate model with either a default model or the model can be specified with options from below. Every model depends on a linear base model where streamflow, \eqn{Q}, is a function of precipitation, \eqn{P}: \eqn{Q = Pa_1 + a_0}. The default model is an analysis of this base linear model with state shifts expected in the intercept, \eqn{a_0}, and standard deviation, \eqn{std}, of the rainfall-runoff relationship. There are additional adjustments to this model with auto-correlations terms, seasonal parameters for a sub-annual analysis, and even evaluation of other state dependent, shifting, parameters. The number of states and assumed error distribution can also be selected. After the model is built, the hydroState model is ready to be fitted with \code{fitModel}
#'
#' @details
#' There are a selection of items to consider when defining the rainfall-runoff relationship and investigating state shifts in this relationship. hydroState provides various options for modelling the rainfall-runoff relationship.
#' \itemize{
#' \item{Data gaps  with \code{input.data}}: When there is missing \code{input.data}, special care was taken to reduce the influence of the missing time periods while making the most of the given data without infilling. For time-periods where the dependent variable, streamflow, is missing, the transition probability for these missing periods is essentially ignored by setting the conditional probability of the missing time-steps equal to one. This results in the time-step after the missing period having the same probability of being in the given state as before the missing period. For time-periods where the independent variable is missing, precipitation, the missing periods are essentially ignored when fitting the models as the likelihood is calculated from only continuous periods with values for the independent variable. Since the log-likelihood is calculated for each continuous period, the sum of the log-likelihoods provides a total log-likelihood to fit models. Again, there is no infilling of missing data.
#' \item{Transform Observations with \code{data.transform}}: Transforms observations to remove heteroscedasticity. Often there is skew within hydrologic data. When defining relationships between observations, this skew results in an unequal variance in the residuals, heteroscedasticity. Transforming observations is often required with observations of streamflow and precipitation. There are several options to transform observations. Since the degree of transformation is not typically known, 'boxcox' is the default. Other options include: 'log', 'burbidge', and of course, 'none' when no transformation is performed.
#' \item{Model Structure with \code{parameters} and \code{seasonal.parameters}}: The structure of the model depends on the \code{parameters}. hydroState simulates runoff, \eqn{Q}, as being in one of finite states, \eqn{i}, at every time-step, \eqn{t}, depending on the distribution of states at prior time steps. This results in a runoff distribution for each state that can vary overtime (\eqn{\widehat{_tQ_i}}). The model defines the relationship that is susceptible to state shifts with precipitation, \eqn{P_t}, as a predictor. This takes the form as a simple linear model \eqn{\widehat{_tQ_i} = f(P_t)}:
#'
#' \eqn{\widehat{_tQ_i} = P_ta_1 + a_0}
#'
#' where \eqn{a_0} and \eqn{a_1} are constant parameters. These parameters and the model error, \eqn{std}, are required parameters for every model. It is possible the relationship contains serial correlation and would be better defined with an auto-regressive term:
#'
#' \eqn{\widehat{_tQ_i} = P_ta_1 + a_0 + AR1\widehat{_{t-1}Q}}
#'
#' where \eqn{AR1} is the lag-1 auto-correlation term. Either, lag-1: \eqn{AR1}, lag-2: \eqn{AR2}, and lag-3: \eqn{AR3} auto-correlation coefficients are an option as additional parameters to better define the rainfall-runoff relationship.
#' For sub-annual analysis, \code{seasonal.parameters} provides the option to assume a sinusoidal function better defines either of the constant parameters or error (\eqn{a_0, a_1, std}) throughout the year, i.e:
#'
#' \eqn{a_0 = a_{0.disp} + a_{0.amp} * sin(2\pi(\frac{M_t}{12} + a_{0.phase}))}
#'
#' where \eqn{M_t} is an integer month at \eqn{t}. Monthly streamflow and precipitation are required as \code{input.data} for the sub-annual analysis.
#'
#' \item{State Dependent Parameters with \code{state.shift.parameters}}: These are state dependent parameters where they are subject to shift in order to better explain the state of streamflow over time. Any or all of the previously chosen parameters can be selected (\eqn{a_0, a_1, std, AR1, AR2, AR3}). The default model evaluates shifts in the rainfall-runoff relationship with \eqn{a_0} \eqn{std} as state dependent parameters.
#' \item{Distribution of the Residuals with \code{error.distribution}}: The distribution of the residuals (error) in the model can be chosen to reduce skew and assist with making models statistically adequate (see \code{plot.residuals}). Either normal: 'normal', truncated normal: 'truc.normal', or gamma: 'gamma' distributions are acceptable. The default is 'truc.normal'. Sub-annual models are restricted to only a 'gamma' distribution.
#' \item{Markov flickering with \code{flickering}}: The form of the Markov model is homogeneous where the transition probabilities are time-invariant. The model be performed with or without flickering between states. Flickering within the Markov model provides a more sensitive analysis. The default is FALSE, no flickering.
#' \item{Number of States with \code{transition.graph}}: The number of possible states in the rainfall-runoff relationship and transition between the states is selected with the transition.graph. The default is a 2-state model in a 2 by 2 matrix with a TRUE transition to and from each state.
#'}
#'
#' @param input.data dataframe of annual, seasonal, or monthly runoff and precipitation observations. Gaps with missing data in either streamflow or precipitation are permitted. Monthly data is required when using \code{seasonal.parameters}.
#' @param data.transform character string with the method of transformation. The default is 'boxcox'. Other options: 'log', 'burbidge', 'none'
#' @param parameters character list of parameters to construct model. Required and default: \code{a0}, \code{a1}, \code{std}. Auto-correlation terms optional: \code{AR1}, \code{AR2}, or \code{AR3}.
#' @param seasonal.parameters character list of one or all parameters (\code{a0}, \code{a1}, \code{std}) defined as a sinusoidal function to represent seasonal variation. Requires monthly or seasonal data. Default is empty list.
#' @param state.shift.parameters character list of one or all parameters (\code{a0}, \code{a1}, \code{std}, \code{AR1}, \code{AR2}, \code{AR3}) able to shift as dependent on state. Default is \code{a0} and \code{std}.
#' @param error.distribution character sting of the distribution in the HMM error. Default is 'truc.normal'. Others include: 'normal' or 'gamma'
#' @param flickering logical T/F. T = allows more sensitive markov flickering between states over time, F = less sensitive and is default.
#' @param transition.graph matrix given the number of states. Default is a 2-state matrix (2 by 2): matrix(TRUE,2,2)
#'
#' @return
#' A built hydroState model object ready to be fitted with \code{fitModel}
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
#' model = buildModel(input.data = streamflow_annual_221201)
#'
#' # OR
#'
#' ## Build annual hydroState model with specified objects

#' # Build hydroState model with: 2-state, normal error distribution,
#' # 1-lag of auto-correlation, and state dependent parameters ('a1', 'std')
#' model = buildModel(input.data = streamflow_annual_221201,
#'                    data.transform = 'boxcox',
#'                    parameters = list('a0','a1','std','AR1'),
#'                    seasonal.parameters = list(),
#'                    state.shift.parameters = list('a1','std'),
#'                    error.distribution = 'normal',
#'                    flickering = FALSE,
#'                    transition.graph = matrix(TRUE,2,2))
#'



buildModel <- function(input.data = data.frame(year=c(), flow=c(), precip=c()),
                       data.transform = 'boxcox',
                       parameters = list('a0','a1','std'),
                       seasonal.parameters = list(),
                       state.shift.parameters = list('a0','std'),
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
  #######################################################

# If no inputs except input.data, just run default analysis
  if(is.null(parameters) && is.null(seasonal.parameters) && is.null(state.shift.parameters)){

    # default transition graph if NULL
    if(is.null(transition.graph))
      transition.graph = matrix(TRUE,2,2)

    # default flickering to FALSE if NULL
    if(is.null(flickering))
      flickering = FALSE

    # default data.transform if NULL
    if(is.null(data.transform)){
      data.transform = select.transform(func = 'boxcox', input.data)
    }else{
      data.transform = select.transform(data.transform, input.data)
    }

    # default state model error.distribution if NULL
    if(is.null(error.distribution))
      error.distribution = 'truc.normal'

    # create state model
    stateModel = select.stateModel(input.data,
                                   parameters = list('a0','a1','std'),
                                   seasonal.parameters = list(),
                                   state.shift.parameters = list('a0','std'),
                                   error.distribution,
                                   transition.graph = transition.graph)

    # create Markov model
    Markov = select.Markov(flickering, transition.graph)

    # build default model
    return(new('hydroState',input.data, data.transform, stateModel, Markov))

  }else if(!is.null(parameters) && is.null(seasonal.parameters) && is.null(state.shift.parameters)){

    # default transition graph if NULL
    if(is.null(transition.graph))
      transition.graph = matrix(TRUE,2,2)

    # default flickering to FALSE if NULL
    if(is.null(flickering))
      flickering = FALSE

    # default data.transform if NULL
    if(is.null(data.transform)){
      data.transform = select.transform(func = 'boxcox', input.data)
    }else{
      data.transform = select.transform(data.transform, input.data)
    }

    # default state model error.distribution if NULL
    if(is.null(error.distribution))
      error.distribution = 'truc.normal'

    # create state model
    stateModel = select.stateModel(input.data,
                                   parameters,
                                   seasonal.parameters = list(),
                                   state.shift.parameters = list('a0','std'),
                                   error.distribution,
                                   transition.graph = transition.graph)

    # create Markov model
    Markov = select.Markov(flickering, transition.graph)

    # build default model
    return(new('hydroState',input.data, data.transform, stateModel, Markov))

  }else if(!is.null(parameters) && !is.null(seasonal.parameters) && is.null(state.shift.parameters)){

    # default transition graph if NULL
    if(is.null(transition.graph))
      transition.graph = matrix(TRUE,2,2)

    # default flickering to FALSE if NULL
    if(is.null(flickering))
      flickering = FALSE

    # default data.transform if NULL
    if(is.null(data.transform)){
      data.transform = select.transform(func = 'boxcox', input.data)
    }else{
      data.transform = select.transform(data.transform, input.data)
    }

    # default state model error.distribution if NULL
    if(is.null(error.distribution))
      error.distribution = 'truc.normal'

    # create state model
    stateModel = select.stateModel(input.data,
                                   parameters,
                                   seasonal.parameters,
                                   state.shift.parameters = list('a0','std'),
                                   error.distribution,
                                   transition.graph = transition.graph)

    # create Markov model
    Markov = select.Markov(flickering, transition.graph)

    # build default model
    return(new('hydroState',input.data, data.transform, stateModel, Markov))

  }else if(!is.null(parameters) && !is.null(seasonal.parameters) && !is.null(state.shift.parameters)){

    # default transition graph if NULL
    if(is.null(transition.graph))
      transition.graph = matrix(TRUE,2,2)

    # default flickering to FALSE if NULL
    if(is.null(flickering))
      flickering = FALSE

    # default data.transform if NULL
    if(is.null(data.transform)){
      data.transform = select.transform(func = 'boxcox', input.data)
    }else{
      data.transform = select.transform(data.transform, input.data)
    }

    # default state model error.distribution if NULL
    if(is.null(error.distribution))
      error.distribution = 'truc.normal'

    # create state model
    stateModel = select.stateModel(input.data,
                                   parameters,
                                   seasonal.parameters,
                                   state.shift.parameters,
                                   error.distribution,
                                   transition.graph = transition.graph)

    # create Markov model
    Markov = select.Markov(flickering, transition.graph)

    # build default model
    return(new('hydroState',input.data, data.transform, stateModel, Markov))

  }else if(is.null(parameters) && is.null(seasonal.parameters) && !is.null(state.shift.parameters)){

    # default transition graph if NULL
    if(is.null(transition.graph))
      transition.graph = matrix(TRUE,2,2)

    # default flickering to FALSE if NULL
    if(is.null(flickering))
      flickering = FALSE

    # default data.transform if NULL
    if(is.null(data.transform)){
      data.transform = select.transform(func = 'boxcox', input.data)
    }else{
      data.transform = select.transform(data.transform, input.data)
    }

    # default state model error.distribution if NULL
    if(is.null(error.distribution))
      error.distribution = 'truc.normal'

    # create state model
    stateModel = select.stateModel(input.data,
                                   parameters = list('a0','a1','std'),
                                   seasonal.parameters = list(),
                                   state.shift.parameters,
                                   error.distribution,
                                   transition.graph = transition.graph)

    # create Markov model
    Markov = select.Markov(flickering, transition.graph)

    # build default model
    return(new('hydroState',input.data, data.transform, stateModel, Markov))

  }else if(is.null(parameters) && !is.null(seasonal.parameters) && !is.null(state.shift.parameters)){

    # default transition graph if NULL
    if(is.null(transition.graph))
      transition.graph = matrix(TRUE,2,2)

    # default flickering to FALSE if NULL
    if(is.null(flickering))
      flickering = FALSE

    # default data.transform if NULL
    if(is.null(data.transform)){
      data.transform = select.transform(func = 'boxcox', input.data)
    }else{
      data.transform = select.transform(data.transform, input.data)
    }

    # default state model error.distribution if NULL
    if(is.null(error.distribution))
      error.distribution = 'truc.normal'

    # create state model
    stateModel = select.stateModel(input.data,
                                   parameters = list('a0','a1','std'),
                                   seasonal.parameters,
                                   state.shift.parameters,
                                   error.distribution,
                                   transition.graph = transition.graph)

    # create Markov model
    Markov = select.Markov(flickering, transition.graph)

    # build default model
    return(new('hydroState',input.data, data.transform, stateModel, Markov))

  }else if(!is.null(parameters) && is.null(seasonal.parameters) && !is.null(state.shift.parameters)){

    # default transition graph if NULL
    if(is.null(transition.graph))
      transition.graph = matrix(TRUE,2,2)

    # default flickering to FALSE if NULL
    if(is.null(flickering))
      flickering = FALSE

    # default data.transform if NULL
    if(is.null(data.transform)){
      data.transform = select.transform(func = 'boxcox', input.data)
    }else{
      data.transform = select.transform(data.transform, input.data)
    }

    # default state model error.distribution if NULL
    if(is.null(error.distribution))
      error.distribution = 'truc.normal'

    # create state model
    stateModel = select.stateModel(input.data,
                                   parameters,
                                   seasonal.parameters = list(),
                                   state.shift.parameters,
                                   error.distribution,
                                   transition.graph = transition.graph)

    # create Markov model
    Markov = select.Markov(flickering, transition.graph)

    # build default model
    return(new('hydroState',input.data, data.transform, stateModel, Markov))

  }else if(is.null(parameters) && !is.null(seasonal.parameters) && is.null(state.shift.parameters)){

    # default transition graph if NULL
    if(is.null(transition.graph))
      transition.graph = matrix(TRUE,2,2)

    # default flickering to FALSE if NULL
    if(is.null(flickering))
      flickering = FALSE

    # default data.transform if NULL
    if(is.null(data.transform)){
      data.transform = select.transform(func = 'boxcox', input.data)
    }else{
      data.transform = select.transform(data.transform, input.data)
    }

    # default state model error.distribution if NULL
    if(is.null(error.distribution))
      error.distribution = 'truc.normal'

    # create state model
    stateModel = select.stateModel(input.data,
                                   parameters = list('a0','a1','std'),
                                   seasonal.parameters,
                                   state.shift.parameters = list('a0','std'),
                                   error.distribution,
                                   transition.graph = transition.graph)

    # create Markov model
    Markov = select.Markov(flickering, transition.graph)

    # build default model
    return(new('hydroState',input.data, data.transform, stateModel, Markov))

  }
}


#' Builds all hydroState models
#'
#' \code{buildModelAll}
#'
#' @description
#' \code{buildModelAll} builds all possible combinations of hydroState models
#'
#' @details
#' All possible combinations of hydroState models are built for each data transformations, auto-correlation lag, and residual distribution from 1 to 3 states for investigating only state changes in the 'a0' and 'std' parameters. Note: annual time-step only
#'
#'
#' @param input.data dataframe of annual, seasonal, or monthly runoff and precipitation observations. Gaps with missing data in either streamflow or precipitation are permitted, and the handling of them is further discussed in \code{select.Markov}. Monthly data is required when using \code{seasonal.parameters} that assumes selected model parameters are better defined with a sinusoidal function.
#' @param ID  character vector of a stream gauge identifier
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
#' # Build all annual models
#' all.annual.models = buildModelAll(input.data = streamflow_annual_221201, ID = '221201')
#'


buildModelAll <- function(input.data = data.frame(year=c(), flow=c(), precip=c()),
                          ID = ''){
  #
  if(is.null(ID)){
    ID = ''
  }

  # Validate
  if(is.character(ID)){

    # If monthly data, sort in ascending order by year and month
    if('month' %in% colnames(input.data)){
      input.data = input.data[order(input.data[,'year'],input.data[,'month']),]


      return(new('hydroState.subAnnual.allModels',ID, input.data, allow.flickering=F))

      # remove month column
      # input.data = input.data[,c("year","flow","precipitation")]
      #
      # input.data = aggregate(input.data[c('flow','precipitation')], by=input.data['year'], sum)
      #
      # message('Note: Monthly data inputted. flow and precipitation summed by year. All models built.')

    # }else if('day' %in% colnames(input.data)){
    #   input.data = input.data[order(input.data[,'year'],input.data[,'month'],input.data[,'day']),]
    #
    #   # remove day and month column
    #   input.data = input.data[,c("year","flow","precipitation")]
    #
    #   input.data = aggregate(input.data[c('flow','precipitation')], by=input.data['year'], sum)
    #
    #   message('Note: Daily data inputted. flow and precipitation summed by year. All models built.')

    }



    return(new('hydroState.allModels',ID, input.data, allow.flickering=F))

  }else{

    stop('Please ensure ID is a character vector. See help')

  }

}

#'Fit hydroState model
#'
#' \code{fitModel}
#'
#' @description
#' \code{fitModel} fits hydrostate model(s) using global optimization by differential evolution \href{https://cran.r-project.org/web/packages/DEoptim/index.html}{DEoptim} library.
#'
#' @details
#' After a hydroState model object is built, the model is ready to be fitted through minimizing the negative log-likelihood function. The likelihood is estimated recursively across each time-step, for each continous period of observations, and the sum of the negative log-likelihood is minimized too calibrate the model parameters. The only required input is the given name of the built hydroState model object. \code{fitModel} works for one built model (\code{buildModel}) or all (\code{buildModelAll}). If fitting all models be sure to install and load the \href{https://cran.r-project.org/web/packages/parallelly/index.html}{parallelly} library. Details on the likelihood function is as follows:
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
#'    For more details, refer to pg. 8-17 in the \href{https://www.science.org/doi/10.1126/science.abd5085#supplementary-materials}{Supplementary Materials} of "Watersheds may not recover from drought".}}
#'  \item{\eqn{\Gamma} is the transition matrix}
#'  \item{\eqn{T} is the number of time-steps.}
#'  }
#'
#' @param model.name name of the built hydroState model or name of the list containing all built models
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
#' @import DEoptim
#' @import methods
#' @import stats
#' @import truncnorm
#' @importFrom utils relist
#'
#' @examples
#'
#' # Load data
#' data(streamflow_annual_221201)
#'
#' ## Build default annual hydroState model
#' model = buildModel(input.data = streamflow_annual_221201)
#'
#' ## Fit built model
#' model = fitModel(model)
#'
#' ## Fit all built models
#' \dontrun{
#'
#' # Load data
#' data(streamflow_annual_221201)
#'
#' ## Build all annual models
#' all.annual.models = buildModelAll(input.data = streamflow_annual_221201, ID = '221201')
#'
#' ## Fit all
#' model = fitModel(all.annual.models)
#'
#' }
#'
#'

fitModel <- function(model.name,
                     pop.size.perParameter = 10,
                     max.generations = 500,
                     doParallel = F){

  # Validate
  if(class(model.name)[1] == "hydroState"){

    if(doParallel == T){

      return(fit(model.name,
                 DEstrategy=3,
                 pop.size.perParameter = pop.size.perParameter,
                 max.generations=max.generations,
                 Domains = NA,
                 reltol=1e-8,
                 steptol=50,
                 print.iterations = 25,
                 use.initial.parameters=F,
                 doParallel = T))

    }else{

      return(fit(model.name,
                 DEstrategy=3,
                 pop.size.perParameter = pop.size.perParameter,
                 max.generations=max.generations,
                 Domains = NA,
                 reltol=1e-8,
                 steptol=50,
                 print.iterations = 25,
                 use.initial.parameters=F,
                 doParallel = F))
    }

  }else if(class(model.name)[1] == "hydroState.allModels"){

    if(doParallel == T){

      return(fit(model.name,
                 DEstrategy=3,
                 pop.size.perParameter = pop.size.perParameter,
                 max.generations=10000,
                 Domains = NA,
                 reltol=1e-8,
                 steptol=50,
                 print.iterations = 25,
                 use.initial.parameters=F,
                 doParallel = T))
    }else{

      return(fit(model.name,
                 DEstrategy=3,
                 pop.size.perParameter = pop.size.perParameter,
                 max.generations=10000,
                 Domains = NA,
                 reltol=1e-8,
                 steptol=50,
                 print.iterations = 25,
                 use.initial.parameters=F,
                 doParallel = F))

    }

  }else if(class(model.name)[1] == "hydroState.subAnnual.allModels"){


    return(fit(model.name,
               DEstrategy=3,
               pop.size.perParameter = pop.size.perParameter,
               max.generations=10000,
               Domains = NA,
               reltol=1e-8,
               steptol=50,
               print.iterations = 25,
               use.initial.parameters=F,
               doParallel = T))

  }else{


    stop('model.name is not an appropriate class. Please ensure input model is a built hydroState model with the class as either of the following: "hydroState", "hydroState.allModels", or "hydroState.subAnnual.allModels"')

  }

}


#'get AIC
#'
#' \code{get.AIC}
#'
#' @description
#' \code{get.AIC} retrieves Akaike information criteria from a fitted model.
#'
#' @details
#' The AIC is the negative log-likelihood of the model plus a penealty for model parameters. This function can be performed on a single model or a selection of models to find the lowest AIC of the set.
#'
#' @param model.name is the name of the fitted \code{hydroState} model object.
#'
#' @return
#' AIC value
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
#' get.AIC()
#'


get.AIC <- function(model.name){


  if(!(class(model.name)[1] %in% c("hydroState","hydroState.allModels", "hydroState.subAnnual.allModels")))
    stop('model.name is not an appropriate class. Please ensure input model is a built hydroState model with the class as either of the following: "hydroState", "hydroState.allModels", or "hydroState.subAnnual.allModels"')


  # Validate
  if(length(model.name@models) >1){

      return(getAIC.bestModel(model.name)$AIC)

  }else{


      return(getAIC(model.name))

  }

}



#'Plot residuals
#'
#' \code{plot.residuals}
#'
#' @description
#' The normal pseudo residuals are plotted for review to check for outliers and validate the fit of the model. It is recommended to ensure the model fit is valid before evaluating results (i.e. \code{plot.states}). Furthermore, to ensure the multi-state model performs better than the one-state model, it is recommended to visually compare \code{plot.resdiuals} of both models.
#'
#' @details
#' \code{plot.residuals} produces five plots to review and validate the fitted hydroState model.
#' \itemize{
#'  \item{A)}{ Time-series of normal-pseudo residuals to ensure the residuals each year are within the confidence intervals.}
#'  \item{B)}{ Auto-correlation function (ACF) of normal-pseudo residuals to ensure there is no serial correlation in residuals. Lag spikes should be below confidence interval at each lag (except 0).}
#'  \item{C)}{ Histogram of uniform-pseudo residuals should show uniform distribution (equal frequency for each residual value)}
#'  \item{D)}{ Histogram of normal-pseudo residuals should show normal distribution centered on zero and with no skew}
#'  \item{E)}{ Quantile-Quantile (Q-Q) plot where normal-pseudo residuals vs. theoretical quantities should align on the diagonal line. The last plot contains the Akaike information criterion (AIC) and Shapiro-Wilk p-value. The AIC is an estimator to determine the most parsimonious, best performing model given the number of parameters. When comparing models, the lowest AIC is the best performing model. Shapiro-Wilks test for normality in the residuals and a p-value greater than 0.05 (chosen alpha level) indicates the residuals are normally distributed; the null hypothesis that the residuals are normally distributed is not rejected.}
#'  }
#'  It is recommended to export the residual plot as a PDF due to it's size. If the R plot windor is too small, two common errors can occur:
#'  \itemize{
#'  \item{"Error in plot.new() : figure margins too large":} reset plot window with "dev.off()", enlarge plot area and re-run \code{plot.residuals}.
#'  \item{"Error in par(op) : invalid value specified for graphical parameter "pin"} if the R plot window is not reset with "dev.off", an additional \code{plot.residuals} attempt will result in this error.
#'  }
#'
#'
#' @param model.name name of the fitted hydroState model object.
#' @param do.pdf option to export residual plots as a pdf. Default is FALSE.
#' @param ID character string of catchment identifier (i.e. gauge ID). Default is NULL. Only recommended when do.pdf = TRUE.
#'
#'
#' @return
#' Plots of residuals to evaluate model fit
#'
#' @keywords plot residuals
#'
#'
#' @export plot.residuals
#' @importFrom stats acf
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @import diagram
#' @import graphics
#'
#' @examples
#' # Load fitted model
#' data(model.annual.fitted.221201)
#'
#' ## Plot residuals
#' plot.residuals(model.name = model.annual.fitted.221201)
#'


plot.residuals <- function(model.name,
                           do.pdf = FALSE,
                           ID = NULL
                             ){

  if(do.pdf != TRUE){
    par(mar=c(4,4,1,1))
    return(check.PseudoResiduals(model.name, do.plot = T))
  }

  if(do.pdf == TRUE){

    pdf(paste(ID,"_Residuals_",class(model.name@QhatModel.object)[1],"_states-",model.name@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
    check.PseudoResiduals(model.name, do.plot = T)
    title(ID)
    dev.off()
    return(check.PseudoResiduals(model.name, do.plot = T))

  }

  if(do.plot != TRUE && do.pdf == TRUE){
    pdf(paste(ID,"_Residuals_",class(model.name@QhatModel.object)[1],"_states-",model.name@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
    check.PseudoResiduals(model.name, do.plot = T)
    title(ID)
    dev.off()
    return(check.PseudoResiduals(model.name, do.plot = F))
  }

  if(do.plot != TRUE && do.pdf != TRUE){
    return(check.PseudoResiduals(model.name, do.plot = F))
  }




}

#'Get residuals
#'
#' \code{get.residuals}
#'
#' @description
#' The normal pseudo residuals are retrieved from the fitted model.
#'
#' @details
#' \code{get.residuals} retrieves residuals from the fitted model and exports them as a data frame.
#'
#'
#' @param model.name name of the fitted hydroState model object.
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
#' get.residuals(model.name = model.annual.fitted.221201)
#'
#'


get.residuals <- function(model.name){

    return(check.PseudoResiduals(model.name, do.plot = F))

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
#' @param model.name name of the fitted hydroState model object.
#' @param initial.year integer with year (YYYY). Default is first year in input.data.
#'
#' @return
#' A fitted hydroState model object with state names for each time-step ready for \code{plot.states}
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
#'                 setInitialYear(model.name = model.annual.fitted.221201,
#'                                initial.year = 1990)
#'
#'


setInitialYear <- function(model.name, initial.year){ #make go to first year of dataframe

  #set state names
  return(setStateNames(model.name, initial.year))


}


#'plot States
#'
#' \code{plot.states}
#'
#' @description
#' \code{plot.states} produces several plots to visualize results of the states overtime. \code{setInitialYear} is required before \code{plot.states}.
#'
#' @details
#' \code{plot.states} produces four figures of the results from the fitted hydroState model. The default produces all four:
#' \itemize{
#'  \item{independent variable}: precipitation
#'  \item{dependent variable and states}: streamflow observations, most likely state, and relative normal state estimate
#'  \item{transformed dependent variable and states}: transformed streamflow observations and most likely state
#'  \item{conditional state probabilities for each state}: probability of hydroState model remaining in given state
#'  }
#'  These are plotted on the same page, and there is an option to export plots as a pdf to the current working directory. There are also options to only plot one of the four plots.
#'
#'
#' @param model.name is the name of the fitted hydroState model object.
#' @param ind.variable option to plot independent variable overtime. Default is TRUE.
#' @param dep.variable option to plot dependent variable and states overtime. Default is TRUE.
#' @param dep.variable.transformed option to plot transformed dependent variable and states overtime. Default is TRUE.
#' @param cond.state.prob option to plot the conditional state probabilities overtime for each state. Default is TRUE.
#' @param do.pdf option to export plots as a pdf. Default is FALSE.
#' @param ID character string of catchment identifier (i.e. gauge ID). Default is NULL. Only recommended when do.pdf = TRUE.
#'
#'
#' @return
#' plots to evaluate rainfall-runoff states overtime along with observations and the conditional probabilities of each state.
#'
#' @keywords plot states results
#'
#'
#' @export plot.states
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @import diagram
#' @import graphics
#' @importFrom utils tail
#' @importFrom zoo na.approx
#'
#' @examples
#' # Load fitted model
#' data(model.annual.fitted.221201)
#'
#' ## Set initial year to set state names
#' model.annual.fitted.221201 =
#'                   setInitialYear(model.name = model.annual.fitted.221201,
#'                   initial.year = 1990)
#'
#' ## Plot all figures
#' plot.states(model.name = model.annual.fitted.221201)
#'
#' ## Plot only dependent variable transformed with markov states
#' plot.states(model.name = model.annual.fitted.221201,
#'              ind.variable = FALSE,
#'              dep.variable = FALSE,
#'              dep.variable.transformed = TRUE,
#'              cond.state.prob = FALSE)
#'


plot.states <- function(model.name,
                        ind.variable = TRUE,
                        dep.variable = TRUE,
                        dep.variable.transformed = TRUE,
                        cond.state.prob = TRUE,
                        do.pdf = FALSE,
                        ID = NULL
                       ){

    # plot everything, default, export as pdf is option too..
    if(ind.variable == TRUE && dep.variable == TRUE && dep.variable.transformed == TRUE && cond.state.prob == TRUE){

      if(do.pdf == TRUE){

        pdf(paste(ID,"_Results_ABCD_",class(model.name@QhatModel.object)[1],"_states-",model.name@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
        viterbi(model.name, do.plot = T, plot.options = c("A","B","C","D"))
        title(ID)
        dev.off()
        return(viterbi(model.name, do.plot = T, plot.options = c("A","B","C","D")))

      } else {

          return(viterbi(model.name, do.plot = T, plot.options = c("A","B","C","D")))

      }

    }

    # plot only A
    if(ind.variable == TRUE && dep.variable != TRUE && dep.variable.transformed != TRUE && cond.state.prob != TRUE){

      if(do.pdf == TRUE){

        pdf(paste(ID,"_Results_A_",class(model.name@QhatModel.object)[1],"_states-",model.name@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
        viterbi(model.name, do.plot = T, plot.options = c("A"))
        title(ID)
        dev.off()
        return(viterbi(model.name, do.plot = T, plot.options = c("A")))

      } else {

            return(viterbi(model.name, do.plot = T, plot.options = c("A")))

        }

      }



    # plot only A and B
    if(ind.variable == TRUE && dep.variable == TRUE && dep.variable.transformed != TRUE && cond.state.prob != TRUE){

      if(do.pdf == TRUE){

        pdf(paste(ID,"_Results_AB_",class(model.name@QhatModel.object)[1],"_states-",model.name@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
        viterbi(model.name, do.plot = T, plot.options = c("A","B"))
        title(ID)
        dev.off()
        return(viterbi(model.name, do.plot = T, plot.options = c("A","B")))

      } else {
        return(viterbi(model.name, do.plot = T, plot.options = c("A","B")))
      }
    }

      # plot only A and B and C
      if(ind.variable == TRUE && dep.variable == TRUE && dep.variable.transformed == TRUE && cond.state.prob != TRUE){

        if(do.pdf == TRUE){

          pdf(paste(ID,"_Results_ABC_",class(model.name@QhatModel.object)[1],"_states-",model.name@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
          viterbi(model.name, do.plot = T, plot.options = c("A","B","C"))
          title(ID)
          dev.off()
          return(viterbi(model.name, do.plot = T, plot.options = c("A","B","C")))

        } else {
          return(viterbi(model.name, do.plot = T, plot.options = c("A","B","C")))
        }
      }

    # plot only A and D
    if(ind.variable == TRUE && dep.variable != TRUE && dep.variable.transformed != TRUE && cond.state.prob == TRUE){

      if(do.pdf == TRUE){

        pdf(paste(ID,"_Results_AD_",class(model.name@QhatModel.object)[1],"_states-",model.name@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
        viterbi(model.name, do.plot = T, plot.options = c("A","D"))
        title(ID)
        dev.off()
        return(viterbi(model.name, do.plot = T, plot.options = c("A","D")))

      } else {
        return(viterbi(model.name, do.plot = T, plot.options = c("A","D")))
      }
    }

    # plot only A and C
    if(ind.variable == TRUE && dep.variable != TRUE && dep.variable.transformed == TRUE && cond.state.prob != TRUE){

      if(do.pdf == TRUE){

        pdf(paste(ID,"_Results_AC_",class(model.name@QhatModel.object)[1],"_states-",model.name@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
        viterbi(model.name, do.plot = T, plot.options = c("A","C"))
        title(ID)
        dev.off()
        return(viterbi(model.name, do.plot = T, plot.options = c("A","C")))

      } else {
        return(viterbi(model.name, do.plot = T, plot.options = c("A","C")))
      }
    }

    # plot only A , C, D
    if(ind.variable == TRUE && dep.variable != TRUE && dep.variable.transformed == TRUE && cond.state.prob == TRUE){

      if(do.pdf == TRUE){

        pdf(paste(ID,"_Results_ACD_",class(model.name@QhatModel.object)[1],"_states-",model.name@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
        viterbi(model.name, do.plot = T, plot.options = c("A","C","D"))
        title(ID)
        dev.off()
        return(viterbi(model.name, do.plot = T, plot.options = c("A","C","D")))

      } else {
        return(viterbi(model.name, do.plot = T, plot.options = c("A","C","D")))
      }
    }

    # plot only B
    if(ind.variable != TRUE && dep.variable == TRUE && dep.variable.transformed != TRUE && cond.state.prob != TRUE){

      if(do.pdf == TRUE){

        pdf(paste(ID,"_Results_B_",class(model.name@QhatModel.object)[1],"_states-",model.name@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
        viterbi(model.name, do.plot = T, plot.options = c("B"))
        title(ID)
        dev.off()
        return(viterbi(model.name, do.plot = T, plot.options = c("B")))

      } else {
        return(viterbi(model.name, do.plot = T, plot.options = c("B")))
      }
    }

    # plot only B and C
    if(ind.variable != TRUE && dep.variable == TRUE && dep.variable.transformed == TRUE && cond.state.prob != TRUE){

      if(do.pdf == TRUE){

        pdf(paste(ID,"_Results_BC_",class(model.name@QhatModel.object)[1],"_states-",model.name@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
        viterbi(model.name, do.plot = T, plot.options = c("B","C"))
        title(ID)
        dev.off()
        return(viterbi(model.name, do.plot = T, plot.options = c("B","C")))

      } else {
        return(viterbi(model.name, do.plot = T, plot.options = c("B","C")))
      }
    }

    # plot only B and D
    if(ind.variable != TRUE && dep.variable == TRUE && dep.variable.transformed != TRUE && cond.state.prob == TRUE){

      if(do.pdf == TRUE){

        pdf(paste(ID,"_Results_BD_",class(model.name@QhatModel.object)[1],"_states-",model.name@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
        viterbi(model.name, do.plot = T, plot.options = c("B","D"))
        title(ID)
        dev.off()
        return(viterbi(model.name, do.plot = T, plot.options = c("B","D")))

      } else {
        return(viterbi(model.name, do.plot = T, plot.options = c("B","D")))
      }
    }

    # plot only B,  C and D
    if(ind.variable != TRUE && dep.variable == TRUE && dep.variable.transformed == TRUE && cond.state.prob == TRUE){

      if(do.pdf == TRUE){

        pdf(paste(ID,"_Results_BCD_",class(model.name@QhatModel.object)[1],"_states-",model.name@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
        viterbi(model.name, do.plot = T, plot.options = c("B","C","D"))
        title(ID)
        dev.off()
        return(viterbi(model.name, do.plot = T, plot.options = c("B","C","D")))

      } else {
        return(viterbi(model.name, do.plot = T, plot.options = c("B","C","D")))
      }
    }

    # plot only  C
    if(ind.variable != TRUE && dep.variable != TRUE && dep.variable.transformed == TRUE && cond.state.prob != TRUE){

      if(do.pdf == TRUE){

          pdf(paste(ID,"_Results_C_",class(model.name@QhatModel.object)[1],"_states-",model.name@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
        viterbi(model.name, do.plot = T, plot.options = c("C"))
        title(ID)
        dev.off()
        return(viterbi(model.name, do.plot = T, plot.options = c("C")))

      } else {
        return(viterbi(model.name, do.plot = T, plot.options = c("C")))
      }
    }

    # plot only  C & D
    if(ind.variable != TRUE && dep.variable != TRUE && dep.variable.transformed == TRUE && cond.state.prob == TRUE){

      if(do.pdf == TRUE){

        pdf(paste(ID,"_Results_CD_",class(model.name@QhatModel.object)[1],"_states-",model.name@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
        viterbi(model.name, do.plot = T, plot.options = c("C","D"))
        title(ID)
        dev.off()
        return(viterbi(model.name, do.plot = T, plot.options = c("C","D")))

      } else {
        return(viterbi(model.name, do.plot = T, plot.options = c("C","D")))
      }
    }

    # plot only   D
    if(ind.variable != TRUE && dep.variable != TRUE && dep.variable.transformed != TRUE && cond.state.prob == TRUE){

      if(do.pdf == TRUE){

        pdf(paste(ID,"_Results_D_",class(model.name@QhatModel.object)[1],"_states-",model.name@QhatModel.object@nStates,".pdf",sep=""), width = 8.5, height = 11)
        viterbi(model.name, do.plot = T, plot.options = c("D"))
        title(ID)
        dev.off()
        return(viterbi(model.name, do.plot = T, plot.options = c("D")))

      } else {
        return(viterbi(model.name, do.plot = T, plot.options = c("D")))
      }
    }


}

#'get states
#'
#' \code{get.states}
#'
#' @description
#' \code{get.states} retrieves results from the fitted \code{hydroState} model.
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
#' @param model.name is the name of the fitted \code{hydroState} model object.
#'
#' @return
#' data frame of results to evaluate the rainfall-runoff states overtime
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
#'                 setInitialYear(model.name = model.annual.fitted.221201,
#'                                initial.year = 1990)
#'
#' ## Get states
#' model.annual.fitted.221201.states =
#'                 get.states(model.name = model.annual.fitted.221201)
#'


get.states <- function(model.name){

      return(viterbi(model.name, do.plot = F))

}



