#' Transforms Observations
#'
#' \code{select.transform}
#'
#' @description
#' Transforms observations to remove heteroscedasticity.
#'
#' @details
#' Often there is skew within hydrologic data. When defining relationships between observations, this skew results in an unequal variance in the residuals, heteroscedasticity. Transforming observations is often required with observations of streamflow, precipitation, and concentration. Qhat provides several options to transform observations. Since the degree of transformation is not typically known, 'boxcox' is the default. Other options include: 'log', 'burbidge', and of course, 'none' when no transformation is performed.
#'
#' For an example of how to transform observations.. see vignette...
#'
#'
#'
#' @param func is the method of transformation. The default is 'boxcox'. Other options: 'log', 'burbidge', 'none'
#' @param input.data dataframe of annual or monthly runoff and precipitation observations.
#'
#' @return
#' A Qhat object with transformed observations ready for \code{buildModel}
#'
#' @keywords Qhat transform
#'
#'
#' @export select.transform
#' @importFrom methods new

select.transform <- function(func = 'boxcox', input.data=data.frame(year=c(), flow=c(), precip=c())){

  #Validate
  func = paste('Qhat.',func,sep='')
  if(func %in% c('Qhat.none','Qhat.log','Qhat.burbidge','Qhat.boxcox')){

    if('month' %in% colnames(input.data)){ # sort in ascending order by year and month
      input.data = input.data[order(input.data[,'year'],input.data[,'month']),]
    }


    return(new(func, input.data))

  }else{

    stop('Please choose an appropriate function: none, log, burbidge, or boxcox')
  }


}

#'
#'
#' \code{select.stateModel}
#'
#' @description
#' \code{select.stateModel} provides various options for constructing the rainfall-runoff relationship. Every \code{stateModel} depends on a linear base model where precipitation is a function of streamflow. The default model is an annual analysis of this base linear model with state shifts expected in the intercept, \eqn{a_0}, and standard deviation, \eqn{std}, of the rainfall-runoff relationship. \code{select.stateModel} provides additional adjustments to this model with auto-correlations terms, seasonal parameters for a sub-annual analysis, and even other or multiple \code{state.shift.parameters}. The number of states and assumed error distribution can also be selected.
#'
#' @details
#' There are a selection of items to consider when defining the rainfall-runoff relationship and investigating state shifts in this relationship. hydroState simulates runoff, \eqn{Q}, as being in one of finite states, \eqn{i}, at every time-step, \eqn{t}, depending on the distribution of states at prior time steps. This results in a runoff distribution for each state that can vary overtime (\eqn{\widehat{_tQ_i}}). The \code{stateModel} defines the relationship that is susceptible to state shifts with precipitation, \eqn{P_t}, as a predictor. This takes the form as a simple linear model \eqn{\widehat{_tQ_i} = f(P_t)}:
#'
#'                    \eqn{\widehat{_tQ_i} = P_ta_0 + a_1}
#'
#' where \eqn{a_0} and \eqn{a_1} are constant parameters. These parameters and the model error, \eqn{std},  are required parameters for every \code{stateModel}. It is possible the relationship contains serial correlation and would be better defined with an auto-regressive term:
#'
#'                \eqn{\widehat{_tQ_i} = P_ta_0 + a_1 + AR1\widehat{_{t-1}Q}}
#'
#' where \eqn{AR1} is the lag-1 auto-correlation term. Either, lag-1: \eqn{AR1}, lag-2: \eqn{AR2}, and lag-3: \eqn{AR3} auto-correlation coefficients are an option as additional parameters to better define the rainfall-runoff relationship.
#' For sub-annual analysis, seasonal.parameters provides the option to assume a sinusoidal function better defines either of the constant parameters or error (\eqn{a_0, a_1, std}) throughout the year, i.e:
#'
#' \eqn{a_0 = a_{0.disp} + a_{0.amp} * sin(2\pi(\frac{M_t}{12} + a_{0.phase}))}
#'
#' where \eqn{M_t} is an integer month at \eqn{t}. Monthly streamflow and precipitation are required as input.data for the sub-annual analysis.
#'
#' Once the model parameters are chosen, the \code{state.shift.parameters} can be selected to investigate shifts based on any or all of the previously chosen parameters (\eqn{a_0, a_1, std, AR1, AR2, AR3}). The selected \code{state.shift.parameters} are state dependent where they are subject to shift in order to better explain the state of streamflow. The default \code{stateModel} evaluates shifts in the rainfall-runoff relationship with \eqn{a_0} as a state dependent parameter.
#'
#' The \code{error.distribution} of the model is chosen as either normal: "normal", truncated normal: "truc.normal", or Gaussian: "gauss" in order to reduce skew. The default is "truc.normal". The number of possible states in the rainfall-runoff relationship and transition between the states is selected with the transition.graph. The default is a 2-state model in a 2 by 2 matrix with a TRUE transition to and from each state.
#'
#' @param input.data dataframe of annual runoff and precipitation observations. For running seasonal models with seasonal.parameters, monthly data is required.
#' @param parameters character list of parameters to construct state model. Required and default: \code{a0}, \code{a1}, \code{std}. Auto-correlation terms optional: \code{AR1}, \code{AR2}, or \code{AR3}.
#' @param seasonal.parameters character list of one or all parameters (\code{a0}, \code{a1}, \code{std}) defined as a sinusoidal function to represent seasonal variation. Requires monthly data. Default is empty list.
#' @param state.shift.parameters character list of one or all parameters (\code{a0}, \code{a1}, \code{std}, \code{AR1}, \code{AR2}, \code{AR3}) able to shift as dependent on state. Default is \code{a_0} and \code{std}.
#' @param error.distribution character name of the distribution in the HMM error. Default is "truc.normal". Others include: "normal" or "gauss"
#' @param transition.graph matrix given the number of states. Default is a 2-state matrix (2 by 2): matrix(TRUE,2,2)
#'
#' @return
#' A \code{stateModel}, QhatModel object, ready to for \code{buildModel}.
#'
#' @keywords stateModel linear-model rainfall-runoff QhatModel
#'
#'
#' @export select.stateModel
#' @importFrom methods new

# Create QhatModel object
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
  }

  # Curate state model
  if('AR1' %in% parameters || 'AR2' %in% parameters || 'AR3' %in% parameters){
    auto.array = c(length(which(parameters == 'AR1')) > 0,
                   length(which(parameters == 'AR2')) > 0,
                   length(which(parameters == 'AR3')) > 0)
    auto.term = paste('AR',which(auto.array == TRUE),sep='')

    # Assume seasonal model if seasona paramaters
    if(length(seasonal.parameters) > 0){

      # Add if statement here if error.distribution for seasonal accepts normal or truc.normal; gamma for now...
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

      # Add if statement here if error.distribution for seasonal accepts normal or truc.normal; gamma for now...
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

        }else{ #gamma does not need use.truncated.dist input...

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





#' \code{select.Markov} selects a hidden Markov model
#'
#' @description
#' select.Markov selects a hidden Markov model.
#'
#' @details
#' There are several hidden Markov models to choose from with different forms. The default markov model is 'annualHomogeneous' where the transition between states only depends on a sequence of observed prior states.
#'
#'
#' For an example of how to.. see vignette...
#'
#' @param func character sting with name of Markov model defining the transition between states. The default is 'annualHomogeneous'. Other options include 'annualHomogeneous.flickering'.
#' @param transition.graph matrix given the number of states. Default is a 2-state matrix (2 by 2): matrix(TRUE,2,2)
#'
#' @return
#' A Markov object with a selected Markov model
#'
#' @keywords Markov HMM homogeneous
#'
#'
#' @export select.Markov
#' @importFrom methods new

select.Markov <- function(func = 'annualHomogeneous',
                        transition.graph = matrix(TRUE,2,2)){

  #Validate
  func = paste('markov.',func,sep='')
  if(func %in% c('markov.annualHomogeneous','markov.annualHomogeneous.flickering')){

    return(new(func, transition.graph))

  }else{

    stop('Please choose an appropriate Markov model. See help')

  }

}

#'Builds hydroState model
#'
#' \code{buildModel}
#'
#' @description
#' \code{buildModel} builds a hydrostate model with either a default \code{stateModel} or the \code{stateModel} can be specified with options from \code{select.stateModel}. After the model is built, the hydroState model is ready to be fitted with \code{fitModel}
#'
#' @details
#' hydroState operates in S4, object oriented programming, and requires three objects to build a hydroState model. Each object can be selected from \code{select.transform}, \code{select.stateModel}, and \code{select.Markov}; however, if no object is defined a default model is built as discussed in the arguments.
#'
#'
#' @param input.data dataframe of annual runoff and precipitation observations. For running seasonal models with \code{seasonal.parameters}, monthly data is required.
#' @param data.transform a \code{Qhat.object} with transformed observations from \code{select.transform}. If blank, the default uses 'boxcox' to transform observations.
#' @param stateModel a \code{QhatModel.object} from \code{select.stateModel}. If blank, the default selects a 2-state 'QhatModel.homo.normal.linear' model with a truncated normal error distribution, and allows the intercept 'a0' and standard deviation 'std' to shift as state dependent parameters.
#' @param Markov a \code{markov.model.object} from \code{select.Markov}. If blank, the default selects a homogeneous Markov model without flickering.
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


buildModel <- function(input.data = data.frame(year=c(), flow=c(), precip=c()),
                       data.transform = NULL,
                       stateModel = NULL,
                       Markov = NULL){


  #Validate input.data
  if(!is.data.frame(input.data)){
    stop("'input.data' is not a dataframe")
  }
  if(!('year' %in% colnames(input.data))){
    stop("'input.data' must contain a 'year' column with an integer of years")
  }

  # If monthly data, sort in ascending order by year and month
  if('month' %in% colnames(input.data)){
    input.data = input.data[order(input.data[,'year'],input.data[,'month']),]
  }

# If no inputs except input.data, just run default annual analysis
  if(is.null(data.transform) && is.null(stateModel) && is.null(Markov)){

    data.transform = select.transform(func = 'boxcox', input.data = input.data)

    stateModel = select.stateModel(input.data = input.data,
                                   parameters = list('a0', 'a1', 'std'),
                                   seasonal.parameters = list(),
                                   state.shift.parameters = list('a0','std'),
                                   error.distribution = "truc.normal",
                                   transition.graph = matrix(TRUE, 2, 2))

    Markov = select.Markov(func = 'annualHomogeneous',
                                     transition.graph = matrix(TRUE, 2, 2))

    if(stateModel@nStates != nrow(Markov@transition.graph)){
      stop("number of possible states in stateModel and 'transition.graph' in Markov are not equal")
    }

    return(new('hydroState',input.data, data.transform, stateModel, Markov))

  }else{

# if no stateModel or Markov still run
    if(is.null(data.transform) && !is.null(stateModel) && !is.null(Markov)){

      data.transform = select.transform(func = 'boxcox', input.data = input.data)

      if(!is(stateModel,"QhatModel")){
        stop("'stateModel' must be a QhatModel object from select.stateModel")
      }
      if(!is(Markov,"markov")){
        stop("'Markov' must be a markov object from select.Markov")
      }

      if(stateModel@nStates != nrow(Markov@transition.graph)){
        stop("number of possible states in stateModel and 'transition.graph' in Markov are not equal")
      }

      if(!identical(stateModel@input.data,input.data)){
        stop("input.data is not the same input.data in the stateModel")
      }

      return(new('hydroState',input.data, data.transform, stateModel, Markov))
    }

# if no transform data and state model, still run
    if(is.null(data.transform) && !is.null(stateModel) && is.null(Markov)){

      data.transform = select.transform(func = 'boxcox', input.data = input.data)

      Markov = select.Markov(func = 'annualHomogeneous',
                                       transition.graph = matrix(TRUE,stateModel@nStates,stateModel@nStates))

      if(!is(stateModel,"QhatModel")){
        stop("'stateModel' must be a QhatModel object from select.stateModel")
      }

      if(stateModel@nStates != nrow(Markov@transition.graph)){
        stop("number of possible states in stateModel and 'transition.graph' in Markov are not equal")
      }

      if(!identical(stateModel@input.data,input.data)){
        stop("input.data is not the same input.data in the stateModel")
      }

      return(new('hydroState',input.data, data.transform, stateModel, Markov))
    }

# if no transform data or state model, still run
    if(!is.null(data.transform) && !is.null(stateModel) && is.null(Markov)){

      if(!is(data.transform,"Qhat")){
        stop("'data.transform' must be a Qhat object from select.transform")
      }
      if(!is(stateModel,"QhatModel")){
        stop("'stateModel' must be a QhatModel object from select.stateModel")
      }

      Markov = select.Markov(func = 'annualHomogeneous',
                                       transition.graph = matrix(TRUE,stateModel@nStates,stateModel@nStates))

      if(stateModel@nStates != nrow(Markov@transition.graph)){
        stop("number of possible states in stateModel and 'transition.graph' in Markov are not equal")
      }

      return(new('hydroState',input.data, data.transform, stateModel, Markov))
    }

# if no markov model, still run
    if(is.null(data.transform) && is.null(stateModel) && !is.null(Markov)){

      if(!is(Markov,"markov")){
        stop("'Markov' must be a markov object from select.Markov")
      }

      data.transform = select.transform(func = 'boxcox', input.data = input.data)

      stateModel = select.stateModel(input.data = input.data,
                                     parameters = list('a0', 'a1', 'std'),
                                     seasonal.parameters = list(),
                                     state.shift.parameters = list('a0','std'),
                                     error.distribution = "truc.normal",
                                     transition.graph = Markov@transition.graph)

      if(stateModel@nStates != nrow(Markov@transition.graph)){
        stop("number of possible states in stateModel and 'transition.graph' in Markov are not equal")
      }


      return(new('hydroState',input.data, data.transform, stateModel, Markov))
    }

    # if no state model, still run
    if(!is.null(data.transform) && is.null(stateModel) && !is.null(Markov)){

      if(!is(data.transform,"Qhat")){
        stop("'data.transform' must be a Qhat object from select.transform")
      }

      if(!is(Markov,"markov")){
        stop("'Markov' must be a markov object from select.Markov")
      }


      stateModel = select.stateModel(input.data = input.data,
                                     parameters = list('a0', 'a1', 'std'),
                                     seasonal.parameters = list(),
                                     state.shift.parameters = list('a0','std'),
                                     error.distribution = "truc.normal",
                                     transition.graph = Markov@transition.graph)

      if(stateModel@nStates != nrow(Markov@transition.graph)){
        stop("number of possible states in stateModel and 'transition.graph' in Markov are not equal")
      }

      if(!identical(data.transform@input.data,input.data)){
        stop("input.data is not the same input.data in the data.transform")
      }


      return(new('hydroState',input.data, data.transform, stateModel, Markov))
    }

  #if no anything, still run with inputs
    if(!is.null(data.transform) && !is.null(stateModel) && !is.null(Markov)){

      if(!is(data.transform,"Qhat")){
        stop("'data.transform' must be a Qhat object from select.transform")
      }

      if(!is(stateModel,"QhatModel")){
        stop("'stateModel' must be a QhatModel object from select.stateModel")
      }
      if(!is(Markov,"markov")){
        stop("'Markov' must be a markov object from select.Markov")
      }
      if(stateModel@nStates != nrow(Markov@transition.graph)){
        stop("number of possible states in stateModel and 'transition.graph' in Markov are not equal")
      }

      if(!identical(stateModel@input.data,input.data)){
        stop("input.data is not the same input.data in the stateModel")
      }

      if(!identical(data.transform@input.data,input.data)){
        stop("input.data is not the same input.data in the data.transform")
      }

      return(new('hydroState',input.data, data.transform, stateModel, Markov))
    }
  }
}


#' Builds all hydroState models # work in progress..
#'
#' \code{buildModelAll}
#'
#' @description
#' \code{buildModelAll} builds all hydroState models
#'
#' @details
#' All hydroState models are build with different QhatModel objects
#'
#' For an example of how to.. see vignette...
#'
#' @param ID  character vector of a stream gauge identifier
#' @param input.data is a dataframe of annual, monthly, or daily, observations of runoff, precipitation, or concentration. For running seasonal models, monthly values are required.
#'
#'
#' @return
#' A list of built hydroState models with every combination of QhatModel, ready to be fitted
#'
#' @keywords hydroState build all
#'
#'
#' @export
#' @importFrom methods new


buildModelAll <- function(ID = '',
                          input.data = data.frame(year=c(), flow=c(), precip=c())){
  # Validate
  if(is.character(ID)){

    # If monthly data, sort in ascending order by year and month
    if('month' %in% colnames(input.data)){
      input.data = input.data[order(input.data[,'year'],input.data[,'month']),]
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
#' \code{fitModel} fits hydrostate model(s) using global optimization by differential evolution \link[DEoptim]{DEoptim}.
#'
#' @details
#' After a hydroState model object is built, the model is ready to be fitted. The only required input is the given name of the built hydroState model object. \code{fitModel} works for one built model (\code{buildModel}) or all (\code{buildModelAll}). If fitting all models be sure to install and load the \href{https://cran.r-project.org/web/packages/parallelly/index.html}{parallelly} library.
#'
#'
#' @param model.name name of the built hydroState model or name of the list containing all built models
#' @param pop.size.perParameter integer that should be greater than or equal to the number of parameters in the model. The default is '10' and is sufficient for all models.
#' @param max.generations integer that will stop the optimizer when set number of generations are reached. The default is '500'.
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
fitModel <- function(model.name = model,
                     pop.size.perParameter = 10,
                     max.generations=500){

  # Validate
  if(class(model.name) %in% c("hydroState", "hydroState.allModels", "hydroState.subAnnual.allModels")){

    # fit model
    return(fit(model.name, pop.size.perParameter = pop.size.perParameter, max.generations=max.generations))

  }else{

    stop('model.name is not an appropriate class. Please ensure input model is a built hydroState model with the class as either of the following: "hydroState", "hydroState.allModels", or "hydroState.subAnnual.allModels"')

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
#' \code{plot.residuals} produces five plots to review and validate the fitted hydroState model. A) Time-series of normal-pseudo residuals to ensure the residuals each year are within the confidence intervals. B) Auto-correlation function (ACF) of normal-pseudo residuals to ensure there is no serial correlation in residuals. Lag spikes should be below confidence interval at each lag (except 0). C) Histogram of uniform-pseudo residuals should show uniform distribution (equal frequency for each residual value) D) Histogram of normal-pseudo residuals should show normal distribution centered on zero and with no skew E) Quantile-Quantile (Q-Q) plot where normal-pseudo residuals vs. theoretical quantities should align on the diagonal line. The last plot contains the Akaike information criterion (AIC) and Shapiro-Wilk p-value. The AIC is an estimator to determine the most parsimonious, best performing model given the number of parameters. When comparing models, the lowest AIC is the best performing model. Shapiro-Wilks test for normality in the residuals and a p-value greater than 0.05 (chosen alpha level) indicates the residuals are normally distributed; the null hypothesis that the residuals are normally distributed is not rejected.
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


plot.residuals <- function(model.name = model,
                           do.pdf = FALSE,
                           ID = NULL
                             ){

  if(do.pdf != TRUE){
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


get.residuals <- function(model.name = model){

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
#' @param initial.year year (YYYY). Default is first year in input.data.
#'
#' @return
#' A fitted hydroState model object with state names for each time-step
#'
#' @keywords state names
#'
#'
#' @export
#'


setInitialYear <- function(model.name = model,
                       initial.year = input.data$year[1]){ #make go to first year of dataframe

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
#' \code{plot.states} produces plots of the results from the fitted hydroState model. The default produces four plots of the A) independent variable, B) dependent variable and states, C) transformed dependent variable and states, and D) conditional state probabilities for each state. These are plotted on the same page, and there is an option to export plots as a pdf to the current working directory. There are also options to only plot one of the four plots.
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
#' plots to evaluate rainfall-runoff states overtime along with observations and the conditional probabilities of each state. The data frame includes several items. These include the time-step (i.e. Year, Month), Viterbi State Number to differentiate states, observations (i.e. stream flow), flow values of the Viterbi state including the 5\% and 95\% confidence intervals, flow values of the normal state including the 5% and 95% confidence intervals, conditional probabilities for each state, and emission density for each state. The Viterbi state flow values are the results, the values of the most likely state at each time-step. The Normal state flow values are the values from the normal state at each time-step. When the most likely state is the Normal state for a time-step, the Viterbi flow state value equals the Normal flow state value. These Normal state values are interpreted as the results for a one-state model, if the multi-state model did not exist. This estimated Normal state can be seen in plot (B) relative to the other states. The Conditional Probabilities of each state show the likelihood of remaining in the given state based on previous states. When the conditional probability is closer to 1, there is a higher probability that hydroState remains in the state for the next time-step. The Emission Density of each state is the result of multiplying the conditional probabilities by the transition probabilities.
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


plot.states <- function(model.name = model,
                        ind.variable = TRUE,
                        dep.variable = TRUE,
                        dep.variable.transformed = TRUE,
                        cond.state.prob = TRUE,
                        do.pdf = FALSE,
                        ID = NULL
                       ){

    # plot everything, default, export as pdf too..
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
#' @details
#' \code{get.states} These results include the time-step (i.e. Year, Month), Viterbi State Number to differentiate states, observations (i.e. stream flow), flow values of the Viterbi state including the 5\% and 95\% confidence intervals, flow values of the normal state including the 5% and 95% confidence intervals, conditional probabilities for each state, and emission density for each state. The Viterbi state flow values are the results, the values of the most likely state at each time-step. The Normal state flow values are the values from the normal state at each time-step. When the most likely state is the Normal state for a time-step, the Viterbi flow state value equals the Normal flow state value. These Normal state values are interpreted as the results for a one-state model, if the multi-state model did not exist. This estimated Normal state can be seen in plot (B) relative to the other states. The Conditional Probabilities of each state show the likelihood of remaining in the given state based on previous states. When the conditional probability is closer to 1, there is a higher probability that hydroState remains in the state for the next time-step. The Emission Density of each state is the result of multiplying the conditional probabilities by the transition probabilities.
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


get.states <- function(model.name = model){

      return(viterbi(model.name, do.plot = F))

}



