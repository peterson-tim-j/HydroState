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

  return(streamflow_monthly.seasonal)
}

#' @keywords internal
select.transform <- function(func = 'boxcox', input.data=data.frame(year=c(), flow=c(), precip=c()), data.transform.constant=1){

  #Validate
  func = paste('Qhat.',func,sep='')
  if(func %in% c('Qhat.none','Qhat.log','Qhat.burbidge','Qhat.boxcox')){

    return(new(func, input.data, constant = data.transform.constant))

  }else{

    stop('Please choose an appropriate function: none, log, burbidge, or boxcox')
  }


}

#' @keywords internal
select.stateModel <- function( input.data,
                               formula.mean = as.formula('q ~ a0 + a1*p'),
                               formula.var = as.formula('e ~ b0'),
                               shift.pars = c('a0', 'b0'),
                               error.distribution = 'trunc.norm',
                               transition.graph=matrix(T,2,2)
                              ) {

  #Validate input.data
  if(!is.data.frame(input.data)){
    stop("'input.data' is not a dataframe")
  }
  if(!('year' %in% colnames(input.data))){
    stop("'input.data' must contain a 'year' column with an integer of years")
  }

  #Validate state.shift.parameters
  if(!is.character(shift.pars)){
    stop("'state.shift.parameters' are not character vector")
  }

  #Validate error.distribution
  if(!is.character(error.distribution)){
    stop("'error.distribution' is not a character string")
  }
  if(!(error.distribution %in% c('normal','trunc.norm','gamma'))){
    stop("'error.distribution' must contain either 'normal', 'truc.normal' or 'gamma'")
  }

  #Validate transition.graph
  if(!is.matrix(transition.graph)){
    stop("'transition.graph' must be a matrix")
  }

  # Check mean formula
  if (!any(all.vars(formula.mean[[2]]) %in% colnames(input.data)))
    stop('The left-hand predictant in the mean formula must be in the input data.frame.')

  # If monthly data, sort in ascending order by year and month
  if('month' %in% colnames(input.data)){
    input.data = input.data[order(input.data[,'year'],input.data[,'month']),]
  }else if('day' %in% colnames(input.data)){
    input.data = input.data[order(input.data[,'year'],input.data[,'month'],input.data[,'day']),]
  }

  # Curate state model
  return(new('QhatModel',
             input.data,
             formula.mean,
             formula.var,
             shift.pars,
             error.distribution,
             transition.graph=transition.graph
             )
         )
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
#' \code{build} builds a hydrostate model object with either a default model or the model can be specified with options from below. Every model depends on a linear base model where streamflow, \eqn{Q}{Q}, is a function of precipitation,  \eqn{P}{P}: \eqn{\hat{Q} = Pa_1+a_0}{Q = P * a1 + a0}. The default model is a variant of this base linear model with state shifts expected in the intercept, \eqn{a_0}{a0}, and standard deviation, \eqn{std}{std}, of the rainfall-runoff relationship. There are additional options to adjust this model with auto-correlation terms and seasonal parameters that can be both independent of state or change with state. The number of states and assumed error distribution can also be selected. After the model is built, the hydroState model is ready to be fitted with \code{\link{fit.hydroState}}
#'
#' @details
#' There are a selection of items to consider when defining the rainfall-runoff relationship and investigating state shifts in this relationship. hydroState provides various options for modelling the rainfall-runoff relationship.
#' \itemize{
#' \item{Data gaps  with \code{input.data}}: When there is missing \code{input.data} in either the dependent variable, streamflow, or independent variable, precipitation, the emissions probability of the missing time-step is set equal to one. This essentially ignores the missing periods. The time step after the missing period has a state probability dependent on the length of the gap. The larger the gap, the closer the state probability gets to approaching a finite probability near zero (as the transition probabilities are recursively multiplied). When the model has auto-correlation terms and there are gaps in the dependent variable, the auto-correlation function restarts at the beginning of each continuous period after the gap. This ignores auto-correlation at the first time steps after the gap. For instance, an 'AR1' model would ignore the contribution of the prior time step for the first (1) observation after the gap.
#' \item{Transform Observations with \code{data.transform}}: Transforms streamflow observations to remove heteroscedasticity. Often there is skew within hydrologic data. When defining relationships between rainfall-runoff, this skew results in an unequal variance in the residuals, heteroscedasticity. Transforming streamflow observations is often required. There are several options to transform observations. Since the degree of transformation is not typically known, \code{boxcox} is the default. Other options include: \code{log}, \code{burbidge}, and of course, \code{none} when no transformation is performed.
#' \item{Model Structure with \code{parameters} and \code{seasonal.parameters}}: The structure of the model depends on the \code{parameters}. hydroState simulates runoff, \eqn{Q}, as being in one of a finite states, \eqn{i}, at every time-step, \eqn{t}, depending on the distribution of states at prior time steps. This results in a runoff distribution for each state that can vary over time (\eqn{\hat{_tQ_i}}). The model defines the relationship that is susceptible to state shifts with precipitation, \eqn{P_t}, as a predictor. This takes the form as a simple linear model \eqn{\hat{_tQ_i} = f(P_t)}:
#'
#'
#' \deqn{
#' \hat{_tQ_i} = P_ta_1 + a_0
#' }
#'
#' where \eqn{a_0} and \eqn{a_1} are constant parameters. These parameters and the model error, \eqn{std}, establish the rainfall-runoff relationship and are required parameters for every built model object. These are the default parameters: \code{c('a0', 'a1', 'std')}.
#' \item{Auto-correlation \code{parameters}}: The relationship may contain serial correlation and would be better defined with an auto-regressive term:
#'
#' \deqn{
#' \hat{_tQ_i} = P_ta_1 + a_0 + AR1\hat{_{t-1}Q}
#' }
#'
#' where \code{AR1} is the lag-1 auto-correlation term. Either, lag-1: \code{AR1}, lag-2: \code{AR2}, and lag-3: \code{AR3} auto-correlation coefficients are an option as additional parameters to better define the rainfall-runoff relationship.
#' \item{Sub-annual analysis with \code{seasonal.parameters}}: Additional options include explaining the seasonal rainfall-runoff relationship with a sinusoidal function that better defines either of the constant parameters or error (\eqn{a_0, a_1, std}) throughout the year, i.e:
#'
#' \deqn{
#' a_0 = a_{0.disp} + a_{0.amp} * sin(2\pi(\frac{M_t}{12} + a_{0.phase}))
#' }
#'
#' where \eqn{M_t} is an integer month at \eqn{t}. Monthly streamflow and precipitation are required as \code{input.data} for the sub-annual analysis.
#'
#' \item{State Dependent Parameters with \code{state.shift.parameters}}: These are state dependent parameters where they are subject to shift in order to better explain the state of streamflow over time. Any or all of the previously chosen parameters can be selected (\code{a_0, a_1, std, AR1, AR2, AR3}). The default model evaluates shifts in the rainfall-runoff relationship with \code{a_0} \code{std} as state dependent parameters.
#' \item{Distribution of the Residuals with \code{error.distribution}}: The distribution of the residuals (error) within a state of the model can be chosen to reduce skew and assist with making models statistically adequate (see \code{plot(pse.residuals = TRUE)}). Either normal: \code{normal}, truncated normal: \code{truc.normal}, or gamma: \code{gamma} distributions are acceptable. These error distribution ensures streamflow is greater than zero \code{Q > 0}, and specifically for \code{truc.normal} greater than or equal to zero \code{Q >= 0}. The default is \code{truc.normal}. Sub-annual models are restricted to only a \code{gamma} distribution.
#' \item{Markov flickering with \code{flickering}}: When flickering is \code{FALSE}, the markov avoids state shifts for very short duration, and hence for a state shift to occur it should last for an extended period. The default is FALSE. If TRUE, flickering between more states is more sensitive. For further explanation on this method, see: \href{https://hess.copernicus.org/articles/7/652/2003/}{Lambert et al., 2003}. The current form of the Markov model is homogeneous where the transition probabilities are time-invariant.
#' \item{Number of States with \code{transition.graph}}: The number of possible states in the rainfall-runoff relationship and transition between the states is selected with the transition.graph. The default is a 2-state model in a 2 by 2 unstructured matrix with a TRUE transition to and from each state (i.e. \code{matrix(TRUE,2,2)}). hydroState accepts 1-state up to 3-states (i.e. for 3-state unstructured transition graph: \code{matrix(TRUE,3,3)}). The unstructured transition graph allows either state to remain in the current state or transition between any state. For the 3-state transition graph, one may want to assume the transitions can only occur in a particular order, as in (Very low -> Low -> Normal->) rather than (Very low <-> Low <-> Normal <-> Very low). Thus, a structured graph is also acceptable (3-state structured: \code{matrix(c(TRUE,TRUE,FALSE,FALSE,TRUE,TRUE,TRUE,FALSE,TRUE),3,3)}). More details in Supplementary Materials of (Peterson TJ, Saft M, Peel MC & John A (2021), Watersheds may not recover from drought, Science, DOI: \doi{10.1126/science.abd5085}).
#'}
#'
#' @param input.data dataframe of annual, seasonal, or monthly runoff and precipitation observations. Gaps with missing data in either streamflow or precipitation are permitted. Monthly data is required when using \code{seasonal.parameters}.
#' @param data.transform character string with the method of transformation. The default is 'boxcox'. Other options: 'log', 'burbidge', 'none'
#' @param data.transform.constant scalar numeric to allow the transform to handle zero flow. Default is 1.
#' @param formula.mean formula for the mean value at each time step. Parameters must be explicit in the
#' formula. Required and default: \code{q ~ a0 + a1 * P}. Any right hand side terms not within \code{input.date} are
#' assumed to be parameters. Hence, errors cannot be thrown if a formula term intended to be a variable is missing
#' from \code{input.date}. Auto-correlation terms optional: \code{AR1}, \code{AR2}, or \code{AR3}. Any or all AR terms
#' can be used, noting that a parameter is expected to be input before the AR term (e.g.\code{a3 * AR1}).
#' @param formula.var formula for the error variance at each time step. Parameters must be explicit in the formula.
#' Required and the default is \code{v ~ b0 * var(q.obs)} where \code{var(q.obs)} is the variance of the
#' observed predictant and is used to allow the b0 parameter to be between zero and one.
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
#' # Build hydroState model with: 2-state, truncated normal error distribution,
#' # 1-lag of auto-correlation, and state dependent parameters ('a0', 'b0')
#' mod <- build(input.data = obs.data,
#'              data.transform = "boxcox",
#'              formula.mean = as.formula(q ~ a0 + a1*p + a2*AR1),
#'              formula.var = as.formula(v ~ b0 * var(q.obs)),
#'              state.shift.parameters = c("a0", "b0"),
#'              error.distribution = "trunc.norm",
#'              transition.graph = matrix(TRUE, 2, 2)
#'              )
#'
build <- function(input.data = data.frame(year=c(), flow=c(), precip=c()),
                  data.transform = 'boxcox',
                  data.transform.constant = 1,
                  formula.mean = as.formula('q ~ a0 + a1*p'),
                  formula.var = as.formula('e ~ b0'),
                  state.shift.parameters =  c('a0', 'b0'),
                  error.distribution = 'trunc.norm',
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
    data.transform = select.transform(func = 'boxcox', input.data, data.transform.constant)
  }else{
    data.transform = select.transform(data.transform, input.data, data.transform.constant)
  }

  # create state model
  stateModel = select.stateModel(input.data,
                                 formula.mean,
                                 formula.var,
                                 state.shift.parameters,
                                 error.distribution,
                                 transition.graph)

  # create Markov model
  Markov = select.Markov(flickering, transition.graph)

  # build default model
  return(new('hydroState',input.data, data.transform, stateModel, Markov))
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
#'   \item{\eqn{\delta} is the initial state distribution, the initial probability of being in each state: \eqn{\delta = \begin{pmatrix} \delta_{1} \\ 1- \delta_{1} \end{pmatrix}}}
#'   \item{\eqn{P(x)} is the \eqn{m} x \eqn{m} diagonal emissions matrix of the probability density for each state using a lower tail truncated Gaussian distribution or a two-parameter Gamma distribution:
#'     \itemize{
#'       \item{\eqn{f_{Gau}(x=\hat{_{obs}q_{t}}; \mu = \hat{_{t}q_{i}}, \sigma = \sigma_{i}, a = 0) = \frac{\phi(\frac{x-\mu}{\sigma})}{\sigma(1-\Phi(\frac{a-\mu}{\sigma}))}}}
#'       \item{\eqn{f_{Gam}(x = \hat{_{obs}q_{t}}; k = \frac{\hat{_{t}q_{i}}^2}{\sigma_{i}^2}, \theta = \frac{\sigma_{i}^2}{\hat{_{t}q_{i}}}) = \frac{x^{k-1}e^{\frac{x}{\theta}}}{\theta^{k}\Gamma(k)}}}
#'       \item{where \eqn{\phi} is the probability density function for the standard normal distribution, \eqn{\Phi} is the cumulative distribution function for the standard normal distribution, \eqn{k} is the shape parameter, \eqn{\theta} is the scale parameter, and \eqn{\Gamma(k)} is the gamma function. For more details, refer to pg. 8–17 in Supplementary Materials of (Peterson TJ, Saft M, Peel MC & John A (2021), Watersheds may not recover from drought, Science, DOI: \doi{10.1126/science.abd5085}).}
#'     }
#'   }
#'   \item{\eqn{\Gamma} is the transition matrix}
#'   \item{\eqn{T} is the number of time-steps.}
#' }
#'
#' @param model built hydroState model object, hydroState.allModels object, or hydroState.subAnnual.allModels object
#' @param pop.size.perParameter integer that should be greater than or equal to the number of parameters in the model. The default is '10' and is sufficient for all models.
#' @param max.generations integer that will stop the optimizer when set number of generations are reached. The default is '500'.
#' @param doParallel TRUE/FALSE to perform fitting in parallel on all computer cores. Default is FALSE
#' @param ... additional options to change the optimization settings: reltol, print.iterations, etc. from the \href{https://cran.r-project.org/package=DEoptim}{DEoptim} library
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
#' \donttest{
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
#' # Fit all (runtime > several hours)
#' \donttest{
#' all.annual.models = fit.hydroState(all.annual.models)
#' }
#'

fit.hydroState <- function(model,
                pop.size.perParameter = 10,
                max.generations = 500,
                doParallel = FALSE,
                ...
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
#' @return AIC value of a single model or a list variable of AIC values for all models
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
#' \code{plot} produces five figures of psuedo residuals OR up to four figures of the results from the fitted hydroState model. When the \code{pse.residuals} is FALSE, the default \code{plot} produces all four markov state plots (ind.variable, dep.variable, dep.variable.transformed, cond.state.prob). Figures are more easily viewed as an exported pdf when directory/file name is given (i.e. \code{file = "flow.state.plots.siteID.pdf"}).
#' \itemize{
#'  \item{psuedo residual figures}
#'    \itemize{
#'     \item{A) Time-series of normal-pseudo residuals to ensure the residuals each year are within the confidence intervals.}
#'     \item{B) Auto-correlation function (ACF) of normal-pseudo residuals to ensure there is minimal serial correlation in residuals. Lag spikes should be below confidence interval at each lag (except 0).}
#'     \item{C) Histogram of uniform-pseudo residuals should show uniform distribution (equal frequency for each residual value)}
#'     \item{D) Histogram of normal-pseudo residuals should show normal distribution centered on zero and with no skew}
#'     \item{E) Quantile-Quantile (Q-Q) plot where normal-pseudo residuals vs. theoretical quantities should align on the diagonal line. The last plot contains the Akaike information criterion (AIC) and Shapiro-Wilk p-value. The AIC is an estimator to determine the most parsimonious, best performing model given the number of parameters. When comparing models, the lowest AIC is the best performing model. Shapiro-Wilks test for normality in the residuals and a p-value greater than 0.05 (chosen alpha level) indicates the residuals are normally distributed; the null hypothesis that the residuals are normally distributed is not rejected.}
#'  }
#'  \item{markov state figures}
#'    \itemize{
#'      \item{A) independent variable: precipitation}
#'      \item{B) dependent variable and states: streamflow observations, most likely state, and relative normal state estimate}
#'      \item{C) transformed dependent variable and states: transformed streamflow observations and most likely state}
#'      \item{D) conditional state probabilities for each state: probability of hydroState model remaining in given state}
#'    }
#'   }
#'
#'   These figures are often large, and below are a few common errors when the plotting window is too small. Exporting the plots as a file is recommend for the pseudo residual figure (\code{file = "file name"}).
#'  \itemize{
#'  \item{"Error in plot.new() : figure margins too large": reset plot window with "dev.off()", enlarge plot area and re-run \code{plot.residuals}.}
#'  \item{"Error in par(op) : invalid value specified for graphical parameter "pin" if the R plot window is not reset with "dev.off", an additional \code{plot.residuals} attempt will result in this error.}
#'  }
#'
#' @param x is the fitted hydroState model object.
#' @param pse.residuals option to plot pseudo residuals. Default is FALSE.
#' @param ind.variable option to plot independent variable over time.
#' @param dep.variable option to plot dependent variable and states over time.
#' @param dep.variable.transformed option to plot transformed dependent variable and states over time.
#' @param cond.state.prob option to plot the conditional state probabilities over time for each state.
#' @param siteID character string of catchment identifier (i.e. gauge ID). Default is NULL. Only recommended when exporting (i.e, file = "filename.pdf").
#' @param file character string of file directory/name to export plots as a pdf: "flow.state.plots.407211.pdf". Default is NULL, no pdf file is exported.
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
#' @importFrom grDevices dev.cur
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
#' ## Plot all markov state figures (default)
#' plot(model.annual.fitted.221201)
#'
#' ## Plot only dependent variable transformed with markov states
#' plot(model.annual.fitted.221201, dep.variable.transformed = TRUE)
#'


plot.hydroState <- function(x, ...,
                       pse.residuals = FALSE,
                       ind.variable = FALSE,
                       dep.variable = FALSE,
                       dep.variable.transformed = FALSE,
                       cond.state.prob = FALSE,
                       siteID = NULL,
                       file = NULL
                       ){

  model = x


  # if .pdf given in file continue else error
    if (!is.null(file)) {

      if(endsWith(file, ".pdf")){

        pdf(file, width = 8.5, height = 11)

        device_opened <- TRUE

      }else{
        stop("Error: file name is not a .pdf file")
      }
    } else {
      device_opened <- FALSE
    }



  # if site ID is NULL
  if (is.null(siteID)) {
    siteID = ""
  }

    if(pse.residuals == TRUE){

      temp.plot = check.PseudoResiduals(model, do.plot = T)


    }else{

      # default is to plot everything
      if(missing(ind.variable) & missing(dep.variable) & missing(dep.variable.transformed) & missing(cond.state.prob)){

        temp.plot = viterbi(model, do.plot = T, plot.options = c("A","B","C","D"))

      }else{

        plot.type = c()
        if(ind.variable) plot.type = c(plot.type, "A")
        if(dep.variable) plot.type = c(plot.type, "B")
        if(dep.variable.transformed) plot.type = c(plot.type, "C")
        if(cond.state.prob) plot.type = c(plot.type, "D")

        temp.plot = viterbi(model, do.plot = T, plot.options = plot.type)
      }
    }

  # if pdf do title
  if (nzchar(siteID) && dev.cur() > 1) {
    try(title(siteID), silent = TRUE)
  }

  # if pdf opened closed
  if (device_opened) dev.off()

  return(temp.plot)


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


      return(viterbi(model, do.plot = F, plot.options = c("A","B","C","D")))

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
#' \donttest{
#' check(model = model.annual.fitted.221201)
#'}


check <- function(model, n.samples = 100000){


  if(class(model)[1] == "hydroState"){

  return(check.viterbi(model, nSamples = n.samples))
  }else{
    stop('model is not a hydroState object')
  }

}



