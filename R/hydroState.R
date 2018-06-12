##' @include abstracts.R
##' @include Qhat.boxcox.R
##' @include Qhatbar.linearRegression.R
##' @include markov.annualHomogeneous.R
##' @export
hydroState <- setClass(
  # Set the name for the class
  "hydroState",

  package='hydroState',

  # Define the slots
  slots = c(
    input.data = "data.frame",
    Qhat.object   = "Qhat",
    Qhatbar.object   = "Qhatbar",
    markov.model.object = "markov",
    calibration.results = "list"
  ),

  # Set the default values for the slots. (optional)
  prototype=list(
    input.data = data.frame(year=c(0),month=c(0),precipitation=c(0),flow=c(0)),
    Qhat.object  = new('Qhat.boxcox',input.data = data.frame(year=c(0),month=c(0),precipitation=c(0),flow=c(0))),
    Qhatbar.object = new('Qhatbar.linearRegression',input.data = data.frame(year=c(0),month=c(0),precipitation=c(0),flow=c(0))),
    markov.model.object = new('markov.annualHomogeneous',transition.graph = matrix(TRUE,2,2)),
    calibration.results = list(value=0,par=0,gradients=0,generations=0,peakgeneration=0,popsize=0,operators=0)
  )
)

# Valid object?
validObject <- function(object) {
  TRUE
}
setValidity("hydroState", validObject)

# Initialise the object.
#setGeneric(name="initialize",def=function(.Object,input.data, Qhat.object, Qhatbar.object, markov.model.object, ...){standardGeneric("initialize")})
setMethod(f="initialize",signature="hydroState",definition=function(.Object, input.data, Qhat.object, Qhatbar.object, markov.model.object)
          {
          .Object@input.data = input.data
          .Object@Qhat.object = Qhat.object
          .Object@Qhatbar.object = Qhatbar.object
          .Object@markov.model.object = markov.model.object
          .Object
          }
)

# Build rhe following methods:
#  - calibrate
#  - maximum.liklihood (it should update Qhat and Qhatbar and then ...)
#  - global.decode
#  - nameStates (hidden? called by calibrate?)
#  - forecast?
#  - plot
#  -

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
            list(
              Qhat = getParameters(.Object@Qhat.object),
              Qhatbar = getParameters(.Object@Qhatbar.object),
              markov = getParameters(.Object@markov.model.object)
            )
          }
)

# Get the full set of model parameters as a list.
#' @exportMethod getParameters.asVector
setGeneric(name="getParameters.asVector",def=function(.Object) {standardGeneric("getParameters.asVector")})
setMethod(f="getParameters.asVector",
          signature="hydroState",
          definition=function(.Object)
          {
            parameters.all = getParameters(.Object);

            return( c(
              unlist(parameters.all$Qhat),
              unlist(parameters.all$Qhatbar),
              unlist(parameters.all$markov)
              ))
          }
)

#' @exportMethod getTransformedParameterBounds.asVector
setGeneric(name="getTransformedParameterBounds.asVector",def=function(.Object) {standardGeneric("getTransformedParameterBounds.asVector")})
setMethod(f="getTransformedParameterBounds.asVector",
          signature="hydroState",
          definition=function(.Object)
          {
            # Get the parameter bounds for each object
            Qhat = getTransformedParameterBounds(.Object@Qhat.object)
            Qhatbar = getTransformedParameterBounds(.Object@Qhatbar.object)
            markov = getTransformedParameterBounds(.Object@markov.model.object)

            lowerBound = c( unlist(Qhat$lower),
                            unlist(Qhatbar$lower),
                            unlist(markov$lower))

            upperBound = c( unlist(Qhat$upper),
                            unlist(Qhatbar$upper),
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
              .Object@Qhat.object <- setParameters(.Object@Qhat.object, parameters$Qhat)
              .Object@Qhatbar.object <- setParameters(.Object@Qhatbar.object, parameters$Qhatbar)
              .Object@markov.model.object <- setParameters(.Object@markov.model.object, parameters$markov)
            } else if(is.numeric(parameters)) {
              .Object@Qhat.object <- setParameters.fromVector(.Object@Qhat.object, parameters$Qhat)
              .Object@Qhatbar.object <- setParameters.fromVector(.Object@Qhatbar.object, parameters$Qhatbar)
              .Object@markov.model.object <- setParameters.fromVector(.Object@markov.model.object, parameters$markov)

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
            Qhat = getParameters(.Object@Qhat.object)
            Qhatbar = getParameters(.Object@Qhatbar.object)
            markov = getParameters(.Object@markov.model.object)

            from =1;
            to = length(Qhat);
            parameters.aslist = relist(parameters[from:to],Qhat)
            .Object@Qhat.object < -setParameters(.Object@Qhat.object, parameters.aslist)

            from =1+length(Qhat);
            to = length(Qhat)+length(Qhatbar);
            parameters.aslist = relist(parameters[from:to],Qhatbar)
            .Object@Qhatbar.object <- setParameters(.Object@Qhatbar.object, parameters.aslist)

            from =1+length(Qhat)+length(Qhatbar);
            to = length(Qhat)+length(Qhatbar)+length(markov);
            parameters.aslist = relist(parameters[from:to],markov)
            .Object@markov.model.object <- setParameters(.Object@markov.model.object, parameters.aslist)

            return(.Object)

          }
)

#' @exportMethod setParameters.fromTransformedVector
setGeneric(name="setParameters.fromTransformedVector",def=function(.Object,parameters) {standardGeneric("setParameters.fromTransformedVector")})
setMethod(f="setParameters.fromTransformedVector",
          signature=c("hydroState","numeric"),
          definition=function(.Object,parameters)
          {
            # Get the parameter structure for each object
            Qhat = getParameters(.Object@Qhat.object)
            Qhatbar = getParameters(.Object@Qhatbar.object)
            markov = getParameters(.Object@markov.model.object)

            # Calc. number parameters per object
            nQhat <- sum(lengths(Qhat))
            nQhatbar <- sum(lengths(Qhatbar))
            nMarkov <- sum(lengths(markov))

            from =1;
            to = nQhat;
            parameters.aslist = relist(parameters[from:to],Qhat)
            .Object@Qhat.object <- setParameters.fromTransformed(.Object@Qhat.object, parameters.aslist)

            from =1+nQhat;
            to = nQhat+nQhatbar;
            parameters.aslist = relist(parameters[from:to],Qhatbar)
            .Object@Qhatbar.object <- setParameters.fromTransformed(.Object@Qhatbar.object, parameters.aslist)

            from =1+nQhat+nQhatbar;
            to = nQhat+nQhatbar+nMarkov;
            parameters.aslist = relist(parameters[from:to],markov)
            .Object@markov.model.object <- setParameters.fromTransformed(.Object@markov.model.object, parameters.aslist)

            return(.Object)
          }
)

# Get the model negative log liklihood.
# #' @exportMethod getNegLogLikelihood
setGeneric(name="getNegLogLikelihood",def=function(.Object, parameters) {standardGeneric("getNegLogLikelihood")})
setMethod(f="getNegLogLikelihood",signature=c(.Object="hydroState"),definition=function(.Object,parameters)
          {
            # Set the parameters.
           .Object <- setParameters(.Object, parameters)

            return(getNegLogLikelihood(.Object))
          }
          )
#setGeneric(name="getNegLogLikelihood",def=function(.Object) {standardGeneric("getNegLogLikelihood")})
setMethod(f="getNegLogLikelihood",signature=c(.Object="hydroState",parameters='missing'),definition=function(.Object)
          {
            # Get the transformed flow, Qhat
            Qhat = getQhat(.Object@Qhat.object, .Object@input.data)

            # Get the transformed flow at the mean precip, Qhatbar.
            Qhatbar = getQhatbar(.Object@Qhatbar.object, Qhat)

            # Add Qhatbar to the input data
            data = cbind.data.frame(.Object@input.data,Qhatbar)

            # Get the markov likelihod and return
            return( -log(getLikelihood(.Object@markov.model.object, data)))
          }
)

#' @exportMethod getNegLogLikelihood.fromTransformedVector
setGeneric(name="getNegLogLikelihood.fromTransformedVector",def=function(parameters,.Object) {standardGeneric("getNegLogLikelihood.fromTransformedVector")})
setMethod(f="getNegLogLikelihood.fromTransformedVector",signature=c(parameters="numeric",.Object="hydroState"),definition=function(parameters,.Object)
          {
            # Set the parameters.
            .Object <- setParameters.fromTransformedVector(.Object, parameters)

            return(getNegLogLikelihood(.Object))
          }
)

#' @exportMethod fit
setGeneric(name="fit",def=function(.Object, pop.size.perParameter=NA,
                                   max.generations=NA,
                                   wait.generations=NA,
                                   Domains=NA,
                                   solution.tolerance=NA,
                                   boundary.enforcement=NA,
                                   print.level=NA,
                                   cluster=NA, ...) {standardGeneric("fit")})
setMethod(f = "fit",signature="hydroState",definition=function(.Object,
                                                             pop.size.perParameter=25,
                                                             max.generations=10000,
                                                             wait.generations=10,
                                                             Domains,
                                                             solution.tolerance=1e-6,
                                                             boundary.enforcement=2,
                                                             print.level=1,
                                                             cluster = "localhost", ...)
  {

    # Setup rgenoud
    Domains <- getTransformedParameterBounds.asVector(.Object)
    nvars <- nrow(Domains)

    # Run optimiser and assign solution to the object.
    calib.results <- genoud(hydroState::getNegLogLikelihood.fromTransformedVector,
                             nvars=nvars,
                             pop.size=pop.size.perParameter*nvars,
                             max.generations=max.generations,
                             wait.generations=wait.generations,
                             hard.generation.limit =T,
                             Domains = Domains,
                             solution.tolerance=solution.tolerance,
                             boundary.enforcement=boundary.enforcement,
                             print.level=print.level,
                             BFGSburnin = -1,
                             optim.method = "Nelder-Mead",
                             .Object=.Object, ...)

    .Object@calibration.results <- calib.results

    # Set calibrated parameters
    .Object <- setParameters.fromTransformedVector(.Object, .Object@calibration.results$par)

    return(.Object)
  }
)

#' @exportMethod setStateNames
setGeneric(name="setStateNames",def=function(.Object, ...) {standardGeneric("setStateNames")})
setMethod(f="setStateNames",signature="hydroState",definition=function(.Object, year.normalFlow)
  {


  }
)

#' @exportMethod plot.fit
setGeneric(name="plot.fit",def=function(.Object) {standardGeneric("plot.fit")})
setMethod(f="plot.fit",signature="hydroState",definition=function(.Object)
  {


  }
)

#' @exportMethod viterbi
setGeneric(name="viterbi",def=function(.Object) {standardGeneric("viterbi")})
setMethod(f="viterbi",signature="hydroState",definition=function(.Object)
  {


  }
)

#' @exportMethod plot.viterbi
setGeneric(name="plot.viterbi",def=function(.Object) {standardGeneric("plot.viterbi")})
setMethod(f="plot.viterbi",signature="hydroState",definition=function(.Object)
  {


  }
)

#' @exportMethod forecast
setGeneric(name="forecast",def=function(.Object, t) {standardGeneric("forecast")})
setMethod(f="forecast",signature="hydroState",definition=function(.Object, t)
  {


  }
)
