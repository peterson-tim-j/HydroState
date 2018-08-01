##' @include abstracts.R parameters.R markov.annualHomogeneous.R
##' @export
markov.annualHomogeneous.flickering <- setClass(
  # Set the name for the class
  "markov.annualHomogeneous.flickering",

  package='hydroState',

  contains=c('markov.annualHomogeneous'),

  # Define the slots
  slots = c(
    allow.flickering = "logical"
  ),

  # Set the default values for the slots. (optional)
  prototype=list(
    allow.flickering = T
  )
)


# Valid object?
validObject <- function(object) {
  return(TRUE)
}
setValidity("markov.annualHomogeneous.flickering", validObject)

# Initialise object
#setGeneric(name="initialize",def=function(.Object,transition.graph) {standardGeneric("initialize")})
setMethod(f="initialize",
          signature="markov.annualHomogeneous.flickering",
          definition=function(.Object, transition.graph, allow.flickering=T)
          {
            .Object@transition.graph <- transition.graph

            # Set te number of transition parameters and a matrxi of indexes to them.
            .Object <- initialize.TransitionParameters(.Object)

            # Set the initial stat probs.
            .Object <- initialize.StateProbParameters(.Object)


            # Check check.state.flickering
            if (is.logical(allow.flickering) && length(allow.flickering)==1) {
              .Object@allow.flickering = allow.flickering
            } else {
              stop('allow.flickering must be a logical variable.')
            }


            validObject(.Object)
            return(.Object)
          }
)

# Get the flickering for each state.
#' @exportMethod getStateFlicker
setGeneric(name="getStateFlicker",def=function(.Object) {standardGeneric("getStateFlicker")})
setMethod(f="getStateFlicker",signature=c("markov.annualHomogeneous.flickering"),definition=function(.Object)
{
  # Get the transition matrix.
  Tprob = getTransitionProbabilities(.Object)

  # Get number of states
  nStates = getNumStates(.Object)

  # Handle one state.
  if (nStates==1)
    return(0)

  # Estimate the flickering from an extenstion of the reference below to >2 states.
  # Specifically, if the sum of the probs of switching from state A to any oter state  and from
  # any other state to A is >=1, the then model DOES NOT persist in state A.
  # Martin F. Lambert, Julian P. Whiting, and Andrew V. Metcalfe, (2003) A non-parametric hidden Markov model for climate state
  # identification, Hydrology and Earth System Sciences, 7(5), 652667
  flickerValue = rep(0,nStates)
  for (i in 1:nStates) {
    # Get prob. of switching out of state i
    ind = 1:nStates
    ind = ind[ind!=i]
    prob.out = sum(Tprob[i,ind])

    # Get prob of switching into state i
    prob.in = sum(Tprob[ind,i])

    flickerValue[i] = (prob.in + prob.out) - 1
  }

  return(flickerValue)

}
)

# Get the log likelihood for the input data.
setMethod(f="getLogLikelihood", signature=c("markov.annualHomogeneous.flickering","data.frame","matrix"),
          definition=function(.Object, data, emission.probs)
          {
          # Check all the emmision probs. are finite.
            if (any(is.infinite(emission.probs)))
              return(Inf)

            # Get number of states
            nStates = getNumStates(.Object)

            # Check required fields exist
            if (!any(names(data)=='Qhat.flow'))
              stop('"data" must be a data.frame with the field "Qhat.flow".')

            # Built filter for non NAs.
            filt <- is.finite(data$Qhat.flow)

            # Handle 1 state model.
            if (nStates==1) {
              return(sum(log(emission.probs[filt])))
            }

            # Get the initial state probabilites
            alpha = getInitialStateProbabilities(.Object)

            # Check the initial state probs are sum to 1 and are all b/w 0 and 1.
            if (abs(sum(alpha)-1)>sqrt(.Machine$double.eps) || any(alpha<0) || any(alpha>1))
              return(Inf)

            # Only accept the transition probs. if the model persists in each state.
            # This is defined as any state having a flicker value >1.
            if (!.Object@allow.flickering && any(getStateFlicker(.Object)>1))
              return(Inf)

            # Only accept the parameters if the forward probabilities from the first to next time step show
            # that the state has not switched state after the first time step.
            P.forward <- getLogForwardProbabilities(.Object, data[1:2,], as.matrix(emission.probs[1:2,], ncol=nStates))
            if (all(!is.na(P.forward))) {
              if (which(max(P.forward[,1])==P.forward[,1]) != which(max(P.forward[,2])==P.forward[,2]))
                return(Inf)
            }

            # Get the transition matrix.
            Tprob = getTransitionProbabilities(.Object)

            # Set Emmision probs 1.0 if no obs data.
            emission.probs[!filt,] = 1

            # Adapted from Zucchini, McDonald and Langrock, 2016, Hidden Markov Models for Time Series, 2nd Edision, Appendix A.
            #---------------------------------------------------------------------

            # Get alpha at first time step and log liklihood
            alpha      <- alpha * emission.probs[1,]
            sumalpha   <- sum(alpha)
            lscale   <- log(sumalpha)
            alpha      <- alpha/sumalpha

            # Loop through 2+ time steps
            for (i in 2:nrow(data)) {
              alpha    <- alpha %*% Tprob * as.vector(emission.probs[i,])
              sumalpha <- sum(alpha)
              lscale <- lscale+log(sumalpha)
              alpha    <- alpha/sumalpha

              # if (!is.finite(lscale))
              #   stop('DBG')
            }

            return(lscale)

          }
)

