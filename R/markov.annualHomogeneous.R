##' @include abstracts.R parameters.R
## @export
markov.annualHomogeneous <- setClass(
  # Set the name for the class
  "markov.annualHomogeneous",

  package='hydroState',

  contains=c('markov'),

  # Define the slots
  slots = c(
    transition.graph = "matrix",
    transition.graph.parameter.index = "matrix",
    parameters = "parameters"

  ),

  # Set the default values for the slots. (optional)
  prototype=list(
    transition.graph = matrix(TRUE,2,2),
    transition.graph.parameter.index = matrix(c(1,-1,2,-1),nrow = 2,ncol=2),
    parameters = new('parameters',c('transition.prob', 'initial.state.prob'),c(1,1))
  )
)


# Valid object?
validObject <- function(object) {
  return(TRUE)
}
setValidity("markov.annualHomogeneous", validObject)

# Initialise object
#setGeneric(name="initialize",def=function(.Object,transition.graph) {standardGeneric("initialize")})
setMethod(f="initialize",
          signature="markov.annualHomogeneous",
          definition=function(.Object, transition.graph)
          {
            .Object@transition.graph <- transition.graph

            # Set te number of transition parameters and a matrxi of indexes to them.
            .Object <- initialize.TransitionParameters(.Object)

            # Set the initial stat probs.
            .Object <- initialize.StateProbParameters(.Object)

            validObject(.Object)
            return(.Object)
          }
)

# Set the number of transition parameters and how they relate to the transition graph.
setGeneric(name="initialize.TransitionParameters",def=function(.Object) {standardGeneric("initialize.TransitionParameters")})
setMethod(f="initialize.TransitionParameters",
          signature="markov.annualHomogeneous",
          definition=function(.Object)
          {

            # Get number of states
            nStates = getNumStates(.Object)

            # Handle 1 state model.
            if (nStates==1) {
              .Object@transition.graph.parameter.index <- matrix(1,1,1)
              .Object@parameters@values$transition.prob <- c()
              return(.Object)
            }

            # Set dimensions of matrix of indexes to transition probs.
            .Object@transition.graph.parameter.index = matrix(0,nStates,nStates)

            # First set the transition parameters to the number of states
            # and define indexes to these parameters.
            .Object@parameters@values$transition.prob = rep(0.5, nStates)
            diag(.Object@transition.graph.parameter.index) = 1:nStates

            # Loop through each row of the graph and assess if there are more than one TRUE
            # entry that is non-diagonal
            for (i in 1:nStates) {
              # Build an index to state transitions from state i excluding the diagonal
              ind <- 1:nStates
              filt <- as.vector(.Object@transition.graph[i,])
              filt[i] <- FALSE
              ind <- ind[filt]

              # Chec there is only one other tranition and asign prob.
              if (length(ind)>1) {
                nParams = length(.Object@parameters@values$transition.prob)
                .Object@parameters@values$transition.prob = c(.Object@parameters@values$transition.prob, rep(0.0, length(ind)-1))
                .Object@transition.graph.parameter.index[i,ind[1:length(ind)-1]] = (nParams+1):(nParams+length(ind)-1);
                .Object@transition.graph.parameter.index[i,ind[length(ind)]] = -1;
              } else {
                .Object@transition.graph.parameter.index[i,ind] = -1;
              }
            }

            return(.Object)
          }
)


setGeneric(name="initialize.StateProbParameters",def=function(.Object) {standardGeneric("initialize.StateProbParameters")})
setMethod(f="initialize.StateProbParameters",
          signature="markov.annualHomogeneous",
          definition=function(.Object)
          {
            # Get number of states
            nStates = getNumStates(.Object)

            # Handle 1 state model.
            if (nStates==1) {
              .Object@parameters@values$initial.state.prob <- c()
              return(.Object)
            }

            # Set dimensions of matrix of indexes to transition probs.
            .Object@parameters@values$initial.state.prob = rep(1/nStates,nStates-1)

            return(.Object)

          }
)

# create a method to assign the parameter
#setGeneric(name="setParameters", def=function(.Object,parameters) {standardGeneric("setParameters")})
setMethod(f="setParameters",
          signature="markov.annualHomogeneous",
          definition=function(.Object,parameters)
          {
            if (length(.Object@parameters$transition.prob)==length(parameters$transition.prob)) {
              .Object@parameters$transition.prob = pmin(pmax(0,parameters$transition.prob),1)
            } else {
              stop('The length of "transition.prob" must equal that within the object.')
            }

            if (length(.Object@parameters$initial.state.prob)==length(parameters$initial.state.prob)) {
              .Object@parameters$initial.state.prob = pmin(pmax(0,parameters$initial.state.prob),1)
            } else {
              stop('The length of "initial.state.prob" must equal that within the object.')
            }
            return(.Object)
          }
)

# Get a vector of initial state probabilies
setMethod(f="getInitialStateProbabilities",
          signature="markov.annualHomogeneous",
          definition=function(.Object)
          {

            # Get number of states
            nStates = getNumStates(.Object)

            # Handle 1 state model.
            if (nStates==1) {
              return(c(1))
            }

            # Get object parameter vector
            parameters = getParameters(.Object@parameters)

            initial.probs = c(parameters$initial.state.prob,0);

            # Limit the probs to between zero and one
            #initial.probs[initial.probs<0] <- 0
            #initial.probs[initial.probs>1] <- 1

            # Importantly, the initial probs must sum to 1.0 and all be between 0 and 1.
            # To achieve this, the first initial.state.prob parameter is taken as a
            # probability. The second parameter is not a probability. It is the fraction
            # of the remaining prob (i.e 1-initial.state.prob[1]) to be assigned to
            # initial prob[2] etc
            Psum <- initial.probs[1];
            if (nStates>2) {
              for (i in seq(2,nStates-1,by=1)) {
              #for (i in 2:(nStates-1)) {
                initial.probs[i] <- (1-Psum)*initial.probs[i]
                Psum <- Psum + initial.probs[i]
              }
            }
            initial.probs[nStates] = 1 - Psum


            return(initial.probs)

          }
)

# Get transition matrix with no input data.
setMethod(f="getTransitionProbabilities",
          signature="markov.annualHomogeneous",
          definition=function(.Object)
          {
            # Get number of states
            nStates = getNumStates(.Object)

            # Handle 1 state model.
            if (nStates==1) {
              return(c(1))
            }

            # Get object parameter vector
            parameters = getParameters(.Object@parameters)

            # Loop though each element of the transition graph to build the matrix
            Tprob = matrix(0.0,nStates,nStates)
            for (i in 1:nStates) {
              Psum <- 0
              k <- 0;
              for (j in 1:nStates) {
                if (.Object@transition.graph.parameter.index[i,j]>0) {
                  k <- k+1

                  if (k==1) {
                    Tprob[i,j] = parameters$transition.prob[.Object@transition.graph.parameter.index[i,j]]
                    Psum = Tprob[i,j]
                  } else {
                    Tprob[i,j] = (1-Psum) * parameters$transition.prob[.Object@transition.graph.parameter.index[i,j]]
                    Psum = Psum + Tprob[i,j]
                  }


                }
              }
            }

            # Calculate the derived probabilities (ie the one -1 value per row of transition.graph.parameter.index)
            for (i in 1:nStates) {
              ind = which(.Object@transition.graph.parameter.index[i,]==-1)
              Tprob[i,ind] = max(min(1, 1-sum(Tprob[i,])),0)
            }

            return(Tprob)
          }
)

# Get the log likelihood for the input data.
setMethod(f="getLogLikelihood", signature=c("markov.annualHomogeneous","data.frame","matrix"),
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

            # Get the transition matrix.
            Tprob = getTransitionProbabilities(.Object)

            # Only accept the transition probs. if the model persists in each state. This is estimated from an extenstion of the
            # reference below to >2 states. Specifically, if the sum of the probs of switching from state A to any oter state  and from
            # any other state to A is >=1, the then model DOES NOT persist in state A.
            # Martin F. Lambert, Julian P. Whiting, and Andrew V. Metcalfe, (2003) A non-parametric hidden Markov model for climate state
            # identification, Hydrology and Earth System Sciences, 7(5), 652667
            for (i in 1:nStates) {
              # Get prob. of switching out of state i
              ind = 1:nStates
              ind = ind[ind!=i]
              prob.out = sum(Tprob[i,ind])

              # Get prob of switching into state i
              prob.in = sum(Tprob[ind,i])

              if ((prob.in + prob.out)>=1) {
                return(Inf)
              }
            }

            # Only accept the parameters if the forward probabilities from the first to next time step show
            # that the state has not switched state after the first time step.
            P.forward <- getLogForwardProbabilities(.Object, data[1:2,], as.matrix(emission.probs[1:2,], ncol=nStates))
            if (all(!is.na(P.forward))) {
              if (any(which(max(P.forward[,1])==P.forward[,1]) != which(max(P.forward[,2])==P.forward[,2])))
                return(Inf)
            }


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

# Get the log forward probabilities for the input data.
# @exportMethod getLogForwardProbabilities
setGeneric(name="getLogForwardProbabilities",def=function(.Object, data, emission.probs) {standardGeneric("getLogForwardProbabilities")})
setMethod(f="getLogForwardProbabilities", signature=c("markov.annualHomogeneous","data.frame","matrix"),
          definition=function(.Object, data, emission.probs)
          {

            # Get number of states
            nStates = getNumStates(.Object)

            # Built filter for non NAs.
            filt <- !is.na(data$Qhat.flow)

            # Handle 1 state model.
            #if (nStates==1) {
            #  return(sum(log(emission.probs[filt])))
            #}

            # Get the initial state probabilites
            alpha = getInitialStateProbabilities(.Object)

            # Get the transition matrix.
            Tprob = getTransitionProbabilities(.Object)

            # Set Emmision probs 1.0 if no obs data.
            emission.probs[!filt,] = 1

            # Adapted from Zucchini, McDonald and Langrock, 2016, Hidden Markov Models for Time Series, 2nd Edision, Appendix A.
            #---------------------------------------------------------------------
            n             <- nrow(data)
            lalpha        <- matrix(NA,nStates,n)
            foo           <- alpha * as.vector(emission.probs[1,])
            sumfoo        <- sum(foo)
            lscale        <- log(sumfoo)
            foo           <- foo/sumfoo
            lalpha[,1]    <- lscale+log(foo)
            for (i in 2:n)
            {
              foo          <- foo %*% Tprob * as.vector(emission.probs[i,])
              sumfoo       <- sum(foo)
              lscale       <- lscale+log(sumfoo)
              foo          <- foo/sumfoo
              lalpha[,i]   <- log(foo)+lscale
            }
            return(lalpha)
            #---------------------------------------------------------------------
          }
)

# Get the log backward probabilities for the input data.
setGeneric(name="getLogBackwardProbabilities",def=function(.Object, data, emission.probs) {standardGeneric("getLogBackwardProbabilities")})
setMethod(f="getLogBackwardProbabilities", signature=c("markov.annualHomogeneous","data.frame","matrix"),
          definition=function(.Object, data, emission.probs)
          {

            # Get number of states
            nStates = getNumStates(.Object)

            # Built filter for non NAs.
            filt <- !is.na(data$Qhat.flow)

            # Handle 1 state model.
            #if (nStates==1) {
            #  return(sum(log(emission.probs[filt])))
            #}

            # Get the initial state probabilites
            alpha = getInitialStateProbabilities(.Object)

            # Get the transition matrix.
            Tprob = getTransitionProbabilities(.Object)

            # Set Emmision probs 1.0 if no obs data.
            emission.probs[!filt,] = 1

            # Adapted from Zucchini, McDonald and Langrock, 2016, Hidden Markov Models for Time Series, 2nd Edision, Appendix A.
            #---------------------------------------------------------------------
            n          <- nrow(data)
            m          <- nStates
            lbeta      <- matrix(NA,m,n)
            lbeta[,n]  <- rep(0,m)
            foo        <- rep(1/m,m)
            lscale     <- log(m)
            for (i in (n-1):1)
            {
              foo        <- Tprob %*% (as.vector(emission.probs[i+1,])*foo)
              lbeta[,i]  <- log(foo)+lscale
              sumfoo     <- sum(foo)
              foo        <- foo/sumfoo
              lscale     <- lscale+log(sumfoo)
            }
            return(lbeta)
            #---------------------------------------------------------------------
          }
)

# Get the conditional probability of a givn Qhat observation at time t given all other observations.
setMethod(f="getConditionalStateProbabilities", signature="markov.annualHomogeneous",
          definition=function(.Object, data, emission.probs)
          {
            # Get number of states
            nStates = getNumStates(.Object)

            # Adapted from Zucchini, McDonald and Langrock, 2016, Hidden Markov Models for Time Series, 2nd Edision, Appendix A.
            #---------------------------------------------------------------------
            n  <- nrow(data)
            m  <- nStates

            la <- getLogForwardProbabilities(.Object, data, emission.probs)
            lb <- getLogBackwardProbabilities(.Object, data, emission.probs)

            c <- max ( la [ , n ])
            llk <- c + log ( sum ( exp ( la [ , n ] - c ) ) )
            stateprobs <- matrix ( NA , ncol =n , nrow = m )
            for ( i in 1: n )
              stateprobs [ , i ] <- exp ( la [ , i ]+ lb [ , i ] - llk )

            return(stateprobs)

          }
)

# Get the conditional probability of a givn Qhat observation at time t given all other observations.
setMethod(f="getConditionalProbabilities", signature="markov.annualHomogeneous",
          definition=function(.Object, data, emission.probs, cumprob.atQhatIncrements)
          {
            # Get number of states
            nStates = getNumStates(.Object)

            # Get the initial state probabilites
            alpha = getInitialStateProbabilities(.Object)

            # Get the transition matrix.
            Tprob = getTransitionProbabilities(.Object)

            # Adapted from Zucchini, McDonald and Langrock, 2016, Hidden Markov Models for Time Series, 2nd Edision, Appendix A.
            #---------------------------------------------------------------------

            n          <- nrow(data)
            m          <- nStates
            nxc       <- dim(cumprob.atQhatIncrements)[3]
            dxc       <- matrix(NA,nrow=nxc,ncol=n)


            la        <- getLogForwardProbabilities(.Object, data, emission.probs)
            lb        <- getLogBackwardProbabilities(.Object, data, emission.probs)

            la        <- cbind(log(alpha),la)
            lafact    <- apply(la,2,max)
            lbfact    <- apply(lb,2,max)
            for (i in 1:n)
            {
              foo      <- (exp(la[,i]-lafact[i]) %*% Tprob)*exp(lb[,i]-lbfact[i])
              foo      <- foo/sum(foo)

              # Note, Zucchini, McDonald and Langrock, 2016 code calculates the probability of a given
              # Xt outside this for-loop because the emmision probs in their model are independent of
              # time. However, here the time-varying means and auto-regressive terms make the emmision probs
              # time-dependent and so here a 3D array is passed to this function with the depth dimension
              # containing the cumulative probs. at pre-defined incrementts of Qhat for a given state and time point.
              Px <- cumprob.atQhatIncrements[i,,]

              dxc[,i]  <- as.vector(foo%*%Px)

            }

            return(dxc)
            #---------------------------------------------------------------------
          }
)

# Get the number of states.
setGeneric(name="getNumStates",def=function(.Object) {standardGeneric("getNumStates")})
setMethod(f="getNumStates",
          signature="markov.annualHomogeneous",
          definition=function(.Object)
          {
            return(ncol(.Object@transition.graph))
          }
)

setMethod(f="generate.sample.states",signature="markov.annualHomogeneous",definition=function(.Object, data)
{

  # Generate a synthtic series from the HMM states.
  # Note, the sampling uses the transition prob. matrix to estimate the long
  # term prob. of being in each state.
  nStates <- getNumStates(.Object)
  nSamples = nrow(data)
  states = rep(0,nSamples)

  # Get the initial state probabilites
  alpha = getInitialStateProbabilities(.Object)

  # Get the transition matrix.
  Tprob = getTransitionProbabilities(.Object)

  states[1] <- sample(1:nStates,1,prob=alpha)
  for (i in 2:nSamples) {
    alpha <- alpha %*% Tprob
    states[i] = sample(1:nStates,1,prob=alpha)
  }

  return(states)
}
)


