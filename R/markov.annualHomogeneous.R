##' @include abstracts.R
##' @export
markov.annualHomogeneous <- setClass(
  # Set the name for the class
  "markov.annualHomogeneous",

  package='hydroState',

  contains=c('markov'),

  # Define the slots
  slots = c(
    transition.graph = "matrix",
    parameters = "list"

  ),

  # Set the default values for the slots. (optional)
  prototype=list(
    transition.graph = matrix(TRUE,2,2),
    parameters = list(mean=rep(1,2), std=rep(1,2),transition.prob = rep(0.5,2), initial.state.prob = 0.5)
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
            validObject(.Object)
            .Object
          }
)


# create a method to assign the parameter
#setGeneric(name="setParameters", def=function(.Object,parameters) {standardGeneric("setParameters")})
setMethod(f="setParameters",
          signature="markov.annualHomogeneous",
          definition=function(.Object,parameters)
          {
            .Object@parameters$mean = parameters$mean
            .Object@parameters$std  = parameters$std
            .Object@parameters$transition.prob = parameters$transition.prob
            .Object@parameters$initial.state.prob = parameters$initial.state.prob
            return(.Object)
          }
)
setMethod(f="setParameters.fromTransformed",
          signature="markov.annualHomogeneous",
          definition=function(.Object,parameters)
          {
            .Object@parameters$mean = parameters$mean
            .Object@parameters$std  = parameters$std
            .Object@parameters$transition.prob = parameters$transition.prob
            .Object@parameters$initial.state.prob = parameters$initial.state.prob
            if (any(is.na(.Object@parameters$transition.prob)))
              stop(print(parameters))
            return(.Object)
          }
)

# create a method to assign the parameter
#setGeneric(name="getParameters",def=function(.Object) {standardGeneric("getParameters")})
setMethod(f="getParameters",
          signature="markov.annualHomogeneous",
          definition=function(.Object)
          {
            return(list(mean=.Object@parameters$mean, std=.Object@parameters$std,transition.prob = .Object@parameters$transition.prob, initial.state.prob = .Object@parameters$initial.state.prob))
          }
)

setMethod(f="getTransformedParameterBounds",
          signature="markov.annualHomogeneous",
          definition=function(.Object)
          {
            parameters = getParameters(.Object)
            lowerBound = list(mean = rep(0.05, length(.Object@parameters$mean)),
                              std = rep(0.05, length(.Object@parameters$std)),
                              transition.prob = rep(0, length(.Object@parameters$transition.prob)),
                              initial.state.prob = rep(0, length(.Object@parameters$initial.state.prob))
                              )
            upperBound = list(mean = rep(2, length(.Object@parameters$mean)),
                              std = rep(2, length(.Object@parameters$std)),
                              transition.prob = rep(1, length(.Object@parameters$transition.prob)),
                              initial.state.prob = rep(1, length(.Object@parameters$initial.state.prob))
            )
            return(list(lower = lowerBound, upper = upperBound))
          }
)



# Get transition matrix with no input data.
#' @exportMethod getTransitionProbMatrix
setGeneric(name="getTransitionProbMatrix",def=function(.Object) {standardGeneric("getTransitionProbMatrix")})
setMethod(f="getTransitionProbMatrix",
          signature="markov.annualHomogeneous",
          definition=function(.Object)
          {
            # Get number of states
            nStates = getNumStates(.Object)

            # Build initial transition matrix
            Tprob = matrix(0.0,nStates,nStates)
            for (i in 1:nStates) {
              Tprob[i,i] =.Object@parameters$transition.prob[i]

              # Get index to state transitions from state i
              ind = which(.Object@transition.graph[i,])

              # Remove diagonal state from the indexes
              ind = ind[-i]

              # Chec there is only one other tranition and asign prob.
              if (length(ind)>1) {
                stop('transition.graph can only have one state to switch to.')
              }
              Tprob[i,ind] = 1-Tprob[i,i];
            }

            return(Tprob)
          }
)

# Get the log likelihood for the input data.
setMethod(f="getLikelihood", signature=c("markov.annualHomogeneous","data.frame"),
          definition=function(.Object, data)
          {

            # Remove NAs form the input data
            data = data[!is.na(data$Qhatbar),]

            # Get number of states
            nStates = getNumStates(.Object)

            #print(unlist(getParameters(.Object)))

            # Build the transition matrix.
            Tprob = getTransitionProbMatrix(.Object)
            #print(paste('...Tprob[1,1]=',Tprob[1,1]))

            # Get the Gaussian mean over time and for each state
            markov.means = getGaussianMean(.Object, data)
            #print(paste('...markov.means[1]=',markov.means[1]))

            # Get the Gaussian mean over time and for each state
            markov.stds = getGaussianStd(.Object, data)
            #print(paste('...markov.stds[1]=',markov.stds[1]))

            # Set initial state probabilities.
            alpha <- matrix(NA,1,2)
            alpha[1] <- .Object@parameters$initial.state.prob
            alpha[2] <- 1 - .Object@parameters$initial.state.prob;

            # Calculate a vector of time step size
            timesteps = getTimesteps(.Object, data)
            timesteps = timesteps$timesteps;


            # Loop through each time step and calculate alpha. At each time step
            # calculate the probability of Qhatbar for each state.
            for (i in 1:nrow(data)) {

              #print(paste('...i=',i))
              #print(paste('timesteps[i]=',timesteps[i]))
              if (is.na(data$Qhatbar[i]) || is.na(timesteps[i]))
                  next

                # Calculate thea probability matrix for the flow Qhatbar and the Gaussian mean and std
                P = diag(dnorm(data$Qhatbar[i],markov.means[i,],markov.stds[i,]))

                #print(paste('...i=',i,', P[1,1]=',P[1,1]))


                # Calc. alpha
                alpha = alpha %*% (Tprob^timesteps[i]) %*% P

                #print(paste('alpha=',alpha))
            }
            return( alpha %*%  matrix(1,nStates,1));
          }
)

# # Get the log lokihood for the input data.
# setMethod(f="getNegLogLikelihood",signature="markov.annualHomogeneous",definition=function(.Object, data)
#           {
#             return( -log(getLikelihood(.Object, data)) )
#           }
# )

# Get the gaussian mean.
#' @exportMethod getGaussianMean
setGeneric(name="getGaussianMean",def=function(.Object, data) {standardGeneric("getGaussianMean")})
setMethod(f="getGaussianMean",
          signature="markov.annualHomogeneous",
          definition=function(.Object, data)
          {
            if (is.numeric(data)) {
              nrows = length(data)
            } else if( is.data.frame(data)){
                nrows = nrow(data)
            }
            # Apply covariates to the Gaussian mean. This can be extdended to allow for
            # serial correlation, seasonality etc.
            return(matrix(.Object@parameters$mean * mean(data$Qhatbar), nrows, length(.Object@parameters$mean),byrow=T))
          }
)

# Get the gaussian standard deviation
#' @exportMethod getGaussianStd
setGeneric(name="getGaussianStd",def=function(.Object, data) {standardGeneric("getGaussianStd")})
setMethod(f="getGaussianStd",
          signature="markov.annualHomogeneous",
          definition=function(.Object, data)
          {

            if (is.numeric(data)) {
              nrows = length(data)
            } else if( is.data.frame(data)){
              nrows = nrow(data)
            }
            # Apply covariates to the Gaussian mean. This can be extdended to allow for
            # serial correlation, seasonality etc.
            return(matrix(.Object@parameters$std * sd(data$Qhatbar), nrows, length(.Object@parameters$std),byrow=T))
          }
)

# Get the number of states.
setGeneric(name="getNumStates",def=function(.Object, Qhatbar) {standardGeneric("getNumStates")})
setMethod(f="getNumStates",
          signature="markov.annualHomogeneous",
          definition=function(.Object)
          {
            return(length(.Object@parameters$mean))
          }
)

# Get the size of each time step AND the units.
setGeneric(name="getTimesteps",def=function(.Object, data) {standardGeneric("getTimesteps")})
setMethod(f="getTimesteps",
          signature=c("markov.annualHomogeneous","data.frame"),
          definition=function(.Object,data)
          {
            timestepType = '(none)'
            nRows = nrow(data)
            if (any(tolower(names(data))=='month')) {
              timestepType = 'monthly'
              timesteps = data$year[2:nRows] - data$year[1:(nRows-1)]
              yearlySteps =  data$year[1:nRows] - data$year[1] + 1;
              timesteps = c(1,yearlySteps*12 + data$month[2:nRows] -  data$month[1:nRows-1])

            } else if (any(tolower(names(data))=='year')) {
              timestepType = 'yearly'

              timesteps = c(1,data$year[2:nRows] - data$year[1:(nRows-1)]);

            } else {
              stop('Only monthly and yearly timesteps are supported. Please delete other time columns from the inpt data.')
            }

            return(list(timestepType=timestepType, timesteps=timesteps))
          }
)

