##' @include abstracts.R QhatModel.homo.normal.linear.R parameters.R
##' @export
QhatModel.homo.normal.linear.AR3 <- setClass(
  # Set the name for the class
  "QhatModel.homo.normal.linear.AR3",

  package='hydroState',

  contains=c('QhatModel.homo.normal.linear'),

  # Set the default values for the slots. (optional)
  prototype=list(
    parameters =  new('parameters',c('mean.a0', 'mean.a1','mean.AR1','mean.AR2','mean.AR3','std.a0'),c(1,1,1,1,1,1))
  )
)

# Valid object?
validObject <- function(object) {
  TRUE

}
setValidity("QhatModel.homo.normal.linear.AR3", validObject)

# Initialise object
#setGeneric(name="initialize",def=function(.Object,input.data){standardGeneric("initialize")})
setMethod("initialize","QhatModel.homo.normal.linear.AR3", function(.Object, input.data, use.truncated.dist=T, transition.graph=matrix(T,2,2),
                                                                    state.dependent.mean.a0=T, state.dependent.mean.a1=F,
                                                                    state.dependent.mean.AR1=F, state.dependent.mean.AR2=F,state.dependent.mean.AR3=F,
                                                                    state.dependent.std.a0=T) {
  .Object@input.data <- input.data
  .Object@use.truncated.dist = use.truncated.dist
  .Object@nStates = ncol(transition.graph)

  parameter.length <- as.numeric(c(state.dependent.mean.a0, state.dependent.mean.a1, state.dependent.mean.AR1, state.dependent.mean.AR2, state.dependent.mean.AR3, state.dependent.std.a0)) * (.Object@nStates-1) + 1

  # Set up model terms for mean and standard deviation.
  .Object@parameters = new('parameters', c('mean.a0', 'mean.a1', 'mean.AR1', 'mean.AR2','mean.AR3','std.a0'), parameter.length)

  validObject(.Object)
  .Object
}
)

setMethod(f="is.stationary",signature=c("QhatModel.homo.normal.linear.AR3"),definition=function(.Object)
{

  # Get object parameter list
  parameters = getParameters(.Object@parameters)

  # Find roots of AR quadratic equations.
  AR1 <- parameters$mean.AR1
  if (length(AR1)==1)
    AR1 = rep(AR1, .Object@nStates)

  AR2 <- parameters$mean.AR2
  if (length(AR2)==1)
    AR2 = rep(AR2, .Object@nStates)

  AR3 <- parameters$mean.AR3
  if (length(AR3)==1)
    AR3 = rep(AR3, .Object@nStates)

  AR.polynomial.coeffs = cbind(rep(1, .Object@nStates), -AR1, -AR2, -AR3)
  AR.polynomial.coeffs = unique(AR.polynomial.coeffs)

  # Find Roots and check if OUTSIDE the unit circle ie stable AR model
  for (i in 1:nrow(AR.polynomial.coeffs)) {
    # get roots
    AR.polynomial.roots = polyroot(as.vector(AR.polynomial.coeffs[i,]))

    # If unit circle, the return false for AR being non-stationary
    if (any(abs(AR.polynomial.roots)<=1))
      return(FALSE)
  }
  return(TRUE)

}
)

# Calculate the transformed flow at the mean annual precip
setMethod(f="getMean",signature=c("QhatModel.homo.normal.linear.AR3","data.frame"),definition=function(.Object, data) {
  return(getMean.AR3(.Object, data))
}
)

# Calculate the transformed flow at the mean annual precip
setGeneric(name="getMean.AR3",def=function(.Object, data) {standardGeneric("getMean.AR3")})
setMethod(f="getMean.AR3",signature=c("QhatModel.homo.normal.linear.AR3","data.frame"),definition=function(.Object, data)
          {
            # Find time points with finite PRECIPITATION  data
            filt <- is.finite(data$Qhat.precipitation)
            data <- data[filt,]

            # Get object parameter list
            parameters = getParameters(.Object@parameters)

            ncols.a1 = length(parameters$mean.a1)
            ncols.a0 = length(parameters$mean.a0)
            ncols.AR1 = length(parameters$mean.AR1)
            ncols.AR2 = length(parameters$mean.AR2)
            ncols.AR3 = length(parameters$mean.AR3)
            nrows = length(data$Qhat.precipitation);
            ncols.max = max(c(ncols.a0 ,ncols.a1,ncols.AR1,ncols.AR2,ncols.AR3))
            if (ncols.max > .Object@nStates)
                stop(paste('The number of parameters for each term of the mean model must must equal 1 or the number of states of ',.Object@nStates))


            # Check the AR1 model is stationary.
            if (!is.stationary(.Object)){
              Qhat.model <- matrix(Inf, length(filt), .Object@nStates)
              return(Qhat.model)
            }

            # Check which terms are uniform for all states and whic terms are unique
            # to each state.
            if (ncols.a0==1 || ncols.a0==.Object@nStates) {
                a0.est = matrix(rep(parameters$mean.a0,each=nrows),nrows,.Object@nStates);
            } else if (ncols.a0<.Object@nStates) {
              stop(paste('The number of parameters for the a0 term of the mean model must must equal 1 or the number of states of ',.Object@nStates))
            }
            if (ncols.a1==1 || ncols.a1==.Object@nStates) {
              a1.est = matrix(rep(parameters$mean.a1,each=nrows),nrows,.Object@nStates);
            } else if (ncols.a1<.Object@nStates) {
              stop(paste('The number of parameters for the a1 term of the mean model must must equal 1 or the number of states of ',.Object@nStates))
            }
            if (ncols.AR1==1 || ncols.AR1 ==.Object@nStates) {
              AR1.est = matrix(rep(parameters$mean.AR1,each=nrows),nrows,.Object@nStates);
            } else if (mean.AR1<.Object@nStates) {
              stop(paste('The number of parameters for the AR1 term of the mean model must must equal 1 or the number of states of ',.Object@nStates))
            }
            if (ncols.AR2==1 || ncols.AR2 ==.Object@nStates) {
              AR2.est = matrix(rep(parameters$mean.AR2,each=nrows),nrows,.Object@nStates);
            } else if (mean.AR2<.Object@nStates) {
              stop(paste('The number of parameters for the AR2 term of the mean model must must equal 1 or the number of states of ',.Object@nStates))
            }
            if (ncols.AR3==1 || ncols.AR3 ==.Object@nStates) {
              AR3.est = matrix(rep(parameters$mean.AR3,each=nrows),nrows,.Object@nStates);
            } else if (mean.AR3<.Object@nStates) {
              stop(paste('The number of parameters for the AR3 term of the mean model must must equal 1 or the number of states of ',.Object@nStates))
            }
            precip.data = matrix(data$Qhat.precipitation,nrows,.Object@nStates);

            # Calculate the mean with the AR1 componants
            a0.est <- 100 * a0.est
            Qhat.model <- precip.data * a1.est + a0.est
            Qhat.model[2,] <- Qhat.model[2,] + Qhat.model[1,] * AR1.est[2,]
            Qhat.model[3,] <- Qhat.model[3,] + Qhat.model[2,] * AR1.est[3,] + Qhat.model[1,] * AR2.est[3,]
            for (i in 4:nrows)
              Qhat.model[i,] <- Qhat.model[i,] + Qhat.model[i-1,] * AR1.est[i,] + Qhat.model[i-2,] * AR2.est[i,] + Qhat.model[i-3,] * AR3.est[i,]

            # NOTE, below is an attempt to improve run time perforamnce of the above loop. The results match from that
            # from the for loop but the run time perforance was very modest (7-14% faster). In the interest of
            # code transparency , the for loop was adopted. The idea comes from https://stackoverflow.com/questions/51312427/vectorizing-an-r-loop-with-backward-dependency
            # ind = 1:3
            # for (i in 1:.Object@nStates)
            #   Qhat.model[,i] <- .Call(stats:::C_rfilter, as.double(Qhat.model[,i]), as.double(c(AR1.est[1,i],AR2.est[1,i],AR3.est[1,i])),
            #                           c(rep(0,3), double(nrows)))[-ind]

            Qhat.model.NAs = matrix(NA,length(filt),.Object@nStates)
            Qhat.model.NAs[filt,] <- Qhat.model

            return(Qhat.model.NAs)
          }
)
