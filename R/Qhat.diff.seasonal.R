##' @include abstracts.R parameters.R
##' @export
Qhat.diff.seasonal <- setClass(
  # Set the name for the class
  "Qhat.diff.seasonal",

  package='hydroState',

  contains=c('Qhat'),

  # Define the slots
  slots = c(
    input.data = "data.frame",
    parameters = 'parameters'
  ),

  # Set the default values for the slots. (optional)
  prototype=list(
    input.data = data.frame(year=c(0),month=c(0),precipitation=c(0)),
    parameters= new('parameters',c(),c())
  )


)

# Initialise object
#setGeneric(name="initialize",def=function(.Object,input.data){standardGeneric("initialize")})
setMethod("initialize","Qhat.diff.seasonal", function(.Object, input.data) {
  .Object@input.data <- input.data
  validObject(.Object)
  .Object
}
)

# Calculate the transformed flow
#setGeneric(name="getQhat",def=function(.Object, data){standardGeneric("getQhat")})
setMethod(f="getQhat",signature=c("Qhat.diff.seasonal",'data.frame'),definition=function(.Object, data)
          {

            if (!is.data.frame(data))
              stop('"data" must be a data.frame.')

            nTimeSteps = 2;
            diff.est <-c(rep(NA,nTimeSteps), diff(data$flow,lag=nTimeSteps))
            diff.est[is.nan(diff.est)] <- NA
            data$Qhat.flow <-diff.est

            diff.est <- c(rep(NA,nTimeSteps), diff(data$precipitation,lag=nTimeSteps))
            diff.est[is.nan(diff.est)] <- NA
            data$Qhat.precipitation <- diff.est
            return(data)
          }
)

# Calculate the transformed flow using the object data
#setGeneric(name="getQhat",def=function(.Object){standardGeneric("getQhat")})
setMethod(f="getQhat",signature="Qhat.diff.seasonal",definition=function(.Object)
          {
             data = .Object@input.data
             return(getQhat(.Object, data))
          }
)


setMethod(f="getQ.backTransformed",signature=c("Qhat.diff.seasonal",'data.frame'),definition=function(.Object, data)
{

  if (!is.data.frame(data))
    stop('"Data" must be a data.frame.')

  # Steps for back transforming integrated flow to flow (note, integrated is as for a ARIMA model, not as in for calculus)
  # Also, the input data$Qhat will be the difference in the flow over 12 months.
  # 1. Find data$flow time steps with a finite value AND with Qhat == NA. These points are the start of a differencing timestep.
  # 2. For each from #1, calc. the cumsum and then add the obs. fow at the start of the period.
  # ind.all = 1:nrow(data)

  # # Find index to rows where both the flow and Qhat are NA.
  # ind.all.na <- which(is.na(data$flow) & is.na(data$Qhat.flow))

  # Find data$flow time steps with a finite value AND with Qhat == NA. These points are the start of a differencing timestep.
  ind_start.of.diff.period <- which(is.finite(data$flow) & is.na(data$Qhat.flow))

  # Remove indexes where there is only one time step between indexes
  #filt = diff(ind_start.of.diff.period, lag=1)>1
  #ind_start.of.diff.period = ind_start.of.diff.period[filt] + 1

  # Initialsie output
  data$flow.modelled <- rep(NA,nrow(data))

  # Loop through each differencing period
  ndiff <- 2
  for (ind in ind_start.of.diff.period) {

    j=ind+ndiff
    data$flow.modelled[ind] <- data$flow[ind]
    while(!is.na(data$Qhat.flow[j])) {
      data$flow.modelled[j] <- data$flow.modelled[j-ndiff] + data$Qhat.flow[j]
      j=j+ndiff
    }

    # # Find index to the end of the differencing period.
    # ind.end <- ind.all.na[ind.all.na > ind]
    # if (length(ind.end)==0) {
    #   ind.end <- nrow(data)
    # } else {
    #   ind.end <- ind.end[1]-1
    # }

    # stop("DBG")

    # # Calc. the cumsum over the differencing period.
    # Qhat.sum = cumsum(data$Qhat.flow[(ind+1):ind.end])

    # Add the obs. fow at the start of the period to get the back trsanformed flow
    # flow.est <- c( data$flow[ind], data$Qhat.flow[(ind+1):ind.end] + data$flow[ind])
    #
    # data$flow.modelled[ind:ind.end] <- flow.est

  }

  return(data)
}
)
