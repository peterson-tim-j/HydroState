##' @include abstracts.R
##' @import methods
## @export
model_parameters <- setClass(
  # Set the name for the class
  "model_parameters",

  package='hydroState',

  slots = c(
    values = 'list',
    bounds = 'matrix',
    transform = 'logical'
  ),

  # Set the default values for the slots. (optional)
  prototype=list(
    values = list(),
    bounds = matrix( c(-Inf, +Inf), 1, 2, dimnames = list('Param.name',c('lower','upper'))),
    transform = F
  )
)

# # Valid object?
# validObject <- function(object) {
#   TRUE
# }
# setValidity("parameters", validObject)

# Initialise object
#setGeneric(name="initialize",def=function(.Object,input.data){standardGeneric("initialize")})
setMethod("initialize","model_parameters", function(.Object, parameter.name, parameter.length) {

  # # Return empty obj if no parameters
  # if (length(parameter.name) ==0 && length(parameter.length)==0) {
  #   .Object@values <- numeric()
  #   .Object@bounds <- numeric()
  #   .Object@transform <- logical()
  # }


  if (length(parameter.length) != length(parameter.name))
    stop('The input "parameter.length" must be a integer vector of equal length to "parameter.name".')

  .Object@values <- vector('list',length(parameter.name))
  names(.Object@values) <- parameter.name

  .Object@transform <- rep(F,length(parameter.name))
  names(.Object@transform) <- parameter.name


  .Object@bounds <-  matrix( 0.0,
                             length(parameter.name), 2,
                             dimnames = list(parameter.name, c('lower','upper')))
  .Object@bounds[,2] <- 1.0

  if (length(parameter.length)>0) {
    for (i in 1:length(parameter.name)) {
      initial.val <- .Object@bounds[i,1] + 0.5*(.Object@bounds[i,2] - .Object@bounds[i,1])
      .Object@values[[i]] = rep(initial.val, parameter.length[i])
    }
  }


  return(.Object)
}
)

# Set/get the NON-TRANSFORMED parameter bounds.
setMethod(f="bounds<-",
          signature=c("model_parameters", value='matrix'),
          definition=function(.Object, value)
          {
            if (nrow(value)==0) {
              .Object@bounds  <- matrix( numeric(), 0, 2, dimnames = list(character(),c('lower','upper')))
              return(.Object)
            }

            parameter.names <- names(.Object@values)

            for (i in 1:length(parameter.names)) {
                .Object@bounds[parameter.names[i], 1] <- value[parameter.names[i], 1]
                .Object@bounds[parameter.names[i], 2] <- value[parameter.names[i], 2]
            }
            return(.Object)
          }
)
setMethod(f="bounds",
          signature=c(.Object="model_parameters", transformed='missing'),
          definition=function(.Object, transformed)
          {
            return(.Object@bounds)
          }
)
setMethod(f="bounds",
          signature=c(.Object="model_parameters",transformed='logical'),
          definition=function(.Object, transformed)
          {
            bnds = .Object@bounds
            if (nrow(bnds)==0)
              return(bnds)

            if (transformed) {
              parameter.names <- names(.Object@values)
              for (i in 1:nrow(bnds)) {
                if (.Object@transform[[parameter.names[i]]])
                  bnds[parameter.names[i], ] <- log10(bnds[parameter.names[i], ])
              }
            }
            return(bnds)
          }
)

setMethod(f="transforms", signature="model_parameters",
          definition=function(.Object)
          {
            return(.Object@transform)
          }
)
setMethod(f="transforms<-", signature = c(.Object="model_parameters", value='logical'),
          definition=function(.Object, value)
          {
            if (length(value)==0)
              return(.Object)

            parameter.names <- names(.Object@values)

            if (length(parameter.names)==0) {
              .Object@transform = logical()
              return(.Object)
            }

            for (i in 1:length(parameter.names))
              .Object@transform[[parameter.names[i]]] <- value[[parameter.names[i]]]

            return(.Object)
          }
)

# Get/Set parameters in real space from a list obj
setMethod(f="parameters", signature = c(.Object="model_parameters", transformed='missing', as_df='missing'),
          definition=function(.Object, transformed, as_df)
          {
            return(.Object@values)
          }
)
setMethod(f="parameters", signature = c(.Object="model_parameters", transformed='logical', as_df='missing'),
          definition=function(.Object, transformed, as_df)
          {
            parameters <- parameters(.Object)
            parameter.name <- names(parameters)
            if (length(parameter.name)==0)
              return(list())

            if (transformed)
              for (i in 1:length(parameter.name))
                if (.Object@transform[[parameter.name[i]]])
                  parameters[[parameter.name[i]]] <- log10(parameters[[parameter.name[i]]])

            return(parameters)
          }
)
setMethod(f="parameters",
          signature = c(.Object="model_parameters", transformed='logical', as_df='logical'),
          definition=function(.Object, transformed, as_df)
          {
            parameters <- parameters(.Object, transformed)
            if (as_df) {
              vals = data.frame(parameter = character(), state = numeric(), value = numeric())
              if (length(parameters) == 0)
                return(vals)

              parameter.name <- names(parameters)
              nstates = max(unlist(lapply(parameters, length)))
              nparams = sum(unlist(lapply(parameters, length)))
              vals = data.frame(parameter = character(), state = numeric(), value = numeric())
              nrow = 0
              for (i in 1:length(parameter.name)) {
                if ( length(parameters[[parameter.name[i]]]) == 1) {
                  nrow = nrow + 1
                  vals[nrow, 'parameter'] = parameter.name[i]
                  vals[nrow, 'state'] = 0
                  vals[nrow, 'value'] = parameters[[ parameter.name[i] ]]
                }
                else {
                  for (j in 1:length( parameters[[ parameter.name[i] ]])) {
                    nrow = nrow + 1
                    vals[nrow, 'parameter'] = parameter.name[i]
                    vals[nrow, 'state'] = j
                    vals[nrow, 'value'] = parameters[[ parameter.name[i] ]][j]
                  }
                }
              }
            }
            return(vals)
          }
)
setMethod(f="parameters<-",
          signature = c(.Object="model_parameters", transformed='missing', value='list'),
          definition=function(.Object, transformed, value)
          {
            if (length(value)==0)
              return(.Object)

            parameter.name <- names(.Object@values)
            for (i in 1:length(parameter.name))
              .Object@values[[parameter.name[i]]] <- value[[parameter.name[i]]]

            return(.Object)
          }
)
setMethod(f="parameters<-",
          signature = c(.Object="model_parameters", transformed='logical', value='list'),
          definition=function(.Object, transformed, value)
          {
            if (length(value)==0)
              return(.Object)

            if (transformed) {
              parameter.name <- names(.Object@values)
              for (i in 1:length(parameter.name))
                if (.Object@transform[parameter.name[i]])
                  value[[parameter.name[i]]] <- 10^value[[parameter.name[i]]]
            }
            parameters(.Object) <- value
            return(.Object)
          }
)

