##' @include abstracts.R
##' @export
parameters <- setClass(
  # Set the name for the class
  "parameters",

  package='hydroState',

  slots = c(
    values = 'list',
    lower.bound = 'numeric',
    upper.bound = 'numeric',
    use.log.transform = 'logical'
  ),

  # Set the default values for the slots. (optional)
  prototype=list(
    values = list(),
    lower.bound = -Inf,
    upper.bound = Inf,
    use.log.transform = F
  )
)

# Valid object?
validObject <- function(object) {
  TRUE

}
setValidity("parameters", validObject)

# Initialise object
#setGeneric(name="initialize",def=function(.Object,input.data){standardGeneric("initialize")})
setMethod("initialize","parameters", function(.Object, parameter.name, parameter.length) {

  if (length(parameter.length) != length(parameter.name))
    stop('The input "parameter.length" must be a integer vector of equal length to "parameter.name".')

  .Object@values <- vector('list',length(parameter.name))
  names(.Object@values) <- parameter.name

  .Object@use.log.transform <- rep(F,length(parameter.name))
  names(.Object@use.log.transform) <- parameter.name
  .Object <- setTransforms(.Object)

  .Object@lower.bound <- rep(-Inf,length(parameter.name))
  .Object@upper.bound <- rep(Inf,length(parameter.name))
  names(.Object@lower.bound) <- parameter.name
  names(.Object@upper.bound) <- parameter.name
  .Object <- setBounds(.Object)

  if (length(parameter.length)>0) {
    for (i in 1:length(parameter.name)) {
      initial.val <- .Object@lower.bound[[i]] + 0.5*(.Object@upper.bound[[i]] - .Object@lower.bound[[i]])
      .Object@values[[i]] = rep(initial.val, parameter.length[i])
    }
  }
  return(.Object)
}
)

# Hard code bounds for each parameter. Importantly, the values are in NON-TRANSFORMED space.
setGeneric(name="setBounds",def=function(.Object) {standardGeneric("setBounds")})
setMethod(f="setBounds",
          signature="parameters",
          definition=function(.Object)
          {

            parameter.names <- names(.Object@values)

            if (length(parameter.names)==0) {
              .Object@lower.bound = numeric()
              .Object@upper.bound = numeric()
              return(.Object)
            }

            for (i in 1:length(parameter.names)) {

              if (parameter.names[i]=='lambda') {
                .Object@lower.bound[[parameter.names[i]]] <- 1e-10;
                .Object@upper.bound[[parameter.names[i]]] <- 1;
              } else if (parameter.names[i]=='lambda.burbidge') {
                .Object@lower.bound[[parameter.names[i]]] <- 1;
                .Object@upper.bound[[parameter.names[i]]] <- 10^3;
              } else if (parameter.names[i]=='mean.a0' || parameter.names[i]=='mean.summer.a0' || parameter.names[i]=='mean.autumn.a0' ||
                         parameter.names[i]=='mean.winter.a0' || parameter.names[i]=='mean.spring.a0') {
                 # .Object@lower.bound[[parameter.names[i]]] <- -1;
                 # .Object@upper.bound[[parameter.names[i]]] <- 1;
                .Object@lower.bound[[parameter.names[i]]] <- -0.5;
                .Object@upper.bound[[parameter.names[i]]] <- 0.5;
              } else if (parameter.names[i]=='mean.a0.amp') {
                .Object@lower.bound[[parameter.names[i]]] <- 0;
                .Object@upper.bound[[parameter.names[i]]] <- 1;
              } else if (parameter.names[i]=='mean.a0.disp') {
                .Object@lower.bound[[parameter.names[i]]] <- -0.5;
                .Object@upper.bound[[parameter.names[i]]] <- 0.5;
              } else if (parameter.names[i]=='mean.a0.phase') {
                .Object@lower.bound[[parameter.names[i]]] <- -1;
                .Object@upper.bound[[parameter.names[i]]] <- 1;
              } else if (parameter.names[i]=='mean.a1' || parameter.names[i]=='mean.summer.a1' || parameter.names[i]=='mean.autumn.a1' ||
                         parameter.names[i]=='mean.winter.a1' || parameter.names[i]=='mean.spring.a1') {
                .Object@lower.bound[[parameter.names[i]]] <- 1e-5;
                .Object@upper.bound[[parameter.names[i]]] <- 1;
              } else if (parameter.names[i]=='mean.trend') {
                .Object@lower.bound[[parameter.names[i]]] <- -0.1;
                .Object@upper.bound[[parameter.names[i]]] <- 0.1;
              } else if (parameter.names[i]=='mean.a1.amp') {
                .Object@lower.bound[[parameter.names[i]]] <-  1e-7;
                .Object@upper.bound[[parameter.names[i]]] <-  1;
              } else if (parameter.names[i]=='mean.a1.disp') {
                .Object@lower.bound[[parameter.names[i]]] <- 1e-5;
                .Object@upper.bound[[parameter.names[i]]] <-  1;
              } else if (parameter.names[i]=='mean.a1.phase') {
                .Object@lower.bound[[parameter.names[i]]] <- -1;
                .Object@upper.bound[[parameter.names[i]]] <-  1;
              } else if (parameter.names[i]=='mean.AR1') {
                .Object@lower.bound[[parameter.names[i]]] <- -1;
                .Object@upper.bound[[parameter.names[i]]] <- 1;
              } else if (parameter.names[i]=='mean.AR2') {
                .Object@lower.bound[[parameter.names[i]]] <- -1;
                .Object@upper.bound[[parameter.names[i]]] <- 1;
              } else if (parameter.names[i]=='mean.AR3') {
                .Object@lower.bound[[parameter.names[i]]] <- -1;
                .Object@upper.bound[[parameter.names[i]]] <- 1;
              } else if (parameter.names[i]=='mean.SAR1') {
                .Object@lower.bound[[parameter.names[i]]] <- -1;
                .Object@upper.bound[[parameter.names[i]]] <- 1;
              } else if (parameter.names[i]=='std.a0' || parameter.names[i]=='std.summer.a0' || parameter.names[i]=='std.autumn.a0' ||
                         parameter.names[i]=='std.winter.a0' || parameter.names[i]=='std.spring.a0') {
                .Object@lower.bound[[parameter.names[i]]] <- 0.05;
                .Object@upper.bound[[parameter.names[i]]] <- 0.75;
              } else if (parameter.names[i]=='std.a0.amp') {
                .Object@lower.bound[[parameter.names[i]]] <-  0;
                .Object@upper.bound[[parameter.names[i]]] <-  1;
              } else if (parameter.names[i]=='std.a0.disp') {
                .Object@lower.bound[[parameter.names[i]]] <- -1;
                .Object@upper.bound[[parameter.names[i]]] <-  1;
              } else if (parameter.names[i]=='std.a0.phase') {
                .Object@lower.bound[[parameter.names[i]]] <- -1;
                .Object@upper.bound[[parameter.names[i]]] <-  1;
              } else if (parameter.names[i]=='shape.a0') {
                .Object@lower.bound[[parameter.names[i]]] <- -10;
                .Object@upper.bound[[parameter.names[i]]] <- 10;
              } else if (parameter.names[i]=='transition.prob') {
                .Object@lower.bound[[parameter.names[i]]] <- 0;
                .Object@upper.bound[[parameter.names[i]]] <- 1;
              } else if (parameter.names[i]=='initial.state.prob') {
                .Object@lower.bound[[parameter.names[i]]] <- 0;
                .Object@upper.bound[[parameter.names[i]]] <- 1;
              } else {
                stop(paste('Default bounds are not defined for the parameter:',parameter.names[i]))
              }

            }

            return(.Object)
          }
)

setGeneric(name="setTransforms",def=function(.Object) {standardGeneric("setTransforms")})
setMethod(f="setTransforms",
          signature="parameters",
          definition=function(.Object)
          {

            parameter.names <- names(.Object@values)

            if (length(parameter.names)==0) {
              .Object@use.log.transform = logical()
              return(.Object)
            }

            for (i in 1:length(parameter.names)) {
              if (parameter.names[i]=='lambda') {
                .Object@use.log.transform[[parameter.names[i]]] <- T
              } else if (parameter.names[i]=='lambda.burbidge') {
                  .Object@use.log.transform[[parameter.names[i]]] <- T
              } else if (parameter.names[i]=='mean.a0' || parameter.names[i]=='mean.summer.a0' || parameter.names[i]=='mean.autumn.a0' ||
                         parameter.names[i]=='mean.winter.a0' || parameter.names[i]=='mean.spring.a0' ||
                         parameter.names[i]=='mean.a0.amp' || parameter.names[i]=='mean.a0.disp' || parameter.names[i]=='mean.a0.phase') {
                .Object@use.log.transform[[parameter.names[i]]] <- F
              } else if (parameter.names[i]=='mean.trend') {
                .Object@use.log.transform[[parameter.names[i]]] <- F
              } else if (parameter.names[i]=='mean.a1' || parameter.names[i]=='mean.summer.a1' || parameter.names[i]=='mean.autumn.a1' ||
                         parameter.names[i]=='mean.winter.a1' || parameter.names[i]=='mean.spring.a1' ||
                         parameter.names[i]=='mean.a1.amp' || parameter.names[i]=='mean.a1.disp') {
                .Object@use.log.transform[[parameter.names[i]]] <- T;
              } else if (parameter.names[i]=='mean.a1.phase') {
                .Object@use.log.transform[[parameter.names[i]]] <- F;
              } else if (parameter.names[i]=='mean.AR1') {
                .Object@use.log.transform[[parameter.names[i]]] <- F;
              } else if (parameter.names[i]=='mean.AR2') {
                .Object@use.log.transform[[parameter.names[i]]] <- F;
              } else if (parameter.names[i]=='mean.AR3') {
                .Object@use.log.transform[[parameter.names[i]]] <- F;
              } else if (parameter.names[i]=='mean.SAR1') {
                .Object@use.log.transform[[parameter.names[i]]] <- F;
              } else if (parameter.names[i]=='std.a0' || parameter.names[i]=='std.summer.a0' || parameter.names[i]=='std.autumn.a0' ||
                         parameter.names[i]=='std.winter.a0' || parameter.names[i]=='std.spring.a0' ||
                         parameter.names[i]=='std.a0.amp' || parameter.names[i]=='std.a0.disp' || parameter.names[i]=='std.a0.phase') {
                .Object@use.log.transform[[parameter.names[i]]] <- F;
              } else if (parameter.names[i]=='shape.a0') {
                .Object@use.log.transform[[parameter.names[i]]] <- F;
              } else if (parameter.names[i]=='transition.prob') {
                .Object@use.log.transform[[parameter.names[i]]] <- F;
              } else if (parameter.names[i]=='initial.state.prob') {
                .Object@use.log.transform[[parameter.names[i]]] <- F;
              } else {
                stop(paste('Default bounds are not defined for the parameter:',parameter.names[i]))
              }

            }

            return(.Object)
          }
)
# create a method to assign the parameter
#setGeneric(name="setParameters",def=function(.Object,parameters) {standardGeneric("setParameters")})
setMethod(f="setParameters",
          signature="parameters",
          definition=function(.Object,parameters)
          {

            parameter.name <- names(.Object@values)

            for (i in 1:length(parameter.name)) {

              if (length(parameters[[parameter.name[i]]]) != length(.Object@values[[parameter.name[i]]])) {
                stop(paste('The length of input parameter "',parameter.name[i],'" does not equal that within the model object.'))
              }
              .Object@values[[parameter.name[i]]] <- parameters[[parameter.name[i]]]

            }

            return(.Object)
          }
)

setMethod(f="setParameters.fromTransformed",
          signature="parameters",
          definition=function(.Object,parameters)
          {

            parameter.name <- names(.Object@values)

            for (i in 1:length(parameter.name)) {
              if (.Object@use.log.transform[[parameter.name[i]]]) {
                parameters[[parameter.name[i]]] <- 10^parameters[[parameter.name[i]]]
              }
            }

            .Object <- setParameters(.Object,parameters)
            return(.Object)
          }
)


# create a method to assign the parameter
#setGeneric(name="getParameters",def=function(.Object,position) {standardGeneric("getParameters")})
setMethod(f="getParameters",
          signature="parameters",
          definition=function(.Object)
          {
            return(.Object@values)
          }
)

# create a method to assign the parameter
#setGeneric(name="getParameters",def=function(.Object,position) {standardGeneric("getParameters")})
setMethod(f="getParameters.transformed",
          signature="parameters",
          definition=function(.Object)
          {

            parameters <- getParameters(.Object)
            parameter.name <- names(parameters)

            if (length(parameter.name)==0) {
              return(list())
            }

            for (i in 1:length(parameter.name)) {
              if (.Object@use.log.transform[[parameter.name[i]]]) {
                parameters[[parameter.name[i]]] <- log10(parameters[[parameter.name[i]]])
              }
            }

            return(parameters)

          }
)


setMethod(f="getBounds.transformed",
          signature="parameters",
          definition=function(.Object)
          {
            lowerBound = vector("list")
            upperBound = vector("list")
            parameter.name <- names(.Object@values)

            if (length(parameter.name)==0) {
              return(list(lower = numeric(), upper = numeric()))
            }

            for (i in 1:length(parameter.name)) {
              upperBound[[parameter.name[i]]] = rep(.Object@upper.bound[[parameter.name[i]]], length(.Object@values[[parameter.name[i]]]))
              lowerBound[[parameter.name[i]]] = rep(.Object@lower.bound[[parameter.name[i]]], length(.Object@values[[parameter.name[i]]]))

              if (.Object@use.log.transform[[parameter.name[i]]]) {
                lowerBound[[parameter.name[i]]] <- log10(lowerBound[[parameter.name[i]]])
                upperBound[[parameter.name[i]]] <- log10(upperBound[[parameter.name[i]]])
              }
              #stop("DBG")
            }
            return(list(lower = unlist(lowerBound), upper = unlist(upperBound)))
          }
)

setMethod(f="getBounds",
          signature="parameters",
          definition=function(.Object)
          {
            lowerBound = .Object@lower.bound
            upperBound = .Object@upper.bound
            parameter.name <- names(.Object@values)

            if (length(parameter.name)==0) {
              return(list(lower = numeric(), upper = numeric()))
            }

            return(list(lower = lowerBound, upper = upperBound))
          }
)
