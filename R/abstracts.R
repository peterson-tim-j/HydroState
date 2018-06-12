##' @export
##'
Qhat <- setClass("Qhat", package='hydroState')

#' @exportMethod setParameters
setGeneric(name="setParameters",def=function(.Object,parameters) {standardGeneric("setParameters")})
setMethod(f="setParameters",signature="Qhat",  definition=function(.Object,parameters){})

#' @exportMethod setParameters.fromTransformed
setGeneric(name="setParameters.fromTransformed",def=function(.Object,parameters) {standardGeneric("setParameters.fromTransformed")})
setMethod(f="setParameters.fromTransformed",signature="Qhat",  definition=function(.Object,parameters){})

#' @exportMethod getParameters
setGeneric(name="getParameters",def=function(.Object) {standardGeneric("getParameters")})
setMethod(f="getParameters",signature="Qhat",  definition=function(.Object){})

#' @exportMethod getTransformedParameterBounds
setGeneric(name="getTransformedParameterBounds",def=function(.Object) {standardGeneric("getTransformedParameterBounds")})
setMethod(f="getTransformedParameterBounds",signature="Qhat",  definition=function(.Object){})

#' @exportMethod getQhat
setGeneric(name="getQhat",def=function(.Object,data) {standardGeneric("getQhat")})
setMethod(f="getQhat",signature="Qhat",  definition=function(.Object, data){})

Qhatbar <- setClass("Qhatbar", package='hydroState')

#' @exportMethod setParameters
setGeneric(name="setParameters",def=function(.Object,parameters) {standardGeneric("setParameters")})
setMethod(f="setParameters",signature="Qhatbar",  definition=function(.Object,parameters){})

#' @exportMethod setParameters.fromTransformed
setGeneric(name="setParameters.fromTransformed",def=function(.Object,parameters) {standardGeneric("setParameters.fromTransformed")})
setMethod(f="setParameters.fromTransformed",signature="Qhatbar",  definition=function(.Object,parameters){})

#' @exportMethod getParameters
setGeneric(name="getParameters",def=function(.Object) {standardGeneric("getParameters")})
setMethod(f="getParameters",signature="Qhatbar",  definition=function(.Object){})

#' @exportMethod getTransformedParameterBounds
setGeneric(name="getTransformedParameterBounds",def=function(.Object) {standardGeneric("getTransformedParameterBounds")})
setMethod(f="getTransformedParameterBounds",signature="Qhatbar",  definition=function(.Object){})

#' @exportMethod getQhatbar
setGeneric(name="getQhatbar",def=function(.Object,data) {standardGeneric("getQhatbar")})
setMethod(f="getQhatbar",signature=c("Qhatbar","numeric"),  definition=function(.Object, data){})

markov <- setClass("markov", package='hydroState')

#' @exportMethod setParameters
setGeneric(name="setParameters",def=function(.Object,parameters) {standardGeneric("setParameters")})
setMethod(f="setParameters",signature="markov",  definition=function(.Object,parameters){})

#' @exportMethod setParameters.fromTransformed
setGeneric(name="setParameters.fromTransformed",def=function(.Object,parameters) {standardGeneric("setParameters.fromTransformed")})
setMethod(f="setParameters.fromTransformed",signature="markov",  definition=function(.Object,parameters){})

#' @exportMethod getParameters
setGeneric(name="getParameters",def=function(.Object) {standardGeneric("getParameters")})
setMethod(f="getParameters",signature="markov",  definition=function(.Object){})

#' @exportMethod getTransformedParameterBounds
setGeneric(name="getTransformedParameterBounds",def=function(.Object) {standardGeneric("getTransformedParameterBounds")})
setMethod(f="getTransformedParameterBounds",signature="markov",  definition=function(.Object){})

#' @exportMethod getNegLogLikelihood
setGeneric(name="getNegLogLikelihood",def=function(.Object, ...) {standardGeneric("getNegLogLikelihood")})

#' @exportMethod getLikelihood
setGeneric(name="getLikelihood",def=function(.Object, data) {standardGeneric("getLikelihood")})
