# @export
##'
##' @import methods

parameters <- setClass("parameters", package='hydroState')
#' @exportMethod parameters
#' @exportMethod parameters<-
setGeneric(name="parameters",   def=function(.Object, transformed, as_df) {standardGeneric("parameters")})
setGeneric(name="parameters<-", def=function(.Object, transformed, value){standardGeneric("parameters<-")})

#' @exportMethod bounds
#' @exportMethod bounds<-
setGeneric(name="bounds",   def=function(.Object, transformed) {standardGeneric("bounds")})
setGeneric(name="bounds<-", def=function(.Object, value){standardGeneric("bounds<-")})

#' @exportMethod transforms
#' @exportMethod transforms<-
setGeneric(name="transforms",  def=function(.Object) {standardGeneric("transforms")})
setGeneric(name="transforms<-",def=function(.Object, value) {standardGeneric("transforms<-")})

Qhat <- setClass("Qhat", package='hydroState')
setGeneric(name="getQhat",def=function(.Object,data) {standardGeneric("getQhat")})
setMethod(f="getQhat",signature="Qhat",  definition=function(.Object, data){})
setGeneric(name="getQ.backTransformed",def=function(.Object,data) {standardGeneric("getQ.backTransformed")})
setMethod(f="getQ.backTransformed",signature="Qhat",  definition=function(.Object, data){})
setGeneric(name="get.zeroFlow",def=function(.Object) {standardGeneric("get.zeroFlow")})
setMethod(f="get.zeroFlow",signature="Qhat",  definition=function(.Object){})

QhatModel <- setClass("QhatModel", package='hydroState')
setGeneric(name="getEmissionDensity",def=function(.Object, data, zero.Flow, cumProb.threshold.Qhat=NA) {standardGeneric("getEmissionDensity")})
setMethod(f="getEmissionDensity",signature="QhatModel",  definition=function(.Object, data, zero.Flow, cumProb.threshold.Qhat=NA){})
setGeneric(name="getDistributionPercentiles",def=function(.Object, data, precentiles) {standardGeneric("getDistributionPercentiles")})
setMethod(f="getDistributionPercentiles",signature="QhatModel",  definition=function(.Object, data, precentiles){})
setGeneric(name="generate.sample.Qhat.fromViterbi",def=function(.Object, data, viterbi.states) {standardGeneric("generate.sample.Qhat.fromViterbi")})
setMethod(f="generate.sample.Qhat.fromViterbi",signature="QhatModel",  definition=function(.Object, data, viterbi.states){})
setGeneric(name="generate.sample.Qhat",def=function(.Object, data, nSamples) {standardGeneric("generate.sample.Qhat")})
setMethod(f="generate.sample.Qhat",signature="QhatModel",  definition=function(.Object, data, nSamples){})
setGeneric(name="getDensityIncrements",def=function(.Object, data, nIncrements) {standardGeneric("getDensityIncrements")})
setMethod(f="getDensityIncrements",signature="QhatModel",definition=function(.Object, data, nIncrements){})

markov <- setClass("markov", package='hydroState')
setGeneric(name="getLogLikelihood", def=function(.Object, data, emission.probs) {standardGeneric("getLogLikelihood")})
setMethod(f="getLogLikelihood",signature="markov",  definition=function(.Object){})
setGeneric(name="getTransitionProbabilities",def=function(.Object) {standardGeneric("getTransitionProbabilities")})
setMethod(f="getTransitionProbabilities",signature="markov",  definition=function(.Object){})
setGeneric(name="getInitialStateProbabilities",def=function(.Object) {standardGeneric("getInitialStateProbabilities")})
setMethod(f="getInitialStateProbabilities",signature="markov",  definition=function(.Object){})
setGeneric(name="generate.sample.states",def=function(.Object, data) {standardGeneric("generate.sample.states")})
setMethod(f="generate.sample.states",signature="markov",  definition=function(.Object, data){})
setGeneric(name="getConditionalStateProbabilities",def=function(.Object, data, emission.probs) {standardGeneric("getConditionalStateProbabilities")})
setMethod(f="getConditionalStateProbabilities",signature="markov",  definition=function(.Object, data, emission.probs){})
setGeneric(name="getConditionalProbabilities",def=function(.Object, data, emission.probs, cumprob.atQhatIncrements) {standardGeneric("getConditionalProbabilities")})
setMethod(f="getConditionalProbabilities",signature="markov",  definition=function(.Object, data, emission.probs, cumprob.atQhatIncrements){})
