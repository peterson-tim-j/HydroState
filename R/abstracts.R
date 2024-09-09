# @export
##'

parameters <- setClass("parameters", package='hydroState')

# @exportMethod setParameters
setGeneric(name="setParameters",def=function(.Object,parameters) {standardGeneric("setParameters")})
setMethod(f="setParameters",signature="parameters",  definition=function(.Object,parameters){})

# @exportMethod setParameters.fromTransformed
setGeneric(name="setParameters.fromTransformed",def=function(.Object,parameters) {standardGeneric("setParameters.fromTransformed")})
setMethod(f="setParameters.fromTransformed",signature="parameters",  definition=function(.Object,parameters){})

# @exportMethod getParameters
setGeneric(name="getParameters",def=function(.Object) {standardGeneric("getParameters")})
setMethod(f="getParameters",signature="parameters",  definition=function(.Object){})

# @exportMethod getParameters.transformed
setGeneric(name="getParameters.transformed",def=function(.Object) {standardGeneric("getParameters.transformed")})
setMethod(f="getParameters.transformed",signature="parameters",  definition=function(.Object){})

# @exportMethod getBounds.transformed
setGeneric(name="getBounds.transformed",def=function(.Object) {standardGeneric("getBounds.transformed")})
setMethod(f="getBounds.transformed",signature="parameters",  definition=function(.Object){})

# @exportMethod getBounds
setGeneric(name="getBounds",def=function(.Object) {standardGeneric("getBounds")})
setMethod(f="getBounds",signature="parameters",  definition=function(.Object){})

Qhat <- setClass("Qhat", package='hydroState')

# @exportMethod getQhat
setGeneric(name="getQhat",def=function(.Object,data) {standardGeneric("getQhat")})
setMethod(f="getQhat",signature="Qhat",  definition=function(.Object, data){})

# @exportMethod getQ.backTransformed
setGeneric(name="getQ.backTransformed",def=function(.Object,data) {standardGeneric("getQ.backTransformed")})
setMethod(f="getQ.backTransformed",signature="Qhat",  definition=function(.Object, data){})

QhatModel <- setClass("QhatModel", package='hydroState')


# @exportMethod getEmissionDensity
setGeneric(name="getEmissionDensity",def=function(.Object, data, cumProb.threshold.Qhat=NA) {standardGeneric("getEmissionDensity")})
setMethod(f="getEmissionDensity",signature="QhatModel",  definition=function(.Object, data, cumProb.threshold.Qhat=NA){})

# @exportMethod getDistributionPercentiles
setGeneric(name="getDistributionPercentiles",def=function(.Object, data, precentiles) {standardGeneric("getDistributionPercentiles")})
setMethod(f="getDistributionPercentiles",signature="QhatModel",  definition=function(.Object, data, precentiles){})

setGeneric(name="generate.sample.Qhat.fromViterbi",def=function(.Object, data, viterbi.states) {standardGeneric("generate.sample.Qhat.fromViterbi")})
setMethod(f="generate.sample.Qhat.fromViterbi",signature="QhatModel",  definition=function(.Object, data, viterbi.states){})

setGeneric(name="generate.sample.Qhat",def=function(.Object, data, nSamples) {standardGeneric("generate.sample.Qhat")})
setMethod(f="generate.sample.Qhat",signature="QhatModel",  definition=function(.Object, data, nSamples){})

setGeneric(name="getDensityIncrements",def=function(.Object, data, nIncrements) {standardGeneric("getDensityIncrements")})
setMethod(f="getDensityIncrements",signature="QhatModel",definition=function(.Object, data, nIncrements){})

# setGeneric(name="getZScore.atObs",def=function(.Object, data) {standardGeneric("getZScore.atObs")})
# setMethod(f="getZScore.atObs",signature="QhatModel",definition=function(.Object, data){})


markov <- setClass("markov", package='hydroState')

# @exportMethod getLogLikelihood
setGeneric(name="getLogLikelihood", def=function(.Object, data, emission.probs) {standardGeneric("getLogLikelihood")})
setMethod(f="getLogLikelihood",signature="markov",  definition=function(.Object){})

# @exportMethod getTransitionProbabilities
setGeneric(name="getTransitionProbabilities",def=function(.Object) {standardGeneric("getTransitionProbabilities")})
setMethod(f="getTransitionProbabilities",signature="markov",  definition=function(.Object){})

# @exportMethod getInitialStateProbabilities
setGeneric(name="getInitialStateProbabilities",def=function(.Object) {standardGeneric("getInitialStateProbabilities")})
setMethod(f="getInitialStateProbabilities",signature="markov",  definition=function(.Object){})

setGeneric(name="generate.sample.states",def=function(.Object, data) {standardGeneric("generate.sample.states")})
setMethod(f="generate.sample.states",signature="markov",  definition=function(.Object, data){})

setGeneric(name="getConditionalStateProbabilities",def=function(.Object, data, emission.probs) {standardGeneric("getConditionalStateProbabilities")})
setMethod(f="getConditionalStateProbabilities",signature="markov",  definition=function(.Object, data, emission.probs){})

#setGeneric(name="getConditionalProbabilities",def=function(.Object, data, emission.probs, zscore.increment.probs) {standardGeneric("getConditionalProbabilities")})
#setMethod(f="getConditionalProbabilities",signature="markov",  definition=function(.Object, data, emission.probs, zscore.increment.probs){})
setGeneric(name="getConditionalProbabilities",def=function(.Object, data, emission.probs, cumprob.atQhatIncrements) {standardGeneric("getConditionalProbabilities")})
setMethod(f="getConditionalProbabilities",signature="markov",  definition=function(.Object, data, emission.probs, cumprob.atQhatIncrements){})
