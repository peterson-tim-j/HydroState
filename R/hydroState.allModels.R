##' @include hydroState.R
hydroState.allModels.allModels <- setClass(
  # Set the name for the class
  "hydroState.allModels",

  package='hydroState',

  # Define the slots
  slots = c(
    one.state.linear.boxcox = 'hydroState',
    one.state.linear.AR1.boxcox = 'hydroState',
    one.state.linear.AR2.boxcox = 'hydroState',
    one.state.linear.AR3.boxcox = 'hydroState',
    one.state.linear.log = 'hydroState',
    one.state.linear.AR1.log = 'hydroState',
    one.state.linear.AR2.log = 'hydroState',
    one.state.linear.AR3.log = 'hydroState',
    two.state.linear.boxcox = 'hydroState',
    two.state.linear.AR1.boxcox = 'hydroState',
    two.state.linear.AR2.boxcox = 'hydroState',
    two.state.linear.AR3.boxcox = 'hydroState',
    two.state.linear.log = 'hydroState',
    two.state.linear.AR1.log = 'hydroState',
    two.state.linear.AR2.log = 'hydroState',
    two.state.linear.AR3.log = 'hydroState',
    three.state.linear.boxcox = 'hydroState',
    three.state.linear.AR1.boxcox = 'hydroState',
    three.state.linear.AR2.boxcox = 'hydroState',
    three.state.linear.AR3.boxcox = 'hydroState',
    three.state.linear.log = 'hydroState',
    three.state.linear.AR1.log = 'hydroState',
    three.state.linear.AR2.log = 'hydroState',
    three.state.linear.AR3.log = 'hydroState',
    three.state.unstructured.linear.boxcox = 'hydroState',
    three.state.unstructured.AR1.boxcox = 'hydroState',
    three.state.unstructured.AR2.boxcox = 'hydroState',
    three.state.unstructured.AR3.boxcox = 'hydroState',
    three.state.unstructured.log = 'hydroState',
    three.state.unstructured.AR1.log = 'hydroState',
    three.state.unstructured.AR2.log = 'hydroState',
    three.state.unstructured.AR3.log = 'hydroState'

  ),

  # Set the default values for the slots. (optional)
  # prototype=list(
  #   input.data = data.frame(year=c(0),month=c(0),precipitation=c(0),flow=c(0)),
  #   Qhat.object  = new('Qhat.boxcox',input.data = data.frame(year=c(0),month=c(0),precipitation=c(0),flow=c(0))),
  #   QhatModel.object = new('QhatModel.homo.normal.linear',input.data = data.frame(year=c(0),month=c(0),precipitation=c(0),flow=c(0))),
  #   markov.model.object = new('markov.annualHomogeneous',transition.graph = matrix(TRUE,2,2)),
  #   calibration.results = list(optim=0,member=0),
  #   state.labels = c('')
  # )
)

# Valid object?
validObject <- function(object) {
  TRUE
}
setValidity("hydroState.allModels", validObject)

# Initialise the object.
#setGeneric(name="initialize",def=function(.Object,input.data, Qhat.object, QhatModel.object, markov.model.object, ...){standardGeneric("initialize")})
setMethod(f="initialize",signature="hydroState.allModels",definition=function(.Object, input.data, Qhat.object, QhatModel.object, markov.model.object)
          {
          # .Object@input.data = input.data
          # .Object@Qhat.object = Qhat.object
          # .Object@QhatModel.object = QhatModel.object
          # .Object@markov.model.object = markov.model.object
          # .Object
          }
)

