##' @include hydroState.allModels.R
##' @export
hydroState.subAnnual.allModels <- setClass(
  # Set the name for the class
  "hydroState.subAnnual.allModels",

  package='hydroState',

  contains=c('hydroState.allModels')

)

# Valid object?
validObject <- function(object) {
  TRUE
}
setValidity("hydroState.allModels", validObject)

# Initialise the object.
#setGeneric(name="initialize",def=function(.Object,input.data, Qhat.object, QhatModel.object, markov.model.object, ...){standardGeneric("initialize")})
setMethod(f="initialize",signature="hydroState.subAnnual.allModels",definition=function(.Object, siteID, input.data, allow.flickering=T, build.3state.models=F)
{

  .Object@siteID <- siteID

  # Define transition graphs
  transition.graph.1State <- matrix(TRUE,1,1)
  transition.graph.2State <- matrix(TRUE,2,2)
  if (build.3state.models) {
    transition.graph.3State.unstructured <- matrix(TRUE,3,3)
    transition.graph.3State <- matrix(TRUE,3,3)
    transition.graph.3State[1,3] <- FALSE
    transition.graph.3State[2,1] <- FALSE
    transition.graph.3State[3,2] <- FALSE
  }

  # Build Qhat flow transform objects.
  Qhat.boxcox = new('Qhat.boxcox', input.data=input.data)
  Qhat.log = new('Qhat.log', input.data=input.data)

  # Check if the input has monthly data. If so, then use the seasonal AR version of the Qhat models.
  if (!any(names(input.data)=='month')) {
    stop('Input data must be sub-annual and contain a column named "month".')
  }
  if (allow.flickering) {
    model.extension.markov = '.flickering'
  } else {
    model.extension.markov = ''
  }

  # Build 1 state Qhat models
  QhatModel.1State.gamma.subAnnual_a0 = new('QhatModel.subAnnual.homo.gamma.linear', input.data=input.data, transition.graph=transition.graph.1State,
                                         state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                         subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=F)
  QhatModel.1State.gamma.subAnnual_a0_a1 = new('QhatModel.subAnnual.homo.gamma.linear', input.data=input.data, transition.graph=transition.graph.1State,
                                         state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                         subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=F)
  QhatModel.1State.gamma.subAnnual_a0_std = new('QhatModel.subAnnual.homo.gamma.linear', input.data=input.data, transition.graph=transition.graph.1State,
                                         state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                         subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=T)
  QhatModel.1State.gamma.subAnnual_a0_a1_std = new('QhatModel.subAnnual.homo.gamma.linear', input.data=input.data, transition.graph=transition.graph.1State,
                                         state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                         subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=T)

  QhatModel.1State.gamma.AR1.subAnnual_a0 = new('QhatModel.subAnnual.homo.gamma.linear.AR1', input.data=input.data, transition.graph=transition.graph.1State,
                                           state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                           subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=F)
  QhatModel.1State.gamma.AR1.subAnnual_a0_a1 = new('QhatModel.subAnnual.homo.gamma.linear.AR1', input.data=input.data, transition.graph=transition.graph.1State,
                                           state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                           subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=F)
  QhatModel.1State.gamma.AR1.subAnnual_a0_std = new('QhatModel.subAnnual.homo.gamma.linear.AR1', input.data=input.data, transition.graph=transition.graph.1State,
                                           state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                           subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=T)
  QhatModel.1State.gamma.AR1.subAnnual_a0_a1_std = new('QhatModel.subAnnual.homo.gamma.linear.AR1', input.data=input.data, transition.graph=transition.graph.1State,
                                            state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                            subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=T)


  QhatModel.1State.gamma.AR2.subAnnual_a0 = new('QhatModel.subAnnual.homo.gamma.linear.AR2', input.data=input.data, transition.graph=transition.graph.1State,
                                           state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                           subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=F)
  QhatModel.1State.gamma.AR2.subAnnual_a0_a1 = new('QhatModel.subAnnual.homo.gamma.linear.AR2', input.data=input.data, transition.graph=transition.graph.1State,
                                           state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                           subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=F)
  QhatModel.1State.gamma.AR2.subAnnual_a0_std = new('QhatModel.subAnnual.homo.gamma.linear.AR2', input.data=input.data, transition.graph=transition.graph.1State,
                                           state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                           subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=T)
  QhatModel.1State.gamma.AR2.subAnnual_a0_a1_std = new('QhatModel.subAnnual.homo.gamma.linear.AR2', input.data=input.data, transition.graph=transition.graph.1State,
                                           state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                           subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=T)

  QhatModel.1State.gamma.AR3.subAnnual_a0 = new('QhatModel.subAnnual.homo.gamma.linear.AR3', input.data=input.data, transition.graph=transition.graph.1State,
                                               state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                               subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=F)
  QhatModel.1State.gamma.AR3.subAnnual_a0_a1 = new('QhatModel.subAnnual.homo.gamma.linear.AR3', input.data=input.data, transition.graph=transition.graph.1State,
                                                  state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                  subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=F)
  QhatModel.1State.gamma.AR3.subAnnual_a0_std = new('QhatModel.subAnnual.homo.gamma.linear.AR3', input.data=input.data, transition.graph=transition.graph.1State,
                                                   state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                   subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=T)
  QhatModel.1State.gamma.AR3.subAnnual_a0_a1_std = new('QhatModel.subAnnual.homo.gamma.linear.AR3', input.data=input.data, transition.graph=transition.graph.1State,
                                                      state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                      subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=T)


  # Build 2 state Qhat models
  QhatModel.2State.gamma.subAnnual_a0 = new('QhatModel.subAnnual.homo.gamma.linear', input.data=input.data, transition.graph=transition.graph.2State,
                                           state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                           subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=F)
  QhatModel.2State.gamma.subAnnual_a0_a1 = new('QhatModel.subAnnual.homo.gamma.linear', input.data=input.data, transition.graph=transition.graph.2State,
                                              state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                              subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=F)
  QhatModel.2State.gamma.subAnnual_a0_std = new('QhatModel.subAnnual.homo.gamma.linear', input.data=input.data, transition.graph=transition.graph.2State,
                                               state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                               subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=T)
  QhatModel.2State.gamma.subAnnual_a0_a1_std = new('QhatModel.subAnnual.homo.gamma.linear', input.data=input.data, transition.graph=transition.graph.2State,
                                                  state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                  subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=T)

  QhatModel.2State.gamma.AR1.subAnnual_a0 = new('QhatModel.subAnnual.homo.gamma.linear.AR1', input.data=input.data, transition.graph=transition.graph.2State,
                                               state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                               subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=F)
  QhatModel.2State.gamma.AR1.subAnnual_a0_a1 = new('QhatModel.subAnnual.homo.gamma.linear.AR1', input.data=input.data, transition.graph=transition.graph.2State,
                                                  state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                  subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=F)
  QhatModel.2State.gamma.AR1.subAnnual_a0_std = new('QhatModel.subAnnual.homo.gamma.linear.AR1', input.data=input.data, transition.graph=transition.graph.2State,
                                                   state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                   subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=T)
  QhatModel.2State.gamma.AR1.subAnnual_a0_a1_std = new('QhatModel.subAnnual.homo.gamma.linear.AR1', input.data=input.data, transition.graph=transition.graph.2State,
                                                      state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                      subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=T)


  QhatModel.2State.gamma.AR2.subAnnual_a0 = new('QhatModel.subAnnual.homo.gamma.linear.AR2', input.data=input.data, transition.graph=transition.graph.2State,
                                               state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                               subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=F)
  QhatModel.2State.gamma.AR2.subAnnual_a0_a1 = new('QhatModel.subAnnual.homo.gamma.linear.AR2', input.data=input.data, transition.graph=transition.graph.2State,
                                                  state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                  subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=F)
  QhatModel.2State.gamma.AR2.subAnnual_a0_std = new('QhatModel.subAnnual.homo.gamma.linear.AR2', input.data=input.data, transition.graph=transition.graph.2State,
                                                   state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                   subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=T)
  QhatModel.2State.gamma.AR2.subAnnual_a0_a1_std = new('QhatModel.subAnnual.homo.gamma.linear.AR2', input.data=input.data, transition.graph=transition.graph.2State,
                                                      state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                      subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=T)

  QhatModel.2State.gamma.AR3.subAnnual_a0 = new('QhatModel.subAnnual.homo.gamma.linear.AR3', input.data=input.data, transition.graph=transition.graph.2State,
                                               state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                               subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=F)
  QhatModel.2State.gamma.AR3.subAnnual_a0_a1 = new('QhatModel.subAnnual.homo.gamma.linear.AR3', input.data=input.data, transition.graph=transition.graph.2State,
                                                  state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                  subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=F)
  QhatModel.2State.gamma.AR3.subAnnual_a0_std = new('QhatModel.subAnnual.homo.gamma.linear.AR3', input.data=input.data, transition.graph=transition.graph.2State,
                                                   state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                   subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=T)
  QhatModel.2State.gamma.AR3.subAnnual_a0_a1_std = new('QhatModel.subAnnual.homo.gamma.linear.AR3', input.data=input.data, transition.graph=transition.graph.2State,
                                                      state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                      subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=T)

  # Build 3 state Qhat models
  if (build.3state.models) {
    QhatModel.3State.gamma.subAnnual_a0 = new('QhatModel.subAnnual.homo.gamma.linear', input.data=input.data, transition.graph=transition.graph.3State,
                                             state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                             subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=F)
    QhatModel.3State.gamma.subAnnual_a0_a1 = new('QhatModel.subAnnual.homo.gamma.linear', input.data=input.data, transition.graph=transition.graph.3State,
                                                state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=F)
    QhatModel.3State.gamma.subAnnual_a0_std = new('QhatModel.subAnnual.homo.gamma.linear', input.data=input.data, transition.graph=transition.graph.3State,
                                                 state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                 subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=T)
    QhatModel.3State.gamma.subAnnual_a0_a1_std = new('QhatModel.subAnnual.homo.gamma.linear', input.data=input.data, transition.graph=transition.graph.3State,
                                                    state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                    subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=T)

    QhatModel.3State.gamma.AR1.subAnnual_a0 = new('QhatModel.subAnnual.homo.gamma.linear.AR1', input.data=input.data, transition.graph=transition.graph.3State,
                                                 state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                 subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=F)
    QhatModel.3State.gamma.AR1.subAnnual_a0_a1 = new('QhatModel.subAnnual.homo.gamma.linear.AR1', input.data=input.data, transition.graph=transition.graph.3State,
                                                    state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                    subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=F)
    QhatModel.3State.gamma.AR1.subAnnual_a0_std = new('QhatModel.subAnnual.homo.gamma.linear.AR1', input.data=input.data, transition.graph=transition.graph.3State,
                                                     state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                     subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=T)
    QhatModel.3State.gamma.AR1.subAnnual_a0_a1_std = new('QhatModel.subAnnual.homo.gamma.linear.AR1', input.data=input.data, transition.graph=transition.graph.3State,
                                                        state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                        subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=T)


    QhatModel.3State.gamma.AR2.subAnnual_a0 = new('QhatModel.subAnnual.homo.gamma.linear.AR2', input.data=input.data, transition.graph=transition.graph.3State,
                                                 state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                 subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=F)
    QhatModel.3State.gamma.AR2.subAnnual_a0_a1 = new('QhatModel.subAnnual.homo.gamma.linear.AR2', input.data=input.data, transition.graph=transition.graph.3State,
                                                    state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                    subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=F)
    QhatModel.3State.gamma.AR2.subAnnual_a0_std = new('QhatModel.subAnnual.homo.gamma.linear.AR2', input.data=input.data, transition.graph=transition.graph.3State,
                                                     state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                     subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=T)
    QhatModel.3State.gamma.AR2.subAnnual_a0_a1_std = new('QhatModel.subAnnual.homo.gamma.linear.AR2', input.data=input.data, transition.graph=transition.graph.3State,
                                                        state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                        subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=T)

    QhatModel.3State.gamma.AR3.subAnnual_a0 = new('QhatModel.subAnnual.homo.gamma.linear.AR3', input.data=input.data, transition.graph=transition.graph.3State,
                                                 state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                 subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=F)
    QhatModel.3State.gamma.AR3.subAnnual_a0_a1 = new('QhatModel.subAnnual.homo.gamma.linear.AR3', input.data=input.data, transition.graph=transition.graph.3State,
                                                    state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                    subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=F)
    QhatModel.3State.gamma.AR3.subAnnual_a0_std = new('QhatModel.subAnnual.homo.gamma.linear.AR3', input.data=input.data, transition.graph=transition.graph.3State,
                                                     state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                     subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=T)
    QhatModel.3State.gamma.AR3.subAnnual_a0_a1_std = new('QhatModel.subAnnual.homo.gamma.linear.AR3', input.data=input.data, transition.graph=transition.graph.3State,
                                                        state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                        subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=T)


    # Build 3 state unstructured Qhat models
    QhatModel.3StateUS.gamma.subAnnual_a0 = new('QhatModel.subAnnual.homo.gamma.linear', input.data=input.data, transition.graph=transition.graph.3State.unstructured,
                                               state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                               subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=F)
    QhatModel.3StateUS.gamma.subAnnual_a0_a1 = new('QhatModel.subAnnual.homo.gamma.linear', input.data=input.data, transition.graph=transition.graph.3State.unstructured,
                                                  state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                  subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=F)
    QhatModel.3StateUS.gamma.subAnnual_a0_std = new('QhatModel.subAnnual.homo.gamma.linear', input.data=input.data, transition.graph=transition.graph.3State.unstructured,
                                                   state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                   subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=T)
    QhatModel.3StateUS.gamma.subAnnual_a0_a1_std = new('QhatModel.subAnnual.homo.gamma.linear', input.data=input.data, transition.graph=transition.graph.3State.unstructured,
                                                      state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                      subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=T)

    QhatModel.3StateUS.gamma.AR1.subAnnual_a0 = new('QhatModel.subAnnual.homo.gamma.linear.AR1', input.data=input.data, transition.graph=transition.graph.3State.unstructured,
                                                   state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                   subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=F)
    QhatModel.3StateUS.gamma.AR1.subAnnual_a0_a1 = new('QhatModel.subAnnual.homo.gamma.linear.AR1', input.data=input.data, transition.graph=transition.graph.3State.unstructured,
                                                      state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                      subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=F)
    QhatModel.3StateUS.gamma.AR1.subAnnual_a0_std = new('QhatModel.subAnnual.homo.gamma.linear.AR1', input.data=input.data, transition.graph=transition.graph.3State.unstructured,
                                                       state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                       subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=T)
    QhatModel.3StateUS.gamma.AR1.subAnnual_a0_a1_std = new('QhatModel.subAnnual.homo.gamma.linear.AR1', input.data=input.data, transition.graph=transition.graph.3State.unstructured,
                                                          state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                          subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=T)


    QhatModel.3StateUS.gamma.AR2.subAnnual_a0 = new('QhatModel.subAnnual.homo.gamma.linear.AR2', input.data=input.data, transition.graph=transition.graph.3State.unstructured,
                                                   state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                   subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=F)
    QhatModel.3StateUS.gamma.AR2.subAnnual_a0_a1 = new('QhatModel.subAnnual.homo.gamma.linear.AR2', input.data=input.data, transition.graph=transition.graph.3State.unstructured,
                                                      state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                      subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=F)
    QhatModel.3StateUS.gamma.AR2.subAnnual_a0_std = new('QhatModel.subAnnual.homo.gamma.linear.AR2', input.data=input.data, transition.graph=transition.graph.3State.unstructured,
                                                       state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                       subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=T)
    QhatModel.3StateUS.gamma.AR2.subAnnual_a0_a1_std = new('QhatModel.subAnnual.homo.gamma.linear.AR2', input.data=input.data, transition.graph=transition.graph.3State.unstructured,
                                                          state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                          subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=T)

    QhatModel.3StateUS.gamma.AR3.subAnnual_a0 = new('QhatModel.subAnnual.homo.gamma.linear.AR3', input.data=input.data, transition.graph=transition.graph.3State.unstructured,
                                                   state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                   subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=F)
    QhatModel.3StateUS.gamma.AR3.subAnnual_a0_a1 = new('QhatModel.subAnnual.homo.gamma.linear.AR3', input.data=input.data, transition.graph=transition.graph.3State.unstructured,
                                                      state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                      subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=F)
    QhatModel.3StateUS.gamma.AR3.subAnnual_a0_std = new('QhatModel.subAnnual.homo.gamma.linear.AR3', input.data=input.data, transition.graph=transition.graph.3State.unstructured,
                                                       state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                       subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=F, subAnnual.dependent.std.a0=T)
    QhatModel.3StateUS.gamma.AR3.subAnnual_a0_a1_std = new('QhatModel.subAnnual.homo.gamma.linear.AR3', input.data=input.data, transition.graph=transition.graph.3State.unstructured,
                                                          state.dependent.mean.a0=T, state.dependent.mean.a1=F, state.dependent.std.a0=T,
                                                          subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T, subAnnual.dependent.std.a0=T)

  }


  # Build Markov model object
  markov.1State = new(paste('markov.annualHomogeneous',model.extension.markov,sep=''), transition.graph=transition.graph.1State)
  markov.2State = new(paste('markov.annualHomogeneous',model.extension.markov,sep=''), transition.graph=transition.graph.2State)
  if (build.3state.models) {
    markov.3State = new(paste('markov.annualHomogeneous',model.extension.markov,sep=''), transition.graph=transition.graph.3State)
    markov.3State.US = new(paste('markov.annualHomogeneous',model.extension.markov,sep=''), transition.graph=transition.graph.3State.unstructured)
  }

  # Build hydrostate models
  .Object@models <- list(
    model.1State.gamma.log.subAnnual_a0 = new('hydroState',       input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.1State.gamma.subAnnual_a0      , markov.model.object=markov.1State),
    model.1State.gamma.log.subAnnual_a0_a1 = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.1State.gamma.subAnnual_a0_a1   , markov.model.object=markov.1State),
    model.1State.gamma.log.subAnnual_a0_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.1State.gamma.subAnnual_a0_std   , markov.model.object=markov.1State),
    model.1State.gamma.log.subAnnual_a0_a1_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.1State.gamma.subAnnual_a0_a1_std   , markov.model.object=markov.1State),

    model.1State.gamma.log.AR1.subAnnual_a0 = new('hydroState',       input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.1State.gamma.AR1.subAnnual_a0      , markov.model.object=markov.1State),
    model.1State.gamma.log.AR1.subAnnual_a0_a1 = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.1State.gamma.AR1.subAnnual_a0_a1   , markov.model.object=markov.1State),
    model.1State.gamma.log.AR1.subAnnual_a0_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.1State.gamma.AR1.subAnnual_a0_std   , markov.model.object=markov.1State),
    model.1State.gamma.log.AR1.subAnnual_a0_a1_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.1State.gamma.AR1.subAnnual_a0_a1_std   , markov.model.object=markov.1State),

    model.1State.gamma.log.AR2.subAnnual_a0 = new('hydroState',       input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.1State.gamma.AR2.subAnnual_a0      , markov.model.object=markov.1State),
    model.1State.gamma.log.AR2.subAnnual_a0_a1 = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.1State.gamma.AR2.subAnnual_a0_a1   , markov.model.object=markov.1State),
    model.1State.gamma.log.AR2.subAnnual_a0_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.1State.gamma.AR2.subAnnual_a0_std   , markov.model.object=markov.1State),
    model.1State.gamma.log.AR2.subAnnual_a0_a1_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.1State.gamma.AR2.subAnnual_a0_a1_std   , markov.model.object=markov.1State),

    model.1State.gamma.log.AR3.subAnnual_a0 = new('hydroState',       input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.1State.gamma.AR3.subAnnual_a0      , markov.model.object=markov.1State),
    model.1State.gamma.log.AR3.subAnnual_a0_a1 = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.1State.gamma.AR3.subAnnual_a0_a1   , markov.model.object=markov.1State),
    model.1State.gamma.log.AR3.subAnnual_a0_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.1State.gamma.AR3.subAnnual_a0_std   , markov.model.object=markov.1State),
    model.1State.gamma.log.AR3.subAnnual_a0_a1_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.1State.gamma.AR3.subAnnual_a0_a1_std   , markov.model.object=markov.1State),

    # 2 state models
    model.2State.gamma.log.subAnnual_a0 = new('hydroState',       input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.2State.gamma.subAnnual_a0      , markov.model.object=markov.2State),
    model.2State.gamma.log.subAnnual_a0_a1 = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.2State.gamma.subAnnual_a0_a1   , markov.model.object=markov.2State),
    model.2State.gamma.log.subAnnual_a0_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.2State.gamma.subAnnual_a0_std   , markov.model.object=markov.2State),
    model.2State.gamma.log.subAnnual_a0_a1_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.2State.gamma.subAnnual_a0_a1_std   , markov.model.object=markov.2State),

    model.2State.gamma.log.AR1.subAnnual_a0 = new('hydroState',       input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.2State.gamma.AR1.subAnnual_a0      , markov.model.object=markov.2State),
    model.2State.gamma.log.AR1.subAnnual_a0_a1 = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.2State.gamma.AR1.subAnnual_a0_a1   , markov.model.object=markov.2State),
    model.2State.gamma.log.AR1.subAnnual_a0_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.2State.gamma.AR1.subAnnual_a0_std   , markov.model.object=markov.2State),
    model.2State.gamma.log.AR1.subAnnual_a0_a1_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.2State.gamma.AR1.subAnnual_a0_a1_std   , markov.model.object=markov.2State),

    model.2State.gamma.log.AR2.subAnnual_a0 = new('hydroState',       input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.2State.gamma.AR2.subAnnual_a0      , markov.model.object=markov.2State),
    model.2State.gamma.log.AR2.subAnnual_a0_a1 = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.2State.gamma.AR2.subAnnual_a0_a1   , markov.model.object=markov.2State),
    model.2State.gamma.log.AR2.subAnnual_a0_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.2State.gamma.AR2.subAnnual_a0_std   , markov.model.object=markov.2State),
    model.2State.gamma.log.AR2.subAnnual_a0_a1_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.2State.gamma.AR2.subAnnual_a0_a1_std   , markov.model.object=markov.2State),

    model.2State.gamma.log.AR3.subAnnual_a0 = new('hydroState',       input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.2State.gamma.AR3.subAnnual_a0      , markov.model.object=markov.2State),
    model.2State.gamma.log.AR3.subAnnual_a0_a1 = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.2State.gamma.AR3.subAnnual_a0_a1   , markov.model.object=markov.2State),
    model.2State.gamma.log.AR3.subAnnual_a0_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.2State.gamma.AR3.subAnnual_a0_std   , markov.model.object=markov.2State),
    model.2State.gamma.log.AR3.subAnnual_a0_a1_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.2State.gamma.AR3.subAnnual_a0_a1_std   , markov.model.object=markov.2State)
  )

  # 3 state models
  if (build.3state.models) {
    .Object@models <- c(.Object@models, list(
      model.3State.gamma.log.subAnnual_a0 = new('hydroState',       input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3State.gamma.subAnnual_a0      , markov.model.object=markov.3State),
      model.3State.gamma.log.subAnnual_a0_a1 = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3State.gamma.subAnnual_a0_a1   , markov.model.object=markov.3State),
      model.3State.gamma.log.subAnnual_a0_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3State.gamma.subAnnual_a0_std   , markov.model.object=markov.3State),
      model.3State.gamma.log.subAnnual_a0_a1_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3State.gamma.subAnnual_a0_a1_std   , markov.model.object=markov.3State),

      model.3State.gamma.log.AR1.subAnnual_a0 = new('hydroState',       input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3State.gamma.AR1.subAnnual_a0      , markov.model.object=markov.3State),
      model.3State.gamma.log.AR1.subAnnual_a0_a1 = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3State.gamma.AR1.subAnnual_a0_a1   , markov.model.object=markov.3State),
      model.3State.gamma.log.AR1.subAnnual_a0_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3State.gamma.AR1.subAnnual_a0_std   , markov.model.object=markov.3State),
      model.3State.gamma.log.AR1.subAnnual_a0_a1_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3State.gamma.AR1.subAnnual_a0_a1_std   , markov.model.object=markov.3State),

      model.3State.gamma.log.AR2.subAnnual_a0 = new('hydroState',       input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3State.gamma.AR2.subAnnual_a0      , markov.model.object=markov.3State),
      model.3State.gamma.log.AR2.subAnnual_a0_a1 = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3State.gamma.AR2.subAnnual_a0_a1   , markov.model.object=markov.3State),
      model.3State.gamma.log.AR2.subAnnual_a0_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3State.gamma.AR2.subAnnual_a0_std   , markov.model.object=markov.3State),
      model.3State.gamma.log.AR2.subAnnual_a0_a1_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3State.gamma.AR2.subAnnual_a0_a1_std   , markov.model.object=markov.3State),

      model.3State.gamma.log.AR3.subAnnual_a0 = new('hydroState',       input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3State.gamma.AR3.subAnnual_a0      , markov.model.object=markov.3State),
      model.3State.gamma.log.AR3.subAnnual_a0_a1 = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3State.gamma.AR3.subAnnual_a0_a1   , markov.model.object=markov.3State),
      model.3State.gamma.log.AR3.subAnnual_a0_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3State.gamma.AR3.subAnnual_a0_std   , markov.model.object=markov.3State),
      model.3State.gamma.log.AR3.subAnnual_a0_a1_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3State.gamma.AR3.subAnnual_a0_a1_std   , markov.model.object=markov.3State),

      # 3 state unstructured models
      model.3StateUS.gamma.log.subAnnual_a0 = new('hydroState',       input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3StateUS.gamma.subAnnual_a0      , markov.model.object=markov.3State.US),
      model.3StateUS.gamma.log.subAnnual_a0_a1 = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3StateUS.gamma.subAnnual_a0_a1   , markov.model.object=markov.3State.US),
      model.3StateUS.gamma.log.subAnnual_a0_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3StateUS.gamma.subAnnual_a0_std   , markov.model.object=markov.3State.US),
      model.3StateUS.gamma.log.subAnnual_a0_a1_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3StateUS.gamma.subAnnual_a0_a1_std   , markov.model.object=markov.3State.US),

      model.3StateUS.gamma.log.AR1.subAnnual_a0 = new('hydroState',       input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3StateUS.gamma.AR1.subAnnual_a0      , markov.model.object=markov.3State.US),
      model.3StateUS.gamma.log.AR1.subAnnual_a0_a1 = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3StateUS.gamma.AR1.subAnnual_a0_a1   , markov.model.object=markov.3State.US),
      model.3StateUS.gamma.log.AR1.subAnnual_a0_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3StateUS.gamma.AR1.subAnnual_a0_std   , markov.model.object=markov.3State.US),
      model.3StateUS.gamma.log.AR1.subAnnual_a0_a1_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3StateUS.gamma.AR1.subAnnual_a0_a1_std   , markov.model.object=markov.3State.US),

      model.3StateUS.gamma.log.AR2.subAnnual_a0 = new('hydroState',       input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3StateUS.gamma.AR2.subAnnual_a0      , markov.model.object=markov.3State.US),
      model.3StateUS.gamma.log.AR2.subAnnual_a0_a1 = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3StateUS.gamma.AR2.subAnnual_a0_a1   , markov.model.object=markov.3State.US),
      model.3StateUS.gamma.log.AR2.subAnnual_a0_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3StateUS.gamma.AR2.subAnnual_a0_std   , markov.model.object=markov.3State.US),
      model.3StateUS.gamma.log.AR2.subAnnual_a0_a1_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3StateUS.gamma.AR2.subAnnual_a0_a1_std   , markov.model.object=markov.3State.US),

      model.3StateUS.gamma.log.AR3.subAnnual_a0 = new('hydroState',       input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3StateUS.gamma.AR3.subAnnual_a0      , markov.model.object=markov.3State.US),
      model.3StateUS.gamma.log.AR3.subAnnual_a0_a1 = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3StateUS.gamma.AR3.subAnnual_a0_a1   , markov.model.object=markov.3State.US),
      model.3StateUS.gamma.log.AR3.subAnnual_a0_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3StateUS.gamma.AR3.subAnnual_a0_std   , markov.model.object=markov.3State.US),
      model.3StateUS.gamma.log.AR3.subAnnual_a0_a1_std = new('hydroState',    input.data=input.data, Qhat.object = Qhat.log, QhatModel.object = QhatModel.3StateUS.gamma.AR3.subAnnual_a0_a1_std   , markov.model.object=markov.3State.US)
    ))
  }

  # Define Refernce models. That is, the calibration obecjive function that this models needs to meet or exceed.
  .Object@calib.reference.model.name <- list(
    model.1State.gamma.log.subAnnual_a0  = '',
    model.1State.gamma.log.subAnnual_a0_a1  = 'model.1State.gamma.log.subAnnual_a0',
    model.1State.gamma.log.subAnnual_a0_std  = 'model.1State.gamma.log.subAnnual_a0',
    model.1State.gamma.log.subAnnual_a0_a1_std  = 'model.1State.gamma.log.subAnnual_a0',

    model.1State.gamma.log.AR1.subAnnual_a0  = 'model.1State.gamma.log.subAnnual_a0',
    model.1State.gamma.log.AR1.subAnnual_a0_a1  = 'model.1State.gamma.log.AR1.subAnnual_a0',
    model.1State.gamma.log.AR1.subAnnual_a0_std  = 'model.1State.gamma.log.AR1.subAnnual_a0',
    model.1State.gamma.log.AR1.subAnnual_a0_a1_std  = 'model.1State.gamma.log.AR1.subAnnual_a0',

    model.1State.gamma.log.AR2.subAnnual_a0  = 'model.1State.gamma.log.AR1.subAnnual_a0',
    model.1State.gamma.log.AR2.subAnnual_a0_a1  = 'model.1State.gamma.log.AR2.subAnnual_a0',
    model.1State.gamma.log.AR2.subAnnual_a0_std  = 'model.1State.gamma.log.AR2.subAnnual_a0',
    model.1State.gamma.log.AR2.subAnnual_a0_a1_std  = 'model.1State.gamma.log.AR2.subAnnual_a0',

    model.1State.gamma.log.AR3.subAnnual_a0  = 'model.1State.gamma.log.AR2.subAnnual_a0',
    model.1State.gamma.log.AR3.subAnnual_a0_a1  = 'model.1State.gamma.log.AR3.subAnnual_a0',
    model.1State.gamma.log.AR3.subAnnual_a0_std  = 'model.1State.gamma.log.AR3.subAnnual_a0',
    model.1State.gamma.log.AR3.subAnnual_a0_a1_std  = 'model.1State.gamma.log.AR3.subAnnual_a0',

    # 2 state models
    model.2State.gamma.log.subAnnual_a0  = 'model.1State.gamma.log.subAnnual_a0',
    model.2State.gamma.log.subAnnual_a0_a1  = 'model.2State.gamma.log.subAnnual_a0',
    model.2State.gamma.log.subAnnual_a0_std  = 'model.2State.gamma.log.subAnnual_a0',
    model.2State.gamma.log.subAnnual_a0_a1_std  = 'model.2State.gamma.log.subAnnual_a0',

    model.2State.gamma.log.AR1.subAnnual_a0  = 'model.2State.gamma.log.subAnnual_a0',
    model.2State.gamma.log.AR1.subAnnual_a0_a1  = 'model.2State.gamma.log.AR1.subAnnual_a0',
    model.2State.gamma.log.AR1.subAnnual_a0_std  = 'model.2State.gamma.log.AR1.subAnnual_a0',
    model.2State.gamma.log.AR1.subAnnual_a0_a1_std  = 'model.2State.gamma.log.AR1.subAnnual_a0',

    model.2State.gamma.log.AR2.subAnnual_a0  = 'model.2State.gamma.log.AR1.subAnnual_a0',
    model.2State.gamma.log.AR2.subAnnual_a0_a1  = 'model.2State.gamma.log.AR2.subAnnual_a0',
    model.2State.gamma.log.AR2.subAnnual_a0_std  = 'model.2State.gamma.log.AR2.subAnnual_a0',
    model.2State.gamma.log.AR2.subAnnual_a0_a1_std  = 'model.2State.gamma.log.AR2.subAnnual_a0',

    model.2State.gamma.log.AR3.subAnnual_a0  = 'model.2State.gamma.log.AR2.subAnnual_a0',
    model.2State.gamma.log.AR3.subAnnual_a0_a1  = 'model.2State.gamma.log.AR3.subAnnual_a0',
    model.2State.gamma.log.AR3.subAnnual_a0_std  = 'model.2State.gamma.log.AR3.subAnnual_a0',
    model.2State.gamma.log.AR3.subAnnual_a0_a1_std  = 'model.2State.gamma.log.AR3.subAnnual_a0'
    )

    # 3 state models
    if (build.3state.models) {
      .Object@calib.reference.model.name = c( list(
        model.3State.gamma.log.subAnnual_a0  = 'model.2State.gamma.log.subAnnual_a0',
        model.3State.gamma.log.subAnnual_a0_a1  = 'model.3State.gamma.log.subAnnual_a0',
        model.3State.gamma.log.subAnnual_a0_std  = 'model.3State.gamma.log.subAnnual_a0',
        model.3State.gamma.log.subAnnual_a0_a1_std  = 'model.3State.gamma.log.subAnnual_a0',

        model.3State.gamma.log.AR1.subAnnual_a0  = 'model.3State.gamma.log.subAnnual_a0',
        model.3State.gamma.log.AR1.subAnnual_a0_a1  = 'model.3State.gamma.log.AR1.subAnnual_a0',
        model.3State.gamma.log.AR1.subAnnual_a0_std  = 'model.3State.gamma.log.AR1.subAnnual_a0',
        model.3State.gamma.log.AR1.subAnnual_a0_a1_std  = 'model.3State.gamma.log.AR1.subAnnual_a0',

        model.3State.gamma.log.AR2.subAnnual_a0  = 'model.3State.gamma.log.AR1.subAnnual_a0',
        model.3State.gamma.log.AR2.subAnnual_a0_a1  = 'model.3State.gamma.log.AR2.subAnnual_a0',
        model.3State.gamma.log.AR2.subAnnual_a0_std  = 'model.3State.gamma.log.AR2.subAnnual_a0',
        model.3State.gamma.log.AR2.subAnnual_a0_a1_std  = 'model.3State.gamma.log.AR2.subAnnual_a0',

        model.3State.gamma.log.AR3.subAnnual_a0  = 'model.3State.gamma.log.AR2.subAnnual_a0',
        model.3State.gamma.log.AR3.subAnnual_a0_a1  = 'model.3State.gamma.log.AR3.subAnnual_a0',
        model.3State.gamma.log.AR3.subAnnual_a0_std  = 'model.3State.gamma.log.AR3.subAnnual_a0',
        model.3State.gamma.log.AR3.subAnnual_a0_a1_std  = 'model.3State.gamma.log.AR3.subAnnual_a0',

        # 3 state unstructured models
        model.3StateUS.gamma.log.subAnnual_a0  = 'model.3State.gamma.log.subAnnual_a0',
        model.3StateUS.gamma.log.subAnnual_a0_a1  = 'model.3StateUS.gamma.log.subAnnual_a0',
        model.3StateUS.gamma.log.subAnnual_a0_std  = 'model.3StateUS.gamma.log.subAnnual_a0',
        model.3StateUS.gamma.log.subAnnual_a0_a1_std  = 'model.3StateUS.gamma.log.subAnnual_a0',

        model.3StateUS.gamma.log.AR1.subAnnual_a0  = 'model.3StateUS.gamma.log.subAnnual_a0',
        model.3StateUS.gamma.log.AR1.subAnnual_a0_a1  = 'model.3StateUS.gamma.log.AR1.subAnnual_a0',
        model.3StateUS.gamma.log.AR1.subAnnual_a0_std  = 'model.3StateUS.gamma.log.AR1.subAnnual_a0',
        model.3StateUS.gamma.log.AR1.subAnnual_a0_a1_std  = 'model.3StateUS.gamma.log.AR1.subAnnual_a0',

        model.3StateUS.gamma.log.AR2.subAnnual_a0  = 'model.3StateUS.gamma.log.AR1.subAnnual_a0',
        model.3StateUS.gamma.log.AR2.subAnnual_a0_a1  = 'model.3StateUS.gamma.log.AR2.subAnnual_a0',
        model.3StateUS.gamma.log.AR2.subAnnual_a0_std  = 'model.3StateUS.gamma.log.AR2.subAnnual_a0',
        model.3StateUS.gamma.log.AR2.subAnnual_a0_a1_std  = 'model.3StateUS.gamma.log.AR2.subAnnual_a0',

        model.3StateUS.gamma.log.AR3.subAnnual_a0  = 'model.3StateUS.gamma.log.AR2.subAnnual_a0',
        model.3StateUS.gamma.log.AR3.subAnnual_a0_a1  = 'model.3StateUS.gamma.log.AR3.subAnnual_a0',
        model.3StateUS.gamma.log.AR3.subAnnual_a0_std  = 'model.3StateUS.gamma.log.AR3.subAnnual_a0',
        model.3StateUS.gamma.log.AR3.subAnnual_a0_a1_std  = 'model.3StateUS.gamma.log.AR3.subAnnual_a0'
    ))
  }

  model.names = names(.Object@models)
  .Object@calib.reference.criteria.met = rep(F,length(model.names))
  names(.Object@calib.reference.criteria.met) <- model.names


  return(.Object)

}
)
