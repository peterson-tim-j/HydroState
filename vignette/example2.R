# Run hydroState
#----------------
rm(list=ls())
library(hydroState)
library(DEoptim)
library(truncnorm)

# Load flow data
flowdata = read.csv('vignette/delwp_monthly.csv',stringsAsFactors=F)

# Extract oine catchment
gaugeID = 221201;
# gaugeID = 226205;
flowdata = flowdata[flowdata$gauge==gaugeID,]

# Convert to format for hydroState
flowdata = data.frame(year = flowdata$year, month=flowdata$month, flow=flowdata$q, precipitation=flowdata$p)

do.Monthly.Analysis=F

if (!do.Monthly.Analysis) {
  # Aggregate monthly data to seasons
  season.months = matrix(0,4,2)
  season.months[1,] = c(12,2)
  season.months[2,] = c(3,5)
  season.months[3,] = c(6,8)
  season.months[4,] = c(9,11)
  nrow=0
  season.end.month.current = 0
  #flowdata.seasonal = data.frame(year=c(), month=c(), flow=c(), precipitation=c(), nmonths=c())
  flowdata.seasonal = matrix(NA,nrow(flowdata),5)
  for (i in 1:nrow(flowdata)) {

    # Get season index.
    if (flowdata$month[i]>=season.months[1,1] || flowdata$month[i]<=season.months[1,2]) {
      season.ind = 1
    } else {
      season.ind = which(flowdata$month[i]>=season.months[,1] & flowdata$month[i]<=season.months[,2])
    }
    season.end.month = season.months[season.ind,2]

    if (season.end.month == season.end.month.current) {
      flowdata.seasonal[nrow,3] = flowdata.seasonal[nrow,3] + flowdata$flow[i]
      flowdata.seasonal[nrow,4] = flowdata.seasonal[nrow,4] + flowdata$precipitation[i]
      flowdata.seasonal[nrow,5] = flowdata.seasonal[nrow,5] + 1
    } else {
      nrow = nrow + 1
      if (season.end.month==2) {
        flowdata.seasonal[nrow,1] = flowdata$year[i]+1
      } else {
        flowdata.seasonal[nrow,1] = flowdata$year[i]
      }
      flowdata.seasonal[nrow,2] = season.end.month
      flowdata.seasonal[nrow,3] = flowdata$flow[i]
      flowdata.seasonal[nrow,4] = flowdata$precipitation[i]
      flowdata.seasonal[nrow,5] = 1

      season.end.month.current = season.end.month
    }
  }

  filt = !is.na(flowdata.seasonal[,1]) & flowdata.seasonal[,5]==3
  flowdata.seasonal = flowdata.seasonal[filt,]
  flowdata.seasonal = data.frame(year=flowdata.seasonal[,1], month=flowdata.seasonal[,2], flow=flowdata.seasonal[,3], precipitation=flowdata.seasonal[,4], nmonths=flowdata.seasonal[,5])

  # Remove the first 5 years if precip to minimise AR1 spin up effects. THIS IS A TEMP FIX!!!
  filt = flowdata.seasonal$year <= (min(flowdata.seasonal$year)+4)
  flowdata.seasonal$flow[filt] = NA

  flowdata = flowdata.seasonal
}

# Build input objects. Note a linear model and a model with first-roder serial correlation is built.
transition.graph=matrix(TRUE,2,2)
Qhat = new('Qhat.log', input.data=flowdata)
QhatModel = new('QhatModel.subAnnual.homo.gamma.linear.AR3', input.data=flowdata, transition.graph=transition.graph,
                state.dependent.mean.a0=T,state.dependent.mean.a1=F, state.dependent.std.a0=T,
                state.dependent.mean.AR1=F, state.dependent.mean.AR2=F, state.dependent.mean.AR3=F,
                subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T,subAnnual.dependent.std.a0=F)

markov = new('markov.annualHomogeneous.flickering', transition.graph=transition.graph, allow.flickering=T)

# Build HydroState object
model = new('hydroState',input.data=flowdata, Qhat.object=Qhat, QhatModel.object=QhatModel, markov.model.object=markov)

# Fit the model
model <- hydroState::fit(model,pop.size.perParameter = 50, parallelType=1, packages = c('hydroState','truncnorm'),parVar=c('model'))

# Name the states
model <- setStateNames(model, 1990)

# Plot Viterbi states
vstates<-viterbi(model,do.plot=T, plot.percentiles = c(0.05, 0.5, 0.95), plot.yearRange=c(1930,2017))

# Check reliability of viterbi statee predictions.
check.viterbi(model, 100000)

# Plot pseduo residuals
check.PseudoResiduals(model)

# Build all seasonal models
all.seasonalModels <- new('hydroState.subAnnual.allModels',as.character(gaugeID), flowdata, allow.flickering=T, build.3state.models=F)
all.seasonalModels <- fit(all.seasonalModels, pop.size.perParameter=10, max.generations=1000, doParallel=T)

