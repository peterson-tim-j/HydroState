# This HydroState examples builds and calibrates one two-state seasonal model and plots the results.
#----------------
rm(list=ls())
library(hydroState)
library(DEoptim)
library(truncnorm)

# Load flow data
data("streamflow_monthly")

# Extract oine catchment
gaugeID = 221201;
streamflow_monthly = streamflow_monthly[streamflow_monthly$gauge==gaugeID,]

# Convert to format for hydroState
streamflow_monthly = data.frame(year = streamflow_monthly$year, month=streamflow_monthly$month, flow=streamflow_monthly$q, precipitation=streamflow_monthly$p)

# Aggregate monthly data to seasonal.
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
  #streamflow_monthly.seasonal = data.frame(year=c(), month=c(), flow=c(), precipitation=c(), nmonths=c())
  streamflow_monthly.seasonal = matrix(NA,nrow(streamflow_monthly),5)
  for (i in 1:nrow(streamflow_monthly)) {

    # Get season index.
    if (streamflow_monthly$month[i]>=season.months[1,1] || streamflow_monthly$month[i]<=season.months[1,2]) {
      season.ind = 1
    } else {
      season.ind = which(streamflow_monthly$month[i]>=season.months[,1] & streamflow_monthly$month[i]<=season.months[,2])
    }
    season.end.month = season.months[season.ind,2]

    if (season.end.month == season.end.month.current) {
      streamflow_monthly.seasonal[nrow,3] = streamflow_monthly.seasonal[nrow,3] + streamflow_monthly$flow[i]
      streamflow_monthly.seasonal[nrow,4] = streamflow_monthly.seasonal[nrow,4] + streamflow_monthly$precipitation[i]
      streamflow_monthly.seasonal[nrow,5] = streamflow_monthly.seasonal[nrow,5] + 1
    } else {
      nrow = nrow + 1
      if (season.end.month==2) {
        streamflow_monthly.seasonal[nrow,1] = streamflow_monthly$year[i]+1
      } else {
        streamflow_monthly.seasonal[nrow,1] = streamflow_monthly$year[i]
      }
      streamflow_monthly.seasonal[nrow,2] = season.end.month
      streamflow_monthly.seasonal[nrow,3] = streamflow_monthly$flow[i]
      streamflow_monthly.seasonal[nrow,4] = streamflow_monthly$precipitation[i]
      streamflow_monthly.seasonal[nrow,5] = 1

      season.end.month.current = season.end.month
    }
  }

  filt = !is.na(streamflow_monthly.seasonal[,1]) & streamflow_monthly.seasonal[,5]==3
  streamflow_monthly.seasonal = streamflow_monthly.seasonal[filt,]
  streamflow_monthly.seasonal = data.frame(year=streamflow_monthly.seasonal[,1], month=streamflow_monthly.seasonal[,2], flow=streamflow_monthly.seasonal[,3], precipitation=streamflow_monthly.seasonal[,4], nmonths=streamflow_monthly.seasonal[,5])

  # Remove the first 5 years if precip to minimise AR1 spin up effects. THIS IS A TEMP FIX!!!
  filt = streamflow_monthly.seasonal$year <= (min(streamflow_monthly.seasonal$year)+4)
  streamflow_monthly.seasonal$flow[filt] = NA

  streamflow_monthly = streamflow_monthly.seasonal
}

# Build input objects. Note a linear model and a model with first-roder serial correlation is built.
transition.graph=matrix(TRUE,2,2)
Qhat = new('Qhat.log', input.data=streamflow_monthly)
QhatModel = new('QhatModel.subAnnual.homo.gamma.linear.AR3', input.data=streamflow_monthly, transition.graph=transition.graph,
                state.dependent.mean.a0=T,state.dependent.mean.a1=F, state.dependent.std.a0=T,
                state.dependent.mean.AR1=F, state.dependent.mean.AR2=F, state.dependent.mean.AR3=F,
                subAnnual.dependent.mean.a0=T, subAnnual.dependent.mean.a1=T,subAnnual.dependent.std.a0=F)

markov = new('markov.annualHomogeneous.flickering', transition.graph=transition.graph, allow.flickering=T)

# Build HydroState object
model = new('hydroState',input.data=streamflow_monthly, Qhat.object=Qhat, QhatModel.object=QhatModel, markov.model.object=markov)

# Fit the model
model <- hydroState::fit(model,pop.size.perParameter = 10, max.generations=500)

# Name the states names with 1990 being defined as a 'normal' runoff year.
model <- setStateNames(model, 1990)

# Plot Viterbi states
vstates<-viterbi(model)

# Plot pseduo residuals
check.PseudoResiduals(model)
