# This HydroState examples builds and calibrates many combinations of one and two state seasonal models
# and then extracts the best model (by AIC) and plots the results.
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

# Build all combinations of seasonal models for this gauge. Note, seasonal flickering between states is allowed.
all.Models <- new('hydroState.subAnnual.allModels',as.character(gaugeID), streamflow_monthly, allow.flickering=T, build.3state.models=F)

# Calibrate (using MLW) each of the models.
all.Models <- fit(all.Models, pop.size.perParameter=10, max.generations=500, doParallel=F)

# Select the best model (byt AIC)
best.model = getAIC.bestModel(all.Models)

# Name the states names with 1990 being defined as a 'norma' runoff year.
best.model <- setStateNames(best.model, 1990)

# Plot Viterbi states
viterbi(best.model)

# Plot pseduo residuals
check.PseudoResiduals(best.model)
