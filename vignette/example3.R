# This HydroState examples builds and calibrates all combinations of one and two state annual models
# and then extracts the best model (by AIC) and plots the results.
#----------------
rm(list=ls())
library(hydroState)
library(DEoptim)
library(truncnorm)

# Load flow data
data(streamflow_annual)

# Extract one catchment
gaugeID = 221201;
streamflow_annual = streamflow_annual[streamflow_annual$gauge==gaugeID,]

# Convert to format for hydroState
streamflow_annual = data.frame(year = streamflow_annual$hy_year, flow=streamflow_annual$q, precipitation=streamflow_annual$p)
filt = streamflow_annual$year>=1900
streamflow_annual = streamflow_annual[filt,]

# Build all combinations of annual models for this gauge.
all.Models <- new('hydroState.allModels',as.character(gaugeID), streamflow_annual, allow.flickering=F)

# Calibrate (using MLW) each of the models.
all.Models <- fit(all.Models, pop.size.perParameter=10, max.generations=500, doParallel=F)

# Select the best model (byt AIC)
best.model = getAIC.bestModel(all.Models)

# Name the states names with 1990 being defined as a 'normal' runoff year.
best.model <- setStateNames(best.model, 1990)

# Plot Viterbi states
viterbi(best.model)

# Plot pseduo residuals
check.PseudoResiduals(best.model)
