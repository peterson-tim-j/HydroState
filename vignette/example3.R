# This HydroState examples builds and calibrates all combinations of one and two state annual models
# and then extracts the best model (by AIC) and plots the results.
#
# NOTE: the calibration settings on line 30 are very basic. In Peterson et al 2021, the following settings were used:
# pop.size.perParameter = 75,max.generations = 10000,reltol=1e-8,steptol=50.
# Uncomment line 31 (and comment out line 30) to use thse more robust settings.
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

# Calibrate each of the models.
# NOTE: comment out line 30 and uncomment line 31 to apply more robust calibration settings.
all.Models <- fit(all.Models, pop.size.perParameter=10, max.generations=500, doParallel=F)
#all.Models <- fit(all.Models,pop.size.perParameter = 75,max.generations = 10000,reltol=1e-8,steptol=50, doParallel=F)

# Select the best model (byt AIC)
all.Models = getAIC.bestModel(all.Models)

# Name the states names with 1990 being defined as a 'normal' runoff year.
# The years after 1990 are the options used if there is no flow data for 1990.
all.Models <- setStateNames(all.Models, c(1990, 1989, 1991, 1988, 1992))

# Plot Viterbi states
viterbi(all.Models)

# Plot pseduo residuals
check.PseudoResiduals(all.Models)
