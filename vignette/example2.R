# Run hydroState
#----------------
rm(list=ls())
library(hydroState)
library(rgenoud)

# Load flow data
flowdata = read.csv('vignette/delwp_hy.csv')
flowdata = data.frame(flowdata)

# Extract oine catchment
gaugeID = 221201;
flowdata = flowdata[flowdata$gauge==gaugeID,]

# Convert to format for hydroState
flowdata = data.frame(year = flowdata$hy_year, flow=flowdata$q, precipitation=flowdata$p)

# Define transition graph for a three state model
transition.graph=matrix(TRUE,3,3)
transition.graph[1,3]=FALSE
transition.graph[2,1]=FALSE
transition.graph[3,2]=FALSE

# Build input object
Qhat = new('Qhat.boxcox', input.data=flowdata)
QhatModel = new('QhatModel.linear', input.data=flowdata, nStates=3)
markov = new('markov.annualHomogeneous', transition.graph=transition.graph)

# Build HydroState object
model_3state = new('hydroState',input.data=flowdata, Qhat.object=Qhat, QhatModel.object=QhatModel, markov.model.object=markov)

# Fit the model
model_3state <- hydroState::fit(model_3state, pop.size.perParameter = 1500, print.level = 1)

