# This example demonstrates the parrallel calibraion of multiple models and
# of 1,2 and 3 states then the selelction of the best mode by AIC.

rm(list=ls())
library(hydroState)
library(parallel)
library(snow)
library(rgenoud)

# Load flow data
flowdata = read.csv('vignette/delwp_hy.csv')
flowdata = data.frame(flowdata)

# Extract oine catchment
gaugeID = 221201;
flowdata = flowdata[flowdata$gauge==gaugeID,]

# Convert to format for hydroState
flowdata = data.frame(year = flowdata$hy_year, flow=flowdata$q, precipitation=flowdata$p)

# Initialsie vector for 3 model objects
models =  vector('list',3)

# Define the one state model
transition.graph=matrix(TRUE,1,1)
Qhat = new('Qhat.boxcox', input.data=flowdata)
QhatModel = new('QhatModel.linear', input.data=flowdata, nStates=1)
markov = new('markov.annualHomogeneous', transition.graph=transition.graph)
models[[1]] <- new('hydroState',input.data=flowdata, Qhat.object=Qhat, QhatModel.object=QhatModel, markov.model.object=markov)

# Define the 2 state model
transition.graph=matrix(TRUE,3,3)
Qhat = new('Qhat.boxcox', input.data=flowdata)
QhatModel = new('QhatModel.linear', input.data=flowdata, nStates=2)
markov <- new('markov.annualHomogeneous', transition.graph=transition.graph)
models[[2]] <- new('hydroState',input.data=flowdata, Qhat.object=Qhat, QhatModel.object=QhatModel, markov.model.object=markov)

# Define a three state model
transition.graph=matrix(TRUE,3,3)
transition.graph[1,3]=FALSE
transition.graph[2,1]=FALSE
transition.graph[3,2]=FALSE
Qhat = new('Qhat.boxcox', input.data=flowdata)
QhatModel = new('QhatModel.linear', input.data=flowdata, nStates=3)
markov <- new('markov.annualHomogeneous', transition.graph=transition.graph)
models[[3]] <- new('hydroState',input.data=flowdata, Qhat.object=Qhat, QhatModel.object=QhatModel, markov.model.object=markov)

# Setup the local cluster for the calibration
nclus = detectCores(all.tests = TRUE, logical = FALSE)
clus <- c(rep("localhost", nclus))
local.cluster <- makeCluster(clus, type = "SOCK")
clusterEvalQ(local.cluster, library(rgenoud))

# Calibrate the three models in parrallel
models <- parLapply(local.cluster, models, function(lst) hydroState::fit(lst, pop.size.perParameter = 1000, gradient.check=F, print.level = 0))
stopCluster(local.cluster)

# Get the Akaike information criterion for each model
models.AIC <- lapply(models,"getAIC")
models.AIC <- unlist(models.AIC)
