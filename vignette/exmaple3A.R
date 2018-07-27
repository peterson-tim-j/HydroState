# Identify the appropriate rgenoud population size for a 3 stat model

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

# Build input object
Qhat = new('Qhat.boxcox', input.data=flowdata)
QhatModel = new('QhatModel.linear', input.data=flowdata, nStates=2)
markov = new('markov.annualHomogeneous', transition.graph=matrix(TRUE,2,2))

# Build HydroState object
model = new('hydroState',input.data=flowdata, Qhat.object=Qhat, QhatModel.object=QhatModel, markov.model.object=markov)

# Set the number of trials per increment in the population size.
nsamples = 10

# Set the population size increments to trial
pop.size.perParameter = seq(100, 4000,by=100)

# Initialse a vector of all models
models =  vector('list',nsamples * length(pop.size.perParameter))

# Setup the local cluster for the calibration
nclus = detectCores(all.tests = TRUE, logical = FALSE)
clus <- c(rep("localhost", nclus))
local.cluster <- makeCluster(clus, type = "SOCK")
clusterEvalQ(local.cluster, library(rgenoud))

# Calibrate each trial for each increment
negLL = matrix(Inf,length(pop.size.perParameter),nsamples)
nGens = matrix(Inf,length(pop.size.perParameter),nsamples)
for (i in 1:length(pop.size.perParameter)) {
  ind = ((i-1)*pop.size.perParameter):(i*pop.size.perParameter)
  for (j in 1:nsamples) {
    models[[ind[j]]] = new('hydroState',input.data=flowdata, Qhat.object=Qhat, QhatModel.object=QhatModel, markov.model.object=markov)
  }

  # Calibrate nsamples models in parrallel
  models[[ind]] <- parLapply(local.cluster, models[[ind]], function(lst) hydroState::fit(lst, pop.size.perParameter = pop.size.perParameter[i], print.level = 0))

  for (j in 1:nsamples) {
    negLL[i,j] = models[[ind[j]]]@calibration.results$value
    nGens[i,j] = models[[ind[j]]]@calibration.results$generations
  }
}
stopCluster(local.cluster)
