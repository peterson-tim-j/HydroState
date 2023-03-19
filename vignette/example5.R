# This HydroState examples builds and calibrates one two-state annual model and plots the results.
#----------------
rm(list=ls())

library(hydroState)
library(DEoptim)
library(truncnorm)

# Load flow data
data("test_monthly") #take 3 deep breathes, this takes 30 seconds to load

# Extract one catchment
gaugeID = 236209;
#test_daily = test_annual[test_annual$gauge==gaugeID,]

# Convert to format for hydroState #flow = TDS, precipitation = flow
test_monthly = data.frame(year = test_monthly$year, month = test_monthly$month, flow = test_monthly$Mean_TDS, precipitation=test_monthly$Mean_Q)
filt = test_monthly$year>=1900
test_monthly = test_monthly[filt,]

# Build input objects.
# Note a 2-state linear model with first-order serial correlation is built. A truncated normal distribution is used for
# each Markov state.
transition.graph=matrix(TRUE,2,2)
Qhat = new('Qhat.log', input.data=test_monthly) #Qhat.log is taking the log of streamflow
QhatModel = new('QhatModel.homo.normal.linear', input.data=test_monthly, transition.graph=transition.graph)
markov = new('markov.annualHomogeneous', transition.graph=transition.graph)

#view output from Chat.log
#getQhat(Chat)
#getQ.backTransformed(Chat,data = )
#getVariance(ChatModel)

# Build HydroState object
model = new('hydroState',input.data=test_monthly, Qhat.object=Qhat, QhatModel.object=QhatModel, markov.model.object=markov)

# Fit the model using DEoptim
model <- hydroState::fit(model,pop.size.perParameter = 10, max.generations=500)

# Name the states names with 1990 being defined as a 'norma' runoff year.
model <- setStateNames(model, 1996)

# Plot Viterbi states
viterbi(model)

# Plot pseduo residuals
check.PseudoResiduals(model)

# Build all seasonal models for this gauge, the claibration each model, select the best (byt AIC) and plot
#-------------------------------------------
all.Models <- new('hydroState.allModels',as.character(gaugeID), test_annual, allow.flickering=F)
all.Models <- fit(all.Models, pop.size.perParameter=5, max.generations=25, doParallel=F)
best.model = getAIC.bestModel(all.Models)
