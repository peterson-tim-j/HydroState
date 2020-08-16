# Run hydroState
#----------------
rm(list=ls())
library(hydroState)
library(DEoptim)
library(truncnorm)

# Load flow data
flowdata = read.csv('vignette/delwp_hy.csv')
flowdata = data.frame(flowdata)

# Extract one catchment
gaugeID = 221201;
flowdata = flowdata[flowdata$gauge==gaugeID,]

# Convert to format for hydroState
flowdata = data.frame(year = flowdata$hy_year, flow=flowdata$q, precipitation=flowdata$p)
filt = flowdata$year>=1900
flowdata = flowdata[filt,]

# flowdata = flowdata[104:nrow(flowdata),]

# Build input objects. Note a linear model and a model with first-roder serial correlation is built.
transition.graph=matrix(TRUE,2,2)
Qhat = new('Qhat.log', input.data=flowdata)
QhatModel = new('QhatModel.homo.normal.linear.AR3', input.data=flowdata, transition.graph=transition.graph)
markov = new('markov.annualHomogeneous', transition.graph=transition.graph)

# Build HydroState object
model = new('hydroState',input.data=flowdata, Qhat.object=Qhat, QhatModel.object=QhatModel, markov.model.object=markov)

# Fit the model
t1 = Sys.time()
model <- hydroState::fit(model,pop.size.perParameter = 10, max.generations=100)
t2 = Sys.time()
t2-t1

# Name the states
model <- setStateNames(model, 1990)

# Plot Viterbi states
viterbi(model)


# Check reliability of viterbi statee predictions.
check.viterbi(model, 100000)

# Plot pseduo residuals
check.PseudoResiduals(model)

# Get resilience index.
drought.resilience.index(model, year.drought.start=1998, year.drought.end=2009,year.postdrought.end=2016)
