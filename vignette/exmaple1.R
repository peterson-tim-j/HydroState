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
Qhatbar = new('Qhatbar.linearRegression', input.data=flowdata)
markov = new('markov.annualHomogeneous', transition.graph=matrix(TRUE,2,2))


# Check markov distn change w parameter chnage.
params.aslist = getParameters(markov)
params.aslist$mean=c(-150,50)
params.aslist$std=c(50,50)
markov <- setParameters(markov,params.aslist);
getTransitionProbMatrix(markov)

# Build HydroState object
model = new('hydroState',input.data=flowdata, Qhat.object=Qhat,Qhatbar.object=Qhatbar, markov.model.object=markov)


# params.aslist = getParameters(model)
# params.asvector = getParameters.asVector(model)
#
# getNegLogLikelihood(model,params.aslist)
#
# getNegLogLikelihood.fromTransformedVector(params.asvector,model)
#
# Domains=getTransformedParameterBounds.asVector(model)
# nvars=nrow(Domains)
# #genoud(hydroState::getNegLogLikelihood.fromTransformedVector,nvars=nvars,Domains = Domains, boundary.enforcement=0,.Object=model)
# est=genoud(hydroState::getNegLogLikelihood.fromTransformedVector,nvars=nvars,pop.size=100,max.generations=10^6,wait.generations=10,Domains = Domains,debug=T,print.level=3,.Object=model)
# est=genoud(hydroState::getNegLogLikelihood.fromTransformedVector,nvars=nvars,pop.size=10^6,max.generations=10^6,wait.generations=10,Domains = Domains,boundary.enforcement=0,.Object=model)

model <- hydroState::fit(model)

getParameters(model)

model@calibration.results
