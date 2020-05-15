#
#           Distribution contrast trees and boosting 
#                   on demographics data
# get data:
#
#   yage = outcome variable (age)
#   xage = predictor variables (other demographics)
#   gbage = training data gradient boosting model for median (yage | xage)
#
load('age_data.R')
#
# load library
#
library(contree)
#
# training and test data subsets
#
dl=1:5000 # training data
dt=5001:8856 # test data
#
# contrast y dist with global y marginal distribution on test data
#
set.seed(5)
zage=yage[sample.int(length(yage))]
treezage=contrast(xage[dt,],yage[dt],zage[dt],tree.size=9,min.node=200)
nodesum(xage[dt,],yage[dt],zage[dt],treezage)
nodeplots(xage[dt,],yage[dt],zage[dt],treezage)
#
# construct residual bootstrap data
#
res=yage-gbage
rbage=gbage+res[sample.int(length(res))]
#
# contrast y dist with residual bootstrap dist.
#
treerbage=contrast(xage[dt,],yage[dt],rbage,tree.size=9,min.node=200)
nodeplots(xage[dt,],yage[dt],rbage[dt],treerbage)
#
# contrast tree model on training data: residual boostrap start
#
mdlrb=modtrast(xage[dl,],yage[dl],rbage[dl])
#
# discrepancy on test data
#
xval(xage[dt,],yage[dt],rbage[dt],mdlrb)
#
# transform test data to solution
#
hrbage=predtrast(xage[dt,],rbage[dt],mdlrb)
#
# contrast y with transformed soln
# 
treehrbage=contrast(xage[dt,],yage[dt],hrbage,tree.size=9,min.node=200)
nodeplots(xage[dt,],yage[dt],hrbage,treerbage)
#
# lack-of-fit contrast curves:
# black = global y marginal distribution
# blue = residual bootstrap distribution 
# red = distribution boosting - residual bootstrap start
#
lofcurve(xage[dt,],yage[dt],zage[dt],treezage)
u=lofcurve(xage[dt,],yage[dt],rbage[dt],treerbage,doplot=F)
lines(u$x,u$y,col='blue')
u=lofcurve(xage[dt,],yage[dt],hrbage,treehrbage,doplot=F)
lines(u$x,u$y,col='red')
#
# CDFs y|x for seleted observations using transformed residual bootrap distribution
#
obs=c(8843,5716,7831,6505,4949,7555,3202,6048,7134)
p=((1:500)-.5)/500
qres=as.numeric(quantile(res,p))
par(mfrow=c(3,3))
for (k in 1:9) { plot(ydist(xage[obs[k],],gbage[obs[k]]+qres,mdlrb),p,type='l',xlim=c(13,100),
   xlab='Age',main=paste('Observation',as.character(obs[k])))
   points(yage[obs[k]],0,col='red')
}
#
# y|x distributions for same observations and model.
#
par(mfrow=c(3,3))
for (k in 1:9) {
   hist(ydist(xage[obs[k],],gbage[obs[k]]+qres,mdlrb),xlim=c(13,100),nclass=10,
   xlab='Age',main=paste('Observation',as.character(obs[k])))
   points(yage[obs[k]],0,col='red')
}
#
par(mfrow=c(1,1))