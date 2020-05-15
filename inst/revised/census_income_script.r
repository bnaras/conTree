
#
#                             Census Income Data
#                      Predict probability salary > $50000
#
# load data:
#   x = training data predictor variables (demographics, financial)
#   y = training data outcome variable (high salary indicator)
#   xt = test data predictor variables
#   yt = test data outcome variables
#   po = training data Gradient Boosting logistic probability estimates
#   pot = test data Gradient Boosting logistic probability estimates
#   prf = training data Random Forest probability estimates
#   prft = test data Random Forest probability estimates
#   psq = training data Gradient Boosting least-squares probability estimates
#   psqt = test data Gradient Boosting least-squares probability estimates
#
load('census_income_data.R')
#
# load contrast procedure
#
library(contree)
#
# test data subsets
#
dx = 1:10000
dxt = 10001:16281
#
# contrast salary indicator with Gradient Boosting logistic probability estimates
#
tree=contrast(xt[dx,],yt[dx],pot[dx],type='prob')
#
# terminal region statistics
#
nodesum(xt[dxt,],yt[dxt],pot[dxt],tree)
#
# terminal region plots
#
nodeplots(xt[dx,],yt[dx],pot[dx],tree)
#
# region boundries
#
treesum(tree,c(7,11))
#
# contrast salary indicator with Gradient Boosting squared-error direct probability estimates
#
tree=contrast(xt[dx,],yt[dx],psqt[dx],type='prob')
nodesum(xt[dxt,],yt[dxt],psqt[dxt],tree)
nodeplots(xt[dxt,],yt[dxt],psqt[dxt],tree)
#
# test data lack-of-fit contrast curves salary indicator with
#   black: Gradient Boosting logistic probability estimates
#   blue: Random Forest probability estimates
#   red: Gradient Boosting least-squares probability estimates
#
tree=contrast(xt[dx,],yt[dx],pot,type='prob',tree.size=50)
lofcurve(xt[dxt,],yt[dxt],pot[dxt],tree)
tree=contrast(xt[dx,],yt[dx],prft[dx,2],type='prob',tree.size=50)
u=lofcurve(xt[dxt,],yt[dxt],prft[dxt,2],tree,doplot=F)
lines(u$x,u$y,col='blue')
tree=contrast(xt[dx,],yt[dx],psqt[dx],type='prob',tree.size=50)
u=lofcurve(xt[dxt,],yt[dxt],psqt[dxt],tree,doplot=F)
lines(u$x,u$y,col='red')
#
# Contrast boost  probability estimates on TRAINING data
#
mdlo=modtrast(x,y,po,type='prob',niter=200)
#
# cross-validate model on test data
#
xval(xt,yt,pot,mdlo)
#
# test data predictions
#
hlot=predtrast(xt,pot,mdlo)
# 
# contrast boost Random Forest model
#
mdrf=modtrast(x,y,prf[,2],type='prob',niter=200)
hrft=predtrast(xt,prft[,2],mdrf)
#
# contrast boost Gradient Boosting squared-error direct model
#
mdsq=modtrast(x,y,psq,type='prob',niter=200)
hlsq=predtrast(xt,psqt,mdsq)
#
# test data lack-of-fit contrast curves salary indicator with 
# contrast boosting solutions
#
tree=contrast(xt[dx,],yt[dx],hrft[dx],type='prob',tree.size=50)
u=lofcurve(xt[dxt,],yt[dxt],hrft[dxt],tree,doplot=F)
plot(u$x,u$y,col='blue',type='l',ylim=c(0,0.04))
tree=contrast(xt[dx,],yt[dx],hlot,type='prob',tree.size=50)
u=lofcurve(xt[dxt,],yt[dxt],hlot[dxt],tree,doplot=F)
lines(u$x,u$y)
tree=contrast(xt[dx,],yt[dx],hlsq[dx],type='prob',tree.size=50)
u=lofcurve(xt[dxt,],yt[dxt],hlsq[dxt],tree,doplot=F)
lines(u$x,u$y,col='red')


 
