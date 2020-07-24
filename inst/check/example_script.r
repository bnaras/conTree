#
# census_income_data
#
library(conTree)
load('census_income_data.Rdata')
dx = 1:10000; dxt = 10001:16281
tree=contrast(xt[dx,],yt[dx],gblt[dx],type='prob')
nodesum(xt[dxt,],yt[dxt],gblt[dxt],tree)
nodeplots(xt[dxt,],yt[dxt],gblt[dxt],tree)
treesum(tree,c(7,29))
nx=getnodes(xt,tree)
plot(gblt[nx==7],gbpt[nx == 7],pch='.',xlab='GBTL',ylab='GBPT')
lines(c(0,1),c(0,1),col='red')
tree=contrast(xt[dx,],yt[dx],gblt,type='prob',tree.size=50)
lofcurve(xt[dxt,],yt[dxt],gblt[dxt],tree)
tree=contrast(xt[dx,],yt[dx],rft[dx,2],type='prob',tree.size=50)
u=lofcurve(xt[dxt,],yt[dxt],rft[dxt,2],tree,doplot=F)
lines(u$x,u$y,col='blue')
tree=contrast(xt[dx,],yt[dx],gbpt[dx],type='prob',tree.size=50)
u=lofcurve(xt[dxt,],yt[dxt],gbpt[dxt],tree,doplot=F)
lines(u$x,u$y,col='red')
modgbl=modtrast(x,y,gbl,type='prob',niter=200)
xval(x,y,gbl,modgbl,col='red')
xval(xt,yt,gblt,modgbl,col='green',doplot='next')
hblt=predtrast(xt,gblt,modgbl)
tree=contrast(xt[dx,],yt[dx],hblt[dx],type='prob')
nodeplots(xt[dxt,],yt[dxt],hblt[dxt],tree)
#
# age_data
#
load('age_data.Rdata')
hist(yage)
dl=1:5000 # training data
dt=5001:8856 # test data
set.seed(5)
zage=yage[sample.int(length(yage))]
treezage=contrast(xage[dt,],yage[dt],zage[dt],tree.size=9,min.node=200)
nodesum(xage[dt,],yage[dt],zage[dt],treezage)
nodeplots(xage[dt,],yage[dt],zage[dt],treezage)
res=yage - gbage
rbage=gbage + res[sample.int(length(res))]
treerbage=contrast(xage[dt,],yage[dt],rbage[dt],tree.size=9,min.node=200)
nodeplots(xage[dt,],yage[dt],rbage[dt],treerbage)
mdlrb=modtrast(xage[dl,],yage[dl],rbage[dl],min.node=200)
xval(xage[dl,],yage[dl],rbage[dl],mdlrb,col='red')
xval(xage[dt,],yage[dt],rbage[dt],mdlrb,col='green',doplot='next')
hrbage=predtrast(xage[dt,],rbage[dt],mdlrb)
treehrbage=contrast(xage[dt,],yage[dt],hrbage,tree.size=9,min.node=200)
nodeplots(xage[dt,],yage[dt],hrbage,treerbage)
lofcurve(xage[dt,],yage[dt],zage[dt],treezage)
u=lofcurve(xage[dt,],yage[dt],rbage[dt],treerbage,doplot=F)
lines(u$x,u$y,col='blue')
u=lofcurve(xage[dt,],yage[dt],hrbage,treehrbage,doplot=F)
lines(u$x,u$y,col='red')
obs=c(8843,5716,7831,6505,4949,7555,3202,6048,7134)
p=((1:500)-0.5)/500; qres=as.numeric(quantile(res,p))
par(mfrow=c(3,3))
for (k in 1:9) {
 plot(ydist(xage[obs[k],],gbage[obs[k]]+qres,mdlrb),p,type='l',xlim=c(13,100),
 xlab='Age',main=paste('Observation',as.character(obs[k])))
 points(yage[obs[k]],0,col='red')
title(paste('Observation',obs[k]))
}
par(mfrow=c(3,3))
for (k in 1:9) {
 hist(ydist(xage[obs[k],],gbage[obs[k]]+qres,mdlrb),xlim=c(13,100),nclass=10,
 xlab='Age',main=paste('Observation',as.character(obs[k])))
 points(yage[obs[k]],0,col='red')
 title(paste('Observation',obs[k]))
}
par(mfrow=c(1,1))
#
# air_quality_data
#
load('air_quality_data.Rdata')
c(mean(yco[zco<0]),mean(yco[zco>0]))
dl=1:6000   # learning
dt=6001:9357 # test
tree=contrast(xco[dl,],yco[dl],zco[dl],mode='twosamp',type='diffmean')
nodesum(xco[dt,],yco[dt],zco[dt],tree)
hist(pr2,nclass=100)
wp=rep(0,length(zco))
wp[zco>0]=1/pr2[zco>0]
wp[zco<0]=1/(1-pr2[zco<0])
wp=length(wp)*wp/sum(wp)
tree=contrast(xco[dl,],yco[dl],zco[dl],w=wp[dl],mode='twosamp',type='diffmean')
nodesum(xco[dt,],yco[dt],zco[dt],tree,w=wp[dt])
nodeplots(xco[dt,],yco[dt],zco[dt],tree,w=wp[dt])
avedisc=rep(0,1000)
set.seed(13)
for (k in 1:1000) {
 zcot=zco[sample.int(length(zco))]
 tre=contrast(xco[dl,],yco[dl],zcot[dl],mode='twosamp',type='diffmean')
 avedisc[k]=nodesum(xco[dt,],yco[dt],zcot[dt],tre)$avecri
 if(k%%10==1) cat('.')
}
cat('\n')
hist(avedisc,xlim=c(0,4))
points(c(2.937557,3.790423),c(0,0),col=c('blue','red'))
treesum(tree,c(11,25))
ycodt=yco[dt]
qqplot(ycodt[zco[dt]<0],ycodt[zco[dt]>0],pch='.')
lines(c(0,200),c(0,200),col='red')
tree=contrast(xco[dl,],yco[dl],zco[dl],w=wp[dl],mode='twosamp',tree.size=9)
nodeplots(xco[dt,],yco[dt],zco[dt],w=wp[dt],tree)
