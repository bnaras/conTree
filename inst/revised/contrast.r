#
#
#                                CONTRAST
#
#                  Find x-regions p(y | x) != p(z | x)
#
#                               (4/26/20)
#
#             Coded and copyright (2020) by Jerome H. Friedman.
#
contrast=function(x,y,z,w=rep(1,nrow(x)),cat.vars=NULL,not.used=NULL,qint=10,
   xmiss=9.0e35,tree.size=10,min.node=500,mode='onesamp',type='dist',pwr=2,
   quant=0.5,nclass=NULL,costs=NULL,cdfsamp=500,verbose=FALSE,
   tree.store=1000000,cat.store=100000,nbump=1,fnodes=0.25,fsamp=1,
   doprint=FALSE) {
   cri='max'; qqtrim=20; n=nrow(x); crts=0;
   for (k in 1:nbump) {
      if(k==1) { s=1:n; yy=y}
      else { s=sample.int(n,as.integer(fsamp*n),replace=TRUE)
         if(is.vector(y)) {yy=y[s]} else { yy=y[s,]}
      }      
      trek=contrastt(x[s,],yy,z[s],w[s],cat.vars,not.used,qint,
         xmiss,tree.size,min.node,cri,mode,type,pwr,qqtrim,quant,nclass,costs,
         cdfsamp,verbose,tree.store,cat.store)
      v=nodesum(x,y,z,trek,doplot=FALSE)
      nf=as.integer(max(1,fnodes*length(v$nodes)))
      crt=sum(v$wt[1:nf]*v$cri[1:nf])/sum(v$wt[1:nf])
      if(doprint) print(c(k,crt))
      if (crt > crts) { crts=crt; trees=trek}
   }
   if(doprint) print(crts)
   invisible(trees)
}
contrastt=function(x,y,z,w=rep(1,nrow(x)),cat.vars=NULL,not.used=NULL,qint=10,
   xmiss=9.0e35,tree.size=10,min.node=500,cri='max',mode='onesamp',
   type='dist',pwr=2,qqtrim=20,quant=0.5,nclass=NULL,costs=NULL,cdfsamp=500,
   verbose=FALSE,tree.store=1000000,cat.store=100000) {
   if(!is.loaded('fcontrast'))
     stop("dyn.load('contrast.so')   linux\n  dyn.load('contrast.dll') windows")
   if(is.data.frame(x)) { p=length(x)}
   else if(is.matrix(x)) { p <- ncol(x)}
   else { stop('x must be a matrix or data frame.')}
   if(is.null(colnames(x))) { varnames=paste('V',as.character(1:p),sep='')}
   else { varnames=colnames(x)}
   if(is.data.frame(x)) { n=length(x[[1]])
      xx=matrix(nrow=n,ncol=p); for (j in 1:p) { xx[,j]=x[[j]]}
      if(is.null(cat.vars)) { lx=as.numeric(sapply(x,is.factor))+1}
      else { lx=rep(1,p)
         iv=getvars(cat.vars,p,lx,varnames); iv=iv[iv!=0]
         if(length(iv)!=0) lx[iv]=2
      }
   }
   else { n <- nrow(x);  xx=x; lx=rep(1,p)
      if(!is.null(cat.vars)) {
         iv=getvars(cat.vars,p,lx,varnames); iv=iv[iv!=0]
         if(length(iv)!=0) lx[iv]=2
      }
   }
   if(length(z)!=n) stop ('x and z dimensions inconsistent.')
   if(is.vector(y)) {
      if(length(y)!=n) stop ('x and y dimensions inconsistent.')   
      if(max(y)==min(y)) stop('all y values are the same.')
   }
   else {
      if(nrow(y)!=n) stop ('x and y dimensions inconsistent.')   
      if(max(y[,1])==min(y[,1])) stop('all y[,1] values are the same.')
      if(max(y[,2])==min(y[,2])) stop('all y[,2] values are the same.')
      if(ncol(y)!=2) stop('y must have two columns.')
   }
   if(length(w)!=n) stop ('x and w dimensions inconsistent.')
   if(!is.null(not.used)) {
      iv=getvars(not.used,p,rep(1,p),varnames); iv=iv[iv!=0]
      if(length(iv)!=0) lx[iv]=0
   }
   if (all(lx == 0)) stop('all predictor variables excluded.')
   xx[is.na(xx)] <- xmiss
   bgstx=max(xx[xx!=xmiss])
   if(xmiss<=bgstx)
      stop(paste('value of xmiss =',xmiss,
      'is smaller than largest training predictor variable value =',bgstx))
   if(any(is.na(y))) { w[is.na(y)] <- 0;
      warning
      ('training response contains NAs - corresponding weights set to 0.')
   }
   if(any(is.na(w))) { w[is.na(wtt)] <- 0;
      warning('training weights contain NAs - zeros substituted.')
   }
   if(any(w<0)) { w[w<0] <- 0;
      warning
      ('training weights contain negative numbers - zeros substituted.')
   }
   if(!is.character(cri)) stop(' cri must be of type chracter.')
   if(!is.character(mode)) stop(' mode must be of type chracter.')
   if(!is.character(type)) stop(' type must be of type chracter.')
   if(!is.null(costs)) {
      if(nrow(costs)!=nclass) stop('nrow(costs) incorrect.')
      if(ncol(costs)!=nclass) stop('ncol(costs) incorrect.')
   }
   tree.size=parchk("tree.size",tree.size,2,n,4)
   tree.store=parchk("tree.store",tree.store,10000,10000000,1000000)
   cat.store=parchk("cat.store",cat.store,10000,10000000,100000)
   qint=parchk('qint',qint,2,20,4)
   min.node=parchk('min.node',min.node,20,n/2,200)
   qqtrim=parchk('qqtrim',qqtrim,10,max(11,n/4),20)
   quant=parchk('quant',quant,0.01,0.99,0.5)
   cdfsamp=parchk('cdfsamp',cdfsamp,100,200000,500)
   if(!is.null(nclass)) nclass=parchk('nclass',nclass,2,100,2)
   v=contrast1(xx,y,z,w,lx,qint,xmiss,tree.size,min.node,cri,mode,type,
      pwr,qqtrim,quant,nclass,costs,cdfsamp,verbose,tree.store,cat.store)
   tree=list(itre=v$itre,rtre=v$rtre,cat=v$cat,kxt=v$kxt,kxc=v$kxc,p=v$p)
   parms=list(cri=cri,mode=mode,type=type,qqtrim=qqtrim,quant=quant,
      nclass=nclass,costs=costs,cdfsamp=cdfsamp,verbose=verbose)
   invisible(list(tree=tree,parms=parms))
}
getvars=function(vars,p,lx,names) {
   lv=length(vars); iv=rep(0,lv);    
   if(!is.character(vars)) {
      for (j in 1:lv) {
         if(vars[j]<1 || vars[j]>p || lx[vars[j]]==0) {
           stop(paste(vars[j],"is not one of the input variables."))
         }
         else { iv[j]=vars[j]}
      }
   }
   else {
      for (j in 1:lv) {  k=(1:p)[names==vars[j]];
         if (length(k)>0) { if(lx[k]>0) { iv[j]=k}}
         if(iv[j]==0) {
            stop(paste(names[lx>0],collapse=' '),'\n',vars[j],
               " is not one of the above input variables.")
         }
      }
   }
   iv
}
parchk = function (ax,x,lx,ux,df) {
   if(x < lx || x > ux) {
      warning(paste('invalid value for',ax,'- default (',df,') used.'))
      df
   }
   else { x}
}
contrast1=function(x,y,z,w,lx,qint,xmiss,tree.size,min.node,cri,
   mode,type,pwr,qqtrim,quant,nclass,costs,cdfsamp,verbose,
   tree.store,cat.store) {
   n=nrow(x); p=ncol(x);
   call=.Fortran('set_miss',arg=as.numeric(xmiss))
   call=.Fortran('set_trm',irg=as.integer(tree.size))
   call=.Fortran('set_ntn',irg=as.integer(min.node))
   call=.Fortran('set_qint',irg=as.integer(qint))
   call=.Fortran('set_pwr',arg=as.numeric(pwr))
   icri=1; if(cri!='max') icri=2
   call=.Fortran('set_cri',irg=as.integer(icri))
   if(mode=='onesamp') {
      kri=1; if(type=='qq') kri=2;
      if(type=='diff') kri=3; if(type=='class') kri=4;
      if(type=='quant') kri=5; if(type=='prob') kri=7;
      if(type=='maxmean') kri=8; if(type=='diffmean') kri=9
      if(kri==1) {
         if(is.matrix(y)) { kri=6
            call=.Fortran('set_samp',irg=as.integer(cdfsamp))
         }
      }
   }
   else {
      if(type=='dist') kri=10; if(type=='qq') kri=11
      if(type=='prob') kri=12; if(type=='diffmean') kri=13
      if(type=='maxmean') kri=14
      if(kri==10) {
         if(is.matrix(y)) { kri=15
            call=.Fortran('set_samp',irg=as.integer(cdfsamp))
         }
      }
   }
   call=.Fortran('set_kri',irg=as.integer(kri),jrg=as.integer(1))
   call=.Fortran('set_qqtrm',irg=as.integer(qqtrim),jrg=as.integer(1))
   call=.Fortran('set_quant',arg=as.numeric(quant))
   ivrb=0; if(verbose) ivrb=1
   call=.Fortran('set_vrb',irg=as.integer(ivrb),jrg=as.integer(1)) 
   if(kri==4) {
      if (is.null(nclass)) nclass=2
      if(is.null(costs)) {
         costs=matrix(rep(1,nclass*nclass),nrow=nclass,ncol=nclass)
         for(k in (1:nclass)) costs[k,k]=0
      }
      call=.Fortran('classin',ient=as.integer(1),nclasssv=as.integer(nclass),
         costssv=as.vector(as.numeric(costs)),nout=integer(1),out=numeric(1))
   }
   if(kri==6 | kri==15) { y2=y[,2]; y=y[,1]} else { y2=rep(0,n)}
   u=.Fortran('fcontrast',no=as.integer(n),ni=as.integer(p),
         x=as.vector(as.numeric(x)),y=as.vector(as.numeric(y)),
         y2=as.vector(as.numeric(y2)),z=as.vector(as.numeric(z)),
         w=as.vector(as.numeric(w)),lx=as.vector(as.integer(lx)),
         mxt=as.integer(tree.store),       
         itre=integer(6*tree.store),rtre=numeric(4*tree.store),
         mxc=as.integer(cat.store),cat=numeric(cat.store),
         ms=as.vector(as.integer(rep(0,2*n*p))),
         isc=as.vector(as.integer(rep(0,n))))
   v=.Fortran('get_stor',kxt=integer(1),kxc=integer(1))
   if (v$kxt > tree.store) stop('tree memory too small.')
   if (v$kxc > cat.store) stop('categorical memory too small.')
   invisible(list(itre=u$itre[1:(6*v$kxt)],rtre=u$rtre[1:(4*v$kxt)],
      cat=u$cat[1:v$kxc],kxt=v$kxt,kxc=v$kxc,p=p))
}
xcheck=function(x,xmiss=9.0e35) {
   if(is.data.frame(x)) { p=length(x); n=length(x[[1]])
      xx=matrix(nrow=n,ncol=p); for (j in 1:p) { xx[,j]=x[[j]]}
   }
   else if(is.matrix(x)) { n <- nrow(x); p <- ncol(x); xx=x;}
   else if(is.vector(x)) { p=length(x); n=1; xx=x}
   else { stop(' x must be a data frame, matrix, or vector.')}
   xx[is.na(xx)] <- xmiss
   call=.Fortran('set_miss',arg=as.numeric(xmiss))
   invisible(xx)
}
prune=function(tree,thr=0.1) {
   if(!is.loaded('prune1'))
      stop("dyn.load('contrast.so')   linux\n  dyn.load('contrast.dll') windows")
   v=tree$tree
   u=.Fortran('prune1',itr=as.vector(as.integer(v$itre)),
         rtr=as.vector(as.numeric(v$rtre)),
         nodes=as.integer(v$kxt),thr=as.numeric(thr),
         itro=integer(6*v$kxt),rtro=numeric(4*v$kxt))
    tr=list(itre=u$itro,rtre=u$rtro,cat=v$cat,
       kxt=v$kxt,kxc=v$kxc,p=v$p)
    invisible(list(tree=tr,parms=tree$parms))
}
crinode=function(tree){
   if(!is.loaded('crinode'))
      stop("dyn.load('contrast.so')   linux\n  dyn.load('contrast.dll') windows")
   u=tree$tree
   v=.Fortran('crinode',itr=as.vector(as.integer(u$itre)),
         rtr=as.vector(as.numeric(u$rtre)),
         mxnodes=as.integer(u$kxt),node=integer(1),nodes=integer(u$kxt),
         cri=numeric(u$kxt),wt=numeric(u$kxt))
   nodes=v$nodes[1:v$node]; cri=v$cri[1:v$node]; wt=v$wt[1:v$node]
   avecri=sum(wt*cri)/sum(wt)
   list(nodes=nodes,cri=cri,wt=wt,avecri=avecri)
}
nodesum=function(x,y,z,tree,w=rep(1,nrow(x)),doplot=FALSE) {
   u=crinode(tree); nds=u$nodes; nd=getnodes(x,tree)
   ln=length(u$nodes); crio=rep(0,ln); wt=rep(0,ln)
   for (j in 1:ln) {
     if (is.vector(y)) { yy=y[nd==nds[j]]} else {yy=y[nd==nds[j],]}
     v=getcri(yy,z[nd==nds[j]],tree,w[nd==nds[j]])
     crio[j]=abs(v$cri); wt[j]=v$wt
   }
   avecri=sum(wt*crio)/sum(wt)
   o=order(-crio)
   if(doplot) {
      opar=par(mfrow=c(2,1)); on.exit(par(opar))
      barplot(crio[o],names=nds[o],xlab='node',ylab='Criterion')
      barplot(wt[o],names=nds[o],xlab='node',ylab='Weight')
   }
   list(nodes=nds[o],cri=crio[o],wt=wt[o],avecri=avecri)
}
getnodes=function(x,tree){
   if(!is.loaded('getnodes1'))
      stop("dyn.load('contrast.so')   linux\n  dyn.load('contrast.dll') windows")
   x=xcheck(x)
   n=nrow(x); p=ncol(x);
   u=tree$tree
   v=.Fortran('getnodes1',no=as.integer(n),ni=as.integer(p),
         x=as.vector(as.numeric(x)),itre=as.vector(as.integer(u$itre)),
         rtre=as.vector(as.numeric(u$rtre)),cat=as.vector(as.numeric(u$cat)),
         nodes=integer(n))
    v$nodes
}
getlims=function(node,tree) {
   if(!is.loaded('getlims'))
      stop("dyn.load('contrast.so')   linux\n  dyn.load('contrast.dll') windows")
   u=tree$tree
   v=.Fortran('getlims',node=as.integer(node),ni=as.integer(u$p),
        itr=as.vector(as.integer(u$itre)),rtr=as.vector(as.numeric(u$rtre)),
        cat=as.vector(as.numeric(u$cat)),nvar=integer(1),jvar=integer(2*1000),
        vlims=numeric(1000),jerr=integer(1))
   if(v$jerr!=0) stop(paste('node',as.character(node),'not terminal.'))
   jvar=v$jvar[1:(2*v$nvar)]; jvar=matrix(jvar,nrow=2); vlims=v$vlims[1:v$nvar]
   if(v$nvar>1) { jvar=jvar[,v$nvar:1]; vlims=vlims[v$nvar:1]}
   list(nvar=v$nvar,jvar=jvar,vlims=vlims)
}
treesum=function(tree,nodes=NULL){
   q=tree$tree; v=crinode(tree);
   if(is.null(nodes)) { nodes=v$nodes}   
   for (k in 1:length(nodes)) {
      cat(paste('node',format(nodes[k],digits=0)))
      u=getlims(nodes[k],tree)      
      cat(paste('  var     dir    split'),'\n')      
      for (j in 1:u$nvar) {
         if (u$jvar[2,j] == 0) {
            if(sign(u$jvar[1,j])<0) {
               cat(paste('         ',format(abs(u$jvar[1,j]),digits=0),
               '     -     ',format(u$vlims[j],digits=2)),'\n')
            }
            else {
              cat(paste('         ',format(abs(u$jvar[1,j]),digits=0),
               '     +     ',format(u$vlims[j],digits=2)),'\n')
            }
         }
         else {
            kp=u$jvar[2,j]; kc=abs(q$cat[kp])
            z=sort(q$cat[(kp+1):(kp+kc)])
            if(u$vlims[j]>0) {
               cat(paste('      cat',format(u$jvar[1,j],digits=0),
                  '     -     ')); cat(z,'\n')
            }
            else {
               cat(paste('      cat',format(u$jvar[1,j],digits=0),
                  '     +     ')); cat(z,'\n')
            }
         }
      }
   }  
}
getcri=function(y,z,tree,w=rep(1,n),cdfsamp=500) {
   if(!is.loaded('andarm'))
      stop("dyn.load('contrast.so')   linux\n  dyn.load('contrast.dll') windows")
   if (is.vector(y)) { n=length(y)} else { n=nrow(y)}
   v=tree$parms; mode=v$mode; type=v$type
   if(mode=='onesamp') {
      kri=1; if(type=='qq') kri=2;
      if(type=='diff') kri=3; if(type=='class') kri=4;
      if(type=='quant') kri=5; if(type=='prob') kri=7;
      if(type=='maxmean') kri=8; if(type=='diffmean') kri=9
      if(kri==1) {
         if(is.matrix(y)) { kri=6
            call=.Fortran('set_samp',irg=as.integer(cdfsamp))
         }
      }
   }
   else {
      if(type=='dist') kri=10; if(type=='qq') kri=11
      if(type=='prob') kri=12; if(type=='diffmean') kri=13
      if(type=='maxmean') kri=14
      if(kri==10) {
         if(is.matrix(y)) { kri=15
            call=.Fortran('set_samp',irg=as.integer(cdfsamp))
         }
      }
   }
   call=.Fortran('set_kri',irg=as.integer(kri),jrg=as.integer(1))
   call=.Fortran('set_qqtrm',irg=as.integer(v$qqtrim),jrg=as.integer(1))
   call=.Fortran('set_quant',arg=as.numeric(v$quant))
   ivrb=0; if(v$verbose) ivrb=1
   call=.Fortran('set_vrb',irg=as.integer(ivrb),jrg=as.integer(1))
   if(kri==4) {
      if (is.null(v$nclass)) v$nclass=2
      if(is.null(v$costs)) {
         costs=matrix(rep(1,v$nclass*v$nclass),nrow=v$nclass,ncol=v$nclass)
         for(k in (1:v$nclass)) costs[k,k]=0
      }
      call=.Fortran('classin',ient=as.integer(1),nclasssv=as.integer(v$nclass),
         costssv=as.vector(as.numeric(costs)),nout=integer(1),out=numeric(1))
   }
   if(kri==6 | kri==15) { y2=y[,2]; y=y[,1]} else { y2=rep(0,n)}
   u=.Fortran('andarm',n=as.integer(n),y=as.vector(as.numeric(y)),
      y2=as.vector(as.numeric(y2)),z=as.vector(as.numeric(z)),
      w=as.vector(as.numeric(w)),dst=numeric(1),sw=numeric(1))
   list(cri=u$dst,wt=u$sw)
}
plotnodes=function(x,y,z,tree,w=rep(1,nrow(x)),nodes=NULL,
   pts='FALSE',xlim=NULL,ylim=NULL) {
   if(!is.loaded('fintcdf1'))
      stop("dyn.load('contrast.so')   linux\n  dyn.load('contrast.dll') windows")
   x=xcheck(x); mode=tree$parms$mode
   if(is.null(nodes)) { v=crinode(tree);
      nodes=v$nodes[1:min(length(v$nodes),9)]
   }
   nplts=length(nodes)
   nd=getnodes(x,tree)
   if (nplts > 2) { 
      nr=trunc(sqrt(nplts))
      if(nr*nr==nplts) { nc=nr}
      else { nr=nr+1; 
         nc=trunc(nplts/nr); if(nc*nr<nplts) nc=nc+1
      }
   }
   if (nplts==2) { nr=2; nc=1}
   if(nplts>=2) { opar=par(mfrow=c(nr,nc)); on.exit(par(opar))}
   for (k in 1:nplts) {
      if(is.vector(y)) {
         zn=z[nd==nodes[k]]; yn=y[nd==nodes[k]]
         if(mode=='onesamp') { p1=zn; p2=yn}
         else {p1=yn[zn<0]; p2=yn[zn>0]}
         if(is.null(xlim)) { xl=c(min(p1),max(p1))} else { xl=xlim}
         if(is.null(ylim)) { yl=c(min(p2),max(p2))} else { yl=ylim}
         qqplot(p1,p2,xlim=xl,ylim=yl,xlab='z',ylab='y',pch='.')
         title(paste('Node',as.character(nodes[k]),':',
            as.character(sum(w[nd==nodes[k]]))))
         lines(c(-1.0e9,1.0e9),c(-1.0e9,1.0e9),col='red')
      }
      else {
         u=y[nd==nodes[k],]; cdfsamp1=1000000
         if(length(unique(c(u[,1],u[,2])))>tree$parms$cdfsamp+2)
            cdfsamp1=tree$parms$cdfsamp
         vrb=0; if(tree$parms$verbose) vrb=1
         cy=cencdf(u,nsamp=cdfsamp1,vrb=vrb)
         u=z[nd==nodes[k]]; cdfsamp1=1000000
         cz=cencdf(cbind(u,u),nsamp=cdfsamp1,vrb=0)
         v=diffcdf(cy$x,cy$y,cz$x,cz$y,xlim=xlim,pts=pts)
         title(paste('Node',as.character(nodes[k]),':',
         as.character(format(v,digits=2))))
      }
   }
   invisible()
}
plotnodes2=function(x,y,z,tree,w=rep(1,nrow(x)),nodes=NULL,
   pts='FALSE',xlim=NULL,ylim=NULL) {
   if(!is.loaded('fintcdf1'))
      stop("dyn.load('contrast.so')   linux\n  dyn.load('contrast.dll') windows")
   x=xcheck(x)
   if(is.null(nodes)) { v=crinode(tree);
      nodes=v$nodes[1:min(length(v$nodes),9)]
   }
   nplts=length(nodes)
   nd=getnodes(x,tree)
   if (nplts > 2) { 
      nr=trunc(sqrt(nplts))
      if(nr*nr==nplts) { nc=nr}
      else { nr=nr+1; 
         nc=trunc(nplts/nr); if(nc*nr<nplts) nc=nc+1
      }
   }
   if (nplts==2) { nr=2; nc=1}
   if(nplts>=2) { opar=par(mfrow=c(nr,nc)); on.exit(par(opar))}
   vrb=0; if(tree$parms$verbose) vrb=1
   for (k in 1:nplts) {      
      u=y[nd==nodes[k],]; v=z[nd==nodes[k]]
      r=v<0; y1=u[r,]; y2=u[!r,]; cdfsamp1=1000000
      if(length(unique(c(y1[,1],y1[,2])))>tree$parms$cdfsamp+2)
         cdfsamp1=tree$parms$cdfsamp        
      cy1=cencdf(y1,nsamp=cdfsamp1,vrb=vrb)
      if(length(unique(c(y2[,1],y2[,2])))>tree$parms$cdfsamp+2)
         cdfsamp1=tree$parms$cdfsamp 
      cy2=cencdf(y2,nsamp=cdfsamp1,vrb=vrb)         
      v=diffcdf(cy1$x,cy1$y,cy2$x,cy2$y,xlim=xlim,pts=pts)
         title(paste('Node',as.character(nodes[k]),':',
         as.character(format(v,digits=2))))
   }
   invisible()
}
plotdiff=function(x,y,z,tree,w=rep(1,nrow(x)),
   nodes=NULL,xlim=NULL,ylim=NULL,span=0.15) {
   if(!is.loaded('fintcdf1'))
      stop("dyn.load('contrast.so')   linux\n  dyn.load('contrast.dll') windows")
   x=xcheck(x)
   if(is.null(xlim)) xlim=c(min(z),max(z))
   if(is.null(ylim)) ylim=c(min(y),max(y))
   if(is.null(nodes)) { v=crinode(tree)
      nodes=v$nodes[1:min(length(v$nodes),9)]
   }
   nplts=length(nodes)
   nd=getnodes(x,tree)
   if (nplts > 2) { 
      nr=trunc(sqrt(nplts))
      if(nr*nr==nplts) { nc=nr}
      else { nr=nr+1; 
         nc=trunc(nplts/nr); if(nc*nr<nplts) nc=nc+1
      }
   }
   if (nplts==2) { nr=2; nc=1}
   if(nplts>=2) { opar=par(mfrow=c(nr,nc)); on.exit(par(opar))}
   for (k in 1:nplts) {
         medplot(z[nd==nodes[k]],y[nd==nodes[k]],xlim=xlim,
            ylim=ylim,xlab='z',ylab='y',span=span)
         title(paste('Node',as.character(nodes[k]),':',
            as.character(sum(w[nd==nodes[k]]))))
         lines(c(-1.0e9,1.0e9),c(-1.0e9,1.0e9),col='blue')
   }
   invisible()
}
showclass=function(x,y,z,tree,w=rep(1,length(y))){
   opar=par(mfrow=c(2,1)); on.exit(par(opar))
   v=nodesum(x,y,z,tree,w); kk=length(v$nodes)
   u=barplot(v$cri,names=v$nodes,xlab='Node',ylab='Class  Risk')
   left=0.5*(u[2]-u[1]); right=0.5*(u[kk]-u[kk-1])
   lines(c(u[1]-left,u[kk]+right),c(v$avecri,v$avecri),col='red')
   barplot(v$wt,names=v$nodes,xlab='node',ylab='Counts')
}
showprobs=function(x,y,z,tree,w=rep(1,length(y)),vlab='Prob ( y = 1 )'){
   x=xcheck(x)
   opar=par(mfrow=c(2,1)); on.exit(par(opar))
   v=nodesum(x,y,z,tree,w); nd=getnodes(x,tree); kk=length(v$nodes)
   plt=matrix(nrow=2,ncol=kk)
   for (k in 1:kk) {
      if(tree$parms$mode=='onesamp') {
         sw=sum(w[nd==v$nodes[k]])
         plt[1,k]=sum(w[nd==v$nodes[k]]*y[nd==v$nodes[k]])/sw
         plt[2,k]=sum(w[nd==v$nodes[k]]*z[nd==v$nodes[k]])/sw
      }
      else { yn=y[nd==v$nodes[k]]; zn=z[nd==v$nodes[k]]; wn=w[nd==v$nodes[k]]
         y1=yn[zn<0]; y2=yn[zn>0]; w1=wn[zn<0]; w2=wn[zn>0];
         plt[1,k]=sum(w1*y1)/sum(w1)
         plt[2,k]=sum(w2*y2)/sum(w2)
      }
   }
   u=barplot(plt,names=v$nodes,xlab='Node',
      ylab=vlab,beside=TRUE,col=c('blue','red'))
   if(tree$parms$mode=='onesamp') {
      barplot(v$wt,names=v$nodes,xlab='node',ylab='Counts',col='green')
   }
   else {
      means=matrix(nrow=kk,ncol=2); nums=means; nodes=v$nodes
      for(k in 1:kk) nums[k,1]=sum(w[nd==nodes[k] & z<0])/sum(w[z<0])
      for(k in 1:kk) nums[k,2]=sum(w[nd==nodes[k] & z>0])/sum(w[z>0])
      barplot(t(nums),beside=T,col=c('blue','red'),
         names=nodes,xlab='Node',ylab='Frequency')
   }
}
prune.seq=function(tree,eps=0.01,plot.it=TRUE) {
   u=crinode(tree); del=eps*u$avecri; z=0;
   n0=length(u$nodes); n=n0
   nodes=rep(0,n); thr=rep(0,n); it=1
   thr[1]=0; nodes[1]=n;
   while(TRUE) {
      z=z+del
      tree=prune(tree,z); u=crinode(tree); n=length(u$nodes)
      if(n<n0) { it=it+1; nodes[it]=n; thr[it]=z}
      n0=n
      if(n<=2) break
   }
   if(plot.it) plot(nodes[1:it],thr[1:it],ylab='Threshold',xlab='Nodes')
   list(thr=thr[1:it],nodes=nodes[1:it])
}
showquants=function(x,y,z,tree,w=rep(1,length(y)),quant=0.5,doplot=TRUE){
   x=xcheck(x)
   opar=par(mfrow=c(2,1)); on.exit(par(opar))
   v=nodesum(x,y,z,tree,w); nd=getnodes(x,tree); kk=length(v$nodes)
   quants=rep(0,kk); num=rep(0,kk)
   for (k in 1:kk) { j=v$nodes[k]; num[k]=sum(nd==j)
      quants[k]=sum(y[nd==j]<=z[nd==j])/num[k]     
   }
   if (doplot) {
      u=barplot(quants-quant,names=v$nodes,xlab='node',ylab=' Probability  Error')
      title(paste(format(quant,digits=2),'- Quantile'))
      barplot(num,names=v$nodes,xlab='node',ylab='Counts')
   }
   sum(num*abs(quants-quant))/sum(num)
}
diffcdf=function(x1,cdf1,x2,cdf2,pts='FALSE',xlab='y',ylab='CDF ( y )',
   xlim=NULL,doplot=TRUE) {
   n1=length(x1); n2=length(x2); n=n1+n2
   z1=rep(0,n1); z2=rep(1,n2); x=c(x1,x2); lab=c(z1,z2)
   o=order(x); lab=lab[o]; diff=0; i1=0; i2=0;
   cdf11=0; cdf22=0
   for (i in 1:n) {
      if(lab[i]==0) { i1=i1+1; cdf11=cdf1[i1]
         diff=diff+abs(cdf11-cdf22);
      }
      else { i2=i2+1; cdf22=cdf2[i2]}
   }
   avediff=diff/i1
   if(doplot) {
      if(is.null(xlim)) xlim=c(min(x1,x2),max(x1,x2))
      if(pts) {
         plot(x1,cdf1,xlab=xlab,ylab=ylab,xlim=xlim,ylim=c(0,1))
      }
      else {
         plot(x1,cdf1,xlab=xlab,ylab=ylab,xlim=xlim,ylim=c(0,1),pch='.')
      }
      lines(x2,cdf2,col='red')
   }
   invisible(avediff)
}
cencdf=function(yin,win=rep(1,n),nsamp=2000,nit=100,thr=1.0e-2,vrb=0,
   xmiss=9.0e35) {
   oldseed=.Random.seed; set.seed(111)
   n=nrow(yin); r=sample.int(n,min(n,nsamp)); y=yin[r,]; w=win[r]
   .Random.seed=oldseed
   u=fintcdf(y,w,nit,thr,xmiss,vrb)
   t=yin[,1]>=yin[,2]; sw=sum(win); st=sum(win[t])/sw; snt=1-st
   b=c(yin[,1],yin[,2]); b=b[b>-xmiss]; b=b[b<xmiss]
   b=sort(unique(b))
   z=xfm2(b,u$x,u$y)
   if(st>0) { v=cdfpoints(b,sort(yin[t,1]),win[t])}
   else { v=0}
   invisible(list(x=b,y=snt*z+st*v))
}
fintcdf=function(y,w=rep(1,n),nit=100,thr=1.0e-2,xmiss=9.0e35,vrb=0) {
   b=c(y[,1],y[,2]); b=b[b>-xmiss]; b=b[b<xmiss]
   b=sort(unique(b)); b=c(b,xmiss); m=length(b)
   yy=y[y[,1]<y[,2],]
   vrb0=vrb
   call=.Fortran('set_vrb',irg=as.integer(vrb))
   z=.Fortran('fintcdf1',n=as.integer(nrow(yy)),y=as.vector(as.numeric(yy)),
     m=as.integer(m),b=as.vector(as.numeric(b)),w=as.vector(as.numeric(w)),
     nit=as.integer(nit),thr=as.numeric(thr/m),cdf=numeric(m),jt=integer(1),
     err=numeric(1))
  call=.Fortran('set_vrb',irg=as.integer(vrb0))
  if(vrb > 0) {
     cat('CDF calc:', format(z$jt,digits=3),'steps',
        ' confac =',format(z$err,digits=4),' ',format(m-1,digits=4),
        'points\n')
  }
  v=1:(m-1)
  invisible(list(x=b[v],y=z$cdf[v]))
}
cdfpoints=function(x,y,w=rep(1,length(y))) {
   z=.Fortran('cdfpoints1',m=as.integer(length(x)),x=as.vector(as.numeric(x)),
     n=as.integer(length(y)),y=as.vector(as.numeric(y)),
     w=as.vector(as.numeric(w)),cdf=numeric(length(x)))
     invisible(z$cdf)
}
xfm2=function(f,b00,b11,efac=7,xmiss=9.0e35) {
   b0=c(-xmiss,b00,xmiss); b1=rep(0,length(b0))
   nb=length(b0)-1; z=.bincode(f,b0); zp1=z+1
   b1[2:nb]=b11
   b0[1]=3*b0[2]-2*b0[3]; b1[1]=3*b1[2]-2*b1[3]
   b0[nb+1]=3*b0[nb]-2*b0[nb-1]; b1[nb+1]=3*b1[nb]-2*b1[nb-1]
   u=f>b0[1] & f<b0[nb+1]; f1=rep(0,length(f))
   if(sum(u)>0) {
      f1[u]=(f[u]-b0[z[u]])*(b1[zp1[u]]-b1[z[u]])/
         (b0[zp1[u]]-b0[z[u]]) + b1[z[u]]
   }
   u=f<=b0[1]
   if(sum(u)>0) {
      f1[u]=extrap(-efac,f[u],b0[z[u]],b0[zp1[u]],b1[z[u]],b1[zp1[u]])
   }
   u=f>=b0[nb+1]
   if(sum(u)>0) {
      f1[u]=extrap(efac,f[u],b0[z[u]],b0[zp1[u]],b1[z[u]],b1[zp1[u]])
   }
   f1
}
nodeplots=function(x,y,z,tree,w=rep(1,nrow(x)),nodes=NULL,
   xlim=NULL,ylim=NULL,pts='FALSE',span=0.15) {
   parms=tree$parms
   if(parms$type=='dist') {
      if(parms$mode=='onesamp' | is.vector(y)) {
         plotnodes(x,y,z,tree,w,nodes,pts,xlim,ylim)
      }
      else { plotnodes2(x,y,z,tree,w,nodes,pts,xlim,ylim)}
      return(invisible())
   }
   if(parms$type=='diff') {
      plotdiff(x,y,z,tree,w,nodes,xlim,ylim,span)
      return(invisible())
   }
   if(parms$type=='class') {
      showclass(x,y,z,tree,w); return(invisible())
   }
   if(parms$type=='prob') {
      showprobs(x,y,z,tree,w); return(invisible())
   }
   if(parms$type=='maxmean') {
      showprobs(x,y,z,tree,w,vlab='Mean'); return(invisible())
   }
   if(parms$type=='diffmean') {
      showprobs(x,y,z,tree,w,vlab='Mean'); return(invisible())
   }
   if(parms$type=='quant') {
      showquants(x,y,z,tree,w,quant=parms$quant,doplot=TRUE)
      return(invisible())
   }
   print('invalid type')
}
modtrast=function(x,y,z,w=rep(1,nrow(x)),cat.vars=NULL,not.used=NULL,qint=10,
   xmiss=9.0e35,tree.size=10,min.node=500,learn.rate=0.1,type='dist',pwr=2,
   quant=0.5,cdfsamp=500,verbose=FALSE,
   tree.store=1000000,cat.store=100000,nbump=1,fnodes=0.25,fsamp=1,
   doprint=FALSE,niter=100,doplot=TRUE,span=0,plot.span=0.15,print.itr=10) {
   if (type=='dist') {
      return(modtrans(x,y,z,w=rep(1,nrow(x)),cat.vars,not.used,qint,
         xmiss,tree.size,min.node,learn.rate,pwr,cdfsamp,verbose,tree.store,
         cat.store,nbump,fnodes,fsamp,doprint,niter,doplot,
         print.itr,span))
   }
   nodes=list(); dels=list(); trees=list(); r=z; acri=rep(0,niter); mode='onesamp'
   for (k in 1:niter) {
      trees[[k]]=contrast(x,y,r,w,cat.vars,not.used,qint,xmiss,tree.size,
         min.node,mode,type,pwr,quant,nclass=NULL,costs=NULL,cdfsamp,verbose,tree.store,
         cat.store,nbump,fnodes,fsamp,doprint)
      acri[k]=nodesum(x,y,r,trees[[k]],w)$avecri
      u=adjnodes(x,y,r,trees[[k]],w,learn.rate)
      r=u$az; dels[[k]]=u$del; nodes[[k]]=u$nodes
      if(k<=10 | k%%print.itr==0)
         cat(paste('Iteration',format(k,digits=3)),' ',format(acri[k],digits=3),"\n")
   }
   if(doplot) {
      if(niter>20 & plot.span>0) {
         medplot(1:k,acri,xlab='Iteration',ylab='Criterion',
            ylim=c(0,max(acri)),span=plot.span)
      }
      else {
         plot(1:k,acri,xlab='Iteration',ylab='Criterion',
            ylim=c(0,max(acri)))
       }
   }
   invisible(list(trees=trees,dels=dels,nodes=nodes,niter=niter))
}      
predtrast=function(x,z,model,num=model$niter) {
      if(model$trees[[1]]$parms$type=='dist') {
         return(predmod(x,z,model,num))
   }
   zo=z; t=model$trees[[1]]$parms$type=='prob'
   if(t) zo=log(zo/(1-zo))
   for(k in 1:num) {
      nd=getnodes(x,model$trees[[k]])
      for(j in 1:length(model$nodes[[k]])) {
         u=nd==model$nodes[[k]][j]
         if(t) { zo[u]=zo[u]+log(model$dels[[k]][j])}
         else { zo[u]=zo[u]+model$dels[[k]][j]}
      }
   }
   if(t) zo=1/(1+exp(-zo))
   zo
}
adjnodes=function(x,y,z,tree,w=rep(1,length(y)),learn.rate=1) {
   t=tree$parms$type
   if(!(t=='diff' | t=='quant' | t=='prob' | t=='maxmean' | t=='diffmean'))
      stop('incorrect tree type')
   u=nodesum(x,y,z,tree); nodes=length(u$nodes); r=z
   v=getnodes(x,tree); del=rep(0,nodes)
   if(t=='quant') { quant=tree$parms$quant}
   else {quant=0.5}
   for (j in 1:nodes) { k=u$nodes[j]
      if(t=='prob') {
         del[j]=(sum(w[v==k]*y[v==k])/sum(w[v==k]*z[v==k]))^learn.rate
         r[v==k]=z[v==k]*del[j]
      }
      else if(t=='maxmean' | t=='diffmean') {
         del[j]=learn.rate*sum(w[v==k]*(y[v==k])-z[v==k])/sum(w[v==k])
         r[v==k]=z[v==k]+del[j]
      }
      else { del[j]=learn.rate*as.numeric(quantile(y[v==k]-z[v==k],quant))
         r[v==k]=z[v==k]+del[j]

      }
   }
   invisible(list(nodes=u$nodes,del=del,az=r))
}
cvtrast=function(x,y,z,model,w=rep(1,length(y)),doplot=TRUE,span=0.5) {
   zo=z; t=model$trees[[1]]$parms$type=='prob'
   if(t) zo=log(zo/(1-zo)); niter=model$niter
   cvcri=rep(0,niter); h=z
   for(k in 1:niter) {
      cvcri[k]=nodesum(x,y,h,model$trees[[k]],w)$avecri
      nd=getnodes(x,model$trees[[k]])
      for(j in 1:length(model$nodes[[k]])) {
         u=nd==model$nodes[[k]][j]
         if(t) { zo[u]=zo[u]+log(model$dels[[k]][j])}
         else { zo[u]=zo[u]+model$dels[[k]][j]}
      }
      if(t) { h=1/(1+exp(-zo))} else { h=zo}      
   }
   if(doplot) {
      if(niter>20 & span>0) {
         scatter.smooth(1:k,cvcri,xlab='Iteration',ylab='Change',
         ylim=c(0,max(cvcri)),span=span)
      }
      else {
         plot(1:k,cvcri,xlab='Iteration',ylab='Change',
            ylim=c(0,max(cvcri)))
       }
   }
   invisible(cvcri)
}
medplot=function(x,y,np=NULL,main=NULL,span=0.15,lnln=FALSE,xlab='x',
   ylab='y',xlim=c(min(x),max(x)),ylim=c(min(y),max(y)),line=FALSE,
   col='red',doplot=TRUE) {
   o=order(x)
   spn=as.integer(length(x)*span); if (spn%%2!=1) spn=spn+1
   if (doplot) {
      p=1:length(x); if(!is.null(np)) p=sample(p,np)
      if (lnln) {
         if (line) {
            plot(x[o],runmed(y[o],spn),type='l',xlab=xlab,ylab=ylab,
              xlim=xlim,ylim=ylim,log='xy',col=col)
         }
         else {
            plot(x[p],y[p],pch='.',xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,log='xy')
         }
      }
      else {
         if (line) {
            plot(x[o],runmed(y[o],spn),xlab=xlab,ylab=ylab,xlim=xlim,
               ylim=ylim,type='l',col=col)
         }
         else { plot(x[p],y[p],pch='.',xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim)}      }
      if(!is.null(main)) title(main)
      if (!line) lines(x[o],runmed(y[o],spn),col=col)
   }
   invisible(list(x=x[o],y=y[o],sy=runmed(y[o],spn)))
}
quantloss=function(y,z,quant=0.5) {
   u=y>z; nu=!u; omq=1-quant
   v=as.numeric(quantile(y,quant))
   risk=omq*sum(abs(y[u]-z[u]))
      +quant*sum(abs(y[nu]-z[nu]))
   u=y>v; nu=!u;
   risk0=omq*sum(abs(y[u]-v))
      +quant*sum(abs(y[nu]-v))
   risk/risk0
}
devexp=function(y,p) {
   ll=mean(y*log(p)+(1-y)*log(1-p))
   u=mean(y); lo=mean(y*log(u)+(1-y)*log(1-u))
   1-(lo-ll)/lo
}
losserr=function(x,y,p,mdl,doplot=TRUE) {
   parms=mdl$tree[[1]]$parms
   z=rep(0,mdl$niter+1)
   if(parms$type=='prob') { z[1]=devexp(y,p)}
   else { z[1]=quantloss(y,p,parms$quant)}
   u=predtrast1(x,p,mdl)
   for (k in 2:(mdl$niter+1)) {
      if(parms$type=='prob') { z[k]=devexp(y,u[,k-1])}
      else { z[k]=quantloss(y,u[,k-1],parms$quant)}
   }
   if(doplot) {
       plot(0:mdl$niter,z,xlab='Iteration',ylab='Unexplained  devience')
   }
   invisible(z)
}
predtrast1=function(x,z,model,num=model$niter) {
   zo=z; t=model$trees[[1]]$parms$type=='prob'
   if(t) zo=log(zo/(1-zo)); zout=matrix(nrow=length(zo),ncol=num)
   for(k in 1:num) {
      nd=getnodes(x,model$trees[[k]])
      for(j in 1:length(model$nodes[[k]])) {
         u=nd==model$nodes[[k]][j]
         if(t) { zo[u]=zo[u]+log(model$dels[[k]][j])}
         else { zo[u]=zo[u]+model$dels[[k]][j]}
         zout[,k]=zo
      }
   }
   if(t) zout[,1:num]=1/(1+exp(-zout[,1:num]))
   invisible(zout)
}
xval =
function(x,y,z,mdl,num=length(mdl$tree),del=10,span=0.15,
   ylab='Average  Discrepancy',doplot=TRUE,doprint=TRUE) {
   parms=mdl$tree[[1]]$parms;
   error=rep(0,num); ntree=rep(0,num); k=0
   for(j in 1:num) {
      if (j==1 | j%%del==0) { k=k+1
         yp=predtrast(x,z,mdl,num=j)
         tree=contrast(x,y,yp,type=parms$type,quant=parms$quant)
         error[k]=crinode(tree)$avecri
         ntree[k]=j
         if(doprint) print(c(ntree[k],error[k]))
      }
   }
   if (doplot) {
      plot(ntree[1:k],error[1:k],xlab='Trees',ylab=ylab,ylim=c(0,max(error[1:k])))
      u=medplot(ntree[1:k],error[1:k],span=span,doplot=F)
      lines(u$x,u$sy,col='red')
   }
   invisible(list(ntree=ntree[1:k],error=error[1:k]))
}
trans=function(y,z,wy=rep(1,ny),wz=rep(1,nz),n=min(ny,nz)) {
   ny=length(y); nz=length(z); n=min(n,ny,nz)
   u=.Fortran('trans',ny=as.integer(ny),y=as.vector(as.numeric(y)),
      wy=as.vector(as.numeric(wy)),nz=as.integer(nz),z=as.vector(as.numeric(z)),
      wz=as.vector(as.numeric(wz)),nt=as.integer(n),t=numeric(2*n+4))
   invisible(list(x=u$t[1:(n+2)],y=u$t[(n+3):(2*n+4)]))
}
xfm=function(y,b0,b1) {
   f=approxfun(b0,b1,rule=2)
   invisible(f(y))
}
untie=function(y) {
   n=length(y)
   v=.Fortran('untie',n=as.integer(n),y=as.vector(as.numeric(y)),u=numeric(n))
   invisible(v$u)
}
modtrans=
function(x,y,z,w=rep(1,nrow(x)),cat.vars=NULL,not.used=NULL,qint=10,
   xmiss=9.0e35,tree.size=10,min.node=500,learn.rate=0.1,pwr=2,cdfsamp=500,verbose=FALSE,
   tree.store=1000000,cat.store=100000,nbump=1,fnodes=0.25,fsamp=1,
   doprint=FALSE,niter=100,doplot=TRUE,print.itr=10,span=0,plot.span=0.15) {
   nodes=list(); trans=list(); trees=list(); r=z; acri=rep(0,niter)
   l=0;
   for (k in 1:niter) {
      trees[[k]]=contrast(x,y,r,w,cat.vars,not.used,qint,xmiss,tree.size,
         min.node,mode='onesamp',type='dist',pwr,quant=0.5,nclass=NULL,costs=NULL,cdfsamp,verbose,tree.store,
         cat.store,nbump,fnodes,fsamp,doprint)
      v=nodesum(x,y,r,trees[[k]],w); nodes[[k]]=v$nodes; lnodes=length(v$nodes)
      acri[k]=v$avecri; nd=getnodes(x,trees[[k]])
      for (i in 1:lnodes) { j=v$nodes[i]; h=nd==j
         u=qqplot(r[h],y[h],plot.it=FALSE); l=l+1
         if(span>0 & span <1) { smo=loess.smooth(u$x,u$y,span=span)
            u$x=smo$x; u$y=smo$y
         }
         else if (span >= 1) { t=rank(u$x)/length(u$x);
            lsf=lsfit(u$x,u$y,t*(1-t))
            u$y=lsf$coefficients[1]+lsf$coefficients[2]*u$x
         }
         u$y=learn.rate*u$y+(1-learn.rate)*u$x
         trans[[l]]=cbind(u$x,u$y)
         r[h]=xfm(r[h],u$x,u$y)
      }
      if(k<10 | k%%print.itr==0)
         cat(paste('Iteration',format(k,digits=3)),' ',format(acri[k],digits=3),"\n")
   }
   if(doplot) {
      if(plot.span>0 & niter>20) {
         plot(1:k,acri,xlab='Iteration',ylab='Criterion',
            ylim=c(0,max(acri)),pch='.')
         u=medplot(1:k,acri,span=plot.span,doplot=F)
         lines(u$x,u$sy,col='red')
      }
      else {
         plot(1:k,acri,xlab='Iteration',ylab='Criterion',
            ylim=c(0,max(acri)))
       }
   }
   invisible(list(trees=trees,trans=trans,nodes=nodes,ty=r,niter=niter,cri=acri[1:k]))
}
predmod=function(x,z,model,num=model$niter) {
   r=z; l=0
   for(k in 1:num) {
      nd=getnodes(x,model$trees[[k]])
      nodes=model$nodes[[k]]
      for(i in 1:length(nodes)) { j=nodes[i]
         l=l+1; u=nd==j; if(sum(u)==0) next
         u1=model$trans[[l]][,1]; u2=model$trans[[l]][,2]
         f=approxfun(u1,u2,rule=2)     
         r[u]=f(r[u])         
      }
   }
   invisible(r)
}
ydist=function(x,z,model,num=model$niter) {
   x=xcheck(x)
   xr=matrix(nrow=length(z),ncol=length(x))
   for (i in 1:nrow(xr)) xr[i,]=x
   h=predmod(xr,z,model,num)
   invisible(h)
}
xvalmod=function(x,y,z,model,num=model$niter,doplot=TRUE) {
   r=z; l=0; cri=rep(0,num+1)
   cri[1]=crinode(contrast(x,y,r))$avecri
   for(k in 1:num) {
      nd=getnodes(x,model$trees[[k]])
      nodes=model$nodes[[k]]
      for(i in 1:length(nodes)) { j=nodes[i]
         l=l+1; u=nd==j; if(sum(u)==0) next      
         r[u]=xfm(r[u],model$trans[[l]][,1],model$trans[[l]][,2])         
      }
      cri[k+1]=crinode(contrast(x,y,r))$avecri
   }
   if(doplot) scatter.smooth(cri,span=0.25,xlab='Trees',ylab='Discrepancy')
   invisible(cri)
}
showdist=function(v,xtrim=NULL,xlim=NULL,xlab='Y | X',main=' ') {
  v=sort(v); v=v[v!=v[1] & v!=v[length(v)]]
  if(!is.null(xlim)) v=v[v>xlim[1] & v<xlim[2]]
  if(is.null(xlim)) xlim=c(min(v),max(v))
  hist(v,xlab=xlab,xlim=xlim,main=main)
  invisible(c(v[1],v[length(v)]))
}
plotpdf=function(x,z,model,num=model$niter,nq=200,span=0.25,xlim=NULL,xlab=NULL) {
   p=((1:nq)-0.5)/nq
   q=as.numeric(quantile(z,p))
   t=as.numeric(quantile(ydist(x,z,model,num),p))
   b=t!=t[1] & t!=t[length(t)]
   p=p[b]; t=t[b]; g=sum(b)
   d=(p[2:g]-p[1:(g-1)])/(t[2:g]-t[1:(g-1)])
   r=loess.smooth(t[2:g],d,type='l',span=span)
   if(is.null(xlim)) xlim=c(min(r$x),max(r$x))
   if(is.null(xlab)) xlab='Y'
   plot(r$x,r$y,type='l',xlim=xlim,xlab=xlab,ylab='PDF')   
   invisible()
}
plotcdf=function(x,z,model,num=model$niter,nq=100,span=0.25,xlim=NULL,xlab=NULL) {
   p=((1:nq)-0.5)/nq
   t=as.numeric(quantile(ydist(x,z,model,num),p))
   r=loess.smooth(t,p,type='l',span=span)
   if(is.null(xlim)) xlim=c(min(r$x),max(r$x))
   if(is.null(xlab)) xlab='Y'
   plot(r$x,r$y,type='l',xlim=xlim,xlab=xlab,ylab='CDF')   
   invisible()
}
plotdist=
function(x,z,model,num=model$niter,type='cdf',nq=100,
   span=0.25,pts=100,xlim=NULL,xlab=NULL,doplot=TRUE) {
   p=((1:nq)-0.5)/nq
   t=as.numeric(quantile(ydist(x,z,model,num),p))
   if (span >0) { r=loess.smooth(t,p,type='l',span=span,evaluation=pts)
      rx=r$x; ry=r$y
   }
   else { rx=t; ry=p}
   if(is.null(xlim)) xlim=c(min(rx),max(rx))
   if(is.null(xlab)) xlab='Y'
   if(type=='cdf') { px=rx; py=ry; ylim=c(0,1)}
   else { nx=length(rx); px=0.5*(rx[1:(nx-1)]+rx[2:nx])
      py=(ry[2:nx]-ry[1:(nx-1)])/(px[2]-px[1])
      ylim=c(0,max(py))
   }
   if(doplot) plot(px,py,type='l',xlim=xlim,ylim=ylim,xlab=xlab,ylab='CDF')
   invisible(list(x=px,y=py))
}
bootcri=
function(x,y,z,w=rep(1,nrow(x)),cat.vars=NULL,not.used=NULL,qint=10,
   xmiss=9.0e35,tree.size=10,min.node=500,mode='onesamp',type='dist',pwr=2,
   quant=0.5,nclass=NULL,costs=NULL,cdfsamp=500,verbose=FALSE,
   tree.store=1000000,cat.store=100000,nbump=100,fnodes=1,fsamp=1,
   doprint=FALSE) {
   cri='max'; qqtrim=20; n=nrow(x); crts=0
   if (length(fnodes)==1) { bcri=rep(0,nbump)}
   else { bcri=matrix(nrow=length(fnodes),ncol=nbump)}
   for (k in 1:nbump) {
      s=sample.int(n,as.integer(fsamp*n),replace=TRUE)
      if(is.vector(y)) {yy=y[s]} else { yy=y[s,]}
      trek=contrastt(x[s,],yy,z[s],w[s],cat.vars,not.used,qint,
         xmiss,tree.size,min.node,cri,mode,type,pwr,qqtrim,quant,nclass,costs,
         cdfsamp,verbose,tree.store,cat.store)
      v=nodesum(x[s,],yy,z[s],trek,doplot=FALSE)
      for (j in 1:length(fnodes)) {
         nf=as.integer(max(1,fnodes[j]*length(v$nodes)))
         crt=sum(v$wt[1:nf]*v$cri[1:nf])/sum(v$wt[1:nf])
         if(is.vector(bcri)) { bcri[k]=crt} else { bcri[j,k]=crt}
      }
      if(doprint) {
         if(is.vector(bcri)) { print(c(k,bcri[k]))}
         else { print(c(k,bcri[,k]))}
      }
   }
   if(is.vector(bcri)) { stds=sqrt(var(bcri))}
   else { stds=rep(0,length(fnodes))
      for(j in 1:length(fnodes)) stds[j]=sqrt(var(bcri[j,]))
   } 
   invisible(list(stds=stds,bcri=bcri))
}
lofcurve=function(x,y,z,tree,w=rep(1,length(y)),doplot=TRUE,ylim=NULL) {
   u=nodesum(x,y,z,tree,w); v=avesum(u$cri,u$wt)
   if(doplot) {
      if(is.null(ylim)) { yl=c(0,max(v$y))} else { yl=ylim}
      plot(v$x/length(y),v$y,ylim=yl,type='l',
         xlab='Fraction of Observations',ylab='Average  Discrepancy')
   }
   invisible(list(x=v$x/length(y),y=v$y))
}
avesum=
function(y,w=rep(1,n)){
   n=length(y); yout=rep(0,n); wout=rep(0,n); scri=0; swt=0
   for (k in 1:n) {
      scri=scri+y[k]*w[k]; swt=swt+w[k]
      yout[k]=scri/swt; wout[k]=swt
   }
   list(x=wout,y=yout)
}
extrap=function(efac,x,x1,x2,y1,y2) {
   if(efac<0){ x=pmax(x,x1+efac*(x2-x1))
      y1+(y2-y1)*(x-x1)/(x2-x1)
   }
   else { x=pmin(x,x2+efac*(x2-x1))
      y2+(y2-y1)*(x-x2)/(x2-x1)}
}

