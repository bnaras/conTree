"
c
c                                CONTRAST
c
c                  Find x-regions p(y | x) != p(z | x)
c                     
c                                (5/21/20)
c
c             Coded and copyright (2020) by Jerome H. Friedman.
c             Modified for R package; search for occurences of: Naras
c
c Build tree:
c
c   call fcontrast (no,ni,x,y,z,w,lx,mxt,itre,rtre,mxc,cat,ms,isc)
c
c input:
c   no = number of training observations (cases)
c   ni = number of input variables (features/attributes)
c   x(no,ni) = input variable training data matrix
c   y(no) = training data y-output values
c   z(no) = contrast data z-output values
c   w(no) = training data observation weights
c   lx(ni) = input variable flags:
c              0 => ignore jth input variable
c      lx(j) = 1 => jth input variable is orderable (real)
c              2 => jth input variable is unorderable (categorical)
c   mxt = 2nd dimension of itr and rtr arrays (set as large as possible)
c   mxc = dimension of cat array (set as large as possible)
c         (see Auxillary output - below)
c
c output:
c   itre(6,mxt) = integer tree storage
c   rtre(4,mxt) = real valued tree storage
c   cat(mxc) = categorical variable split storage
c
c scratch workspace:
c
c   ms(no,ni,2) : integer
c   isc(no) : integer
c
c
c Bottom-up tree pruning:
c
c   call prune (itre,rtre,mxt,thr)
c 
c   input:
c     itre(6,mxt),rtre(4,mxt)  = output from omnilog
c     mxt = input to omnilog
c     thr = split remove pruning threshold
c   
c   Output:
c     itre(6,mxt),rtre(4,mxt)  = resulting imnilog pruned tree.
c
c
c Terminal node (region) summary:
c
c   call crinode (itre,rtre,mxt,nunnode,nodes,cri,wt)
c
c   input:
c     itre(6,mxt),rtre(4,mxt)  = output from omnilog
c     mxt = input to omnilog
c
c   output:
c      numnode = number of terminal nodes
c      nodes(numnode) = terminal node identity
c      cri(numnodes)) = lack-of-fit criterion value for each terminal node
c      wt(numnodes)) = sum of weights for each terminal node
c  
c Identify terminal nodes containing specified observations (x-values):
c
c   call getnodes (no,ni,x,itre,rtre,cat,nodes)
c
c   input:
c      no = number of observations (x-values)
c      ni = number of predictor variables (input to omnilog)
c      x(no,ni) = specified observation x-values
c      itre,rtre,cat = output from omnilog
c
c   output:
c      nodes(no) = respective terminal node identities 
c      
c Obtain predictor variable boundaries for specified terminal nodes:
c
c   call getlims(node,ni,itr,rtr,cat,nvar,jvar,vlims,jerr)
c
c   input:
c      node = terminal node identity (output from getnodes)
c      ni = number of predictor variables (input to omnilog)
c      itr,rtr,cat = output from ominlog
c
c   output:
c      nvar = number of variables split on path to node
c      jvar(nvar,2) = variable identities  and type
c         jvar(j,2)=0 => jvar(1,j) = numeric splitting variable
c                       (<0 => upper bound, >0 => lower bound)
c         jvar(j,2)!=0 => categorical splitting variable pointer to cat array
c      vlims(nvar) = split point (numeric split)
c                  = pointer to cat array (categorical split)
c
c
c Auxillary output (values can be obtained with following call):
c
c   call get_stor (kxt,kxc): actual storage used.
c      kxt = minimum value required for 2nd dimension of itr and rtr arrays
c      kxc = minimum value required for dimension of cat array
c      [if these values are larger than (or close to) the corresponding input
c       values (mxt and mxc - see above) then a storage error may have
c       occurred.]
c
c
c
c Additional procedure parameters (current value can be changed with
c   the following calls):
c
c
c   call set_trm (nterm): maximum tree size
c      nterm = maximum number of terminal nodes in each tree (default = 10)
c
c
c   call set_miss (xmiss): missing value flag
c      All missing input values must be set to the value xmiss, which must
c      be larger than any possible non-missing input data value.
c      (default = 9.0e30)
c
c
c   call set_ntn (ntn): terminal node size parameter
c      ntn = minimum number of (non-zero weight) training observations in
c            each terminal node (default =500)
c
c   call set_cri (icri): predictor variable selection criterion
c      icri=1 => maximum of value on both sides (default)
c      icri=2 => absolute difference between values on both sides
c
c   call set_kri (kri): splitting criterion
c      kri = 1 => weighted CDF difference criterion (default))
c      kri = 2 => quantile difference criterion
c      kri = 3 => absolute value difference criterion
c      kri = 4 => classification predicted difference
c      kri = 5 => quantile predicted difference
c
c   call set_qqtrm (itrim): quantile difference trimming
c      itrim = number of extreme valued observations trimmed
c         (default = 20, kri=2 only)
c
c
c Description of arrays:
c
c   ms(no,ni) = permutation vector sorting each input variable
c               in ascending order for each terminal node.
c   itre(6,mxt): integer tree storage
c       (1,.) = split variable ( < 0 => categorical)
c       (2,.) = left son
c       (3,.) = right son
c       (4,.) = parent ( < 0 => terminal node)
c       (5,.) = points to first node obs in ms array
c       (6,.) = points to last node obs in ms array
c
c   rtre(4,mxt): real tree storage
c       (1,.) = split point [itre(1,.) < 0 => pointer (kp) to cat array]
c       (2,.) = split improvement ( < 0 => missing value pre-split)
c       (3,.) = lack-of-fit value in node
c       (4,.) = sum w(i), i in node
c
c   cat(mxc) = categorical split storage
c      (kp) = abs-value => number of cat values defining split
c      (kp+1 ... kp + abs(cat(kp))) = set of values defining split
c         cat(kp) < 0 => go left if x in set of values, else go right
c         cat(kp) > 0 => go right if x in set of values, else go left
c
c
"
subroutine fcontrast (no,ni,x,y,y2,z,w,lx,mxt,itre,rtre,mxc,cat,ms,isc);
%fortran
!DEC$ ATTRIBUTES DLLEXPORT :: fcontrast
%mortran
integer lx(ni),ms(no,ni,2),itre(6,mxt),isc(no);
real x(no,ni),y(no),y2(no),z(no),w(no),rtre(4,mxt),cat(mxc);
save nodes,nct,kstor,lstor;
data mxtrm /10/;
call dosort(no,ni,x,w,ms,nu,isc);
lstor=0;
itre(4,1)=0; itre(5,1)=1; itre(6,1)=nu; itre(2,1)=2; itre(3,1)=3;
call andarm(no,y,y2,z,w,rtre(2,1),rtre(4,1));
nct=1;
call split7(no,ni,x,y,y2,z,w,lx,ms,1,nu,itre(1,1),rtre(1,1),
   cri1,cri2,w1,w2,rtre(3,1),kct,cat(nct));
if itre(1,1).eq.0 < itre(4,1)=-9999; return;>
if itre(1,1).lt.0 < rtre(1,1)=nct; nct=nct+kct; if(nct.gt.mxc) lstor=1;>
/rtre(2,2),rtre(2,3)/=0.0; rtre(3,2)=cri1; rtre(3,3)=cri2;
/itre(4,2),itre(4,3)/=-1; rtre(4,2)=w1; rtre(4,3)=w2;
rtre(2,1)=sign(max(0.0,rtre(3,1)-rtre(2,1)),rtre(3,1));
if(lstor.ne.0) return;
nodes=3; ntrm=0; if(rtre(2,1).gt.0.0) ntrm=1; ktrg=1; kstor=0;
until ntrm.ge.mxtrm < jt=itre(1,ktrg); st=rtre(1,ktrg);
   k5=itre(5,ktrg); k6=itre(6,ktrg);
"Naras fix: explicit conversion to integer"
"   if jt.lt.0 < ju=-jt; kp=st+0.1; np=abs(cat(kp))+0.1;>"
if jt.lt.0 < ju=-jt; kp=int(st+0.1); np=int(abs(cat(kp))+0.1);>
   <j=1,ni; if(ni.gt.1.and.j.eq.jt) next; kl=k5-1; kr=k6+1;
      <i=k5,k6; l=ms(i,j,1);
         if jt.gt.0 <
            if x(l,jt).lt.st < kl=kl+1; isc(kl)=l;>
            else < kr=kr-1; isc(kr)=l;>
         >
         else < in=0;
            <im=1,np; if(x(l,ju).ne.cat(kp+im)) next; in=1; exit;>
            if cat(kp).gt.0.0 <
               if in.eq.0 < kl=kl+1; isc(kl)=l;>
               else < kr=kr-1; isc(kr)=l;>
            >
            else <
               if in.ne.0 < kl=kl+1; isc(kl)=l;>
               else < kr=kr-1; isc(kr)=l;>
            >
         >
      >
      <i=k5,kl;  ms(i,j,1)=isc(i);> <i=kr,k6; ms(i,j,1)=isc(k6-i+kr);>
   >
   itre(5,itre(2,ktrg))=k5; itre(6,itre(2,ktrg))=kl;
   itre(5,itre(3,ktrg))=kr; itre(6,itre(3,ktrg))=k6;
   nde=itre(2,ktrg);
   call split7(no,ni,x,y,y2,z,w,lx,ms,itre(5,nde),itre(6,nde),itre(1,nde),
      rtre(1,nde),cri1,cri2,w1,w2,cri,kct,cat(nct));
   if itre(1,nde).lt.0 < rtre(1,nde)=nct; nct=nct+kct; if(nct.gt.mxc) lstor=1;>
   if nodes+2.gt.mxt < kstor=1;>
   else < itre(2,nde)=nodes+1; itre(3,nde)=nodes+2;
      l=nodes+1; /rtre(2,l),rtre(2,l+1)/=0.0; rtre(3,l)=cri1; rtre(3,l+1)=cri2;
      rtre(4,l)=w1; rtre(4,l+1)=w2; /itre(4,l),itre(4,l+1)/=-nde;
      rtre(2,nde)=sign(max(0.0,cri-rtre(3,nde)),cri);
   >
   nde=itre(3,ktrg);
   call split7(no,ni,x,y,y2,z,w,lx,ms,itre(5,nde),itre(6,nde),itre(1,nde),
      rtre(1,nde),cri1,cri2,w1,w2,cri,kct,cat(nct));
   if itre(1,nde).lt.0 < rtre(1,nde)=nct; nct=nct+kct; if(nct.gt.mxc) lstor=1;>
   if nodes+4.gt.mxt < kstor=1;>
   else < itre(2,nde)=nodes+3; itre(3,nde)=nodes+4;
      l=nodes+3; /rtre(2,l),rtre(2,l+1)/=0.0; rtre(3,l)=cri1; rtre(3,l+1)=cri2;
      rtre(4,l)=w1; rtre(4,l+1)=w2; /itre(4,l),itre(4,l+1)/=-nde;
      rtre(2,nde)=sign(max(0.0,cri-rtre(3,nde)),cri);
      itre(4,ktrg)=-itre(4,ktrg); rsv=rtre(2,ktrg); crix=0.0;     
      <k=1,nodes; if(itre(4,k).ge.0) next; if(abs(rtre(2,k)).le.crix) next;
         if(itre(1,k).eq.0) next; crix=abs(rtre(2,k)); ktrg=k;
      >
   >
   if(crix.le.0.0.or.kstor.ne.0.or.lstor.ne.0) return;
   nodes=nodes+4; if(rsv.gt.0.0) ntrm=ntrm+1;
>
return;
entry set_trm(irg);
%fortran
!DEC$ ATTRIBUTES DLLEXPORT :: set_trm
%mortran
mxtrm=max0(irg,2); return;
entry get_stor(irg1,irg2);
%fortran
!DEC$ ATTRIBUTES DLLEXPORT :: get_stor
%mortran
irg1=nodes; irg2=nct-1; return;
entry get_err(irg1,irg2); irg1=kstor; irg2=lstor; return;
end;
subroutine ans (x,itre,rtre,cat,yh);
integer itre(6,*); real x(*),rtre(4,*),cat(*);
k=1;
until itre(4,k).lt.0 <
   if itre(1,k).gt.0 <
      if x(itre(1,k)).lt.rtre(1,k) < k=itre(2,k);> else < k=itre(3,k);>
   >
"Naras fix: explicit conversion to integer"
"   else < j=-itre(1,k); kp=rtre(1,k)+0.1; in=0; np=abs(cat(kp))+0.1;"
   else < j=-itre(1,k); kp=int(rtre(1,k)+0.1); in=0; np=int(abs(cat(kp))+0.1);
      <i=1,np; if(x(j).ne.cat(kp+i)) next; in=1; exit;>
      if cat(kp).gt.0.0 < if in.eq.0 < k=itre(2,k);> else < k=itre(3,k);>>
      else < if in.ne.0 < k=itre(2,k);> else < k=itre(3,k);>>
   >
>
yh=rtre(3,k);
return;
end;
subroutine dosort(no,ni,x,w,ms,nu,isc);
integer ms(no,ni,*),isc(no); real x(no,ni),w(no);
data new /1/;
"Naras fix: gfortran says isc is unused; so use it benignly."
isc = isc + 0;
if new.ne.0 <
   <j=1,ni; <i=1,no; ms(i,j,1)=i;> call psort8(x(1,j),ms(1,j,1),1,no);>
   <j=1,ni; <i=1,no; ms(i,j,2)=ms(i,j,1);> nu=0;
      <i=1,no; if(w(ms(i,j,2)).le.0.0) next; nu=nu+1; ms(nu,j,1)=ms(i,j,2);>
   >
   return;
>
<j=1,ni; nu=0;
   <i=1,no; if(w(ms(i,j,2)).le.0.0) next; nu=nu+1; ms(nu,j,1)=ms(i,j,2);>
>
return;
entry set_new(irg); new=irg; return;
end;
subroutine split7 (no,ni,x,y,y2,z,w,lx,m,m1,m2,
   jt,sp,cri1,cri2,w1,w2,crm,kct,cat);
parameter(maxcat=1000, big=9.9e35);
integer lx(ni),m(no,ni);
real x(no,ni),y(no),y2(no),z(no),w(no),cat(*),tcat(maxcat);
data xmiss,ntn,pwr /9.0e35,500,2/;
crm=0.0; jt=0; if(m2-m1+1.lt.2*ntn) return;
<j=1,ni; if(x(m(m1,j),j).ge.x(m(m2,j),j).or.lx(j).eq.0) next;
   if lx(j).eq.1 <
      call eav(x(1,j),y,y2,z,w,m(1,j),m1,m2,ntn,pwr,xmiss,tp,
         cril,criu,wl,wu,cri);
      if(cri.eq.-xmiss) next;
      if abs(cri).ge.crm < crm=abs(cri); cri1=cril; cri2=criu;
         w1=wl; w2=wu; jt=j; sp=tp; crx=cri;
      >
   >
   else <
      call ceav(x(1,j),y,y2,z,w,m(1,j),m1,m2,ntn,pwr,lct,
         tcat,cril,criu,wl,wu,cri); 
      if cri.ge.crm < crm=cri; cri1=cril; cri2=criu; 
         w1=wl; w2=wu; jt=-j; kct=lct; <k=1,kct; cat(k)=tcat(k);>
      >
   >
>
if(jt.gt.0) crm=crx;
return;
entry set_miss(arg);
%fortran
!DEC$ ATTRIBUTES DLLEXPORT :: set_miss
%mortran
xmiss=arg; return;
entry set_ntn(irg);
%fortran
!DEC$ ATTRIBUTES DLLEXPORT :: set_ntn
%mortran
ntn=max0(irg,3); return;
entry set_pwr(arg);
%fortran
!DEC$ ATTRIBUTES DLLEXPORT :: set_pwr
%mortran
pwr=arg; return;
end;
subroutine eav (x,y,y2,z,w,m,m1,m2,ntn,pwr,xmiss,tp,
   cri1s,cri2s,w1s,w2s,cri);
integer m(*); real x(*),y(*),y2(*),z(*),w(*);
data nint,icri /4,1/;
xl=x(m(m1)); if xl.ge.xmiss < cri=-xmiss; return;> nu=m2;
k=m(nu); until x(k).lt.xmiss < nu=nu-1; k=m(nu);>
/crx,crimx/=0.0;
if nu.lt.m2 <   
   call andarm(nu-m1+1,y(m(m1:nu)),y2(m(m1:nu)),
      z(m(m1:nu)),w(m(m1:nu)),cri1xs,w1x);
   call andarm(m2-nu,y(m((nu+1):m2)),y2(m((nu+1):m2)),
      z(m((nu+1):m2)),w(m((nu+1):m2)),cri2xs,w2x);
   crx=max(cri1xs,cri2xs);
   crimx=(crx**pwr)*float(nu-m1+1)*float(m2-nu)/float(m2-m1+1)**2;
>
"Naras fix: explicit conversion to integer"
"rq=(nu-m1+1); mq=rq/nint; kq=1; crim=-xmiss; cris=-xmiss;"
rq=(nu-m1+1); mq=int(rq/nint); kq=1; crim=-xmiss; cris=-xmiss;
<i=nu,m1+1,-1; k=m(i);
   tq=0.5*(x(k)+x(m(i-1)));
   if(tq.le.x(m(i-1))) next; if(tq.ge.x(k)) next;
   if(i-m1.lt.ntn.or.nu-i+1.lt.ntn) next;
   if i.lt.nu-kq*mq+1 < kq=kq+1;
      call andarm(i-m1,y(m(m1:(i-1))),y2(m(m1:(i-1))),
         z(m(m1:(i-1))),w(m(m1:(i-1))),cri1,w1);
      call andarm(nu-i+1,y(m(i:nu)),y2(m(i:nu)),
         z(m(i:nu)),w(m(i:nu)),cri2,w2);
      if icri.eq.1 < cri=max(cri1,cri2);>
      else < cri=abs(cri1-cri2);>
      crit=(cri**pwr)*float(i-m1)*float(nu-i+1)/float(nu-m1+1)**2;
      if(crit.lt.crim) next;
      crim=crit; tp=tq;
      cris=cri; w1s=w1; w2s=w2; cri1s=cri1; cri2s=cri2;
   >
>
cri=cris;
if(x(m(m2)).lt.xmiss) return;
tp=xmiss; cri1s=cri1xs;
cri2s=cri2xs; w1s=w1x; w2s=w2x;
if crx.gt.cri < cri=crx; return;>
cri=-cri;
return;
entry set_qint(irg);
%fortran
!DEC$ ATTRIBUTES DLLEXPORT :: set_qint
%mortran
nint=irg; return;
entry set_cri(irg);
%fortran
!DEC$ ATTRIBUTES DLLEXPORT :: set_cri
%mortran
icri=irg; return;
end;
subroutine ceav (x,y,y2,z,w,m,m1,m2,ntn,pwr,kct,cat,cri1,cri2,w1,w2,cri);
parameter(maxcat=1000);
integer m(*),mt(maxcat);
real x(*),y(*),y2(*),z(*),w(*),cat(*),v(maxcat,3);
fntn=ntn; l=0; i1=m1;
<i=m1,m2-1; k=m(i); if(x(m(i+1)).le.x(k)) next;
   l=l+1; v(l,1)=x(k); i2=i-1;
   call andarm(i2-i1+1,y(m(i1:i2)),y2(m(i1:i2)),
      z(m(i1:i2)),w(m(i1:i2)),v(l,2),v(l,3));
   i1=i;
>
k=m(m2); l=l+1; v(l,1)=x(k);
call andarm(m2-i1+1,y(m(i1:m2)),y2(m(i1:m2)),z(m(i1:m2)),w(m(i1:m2)),
   v(l,2),v(l,3));
<i=1,l; mt(i)=i;> call psort8(v(1:l,2),mt,1,l);
<j=1,l; v(j,2)=v(j,2)*v(j,3);>
/sl,wl,cri,scri/=0.0; sr=sum(v(1:l,2)); wr=sum(v(1:l,3)); kct=0;
<i=1,l-1; k=mt(i); sl=sl+v(k,2); sr=sr-v(k,2);
   wl=wl+v(k,3); wr=wr-v(k,3);
   if(wl.lt.fntn.or.wr.lt.fntn) next;
   c=wr*wl*max(sr/wr,sl/wl)**pwr; if(c.le.cri) next;
   cri=c; kct=i; cri1=sl/wl; cri2=sr/wr; w1=wl; w2=wr; scri=max(cri1,cri2);
>
if kct.eq.0 < cri=0.0; return;>
cat(1)=-kct; <i=1,kct; cat(i+1)=v(mt(i),1);>
cri=scri; kct=kct+1;
return;
end;
subroutine andarm(n,y,y2,z,w,dst,sw);
%fortran
!DEC$ ATTRIBUTES DLLEXPORT :: andarm
%mortran
real y(n),y2(n),z(n),w(n);
call set_kri(kri,2);
if kri.eq.1000 < call rfcall(n, y, z, w, dst); sw=sum(w);>
elseif kri.eq.1 < call andarm1(n,y,z,w,dst,sw);>
elseif kri.eq.2 < call andarm2(n,y,z,w,dst,sw);>
elseif kri.eq.3 < call andarm3(n,y,z,w,dst,sw);>
elseif kri.eq.4 < call andarm4(n,y,z,w,dst,sw);>
elseif kri.eq.5 < call andarm5(n,y,z,w,dst,sw);>
elseif kri.eq.6 < call andarm6(n,y,y2,z,w,dst,sw);>
elseif kri.eq.7 < call andarm7(n,y,z,w,dst,sw);>
elseif kri.eq.8 < call andarm8(n,y,z,w,dst,sw);>
elseif kri.eq.9 < call andarm7(n,y,z,w,dst,sw);>
elseif kri.eq.10 < call andarm10(n,y,z,w,dst,sw);>
"Naras fix: andarm11 never uses args n,y,z,w! Orig followed by changed"
" elseif kri.eq.11 < call andarm11(n,y,z,w,dst,sw);> "
elseif kri.eq.11 < call andarm11(dst,sw);>
elseif kri.eq.12 < call andarm12(n,y,z,w,dst,sw);>
elseif kri.eq.13 < call andarm12(n,y,z,w,dst,sw);>
elseif kri.eq.14 < call andarm14(n,y,z,w,dst,sw);>
else < call andarm15(n,y,y2,z,w,dst,sw);>
return;
end;
subroutine set_kri(irg,jrg);
%fortran
!DEC$ ATTRIBUTES DLLEXPORT :: set_kri
%mortran
save kri;
if jrg.eq.1 < kri=irg; return;>                                           
irg=kri;                                       
return;                                                            
end;
%mortran       
"Naras fix: andarm11 never uses args n,y,z,w! Orig followed by changed"
" subroutine andarm11(n,y,z,w,dst,sw);"
subroutine andarm11(dst,sw);
/dst,sw/=0.0; return; end;
subroutine andarm2(n,y,z,w,dst,sw);
parameter(nmin=100);
real y(n),z(n),w(n); integer my(n),mz(n);
call set_qqtrm(itrm,2);
if n.lt.nmin < dst=0.0; sw=sum(w); return;> 
if n.lt.2*itrm < /dst,sw/=0.0; return;>
<i=1,n; my(i)=i;> call psort8(y,my,1,n);
<i=1,n; mz(i)=i;> call psort8(z,mz,1,n);
/dst,sw1/=0.0;
<i=itrm+1,n-itrm; sw1=sw1+w(my(i));
   dst=dst+w(my(i))*abs(y(my(i))-z(mz(i)));
>
dst=dst/sw1; sw=sum(w);
return;
end;
subroutine set_qqtrm(irg,jrg);
%fortran
!DEC$ ATTRIBUTES DLLEXPORT :: set_qqtrm
%mortran
save itrm;
if jrg.eq.1 <itrm=irg; return;>
irg=itrm;
return;
end;
subroutine andarm1(n,y,z,w,dst,sw);
parameter(eps=1.0e-5,nmin=100);
real y(n),z(n),w(n),q(2*n),w2(2*n); integer m(2*n),iq(2*n);
if n.lt.nmin < dst=0.0; sw=sum(w); return;> 
n2=2*n;
<i=1,n; q(i)=y(i); iq(i)=0; q(i+n)=z(i);
   iq(i+n)=1; w2(i)=w(i); w2(i+n)=w(i);
>
<i=1,n2; m(i)=i;> call psort8(q,m,1,n2);
/sw,tw,dst/=0.0; sw2=2.0*sum(w);
<i=1,n2; k=m(i);
   if iq(k).eq.0 < sw=sw+w2(k);> else < tw=tw+w2(k);>
   pw=(sw+tw)*(sw2-sw-tw);
   pw=max(eps,pw);
   dst=dst+abs(sw-tw)/sqrt(pw);
>
dst=dst/n;
return;
end;
subroutine andarm3(n,y,z,w,dst,sw);
real y(n),z(n),w(n);
sw=sum(w); dst=dot_product(w,abs(y-z))/sw;
return;
end;
subroutine andarm7(n,y,z,w,dst,sw);
parameter(nmin=20);
real y(n),z(n),w(n);
if n.lt.nmin < dst=0.0; sw=sum(w); return;> 
sw=sum(w); dst=abs(dot_product(w,y)/sw-dot_product(w,z)/sw);
return;
end;
subroutine andarm12(n,y,z,w,dst,sw);
parameter(fmin=20);
real y(n),z(n),w(n);
if n.lt.2*int(fmin) < dst=0.0; sw=sum(w); return;>
/dst1,dst2,sw1,sw2/=0.0;
<i=1,n;
   if z(i).lt.0.0 < sw1=sw1+w(i); dst1=dst1+w(i)*y(i);>
   else < sw2=sw2+w(i); dst2=dst2+w(i)*y(i);>
>
sw=sum(w);
if n*sw1/sw.lt.fmin.or.n*sw2/sw.lt.fmin < dst=0.0; return;>
dst=abs(dst2/sw2-dst1/sw1);
return;
end;
subroutine andarm14(n,y,z,w,dst,sw);
parameter(fmin=20,sml=-1.0e20);
real y(n),z(n),w(n);
if n.lt.2*int(fmin) < dst=sml; sw=sum(w); return;>
/dst1,dst2,sw1,sw2/=0.0;
<i=1,n;
   if z(i).lt.0.0 < sw1=sw1+w(i); dst1=dst1+w(i)*y(i);>
   else < sw2=sw2+w(i); dst2=dst2+w(i)*y(i);>
>
sw=sum(w);
if n*sw1/sw.lt.fmin.or.n*sw2/sw.lt.fmin < dst=sml; return;>
dst=dst2/sw2-dst1/sw1;
return;
end;
subroutine andarm8(n,y,z,w,dst,sw);
parameter(nmin=20,sml=-1.0e20);
real y(n),z(n),w(n);
if n.lt.nmin < dst=sml; sw=sum(w); return;> 
sw=sum(w); dst=dot_product(w,y)/sw-dot_product(w,z)/sw;
return;
end;
subroutine andarm4(n,y,z,w,dst,sw);
parameter(maxclass2=10000,nmin=100,idum=2);
real y(n),z(n),w(n),out(maxclass2),dum(2,2);
%fortran
      real, dimension (:,:), allocatable :: costs
%mortran
if n.lt.nmin < dst=0.0; sw=sum(w); return;> 
call classin(2,idum,dum,nclass,out);
%fortran
      allocate(costs(1:nclass,1:nclass),stat=jerr);
%mortran
call reorg(2,nclass,out,costs);
dst=0.0;
"Naras fix: explicit conversion to integer"
"<i=1,n; ky=y(i)+0.1; kz=z(i)+0.1; "
<i=1,n; ky=int(y(i)+0.1); kz=int(z(i)+0.1); 
   dst=dst+w(i)*costs(ky,kz);
>
sw=sum(w); dst=dst/sw;
return;
end;
subroutine andarm5(n,y,z,w,dst,sw);
parameter(nmin=50);
real y(n),z(n),w(n);
data qntl /0.5/;
if n.lt.nmin < dst=0.0; sw=sum(w); return;> 
sw=sum(w); up=0.0;
<i=1,n; if(y(i).le.z(i)) up=up+w(i);>
dst=abs(up/sw-qntl);
return;
entry set_quant(arg);
%fortran
!DEC$ ATTRIBUTES DLLEXPORT :: set_quant
%mortran
qntl=arg; return;
end;
subroutine andarm6(n,y,y2,z,w,dst,sw);
parameter(nit=100,xmiss=9.0e35,thr=1.0e-2,nmin=100);
real y(n),y2(n),z(n),w(n),yy(n,2);
if n.lt.nmin < dst=0.0; sw=sum(w); return;> 
yy(:,1)=y; yy(:,2)=y2;
call cendst(n,yy,z,w,nit,thr,xmiss,dst,sw);
sw=sum(w);
return;
end;
subroutine andarm15(n,y,y2,z,w,dst,sw);
parameter(nit=100,xmiss=9.0e35,thr=1.0e-2,nmin=100);
real y(n),y2(n),z(n),w(n),yy(n,2);
if n.lt.nmin < dst=0.0; sw=sum(w); return;> 
yy(:,1)=y; yy(:,2)=y2;
call cendst1(n,yy,z,w,nit,thr,xmiss,dst,sw);
sw=sum(w);
return;
end;
subroutine andarm10(n,y,z,w,dst,sw);
parameter(eps=1.0e-5,nmin=100);
real y(n),z(n),w(n); integer m(n);
sw=sum(w); if n.lt.nmin < dst=0.0; return;>
/sw1,sw2/=0.0;
<i=1,n; m(i)=i;
   if z(i).lt.0.0 < sw1=sw1+w(i);> else < sw2=sw2+w(i);>
>
call psort8(y,m,1,n);
/s1,s2,s,dst/=0.0;
<i=1,n; k=m(i); s=s+w(k);
   if z(k).lt.0 < s1=s1+w(k)/sw1;>
   else < s2=s2+w(k)/sw2;>
   pw=s*(sw-s); pw=max(eps,pw);
   dst=dst+abs(s1-s2)/sqrt(pw);
>
return;
end;
%fortran
      subroutine stput (iseed)
      real x(*)
      data i /987654321/
      i=iseed
      return
      entry rget (x,n)
      do 1 j=1,n
c Naras fix: explicit conversion of nq to integer
c      i=mod(i*16807.0,2147483647.0)	   
      i=int(mod(i*16807.0,2147483647.0))
      u=i
      u=u*.465661287d-9
c Naras fix: gcc fortran warns about label for statement in do loop
c So put label on a separate continue statement
      x(j) = u
    1 continue
      return
      entry stget (irg)
      irg=i
      return
      end
%mortran
subroutine cendst(n,y,z,w,nit,thr,xmiss,dst,sw);
parameter(eps=0.1);
real y(n,2),z(n),w(n),b(2*n+1),q(3*n),cdf(3*n),r(n);
integer iq(3*n),mm(3*n),mz(n);
data nsamp /500/;
"Naras fix: gfortran warns argument sw is never used, so use it benignly."
sw = sw + 0;
n2=2*n; n3=3*n; sw=sum(w);
<i=1,n; mz(i)=i;> call psort8(z,mz,1,n);
"Naras fix: explicit conversion of nq to integer"
"nq=0.25*n; teps=(z(mz(n-nq))-z(mz(nq)))*eps;"
nq=int(0.25*n); teps=(z(mz(n-nq))-z(mz(nq)))*eps;
<i=1,n; if(y(i,2)-y(i,1).ge.teps) next;
   y(i,1)=y(i,1)-teps; y(i,2)=y(i,2)+teps;
>
<i=1,n; b(i)=y(i,1); b(i+n)=y(i,2);>
m=0;
<i=1,n2; if(b(i).le.-xmiss) next; if(b(i).ge.xmiss) next;
   m=m+1; b(m)=b(i);
>
call unique(m,b,nu);
if nu.gt.nsamp < call rget(r,nsamp);
   <i=1,nsamp; r(i)=b(int(nu*r(i))+1);>
   nu=nsamp; b(1:nu)=r(1:nu);
   call sort(b,nu);
>
m=nu+1; b(m)=xmiss;
call fintcdf1(n,y,m,b,w,nit,thr/m,cdf,jt,err);
m=m-1;
<i=1,m; q(i)=b(i); iq(i)=0;>
<i=1,n; mz(i)=i;> call psort8(z,mz,1,n);
mpn=m+n; k=0;
<i=m+1,mpn; k=k+1; q(i)=z(mz(k)); iq(i)=1;>
<i=1,mpn; mm(i)=i;> call psort8(q,mm,1,mpn);
/ycdf,zcdf,dst,spw/=0.0; /ny,nz/=0;
<i=1,mpn; k=mm(i);
   if iq(k).eq.0 < ny=ny+1; ycdf=cdf(ny);
      pw=float(i)*float(mpn-i)/float(mpn)**2;
      pw=max(eps,pw); pw=1.0/sqrt(pw); spw=spw+pw;
      dst=dst+pw*abs(ycdf-zcdf);
   >
   else < nz=nz+1; zcdf=zcdf+w(nz)/sw;>
>
dst=dst/spw;
return;
entry set_samp(irg);
%fortran
!DEC$ ATTRIBUTES DLLEXPORT :: set_samp
%mortran
nsamp=irg; call set_samp1(nsamp); return;
end;
subroutine cendst1(n,y,z,w,nit,thr,xmiss,dst,sw);
"Naras fix: POSSIBLY WRONG! teps is not defined so add parameter"
parameter(teps=0.01);
real y(n,2),z(n),w(n),b(2*n+1),cdf1(3*n),cdf2(3*n),r(n);
real y1(n,2),y2(n,2),w1(n),w2(n);
data nsamp /500/;
/n1,n2/=0.0;
"Naras fix: gfortran says sw is never used; use it benignly!"
sw = sw + 0;
<i=1,n;
   if(y(i,1).le.-xmiss) next; if(y(i,2).ge.xmiss) next;
   if(y(i,2)-y(i,1).ge.teps) next;
   y(i,1)=y(i,1)-teps; y(i,2)=y(i,2)+teps;
>
<i=1,n;
   if z(i).lt.0.0 < n1=n1+1; y1(n1,:)=y(i,:); w1(n1)=w(i);>
   else < n2=n2+1; y2(n2,:)=y(i,:); w2(n2)=w(i);>
>
<i=1,n; b(i)=y(i,1); b(i+n)=y(i,2);>
m=0;
<i=1,n2; if(b(i).le.-xmiss) next; if(b(i).ge.xmiss) next;
   m=m+1; b(m)=b(i);
>
call unique(m,b,nu);
if nu.gt.nsamp < call rget(r,nsamp);
   <i=1,nsamp; r(i)=b(int(nu*r(i))+1);>
   nu=nsamp; b(1:nu)=r(1:nu);
   call sort(b,nu);
>
m=nu+1; b(m)=xmiss;
call getcdf1(n1,y1,w1,nit,thr,xmiss,nsamp,m,b,cdf1,sw1);
call getcdf1(n2,y2,w2,nit,thr,xmiss,nsamp,m,b,cdf2,sw2);
call diffcdf(m,cdf1,cdf2,dst);
return;
entry set_samp1(irg);
nsamp=irg; return;
end;
subroutine getcdf1(n,y,w,nit,thr,xmiss,nsamp,m,b,cdf,sw);
parameter(teps=0.01);
"Naras fix: gfortran says z,q,r,iq,mm,mz are never used!"
"real y(n,2),z(n),w(n),b(2*n+1),q(3*n),cdf(3*n),r(n);"
real y(n,2),w(n),b(2*n+1),cdf(3*n);
"integer iq(3*n),mm(3*n),mz(n);"
"Naras fix: gfortran says xmiss, nsamp is never used; so use it benignly!"
xmiss = xmiss + 0;
nsamp = nsamp + 0;
n2=2*n; sw=sum(w);
call fintcdf1(n,y,m,b,w,nit,thr/m,cdf,jt,err);
m=m-1;
return;
end;
subroutine diffcdf(m,cdf1,cdf2,dst);
real cdf1(m),cdf2(m);
f12=sqrt(float(m)); dst=0.0;
<i=1,m;  
   dst=dst+abs(cdf1(i)-cdf2(i))/sqrt(float(i)*float(m-i+1));
>
dst=f12*dst/m;
return;
end;
subroutine unique(n,y,nu);
real y(n),yu(n); integer m(n);
<i=1,n; m(i)=i;> call psort8(y,m,1,n);
nu=1; yu(1)=y(m(1));
<i=2,n; if(y(m(i-1)).ge.y(m(i))) next;
   nu=nu+1; yu(nu)=y(m(i));
>
y(1:nu)=yu(1:nu);
return;
end;
subroutine fintcdf1(n,y,m,b,w1,nit,thr,cdf,jt,err);
%fortran
!DEC$ ATTRIBUTES DLLEXPORT :: fintcdf1
%mortran
real y(n,2),b(m),w(n),w1(n),p(m),pij(n,m),ps(m),cdf(m);
real djunk(1); integer ijunk(1);
%fortran
      integer, dimension (:), allocatable :: ic,jc,kc1,kc2,lc
      real, dimension (:), allocatable :: smo
%mortran
call set_vrb(ivrb,2);
w=w1/sum(w1); p=1.0/m;
nt=0;
<i=1,n; <k=1,m; if(y(i,1).ge.b(k)) next;
   if(y(i,2).lt.b(k)) next; nt=nt+1;
>>
allocate(ic(1:(n+1)),stat=jerr);
allocate(jc(1:nt),stat=ierr); jerr=jerr+ierr;
allocate(kc1(1:m),stat=ierr); jerr=jerr+ierr;
allocate(kc2(1:m),stat=ierr); jerr=jerr+ierr;
allocate(smo(1:m),stat=ierr); jerr=jerr+ierr;
allocate(lc(1:nt),stat=ierr); jerr=jerr+ierr;
if jerr.ne.0 < err=8888.0; return;>
nt=0; ic(1)=1;
<i=1,n;
   <k=1,m; if(y(i,1).ge.b(k)) next; if(y(i,2).lt.b(k)) next;
      nt=nt+1; jc(nt)=k;
   >
   ic(i+1)=nt+1;
>
nt=0;
<j=1,m; kc1(j)=nt+1;
   <i=1,n;
      if(y(i,1).ge.b(j)) next; if(y(i,2).lt.b(j)) next;
      nt=nt+1; lc(nt)=i;
   >
   kc2(j)=nt;
>
"Naras Fix: R does not allow FORTRAN output; so use R function"
" if(ivrb.gt.0) write(6,'(''CDF iterations'',$)'); "
if(ivrb.gt.0) < call intpr('CDF iterations', -1, ijunk, 0);>
<it=1,nit; jt=it; ps=p;
   <j=1,m; pij(:,j)=0.0;
      <ii=kc1(j),kc2(j); i=lc(ii);
         s=sum(p(jc(ic(i):(ic(i+1)-1))));
         if s.le.0.0 < err=-7777.0; return;>
         pij(i,j)=w(i)*p(j)/s
      >
      p(j)=sum(pij(:,j));
   >
   if m.gt.100 <
      smo(1)=(2.0*p(1)+p(2))/3.0; smo(m)=(2.0*p(m)+p(m-1))/3.0;   
      smo(2)=0.25*(p(1)+2.0*p(2)+p(3)); smo(m-1)=0.25*(p(m)+2.0*p(m-1)+p(m-2));
      <j=3,m-2; smo(j)=(p(j-2)+2.0*p(j-1)+3.0*p(j)+2.0*p(j+1)+p(j+2))/9.0;>
      p=smo;
   >
   err=sum(abs(p-ps))/m;
   if kbad(err).gt.0 < err=7777.0; return>
   if(err.lt.thr)  exit;
"Naras Fix: R does not allow FORTRAN output; so use R function"
"   if(ivrb.gt.0) write(6,'(''.'',$)');"
   if(ivrb.gt.0) < call intpr('.', -1, ijunk, 0);>
>
cdf(1)=p(1); <j=2,m; cdf(j)=cdf(j-1)+p(j);>
"Naras Fix: R does not allow FORTRAN output; so use R function"
" if ivrb.gt.0 < <w> err; (g10.2);>"
if ivrb.gt.0 < djunk(1) = err; call dblepr('Err = ', -1, djunk, 1);>
return;
end;
subroutine set_vrb(irg,jrg);
%fortran
!DEC$ ATTRIBUTES DLLEXPORT :: set_vrb
%mortran
save ivrb;
if jrg.eq.1 < ivrb=irg; return;>
irg=ivrb;
return;
end;
subroutine cdfpoints1(m,x,n,y,w,cdf);
%fortran
!DEC$ ATTRIBUTES DLLEXPORT :: cdfpoints1
%mortran
real x(m),y(n),w(n),cdf(m);
i=1; j=0; sw=0.0;
loop < j=j+1; if(j.gt.m) go to :out2:;
   until y(i).gt.x(j) < sw=sw+w(i); i=i+1;> until i.gt.n;
   if i.gt.n <
      <k=j,m; cdf(k)=sw;> go to :out1:;
   >
   cdf(j)=sw;
>
:out1: continue;
:out2: continue;
cdf=cdf/sum(w);
return;
end;
subroutine sort(x,n);
real x(n),z(n); integer m(n);
<i=1,n; m(i)=i> z=x;
call psort8(z,m,1,n);
<i=1,n; x(i)=z(m(i));>
return;
end;
function kbad(u);
kbad=0;
if(isnan(u).or.abs(u).ge.abs(huge(u))) kbad=1;
return;
end;
subroutine classin(ient,nclasssv,costssv,nout,out);
%fortran
!DEC$ ATTRIBUTES DLLEXPORT :: classin
%mortran
real costssv(nclasssv,nclasssv),out(1);
%fortran
      real, dimension (:), allocatable :: costs
%mortran
save costs,nclass;
nq=nclasssv*nclasssv;
allocate(costs(1:nq),stat=jerr);
if ient.eq.1 < nclass=nclasssv; 
   call reorg(1,nclass,costs,costssv);
   nout=1; out=1.0;
>
else < nout=nclass; call reorg(2,nclass,costs,out);>
return;
end;
subroutine reorg(ient,n,a,b);
real a(n*n),b(n,n);
i=0;
if ient.eq.2 <<k=1,n; <j=1,n; i=i+1; b(j,k)=a(i);>>>
else <<k=1,n; <j=1,n; i=i+1; a(i)=b(j,k);>>>
return;
end; 
subroutine crinode (itr,rtr,mxnodes,node,nodes,cri,wt);
%fortran
!DEC$ ATTRIBUTES DLLEXPORT :: crinode
%mortran
integer itr(6,*),nodes(mxnodes),m(mxnodes),ic(mxnodes);
real rtr(4,*),cri(mxnodes),wt(mxnodes),sc(mxnodes,2);
node=0; k=itr(2,1);
loop < if itr(4,k).ge.0 < k=itr(2,k); next;> i1=itr(5,k); i2=itr(6,k);
   node=node+1; if(node.gt.mxnodes) return;
   nodes(node)=k; cri(node)=rtr(3,k); wt(node)=rtr(4,k);
   until k.eq.itr(2,iabs(itr(4,k))) < k=iabs(itr(4,k));> until k.eq.1;
   if(k.eq.1) exit; k=itr(3,iabs(itr(4,k)));
>
<k=1,node; m(k)=k;> call psort8(-cri,m,1,node);
<i=1,node; ic(i)=nodes(m(i)); sc(i,1)=cri(m(i)); sc(i,2)=wt(m(i));>
<i=1,node; nodes(i)=ic(i); cri(i)=sc(i,1); wt(i)=sc(i,2);> 
return;
end;
subroutine prune1 (itr,rtr,nodes,thr,itro,rtro);
%fortran
!DEC$ ATTRIBUTES DLLEXPORT :: prune1
%mortran
integer itr(6,nodes); real rtr(4,nodes);
integer itro(6,nodes); real rtro(4,nodes);
call prune(itr,rtr,nodes,thr);
itro=itr; rtro=rtr;
return;
end;
subroutine prune (itr,rtr,nodes,thr);
integer itr(6,nodes); real rtr(4,nodes);
loop < nch=0;
   <k=1,nodes; if(itr(4,k).le.0) next;
      nl=itr(2,k); nr=itr(3,k);
      if (itr(4,nl).ge.0) next;
      if (itr(4,nr).ge.0) next;
      if(max(rtr(3,nl),rtr(3,nr)).gt.rtr(3,k)+thr) next;
      itr(4,k)=-itr(4,k); nch=nch+1;
   >
> until nch.eq.0;
return;
end;
subroutine getnode (x,itre,rtre,cat,node);
integer itre(6,*); real x(*),rtre(4,*),cat(*);
k=1;
until itre(4,k).lt.0 <
   if itre(1,k).gt.0 <
      if x(itre(1,k)).lt.rtre(1,k) < k=itre(2,k);> else < k=itre(3,k);>
   >
"Naras fix: explicit conversion to integer"
"   else < j=-itre(1,k); kp=rtre(1,k)+0.1; in=0; np=abs(cat(kp))+0.1;"
   else < j=-itre(1,k); kp=int(rtre(1,k)+0.1); in=0; np=int(abs(cat(kp))+0.1);
      <i=1,np; if(x(j).ne.cat(kp+i)) next; in=1; exit;>
      if cat(kp).gt.0.0 < if in.eq.0 < k=itre(2,k);> else < k=itre(3,k);>>
      else < if in.ne.0 < k=itre(2,k);> else < k=itre(3,k);>>
   > 
>
node=k;
return;
end;
subroutine getnodes1 (no,ni,x,itre,rtre,cat,nodes);
%fortran
!DEC$ ATTRIBUTES DLLEXPORT :: getnodes1
%mortran
integer itre(6,*),nodes(no); real x(no,ni),rtre(4,*),cat(*);
<i=1,no; call getnode (x(i,:),itre,rtre,cat,node); nodes(i)=node;>
return;
end;
subroutine getlims(node,ni,itr,rtr,cat,nvar,jvar,vlims,jerr);
%fortran
!DEC$ ATTRIBUTES DLLEXPORT :: getlims
%mortran
integer itr(6,*),jvar(2,*); real rtr(4,*),cat(*),vlims(*);
"Naras fix: gfortran says ni is unused; so use it benignly!"
ni = ni + 0;
jerr=0; if itr(4,node).ge.0 < jerr=1; return;>
nvar=0;
k=node;
loop < nvar=nvar+1; kpp=abs(itr(4,k));
   if itr(1,kpp).gt.0 < jvar(2,nvar)=0;
      if itr(2,kpp).eq.k < jvar(1,nvar)=-itr(1,kpp);>
      else < jvar(1,nvar)=itr(1,kpp);>
      vlims(nvar)=rtr(1,kpp);
   >
   else <
      if k.eq.itr(2,kpp) < sgn=-1.0;> else < sgn=1.0;>
"Naras fix: explicit conversion to integer"
"      jvar(1,nvar)=-itr(1,kpp); kp=rtr(1,kpp)+0.1; jvar(2,nvar)=kp;"
jvar(1,nvar)=-itr(1,kpp); kp=int(rtr(1,kpp)+0.1); jvar(2,nvar)=kp;
      vlims(nvar)=sgn*abs(cat(kp));
   > 
   k=kpp;
> until kpp.eq.1;
return;
end;
subroutine trans(ny,y,wy,nz,z,wz,nt,t);
%fortran
!DEC$ ATTRIBUTES DLLEXPORT :: trans
%mortran
real y(ny),wy(ny),z(nz),wz(nz),t(nt+2,2),u(max(ny,nz));
real p(nt);
integer m(max(ny,nz));
<i=1,ny; m(i)=i; u(i)=y(i);> call psort8(u,m,1,ny);
<i=1,ny; y(i)=u(m(i));> u=wy; <i=1,ny; wy(i)=u(m(i));>
<i=1,nz; m(i)=i; u(i)=z(i);> call psort8(u,m,1,nz);
<i=1,nz; z(i)=u(m(i));>  u=wz; <i=1,nz; wz(i)=u(m(i));>
<i=1,nt; p(i)=(i-0.5)/float(nt);>
call untie(ny,y,u);
call qntl(ny,u,wy,nt,p,t(:,1));
call untie(nz,z,u);
call qntl(nz,u,wz,nt,p,t(:,2));
return;
end;
subroutine qntl(n,y,w,nq,p,q);
real y(n),w(n),p(nq),q(nq+2);
sw=sum(w); k=1; ff=w(1); q(1)=y(1); q(nq+2)=y(n);
<i=2,n; ff=ff+w(i); pp=ff/sw; if(pp.lt.p(k)) next;
   k=k+1; q(k)=0.5*(y(i)+y(i-1));
   if(k.ge.nq) exit;
>
q(nq+1)=0.5*(q(nq+2)+q(nq));
return;
end;
subroutine untie(n,y,u);
%fortran
!DEC$ ATTRIBUTES DLLEXPORT :: untie
%mortran
real y(n),u(n);
i=1; k=0;
until i.ge.n < 
   if y(i+1).gt.y(i) < k=k+1; u(k)=y(i); i=i+1; next;>
   i1=i; until y(i+1).gt.y(i) < i=i+1;> until i.ge.n; i2=i;
   if i1.le.1 < a=y(i1+1); b=y(i2+1); u(1)=y(1); k=1;  
      <j=i1+1,i2; k=k+1; u(k)=a+(b-a)*(j-i1)/(i2-i1+1);>
      i=i2+1;
   >
   elseif i2.ge.n <
      a=y(i1-1); b=(y(i2)-a)/(i2-i1+1);
      <j=i1,i2; k=k+1; u(k)=a+b*(j-i1+1);>
   >
   else < a=y(i1-1); b=y(i2);
      <j=i1,i2; k=k+1; u(k)=a+(b-a)*(j-i1+1)/(i2-i1+1);>
      i=i+1;
   >
>
if k.lt.n < k=k+1; u(k)=y(n);>
return;
end;
%fortran
      subroutine psort8 (v,a,ii,jj)
c
c     puts into a the permutation vector which sorts v into
c     increasing order. the array v is not modified.
c     only elements from ii to jj are considered.
c     arrays iu(k) and il(k) permit sorting up to 2**(k+1)-1 elements
c
c     this is a modification of cacm algorithm #347 by r. c. singleton,
c     which is a modified hoare quicksort.
c
      dimension a(jj),v(jj),iu(20),il(20)
      integer t,tt
      integer a
      real v
      m=1
      i=ii
      j=jj
 10   if (i.ge.j) go to 80
 20   k=i
      ij=(j+i)/2
      t=a(ij)
      vt=v(t)
      if (v(a(i)).le.vt) go to 30
      a(ij)=a(i)
      a(i)=t
      t=a(ij)
      vt=v(t)
 30   l=j
      if (v(a(j)).ge.vt) go to 50
      a(ij)=a(j)
      a(j)=t
      t=a(ij)
      vt=v(t)
      if (v(a(i)).le.vt) go to 50
      a(ij)=a(i)
      a(i)=t
      t=a(ij)
      vt=v(t)
      go to 50
 40   a(l)=a(k)
      a(k)=tt
 50   l=l-1
      if (v(a(l)).gt.vt) go to 50
      tt=a(l)
      vtt=v(tt)
 60   k=k+1
      if (v(a(k)).lt.vt) go to 60
      if (k.le.l) go to 40
      if (l-i.le.j-k) go to 70
      il(m)=i
      iu(m)=l
      i=k
      m=m+1
      go to 90
 70   il(m)=k
      iu(m)=j
      j=l
      m=m+1
      go to 90
 80   m=m-1
      if (m.eq.0) return
      i=il(m)
      j=iu(m)
 90   if (j-i.gt.10) go to 20
      if (i.eq.ii) go to 10
      i=i-1
 100  i=i+1
      if (i.eq.j) go to 80
      t=a(i+1)
      vt=v(t)
      if (v(a(i)).le.vt) go to 100
      k=i
 110  a(k+1)=a(k)
      k=k-1
      if (vt.lt.v(a(k))) go to 110
      a(k+1)=t
      go to 100
      end      
%%
