c     mortran 2.0     (version of 7/04/75 mod 7/4/87 (ajc))             
      subroutine fcontrast (no,ni,x,y,y2,z,w,lx,mxt,itre,rtre,mxc,cat,ms
     *,isc)
      implicit double precision(a-h,o-z)                                
!DEC$ ATTRIBUTES DLLEXPORT :: fcontrast                                 
      integer lx(ni),ms(no,ni,2),itre(6,mxt),isc(no)                    
      double precision x(no,ni),y(no),y2(no),z(no),w(no),rtre(4,mxt),cat
     *(mxc)
      save nodes,nct,kstor,lstor                                        
      data mxtrm /10/                                                   
      call dosort(no,ni,x,w,ms,nu,isc)                                  
      lstor=0                                                           
      itre(4,1)=0                                                       
      itre(5,1)=1                                                       
      itre(6,1)=nu                                                      
      itre(2,1)=2                                                       
      itre(3,1)=3                                                       
      call andarm(no,y,y2,z,w,rtre(2,1),rtre(4,1))                      
      nct=1                                                             
      call split7(no,ni,x,y,y2,z,w,lx,ms,1,nu,itre(1,1),rtre(1,1),  cri1
     *,cri2,w1,w2,rtre(3,1),kct,cat(nct))
      if(itre(1,1) .ne. 0)goto 10021                                    
      itre(4,1)=-9999                                                   
      return                                                            
10021 continue                                                          
      if(itre(1,1) .ge. 0)goto 10041                                    
      rtre(1,1)=nct                                                     
      nct=nct+kct                                                       
      if(nct.gt.mxc) lstor=1                                            
10041 continue                                                          
      rtre(2,2)=0.0                                                     
      rtre(2,3)=rtre(2,2)                                               
      rtre(3,2)=cri1                                                    
      rtre(3,3)=cri2                                                    
      itre(4,2)=-1                                                      
      itre(4,3)=itre(4,2)                                               
      rtre(4,2)=w1                                                      
      rtre(4,3)=w2                                                      
      rtre(2,1)=sign(max(0.0,rtre(3,1)-rtre(2,1)),rtre(3,1))            
      if(lstor.ne.0) return                                             
      nodes=3                                                           
      ntrm=0                                                            
      if(rtre(2,1).gt.0.0) ntrm=1                                       
      ktrg=1                                                            
      kstor=0                                                           
10050 continue                                                          
10051 if(ntrm.ge.mxtrm)goto 10052                                       
      jt=itre(1,ktrg)                                                   
      st=rtre(1,ktrg)                                                   
      k5=itre(5,ktrg)                                                   
      k6=itre(6,ktrg)                                                   
      if(jt .ge. 0)goto 10071                                           
      ju=-jt                                                            
      kp=st+0.1                                                         
      np=abs(cat(kp))+0.1                                               
10071 continue                                                          
10080 do 10081 j=1,ni                                                   
      if(ni.gt.1.and.j.eq.jt)goto 10081                                 
      kl=k5-1                                                           
      kr=k6+1                                                           
10090 do 10091 i=k5,k6                                                  
      l=ms(i,j,1)                                                       
      if(jt .le. 0)goto 10111                                           
      if(x(l,jt) .ge. st)goto 10131                                     
      kl=kl+1                                                           
      isc(kl)=l                                                         
      goto 10141                                                        
10131 continue                                                          
      kr=kr-1                                                           
      isc(kr)=l                                                         
10141 continue                                                          
10121 continue                                                          
      goto 10151                                                        
10111 continue                                                          
      in=0                                                              
10160 do 10161 im=1,np                                                  
      if(x(l,ju).ne.cat(kp+im))goto 10161                               
      in=1                                                              
      goto 10162                                                        
10161 continue                                                          
10162 continue                                                          
      if(cat(kp) .le. 0.0)goto 10181                                    
      if(in .ne. 0)goto 10201                                           
      kl=kl+1                                                           
      isc(kl)=l                                                         
      goto 10211                                                        
10201 continue                                                          
      kr=kr-1                                                           
      isc(kr)=l                                                         
10211 continue                                                          
10191 continue                                                          
      goto 10221                                                        
10181 continue                                                          
      if(in .eq. 0)goto 10241                                           
      kl=kl+1                                                           
      isc(kl)=l                                                         
      goto 10251                                                        
10241 continue                                                          
      kr=kr-1                                                           
      isc(kr)=l                                                         
10251 continue                                                          
10231 continue                                                          
10221 continue                                                          
10171 continue                                                          
10151 continue                                                          
10101 continue                                                          
10091 continue                                                          
10092 continue                                                          
10260 do 10261 i=k5,kl                                                  
      ms(i,j,1)=isc(i)                                                  
10261 continue                                                          
10262 continue                                                          
10270 do 10271 i=kr,k6                                                  
      ms(i,j,1)=isc(k6-i+kr)                                            
10271 continue                                                          
10272 continue                                                          
10081 continue                                                          
10082 continue                                                          
      itre(5,itre(2,ktrg))=k5                                           
      itre(6,itre(2,ktrg))=kl                                           
      itre(5,itre(3,ktrg))=kr                                           
      itre(6,itre(3,ktrg))=k6                                           
      nde=itre(2,ktrg)                                                  
      call split7(no,ni,x,y,y2,z,w,lx,ms,itre(5,nde),itre(6,nde),itre(1,
     *nde),  rtre(1,nde),cri1,cri2,w1,w2,cri,kct,cat(nct))
      if(itre(1,nde) .ge. 0)goto 10291                                  
      rtre(1,nde)=nct                                                   
      nct=nct+kct                                                       
      if(nct.gt.mxc) lstor=1                                            
10291 continue                                                          
      if(nodes+2 .le. mxt)goto 10311                                    
      kstor=1                                                           
      goto 10321                                                        
10311 continue                                                          
      itre(2,nde)=nodes+1                                               
      itre(3,nde)=nodes+2                                               
      l=nodes+1                                                         
      rtre(2,l)=0.0                                                     
      rtre(2,l+1)=rtre(2,l)                                             
      rtre(3,l)=cri1                                                    
      rtre(3,l+1)=cri2                                                  
      rtre(4,l)=w1                                                      
      rtre(4,l+1)=w2                                                    
      itre(4,l)=-nde                                                    
      itre(4,l+1)=itre(4,l)                                             
      rtre(2,nde)=sign(max(0.0,cri-rtre(3,nde)),cri)                    
10321 continue                                                          
10301 continue                                                          
      nde=itre(3,ktrg)                                                  
      call split7(no,ni,x,y,y2,z,w,lx,ms,itre(5,nde),itre(6,nde),itre(1,
     *nde),  rtre(1,nde),cri1,cri2,w1,w2,cri,kct,cat(nct))
      if(itre(1,nde) .ge. 0)goto 10341                                  
      rtre(1,nde)=nct                                                   
      nct=nct+kct                                                       
      if(nct.gt.mxc) lstor=1                                            
10341 continue                                                          
      if(nodes+4 .le. mxt)goto 10361                                    
      kstor=1                                                           
      goto 10371                                                        
10361 continue                                                          
      itre(2,nde)=nodes+3                                               
      itre(3,nde)=nodes+4                                               
      l=nodes+3                                                         
      rtre(2,l)=0.0                                                     
      rtre(2,l+1)=rtre(2,l)                                             
      rtre(3,l)=cri1                                                    
      rtre(3,l+1)=cri2                                                  
      rtre(4,l)=w1                                                      
      rtre(4,l+1)=w2                                                    
      itre(4,l)=-nde                                                    
      itre(4,l+1)=itre(4,l)                                             
      rtre(2,nde)=sign(max(0.0,cri-rtre(3,nde)),cri)                    
      itre(4,ktrg)=-itre(4,ktrg)                                        
      rsv=rtre(2,ktrg)                                                  
      crix=0.0                                                          
10380 do 10381 k=1,nodes                                                
      if(itre(4,k).ge.0)goto 10381                                      
      if(abs(rtre(2,k)).le.crix)goto 10381                              
      if(itre(1,k).eq.0)goto 10381                                      
      crix=abs(rtre(2,k))                                               
      ktrg=k                                                            
10381 continue                                                          
10382 continue                                                          
10371 continue                                                          
10351 continue                                                          
      if(crix.le.0.0.or.kstor.ne.0.or.lstor.ne.0) return                
      nodes=nodes+4                                                     
      if(rsv.gt.0.0) ntrm=ntrm+1                                        
      goto 10051                                                        
10052 continue                                                          
      return                                                            
      entry set_trm(irg)                                                
!DEC$ ATTRIBUTES DLLEXPORT :: set_trm                                   
      mxtrm=max0(irg,2)                                                 
      return                                                            
      entry get_stor(irg1,irg2)                                         
!DEC$ ATTRIBUTES DLLEXPORT :: get_stor                                  
      irg1=nodes                                                        
      irg2=nct-1                                                        
      return                                                            
      entry get_err(irg1,irg2)                                          
      irg1=kstor                                                        
      irg2=lstor                                                        
      return                                                            
      end                                                               
      subroutine ans (x,itre,rtre,cat,yh)                               
      implicit double precision(a-h,o-z)                                
      integer itre(6,*)                                                 
      double precision x(*),rtre(4,*),cat(*)                            
      k=1                                                               
10390 continue                                                          
10391 if(itre(4,k).lt.0)goto 10392                                      
      if(itre(1,k) .le. 0)goto 10411                                    
      if(x(itre(1,k)) .ge. rtre(1,k))goto 10431                         
      k=itre(2,k)                                                       
      goto 10441                                                        
10431 continue                                                          
      k=itre(3,k)                                                       
10441 continue                                                          
10421 continue                                                          
      goto 10451                                                        
10411 continue                                                          
      j=-itre(1,k)                                                      
      kp=rtre(1,k)+0.1                                                  
      in=0                                                              
      np=abs(cat(kp))+0.1                                               
10460 do 10461 i=1,np                                                   
      if(x(j).ne.cat(kp+i))goto 10461                                   
      in=1                                                              
      goto 10462                                                        
10461 continue                                                          
10462 continue                                                          
      if(cat(kp) .le. 0.0)goto 10481                                    
      if(in .ne. 0)goto 10501                                           
      k=itre(2,k)                                                       
      goto 10511                                                        
10501 continue                                                          
      k=itre(3,k)                                                       
10511 continue                                                          
10491 continue                                                          
      goto 10521                                                        
10481 continue                                                          
      if(in .eq. 0)goto 10541                                           
      k=itre(2,k)                                                       
      goto 10551                                                        
10541 continue                                                          
      k=itre(3,k)                                                       
10551 continue                                                          
10531 continue                                                          
10521 continue                                                          
10471 continue                                                          
10451 continue                                                          
10401 continue                                                          
      goto 10391                                                        
10392 continue                                                          
      yh=rtre(3,k)                                                      
      return                                                            
      end                                                               
      subroutine dosort(no,ni,x,w,ms,nu,isc)                            
      implicit double precision(a-h,o-z)                                
      integer ms(no,ni,*),isc(no)                                       
      double precision x(no,ni),w(no)                                   
      data new /1/                                                      
      if(new .eq. 0)goto 10571                                          
10580 do 10581 j=1,ni                                                   
10590 do 10591 i=1,no                                                   
      ms(i,j,1)=i                                                       
10591 continue                                                          
10592 continue                                                          
      call psort8(x(1,j),ms(1,j,1),1,no)                                
10581 continue                                                          
10582 continue                                                          
10600 do 10601 j=1,ni                                                   
10610 do 10611 i=1,no                                                   
      ms(i,j,2)=ms(i,j,1)                                               
10611 continue                                                          
10612 continue                                                          
      nu=0                                                              
10620 do 10621 i=1,no                                                   
      if(w(ms(i,j,2)).le.0.0)goto 10621                                 
      nu=nu+1                                                           
      ms(nu,j,1)=ms(i,j,2)                                              
10621 continue                                                          
10622 continue                                                          
10601 continue                                                          
10602 continue                                                          
      return                                                            
10571 continue                                                          
10630 do 10631 j=1,ni                                                   
      nu=0                                                              
10640 do 10641 i=1,no                                                   
      if(w(ms(i,j,2)).le.0.0)goto 10641                                 
      nu=nu+1                                                           
      ms(nu,j,1)=ms(i,j,2)                                              
10641 continue                                                          
10642 continue                                                          
10631 continue                                                          
10632 continue                                                          
      return                                                            
      entry set_new(irg)                                                
      new=irg                                                           
      return                                                            
      end                                                               
      subroutine split7 (no,ni,x,y,y2,z,w,lx,m,m1,m2,  jt,sp,cri1,cri2,w
     *1,w2,crm,kct,cat)
      implicit double precision(a-h,o-z)                                
      parameter(maxcat=1000, big=9.9e35)                                
      integer lx(ni),m(no,ni)                                           
      double precision x(no,ni),y(no),y2(no),z(no),w(no),cat(*),tcat(max
     *cat)
      data xmiss,ntn,pwr /9.0e35,500,2/                                 
      crm=0.0                                                           
      jt=0                                                              
      if(m2-m1+1.lt.2*ntn) return                                       
10650 do 10651 j=1,ni                                                   
      if(x(m(m1,j),j).ge.x(m(m2,j),j).or.lx(j).eq.0)goto 10651          
      if(lx(j) .ne. 1)goto 10671                                        
      call eav(x(1,j),y,y2,z,w,m(1,j),m1,m2,ntn,pwr,xmiss,tp,  cril,criu
     *,wl,wu,cri)
      if(cri.eq.-xmiss)goto 10651                                       
      if(abs(cri) .lt. crm)goto 10691                                   
      crm=abs(cri)                                                      
      cri1=cril                                                         
      cri2=criu                                                         
      w1=wl                                                             
      w2=wu                                                             
      jt=j                                                              
      sp=tp                                                             
      crx=cri                                                           
10691 continue                                                          
      goto 10701                                                        
10671 continue                                                          
      call ceav(x(1,j),y,y2,z,w,m(1,j),m1,m2,ntn,pwr,lct,  tcat,cril,cri
     *u,wl,wu,cri)
      if(cri .lt. crm)goto 10721                                        
      crm=cri                                                           
      cri1=cril                                                         
      cri2=criu                                                         
      w1=wl                                                             
      w2=wu                                                             
      jt=-j                                                             
      kct=lct                                                           
10730 do 10731 k=1,kct                                                  
      cat(k)=tcat(k)                                                    
10731 continue                                                          
10732 continue                                                          
10721 continue                                                          
10701 continue                                                          
10661 continue                                                          
10651 continue                                                          
10652 continue                                                          
      if(jt.gt.0) crm=crx                                               
      return                                                            
      entry set_miss(arg)                                               
!DEC$ ATTRIBUTES DLLEXPORT :: set_miss                                  
      xmiss=arg                                                         
      return                                                            
      entry set_ntn(irg)                                                
!DEC$ ATTRIBUTES DLLEXPORT :: set_ntn                                   
      ntn=max0(irg,3)                                                   
      return                                                            
      entry set_pwr(arg)                                                
!DEC$ ATTRIBUTES DLLEXPORT :: set_pwr                                   
      pwr=arg                                                           
      return                                                            
      end                                                               
      subroutine eav (x,y,y2,z,w,m,m1,m2,ntn,pwr,xmiss,tp,  cri1s,cri2s,
     *w1s,w2s,cri)
      implicit double precision(a-h,o-z)                                
      integer m(*)                                                      
      double precision x(*),y(*),y2(*),z(*),w(*)                        
      data nint,icri /4,1/                                              
      xl=x(m(m1))                                                       
      if(xl .lt. xmiss)goto 10751                                       
      cri=-xmiss                                                        
      return                                                            
10751 continue                                                          
      nu=m2                                                             
      k=m(nu)                                                           
10760 continue                                                          
10761 if(x(k).lt.xmiss)goto 10762                                       
      nu=nu-1                                                           
      k=m(nu)                                                           
      goto 10761                                                        
10762 continue                                                          
      crx=0.0                                                           
      crimx=crx                                                         
      if(nu .ge. m2)goto 10781                                          
      call andarm(nu-m1+1,y(m(m1:nu)),y2(m(m1:nu)),  z(m(m1:nu)),w(m(m1:
     *nu)),cri1xs,w1x)
      call andarm(m2-nu,y(m((nu+1):m2)),y2(m((nu+1):m2)),  z(m((nu+1):m2
     *)),w(m((nu+1):m2)),cri2xs,w2x)
      crx=max(cri1xs,cri2xs)                                            
      crimx=(crx**pwr)*float(nu-m1+1)*float(m2-nu)/float(m2-m1+1)**2    
10781 continue                                                          
      rq=(nu-m1+1)                                                      
      mq=rq/nint                                                        
      kq=1                                                              
      crim=-xmiss                                                       
      cris=-xmiss                                                       
10790 do 10791 i=nu,m1+1,-1                                             
      k=m(i)                                                            
      tq=0.5*(x(k)+x(m(i-1)))                                           
      if(tq.le.x(m(i-1)))goto 10791                                     
      if(tq.ge.x(k))goto 10791                                          
      if(i-m1.lt.ntn.or.nu-i+1.lt.ntn)goto 10791                        
      if(i .ge. nu-kq*mq+1)goto 10811                                   
      kq=kq+1                                                           
      call andarm(i-m1,y(m(m1:(i-1))),y2(m(m1:(i-1))),  z(m(m1:(i-1))),w
     *(m(m1:(i-1))),cri1,w1)
      call andarm(nu-i+1,y(m(i:nu)),y2(m(i:nu)),  z(m(i:nu)),w(m(i:nu)),
     *cri2,w2)
      if(icri .ne. 1)goto 10831                                         
      cri=max(cri1,cri2)                                                
      goto 10841                                                        
10831 continue                                                          
      cri=abs(cri1-cri2)                                                
10841 continue                                                          
10821 continue                                                          
      crit=(cri**pwr)*float(i-m1)*float(nu-i+1)/float(nu-m1+1)**2       
      if(crit.lt.crim)goto 10791                                        
      crim=crit                                                         
      tp=tq                                                             
      cris=cri                                                          
      w1s=w1                                                            
      w2s=w2                                                            
      cri1s=cri1                                                        
      cri2s=cri2                                                        
10811 continue                                                          
10791 continue                                                          
10792 continue                                                          
      cri=cris                                                          
      if(x(m(m2)).lt.xmiss) return                                      
      tp=xmiss                                                          
      cri1s=cri1xs                                                      
      cri2s=cri2xs                                                      
      w1s=w1x                                                           
      w2s=w2x                                                           
      if(crx .le. cri)goto 10861                                        
      cri=crx                                                           
      return                                                            
10861 continue                                                          
      cri=-cri                                                          
      return                                                            
      entry set_qint(irg)                                               
!DEC$ ATTRIBUTES DLLEXPORT :: set_qint                                  
      nint=irg                                                          
      return                                                            
      entry set_cri(irg)                                                
!DEC$ ATTRIBUTES DLLEXPORT :: set_cri                                   
      icri=irg                                                          
      return                                                            
      end                                                               
      subroutine ceav (x,y,y2,z,w,m,m1,m2,ntn,pwr,kct,cat,cri1,cri2,w1,w
     *2,cri)
      implicit double precision(a-h,o-z)                                
      parameter(maxcat=1000)                                            
      integer m(*),mt(maxcat)                                           
      double precision x(*),y(*),y2(*),z(*),w(*),cat(*),v(maxcat,3)     
      fntn=ntn                                                          
      l=0                                                               
      i1=m1                                                             
10870 do 10871 i=m1,m2-1                                                
      k=m(i)                                                            
      if(x(m(i+1)).le.x(k))goto 10871                                   
      l=l+1                                                             
      v(l,1)=x(k)                                                       
      i2=i-1                                                            
      call andarm(i2-i1+1,y(m(i1:i2)),y2(m(i1:i2)),  z(m(i1:i2)),w(m(i1:
     *i2)),v(l,2),v(l,3))
      i1=i                                                              
10871 continue                                                          
10872 continue                                                          
      k=m(m2)                                                           
      l=l+1                                                             
      v(l,1)=x(k)                                                       
      call andarm(m2-i1+1,y(m(i1:m2)),y2(m(i1:m2)),z(m(i1:m2)),w(m(i1:m2
     *)),  v(l,2),v(l,3))
10880 do 10881 i=1,l                                                    
      mt(i)=i                                                           
10881 continue                                                          
10882 continue                                                          
      call psort8(v(1:l,2),mt,1,l)                                      
10890 do 10891 j=1,l                                                    
      v(j,2)=v(j,2)*v(j,3)                                              
10891 continue                                                          
10892 continue                                                          
      sl=0.0                                                            
      wl=sl                                                             
      cri=wl                                                            
      scri=cri                                                          
      sr=sum(v(1:l,2))                                                  
      wr=sum(v(1:l,3))                                                  
      kct=0                                                             
10900 do 10901 i=1,l-1                                                  
      k=mt(i)                                                           
      sl=sl+v(k,2)                                                      
      sr=sr-v(k,2)                                                      
      wl=wl+v(k,3)                                                      
      wr=wr-v(k,3)                                                      
      if(wl.lt.fntn.or.wr.lt.fntn)goto 10901                            
      c=wr*wl*max(sr/wr,sl/wl)**pwr                                     
      if(c.le.cri)goto 10901                                            
      cri=c                                                             
      kct=i                                                             
      cri1=sl/wl                                                        
      cri2=sr/wr                                                        
      w1=wl                                                             
      w2=wr                                                             
      scri=max(cri1,cri2)                                               
10901 continue                                                          
10902 continue                                                          
      if(kct .ne. 0)goto 10921                                          
      cri=0.0                                                           
      return                                                            
10921 continue                                                          
      cat(1)=-kct                                                       
10930 do 10931 i=1,kct                                                  
      cat(i+1)=v(mt(i),1)                                               
10931 continue                                                          
10932 continue                                                          
      cri=scri                                                          
      kct=kct+1                                                         
      return                                                            
      end                                                               
      subroutine andarm(n,y,y2,z,w,dst,sw)                              
      implicit double precision(a-h,o-z)                                
!DEC$ ATTRIBUTES DLLEXPORT :: andarm                                    
      double precision y(n),y2(n),z(n),w(n)                             
      call set_kri(kri,2)                                               
      if(kri .ne. 1)goto 10951                                          
      call andarm1(n,y,z,w,dst,sw)                                      
      goto 10941                                                        
10951 if(kri .ne. 2)goto 10961                                          
      call andarm2(n,y,z,w,dst,sw)                                      
      goto 10941                                                        
10961 if(kri .ne. 3)goto 10971                                          
      call andarm3(n,y,z,w,dst,sw)                                      
      goto 10941                                                        
10971 if(kri .ne. 4)goto 10981                                          
      call andarm4(n,y,z,w,dst,sw)                                      
      goto 10941                                                        
10981 if(kri .ne. 5)goto 10991                                          
      call andarm5(n,y,z,w,dst,sw)                                      
      goto 10941                                                        
10991 if(kri .ne. 6)goto 11001                                          
      call andarm6(n,y,y2,z,w,dst,sw)                                   
      goto 10941                                                        
11001 if(kri .ne. 7)goto 11011                                          
      call andarm7(n,y,z,w,dst,sw)                                      
      goto 10941                                                        
11011 if(kri .ne. 8)goto 11021                                          
      call andarm8(n,y,z,w,dst,sw)                                      
      goto 10941                                                        
11021 if(kri .ne. 9)goto 11031                                          
      call andarm7(n,y,z,w,dst,sw)                                      
      goto 10941                                                        
11031 if(kri .ne. 10)goto 11041                                         
      call andarm10(n,y,z,w,dst,sw)                                     
      goto 10941                                                        
11041 if(kri .ne. 11)goto 11051                                         
      call andarm11(dst,sw)                                     
      goto 10941                                                        
11051 if(kri .ne. 12)goto 11061                                         
      call andarm12(n,y,z,w,dst,sw)                                     
      goto 10941                                                        
11061 if(kri .ne. 13)goto 11071                                         
      call andarm12(n,y,z,w,dst,sw)                                     
      goto 10941                                                        
11071 if(kri .ne. 14)goto 11081                                         
      call andarm14(n,y,z,w,dst,sw)                                     
      goto 11091                                                        
11081 continue                                                          
      call andarm15(n,y,y2,z,w,dst,sw)                                  
11091 continue                                                          
10941 continue                                                          
      return                                                            
      end                                                               
      subroutine set_kri(irg,jrg)                                       
      implicit double precision(a-h,o-z)                                
!DEC$ ATTRIBUTES DLLEXPORT :: set_kri                                   
      save kri                                                          
      if(jrg .ne. 1)goto 11111                                          
      kri=irg                                                           
      return                                                            
11111 continue                                                          
      irg=kri                                                           
      return                                                            
      end                                                               
      subroutine andarm11(dst,sw)                               
      implicit double precision(a-h,o-z)                                
      dst=0.0                                                           
      sw=dst                                                            
      return                                                            
      end                                                               
      subroutine andarm2(n,y,z,w,dst,sw)                                
      implicit double precision(a-h,o-z)                                
      parameter(nmin=100)                                               
      double precision y(n),z(n),w(n)                                   
      integer my(n),mz(n)                                               
      call set_qqtrm(itrm,2)                                            
      if(n .ge. nmin)goto 11131                                         
      dst=0.0                                                           
      sw=sum(w)                                                         
      return                                                            
11131 continue                                                          
      if(n .ge. 2*itrm)goto 11151                                       
      dst=0.0                                                           
      sw=dst                                                            
      return                                                            
11151 continue                                                          
11160 do 11161 i=1,n                                                    
      my(i)=i                                                           
11161 continue                                                          
11162 continue                                                          
      call psort8(y,my,1,n)                                             
11170 do 11171 i=1,n                                                    
      mz(i)=i                                                           
11171 continue                                                          
11172 continue                                                          
      call psort8(z,mz,1,n)                                             
      dst=0.0                                                           
      sw1=dst                                                           
11180 do 11181 i=itrm+1,n-itrm                                          
      sw1=sw1+w(my(i))                                                  
      dst=dst+w(my(i))*abs(y(my(i))-z(mz(i)))                           
11181 continue                                                          
11182 continue                                                          
      dst=dst/sw1                                                       
      sw=sum(w)                                                         
      return                                                            
      end                                                               
      subroutine set_qqtrm(irg,jrg)                                     
      implicit double precision(a-h,o-z)                                
!DEC$ ATTRIBUTES DLLEXPORT :: set_qqtrm                                 
      save itrm                                                         
      if(jrg .ne. 1)goto 11201                                          
      itrm=irg                                                          
      return                                                            
11201 continue                                                          
      irg=itrm                                                          
      return                                                            
      end                                                               
      subroutine andarm1(n,y,z,w,dst,sw)                                
      implicit double precision(a-h,o-z)                                
      parameter(eps=1.0e-5,nmin=100)                                    
      double precision y(n),z(n),w(n),q(2*n),w2(2*n)                    
      integer m(2*n),iq(2*n)                                            
      if(n .ge. nmin)goto 11221                                         
      dst=0.0                                                           
      sw=sum(w)                                                         
      return                                                            
11221 continue                                                          
      n2=2*n                                                            
11230 do 11231 i=1,n                                                    
      q(i)=y(i)                                                         
      iq(i)=0                                                           
      q(i+n)=z(i)                                                       
      iq(i+n)=1                                                         
      w2(i)=w(i)                                                        
      w2(i+n)=w(i)                                                      
11231 continue                                                          
11232 continue                                                          
11240 do 11241 i=1,n2                                                   
      m(i)=i                                                            
11241 continue                                                          
11242 continue                                                          
      call psort8(q,m,1,n2)                                             
      sw=0.0                                                            
      tw=sw                                                             
      dst=tw                                                            
      sw2=2.0*sum(w)                                                    
11250 do 11251 i=1,n2                                                   
      k=m(i)                                                            
      if(iq(k) .ne. 0)goto 11271                                        
      sw=sw+w2(k)                                                       
      goto 11281                                                        
11271 continue                                                          
      tw=tw+w2(k)                                                       
11281 continue                                                          
11261 continue                                                          
      pw=(sw+tw)*(sw2-sw-tw)                                            
      pw=max(eps,pw)                                                    
      dst=dst+abs(sw-tw)/sqrt(pw)                                       
11251 continue                                                          
11252 continue                                                          
      dst=dst/n                                                         
      return                                                            
      end                                                               
      subroutine andarm3(n,y,z,w,dst,sw)                                
      implicit double precision(a-h,o-z)                                
      double precision y(n),z(n),w(n)                                   
      sw=sum(w)                                                         
      dst=dot_product(w,abs(y-z))/sw                                    
      return                                                            
      end                                                               
      subroutine andarm7(n,y,z,w,dst,sw)                                
      implicit double precision(a-h,o-z)                                
      parameter(nmin=20)                                                
      double precision y(n),z(n),w(n)                                   
      if(n .ge. nmin)goto 11301                                         
      dst=0.0                                                           
      sw=sum(w)                                                         
      return                                                            
11301 continue                                                          
      sw=sum(w)                                                         
      dst=abs(dot_product(w,y)/sw-dot_product(w,z)/sw)                  
      return                                                            
      end                                                               
      subroutine andarm12(n,y,z,w,dst,sw)                               
      implicit double precision(a-h,o-z)                                
      parameter(fmin=20)                                                
      double precision y(n),z(n),w(n)                                   
      if(n .ge. 2*int(fmin))goto 11321                                  
      dst=0.0                                                           
      sw=sum(w)                                                         
      return                                                            
11321 continue                                                          
      dst1=0.0                                                          
      dst2=dst1                                                         
      sw1=dst2                                                          
      sw2=sw1                                                           
11330 do 11331 i=1,n                                                    
      if(z(i) .ge. 0.0)goto 11351                                       
      sw1=sw1+w(i)                                                      
      dst1=dst1+w(i)*y(i)                                               
      goto 11361                                                        
11351 continue                                                          
      sw2=sw2+w(i)                                                      
      dst2=dst2+w(i)*y(i)                                               
11361 continue                                                          
11341 continue                                                          
11331 continue                                                          
11332 continue                                                          
      sw=sum(w)                                                         
      if((n*sw1/sw .ge. fmin) .and. (n*sw2/sw .ge. fmin))goto 11381     
      dst=0.0                                                           
      return                                                            
11381 continue                                                          
      dst=abs(dst2/sw2-dst1/sw1)                                        
      return                                                            
      end                                                               
      subroutine andarm14(n,y,z,w,dst,sw)                               
      implicit double precision(a-h,o-z)                                
      parameter(fmin=20,sml=-1.0e20)                                    
      double precision y(n),z(n),w(n)                                   
      if(n .ge. 2*int(fmin))goto 11401                                  
      dst=sml                                                           
      sw=sum(w)                                                         
      return                                                            
11401 continue                                                          
      dst1=0.0                                                          
      dst2=dst1                                                         
      sw1=dst2                                                          
      sw2=sw1                                                           
11410 do 11411 i=1,n                                                    
      if(z(i) .ge. 0.0)goto 11431                                       
      sw1=sw1+w(i)                                                      
      dst1=dst1+w(i)*y(i)                                               
      goto 11441                                                        
11431 continue                                                          
      sw2=sw2+w(i)                                                      
      dst2=dst2+w(i)*y(i)                                               
11441 continue                                                          
11421 continue                                                          
11411 continue                                                          
11412 continue                                                          
      sw=sum(w)                                                         
      if((n*sw1/sw .ge. fmin) .and. (n*sw2/sw .ge. fmin))goto 11461     
      dst=sml                                                           
      return                                                            
11461 continue                                                          
      dst=dst2/sw2-dst1/sw1                                             
      return                                                            
      end                                                               
      subroutine andarm8(n,y,z,w,dst,sw)                                
      implicit double precision(a-h,o-z)                                
      parameter(nmin=20,sml=-1.0e20)                                    
      double precision y(n),z(n),w(n)                                   
      if(n .ge. nmin)goto 11481                                         
      dst=sml                                                           
      sw=sum(w)                                                         
      return                                                            
11481 continue                                                          
      sw=sum(w)                                                         
      dst=dot_product(w,y)/sw-dot_product(w,z)/sw                       
      return                                                            
      end                                                               
      subroutine andarm4(n,y,z,w,dst,sw)                                
      implicit double precision(a-h,o-z)                                
      parameter(maxclass2=10000,nmin=100)                               
      double precision y(n),z(n),w(n),out(maxclass2)                    
      double precision, dimension (:,:), allocatable :: costs           
      if(n .ge. nmin)goto 11501                                         
      dst=0.0                                                           
      sw=sum(w)                                                         
      return                                                            
11501 continue                                                          
      allocate(costs(1:nclass,1:nclass),stat=jerr);                     
      call classin(2,idum,costs,nclass,out)                               
      call reorg(2,nclass,out,costs)                                    
      dst=0.0                                                           
11510 do 11511 i=1,n                                                    
      ky=y(i)+0.1                                                       
      kz=z(i)+0.1                                                       
      dst=dst+w(i)*costs(ky,kz)                                         
11511 continue                                                          
11512 continue                                                          
      sw=sum(w)                                                         
      dst=dst/sw                                                        
      return                                                            
      end                                                               
      subroutine andarm5(n,y,z,w,dst,sw)                                
      implicit double precision(a-h,o-z)                                
      parameter(nmin=50)                                                
      double precision y(n),z(n),w(n)                                   
      data qntl /0.5/                                                   
      if(n .ge. nmin)goto 11531                                         
      dst=0.0                                                           
      sw=sum(w)                                                         
      return                                                            
11531 continue                                                          
      sw=sum(w)                                                         
      up=0.0                                                            
11540 do 11541 i=1,n                                                    
      if(y(i).le.z(i)) up=up+w(i)                                       
11541 continue                                                          
11542 continue                                                          
      dst=abs(up/sw-qntl)                                               
      return                                                            
      entry set_quant(arg)                                              
!DEC$ ATTRIBUTES DLLEXPORT :: set_quant                                 
      qntl=arg                                                          
      return                                                            
      end                                                               
      subroutine andarm6(n,y,y2,z,w,dst,sw)                             
      implicit double precision(a-h,o-z)                                
      parameter(nit=100,xmiss=9.0e35,thr=1.0e-2,nmin=100)               
      double precision y(n),y2(n),z(n),w(n),yy(n,2)                     
      if(n .ge. nmin)goto 11561                                         
      dst=0.0                                                           
      sw=sum(w)                                                         
      return                                                            
11561 continue                                                          
      yy(:,1)=y                                                         
      yy(:,2)=y2                                                        
      call cendst(n,yy,z,w,nit,thr,xmiss,dst,sw)                        
      sw=sum(w)                                                         
      return                                                            
      end                                                               
      subroutine andarm15(n,y,y2,z,w,dst,sw)                            
      implicit double precision(a-h,o-z)                                
      parameter(nit=100,xmiss=9.0e35,thr=1.0e-2,nmin=100)               
      double precision y(n),y2(n),z(n),w(n),yy(n,2)                     
      if(n .ge. nmin)goto 11581                                         
      dst=0.0                                                           
      sw=sum(w)                                                         
      return                                                            
11581 continue                                                          
      yy(:,1)=y                                                         
      yy(:,2)=y2                                                        
      call cendst1(n,yy,z,w,nit,thr,xmiss,dst,sw)                       
      sw=sum(w)                                                         
      return                                                            
      end                                                               
      subroutine andarm10(n,y,z,w,dst,sw)                               
      implicit double precision(a-h,o-z)                                
      parameter(eps=1.0e-5,nmin=100)                                    
      double precision y(n),z(n),w(n)                                   
      integer m(n)                                                      
      sw=sum(w)                                                         
      if(n .ge. nmin)goto 11601                                         
      dst=0.0                                                           
      return                                                            
11601 continue                                                          
      sw1=0.0                                                           
      sw2=sw1                                                           
11610 do 11611 i=1,n                                                    
      m(i)=i                                                            
      if(z(i) .ge. 0.0)goto 11631                                       
      sw1=sw1+w(i)                                                      
      goto 11641                                                        
11631 continue                                                          
      sw2=sw2+w(i)                                                      
11641 continue                                                          
11621 continue                                                          
11611 continue                                                          
11612 continue                                                          
      call psort8(y,m,1,n)                                              
      s1=0.0                                                            
      s2=s1                                                             
      s=s2                                                              
      dst=s                                                             
11650 do 11651 i=1,n                                                    
      k=m(i)                                                            
      s=s+w(k)                                                          
      if(z(k) .ge. 0)goto 11671                                         
      s1=s1+w(k)/sw1                                                    
      goto 11681                                                        
11671 continue                                                          
      s2=s2+w(k)/sw2                                                    
11681 continue                                                          
11661 continue                                                          
      pw=s*(sw-s)                                                       
      pw=max(eps,pw)                                                    
      dst=dst+abs(s1-s2)/sqrt(pw)                                       
11651 continue                                                          
11652 continue                                                          
      return                                                            
      end                                                               
      subroutine stput (iseed)                                          
      implicit double precision(a-h,o-z)                                
      double precision x(*)                                             
      data i /987654321/                                                
      i=iseed                                                           
      return                                                            
      entry rget (x,n)                                                  
      do 1 j=1,n                                                        
      i=mod(i*16807.0,2147483647.0)                                     
      u=i                                                               
      u=u*.465661287d-9                                                 
      x(j)=u
 1    continue
      return                                                            
      entry stget (irg)                                                 
      irg=i                                                             
      return                                                            
      end                                                               
      subroutine cendst(n,y,z,w,nit,thr,xmiss,dst,sw)                   
      implicit double precision(a-h,o-z)                                
      parameter(eps=0.1)                                                
      double precision y(n,2),z(n),w(n),b(2*n+1),q(3*n),cdf(3*n),r(n)   
      integer iq(3*n),mm(3*n),mz(n)                                     
      data nsamp /500/                                                  
      n2=2*n                                                            
      n3=3*n                                                            
      sw=sum(w)                                                         
11690 do 11691 i=1,n                                                    
      mz(i)=i                                                           
11691 continue                                                          
11692 continue                                                          
      call psort8(z,mz,1,n)                                             
      nq=0.25*n                                                         
      teps=(z(mz(n-nq))-z(mz(nq)))*eps                                  
11700 do 11701 i=1,n                                                    
      if(y(i,2)-y(i,1).ge.teps)goto 11701                               
      y(i,1)=y(i,1)-teps                                                
      y(i,2)=y(i,2)+teps                                                
11701 continue                                                          
11702 continue                                                          
11710 do 11711 i=1,n                                                    
      b(i)=y(i,1)                                                       
      b(i+n)=y(i,2)                                                     
11711 continue                                                          
11712 continue                                                          
      m=0                                                               
11720 do 11721 i=1,n2                                                   
      if(b(i).le.-xmiss)goto 11721                                      
      if(b(i).ge.xmiss)goto 11721                                       
      m=m+1                                                             
      b(m)=b(i)                                                         
11721 continue                                                          
11722 continue                                                          
      call unique(m,b,nu)                                               
      if(nu .le. nsamp)goto 11741                                       
      call rget(r,nsamp)                                                
11750 do 11751 i=1,nsamp                                                
      r(i)=b(int(nu*r(i))+1)                                            
11751 continue                                                          
11752 continue                                                          
      nu=nsamp                                                          
      b(1:nu)=r(1:nu)                                                   
      call sort(b,nu)                                                   
11741 continue                                                          
      m=nu+1                                                            
      b(m)=xmiss                                                        
      call fintcdf1(n,y,m,b,w,nit,thr/m,cdf,jt,err)                     
      m=m-1                                                             
11760 do 11761 i=1,m                                                    
      q(i)=b(i)                                                         
      iq(i)=0                                                           
11761 continue                                                          
11762 continue                                                          
11770 do 11771 i=1,n                                                    
      mz(i)=i                                                           
11771 continue                                                          
11772 continue                                                          
      call psort8(z,mz,1,n)                                             
      mpn=m+n                                                           
      k=0                                                               
11780 do 11781 i=m+1,mpn                                                
      k=k+1                                                             
      q(i)=z(mz(k))                                                     
      iq(i)=1                                                           
11781 continue                                                          
11782 continue                                                          
11790 do 11791 i=1,mpn                                                  
      mm(i)=i                                                           
11791 continue                                                          
11792 continue                                                          
      call psort8(q,mm,1,mpn)                                           
      ycdf=0.0                                                          
      zcdf=ycdf                                                         
      dst=zcdf                                                          
      spw=dst                                                           
      ny=0                                                              
      nz=ny                                                             
11800 do 11801 i=1,mpn                                                  
      k=mm(i)                                                           
      if(iq(k) .ne. 0)goto 11821                                        
      ny=ny+1                                                           
      ycdf=cdf(ny)                                                      
      pw=float(i)*float(mpn-i)/float(mpn)**2                            
      pw=max(eps,pw)                                                    
      pw=1.0/sqrt(pw)                                                   
      spw=spw+pw                                                        
      dst=dst+pw*abs(ycdf-zcdf)                                         
      goto 11831                                                        
11821 continue                                                          
      nz=nz+1                                                           
      zcdf=zcdf+w(nz)/sw                                                
11831 continue                                                          
11811 continue                                                          
11801 continue                                                          
11802 continue                                                          
      dst=dst/spw                                                       
      return                                                            
      entry set_samp(irg)                                               
!DEC$ ATTRIBUTES DLLEXPORT :: set_samp                                  
      nsamp=irg                                                         
      call set_samp1(nsamp)                                             
      return                                                            
      end                                                               
      subroutine cendst1(n,y,z,w,nit,thr,xmiss,dst,sw)                  
      implicit double precision(a-h,o-z)                                
      double precision y(n,2),z(n),w(n),b(2*n+1),cdf1(3*n),cdf2(3*n),r(n
     *)
      double precision y1(n,2),y2(n,2),w1(n),w2(n)                      
      data nsamp /500/                                                  
      n1=0.0                                                            
      n2=n1                                                             
11840 do 11841 i=1,n                                                    
      if(y(i,1).le.-xmiss)goto 11841                                    
      if(y(i,2).ge.xmiss)goto 11841                                     
      if(y(i,2)-y(i,1).ge.teps)goto 11841                               
      y(i,1)=y(i,1)-teps                                                
      y(i,2)=y(i,2)+teps                                                
11841 continue                                                          
11842 continue                                                          
11850 do 11851 i=1,n                                                    
      if(z(i) .ge. 0.0)goto 11871                                       
      n1=n1+1                                                           
      y1(n1,:)=y(i,:)                                                   
      w1(n1)=w(i)                                                       
      goto 11881                                                        
11871 continue                                                          
      n2=n2+1                                                           
      y2(n2,:)=y(i,:)                                                   
      w2(n2)=w(i)                                                       
11881 continue                                                          
11861 continue                                                          
11851 continue                                                          
11852 continue                                                          
11890 do 11891 i=1,n                                                    
      b(i)=y(i,1)                                                       
      b(i+n)=y(i,2)                                                     
11891 continue                                                          
11892 continue                                                          
      m=0                                                               
11900 do 11901 i=1,n2                                                   
      if(b(i).le.-xmiss)goto 11901                                      
      if(b(i).ge.xmiss)goto 11901                                       
      m=m+1                                                             
      b(m)=b(i)                                                         
11901 continue                                                          
11902 continue                                                          
      call unique(m,b,nu)                                               
      if(nu .le. nsamp)goto 11921                                       
      call rget(r,nsamp)                                                
11930 do 11931 i=1,nsamp                                                
      r(i)=b(int(nu*r(i))+1)                                            
11931 continue                                                          
11932 continue                                                          
      nu=nsamp                                                          
      b(1:nu)=r(1:nu)                                                   
      call sort(b,nu)                                                   
11921 continue                                                          
      m=nu+1                                                            
      b(m)=xmiss                                                        
      call getcdf1(n1,y1,w1,nit,thr,xmiss,nsamp,m,b,cdf1,sw1)           
      call getcdf1(n2,y2,w2,nit,thr,xmiss,nsamp,m,b,cdf2,sw2)           
      call diffcdf(m,cdf1,cdf2,dst)                                     
      return                                                            
      entry set_samp1(irg)                                              
      nsamp=irg                                                         
      return                                                            
      end                                                               
      subroutine getcdf1(n,y,w,nit,thr,xmiss,nsamp,m,b,cdf,sw)          
      implicit double precision(a-h,o-z)                                
      parameter(teps=0.01)                                              
      double precision y(n,2),z(n),w(n),b(2*n+1),q(3*n),cdf(3*n),r(n)   
      integer iq(3*n),mm(3*n),mz(n)                                     
      n2=2*n                                                            
      sw=sum(w)                                                         
      call fintcdf1(n,y,m,b,w,nit,thr/m,cdf,jt,err)                     
      m=m-1                                                             
      return                                                            
      end                                                               
      subroutine diffcdf(m,cdf1,cdf2,dst)                               
      implicit double precision(a-h,o-z)                                
      double precision cdf1(m),cdf2(m)                                  
      f12=sqrt(float(m))                                                
      dst=0.0                                                           
11940 do 11941 i=1,m                                                    
      dst=dst+abs(cdf1(i)-cdf2(i))/sqrt(float(i)*float(m-i+1))          
11941 continue                                                          
11942 continue                                                          
      dst=f12*dst/m                                                     
      return                                                            
      end                                                               
      subroutine unique(n,y,nu)                                         
      implicit double precision(a-h,o-z)                                
      double precision y(n),yu(n)                                       
      integer m(n)                                                      
11950 do 11951 i=1,n                                                    
      m(i)=i                                                            
11951 continue                                                          
11952 continue                                                          
      call psort8(y,m,1,n)                                              
      nu=1                                                              
      yu(1)=y(m(1))                                                     
11960 do 11961 i=2,n                                                    
      if(y(m(i-1)).ge.y(m(i)))goto 11961                                
      nu=nu+1                                                           
      yu(nu)=y(m(i))                                                    
11961 continue                                                          
11962 continue                                                          
      y(1:nu)=yu(1:nu)                                                  
      return                                                            
      end                                                               
      subroutine fintcdf1(n,y,m,b,w1,nit,thr,cdf,jt,err)                
      implicit double precision(a-h,o-z)                                
!DEC$ ATTRIBUTES DLLEXPORT :: fintcdf1                                  
      double precision y(n,2),b(m),w(n),w1(n),p(m),pij(n,m),ps(m),cdf(m)
      integer, dimension (:), allocatable :: ic,jc,kc1,kc2,lc           
      double precision, dimension (:), allocatable :: smo               
      call set_vrb(ivrb,2)                                              
      w=w1/sum(w1)                                                      
      p=1.0/m                                                           
      nt=0                                                              
11970 do 11971 i=1,n                                                    
11980 do 11981 k=1,m                                                    
      if(y(i,1).ge.b(k))goto 11981                                      
      if(y(i,2).lt.b(k))goto 11981                                      
      nt=nt+1                                                           
11981 continue                                                          
11982 continue                                                          
11971 continue                                                          
11972 continue                                                          
      allocate(ic(1:(n+1)),stat=jerr)                                   
      if(jerr.ne.0) return                                              
      allocate(jc(1:nt),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(kc1(1:m),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(kc2(1:m),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(smo(1:m),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(lc(1:nt),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      if(jerr .eq. 0)goto 12001                                         
      err=8888.0                                                        
      return                                                            
12001 continue                                                          
      nt=0                                                              
      ic(1)=1                                                           
12010 do 12011 i=1,n                                                    
12020 do 12021 k=1,m                                                    
      if(y(i,1).ge.b(k))goto 12021                                      
      if(y(i,2).lt.b(k))goto 12021                                      
      nt=nt+1                                                           
      jc(nt)=k                                                          
12021 continue                                                          
12022 continue                                                          
      ic(i+1)=nt+1                                                      
12011 continue                                                          
12012 continue                                                          
      nt=0                                                              
12030 do 12031 j=1,m                                                    
      kc1(j)=nt+1                                                       
12040 do 12041 i=1,n                                                    
      if(y(i,1).ge.b(j))goto 12041                                      
      if(y(i,2).lt.b(j))goto 12041                                      
      nt=nt+1                                                           
      lc(nt)=i                                                          
12041 continue                                                          
12042 continue                                                          
      kc2(j)=nt                                                         
12031 continue                                                          
12032 continue                                                          
      if(ivrb.gt.0) call labelpr('CDF iterations', -1)
c     if(ivrb.gt.0) write(6,'(''CDF iterations'',$)')                   
12050 do 12051 it=1,nit                                                 
      jt=it                                                             
      ps=p                                                              
12060 do 12061 j=1,m                                                    
      pij(:,j)=0.0                                                      
12070 do 12071 ii=kc1(j),kc2(j)                                         
      i=lc(ii)                                                          
      s=sum(p(jc(ic(i):(ic(i+1)-1))))                                   
      if(s .gt. 0.0)goto 12091                                          
      err=-7777.0                                                       
      return                                                            
12091 continue                                                          
      pij(i,j)=w(i)*p(j)/s                                              
12071 continue                                                          
12072 continue                                                          
      p(j)=sum(pij(:,j))                                                
12061 continue                                                          
12062 continue                                                          
      if(m .le. 100)goto 12111                                          
      smo(1)=(2.0*p(1)+p(2))/3.0                                        
      smo(m)=(2.0*p(m)+p(m-1))/3.0                                      
      smo(2)=0.25*(p(1)+2.0*p(2)+p(3))                                  
      smo(m-1)=0.25*(p(m)+2.0*p(m-1)+p(m-2))                            
12120 do 12121 j=3,m-2                                                  
      smo(j)=(p(j-2)+2.0*p(j-1)+3.0*p(j)+2.0*p(j+1)+p(j+2))/9.0         
12121 continue                                                          
12122 continue                                                          
      p=smo                                                             
12111 continue                                                          
      err=sum(abs(p-ps))/m                                              
      if(kbad(err) .le. 0)goto 12141                                    
      err=7777.0                                                        
      return                                                            
12141 continue                                                          
      if(err.lt.thr)goto 12052
      call labelpr('.', 1)
c     if(ivrb.gt.0) write(6,'(''.'',$)')                                
12051 continue                                                          
12052 continue                                                          
      cdf(1)=p(1)                                                       
12150 do 12151 j=2,m                                                    
      cdf(j)=cdf(j-1)+p(j)                                              
12151 continue                                                          
12152 continue                                                          
      if(ivrb .le. 0)goto 12171                                         
      call dblepr1('Err = ', -1, err)
c$$$  write(6,12180)err                                                 
c$$$12180 format (g10.2)                                                    
12171 continue                                                          
      return                                                            
      end                                                               
      subroutine set_vrb(irg,jrg)                                       
      implicit double precision(a-h,o-z)                                
!DEC$ ATTRIBUTES DLLEXPORT :: set_vrb                                   
      save ivrb                                                         
      if(jrg .ne. 1)goto 12201                                          
      ivrb=irg                                                          
      return                                                            
12201 continue                                                          
      irg=ivrb                                                          
      return                                                            
      end                                                               
      subroutine cdfpoints1(m,x,n,y,w,cdf)                              
      implicit double precision(a-h,o-z)                                
!DEC$ ATTRIBUTES DLLEXPORT :: cdfpoints1                                
      double precision x(m),y(n),w(n),cdf(m)                            
      i=1                                                               
      j=0                                                               
      sw=0.0                                                            
12210 continue                                                          
12211 continue                                                          
      j=j+1                                                             
      if(j.gt.m) go to 12220                                            
12230 continue                                                          
12231 if(y(i).gt.x(j))goto 12232                                        
      sw=sw+w(i)                                                        
      i=i+1                                                             
      if(i.gt.n)goto 12232                                              
      goto 12231                                                        
12232 continue                                                          
      if(i .le. n)goto 12251                                            
12260 do 12261 k=j,m                                                    
      cdf(k)=sw                                                         
12261 continue                                                          
12262 continue                                                          
      go to 12270                                                       
12251 continue                                                          
      cdf(j)=sw                                                         
      goto 12211                                                        
12212 continue                                                          
12270 continue                                                          
12220 continue                                                          
      cdf=cdf/sum(w)                                                    
      return                                                            
      end                                                               
      subroutine sort(x,n)                                              
      implicit double precision(a-h,o-z)                                
      double precision x(n),z(n)                                        
      integer m(n)                                                      
12280 do 12281 i=1,n                                                    
      m(i)=i                                                            
12281 continue                                                          
12282 continue                                                          
      z=x                                                               
      call psort8(z,m,1,n)                                              
12290 do 12291 i=1,n                                                    
      x(i)=z(m(i))                                                      
12291 continue                                                          
12292 continue                                                          
      return                                                            
      end                                                               
      function kbad(u)                                                  
      implicit double precision(a-h,o-z)                                
      kbad=0                                                            
      if(isnan(u).or.abs(u).ge.abs(huge(u))) kbad=1                     
      return                                                            
      end                                                               
      subroutine classin(ient,nclasssv,costssv,nout,out)                
      implicit double precision(a-h,o-z)                                
!DEC$ ATTRIBUTES DLLEXPORT :: classin                                   
      double precision costssv(nclasssv,nclasssv),out(1)                
      double precision, dimension (:), allocatable :: costs             
      save costs,nclass                                                 
      nq=nclasssv*nclasssv                                              
      allocate(costs(1:nq),stat=jerr)                                   
      if(jerr.ne.0) return                                              
      if(ient .ne. 1)goto 12311                                         
      nclass=nclasssv                                                   
      call reorg(1,nclass,costs,costssv)                                
      nout=1                                                            
      out=1.0                                                           
      goto 12321                                                        
12311 continue                                                          
      nout=nclass                                                       
      call reorg(2,nclass,costs,out)                                    
12321 continue                                                          
12301 continue                                                          
      return                                                            
      end                                                               
      subroutine reorg(ient,n,a,b)                                      
      implicit double precision(a-h,o-z)                                
      double precision a(n*n),b(n,n)                                    
      i=0                                                               
      if(ient .ne. 2)goto 12341                                         
12350 do 12351 k=1,n                                                    
12360 do 12361 j=1,n                                                    
      i=i+1                                                             
      b(j,k)=a(i)                                                       
12361 continue                                                          
12362 continue                                                          
12351 continue                                                          
12352 continue                                                          
      goto 12371                                                        
12341 continue                                                          
12380 do 12381 k=1,n                                                    
12390 do 12391 j=1,n                                                    
      i=i+1                                                             
      a(i)=b(j,k)                                                       
12391 continue                                                          
12392 continue                                                          
12381 continue                                                          
12382 continue                                                          
12371 continue                                                          
12331 continue                                                          
      return                                                            
      end                                                               
      subroutine crinode (itr,rtr,mxnodes,node,nodes,cri,wt)            
      implicit double precision(a-h,o-z)                                
!DEC$ ATTRIBUTES DLLEXPORT :: crinode                                   
      integer itr(6,*),nodes(mxnodes),m(mxnodes),ic(mxnodes)            
      double precision rtr(4,*),cri(mxnodes),wt(mxnodes),sc(mxnodes,2)  
      node=0                                                            
      k=itr(2,1)                                                        
12400 continue                                                          
12401 continue                                                          
      if(itr(4,k) .lt. 0)goto 12421                                     
      k=itr(2,k)                                                        
      goto 12401                                                        
12421 continue                                                          
      i1=itr(5,k)                                                       
      i2=itr(6,k)                                                       
      node=node+1                                                       
      if(node.gt.mxnodes) return                                        
      nodes(node)=k                                                     
      cri(node)=rtr(3,k)                                                
      wt(node)=rtr(4,k)                                                 
12430 continue                                                          
12431 if(k.eq.itr(2,iabs(itr(4,k))))goto 12432                          
      k=iabs(itr(4,k))                                                  
      if(k.eq.1)goto 12432                                              
      goto 12431                                                        
12432 continue                                                          
      if(k.eq.1)goto 12402                                              
      k=itr(3,iabs(itr(4,k)))                                           
      goto 12401                                                        
12402 continue                                                          
12440 do 12441 k=1,node                                                 
      m(k)=k                                                            
12441 continue                                                          
12442 continue                                                          
      call psort8(-cri,m,1,node)                                        
12450 do 12451 i=1,node                                                 
      ic(i)=nodes(m(i))                                                 
      sc(i,1)=cri(m(i))                                                 
      sc(i,2)=wt(m(i))                                                  
12451 continue                                                          
12452 continue                                                          
12460 do 12461 i=1,node                                                 
      nodes(i)=ic(i)                                                    
      cri(i)=sc(i,1)                                                    
      wt(i)=sc(i,2)                                                     
12461 continue                                                          
12462 continue                                                          
      return                                                            
      end                                                               
      subroutine prune1 (itr,rtr,nodes,thr,itro,rtro)                   
      implicit double precision(a-h,o-z)                                
!DEC$ ATTRIBUTES DLLEXPORT :: prune1                                    
      integer itr(6,nodes)                                              
      double precision rtr(4,nodes)                                     
      integer itro(6,nodes)                                             
      double precision rtro(4,nodes)                                    
      call prune(itr,rtr,nodes,thr)                                     
      itro=itr                                                          
      rtro=rtr                                                          
      return                                                            
      end                                                               
      subroutine prune (itr,rtr,nodes,thr)                              
      implicit double precision(a-h,o-z)                                
      integer itr(6,nodes)                                              
      double precision rtr(4,nodes)                                     
12470 continue                                                          
12471 continue                                                          
      nch=0                                                             
12480 do 12481 k=1,nodes                                                
      if(itr(4,k).le.0)goto 12481                                       
      nl=itr(2,k)                                                       
      nr=itr(3,k)                                                       
      if(itr(4,nl).ge.0)goto 12481                                      
      if(itr(4,nr).ge.0)goto 12481                                      
      if(max(rtr(3,nl),rtr(3,nr)).gt.rtr(3,k)+thr)goto 12481            
      itr(4,k)=-itr(4,k)                                                
      nch=nch+1                                                         
12481 continue                                                          
12482 continue                                                          
      if(nch.eq.0)goto 12472                                            
      goto 12471                                                        
12472 continue                                                          
      return                                                            
      end                                                               
      subroutine getnode (x,itre,rtre,cat,node)                         
      implicit double precision(a-h,o-z)                                
      integer itre(6,*)                                                 
      double precision x(*),rtre(4,*),cat(*)                            
      k=1                                                               
12490 continue                                                          
12491 if(itre(4,k).lt.0)goto 12492                                      
      if(itre(1,k) .le. 0)goto 12511                                    
      if(x(itre(1,k)) .ge. rtre(1,k))goto 12531                         
      k=itre(2,k)                                                       
      goto 12541                                                        
12531 continue                                                          
      k=itre(3,k)                                                       
12541 continue                                                          
12521 continue                                                          
      goto 12551                                                        
12511 continue                                                          
      j=-itre(1,k)                                                      
      kp=rtre(1,k)+0.1                                                  
      in=0                                                              
      np=abs(cat(kp))+0.1                                               
12560 do 12561 i=1,np                                                   
      if(x(j).ne.cat(kp+i))goto 12561                                   
      in=1                                                              
      goto 12562                                                        
12561 continue                                                          
12562 continue                                                          
      if(cat(kp) .le. 0.0)goto 12581                                    
      if(in .ne. 0)goto 12601                                           
      k=itre(2,k)                                                       
      goto 12611                                                        
12601 continue                                                          
      k=itre(3,k)                                                       
12611 continue                                                          
12591 continue                                                          
      goto 12621                                                        
12581 continue                                                          
      if(in .eq. 0)goto 12641                                           
      k=itre(2,k)                                                       
      goto 12651                                                        
12641 continue                                                          
      k=itre(3,k)                                                       
12651 continue                                                          
12631 continue                                                          
12621 continue                                                          
12571 continue                                                          
12551 continue                                                          
12501 continue                                                          
      goto 12491                                                        
12492 continue                                                          
      node=k                                                            
      return                                                            
      end                                                               
      subroutine getnodes1 (no,ni,x,itre,rtre,cat,nodes)                
      implicit double precision(a-h,o-z)                                
!DEC$ ATTRIBUTES DLLEXPORT :: getnodes1                                 
      integer itre(6,*),nodes(no)                                       
      double precision x(no,ni),rtre(4,*),cat(*)                        
12660 do 12661 i=1,no                                                   
      call getnode (x(i,:),itre,rtre,cat,node)                          
      nodes(i)=node                                                     
12661 continue                                                          
12662 continue                                                          
      return                                                            
      end                                                               
      subroutine getlims(node,ni,itr,rtr,cat,nvar,jvar,vlims,jerr)      
      implicit double precision(a-h,o-z)                                
!DEC$ ATTRIBUTES DLLEXPORT :: getlims                                   
      integer itr(6,*),jvar(2,*)                                        
      double precision rtr(4,*),cat(*),vlims(*)                         
      jerr=0                                                            
      if(itr(4,node) .lt. 0)goto 12681                                  
      jerr=1                                                            
      return                                                            
12681 continue                                                          
      nvar=0                                                            
      k=node                                                            
12690 continue                                                          
12691 continue                                                          
      nvar=nvar+1                                                       
      kpp=abs(itr(4,k))                                                 
      if(itr(1,kpp) .le. 0)goto 12711                                   
      jvar(2,nvar)=0                                                    
      if(itr(2,kpp) .ne. k)goto 12731                                   
      jvar(1,nvar)=-itr(1,kpp)                                          
      goto 12741                                                        
12731 continue                                                          
      jvar(1,nvar)=itr(1,kpp)                                           
12741 continue                                                          
12721 continue                                                          
      vlims(nvar)=rtr(1,kpp)                                            
      goto 12751                                                        
12711 continue                                                          
      if(k .ne. itr(2,kpp))goto 12771                                   
      sgn=-1.0                                                          
      goto 12781                                                        
12771 continue                                                          
      sgn=1.0                                                           
12781 continue                                                          
12761 continue                                                          
      jvar(1,nvar)=-itr(1,kpp)                                          
      kp=rtr(1,kpp)+0.1                                                 
      jvar(2,nvar)=kp                                                   
      vlims(nvar)=sgn*abs(cat(kp))                                      
12751 continue                                                          
12701 continue                                                          
      k=kpp                                                             
      if(kpp.eq.1)goto 12692                                            
      goto 12691                                                        
12692 continue                                                          
      return                                                            
      end                                                               
      subroutine trans(ny,y,wy,nz,z,wz,nt,t)                            
      implicit double precision(a-h,o-z)                                
!DEC$ ATTRIBUTES DLLEXPORT :: trans                                     
      double precision y(ny),wy(ny),z(nz),wz(nz),t(nt+2,2),u(max(ny,nz))
      double precision p(nt)                                            
      integer m(max(ny,nz))                                             
12790 do 12791 i=1,ny                                                   
      m(i)=i                                                            
      u(i)=y(i)                                                         
12791 continue                                                          
12792 continue                                                          
      call psort8(u,m,1,ny)                                             
12800 do 12801 i=1,ny                                                   
      y(i)=u(m(i))                                                      
12801 continue                                                          
12802 continue                                                          
      u=wy                                                              
12810 do 12811 i=1,ny                                                   
      wy(i)=u(m(i))                                                     
12811 continue                                                          
12812 continue                                                          
12820 do 12821 i=1,nz                                                   
      m(i)=i                                                            
      u(i)=z(i)                                                         
12821 continue                                                          
12822 continue                                                          
      call psort8(u,m,1,nz)                                             
12830 do 12831 i=1,nz                                                   
      z(i)=u(m(i))                                                      
12831 continue                                                          
12832 continue                                                          
      u=wz                                                              
12840 do 12841 i=1,nz                                                   
      wz(i)=u(m(i))                                                     
12841 continue                                                          
12842 continue                                                          
12850 do 12851 i=1,nt                                                   
      p(i)=(i-0.5)/float(nt)                                            
12851 continue                                                          
12852 continue                                                          
      call untie(ny,y,u)                                                
      call qntl(ny,u,wy,nt,p,t(:,1))                                    
      call untie(nz,z,u)                                                
      call qntl(nz,u,wz,nt,p,t(:,2))                                    
      return                                                            
      end                                                               
      subroutine qntl(n,y,w,nq,p,q)                                     
      implicit double precision(a-h,o-z)                                
      double precision y(n),w(n),p(nq),q(nq+2)                          
      sw=sum(w)                                                         
      k=1                                                               
      ff=w(1)                                                           
      q(1)=y(1)                                                         
      q(nq+2)=y(n)                                                      
12860 do 12861 i=2,n                                                    
      ff=ff+w(i)                                                        
      pp=ff/sw                                                          
      if(pp.lt.p(k))goto 12861                                          
      k=k+1                                                             
      q(k)=0.5*(y(i)+y(i-1))                                            
      if(k.ge.nq)goto 12862                                             
12861 continue                                                          
12862 continue                                                          
      q(nq+1)=0.5*(q(nq+2)+q(nq))                                       
      return                                                            
      end                                                               
      subroutine untie(n,y,u)                                           
      implicit double precision(a-h,o-z)                                
!DEC$ ATTRIBUTES DLLEXPORT :: untie                                     
      double precision y(n),u(n)                                        
      i=1                                                               
      k=0                                                               
12870 continue                                                          
12871 if(i.ge.n)goto 12872                                              
      if(y(i+1) .le. y(i))goto 12891                                    
      k=k+1                                                             
      u(k)=y(i)                                                         
      i=i+1                                                             
      goto 12871                                                        
12891 continue                                                          
      i1=i                                                              
12900 continue                                                          
12901 if(y(i+1).gt.y(i))goto 12902                                      
      i=i+1                                                             
      if(i.ge.n)goto 12902                                              
      goto 12901                                                        
12902 continue                                                          
      i2=i                                                              
      if(i1 .gt. 1)goto 12921                                           
      a=y(i1+1)                                                         
      b=y(i2+1)                                                         
      u(1)=y(1)                                                         
      k=1                                                               
12930 do 12931 j=i1+1,i2                                                
      k=k+1                                                             
      u(k)=a+(b-a)*(j-i1)/(i2-i1+1)                                     
12931 continue                                                          
12932 continue                                                          
      i=i2+1                                                            
      goto 12911                                                        
12921 if(i2 .lt. n)goto 12941                                           
      a=y(i1-1)                                                         
      b=(y(i2)-a)/(i2-i1+1)                                             
12950 do 12951 j=i1,i2                                                  
      k=k+1                                                             
      u(k)=a+b*(j-i1+1)                                                 
12951 continue                                                          
12952 continue                                                          
      goto 12961                                                        
12941 continue                                                          
      a=y(i1-1)                                                         
      b=y(i2)                                                           
12970 do 12971 j=i1,i2                                                  
      k=k+1                                                             
      u(k)=a+(b-a)*(j-i1+1)/(i2-i1+1)                                   
12971 continue                                                          
12972 continue                                                          
      i=i+1                                                             
12961 continue                                                          
12911 continue                                                          
      goto 12871                                                        
12872 continue                                                          
      if(k .ge. n)goto 12991                                            
      k=k+1                                                             
      u(k)=y(n)                                                         
12991 continue                                                          
      return                                                            
      end                                                               
      subroutine psort8 (v,a,ii,jj)                                     
      implicit double precision(a-h,o-z)                                
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
      double precision v                                                
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
