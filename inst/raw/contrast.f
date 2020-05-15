c     mortran 2.0     (version of 7/04/75 mod 7/4/87 (ajc))              
      subroutine fcontrast (no,ni,x,y,y2,z,w,lx,mxt,itre,rtre,mxc,cat,ms    258 
     *,isc)
!DEC$ ATTRIBUTES DLLEXPORT :: fcontrast                                         
      integer lx(ni),ms(no,ni,2),itre(6,mxt),isc(no)                        262
      real x(no,ni),y(no),y2(no),z(no),w(no),rtre(4,mxt),cat(mxc)           263
      save nodes,nct,kstor,lstor                                            264
      data mxtrm /10/                                                       265
      call dosort(no,ni,x,w,ms,nu,isc)                                      266
      lstor=0                                                               267
      itre(4,1)=0                                                           267
      itre(5,1)=1                                                           267
      itre(6,1)=nu                                                          267
      itre(2,1)=2                                                           267
      itre(3,1)=3                                                           268
      call andarm(no,y,y2,z,w,rtre(2,1),rtre(4,1))                          269
      nct=1                                                                 270
      call split7(no,ni,x,y,y2,z,w,lx,ms,1,nu,itre(1,1),rtre(1,1),  cri1    272 
     *,cri2,w1,w2,rtre(3,1),kct,cat(nct))
      if(itre(1,1) .ne. 0)goto 10021                                        272
      itre(4,1)=-9999                                                       272
      return                                                                272
10021 continue                                                              273
      if(itre(1,1) .ge. 0)goto 10041                                        273
      rtre(1,1)=nct                                                         273
      nct=nct+kct                                                           273
      if(nct.gt.mxc) lstor=1                                                273
10041 continue                                                              274
      rtre(2,2)=0.0                                                         274
      rtre(2,3)=rtre(2,2)                                                   274
      rtre(3,2)=cri1                                                        274
      rtre(3,3)=cri2                                                        275
      itre(4,2)=-1                                                          275
      itre(4,3)=itre(4,2)                                                   275
      rtre(4,2)=w1                                                          275
      rtre(4,3)=w2                                                          276
      rtre(2,1)=sign(max(0.0,rtre(3,1)-rtre(2,1)),rtre(3,1))                277
      if(lstor.ne.0) return                                                 278
      nodes=3                                                               278
      ntrm=0                                                                278
      if(rtre(2,1).gt.0.0) ntrm=1                                           278
      ktrg=1                                                                278
      kstor=0                                                               279
10050 continue                                                              279
10051 if(ntrm.ge.mxtrm)goto 10052                                           279
      jt=itre(1,ktrg)                                                       279
      st=rtre(1,ktrg)                                                       280
      k5=itre(5,ktrg)                                                       280
      k6=itre(6,ktrg)                                                       281
      if(jt .ge. 0)goto 10071                                               281
      ju=-jt                                                                281
      kp=st+0.1                                                             281
      np=abs(cat(kp))+0.1                                                   281
10071 continue                                                              282
10080 do 10081 j=1,ni                                                       282
      if(ni.gt.1.and.j.eq.jt)goto 10081                                     282
      kl=k5-1                                                               282
      kr=k6+1                                                               283
10090 do 10091 i=k5,k6                                                      283
      l=ms(i,j,1)                                                           284
      if(jt .le. 0)goto 10111                                               285
      if(x(l,jt) .ge. st)goto 10131                                         285
      kl=kl+1                                                               285
      isc(kl)=l                                                             285
      goto 10141                                                            286
10131 continue                                                              286
      kr=kr-1                                                               286
      isc(kr)=l                                                             286
10141 continue                                                              287
10121 continue                                                              287
      goto 10151                                                            288
10111 continue                                                              288
      in=0                                                                  289
10160 do 10161 im=1,np                                                      289
      if(x(l,ju).ne.cat(kp+im))goto 10161                                   289
      in=1                                                                  289
      goto 10162                                                            289
10161 continue                                                              290
10162 continue                                                              290
      if(cat(kp) .le. 0.0)goto 10181                                        291
      if(in .ne. 0)goto 10201                                               291
      kl=kl+1                                                               291
      isc(kl)=l                                                             291
      goto 10211                                                            292
10201 continue                                                              292
      kr=kr-1                                                               292
      isc(kr)=l                                                             292
10211 continue                                                              293
10191 continue                                                              293
      goto 10221                                                            294
10181 continue                                                              295
      if(in .eq. 0)goto 10241                                               295
      kl=kl+1                                                               295
      isc(kl)=l                                                             295
      goto 10251                                                            296
10241 continue                                                              296
      kr=kr-1                                                               296
      isc(kr)=l                                                             296
10251 continue                                                              297
10231 continue                                                              297
10221 continue                                                              298
10171 continue                                                              298
10151 continue                                                              299
10101 continue                                                              299
10091 continue                                                              300
10092 continue                                                              300
10260 do 10261 i=k5,kl                                                      300
      ms(i,j,1)=isc(i)                                                      300
10261 continue                                                              300
10262 continue                                                              300
10270 do 10271 i=kr,k6                                                      300
      ms(i,j,1)=isc(k6-i+kr)                                                300
10271 continue                                                              301
10272 continue                                                              301
10081 continue                                                              302
10082 continue                                                              302
      itre(5,itre(2,ktrg))=k5                                               302
      itre(6,itre(2,ktrg))=kl                                               303
      itre(5,itre(3,ktrg))=kr                                               303
      itre(6,itre(3,ktrg))=k6                                               304
      nde=itre(2,ktrg)                                                      305
      call split7(no,ni,x,y,y2,z,w,lx,ms,itre(5,nde),itre(6,nde),itre(1,    307 
     *nde),  rtre(1,nde),cri1,cri2,w1,w2,cri,kct,cat(nct))
      if(itre(1,nde) .ge. 0)goto 10291                                      307
      rtre(1,nde)=nct                                                       307
      nct=nct+kct                                                           307
      if(nct.gt.mxc) lstor=1                                                307
10291 continue                                                              308
      if(nodes+2 .le. mxt)goto 10311                                        308
      kstor=1                                                               308
      goto 10321                                                            309
10311 continue                                                              309
      itre(2,nde)=nodes+1                                                   309
      itre(3,nde)=nodes+2                                                   310
      l=nodes+1                                                             310
      rtre(2,l)=0.0                                                         310
      rtre(2,l+1)=rtre(2,l)                                                 310
      rtre(3,l)=cri1                                                        310
      rtre(3,l+1)=cri2                                                      311
      rtre(4,l)=w1                                                          311
      rtre(4,l+1)=w2                                                        311
      itre(4,l)=-nde                                                        311
      itre(4,l+1)=itre(4,l)                                                 312
      rtre(2,nde)=sign(max(0.0,cri-rtre(3,nde)),cri)                        313
10321 continue                                                              314
10301 continue                                                              314
      nde=itre(3,ktrg)                                                      315
      call split7(no,ni,x,y,y2,z,w,lx,ms,itre(5,nde),itre(6,nde),itre(1,    317 
     *nde),  rtre(1,nde),cri1,cri2,w1,w2,cri,kct,cat(nct))
      if(itre(1,nde) .ge. 0)goto 10341                                      317
      rtre(1,nde)=nct                                                       317
      nct=nct+kct                                                           317
      if(nct.gt.mxc) lstor=1                                                317
10341 continue                                                              318
      if(nodes+4 .le. mxt)goto 10361                                        318
      kstor=1                                                               318
      goto 10371                                                            319
10361 continue                                                              319
      itre(2,nde)=nodes+3                                                   319
      itre(3,nde)=nodes+4                                                   320
      l=nodes+3                                                             320
      rtre(2,l)=0.0                                                         320
      rtre(2,l+1)=rtre(2,l)                                                 320
      rtre(3,l)=cri1                                                        320
      rtre(3,l+1)=cri2                                                      321
      rtre(4,l)=w1                                                          321
      rtre(4,l+1)=w2                                                        321
      itre(4,l)=-nde                                                        321
      itre(4,l+1)=itre(4,l)                                                 322
      rtre(2,nde)=sign(max(0.0,cri-rtre(3,nde)),cri)                        323
      itre(4,ktrg)=-itre(4,ktrg)                                            323
      rsv=rtre(2,ktrg)                                                      323
      crix=0.0                                                              324
10380 do 10381 k=1,nodes                                                    324
      if(itre(4,k).ge.0)goto 10381                                          324
      if(abs(rtre(2,k)).le.crix)goto 10381                                  325
      if(itre(1,k).eq.0)goto 10381                                          325
      crix=abs(rtre(2,k))                                                   325
      ktrg=k                                                                326
10381 continue                                                              327
10382 continue                                                              327
10371 continue                                                              328
10351 continue                                                              328
      if(crix.le.0.0.or.kstor.ne.0.or.lstor.ne.0) return                    329
      nodes=nodes+4                                                         329
      if(rsv.gt.0.0) ntrm=ntrm+1                                            330
      goto 10051                                                            331
10052 continue                                                              331
      return                                                                332
      entry set_trm(irg)                                                    333
!DEC$ ATTRIBUTES DLLEXPORT :: set_trm                                           
      mxtrm=max0(irg,2)                                                     336
      return                                                                337
      entry get_stor(irg1,irg2)                                             338
!DEC$ ATTRIBUTES DLLEXPORT :: get_stor                                          
      irg1=nodes                                                            341
      irg2=nct-1                                                            341
      return                                                                342
      entry get_err(irg1,irg2)                                              342
      irg1=kstor                                                            342
      irg2=lstor                                                            342
      return                                                                343
      end                                                                   344
      subroutine ans (x,itre,rtre,cat,yh)                                   345
      integer itre(6,*)                                                     345
      real x(*),rtre(4,*),cat(*)                                            346
      k=1                                                                   347
10390 continue                                                              347
10391 if(itre(4,k).lt.0)goto 10392                                          347
      if(itre(1,k) .le. 0)goto 10411                                        349
      if(x(itre(1,k)) .ge. rtre(1,k))goto 10431                             349
      k=itre(2,k)                                                           349
      goto 10441                                                            349
10431 continue                                                              349
      k=itre(3,k)                                                           349
10441 continue                                                              350
10421 continue                                                              350
      goto 10451                                                            351
10411 continue                                                              351
      j=-itre(1,k)                                                          351
      kp=rtre(1,k)+0.1                                                      351
      in=0                                                                  351
      np=abs(cat(kp))+0.1                                                   352
10460 do 10461 i=1,np                                                       352
      if(x(j).ne.cat(kp+i))goto 10461                                       352
      in=1                                                                  352
      goto 10462                                                            352
10461 continue                                                              353
10462 continue                                                              353
      if(cat(kp) .le. 0.0)goto 10481                                        353
      if(in .ne. 0)goto 10501                                               353
      k=itre(2,k)                                                           353
      goto 10511                                                            353
10501 continue                                                              353
      k=itre(3,k)                                                           353
10511 continue                                                              353
10491 continue                                                              353
      goto 10521                                                            354
10481 continue                                                              354
      if(in .eq. 0)goto 10541                                               354
      k=itre(2,k)                                                           354
      goto 10551                                                            354
10541 continue                                                              354
      k=itre(3,k)                                                           354
10551 continue                                                              354
10531 continue                                                              354
10521 continue                                                              355
10471 continue                                                              355
10451 continue                                                              356
10401 continue                                                              356
      goto 10391                                                            357
10392 continue                                                              357
      yh=rtre(3,k)                                                          358
      return                                                                359
      end                                                                   360
      subroutine dosort(no,ni,x,w,ms,nu,isc)                                361
      integer ms(no,ni,*),isc(no)                                           361
      real x(no,ni),w(no)                                                   362
      data new /1/                                                          363
      if(new .eq. 0)goto 10571                                              364
10580 do 10581 j=1,ni                                                       364
10590 do 10591 i=1,no                                                       364
      ms(i,j,1)=i                                                           364
10591 continue                                                              364
10592 continue                                                              364
      call psort8(x(1,j),ms(1,j,1),1,no)                                    364
10581 continue                                                              365
10582 continue                                                              365
10600 do 10601 j=1,ni                                                       365
10610 do 10611 i=1,no                                                       365
      ms(i,j,2)=ms(i,j,1)                                                   365
10611 continue                                                              365
10612 continue                                                              365
      nu=0                                                                  366
10620 do 10621 i=1,no                                                       366
      if(w(ms(i,j,2)).le.0.0)goto 10621                                     366
      nu=nu+1                                                               366
      ms(nu,j,1)=ms(i,j,2)                                                  366
10621 continue                                                              367
10622 continue                                                              367
10601 continue                                                              368
10602 continue                                                              368
      return                                                                369
10571 continue                                                              370
10630 do 10631 j=1,ni                                                       370
      nu=0                                                                  371
10640 do 10641 i=1,no                                                       371
      if(w(ms(i,j,2)).le.0.0)goto 10641                                     371
      nu=nu+1                                                               371
      ms(nu,j,1)=ms(i,j,2)                                                  371
10641 continue                                                              372
10642 continue                                                              372
10631 continue                                                              373
10632 continue                                                              373
      return                                                                374
      entry set_new(irg)                                                    374
      new=irg                                                               374
      return                                                                375
      end                                                                   376
      subroutine split7 (no,ni,x,y,y2,z,w,lx,m,m1,m2,  jt,sp,cri1,cri2,w    378 
     *1,w2,crm,kct,cat)
      parameter(maxcat=1000, big=9.9e35)                                    379
      integer lx(ni),m(no,ni)                                               380
      real x(no,ni),y(no),y2(no),z(no),w(no),cat(*),tcat(maxcat)            381
      data xmiss,ntn,pwr /9.0e35,500,2/                                     382
      crm=0.0                                                               382
      jt=0                                                                  382
      if(m2-m1+1.lt.2*ntn) return                                           383
10650 do 10651 j=1,ni                                                       383
      if(x(m(m1,j),j).ge.x(m(m2,j),j).or.lx(j).eq.0)goto 10651              384
      if(lx(j) .ne. 1)goto 10671                                            385
      call eav(x(1,j),y,y2,z,w,m(1,j),m1,m2,ntn,pwr,xmiss,tp,  cril,criu    387 
     *,wl,wu,cri)
      if(cri.eq.-xmiss)goto 10651                                           388
      if(abs(cri) .lt. crm)goto 10691                                       388
      crm=abs(cri)                                                          388
      cri1=cril                                                             388
      cri2=criu                                                             389
      w1=wl                                                                 389
      w2=wu                                                                 389
      jt=j                                                                  389
      sp=tp                                                                 389
      crx=cri                                                               390
10691 continue                                                              391
      goto 10701                                                            392
10671 continue                                                              393
      call ceav(x(1,j),y,y2,z,w,m(1,j),m1,m2,ntn,pwr,lct,  tcat,cril,cri    395 
     *u,wl,wu,cri)
      if(cri .lt. crm)goto 10721                                            395
      crm=cri                                                               395
      cri1=cril                                                             395
      cri2=criu                                                             396
      w1=wl                                                                 396
      w2=wu                                                                 396
      jt=-j                                                                 396
      kct=lct                                                               396
10730 do 10731 k=1,kct                                                      396
      cat(k)=tcat(k)                                                        396
10731 continue                                                              397
10732 continue                                                              397
10721 continue                                                              398
10701 continue                                                              399
10661 continue                                                              399
10651 continue                                                              400
10652 continue                                                              400
      if(jt.gt.0) crm=crx                                                   401
      return                                                                402
      entry set_miss(arg)                                                   403
!DEC$ ATTRIBUTES DLLEXPORT :: set_miss                                          
      xmiss=arg                                                             406
      return                                                                407
      entry set_ntn(irg)                                                    408
!DEC$ ATTRIBUTES DLLEXPORT :: set_ntn                                           
      ntn=max0(irg,3)                                                       411
      return                                                                412
      entry set_pwr(arg)                                                    413
!DEC$ ATTRIBUTES DLLEXPORT :: set_pwr                                           
      pwr=arg                                                               416
      return                                                                417
      end                                                                   418
      subroutine eav (x,y,y2,z,w,m,m1,m2,ntn,pwr,xmiss,tp,  cri1s,cri2s,    420 
     *w1s,w2s,cri)
      integer m(*)                                                          420
      real x(*),y(*),y2(*),z(*),w(*)                                        421
      data nint,icri /4,1/                                                  422
      xl=x(m(m1))                                                           422
      if(xl .lt. xmiss)goto 10751                                           422
      cri=-xmiss                                                            422
      return                                                                422
10751 continue                                                              422
      nu=m2                                                                 423
      k=m(nu)                                                               423
10760 continue                                                              423
10761 if(x(k).lt.xmiss)goto 10762                                           423
      nu=nu-1                                                               423
      k=m(nu)                                                               423
      goto 10761                                                            424
10762 continue                                                              424
      crx=0.0                                                               424
      crimx=crx                                                             425
      if(nu .ge. m2)goto 10781                                              426
      call andarm(nu-m1+1,y(m(m1:nu)),y2(m(m1:nu)),  z(m(m1:nu)),w(m(m1:    428 
     *nu)),cri1xs,w1x)
      call andarm(m2-nu,y(m((nu+1):m2)),y2(m((nu+1):m2)),  z(m((nu+1):m2    430 
     *)),w(m((nu+1):m2)),cri2xs,w2x)
      crx=max(cri1xs,cri2xs)                                                431
      crimx=(crx**pwr)*float(nu-m1+1)*float(m2-nu)/float(m2-m1+1)**2        432
10781 continue                                                              433
      rq=(nu-m1+1)                                                          433
      mq=rq/nint                                                            433
      kq=1                                                                  433
      crim=-xmiss                                                           433
      cris=-xmiss                                                           434
10790 do 10791 i=nu,m1+1,-1                                                 434
      k=m(i)                                                                435
      tq=0.5*(x(k)+x(m(i-1)))                                               436
      if(tq.le.x(m(i-1)))goto 10791                                         436
      if(tq.ge.x(k))goto 10791                                              437
      if(i-m1.lt.ntn.or.nu-i+1.lt.ntn)goto 10791                            438
      if(i .ge. nu-kq*mq+1)goto 10811                                       438
      kq=kq+1                                                               439
      call andarm(i-m1,y(m(m1:(i-1))),y2(m(m1:(i-1))),  z(m(m1:(i-1))),w    441 
     *(m(m1:(i-1))),cri1,w1)
      call andarm(nu-i+1,y(m(i:nu)),y2(m(i:nu)),  z(m(i:nu)),w(m(i:nu)),    443 
     *cri2,w2)
      if(icri .ne. 1)goto 10831                                             443
      cri=max(cri1,cri2)                                                    443
      goto 10841                                                            444
10831 continue                                                              444
      cri=abs(cri1-cri2)                                                    444
10841 continue                                                              445
10821 continue                                                              445
      crit=(cri**pwr)*float(i-m1)*float(nu-i+1)/float(nu-m1+1)**2           446
      if(crit.lt.crim)goto 10791                                            447
      crim=crit                                                             447
      tp=tq                                                                 448
      cris=cri                                                              448
      w1s=w1                                                                448
      w2s=w2                                                                448
      cri1s=cri1                                                            448
      cri2s=cri2                                                            449
10811 continue                                                              450
10791 continue                                                              451
10792 continue                                                              451
      cri=cris                                                              452
      if(x(m(m2)).lt.xmiss) return                                          453
      tp=xmiss                                                              453
      cri1s=cri1xs                                                          454
      cri2s=cri2xs                                                          454
      w1s=w1x                                                               454
      w2s=w2x                                                               455
      if(crx .le. cri)goto 10861                                            455
      cri=crx                                                               455
      return                                                                455
10861 continue                                                              456
      cri=-cri                                                              457
      return                                                                458
      entry set_qint(irg)                                                   459
!DEC$ ATTRIBUTES DLLEXPORT :: set_qint                                          
      nint=irg                                                              462
      return                                                                463
      entry set_cri(irg)                                                    464
!DEC$ ATTRIBUTES DLLEXPORT :: set_cri                                           
      icri=irg                                                              467
      return                                                                468
      end                                                                   469
      subroutine ceav (x,y,y2,z,w,m,m1,m2,ntn,pwr,kct,cat,cri1,cri2,w1,w    470 
     *2,cri)
      parameter(maxcat=1000)                                                471
      integer m(*),mt(maxcat)                                               472
      real x(*),y(*),y2(*),z(*),w(*),cat(*),v(maxcat,3)                     473
      fntn=ntn                                                              473
      l=0                                                                   473
      i1=m1                                                                 474
10870 do 10871 i=m1,m2-1                                                    474
      k=m(i)                                                                474
      if(x(m(i+1)).le.x(k))goto 10871                                       475
      l=l+1                                                                 475
      v(l,1)=x(k)                                                           475
      i2=i-1                                                                476
      call andarm(i2-i1+1,y(m(i1:i2)),y2(m(i1:i2)),  z(m(i1:i2)),w(m(i1:    478 
     *i2)),v(l,2),v(l,3))
      i1=i                                                                  479
10871 continue                                                              480
10872 continue                                                              480
      k=m(m2)                                                               480
      l=l+1                                                                 480
      v(l,1)=x(k)                                                           481
      call andarm(m2-i1+1,y(m(i1:m2)),y2(m(i1:m2)),z(m(i1:m2)),w(m(i1:m2    483 
     *)),  v(l,2),v(l,3))
10880 do 10881 i=1,l                                                        483
      mt(i)=i                                                               483
10881 continue                                                              483
10882 continue                                                              483
      call psort8(v(1:l,2),mt,1,l)                                          484
10890 do 10891 j=1,l                                                        484
      v(j,2)=v(j,2)*v(j,3)                                                  484
10891 continue                                                              485
10892 continue                                                              485
      sl=0.0                                                                485
      wl=sl                                                                 485
      cri=wl                                                                485
      scri=cri                                                              485
      sr=sum(v(1:l,2))                                                      485
      wr=sum(v(1:l,3))                                                      485
      kct=0                                                                 486
10900 do 10901 i=1,l-1                                                      486
      k=mt(i)                                                               486
      sl=sl+v(k,2)                                                          486
      sr=sr-v(k,2)                                                          487
      wl=wl+v(k,3)                                                          487
      wr=wr-v(k,3)                                                          488
      if(wl.lt.fntn.or.wr.lt.fntn)goto 10901                                489
      c=wr*wl*max(sr/wr,sl/wl)**pwr                                         489
      if(c.le.cri)goto 10901                                                490
      cri=c                                                                 490
      kct=i                                                                 490
      cri1=sl/wl                                                            490
      cri2=sr/wr                                                            490
      w1=wl                                                                 490
      w2=wr                                                                 490
      scri=max(cri1,cri2)                                                   491
10901 continue                                                              492
10902 continue                                                              492
      if(kct .ne. 0)goto 10921                                              492
      cri=0.0                                                               492
      return                                                                492
10921 continue                                                              493
      cat(1)=-kct                                                           493
10930 do 10931 i=1,kct                                                      493
      cat(i+1)=v(mt(i),1)                                                   493
10931 continue                                                              494
10932 continue                                                              494
      cri=scri                                                              494
      kct=kct+1                                                             495
      return                                                                496
      end                                                                   497
      subroutine andarm(n,y,y2,z,w,dst,sw)                                  498
!DEC$ ATTRIBUTES DLLEXPORT :: andarm                                            
      real y(n),y2(n),z(n),w(n)                                             502
      call set_kri(kri,2)                                                   503
      if(kri .ne. 1)goto 10951                                              503
      call andarm1(n,y,z,w,dst,sw)                                          503
      goto 10941                                                            504
10951 if(kri .ne. 2)goto 10961                                              504
      call andarm2(n,y,z,w,dst,sw)                                          504
      goto 10941                                                            505
10961 if(kri .ne. 3)goto 10971                                              505
      call andarm3(n,y,z,w,dst,sw)                                          505
      goto 10941                                                            506
10971 if(kri .ne. 4)goto 10981                                              506
      call andarm4(n,y,z,w,dst,sw)                                          506
      goto 10941                                                            507
10981 if(kri .ne. 5)goto 10991                                              507
      call andarm5(n,y,z,w,dst,sw)                                          507
      goto 10941                                                            508
10991 if(kri .ne. 6)goto 11001                                              508
      call andarm6(n,y,y2,z,w,dst,sw)                                       508
      goto 10941                                                            509
11001 if(kri .ne. 7)goto 11011                                              509
      call andarm7(n,y,z,w,dst,sw)                                          509
      goto 10941                                                            510
11011 if(kri .ne. 8)goto 11021                                              510
      call andarm8(n,y,z,w,dst,sw)                                          510
      goto 10941                                                            511
11021 if(kri .ne. 9)goto 11031                                              511
      call andarm7(n,y,z,w,dst,sw)                                          511
      goto 10941                                                            512
11031 if(kri .ne. 10)goto 11041                                             512
      call andarm10(n,y,z,w,dst,sw)                                         512
      goto 10941                                                            513
11041 if(kri .ne. 11)goto 11051                                             513
      call andarm11(n,y,z,w,dst,sw)                                         513
      goto 10941                                                            514
11051 if(kri .ne. 12)goto 11061                                             514
      call andarm12(n,y,z,w,dst,sw)                                         514
      goto 10941                                                            515
11061 if(kri .ne. 13)goto 11071                                             515
      call andarm12(n,y,z,w,dst,sw)                                         515
      goto 10941                                                            516
11071 if(kri .ne. 14)goto 11081                                             516
      call andarm14(n,y,z,w,dst,sw)                                         516
      goto 11091                                                            517
11081 continue                                                              517
      call andarm15(n,y,y2,z,w,dst,sw)                                      517
11091 continue                                                              518
10941 continue                                                              518
      return                                                                519
      end                                                                   520
      subroutine set_kri(irg,jrg)                                           521
!DEC$ ATTRIBUTES DLLEXPORT :: set_kri                                           
      save kri                                                              525
      if(jrg .ne. 1)goto 11111                                              525
      kri=irg                                                               525
      return                                                                525
11111 continue                                                              526
      irg=kri                                                               527
      return                                                                528
      end                                                                   529
      subroutine andarm11(n,y,z,w,dst,sw)                                   531
      dst=0.0                                                               531
      sw=dst                                                                531
      return                                                                531
      end                                                                   532
      subroutine andarm2(n,y,z,w,dst,sw)                                    533
      parameter(nmin=100)                                                   534
      real y(n),z(n),w(n)                                                   534
      integer my(n),mz(n)                                                   535
      call set_qqtrm(itrm,2)                                                536
      if(n .ge. nmin)goto 11131                                             536
      dst=0.0                                                               536
      sw=sum(w)                                                             536
      return                                                                536
11131 continue                                                              537
      if(n .ge. 2*itrm)goto 11151                                           537
      dst=0.0                                                               537
      sw=dst                                                                537
      return                                                                537
11151 continue                                                              538
11160 do 11161 i=1,n                                                        538
      my(i)=i                                                               538
11161 continue                                                              538
11162 continue                                                              538
      call psort8(y,my,1,n)                                                 539
11170 do 11171 i=1,n                                                        539
      mz(i)=i                                                               539
11171 continue                                                              539
11172 continue                                                              539
      call psort8(z,mz,1,n)                                                 540
      dst=0.0                                                               540
      sw1=dst                                                               541
11180 do 11181 i=itrm+1,n-itrm                                              541
      sw1=sw1+w(my(i))                                                      542
      dst=dst+w(my(i))*abs(y(my(i))-z(mz(i)))                               543
11181 continue                                                              544
11182 continue                                                              544
      dst=dst/sw1                                                           544
      sw=sum(w)                                                             545
      return                                                                546
      end                                                                   547
      subroutine set_qqtrm(irg,jrg)                                         548
!DEC$ ATTRIBUTES DLLEXPORT :: set_qqtrm                                         
      save itrm                                                             552
      if(jrg .ne. 1)goto 11201                                              552
      itrm=irg                                                              552
      return                                                                552
11201 continue                                                              553
      irg=itrm                                                              554
      return                                                                555
      end                                                                   556
      subroutine andarm1(n,y,z,w,dst,sw)                                    557
      parameter(eps=1.0e-5,nmin=100)                                        558
      real y(n),z(n),w(n),q(2*n),w2(2*n)                                    558
      integer m(2*n),iq(2*n)                                                559
      if(n .ge. nmin)goto 11221                                             559
      dst=0.0                                                               559
      sw=sum(w)                                                             559
      return                                                                559
11221 continue                                                              560
      n2=2*n                                                                561
11230 do 11231 i=1,n                                                        561
      q(i)=y(i)                                                             561
      iq(i)=0                                                               561
      q(i+n)=z(i)                                                           562
      iq(i+n)=1                                                             562
      w2(i)=w(i)                                                            562
      w2(i+n)=w(i)                                                          563
11231 continue                                                              564
11232 continue                                                              564
11240 do 11241 i=1,n2                                                       564
      m(i)=i                                                                564
11241 continue                                                              564
11242 continue                                                              564
      call psort8(q,m,1,n2)                                                 565
      sw=0.0                                                                565
      tw=sw                                                                 565
      dst=tw                                                                565
      sw2=2.0*sum(w)                                                        566
11250 do 11251 i=1,n2                                                       566
      k=m(i)                                                                567
      if(iq(k) .ne. 0)goto 11271                                            567
      sw=sw+w2(k)                                                           567
      goto 11281                                                            567
11271 continue                                                              567
      tw=tw+w2(k)                                                           567
11281 continue                                                              568
11261 continue                                                              568
      pw=(sw+tw)*(sw2-sw-tw)                                                569
      pw=max(eps,pw)                                                        570
      dst=dst+abs(sw-tw)/sqrt(pw)                                           571
11251 continue                                                              572
11252 continue                                                              572
      dst=dst/n                                                             573
      return                                                                574
      end                                                                   575
      subroutine andarm3(n,y,z,w,dst,sw)                                    576
      real y(n),z(n),w(n)                                                   577
      sw=sum(w)                                                             577
      dst=dot_product(w,abs(y-z))/sw                                        578
      return                                                                579
      end                                                                   580
      subroutine andarm7(n,y,z,w,dst,sw)                                    581
      parameter(nmin=20)                                                    582
      real y(n),z(n),w(n)                                                   583
      if(n .ge. nmin)goto 11301                                             583
      dst=0.0                                                               583
      sw=sum(w)                                                             583
      return                                                                583
11301 continue                                                              584
      sw=sum(w)                                                             584
      dst=abs(dot_product(w,y)/sw-dot_product(w,z)/sw)                      585
      return                                                                586
      end                                                                   587
      subroutine andarm12(n,y,z,w,dst,sw)                                   588
      parameter(fmin=20)                                                    589
      real y(n),z(n),w(n)                                                   590
      if(n .ge. 2*int(fmin))goto 11321                                      590
      dst=0.0                                                               590
      sw=sum(w)                                                             590
      return                                                                590
11321 continue                                                              591
      dst1=0.0                                                              591
      dst2=dst1                                                             591
      sw1=dst2                                                              591
      sw2=sw1                                                               592
11330 do 11331 i=1,n                                                        593
      if(z(i) .ge. 0.0)goto 11351                                           593
      sw1=sw1+w(i)                                                          593
      dst1=dst1+w(i)*y(i)                                                   593
      goto 11361                                                            594
11351 continue                                                              594
      sw2=sw2+w(i)                                                          594
      dst2=dst2+w(i)*y(i)                                                   594
11361 continue                                                              595
11341 continue                                                              595
11331 continue                                                              596
11332 continue                                                              596
      sw=sum(w)                                                             597
      if((n*sw1/sw .ge. fmin) .and. (n*sw2/sw .ge. fmin))goto 11381         597
      dst=0.0                                                               597
      return                                                                597
11381 continue                                                              598
      dst=abs(dst2/sw2-dst1/sw1)                                            599
      return                                                                600
      end                                                                   601
      subroutine andarm14(n,y,z,w,dst,sw)                                   602
      parameter(fmin=20,sml=-1.0e20)                                        603
      real y(n),z(n),w(n)                                                   604
      if(n .ge. 2*int(fmin))goto 11401                                      604
      dst=sml                                                               604
      sw=sum(w)                                                             604
      return                                                                604
11401 continue                                                              605
      dst1=0.0                                                              605
      dst2=dst1                                                             605
      sw1=dst2                                                              605
      sw2=sw1                                                               606
11410 do 11411 i=1,n                                                        607
      if(z(i) .ge. 0.0)goto 11431                                           607
      sw1=sw1+w(i)                                                          607
      dst1=dst1+w(i)*y(i)                                                   607
      goto 11441                                                            608
11431 continue                                                              608
      sw2=sw2+w(i)                                                          608
      dst2=dst2+w(i)*y(i)                                                   608
11441 continue                                                              609
11421 continue                                                              609
11411 continue                                                              610
11412 continue                                                              610
      sw=sum(w)                                                             611
      if((n*sw1/sw .ge. fmin) .and. (n*sw2/sw .ge. fmin))goto 11461         611
      dst=sml                                                               611
      return                                                                611
11461 continue                                                              612
      dst=dst2/sw2-dst1/sw1                                                 613
      return                                                                614
      end                                                                   615
      subroutine andarm8(n,y,z,w,dst,sw)                                    616
      parameter(nmin=20,sml=-1.0e20)                                        617
      real y(n),z(n),w(n)                                                   618
      if(n .ge. nmin)goto 11481                                             618
      dst=sml                                                               618
      sw=sum(w)                                                             618
      return                                                                618
11481 continue                                                              619
      sw=sum(w)                                                             619
      dst=dot_product(w,y)/sw-dot_product(w,z)/sw                           620
      return                                                                621
      end                                                                   622
      subroutine andarm4(n,y,z,w,dst,sw)                                    623
      parameter(maxclass2=10000,nmin=100)                                   624
      real y(n),z(n),w(n),out(maxclass2)                                    625
      real, dimension (:,:), allocatable :: costs                               
      if(n .ge. nmin)goto 11501                                             628
      dst=0.0                                                               628
      sw=sum(w)                                                             628
      return                                                                628
11501 continue                                                              629
      call classin(2,idum,dum,nclass,out)                                   630
      allocate(costs(1:nclass,1:nclass),stat=jerr);                             
      call reorg(2,nclass,out,costs)                                        634
      dst=0.0                                                               635
11510 do 11511 i=1,n                                                        635
      ky=y(i)+0.1                                                           635
      kz=z(i)+0.1                                                           636
      dst=dst+w(i)*costs(ky,kz)                                             637
11511 continue                                                              638
11512 continue                                                              638
      sw=sum(w)                                                             638
      dst=dst/sw                                                            639
      return                                                                640
      end                                                                   641
      subroutine andarm5(n,y,z,w,dst,sw)                                    642
      parameter(nmin=50)                                                    643
      real y(n),z(n),w(n)                                                   644
      data qntl /0.5/                                                       645
      if(n .ge. nmin)goto 11531                                             645
      dst=0.0                                                               645
      sw=sum(w)                                                             645
      return                                                                645
11531 continue                                                              646
      sw=sum(w)                                                             646
      up=0.0                                                                647
11540 do 11541 i=1,n                                                        647
      if(y(i).le.z(i)) up=up+w(i)                                           647
11541 continue                                                              648
11542 continue                                                              648
      dst=abs(up/sw-qntl)                                                   649
      return                                                                650
      entry set_quant(arg)                                                  651
!DEC$ ATTRIBUTES DLLEXPORT :: set_quant                                         
      qntl=arg                                                              654
      return                                                                655
      end                                                                   656
      subroutine andarm6(n,y,y2,z,w,dst,sw)                                 657
      parameter(nit=100,xmiss=9.0e35,thr=1.0e-2,nmin=100)                   658
      real y(n),y2(n),z(n),w(n),yy(n,2)                                     659
      if(n .ge. nmin)goto 11561                                             659
      dst=0.0                                                               659
      sw=sum(w)                                                             659
      return                                                                659
11561 continue                                                              660
      yy(:,1)=y                                                             660
      yy(:,2)=y2                                                            661
      call cendst(n,yy,z,w,nit,thr,xmiss,dst,sw)                            662
      sw=sum(w)                                                             663
      return                                                                664
      end                                                                   665
      subroutine andarm15(n,y,y2,z,w,dst,sw)                                666
      parameter(nit=100,xmiss=9.0e35,thr=1.0e-2,nmin=100)                   667
      real y(n),y2(n),z(n),w(n),yy(n,2)                                     668
      if(n .ge. nmin)goto 11581                                             668
      dst=0.0                                                               668
      sw=sum(w)                                                             668
      return                                                                668
11581 continue                                                              669
      yy(:,1)=y                                                             669
      yy(:,2)=y2                                                            670
      call cendst1(n,yy,z,w,nit,thr,xmiss,dst,sw)                           671
      sw=sum(w)                                                             672
      return                                                                673
      end                                                                   674
      subroutine andarm10(n,y,z,w,dst,sw)                                   675
      parameter(eps=1.0e-5,nmin=100)                                        676
      real y(n),z(n),w(n)                                                   676
      integer m(n)                                                          677
      sw=sum(w)                                                             677
      if(n .ge. nmin)goto 11601                                             677
      dst=0.0                                                               677
      return                                                                677
11601 continue                                                              678
      sw1=0.0                                                               678
      sw2=sw1                                                               679
11610 do 11611 i=1,n                                                        679
      m(i)=i                                                                680
      if(z(i) .ge. 0.0)goto 11631                                           680
      sw1=sw1+w(i)                                                          680
      goto 11641                                                            680
11631 continue                                                              680
      sw2=sw2+w(i)                                                          680
11641 continue                                                              681
11621 continue                                                              681
11611 continue                                                              682
11612 continue                                                              682
      call psort8(y,m,1,n)                                                  683
      s1=0.0                                                                683
      s2=s1                                                                 683
      s=s2                                                                  683
      dst=s                                                                 684
11650 do 11651 i=1,n                                                        684
      k=m(i)                                                                684
      s=s+w(k)                                                              685
      if(z(k) .ge. 0)goto 11671                                             685
      s1=s1+w(k)/sw1                                                        685
      goto 11681                                                            686
11671 continue                                                              686
      s2=s2+w(k)/sw2                                                        686
11681 continue                                                              687
11661 continue                                                              687
      pw=s*(sw-s)                                                           687
      pw=max(eps,pw)                                                        688
      dst=dst+abs(s1-s2)/sqrt(pw)                                           689
11651 continue                                                              690
11652 continue                                                              690
      return                                                                691
      end                                                                   692
      subroutine stput (iseed)                                                  
      real x(*)                                                                 
      data i /987654321/                                                        
      i=iseed                                                                   
      return                                                                    
      entry rget (x,n)                                                          
      do 1 j=1,n                                                                
      i=mod(i*16807.0,2147483647.0)                                             
      u=i                                                                       
      u=u*.465661287d-9                                                         
    1 x(j)=u                                                                    
      return                                                                    
      entry stget (irg)                                                         
      irg=i                                                                     
      return                                                                    
      end                                                                       
      subroutine cendst(n,y,z,w,nit,thr,xmiss,dst,sw)                       711
      parameter(eps=0.1)                                                    712
      real y(n,2),z(n),w(n),b(2*n+1),q(3*n),cdf(3*n),r(n)                   713
      integer iq(3*n),mm(3*n),mz(n)                                         714
      data nsamp /500/                                                      715
      n2=2*n                                                                715
      n3=3*n                                                                715
      sw=sum(w)                                                             716
11690 do 11691 i=1,n                                                        716
      mz(i)=i                                                               716
11691 continue                                                              716
11692 continue                                                              716
      call psort8(z,mz,1,n)                                                 717
      nq=0.25*n                                                             717
      teps=(z(mz(n-nq))-z(mz(nq)))*eps                                      718
11700 do 11701 i=1,n                                                        718
      if(y(i,2)-y(i,1).ge.teps)goto 11701                                   719
      y(i,1)=y(i,1)-teps                                                    719
      y(i,2)=y(i,2)+teps                                                    720
11701 continue                                                              721
11702 continue                                                              721
11710 do 11711 i=1,n                                                        721
      b(i)=y(i,1)                                                           721
      b(i+n)=y(i,2)                                                         721
11711 continue                                                              722
11712 continue                                                              722
      m=0                                                                   723
11720 do 11721 i=1,n2                                                       723
      if(b(i).le.-xmiss)goto 11721                                          723
      if(b(i).ge.xmiss)goto 11721                                           724
      m=m+1                                                                 724
      b(m)=b(i)                                                             725
11721 continue                                                              726
11722 continue                                                              726
      call unique(m,b,nu)                                                   727
      if(nu .le. nsamp)goto 11741                                           727
      call rget(r,nsamp)                                                    728
11750 do 11751 i=1,nsamp                                                    728
      r(i)=b(int(nu*r(i))+1)                                                728
11751 continue                                                              729
11752 continue                                                              729
      nu=nsamp                                                              729
      b(1:nu)=r(1:nu)                                                       730
      call sort(b,nu)                                                       731
11741 continue                                                              732
      m=nu+1                                                                732
      b(m)=xmiss                                                            733
      call fintcdf1(n,y,m,b,w,nit,thr/m,cdf,jt,err)                         734
      m=m-1                                                                 735
11760 do 11761 i=1,m                                                        735
      q(i)=b(i)                                                             735
      iq(i)=0                                                               735
11761 continue                                                              736
11762 continue                                                              736
11770 do 11771 i=1,n                                                        736
      mz(i)=i                                                               736
11771 continue                                                              736
11772 continue                                                              736
      call psort8(z,mz,1,n)                                                 737
      mpn=m+n                                                               737
      k=0                                                                   738
11780 do 11781 i=m+1,mpn                                                    738
      k=k+1                                                                 738
      q(i)=z(mz(k))                                                         738
      iq(i)=1                                                               738
11781 continue                                                              739
11782 continue                                                              739
11790 do 11791 i=1,mpn                                                      739
      mm(i)=i                                                               739
11791 continue                                                              739
11792 continue                                                              739
      call psort8(q,mm,1,mpn)                                               740
      ycdf=0.0                                                              740
      zcdf=ycdf                                                             740
      dst=zcdf                                                              740
      spw=dst                                                               740
      ny=0                                                                  740
      nz=ny                                                                 741
11800 do 11801 i=1,mpn                                                      741
      k=mm(i)                                                               742
      if(iq(k) .ne. 0)goto 11821                                            742
      ny=ny+1                                                               742
      ycdf=cdf(ny)                                                          743
      pw=float(i)*float(mpn-i)/float(mpn)**2                                744
      pw=max(eps,pw)                                                        744
      pw=1.0/sqrt(pw)                                                       744
      spw=spw+pw                                                            745
      dst=dst+pw*abs(ycdf-zcdf)                                             746
      goto 11831                                                            747
11821 continue                                                              747
      nz=nz+1                                                               747
      zcdf=zcdf+w(nz)/sw                                                    747
11831 continue                                                              748
11811 continue                                                              748
11801 continue                                                              749
11802 continue                                                              749
      dst=dst/spw                                                           750
      return                                                                751
      entry set_samp(irg)                                                   752
!DEC$ ATTRIBUTES DLLEXPORT :: set_samp                                          
      nsamp=irg                                                             755
      call set_samp1(nsamp)                                                 755
      return                                                                756
      end                                                                   757
      subroutine cendst1(n,y,z,w,nit,thr,xmiss,dst,sw)                      758
      real y(n,2),z(n),w(n),b(2*n+1),cdf1(3*n),cdf2(3*n),r(n)               759
      real y1(n,2),y2(n,2),w1(n),w2(n)                                      760
      data nsamp /500/                                                      761
      n1=0.0                                                                761
      n2=n1                                                                 762
11840 do 11841 i=1,n                                                        763
      if(y(i,1).le.-xmiss)goto 11841                                        763
      if(y(i,2).ge.xmiss)goto 11841                                         764
      if(y(i,2)-y(i,1).ge.teps)goto 11841                                   765
      y(i,1)=y(i,1)-teps                                                    765
      y(i,2)=y(i,2)+teps                                                    766
11841 continue                                                              767
11842 continue                                                              767
11850 do 11851 i=1,n                                                        768
      if(z(i) .ge. 0.0)goto 11871                                           768
      n1=n1+1                                                               768
      y1(n1,:)=y(i,:)                                                       768
      w1(n1)=w(i)                                                           768
      goto 11881                                                            769
11871 continue                                                              769
      n2=n2+1                                                               769
      y2(n2,:)=y(i,:)                                                       769
      w2(n2)=w(i)                                                           769
11881 continue                                                              770
11861 continue                                                              770
11851 continue                                                              771
11852 continue                                                              771
11890 do 11891 i=1,n                                                        771
      b(i)=y(i,1)                                                           771
      b(i+n)=y(i,2)                                                         771
11891 continue                                                              772
11892 continue                                                              772
      m=0                                                                   773
11900 do 11901 i=1,n2                                                       773
      if(b(i).le.-xmiss)goto 11901                                          773
      if(b(i).ge.xmiss)goto 11901                                           774
      m=m+1                                                                 774
      b(m)=b(i)                                                             775
11901 continue                                                              776
11902 continue                                                              776
      call unique(m,b,nu)                                                   777
      if(nu .le. nsamp)goto 11921                                           777
      call rget(r,nsamp)                                                    778
11930 do 11931 i=1,nsamp                                                    778
      r(i)=b(int(nu*r(i))+1)                                                778
11931 continue                                                              779
11932 continue                                                              779
      nu=nsamp                                                              779
      b(1:nu)=r(1:nu)                                                       780
      call sort(b,nu)                                                       781
11921 continue                                                              782
      m=nu+1                                                                782
      b(m)=xmiss                                                            783
      call getcdf1(n1,y1,w1,nit,thr,xmiss,nsamp,m,b,cdf1,sw1)               784
      call getcdf1(n2,y2,w2,nit,thr,xmiss,nsamp,m,b,cdf2,sw2)               785
      call diffcdf(m,cdf1,cdf2,dst)                                         786
      return                                                                787
      entry set_samp1(irg)                                                  788
      nsamp=irg                                                             788
      return                                                                789
      end                                                                   790
      subroutine getcdf1(n,y,w,nit,thr,xmiss,nsamp,m,b,cdf,sw)              791
      parameter(teps=0.01)                                                  792
      real y(n,2),z(n),w(n),b(2*n+1),q(3*n),cdf(3*n),r(n)                   793
      integer iq(3*n),mm(3*n),mz(n)                                         794
      n2=2*n                                                                794
      sw=sum(w)                                                             795
      call fintcdf1(n,y,m,b,w,nit,thr/m,cdf,jt,err)                         796
      m=m-1                                                                 797
      return                                                                798
      end                                                                   799
      subroutine diffcdf(m,cdf1,cdf2,dst)                                   800
      real cdf1(m),cdf2(m)                                                  801
      f12=sqrt(float(m))                                                    801
      dst=0.0                                                               802
11940 do 11941 i=1,m                                                        803
      dst=dst+abs(cdf1(i)-cdf2(i))/sqrt(float(i)*float(m-i+1))              804
11941 continue                                                              805
11942 continue                                                              805
      dst=f12*dst/m                                                         806
      return                                                                807
      end                                                                   808
      subroutine unique(n,y,nu)                                             809
      real y(n),yu(n)                                                       809
      integer m(n)                                                          810
11950 do 11951 i=1,n                                                        810
      m(i)=i                                                                810
11951 continue                                                              810
11952 continue                                                              810
      call psort8(y,m,1,n)                                                  811
      nu=1                                                                  811
      yu(1)=y(m(1))                                                         812
11960 do 11961 i=2,n                                                        812
      if(y(m(i-1)).ge.y(m(i)))goto 11961                                    813
      nu=nu+1                                                               813
      yu(nu)=y(m(i))                                                        814
11961 continue                                                              815
11962 continue                                                              815
      y(1:nu)=yu(1:nu)                                                      816
      return                                                                817
      end                                                                   818
      subroutine fintcdf1(n,y,m,b,w1,nit,thr,cdf,jt,err)                    819
!DEC$ ATTRIBUTES DLLEXPORT :: fintcdf1                                          
      real y(n,2),b(m),w(n),w1(n),p(m),pij(n,m),ps(m),cdf(m)                823
      integer, dimension (:), allocatable :: ic,jc,kc1,kc2,lc                   
      real, dimension (:), allocatable :: smo                                   
      call set_vrb(ivrb,2)                                                  828
      w=w1/sum(w1)                                                          828
      p=1.0/m                                                               829
      nt=0                                                                  830
11970 do 11971 i=1,n                                                        830
11980 do 11981 k=1,m                                                        830
      if(y(i,1).ge.b(k))goto 11981                                          831
      if(y(i,2).lt.b(k))goto 11981                                          831
      nt=nt+1                                                               832
11981 continue                                                              832
11982 continue                                                              832
11971 continue                                                              833
11972 continue                                                              833
      allocate(ic(1:(n+1)),stat=jerr)                                       834
      allocate(jc(1:nt),stat=ierr)                                          834
      jerr=jerr+ierr                                                        835
      allocate(kc1(1:m),stat=ierr)                                          835
      jerr=jerr+ierr                                                        836
      allocate(kc2(1:m),stat=ierr)                                          836
      jerr=jerr+ierr                                                        837
      allocate(smo(1:m),stat=ierr)                                          837
      jerr=jerr+ierr                                                        838
      allocate(lc(1:nt),stat=ierr)                                          838
      jerr=jerr+ierr                                                        839
      if(jerr .eq. 0)goto 12001                                             839
      err=8888.0                                                            839
      return                                                                839
12001 continue                                                              840
      nt=0                                                                  840
      ic(1)=1                                                               841
12010 do 12011 i=1,n                                                        842
12020 do 12021 k=1,m                                                        842
      if(y(i,1).ge.b(k))goto 12021                                          842
      if(y(i,2).lt.b(k))goto 12021                                          843
      nt=nt+1                                                               843
      jc(nt)=k                                                              844
12021 continue                                                              845
12022 continue                                                              845
      ic(i+1)=nt+1                                                          846
12011 continue                                                              847
12012 continue                                                              847
      nt=0                                                                  848
12030 do 12031 j=1,m                                                        848
      kc1(j)=nt+1                                                           849
12040 do 12041 i=1,n                                                        850
      if(y(i,1).ge.b(j))goto 12041                                          850
      if(y(i,2).lt.b(j))goto 12041                                          851
      nt=nt+1                                                               851
      lc(nt)=i                                                              852
12041 continue                                                              853
12042 continue                                                              853
      kc2(j)=nt                                                             854
12031 continue                                                              855
12032 continue                                                              855
      if(ivrb.gt.0) write(6,'(''CDF iterations'',$)')                       856
12050 do 12051 it=1,nit                                                     856
      jt=it                                                                 856
      ps=p                                                                  857
12060 do 12061 j=1,m                                                        857
      pij(:,j)=0.0                                                          858
12070 do 12071 ii=kc1(j),kc2(j)                                             858
      i=lc(ii)                                                              859
      s=sum(p(jc(ic(i):(ic(i+1)-1))))                                       860
      if(s .gt. 0.0)goto 12091                                              860
      err=-7777.0                                                           860
      return                                                                860
12091 continue                                                              861
      pij(i,j)=w(i)*p(j)/s                                                  863
12071 continue                                                              863
12072 continue                                                              863
      p(j)=sum(pij(:,j))                                                    864
12061 continue                                                              865
12062 continue                                                              865
      if(m .le. 100)goto 12111                                              866
      smo(1)=(2.0*p(1)+p(2))/3.0                                            866
      smo(m)=(2.0*p(m)+p(m-1))/3.0                                          867
      smo(2)=0.25*(p(1)+2.0*p(2)+p(3))                                      867
      smo(m-1)=0.25*(p(m)+2.0*p(m-1)+p(m-2))                                868
12120 do 12121 j=3,m-2                                                      868
      smo(j)=(p(j-2)+2.0*p(j-1)+3.0*p(j)+2.0*p(j+1)+p(j+2))/9.0             868
12121 continue                                                              869
12122 continue                                                              869
      p=smo                                                                 870
12111 continue                                                              871
      err=sum(abs(p-ps))/m                                                  872
      if(kbad(err) .le. 0)goto 12141                                        872
      err=7777.0                                                            872
      return                                                                873
12141 continue                                                              873
      if(err.lt.thr)goto 12052                                              874
      if(ivrb.gt.0) write(6,'(''.'',$)')                                    875
12051 continue                                                              876
12052 continue                                                              876
      cdf(1)=p(1)                                                           876
12150 do 12151 j=2,m                                                        876
      cdf(j)=cdf(j-1)+p(j)                                                  876
12151 continue                                                              877
12152 continue                                                              877
      if(ivrb .le. 0)goto 12171                                             877
      write(6,12180)err                                                     877
12180 format (g10.2)                                                        877
12171 continue                                                              878
      return                                                                879
      end                                                                   880
      subroutine set_vrb(irg,jrg)                                           881
!DEC$ ATTRIBUTES DLLEXPORT :: set_vrb                                           
      save ivrb                                                             885
      if(jrg .ne. 1)goto 12201                                              885
      ivrb=irg                                                              885
      return                                                                885
12201 continue                                                              886
      irg=ivrb                                                              887
      return                                                                888
      end                                                                   889
      subroutine cdfpoints1(m,x,n,y,w,cdf)                                  890
!DEC$ ATTRIBUTES DLLEXPORT :: cdfpoints1                                        
      real x(m),y(n),w(n),cdf(m)                                            894
      i=1                                                                   894
      j=0                                                                   894
      sw=0.0                                                                895
12210 continue                                                              895
12211 continue                                                              895
      j=j+1                                                                 895
      if(j.gt.m) go to 12220                                                896
12230 continue                                                              896
12231 if(y(i).gt.x(j))goto 12232                                            896
      sw=sw+w(i)                                                            896
      i=i+1                                                                 896
      if(i.gt.n)goto 12232                                                  896
      goto 12231                                                            897
12232 continue                                                              897
      if(i .le. n)goto 12251                                                898
12260 do 12261 k=j,m                                                        898
      cdf(k)=sw                                                             898
12261 continue                                                              898
12262 continue                                                              898
      go to 12270                                                           899
12251 continue                                                              900
      cdf(j)=sw                                                             901
      goto 12211                                                            902
12212 continue                                                              902
12270 continue                                                              903
12220 continue                                                              904
      cdf=cdf/sum(w)                                                        905
      return                                                                906
      end                                                                   907
      subroutine sort(x,n)                                                  908
      real x(n),z(n)                                                        908
      integer m(n)                                                          909
12280 do 12281 i=1,n                                                        909
      m(i)=i                                                                909
12281 continue                                                              909
12282 continue                                                              909
      z=x                                                                   910
      call psort8(z,m,1,n)                                                  911
12290 do 12291 i=1,n                                                        911
      x(i)=z(m(i))                                                          911
12291 continue                                                              912
12292 continue                                                              912
      return                                                                913
      end                                                                   914
      function kbad(u)                                                      915
      kbad=0                                                                916
      if(isnan(u).or.abs(u).ge.abs(huge(u))) kbad=1                         917
      return                                                                918
      end                                                                   919
      subroutine classin(ient,nclasssv,costssv,nout,out)                    920
!DEC$ ATTRIBUTES DLLEXPORT :: classin                                           
      real costssv(nclasssv,nclasssv),out(1)                                924
      real, dimension (:), allocatable :: costs                                 
      save costs,nclass                                                     928
      nq=nclasssv*nclasssv                                                  929
      allocate(costs(1:nq),stat=jerr)                                       930
      if(ient .ne. 1)goto 12311                                             930
      nclass=nclasssv                                                       931
      call reorg(1,nclass,costs,costssv)                                    932
      nout=1                                                                932
      out=1.0                                                               933
      goto 12321                                                            934
12311 continue                                                              934
      nout=nclass                                                           934
      call reorg(2,nclass,costs,out)                                        934
12321 continue                                                              935
12301 continue                                                              935
      return                                                                936
      end                                                                   937
      subroutine reorg(ient,n,a,b)                                          938
      real a(n*n),b(n,n)                                                    939
      i=0                                                                   940
      if(ient .ne. 2)goto 12341                                             940
12350 do 12351 k=1,n                                                        940
12360 do 12361 j=1,n                                                        940
      i=i+1                                                                 940
      b(j,k)=a(i)                                                           940
12361 continue                                                              940
12362 continue                                                              940
12351 continue                                                              940
12352 continue                                                              940
      goto 12371                                                            941
12341 continue                                                              941
12380 do 12381 k=1,n                                                        941
12390 do 12391 j=1,n                                                        941
      i=i+1                                                                 941
      a(i)=b(j,k)                                                           941
12391 continue                                                              941
12392 continue                                                              941
12381 continue                                                              941
12382 continue                                                              941
12371 continue                                                              942
12331 continue                                                              942
      return                                                                943
      end                                                                   944
      subroutine crinode (itr,rtr,mxnodes,node,nodes,cri,wt)                945
!DEC$ ATTRIBUTES DLLEXPORT :: crinode                                           
      integer itr(6,*),nodes(mxnodes),m(mxnodes),ic(mxnodes)                949
      real rtr(4,*),cri(mxnodes),wt(mxnodes),sc(mxnodes,2)                  950
      node=0                                                                950
      k=itr(2,1)                                                            951
12400 continue                                                              951
12401 continue                                                              951
      if(itr(4,k) .lt. 0)goto 12421                                         951
      k=itr(2,k)                                                            951
      goto 12401                                                            951
12421 continue                                                              951
      i1=itr(5,k)                                                           951
      i2=itr(6,k)                                                           952
      node=node+1                                                           952
      if(node.gt.mxnodes) return                                            953
      nodes(node)=k                                                         953
      cri(node)=rtr(3,k)                                                    953
      wt(node)=rtr(4,k)                                                     954
12430 continue                                                              954
12431 if(k.eq.itr(2,iabs(itr(4,k))))goto 12432                              954
      k=iabs(itr(4,k))                                                      954
      if(k.eq.1)goto 12432                                                  954
      goto 12431                                                            955
12432 continue                                                              955
      if(k.eq.1)goto 12402                                                  955
      k=itr(3,iabs(itr(4,k)))                                               956
      goto 12401                                                            957
12402 continue                                                              957
12440 do 12441 k=1,node                                                     957
      m(k)=k                                                                957
12441 continue                                                              957
12442 continue                                                              957
      call psort8(-cri,m,1,node)                                            958
12450 do 12451 i=1,node                                                     958
      ic(i)=nodes(m(i))                                                     958
      sc(i,1)=cri(m(i))                                                     958
      sc(i,2)=wt(m(i))                                                      958
12451 continue                                                              959
12452 continue                                                              959
12460 do 12461 i=1,node                                                     959
      nodes(i)=ic(i)                                                        959
      cri(i)=sc(i,1)                                                        959
      wt(i)=sc(i,2)                                                         959
12461 continue                                                              960
12462 continue                                                              960
      return                                                                961
      end                                                                   962
      subroutine prune1 (itr,rtr,nodes,thr,itro,rtro)                       963
!DEC$ ATTRIBUTES DLLEXPORT :: prune1                                            
      integer itr(6,nodes)                                                  966
      real rtr(4,nodes)                                                     967
      integer itro(6,nodes)                                                 967
      real rtro(4,nodes)                                                    968
      call prune(itr,rtr,nodes,thr)                                         969
      itro=itr                                                              969
      rtro=rtr                                                              970
      return                                                                971
      end                                                                   972
      subroutine prune (itr,rtr,nodes,thr)                                  973
      integer itr(6,nodes)                                                  973
      real rtr(4,nodes)                                                     974
12470 continue                                                              974
12471 continue                                                              974
      nch=0                                                                 975
12480 do 12481 k=1,nodes                                                    975
      if(itr(4,k).le.0)goto 12481                                           976
      nl=itr(2,k)                                                           976
      nr=itr(3,k)                                                           977
      if(itr(4,nl).ge.0)goto 12481                                          978
      if(itr(4,nr).ge.0)goto 12481                                          979
      if(max(rtr(3,nl),rtr(3,nr)).gt.rtr(3,k)+thr)goto 12481                980
      itr(4,k)=-itr(4,k)                                                    980
      nch=nch+1                                                             981
12481 continue                                                              982
12482 continue                                                              982
      if(nch.eq.0)goto 12472                                                982
      goto 12471                                                            983
12472 continue                                                              983
      return                                                                984
      end                                                                   985
      subroutine getnode (x,itre,rtre,cat,node)                             986
      integer itre(6,*)                                                     986
      real x(*),rtre(4,*),cat(*)                                            987
      k=1                                                                   988
12490 continue                                                              988
12491 if(itre(4,k).lt.0)goto 12492                                          988
      if(itre(1,k) .le. 0)goto 12511                                        990
      if(x(itre(1,k)) .ge. rtre(1,k))goto 12531                             990
      k=itre(2,k)                                                           990
      goto 12541                                                            990
12531 continue                                                              990
      k=itre(3,k)                                                           990
12541 continue                                                              991
12521 continue                                                              991
      goto 12551                                                            992
12511 continue                                                              992
      j=-itre(1,k)                                                          992
      kp=rtre(1,k)+0.1                                                      992
      in=0                                                                  992
      np=abs(cat(kp))+0.1                                                   993
12560 do 12561 i=1,np                                                       993
      if(x(j).ne.cat(kp+i))goto 12561                                       993
      in=1                                                                  993
      goto 12562                                                            993
12561 continue                                                              994
12562 continue                                                              994
      if(cat(kp) .le. 0.0)goto 12581                                        994
      if(in .ne. 0)goto 12601                                               994
      k=itre(2,k)                                                           994
      goto 12611                                                            994
12601 continue                                                              994
      k=itre(3,k)                                                           994
12611 continue                                                              994
12591 continue                                                              994
      goto 12621                                                            995
12581 continue                                                              995
      if(in .eq. 0)goto 12641                                               995
      k=itre(2,k)                                                           995
      goto 12651                                                            995
12641 continue                                                              995
      k=itre(3,k)                                                           995
12651 continue                                                              995
12631 continue                                                              995
12621 continue                                                              996
12571 continue                                                              996
12551 continue                                                              997
12501 continue                                                              997
      goto 12491                                                            998
12492 continue                                                              998
      node=k                                                                999
      return                                                               1000
      end                                                                  1001
      subroutine getnodes1 (no,ni,x,itre,rtre,cat,nodes)                   1002
!DEC$ ATTRIBUTES DLLEXPORT :: getnodes1                                         
      integer itre(6,*),nodes(no)                                          1005
      real x(no,ni),rtre(4,*),cat(*)                                       1006
12660 do 12661 i=1,no                                                      1006
      call getnode (x(i,:),itre,rtre,cat,node)                             1006
      nodes(i)=node                                                        1006
12661 continue                                                             1007
12662 continue                                                             1007
      return                                                               1008
      end                                                                  1009
      subroutine getlims(node,ni,itr,rtr,cat,nvar,jvar,vlims,jerr)         1010
!DEC$ ATTRIBUTES DLLEXPORT :: getlims                                           
      integer itr(6,*),jvar(2,*)                                           1013
      real rtr(4,*),cat(*),vlims(*)                                        1014
      jerr=0                                                               1014
      if(itr(4,node) .lt. 0)goto 12681                                     1014
      jerr=1                                                               1014
      return                                                               1014
12681 continue                                                             1015
      nvar=0                                                               1016
      k=node                                                               1017
12690 continue                                                             1017
12691 continue                                                             1017
      nvar=nvar+1                                                          1017
      kpp=abs(itr(4,k))                                                    1018
      if(itr(1,kpp) .le. 0)goto 12711                                      1018
      jvar(2,nvar)=0                                                       1019
      if(itr(2,kpp) .ne. k)goto 12731                                      1019
      jvar(1,nvar)=-itr(1,kpp)                                             1019
      goto 12741                                                           1020
12731 continue                                                             1020
      jvar(1,nvar)=itr(1,kpp)                                              1020
12741 continue                                                             1021
12721 continue                                                             1021
      vlims(nvar)=rtr(1,kpp)                                               1022
      goto 12751                                                           1023
12711 continue                                                             1024
      if(k .ne. itr(2,kpp))goto 12771                                      1024
      sgn=-1.0                                                             1024
      goto 12781                                                           1024
12771 continue                                                             1024
      sgn=1.0                                                              1024
12781 continue                                                             1025
12761 continue                                                             1025
      jvar(1,nvar)=-itr(1,kpp)                                             1025
      kp=rtr(1,kpp)+0.1                                                    1025
      jvar(2,nvar)=kp                                                      1026
      vlims(nvar)=sgn*abs(cat(kp))                                         1027
12751 continue                                                             1028
12701 continue                                                             1028
      k=kpp                                                                1029
      if(kpp.eq.1)goto 12692                                               1029
      goto 12691                                                           1030
12692 continue                                                             1030
      return                                                               1031
      end                                                                  1032
      subroutine trans(ny,y,wy,nz,z,wz,nt,t)                               1033
!DEC$ ATTRIBUTES DLLEXPORT :: trans                                             
      real y(ny),wy(ny),z(nz),wz(nz),t(nt+2,2),u(max(ny,nz))               1037
      real p(nt)                                                           1038
      integer m(max(ny,nz))                                                1039
12790 do 12791 i=1,ny                                                      1039
      m(i)=i                                                               1039
      u(i)=y(i)                                                            1039
12791 continue                                                             1039
12792 continue                                                             1039
      call psort8(u,m,1,ny)                                                1040
12800 do 12801 i=1,ny                                                      1040
      y(i)=u(m(i))                                                         1040
12801 continue                                                             1040
12802 continue                                                             1040
      u=wy                                                                 1040
12810 do 12811 i=1,ny                                                      1040
      wy(i)=u(m(i))                                                        1040
12811 continue                                                             1041
12812 continue                                                             1041
12820 do 12821 i=1,nz                                                      1041
      m(i)=i                                                               1041
      u(i)=z(i)                                                            1041
12821 continue                                                             1041
12822 continue                                                             1041
      call psort8(u,m,1,nz)                                                1042
12830 do 12831 i=1,nz                                                      1042
      z(i)=u(m(i))                                                         1042
12831 continue                                                             1042
12832 continue                                                             1042
      u=wz                                                                 1042
12840 do 12841 i=1,nz                                                      1042
      wz(i)=u(m(i))                                                        1042
12841 continue                                                             1043
12842 continue                                                             1043
12850 do 12851 i=1,nt                                                      1043
      p(i)=(i-0.5)/float(nt)                                               1043
12851 continue                                                             1044
12852 continue                                                             1044
      call untie(ny,y,u)                                                   1045
      call qntl(ny,u,wy,nt,p,t(:,1))                                       1046
      call untie(nz,z,u)                                                   1047
      call qntl(nz,u,wz,nt,p,t(:,2))                                       1048
      return                                                               1049
      end                                                                  1050
      subroutine qntl(n,y,w,nq,p,q)                                        1051
      real y(n),w(n),p(nq),q(nq+2)                                         1052
      sw=sum(w)                                                            1052
      k=1                                                                  1052
      ff=w(1)                                                              1052
      q(1)=y(1)                                                            1052
      q(nq+2)=y(n)                                                         1053
12860 do 12861 i=2,n                                                       1053
      ff=ff+w(i)                                                           1053
      pp=ff/sw                                                             1053
      if(pp.lt.p(k))goto 12861                                             1054
      k=k+1                                                                1054
      q(k)=0.5*(y(i)+y(i-1))                                               1055
      if(k.ge.nq)goto 12862                                                1056
12861 continue                                                             1057
12862 continue                                                             1057
      q(nq+1)=0.5*(q(nq+2)+q(nq))                                          1058
      return                                                               1059
      end                                                                  1060
      subroutine untie(n,y,u)                                              1061
!DEC$ ATTRIBUTES DLLEXPORT :: untie                                             
      real y(n),u(n)                                                       1065
      i=1                                                                  1065
      k=0                                                                  1066
12870 continue                                                             1066
12871 if(i.ge.n)goto 12872                                                 1066
      if(y(i+1) .le. y(i))goto 12891                                       1067
      k=k+1                                                                1067
      u(k)=y(i)                                                            1067
      i=i+1                                                                1067
      goto 12871                                                           1067
12891 continue                                                             1068
      i1=i                                                                 1068
12900 continue                                                             1068
12901 if(y(i+1).gt.y(i))goto 12902                                         1068
      i=i+1                                                                1068
      if(i.ge.n)goto 12902                                                 1068
      goto 12901                                                           1068
12902 continue                                                             1068
      i2=i                                                                 1069
      if(i1 .gt. 1)goto 12921                                              1069
      a=y(i1+1)                                                            1069
      b=y(i2+1)                                                            1069
      u(1)=y(1)                                                            1069
      k=1                                                                  1070
12930 do 12931 j=i1+1,i2                                                   1070
      k=k+1                                                                1070
      u(k)=a+(b-a)*(j-i1)/(i2-i1+1)                                        1070
12931 continue                                                             1071
12932 continue                                                             1071
      i=i2+1                                                               1072
      goto 12911                                                           1073
12921 if(i2 .lt. n)goto 12941                                              1074
      a=y(i1-1)                                                            1074
      b=(y(i2)-a)/(i2-i1+1)                                                1075
12950 do 12951 j=i1,i2                                                     1075
      k=k+1                                                                1075
      u(k)=a+b*(j-i1+1)                                                    1075
12951 continue                                                             1076
12952 continue                                                             1076
      goto 12961                                                           1077
12941 continue                                                             1077
      a=y(i1-1)                                                            1077
      b=y(i2)                                                              1078
12970 do 12971 j=i1,i2                                                     1078
      k=k+1                                                                1078
      u(k)=a+(b-a)*(j-i1+1)/(i2-i1+1)                                      1078
12971 continue                                                             1079
12972 continue                                                             1079
      i=i+1                                                                1080
12961 continue                                                             1081
12911 continue                                                             1081
      goto 12871                                                           1082
12872 continue                                                             1082
      if(k .ge. n)goto 12991                                               1082
      k=k+1                                                                1082
      u(k)=y(n)                                                            1082
12991 continue                                                             1083
      return                                                               1084
      end                                                                  1085
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
