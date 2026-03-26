        subroutine subaeroretwat(iband1,iband3,xts,xtv,xfi,pres,uoz,uwv,erelc,troatm,
     c       tpres,aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,iv,raot,residual,nit,snext,ierr,
     s       ifast,tbtgo,tbroatm,tbttatmg,tbsatm,iaots,eps)
        IMPLICIT NONE
c driver for the atmopsheric correction program
c read up the look up table and perform correction
c look up table variable
c   Aerosol model (string) and length of the string
c local variable
        real xts,xtv,xfi,raot550nm,uoz,uwv,pres
        real residual
        integer ib,iv
        real roslamb,ros1,ros3,raot1,raot2
        integer flagn
       
       real tauray(16)
       real oztransa(16)
       real wvtransa(16),wvtransb(16),wvtransc(16)
       real ogtransa0(16),ogtransa1(16)
       real ogtransb0(16),ogtransb1(16)
       real ogtransc0(16),ogtransc1(16)      
       real rolutt(16,7,22,8000),nbfi(22,20),tsmax(22,20),tsmin(22,20)
       real ttv(22,20),nbfic(22,20)
       real tts(22)
       integer indts(22)
       real transt(16,7,22,22)
       real sphalbt(16,7,22),normext(16,7,22)
       real xtsmin,xtsstep,xtvmin,xtvstep
       real aot550nm(22)
       real next,snext
       real tpres(7)
       integer retval
       integer iaot
        real tgo,roatm,ttatmg,satm,xrorayp
        real pratio
        real troatm(16)
        real residual1,residual2
        real raot
        real erelc(16)
        integer iband1,iband3,nit
        character*80 err_msg
        integer nbval,testth,ierr
        real tth(16)
        integer iverbose
        real xa,xb,xc,xd,xe,xf,coefa,coefb,raotmin,residualm,raotsaved
        integer ifast,iaots
        real tbtgo(16,22),tbroatm(16,22),tbttatmg(16,22),tbsatm(16,22)
        integer iaot1,iaot2
        real eps
           
c        eps=1.
        iverbose=iv
c correct band 3 and band 1 with increasing AOT (using pre till ratio is equal to erelc(3)
        iaot=iaots
c        iaots=1
        pratio=erelc(iband3)/erelc(iband1)
c pratio is the targeted ratio between the surface reflectance in band 4 (ros4) and 2 (ros2)
        residual1=2000.
        residual2=1000.
        iaot2=1
        iaot1=1
        raot2=1.E-06
        raot1=1.E-05
        ros1=1.0
        ros3=1.0
        raot1=0.0001
        flagn=0
c ESPA these are the Landsat water bands (error) vs. the Sentinel bands 
c ORIG
c        tth(1)=1.E-03
c        tth(2)=1.E-03
c        tth(4)=1.E-03
c        tth(5)=1.E-03
c        tth(7)=1.E-04

c ESPA
         tth(1)=1.E-03
         tth(4)=1.E-03
         tth(9)=1.E-03
         tth(13)=1.E-04
        nit=0
        residual=0
        iaot=iaots
        raot550nm=aot550nm(iaot)
        ib=iband1
        testth=0
        ierr=0

        if (ifast.eq.1) then
             roslamb=troatm(ib)/tbtgo(ib,iaot)
              roslamb=roslamb-tbroatm(ib,iaot)
              roslamb=roslamb/tbttatmg(ib,iaot)
              roslamb=roslamb/(1.+tbsatm(ib,iaot)*roslamb)
        else        
        call atmcorlamb2(xts,xtv,xfi,raot550nm,ib,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       troatm(ib),roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval,eps,iverbose)
        endif
     
        if ((roslamb-tth(iband1)).lt.0.) testth=1
        ros1=roslamb
        
c        if (iverbose.eq.1) then
c        write(6,*) "water1 band1 ",ib
c        write(6,*) "      tau troatm1 ros1 eps ",raot550nm,troatm(ib),ros1,eps
c        write(6,*) "      iaot ",iaot
c        endif
        
        nbval=0
        residual=0.
        do ib=1,16
        if (erelc(ib).gt.0.) then
        
        if (ifast.eq.1) then
             roslamb=troatm(ib)/tbtgo(ib,iaot)
              roslamb=roslamb-tbroatm(ib,iaot)
              roslamb=roslamb/tbttatmg(ib,iaot)
              roslamb=roslamb/(1.+tbsatm(ib,iaot)*roslamb)
        else        
        call atmcorlamb2(xts,xtv,xfi,raot550nm,ib,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       troatm(ib),roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval,eps,iverbose)
        endif
        residual=residual+roslamb*roslamb
        nbval=nbval+1
        if ((roslamb-tth(ib)).lt.0.) testth=1
        endif
        enddo
        residual=sqrt(residual)/nbval
        
        iaot=iaot+1
        do while ((iaot.le.22).and.(residual.lt.residual1).and.(testth.ne.1))
        residual2=residual1
        iaot2=iaot1
        raot2=raot1
        residual1=residual
        raot1=raot550nm
        iaot1=iaot
        raot550nm=aot550nm(iaot)
        ib=iband1
        testth=0
        if (ifast.eq.1) then
             roslamb=troatm(ib)/tbtgo(ib,iaot)
              roslamb=roslamb-tbroatm(ib,iaot)
              roslamb=roslamb/tbttatmg(ib,iaot)
              roslamb=roslamb/(1.+tbsatm(ib,iaot)*roslamb)
        else        
        call atmcorlamb2(xts,xtv,xfi,raot550nm,ib,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       troatm(ib),roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval,eps,iverbose)
        endif
        if ((roslamb-tth(iband1)).lt.0.) testth=1
        ros1=roslamb
        nbval=0
        residual=0.
        do ib=1,16
        if ((erelc(ib).gt.0.)) then
        
        if (ifast.eq.1) then
             roslamb=troatm(ib)/tbtgo(ib,iaot)
              roslamb=roslamb-tbroatm(ib,iaot)
              roslamb=roslamb/tbttatmg(ib,iaot)
              roslamb=roslamb/(1.+tbsatm(ib,iaot)*roslamb)
        else        
        call atmcorlamb2(xts,xtv,xfi,raot550nm,ib,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       troatm(ib),roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval,eps,iverbose)
             endif
        if ((roslamb-tth(ib)).lt.0.) testth=1
        residual=residual+(roslamb*roslamb)
        nbval=nbval+1
        endif
        enddo

        residual=sqrt(residual)/nbval
        iaot=iaot+1
        enddo
        
        if ((iaot.eq.2)) then
        raot=raot550nm
        goto 999
        endif
        

        raot=raot550nm
        raotsaved=raot
        xa=(raot1*raot1)-(raot*raot)
        xd=(raot2*raot2-raot*raot)
        xb=(raot1-raot)
        xe=(raot2-raot)
        xc=residual1-residual
        xf=residual2-residual
        coefa=(xc*xe-xb*xf)/(xa*xe-xb*xd)
        coefb=(xa*xf-xc*xd)/(xa*xe-xb*xd)
        raotmin=-coefb/(2*coefa)

        if ((raotmin.lt.0.01).or.(raotmin.gt.4.0)) then
        raotmin=raot
        endif

        raot550nm=raotmin
        ib=iband1
        testth=0
        call atmcorlamb2(xts,xtv,xfi,raot550nm,ib,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       troatm(ib),roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval,eps,iverbose)
        if ((roslamb-tth(iband1)).lt.0.) testth=1
        ros1=roslamb

        nbval=0
        residualm=0.
        do ib=1,16
        if ((erelc(ib).gt.0.)) then
        call atmcorlamb2(xts,xtv,xfi,raot550nm,ib,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       troatm(ib),roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval,eps,iverbose)
        if ((roslamb-tth(ib)).lt.0.) testth=1
        residualm=residualm+roslamb*roslamb
        nbval=nbval+1
        endif
        enddo
        residualm=sqrt(residualm)/nbval
        raot=raot550nm
        
        if (residualm.gt.residual) then
            residualm=residual
            raot=raotsaved
            endif
        if (residualm.gt.residual1) then
            residualm=residual1
            raot=raot1
            endif
        if (residualm.gt.residual2) then
            residualm=residual2
            raot=raot2
            endif
        
        residual=residualm
        
        if ((iaot.eq.2)) then
        ierr=1
        iaots=1
        goto 999
        endif
        
        iaots=max((iaot2-3),1)        
        
999      if (ierr.eq.1) write(6,*) "problem to be solved"
        
        if (testth.eq.1) ierr=1
          return
        end
