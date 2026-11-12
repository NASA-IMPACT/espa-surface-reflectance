       subroutine invaero(xts,xtv,xfi,aot550nm,rolutt,
     s         pres,tpres,  
     s            transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s	     trotoa,erelc,ib1,ib2,raot550nm,roslamb1,residual,snext)
     
       real xts,xtv,xfi
       real aot550nm(22),raot550nm
       integer ib
       real rolutt(16,7,22,8000),nbfi(22,20),tsmax(22,20),tsmin(22,20)
       real transt(16,7,22,22),sphalbt(16,7,22),normext(16,7,22)
       real ttv(22,20),nbfic(22,20)
       real tts(22)
       real uoz,uwv
       real pres,tpres(7)
       real oztransa(16)
       real wvtransa(16),wvtransb(16),wvtransc(16)
       real ogtransa0(16),ogtransa1(16)
       real ogtransb0(16),ogtransb1(16)
       real ogtransc0(16),ogtransc1(16)
       real tauray(16)
       real roslamb,roatm,rotoa,next,snext
       real fac,pi
       real trotoa(16),erelc(16)
       integer indts(22)
       integer ib1,ib2
c        ib1=7
c	ib2=3
        iaot=1
        deltasrp=0.0      
  111    call atmcorlamb(xts,xtv,xfi,aot550nm(iaot),ib1,
     s       pres,tpres,aot550nm,
     s       rolutt,
     s            transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s	     trotoa(ib1),roslamb1,next)
C
       
       call atmcorlamb(xts,xtv,xfi,aot550nm(iaot),ib2,
     s          pres,tpres,aot550nm,
     s           rolutt,
     s            transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s	     trotoa(ib2),roslamb2,next)
       
       deltasr=roslamb2*erelc(ib1)-roslamb1*erelc(ib2)
       if ((deltasr.ge.0.).and.(iaot.le.21)) then
           iaot=iaot+1
	   deltasrp=deltasr
	   goto 111
	   endif
	if (iaot.ne.1) then
           slope=aot550nm(iaot)-aot550nm(iaot-1)
           slope=slope/(deltasr-deltasrp)
           raot550nm=aot550nm(iaot-1)-deltasrp*slope
           if (raot550nm.lt.0.01) raot550nm=0.01
	   else
	   raot550nm=0.01
	   endif
c       sanity check on upper limit of raot550nm, 22-JAN-07	   
	if (raot550nm.gt.5.0) then
	    raot550nm=4.99
	endif
c       write(6,*) "raot550nm ",raot550nm        
C       
C  bug May,30,2006 needs to recompute roslamb1 with new
C  aot before proceeding with residual computation
       call atmcorlamb(xts,xtv,xfi,raot550nm,ib1,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s	     trotoa(ib1),roslamb1,next) 

C   Added residual calculation to invaero (05-JUL-05)
  40   residual=0.0
       nband=0
       do ib=1,16 
       if (erelc(ib).gt.0.) then
       call atmcorlamb(xts,xtv,xfi,raot550nm,ib,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s	     trotoa(ib),roslamb,next)       
c May,30,2006 fix the computation of residual to have it in quadratic means versus
c absolute value average  
        if (ib.eq.ib3) snext=next   
       residual=residual+(roslamb*erelc(ib1)-roslamb1*erelc(ib))
     s  *(roslamb*erelc(ib1)-roslamb1*erelc(ib))
       nband=nband+1
       endif
       enddo
       residual=sqrt(residual)/(nband-1)
       
       return
       end
       
       subroutine invaeroocean(xts,xtv,xfi,aot550nm,rolutt,
     s         pres,tpres,  
     s            transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s	     trotoa,erelc,ib1,ib2,aot2,roslamb1,residual,angexp)
     
       real xts,xtv,xfi
       real aot550nm(22),raot550nm
       integer ib
       real rolutt(16,7,22,8000),nbfi(22,20),tsmax(22,20),tsmin(22,20)
       real transt(16,7,22,22),sphalbt(16,7,22)
       real ttv(22,20),nbfic(22,20)
       real tts(22)
       real uoz,uwv
       real pres,tpres(7)
       real oztransa(16)
       real wvtransa(16),wvtransb(16),wvtransc(16)
       real ogtransa0(16),ogtransa1(16)
       real ogtransb0(16),ogtransb1(16)
       real ogtransc0(16),ogtransc1(16)
       real tauray(16)
       real roslamb,roatm,rotoa
       real fac,pi
       real trotoa(16),erelc(16)
       integer indts(22)
       integer ib1,ib2
       real aot5,aot2
       real wave(16)
       real angexp
       character*80 err_msg
       integer retval
       data wave/670.0,870.0,480.,550.,1020.,1670.,2130.,
     s           412.,443.,490.,530.,545.,663.,670.,740.,865./
        iaot=1
        deltasrp=0.0      
	angexp=0.
  111    call atmcorocea2(xts,xtv,xfi,aot550nm(iaot),2,
     s       pres,tpres,aot550nm,
     s       rolutt,
     s            transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s	     trotoa(2),roslamb2,angexp,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval)
     
C do retrieval in band 5
C
       deltasr=roslamb2
       if ((deltasr.ge.0.).and.(iaot.le.21)) then
           iaot=iaot+1
	   deltasrp=deltasr
	   goto 111
	   endif
	if (iaot.ne.1) then
           slope=aot550nm(iaot)-aot550nm(iaot-1)
           slope=slope/(deltasr-deltasrp)
           raot550nm=aot550nm(iaot-1)-deltasrp*slope
           if (raot550nm.lt.0.01) raot550nm=0.01
	   else
	   raot550nm=0.01
	   endif
C        write(6,*) "raot550nm band 2",raot550nm        
	aot2=raot550nm
C       
        iaot=1
        deltasrp=0.0      
  112   angexp=0.
         call atmcorocea2(xts,xtv,xfi,aot550nm(iaot),5,
     s       pres,tpres,aot550nm,
     s       rolutt,
     s            transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s	     trotoa(5),roslamb2,angexp,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval)
C do retrieval in band 2
C
       deltasr=roslamb2
       if ((deltasr.ge.0.).and.(iaot.le.21)) then
           iaot=iaot+1
	   deltasrp=deltasr
	   goto 112
	   endif
	if (iaot.ne.1) then
           slope=aot550nm(iaot)-aot550nm(iaot-1)
           slope=slope/(deltasr-deltasrp)
           raot550nm=aot550nm(iaot-1)-deltasrp*slope
           if (raot550nm.lt.0.01) raot550nm=0.01
	   else
	   raot550nm=0.01
	   endif
C        write(6,*) "raot550nm band 5",raot550nm
	aot5=raot550nm  
	      
       angexp=log(aot2/aot5)/log(870./1020.)
c threshold on angexp for Urban clean model
       if (angexp.lt.-2.) angexp=-2.
       if (angexp.gt.1.0) angexp=1.0       
       
C       write(6,*) "angexp ",angexp
       residual=0.0
       nband=0
       do ib=1,16 
       if (erelc(ib).gt.-1.0) then
        call atmcorocea2(xts,xtv,xfi,aot2,ib,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s	     trotoa(ib),roslamb,angexp,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval)   
       if (erelc(ib).gt.0.0) then
       residual=residual+sqrt((roslamb*erelc(ib1)-roslamb1*erelc(ib))
     s  *(roslamb*erelc(ib1)-roslamb1*erelc(ib)))
       else
       residual=residual+sqrt((roslamb*roslamb))
       endif
       nband=nband+1
       endif
       enddo
       residual=residual/(nband-1)
       
       return
       end
       
       
       subroutine atmcorocea2(xts,xtv,xfi,aot2,ib,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s	     rotoa,roslamb,angexp,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval)
     
       parameter (fac = 0.017453293)
       real xts,xtv,xfi
       real aot550nm(22),raot550nm
       integer ib
       real rolutt(16,7,22,8000),nbfi(22,20),tsmax(22,20),tsmin(22,20)
       real transt(16,7,22,22),sphalbt(16,7,22)
       real ttv(22,20),nbfic(22,20)
       real tts(22)
       real uoz,uwv
       real oztransa(16)
       real wvtransa(16),wvtransb(16),wvtransc(16)
       real ogtransa0(16),ogtransa1(16)
       real ogtransb0(16),ogtransb1(16)
       real ogtransc0(16),ogtransc1(16)
       real tauray(16)
       real roslamb,roatm,rotoa
C       real fac,pi
       integer indts(22)
       real tgo,ttatmg,satm,next
       real tpres(7),pres
       real aot5,aot2
       real wave(16)
       real angexp
       character*80 err_msg
       integer retval
        data wave/670.0,870.0,480.,550.,1020.,1670.,2130.,
     s           412.,443.,490.,530.,545.,663.,670.,740.,865./
     
C       pi=acos(0.)*2.
C       fac=pi/180.
       
       if ((ib.ne.7).and.(ib.ne.6)) then
       raot550nm=aot2*exp(log(wave(ib)/870.)*angexp)
       else
       raot550nm=0.01
       endif
C
C     New, 21-AUG-06: reject raot550nm values above aot550nm(22), which is 5.000 (just for ocean)
C       
       if (raot550nm .gt. aot550nm(22)) then
       raot550nm=aot550nm(22)
       endif
       if (raot550nm .lt. aot550nm(1)) then
       raot550nm=aot550nm(1)
       endif
    
C       write(6,*)raot550nm, angexp, aot2*exp(log(wave(ib)/870.)*angexp)
       call comproatm(xts,xtv,xfi,raot550nm,ib,pres,tpres,
     s      aot550nm,rolutt,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s			roatm,err_msg,retval)
     

     
       call comptrans(xts,raot550nm,ib,pres,tpres,aot550nm,transt,
     s             xtsstep,xtsmin,tts,  
     s                      xtts,err_msg,retval)
       call comptrans(xtv,raot550nm,ib,pres,tpres,aot550nm,transt,
     s             xtvstep,xtvmin,tts,  
     s                      xttv,err_msg,retval)
C Compute total transmission (product downward by  upward)
       ttatm=xtts*xttv
    
C Compute SPHERICAL ALBEDO
       call compsalb(raot550nm,ib,pres,tpres,aot550nm,sphalbt,normext,
     s                      satm)
     
       call comptg(ib,xts,xtv,uoz,uwv,pres,
     s ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s ogtransc0,ogtransc1,
     s wvtransa,wvtransb,wvtransc,oztransa,  
     s            tgoz,tgwv,tgwvhalf,tgog)
     
	
C compute rayleigh component (intrinsic reflectance, at p0)
        xtaur=tauray(ib)
        xphi=xfi
        xmus=cos(xts*fac)
        xmuv=cos(xtv*fac)
C compute rayleigh component (intrinsic reflectance, at p=pres)
        xtaur=tauray(ib)*pres/1013.0
c	write(6,*) "for chand ",xtaur,xphi,xmus,xmuv,fac
        call local_chand(xphi,xmus,xmuv,xtaur,xrorayp)
	

C Perform atmospheric correction 
        roslamb=rotoa/tgog/tgoz
        roslamb=(roslamb-(roatm-xrorayp)*tgwvhalf-xrorayp)
c	write(6,*) "roatm ",(roatm-xrorayp)*tgwvhalf+xrorayp
c	write(6,*) "ttatm ",ttatm
c	write(6,*) "satm ",satm
        roslamb=roslamb/ttatm/tgwv
        roslamb=roslamb/(1.+satm*roslamb)
c        write(6,*) "roslamb: ", roslamb
C  New, 07-JUL-05 from atmcorlamb2:
	tgo=(tgog*tgoz)
        roatm=(roatm-xrorayp)*tgwvhalf+xrorayp
        ttatmg=(ttatm*tgwv)

	return
	end
       
       subroutine atmcorlamb2(xts,xtv,xfi,raot550nm,ib,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s	     rotoa,roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval,eps)
C
C  07-JUL-05: this routine returns variables for calculating roslamb.
C     
       parameter (fac = 0.017453293)
       real xts,xtv,xfi
       real aot550nm(22),raot550nm,mraot550nm
       integer ib
       real rolutt(16,7,22,8000),nbfi(22,20),tsmax(22,20),tsmin(22,20)
       real transt(16,7,22,22),sphalbt(16,7,22),normext(16,7,22)
       real ttv(22,20),nbfic(22,20)
       real tts(22)
       real uoz,uwv
       real oztransa(16)
       real wvtransa(16),wvtransb(16),wvtransc(16)
       real ogtransa0(16),ogtransa1(16)
       real ogtransb0(16),ogtransb1(16)
       real ogtransc0(16),ogtransc1(16)
       real tauray(16)
       real lambda(16)
       real roslamb,roatm,rotoa,next
C       real fac,pi
       integer indts(22)
       real tgo,ttatmg,satm
       real tpres(7),pres
       character*80 err_msg
       integer retval
       data lambda/0.443,0.480,0.585,0.655,0.865,1.61,2.2,
     s  4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0/ 
        real eps
C       pi=acos(0.)*2.
C       fac=pi/180.
       
c       eps=1.
       if  (eps.lt.0.) then
       mraot550nm=raot550nm
       else
       if (ib.le.7) then
       mraot550nm=(raot550nm/normext(ib,1,4))*((lambda(ib)/0.55)**(-eps))
       else
       mraot550nm=raot550nm
       endif
       endif
       retval=0
       call comproatm(xts,xtv,xfi,mraot550nm,ib,pres,tpres,
     s      aot550nm,rolutt,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s			roatm,err_msg,retval)
     
       if (retval .ne. 0) then
          return
       endif
     
       retval=0
       call comptrans(xts,mraot550nm,ib,pres,tpres,aot550nm,transt,
     s             xtsstep,xtsmin,tts,  
     s                      xtts,err_msg,retval)
       if (retval .ne. 0) then
          return
       endif
       
       retval=0
       call comptrans(xtv,mraot550nm,ib,pres,tpres,aot550nm,transt,
     s             xtvstep,xtvmin,tts,  
     s                      xttv,err_msg,retval)
       if (retval .ne. 0) then
          return
       endif
     
C Compute total transmission (product downward by  upward)
       ttatm=xtts*xttv
    
C Compute SPHERICAL ALBEDO
       call compsalb(mraot550nm,ib,pres,tpres,aot550nm,sphalbt,normext,
     s                      satm,next)
     
       call comptg(ib,xts,xtv,uoz,uwv,pres,
     s ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s ogtransc0,ogtransc1,
     s wvtransa,wvtransb,wvtransc,oztransa,  
     s            tgoz,tgwv,tgwvhalf,tgog)
     
	
C compute rayleigh component (intrinsic reflectance, at p0)
        xtaur=tauray(ib)
        xphi=xfi
        xmus=cos(xts*fac)
        xmuv=cos(xtv*fac)
C compute rayleigh component (intrinsic reflectance, at p=pres)
        xtaur=tauray(ib)*pres/1013.0
c	write(6,*) "for local_chand ",xtaur,xphi,xmus,xmuv,fac
        call local_chand(xphi,xmus,xmuv,xtaur,xrorayp)
	

C Perform atmospheric correction 
        roslamb=rotoa/(tgog*tgoz)
        roslamb=(roslamb-(roatm-xrorayp)*tgwvhalf-xrorayp)
c	write(6,*) "roatm ",(roatm-xrorayp)*tgwvhalf+xrorayp
c	write(6,*) "ttatm ",ttatm      
c	write(6,*) "satm ",satm
        roslamb=roslamb/(ttatm*tgwv)
        roslamb=roslamb/(1.+satm*roslamb)
c        write(6,*) "roslamb: ", roslamb
C
C  New, 07-JUL-05:
	tgo=(tgog*tgoz)
        roatm=(roatm-xrorayp)*tgwvhalf+xrorayp
        ttatmg=(ttatm*tgwv)
C	
	return
	end

       subroutine raycorlamb2(xts,xtv,xfi,ib,pres,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s	     rotoa,roslamb,tgo,roatm,ttatmg,satm,xrorayp)
C
C  07-JUL-05: this routine returns variables for calculating roslamb.
C     
       parameter (fac = 0.017453293)
       real xts,xtv,xfi
       integer ib
       real uoz,uwv
       real oztransa(16)
       real wvtransa(16),wvtransb(16),wvtransc(16)
       real ogtransa0(16),ogtransa1(16)
       real ogtransb0(16),ogtransb1(16)
       real ogtransc0(16),ogtransc1(16)
       real tauray(16)
       real roslamb,roatm,rotoa
C       real fac,pi
       integer indts(22)
       real tgo,ttatmg,satm
       real pres
C       pi=acos(0.)*2.
C       fac=pi/180.
       
     
       xtaur=tauray(ib)*pres/1013.0
       xmus=cos(xts*fac)
       xmuv=cos(xtv*fac)
       call comptransray(xtaur,xmus,xtts)
       call comptransray(xtaur,xmuv,xttv)

C Compute total transmission (product downward by  upward)
       ttatm=xtts*xttv
    
C Compute SPHERICAL ALBEDO
       call local_csalbr(xtaur,satm)

       call comptg(ib,xts,xtv,uoz,uwv,pres,
     s ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s ogtransc0,ogtransc1,
     s wvtransa,wvtransb,wvtransc,oztransa,  
     s            tgoz,tgwv,tgwvhalf,tgog)
     
	
C compute rayleigh component (intrinsic reflectance, at p0)
        xtaur=tauray(ib)
        xphi=xfi
C        xmus=cos(xts*fac)   calculated already
C        xmuv=cos(xtv*fac)
C compute rayleigh component (intrinsic reflectance, at p=pres)
        xtaur=tauray(ib)*pres/1013.0
c	write(6,*) "for local_chand ",xtaur,xphi,xmus,xmuv,fac
        call local_chand(xphi,xmus,xmuv,xtaur,xrorayp)
	

C Perform atmospheric correction 
        roatm=xrorayp
        roslamb=rotoa/(tgog*tgoz)
        roslamb=(roslamb-(roatm-xrorayp)*tgwvhalf-xrorayp)
c	write(6,*) "roatm ",(roatm-xrorayp)*tgwvhalf+xrorayp
c	write(6,*) "ttatm ",ttatm      
c	write(6,*) "satm ",satm
        roslamb=roslamb/(ttatm*tgwv)
        roslamb=roslamb/(1.+satm*roslamb)
c        write(6,*) "roslamb: ", roslamb
C
C  New, 07-JUL-05:
	tgo=(tgog*tgoz)
        roatm=(roatm-xrorayp)*tgwvhalf+xrorayp
        ttatmg=(ttatm*tgwv)
C	
	return
	end

       subroutine atmcorlamb(xts,xtv,xfi,raot550nm,ib,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s	     rotoa,roslamb,next)
     
       parameter (fac = 0.017453293)
       real xts,xtv,xfi
       real aot550nm(22),raot550nm
       integer ib
       real rolutt(16,7,22,8000),nbfi(22,20),tsmax(22,20),tsmin(22,20)
       real transt(16,7,22,22),sphalbt(16,7,22),normext(16,7,22)
       real ttv(22,20),nbfic(22,20)
       real tts(22)
       real uoz,uwv
       real oztransa(16)
       real wvtransa(16),wvtransb(16),wvtransc(16)
       real ogtransa0(16),ogtransa1(16)
       real ogtransb0(16),ogtransb1(16)
       real ogtransc0(16),ogtransc1(16)
       real tauray(16)
       real roslamb,roatm,rotoa,next
C       real fac,pi
       integer indts(22)
       real tpres(7),pres
       character*80 err_msg
       integer retval
C       pi=acos(0.)*2.
C       fac=pi/180.

       call comproatm(xts,xtv,xfi,raot550nm,ib,pres,tpres,
     s      aot550nm,rolutt,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s			roatm,err_msg,retval)
     

     
       call comptrans(xts,raot550nm,ib,pres,tpres,aot550nm,transt,
     s             xtsstep,xtsmin,tts,  
     s                      xtts,err_msg,retval)
       call comptrans(xtv,raot550nm,ib,pres,tpres,aot550nm,transt,
     s             xtvstep,xtvmin,tts,  
     s                      xttv,err_msg,retval)
C Compute total transmission (product downward by  upward)
       ttatm=xtts*xttv
    
C Compute SPHERICAL ALBEDO
       call compsalb(raot550nm,ib,pres,tpres,aot550nm,sphalbt,normext,
     s                      satm,next)
     
       call comptg(ib,xts,xtv,uoz,uwv,pres,
     s ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s ogtransc0,ogtransc1,
     s wvtransa,wvtransb,wvtransc,oztransa,  
     s            tgoz,tgwv,tgwvhalf,tgog)
     
	
C compute rayleigh component (intrinsic reflectance, at p0)
        xtaur=tauray(ib)
        xphi=xfi
        xmus=cos(xts*fac)
        xmuv=cos(xtv*fac)
C compute rayleigh component (intrinsic reflectance, at p=pres)
        xtaur=tauray(ib)*pres/1013.0
c	write(6,*) "for local_chand ",xtaur,xphi,xmus,xmuv,fac
        call local_chand(xphi,xmus,xmuv,xtaur,xrorayp)
	

C Perform atmospheric correction 
        roslamb=rotoa/(tgog*tgoz)
        roslamb=(roslamb-(roatm-xrorayp)*tgwvhalf-xrorayp)
c	write(6,*) "roatm ",(roatm-xrorayp)*tgwvhalf+xrorayp
c	write(6,*) "ttatm ",ttatm
c	write(6,*) "satm ",satm
        roslamb=roslamb/(ttatm*tgwv)
        roslamb=roslamb/(1.+satm*roslamb)
c        write(6,*) "roslamb: ", roslamb

	return
	end


      subroutine local_csalbr(xtau,xalb)
      real xtau,xalb,fintexp3
      xalb=(3*xtau-fintexp3(xtau)*(4+2*xtau)+2*exp(-xtau))
      xalb=xalb/(4.+3*xtau)
      return
      end
      real function fintexp3(xtau)
      real xx,xtau,fintexp1
      xx=(exp(-xtau)*(1.-xtau)+xtau*xtau*fintexp1(xtau))/2.
      fintexp3=xx
      return
      end
      real function fintexp1(xtau)
c accuracy 2e-07... for 0<xtau<1
      real xx,a(0:5),xtau,xftau
      integer i
      data (a(i),i=0,5) /-.57721566,0.99999193,-0.24991055,
     c                  0.05519968,-0.00976004,0.00107857/
      xx=a(0)
      xftau=1.
      do i=1,5
      xftau=xftau*xtau
      xx=xx+a(i)*xftau
      enddo
      fintexp1=xx-log(xtau)
      return
      end

	subroutine local_chand (xphi,xmuv,xmus,xtau
     s			,xrray)
c input parameters: xphi,xmus,xmuv,xtau
c xphi: azimuthal difference between sun and observation (xphi=0,
c in backscattering) and expressed in degree (0.:360.)
c xmus: cosine of the sun zenith angle
c xmuv: cosine of the observation zenith angle
c xtau: molecular optical depth
c output parameter: xrray : molecular reflectance (0.:1.)
c constant : xdep: depolarization factor (0.0279)
        parameter (fac = 0.017453293)
	real xdep,pl(10)
	real fs0,fs1,fs2
	real as0(10),as1(2),as2(2)
        real xphi,xmus,xmuv,xtau,xrray,pi,phios,xcosf1,xcosf2
        real xcosf3,xbeta2,xfd,xph1,xph2,xph3,xitm, xp1, xp2, xp3
        real cfonc1,cfonc2,cfonc3,xlntau,xitot1,xitot2,xitot3
        integer i
	data (as0(i),i=1,10) /.33243832,-6.777104e-02,.16285370
     s	,1.577425e-03,-.30924818,-1.240906e-02,-.10324388
     s	,3.241678e-02,.11493334,-3.503695e-02/
	data (as1(i),i=1,2) /.19666292, -5.439061e-02/
	data (as2(i),i=1,2) /.14545937,-2.910845e-02/
C	pi=3.1415927
C	fac=pi/180.
	phios=180.-xphi
	xcosf1=1.
	xcosf2=cos(phios*fac)
	xcosf3=cos(2*phios*fac)
	xbeta2=0.5
	xdep=0.0279
	xfd=xdep/(2-xdep)
	xfd=(1-xfd)/(1+2*xfd)
	xph1=1+(3*xmus*xmus-1)*(3*xmuv*xmuv-1)*xfd/8.
	xph2=-xmus*xmuv*sqrt(1-xmus*xmus)*sqrt(1-xmuv*xmuv)
	xph2=xph2*xfd*xbeta2*1.5
	xph3=(1-xmus*xmus)*(1-xmuv*xmuv)
	xph3=xph3*xfd*xbeta2*0.375
	xitm=(1-exp(-xtau*(1/xmus+1/xmuv)))*xmus/(4*(xmus+xmuv))
	xp1=xph1*xitm
	xp2=xph2*xitm
	xp3=xph3*xitm
	xitm=(1-exp(-xtau/xmus))*(1-exp(-xtau/xmuv))
	cfonc1=xph1*xitm
	cfonc2=xph2*xitm
	cfonc3=xph3*xitm
	xlntau=log(xtau)
	pl(1)=1.
	pl(2)=xlntau
	pl(3)=xmus+xmuv
	pl(4)=xlntau*pl(3)
	pl(5)=xmus*xmuv
	pl(6)=xlntau*pl(5)
	pl(7)=xmus*xmus+xmuv*xmuv
	pl(8)=xlntau*pl(7)
	pl(9)=xmus*xmus*xmuv*xmuv
	pl(10)=xlntau*pl(9)
	fs0=0.
	do i=1,10
	fs0=fs0+pl(i)*as0(i)
	enddo
	fs1=pl(1)*as1(1)+pl(2)*as1(2)
	fs2=pl(1)*as2(1)+pl(2)*as2(2)
	xitot1=xp1+cfonc1*fs0*xmus
	xitot2=xp2+cfonc2*fs1*xmus
	xitot3=xp3+cfonc3*fs2*xmus
	xrray=xitot1*xcosf1
	xrray=xrray+xitot2*xcosf2*2
	xrray=xrray+xitot3*xcosf3*2
	xrray=xrray/xmus
	return
	end
	
       subroutine comptransray(xtaur,xmus,ttray)
       real xtaur,xmus,ttray,ddiftt,ddirtt
       
       ddiftt=(2./3.+xmus)+(2./3.-xmus)*exp(-xtaur/xmus)
       ddiftt=ddiftt/((4./3.)+xtaur)-exp(-xtaur/xmus)
       ddirtt=exp(-xtaur/xmus)
       ttray=ddirtt+ddiftt
       return
       end

       subroutine comptg(ib,xts,xtv,uoz,uwv,pres,
     s ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s ogtransc0,ogtransc1,
     s wvtransa,wvtransb,wvtransc,oztransa,  
     s            tgoz,tgwv,tgwvhalf,tgog)
       parameter (fac = 0.017453293)
       integer ib
       real xts,xtv,uoz,uwv,pres,tgoz,tgwv,tgwvhalf,tgog
       real ogtransa0(16),ogtransa1(16),ogtransb0(16),ogtransb1(16)
       real ogtransc0(16),ogtransc1(16)
       real wvtransa(16),wvtransb(16),wvtransc(16)
       real oztransa(16)
       real pi,m
       
C       pi=atan(1.)*4.
C       m=1./cos(xts*pi/180.0)+1./cos(xtv*pi/180.0)
       m=1./cos(xts*fac)+1./cos(xtv*fac)
       tgoz=exp(oztransa(ib)*m*uoz)
	
C compute water vapor transmission
        a=wvtransa(ib)
        b=wvtransb(ib)
c        c=wvtransc(ib)
	x=m*uwv
        if (x.gt.1.E-06) then
        tgwv=exp(-a*(exp(log(x)*b)))
        else
        tgwv=1.
        endif
C compute water vapor transmission half the content
	x=m*uwv/2.
        if (x.gt.1.E-06) then
        tgwvhalf=exp(-a*(exp(log(x)*b)))
        else
        tgwvhalf=1.
        endif
C compute other gases transmission
        a1=ogtransa1(ib)
	b0=ogtransb0(ib)
	b1=ogtransb1(ib)
	p=pres/1013.
c p is the pressure in atmosphere	
        tgog=-(a1*p)*(m**(exp(-(b0+b1*p))))
        tgog=exp(tgog)
        return
	end


       subroutine compsalb(raot550nm,ib,pres,tpres,aot550nm,sphalbt,normext,
     s                      satm,next)
       real aot550nm(22),raot550nm,sphalbt(16,7,22),normext(16,7,22)
       real xtiaot1,xtiaot2
       integer iaot1,iaot2
       integer ip1,ip2,ip
       real satm1,satm2,satm,next1,next2,next
       real tpres(7)
       real pres,dpres
       real deltaaot
       integer ib
       
        ip1=1
        do ip=1,7
        if (pres.lt.tpres(ip)) then
	     ip1=ip
	     endif
	enddo   
	if (ip1.eq.7) ip1=6  
        ip2=ip1+1	
      
       
        iaot1=1
        do iaot=1,22
        if (raot550nm.gt.aot550nm(iaot)) then
	     iaot1=iaot
	     endif
	enddo   
	if (iaot1.eq.22) iaot1=21
       iaot2=iaot1+1	
       
                 
       xtiaot1=sphalbt(ib,ip1,iaot1)
       xtiaot2=sphalbt(ib,ip1,iaot2)
       deltaaot=(raot550nm-aot550nm(iaot1))
       deltaaot=deltaaot/(aot550nm(iaot2)-aot550nm(iaot1))
       satm1=xtiaot1+(xtiaot2-xtiaot1)*deltaaot
       
       xtiaot1=sphalbt(ib,ip2,iaot1)
       xtiaot2=sphalbt(ib,ip2,iaot2)
       deltaaot=(raot550nm-aot550nm(iaot1))
       deltaaot=deltaaot/(aot550nm(iaot2)-aot550nm(iaot1))
       satm2=xtiaot1+(xtiaot2-xtiaot1)*deltaaot
       
       dpres=(pres-tpres(ip1))/(tpres(ip2)-tpres(ip1))
       satm=satm1+(satm2-satm1)*dpres
       
       xtiaot1=normext(ib,ip1,iaot1)
       xtiaot2=normext(ib,ip1,iaot2)
       deltaaot=(raot550nm-aot550nm(iaot1))
       deltaaot=deltaaot/(aot550nm(iaot2)-aot550nm(iaot1))
       next1=xtiaot1+(xtiaot2-xtiaot1)*deltaaot
       
       xtiaot1=normext(ib,ip2,iaot1)
       xtiaot2=normext(ib,ip2,iaot2)
       deltaaot=(raot550nm-aot550nm(iaot1))
       deltaaot=deltaaot/(aot550nm(iaot2)-aot550nm(iaot1))
       next2=xtiaot1+(xtiaot2-xtiaot1)*deltaaot
       
       dpres=(pres-tpres(ip1))/(tpres(ip2)-tpres(ip1))
       next=next1+(next2-next1)*dpres
       
       return
       end
       
       
       subroutine comptrans(xts,raot550nm,ib,pres,tpres,aot550nm,transt,
     s             xtsstep,xtsmin,tts,  
     s                      xtts,err_msg,retval)
       real xts,xtts,tts(22)
       real aot550nm(22),raot550nm,transt(16,7,22,22)
       real xtiaot1,xtiaot2
       integer iaot1,iaot2
       integer ip1,ip2,ip
       real xtts1,xtts2
       real tpres(7)
       real pres,dpres
       real xtsstep,xtsmin,xmts,xtranst
       character*80 err_msg
       integer retval

       
        ip1=1
        do ip=1,7
        if (pres.lt.tpres(ip)) then
	     ip1=ip
	     endif
	enddo   
	if (ip1.eq.7) ip1=6  
        ip2=ip1+1	

        iaot1=1
        do iaot=1,22
        if (raot550nm.gt.aot550nm(iaot)) then
	     iaot1=iaot
	     endif
	enddo   
	
	if (iaot1.eq.22) iaot1=21 
       if (xts.le.xtsmin) then
          its=1
	  else
          its=int((xts-xtsmin)/xtsstep)+1
       endif
       
       if (its.gt.20) then
C           write(6,*) "Xts is too large",xts
C           stop
            write(err_msg,*) "Xts is too large",xts," "//char(0)
            retval=2
            return
       endif
!
! 18-JUL-05 Changes to save some time.
!       
!       xmts=(xts-tts(its))/4.0
!       xtiaot1=transt(ib,ip1,iaot1,its)+(transt(ib,ip1,iaot1,its+1)-
!     &   transt(ib,ip1,iaot1,its))*xmts
!       iaot2=iaot1+1	          
!       xtiaot2=transt(ib,ip1,iaot2,its)+(transt(ib,ip1,iaot2,its+1)-
!     &   transt(ib,ip1,iaot2,its))*xmts
!       deltaaot=(raot550nm-aot550nm(iaot1))
!       deltaaot=deltaaot/(aot550nm(iaot2)-aot550nm(iaot1))
!       xtts1=xtiaot1+(xtiaot2-xtiaot1)*deltaaot
!       
!       xtiaot1=transt(ib,ip2,iaot1,its)+(transt(ib,ip2,iaot1,its+1)-
!     &   transt(ib,ip2,iaot1,its))*xmts
!c       iaot2=iaot1+1	          
!       xtiaot2=transt(ib,ip2,iaot2,its)+(transt(ib,ip2,iaot2,its+1)-
!     &   transt(ib,ip2,iaot2,its))*xmts
!c       deltaaot=(raot550nm-aot550nm(iaot1))
!c       deltaaot=deltaaot/(aot550nm(iaot2)-aot550nm(iaot1))
!       xtts2=xtiaot1+(xtiaot2-xtiaot1)*deltaaot

       xmts=(xts-tts(its))/4.0
       xtranst=transt(ib,ip1,iaot1,its)
       xtiaot1=xtranst+(transt(ib,ip1,iaot1,its+1)-
     &   xtranst)*xmts
       iaot2=iaot1+1	          
       xtranst=transt(ib,ip1,iaot2,its)
       xtiaot2=xtranst+(transt(ib,ip1,iaot2,its+1)-
     &   xtranst)*xmts
       deltaaot=(raot550nm-aot550nm(iaot1))
       deltaaot=deltaaot/(aot550nm(iaot2)-aot550nm(iaot1))
       xtts1=xtiaot1+(xtiaot2-xtiaot1)*deltaaot
       
       xtranst=transt(ib,ip2,iaot1,its)
       xtiaot1=xtranst+(transt(ib,ip2,iaot1,its+1)-
     &   xtranst)*xmts
c       iaot2=iaot1+1	          
       xtranst=transt(ib,ip2,iaot2,its)
       xtiaot2=xtranst+(transt(ib,ip2,iaot2,its+1)-
     &   xtranst)*xmts
c       deltaaot=(raot550nm-aot550nm(iaot1))
c       deltaaot=deltaaot/(aot550nm(iaot2)-aot550nm(iaot1))
       xtts2=xtiaot1+(xtiaot2-xtiaot1)*deltaaot


       dpres=(pres-tpres(ip1))/(tpres(ip2)-tpres(ip1))
       xtts=xtts1+(xtts2-xtts1)*dpres
        
       return
       end
       

       subroutine comproatm(xts,xtv,xfi,raot550nm,ib,pres,tpres,
     s       aot550nm,rolutt,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s			roatm,err_msg,retval)
       parameter (fac = 0.017453293)
       real nbfic(22,20),tts(22),xtsmin,xtsstep,xtvmin,xtvstep,pi
       real rolutt(16,7,22,8000),nbfi(22,20),tsmax(22,20),tsmin(22,20)
       real ttv(22,20)
       integer indts(22)
       real xts,xtv,xfi,xmus,xmuv,scaa,xtsmax
       real aot550nm(22),raot550nm       
       integer its,itv,ita1,iaot,iaot2
       integer ip1,ip2,ip
       real rop1,rop2
       real tpres(7)
       real pres,dpres
       character*80 err_msg
       integer retval
       real logaot550nm(22)
C       data aot550nm /0.01,0.05,0.10,0.15,0.20,0.30,0.40,0.60,
C     &                0.80,1.00,1.20,1.40,1.60,1.80,2.00,
C     &                2.30,2.60,3.00,3.50,4.00,4.50,5.00/
       data logaot550nm /-4.605170186,-2.995732274,-2.302585093,
     &                   -1.897119985,-1.609437912,-1.203972804,
     &                   -0.916290732,-0.510825624,-0.223143551,
     &                    0.000000000,0.182321557,0.336472237,
     &                    0.470003629,0.587786665,0.693157181,
     &                    0.832909123,0.955511445,1.098612289,
     &                    1.252762969,1.386294361,1.504077397,
     &                    1.609437912/
       real nbfic1, nbfic2, nbfic3, nbfic4
       real nbfi1, nbfi2, nbfi3, nbfi4
       integer get_value, gt
       real result

        ip1=1
        do ip=1,7
        if (pres.lt.tpres(ip)) then
	     ip1=ip
	     endif
	enddo   
	if (ip1.eq.7) ip1=6  
        ip2=ip1+1	
       
       xtsmin=0
       xtsstep=4.0
       xtvmin=2.84090
       xtvstep=(6.52107-2.84090)
C       pi=acos(0.0)*2.0
       if (xtv.le.xtvmin) then 
           itv=1
	   else
	   itv=int((xtv-xtvmin)/xtvstep)+2
       endif
       its=int((xts-xtsmin)/xtsstep)+1
       if (its.gt.20) then
C          write(6,*) "Xts is too large",xts
C          stop
          write(err_msg,*) "Xts is too large",xts," "//char(0)
          retval = 2
          return
       endif
C       fac=pi/180.
       xmuv=cos(xtv*fac)  
       xmus=cos(xts*fac)
       cscaa=-xmus*xmuv-cos(xfi*fac)
     S  *sqrt(1.-xmus*xmus)*sqrt(1.-xmuv*xmuv)
       scaa=acos(cscaa)/fac
c        write(6,*) xts,tts(its),tts(its+1),xtv,ttv(its,itv),
c     S  ttv(its,itv+1),scaa
C determine lower and uper limit in aot
        iaot1=1
        do iaot=1,22
        if (raot550nm.gt.aot550nm(iaot)) then
	     iaot1=iaot
	     endif
	enddo   
	if (iaot1.eq.22) iaot1=21 
	
	nbfic1=nbfic(its,itv)
	nbfi1=nbfi(its,itv)
	nbfic2=nbfic(its+1,itv)
	nbfi2=nbfi(its+1,itv)
	nbfic3=nbfic(its,itv+1)
	nbfi3=nbfi(its,itv+1)
	nbfic4=nbfic(its+1,itv+1)
	nbfi4=nbfi(its+1,itv+1)
	
C compute for ip1
C Compute for iaot1	          
C interpolate point 1 (its,itv) vs scattering angle 
       xtsmax=tsmax(its,itv)
       if ((its.ne.1).and.(itv.ne.1)) then
          isca=int((xtsmax-scaa)/4.0)+1
	  if (isca.le.0) isca=1
          if ((isca+1).lt.nbfi(its,itv)) then
!             sca1=tsmax(its,itv)-(isca-1)*4.0
!             sca2=tsmax(its,itv)-isca*4.0
             sca1=xtsmax-(isca-1)*4.0
             sca2=xtsmax-isca*4.0
             else
 	     isca=nbfi(its,itv)-1
!             sca1=tsmax(its,itv)-(isca-1)*4.0
             sca1=xtsmax-(isca-1)*4.0
	     sca2=tsmin(its,itv)
	     endif
!	   iindex=indts(its)+nbfic(its,itv)-nbfi(its,itv)+isca
	   iindex=indts(its)+nbfic1-nbfi1+isca
C	   
C	   gt = get_value(ib,ip1,iaot1,iindex,result)
C	   write(6,*)result, rolutt(ib,ip1,iaot1,iindex)
C
	   roinf=rolutt(ib,ip1,iaot1,iindex)
	   rosup=rolutt(ib,ip1,iaot1,iindex+1)
	   ro1=roinf+(rosup-roinf)*(scaa-sca1)/(sca2-sca1)
           else
!           sca1=tsmax(its,itv)
!           sca2=tsmax(its,itv)
           sca1=xtsmax
           sca2=xtsmax
!	   iindex=indts(its)+nbfic(its,itv)-nbfi(its,itv)+1
	   iindex=indts(its)+nbfic1-nbfi1+1
	   roinf=rolutt(ib,ip1,iaot1,iindex)
!	   rosup=rolutt(ib,ip1,iaot1,iindex)
!	   ro1=rolutt(ib,ip1,iaot1,iindex)
	   rosup=roinf
	   ro1=roinf
           endif
C        write(6,*) "Point its,itv sca,sca1,sca2,ro..",scaa,sca1,sca2,ro1,
C     S             roinf,rosup
C interpolate point 2 (its+1,itv) vs scattering angle 
       xtsmax=tsmax(its+1,itv)
       if ((itv.ne.1)) then
          isca=int((xtsmax-scaa)/4.0)+1
	  if (isca.le.0) isca=1
          if ((isca+1).lt.nbfi(its+1,itv)) then
!             sca1=tsmax(its+1,itv)-(isca-1)*4.0
!             sca2=tsmax(its+1,itv)-isca*4.0
             sca1=xtsmax-(isca-1)*4.0
             sca2=xtsmax-isca*4.0
             else
 	     isca=nbfi(its+1,itv)-1
!             sca1=tsmax(its+1,itv)-(isca-1)*4.0
             sca1=xtsmax-(isca-1)*4.0
	     sca2=tsmin(its+1,itv)
	     endif
!	   iindex=indts(its+1)+nbfic(its+1,itv)-nbfi(its+1,itv)+isca
	   iindex=indts(its+1)+nbfic2-nbfi2+isca
	   roinf=rolutt(ib,ip1,iaot1,iindex)
	   rosup=rolutt(ib,ip1,iaot1,iindex+1)
	   ro2=roinf+(rosup-roinf)*(scaa-sca1)/(sca2-sca1)
           else
!           sca1=tsmax(its+1,itv)
!           sca2=tsmax(its+1,itv)
           sca1=xtsmax
           sca2=xtsmax
!	   iindex=indts(its+1)+nbfic(its+1,itv)-nbfi(its+1,itv)+1
	   iindex=indts(its+1)+nbfic2-nbfi2+1
	   roinf=rolutt(ib,ip1,iaot1,iindex)
!	   rosup=rolutt(ib,ip1,iaot1,iindex)
!	   ro2=rolutt(ib,ip1,iaot1,iindex)
	   rosup=roinf
	   ro2=roinf
           endif
c       write(6,*) "Point its+1,itv sca,sca1,sca2,ro..",scaa,sca1,sca2,
c     S            ro2,roinf,rosup
C interpolate point 3 (its,itv+1) vs scattering angle 
       xtsmax=tsmax(its,itv+1)
       if ((its.ne.1)) then
          isca=int((xtsmax-scaa)/4.0)+1
	  if (isca.le.0) isca=1
          if ((isca+1).lt.nbfi(its,itv+1)) then
!             sca1=tsmax(its,itv+1)-(isca-1)*4.0
!             sca2=tsmax(its,itv+1)-isca*4.0
             sca1=xtsmax-(isca-1)*4.0
             sca2=xtsmax-isca*4.0
             else
	     isca=nbfi(its,itv+1)-1
!             sca1=tsmax(its,itv+1)-(isca-1)*4.0
             sca1=xtsmax-(isca-1)*4.0
	     sca2=tsmin(its,itv+1)
	     endif
!           iindex=indts(its)+nbfic(its,itv+1)-nbfi(its,itv+1)+isca
           iindex=indts(its)+nbfic3-nbfi3+isca
	   roinf=rolutt(ib,ip1,iaot1,iindex)
	   rosup=rolutt(ib,ip1,iaot1,iindex+1)
	   ro3=roinf+(rosup-roinf)*(scaa-sca1)/(sca2-sca1)
           else
!           sca1=tsmax(its,itv+1)
!           sca2=tsmax(its,itv+1)
           sca1=xtsmax
           sca2=xtsmax
!	   iindex=indts(its)+nbfic(its,itv+1)-nbfi(its,itv+1)+1
	   iindex=indts(its)+nbfic3-nbfi3+1
	   roinf=rolutt(ib,ip1,iaot1,iindex)
!	   rosup=rolutt(ib,ip1,iaot1,iindex)
!	   ro3=rolutt(ib,ip1,iaot1,iindex)
	   rosup=roinf
	   ro3=roinf
           endif
c       write(6,*) "Point its,itv+1 sca,sca1,sca2,ro..",scaa,sca1,sca2,
c     S            ro3,roinf,rosup
C interpolate point 4 (its+1,itv+1) vs scattering angle 
          xtsmax=tsmax(its+1,itv+1)
          isca=int((xtsmax-scaa)/4.0)+1
	  if (isca.le.0) isca=1
          if ((isca+1).lt.nbfi(its+1,itv+1)) then
!             sca1=tsmax(its+1,itv+1)-(isca-1)*4.0
!             sca2=tsmax(its+1,itv+1)-isca*4.0
             sca1=xtsmax-(isca-1)*4.0
             sca2=xtsmax-isca*4.0
             else
	     isca=nbfi(its+1,itv+1)-1
!             sca1=tsmax(its+1,itv+1)-(isca-1)*4.0
             sca1=xtsmax-(isca-1)*4.0
	     sca2=tsmin(its+1,itv+1)
	     endif
!           iindex=indts(its+1)+nbfic(its+1,itv+1)-nbfi(its+1,itv+1)+isca
           iindex=indts(its+1)+nbfic4-nbfi4+isca
	   roinf=rolutt(ib,ip1,iaot1,iindex)
	   rosup=rolutt(ib,ip1,iaot1,iindex+1)
	   ro4=roinf+(rosup-roinf)*(scaa-sca1)/(sca2-sca1)
c       write(6,*) "Point its+1,itv+1 sca,sca1,sca2,ro..",scaa,sca1,sca2,
c     S            ro4,roinf,rosup
     
       t=(tts(its+1)-xts)/(tts(its+1)-tts(its))
       u=(ttv(its,itv+1)-xtv)/(ttv(its,itv+1)-ttv(its,itv))
       roiaot1=ro1*t*u+ro2*u*(1.-t)+ro3*(1.-u)*t+ro4*(1.-u)*(1.-t)
c       write(6,*) "final interpolated value xts,xtv,fi,scaa,ro..",xts,
c     S            xtv,xfi,scaa,ro
c        write(6,*)   "ro at aot=",aot550nm(iaot1), roiaot1    
C Compute for iaot2
        iaot2=iaot1+1	          
C interpolate point 1 (its,itv) vs scattering angle 
       xtsmax=tsmax(its,itv)
       if ((its.ne.1).and.(itv.ne.1)) then
          isca=int((xtsmax-scaa)/4.0)+1
	  if (isca.le.0) isca=1
          if ((isca+1).lt.nbfi(its,itv)) then
!             sca1=tsmax(its,itv)-(isca-1)*4.0
!             sca2=tsmax(its,itv)-isca*4.0
             sca1=xtsmax-(isca-1)*4.0
             sca2=xtsmax-isca*4.0
             else
 	     isca=nbfi(its,itv)-1
!             sca1=tsmax(its,itv)-(isca-1)*4.0
             sca1=xtsmax-(isca-1)*4.0
	     sca2=tsmin(its,itv)
	     endif
!	   iindex=indts(its)+nbfic(its,itv)-nbfi(its,itv)+isca
	   iindex=indts(its)+nbfic1-nbfi1+isca
	   roinf=rolutt(ib,ip1,iaot2,iindex)
	   rosup=rolutt(ib,ip1,iaot2,iindex+1)
	   ro1=roinf+(rosup-roinf)*(scaa-sca1)/(sca2-sca1)
           else
!           sca1=tsmax(its,itv)
!           sca2=tsmax(its,itv)
           sca1=xtsmax
           sca2=xtsmax
!	   iindex=indts(its)+nbfic(its,itv)-nbfi(its,itv)+1
	   iindex=indts(its)+nbfic1-nbfi1+1
	   roinf=rolutt(ib,ip1,iaot2,iindex)
!	   rosup=rolutt(ib,ip1,iaot2,iindex)
!	   ro1=rolutt(ib,ip1,iaot2,iindex)
	   rosup=roinf
	   ro1=roinf
           endif
c       write(6,*) "Point its,itv sca,sca1,sca2,ro..",scaa,sca1,sca2,ro1,
c     S             roinf,rosup
C interpolate point 2 (its+1,itv) vs scattering angle 
       xtsmax=tsmax(its+1,itv)
       if ((itv.ne.1)) then
          isca=int((xtsmax-scaa)/4.0)+1
	  if (isca.le.0) isca=1
          if ((isca+1).lt.nbfi(its+1,itv)) then
!             sca1=tsmax(its+1,itv)-(isca-1)*4.0
!             sca2=tsmax(its+1,itv)-isca*4.0
             sca1=xtsmax-(isca-1)*4.0
             sca2=xtsmax-isca*4.0
             else
 	     isca=nbfi(its+1,itv)-1
!             sca1=tsmax(its+1,itv)-(isca-1)*4.0
             sca1=xtsmax-(isca-1)*4.0
	     sca2=tsmin(its+1,itv)
	     endif
!	   iindex=indts(its+1)+nbfic(its+1,itv)-nbfi(its+1,itv)+isca
	   iindex=indts(its+1)+nbfic2-nbfi2+isca
	   roinf=rolutt(ib,ip1,iaot2,iindex)
	   rosup=rolutt(ib,ip1,iaot2,iindex+1)
	   ro2=roinf+(rosup-roinf)*(scaa-sca1)/(sca2-sca1)
           else
!           sca1=tsmax(its+1,itv)
!           sca2=tsmax(its+1,itv)
           sca1=xtsmax
           sca2=xtsmax
!	   iindex=indts(its+1)+nbfic(its+1,itv)-nbfi(its+1,itv)+1
	   iindex=indts(its+1)+nbfic2-nbfi2+1
	   roinf=rolutt(ib,ip1,iaot2,iindex)
!	   rosup=rolutt(ib,ip1,iaot2,iindex)
!	   ro2=rolutt(ib,ip1,iaot2,iindex)
	   rosup=roinf
	   ro2=roinf
           endif
c       write(6,*) "Point its+1,itv sca,sca1,sca2,ro..",scaa,sca1,sca2,
c     S            ro2,roinf,rosup
C interpolate point 3 (its,itv+1) vs scattering angle 
       xtsmax=tsmax(its,itv+1)
       if ((its.ne.1)) then
          isca=int((xtsmax-scaa)/4.0)+1
	  if (isca.le.0) isca=1
          if ((isca+1).lt.nbfi(its,itv+1)) then
!             sca1=tsmax(its,itv+1)-(isca-1)*4.0
!             sca2=tsmax(its,itv+1)-isca*4.0
             sca1=xtsmax-(isca-1)*4.0
             sca2=xtsmax-isca*4.0
             else
	     isca=nbfi(its,itv+1)-1
!             sca1=tsmax(its,itv+1)-(isca-1)*4.0
             sca1=xtsmax-(isca-1)*4.0
	     sca2=tsmin(its,itv+1)
	     endif
!           iindex=indts(its)+nbfic(its,itv+1)-nbfi(its,itv+1)+isca
           iindex=indts(its)+nbfic3-nbfi3+isca
	   roinf=rolutt(ib,ip1,iaot2,iindex)
	   rosup=rolutt(ib,ip1,iaot2,iindex+1)
	   ro3=roinf+(rosup-roinf)*(scaa-sca1)/(sca2-sca1)
           else
!           sca1=tsmax(its,itv+1)
!           sca2=tsmax(its,itv+1)
           sca1=xtsmax
           sca2=xtsmax
!	   iindex=indts(its)+nbfic(its,itv+1)-nbfi(its,itv+1)+1
	   iindex=indts(its)+nbfic3-nbfi3+1
	   roinf=rolutt(ib,ip1,iaot2,iindex)
!	   rosup=rolutt(ib,ip1,iaot2,iindex)
!	   ro3=rolutt(ib,ip1,iaot2,iindex)
	   rosup=roinf
	   ro3=roinf
           endif
c       write(6,*) "Point its,itv+1 sca,sca1,sca2,ro..",scaa,sca1,sca2,
c     S            ro3,roinf,rosup
C interpolate point 4 (its+1,itv+1) vs scattering angle 
          xtsmax=tsmax(its+1,itv+1)
          isca=int((xtsmax-scaa)/4.0)+1
	  if (isca.le.0) isca=1
          if ((isca+1).lt.nbfi(its+1,itv+1)) then
!             sca1=tsmax(its+1,itv+1)-(isca-1)*4.0
!             sca2=tsmax(its+1,itv+1)-isca*4.0
             sca1=xtsmax-(isca-1)*4.0
             sca2=xtsmax-isca*4.0
             else
	     isca=nbfi(its+1,itv+1)-1
!             sca1=tsmax(its+1,itv+1)-(isca-1)*4.0
             sca1=xtsmax-(isca-1)*4.0
	     sca2=tsmin(its+1,itv+1)
	     endif
!           iindex=indts(its+1)+nbfic(its+1,itv+1)-nbfi(its+1,itv+1)+isca
           iindex=indts(its+1)+nbfic4-nbfi4+isca
	   roinf=rolutt(ib,ip1,iaot2,iindex)
	   rosup=rolutt(ib,ip1,iaot2,iindex+1)
	   ro4=roinf+(rosup-roinf)*(scaa-sca1)/(sca2-sca1)
c       write(6,*) "Point its+1,itv+1 sca,sca1,sca2,ro..",scaa,sca1,sca2,
c    S            ro4,roinf,rosup
     
       t=(tts(its+1)-xts)/(tts(its+1)-tts(its))
       u=(ttv(its,itv+1)-xtv)/(ttv(its,itv+1)-ttv(its,itv))
       roiaot2=ro1*t*u+ro2*u*(1.-t)+ro3*(1.-u)*t+ro4*(1.-u)*(1.-t)
c       write(6,*) "final interpolated value xts,xtv,fi,scaa,ro..",xts,
c     S            xtv,xfi,scaa,ro
c       write(6,*)   "ro at aot=",aot550nm(iaot2), roiaot2 
c interpolation as log of tau       
!       deltaaot=log(aot550nm(iaot2))-log(aot550nm(iaot1))
!       deltaaot=(log(raot550nm)-log(aot550nm(iaot1)))/deltaaot
       deltaaot=logaot550nm(iaot2)-logaot550nm(iaot1)
       deltaaot=(log(raot550nm)-logaot550nm(iaot1))/deltaaot
       ro=roiaot1+(roiaot2-roiaot1)*deltaaot
       rop1=ro
c       write(6,*) " ro at aot= ",raot550nm," at pres ",tpres(ip1),
c     s  "for band ",ib,"="
c     S  ,rop1
C compute for ip2
C Compute for iaot1	          
C interpolate point 1 (its,itv) vs scattering angle 
       xtsmax=tsmax(its,itv)
       if ((its.ne.1).and.(itv.ne.1)) then
          isca=int((xtsmax-scaa)/4.0)+1
	  if (isca.le.0) isca=1
          if ((isca+1).lt.nbfi(its,itv)) then
!             sca1=tsmax(its,itv)-(isca-1)*4.0
!             sca2=tsmax(its,itv)-isca*4.0
             sca1=xtsmax-(isca-1)*4.0
             sca2=xtsmax-isca*4.0
             else
 	     isca=nbfi(its,itv)-1
!             sca1=tsmax(its,itv)-(isca-1)*4.0
             sca1=xtsmax-(isca-1)*4.0
	     sca2=tsmin(its,itv)
	     endif
!	   iindex=indts(its)+nbfic(its,itv)-nbfi(its,itv)+isca
	   iindex=indts(its)+nbfic1-nbfi1+isca
	   roinf=rolutt(ib,ip2,iaot1,iindex)
	   rosup=rolutt(ib,ip2,iaot1,iindex+1)
	   ro1=roinf+(rosup-roinf)*(scaa-sca1)/(sca2-sca1)
           else
!           sca1=tsmax(its,itv)
!           sca2=tsmax(its,itv)
           sca1=xtsmax
           sca2=xtsmax
!	   iindex=indts(its)+nbfic(its,itv)-nbfi(its,itv)+1
	   iindex=indts(its)+nbfic1-nbfi1+1
	   roinf=rolutt(ib,ip2,iaot1,iindex)
!	   rosup=rolutt(ib,ip2,iaot1,iindex)
!	   ro1=rolutt(ib,ip2,iaot1,iindex)
	   rosup=roinf
	   ro1=roinf
           endif
C        write(6,*) "Point its,itv sca,sca1,sca2,ro..",scaa,sca1,sca2,ro1,
C     S             roinf,rosup
C interpolate point 2 (its+1,itv) vs scattering angle 
       xtsmax=tsmax(its+1,itv)
       if ((itv.ne.1)) then
          isca=int((xtsmax-scaa)/4.0)+1
	  if (isca.le.0) isca=1
          if ((isca+1).lt.nbfi(its+1,itv)) then
!             sca1=tsmax(its+1,itv)-(isca-1)*4.0
!             sca2=tsmax(its+1,itv)-isca*4.0
             sca1=xtsmax-(isca-1)*4.0
             sca2=xtsmax-isca*4.0
             else
 	     isca=nbfi(its+1,itv)-1
!             sca1=tsmax(its+1,itv)-(isca-1)*4.0
             sca1=xtsmax-(isca-1)*4.0
	     sca2=tsmin(its+1,itv)
	     endif
!	   iindex=indts(its+1)+nbfic(its+1,itv)-nbfi(its+1,itv)+isca
	   iindex=indts(its+1)+nbfic2-nbfi2+isca
	   roinf=rolutt(ib,ip2,iaot1,iindex)
	   rosup=rolutt(ib,ip2,iaot1,iindex+1)
	   ro2=roinf+(rosup-roinf)*(scaa-sca1)/(sca2-sca1)
           else
!           sca1=tsmax(its+1,itv)
!           sca2=tsmax(its+1,itv)
           sca1=xtsmax
           sca2=xtsmax
!	   iindex=indts(its+1)+nbfic(its+1,itv)-nbfi(its+1,itv)+1
	   iindex=indts(its+1)+nbfic2-nbfi2+1
	   roinf=rolutt(ib,ip2,iaot1,iindex)
!	   rosup=rolutt(ib,ip2,iaot1,iindex)
!	   ro2=rolutt(ib,ip2,iaot1,iindex)
	   rosup=roinf
	   ro2=roinf
           endif
c       write(6,*) "Point its+1,itv sca,sca1,sca2,ro..",scaa,sca1,sca2,
c     S            ro2,roinf,rosup
C interpolate point 3 (its,itv+1) vs scattering angle 
       xtsmax=tsmax(its,itv+1)
       if ((its.ne.1)) then
          isca=int((xtsmax-scaa)/4.0)+1
	  if (isca.le.0) isca=1
          if ((isca+1).lt.nbfi(its,itv+1)) then
!             sca1=tsmax(its,itv+1)-(isca-1)*4.0
!             sca2=tsmax(its,itv+1)-isca*4.0
             sca1=xtsmax-(isca-1)*4.0
             sca2=xtsmax-isca*4.0
             else
	     isca=nbfi(its,itv+1)-1
!             sca1=tsmax(its,itv+1)-(isca-1)*4.0
             sca1=xtsmax-(isca-1)*4.0
	     sca2=tsmin(its,itv+1)
	     endif
!           iindex=indts(its)+nbfic(its,itv+1)-nbfi(its,itv+1)+isca
           iindex=indts(its)+nbfic3-nbfi3+isca
	   roinf=rolutt(ib,ip2,iaot1,iindex)
	   rosup=rolutt(ib,ip2,iaot1,iindex+1)
	   ro3=roinf+(rosup-roinf)*(scaa-sca1)/(sca2-sca1)
           else
!           sca1=tsmax(its,itv+1)
!           sca2=tsmax(its,itv+1)
           sca1=xtsmax
           sca2=xtsmax
!	   iindex=indts(its)+nbfic(its,itv+1)-nbfi(its,itv+1)+1
	   iindex=indts(its)+nbfic3-nbfi3+1
	   roinf=rolutt(ib,ip2,iaot1,iindex)
!	   rosup=rolutt(ib,ip2,iaot1,iindex)
!	   ro3=rolutt(ib,ip2,iaot1,iindex)
	   rosup=roinf
	   ro3=roinf
           endif
c       write(6,*) "Point its,itv+1 sca,sca1,sca2,ro..",scaa,sca1,sca2,
c     S            ro3,roinf,rosup
C interpolate point 4 (its+1,itv+1) vs scattering angle 
          xtsmax=tsmax(its+1,itv+1)
          isca=int((xtsmax-scaa)/4.0)+1
	  if (isca.le.0) isca=1
          if ((isca+1).lt.nbfi(its+1,itv+1)) then
!             sca1=tsmax(its+1,itv+1)-(isca-1)*4.0
!             sca2=tsmax(its+1,itv+1)-isca*4.0
             sca1=xtsmax-(isca-1)*4.0
             sca2=xtsmax-isca*4.0
             else
	     isca=nbfi(its+1,itv+1)-1
!             sca1=tsmax(its+1,itv+1)-(isca-1)*4.0
             sca1=xtsmax-(isca-1)*4.0
	     sca2=tsmin(its+1,itv+1)
	     endif
!           iindex=indts(its+1)+nbfic(its+1,itv+1)-nbfi(its+1,itv+1)+isca
           iindex=indts(its+1)+nbfic4-nbfi4+isca
	   roinf=rolutt(ib,ip2,iaot1,iindex)
	   rosup=rolutt(ib,ip2,iaot1,iindex+1)
	   ro4=roinf+(rosup-roinf)*(scaa-sca1)/(sca2-sca1)
c       write(6,*) "Point its+1,itv+1 sca,sca1,sca2,ro..",scaa,sca1,sca2,
c     S            ro4,roinf,rosup
     
       t=(tts(its+1)-xts)/(tts(its+1)-tts(its))
       u=(ttv(its,itv+1)-xtv)/(ttv(its,itv+1)-ttv(its,itv))
       roiaot1=ro1*t*u+ro2*u*(1.-t)+ro3*(1.-u)*t+ro4*(1.-u)*(1.-t)
c       write(6,*) "final interpolated value xts,xtv,fi,scaa,ro..",xts,
c     S            xtv,xfi,scaa,ro
c        write(6,*)   "ro at aot=",aot550nm(iaot1), roiaot1    
C Compute for iaot2
        iaot2=iaot1+1	          
C interpolate point 1 (its,itv) vs scattering angle 
       xtsmax=tsmax(its,itv)
       if ((its.ne.1).and.(itv.ne.1)) then
          isca=int((xtsmax-scaa)/4.0)+1
	  if (isca.le.0) isca=1
          if ((isca+1).lt.nbfi(its,itv)) then
!             sca1=tsmax(its,itv)-(isca-1)*4.0
!             sca2=tsmax(its,itv)-isca*4.0
             sca1=xtsmax-(isca-1)*4.0
             sca2=xtsmax-isca*4.0
             else
 	     isca=nbfi(its,itv)-1
!             sca1=tsmax(its,itv)-(isca-1)*4.0
             sca1=xtsmax-(isca-1)*4.0
	     sca2=tsmin(its,itv)
	     endif
!	   iindex=indts(its)+nbfic(its,itv)-nbfi(its,itv)+isca
	   iindex=indts(its)+nbfic1-nbfi1+isca
	   roinf=rolutt(ib,ip2,iaot2,iindex)
	   rosup=rolutt(ib,ip2,iaot2,iindex+1)
	   ro1=roinf+(rosup-roinf)*(scaa-sca1)/(sca2-sca1)
           else
!           sca1=tsmax(its,itv)
!           sca2=tsmax(its,itv)
           sca1=xtsmax
           sca2=xtsmax
!	   iindex=indts(its)+nbfic(its,itv)-nbfi(its,itv)+1
	   iindex=indts(its)+nbfic1-nbfi1+1
	   roinf=rolutt(ib,ip2,iaot2,iindex)
!	   rosup=rolutt(ib,ip2,iaot2,iindex)
!	   ro1=rolutt(ib,ip2,iaot2,iindex)
	   rosup=roinf
	   ro1=roinf
           endif
c       write(6,*) "Point its,itv sca,sca1,sca2,ro..",scaa,sca1,sca2,ro1,
c     S             roinf,rosup
C interpolate point 2 (its+1,itv) vs scattering angle 
       xtsmax=tsmax(its+1,itv)
       if ((itv.ne.1)) then
          isca=int((xtsmax-scaa)/4.0)+1
	  if (isca.le.0) isca=1
          if ((isca+1).lt.nbfi(its+1,itv)) then
!             sca1=tsmax(its+1,itv)-(isca-1)*4.0
!             sca2=tsmax(its+1,itv)-isca*4.0
             sca1=xtsmax-(isca-1)*4.0
             sca2=xtsmax-isca*4.0
             else
 	     isca=nbfi(its+1,itv)-1
!             sca1=tsmax(its+1,itv)-(isca-1)*4.0
             sca1=xtsmax-(isca-1)*4.0
	     sca2=tsmin(its+1,itv)
	     endif
!	   iindex=indts(its+1)+nbfic(its+1,itv)-nbfi(its+1,itv)+isca
	   iindex=indts(its+1)+nbfic2-nbfi2+isca
	   roinf=rolutt(ib,ip2,iaot2,iindex)
	   rosup=rolutt(ib,ip2,iaot2,iindex+1)
	   ro2=roinf+(rosup-roinf)*(scaa-sca1)/(sca2-sca1)
           else
!           sca1=tsmax(its+1,itv)
!           sca2=tsmax(its+1,itv)
           sca1=xtsmax
           sca2=xtsmax
!	   iindex=indts(its+1)+nbfic(its+1,itv)-nbfi(its+1,itv)+1
	   iindex=indts(its+1)+nbfic2-nbfi2+1
	   roinf=rolutt(ib,ip2,iaot2,iindex)
!	   rosup=rolutt(ib,ip2,iaot2,iindex)
!	   ro2=rolutt(ib,ip2,iaot2,iindex)
	   rosup=roinf
	   ro2=roinf
           endif
c       write(6,*) "Point its+1,itv sca,sca1,sca2,ro..",scaa,sca1,sca2,
c     S            ro2,roinf,rosup
C interpolate point 3 (its,itv+1) vs scattering angle 
       xtsmax=tsmax(its,itv+1)
       if ((its.ne.1)) then
          isca=int((xtsmax-scaa)/4.0)+1
	  if (isca.le.0) isca=1
          if ((isca+1).lt.nbfi(its,itv+1)) then
!             sca1=tsmax(its,itv+1)-(isca-1)*4.0
!             sca2=tsmax(its,itv+1)-isca*4.0
             sca1=xtsmax-(isca-1)*4.0
             sca2=xtsmax-isca*4.0
             else
	     isca=nbfi(its,itv+1)-1
!             sca1=tsmax(its,itv+1)-(isca-1)*4.0
             sca1=xtsmax-(isca-1)*4.0
	     sca2=tsmin(its,itv+1)
	     endif
!           iindex=indts(its)+nbfic(its,itv+1)-nbfi(its,itv+1)+isca
           iindex=indts(its)+nbfic3-nbfi3+isca
	   roinf=rolutt(ib,ip2,iaot2,iindex)
	   rosup=rolutt(ib,ip2,iaot2,iindex+1)
	   ro3=roinf+(rosup-roinf)*(scaa-sca1)/(sca2-sca1)
           else
!           sca1=tsmax(its,itv+1)
!           sca2=tsmax(its,itv+1)
           sca1=xtsmax
           sca2=xtsmax
!	   iindex=indts(its)+nbfic(its,itv+1)-nbfi(its,itv+1)+1
	   iindex=indts(its)+nbfic3-nbfi3+1
	   roinf=rolutt(ib,ip2,iaot2,iindex)
!	   rosup=rolutt(ib,ip2,iaot2,iindex)
!	   ro3=rolutt(ib,ip2,iaot2,iindex)
	   rosup=roinf
	   ro3=roinf
           endif
c       write(6,*) "Point its,itv+1 sca,sca1,sca2,ro..",scaa,sca1,sca2,
c     S            ro3,roinf,rosup
C interpolate point 4 (its+1,itv+1) vs scattering angle 
          xtsmax=tsmax(its+1,itv+1)
          isca=int((xtsmax-scaa)/4.0)+1
	  if (isca.le.0) isca=1
          if ((isca+1).lt.nbfi(its+1,itv+1)) then
!             sca1=tsmax(its+1,itv+1)-(isca-1)*4.0
!             sca2=tsmax(its+1,itv+1)-isca*4.0
             sca1=xtsmax-(isca-1)*4.0
             sca2=xtsmax-isca*4.0
             else
	     isca=nbfi(its+1,itv+1)-1
!             sca1=tsmax(its+1,itv+1)-(isca-1)*4.0
             sca1=xtsmax-(isca-1)*4.0
	     sca2=tsmin(its+1,itv+1)
	     endif
!           iindex=indts(its+1)+nbfic(its+1,itv+1)-nbfi(its+1,itv+1)+isca
           iindex=indts(its+1)+nbfic4-nbfi4+isca
	   roinf=rolutt(ib,ip2,iaot2,iindex)
	   rosup=rolutt(ib,ip2,iaot2,iindex+1)
	   ro4=roinf+(rosup-roinf)*(scaa-sca1)/(sca2-sca1)
c       write(6,*) "Point its+1,itv+1 sca,sca1,sca2,ro..",scaa,sca1,sca2,
c    S            ro4,roinf,rosup
     
       t=(tts(its+1)-xts)/(tts(its+1)-tts(its))
       u=(ttv(its,itv+1)-xtv)/(ttv(its,itv+1)-ttv(its,itv))
       roiaot2=ro1*t*u+ro2*u*(1.-t)+ro3*(1.-u)*t+ro4*(1.-u)*(1.-t)
c       write(6,*) "final interpolated value xts,xtv,fi,scaa,ro..",xts,
c     S            xtv,xfi,scaa,ro
c       write(6,*)   "ro at aot=",aot550nm(iaot2), roiaot2 
c interpolation as log of tau       
!       deltaaot=log(aot550nm(iaot2))-log(aot550nm(iaot1))
!       deltaaot=(log(raot550nm)-log(aot550nm(iaot1)))/deltaaot
       deltaaot=logaot550nm(iaot2)-logaot550nm(iaot1)
       deltaaot=(log(raot550nm)-logaot550nm(iaot1))/deltaaot
       ro=roiaot1+(roiaot2-roiaot1)*deltaaot
       rop2=ro
c       write(6,*) " ro at aot= ",raot550nm, "at pres= ",tpres(ip2),
c     s  "for band ",ib,"="
c     S  ,rop2
     
       dpres=(pres-tpres(ip1))/(tpres(ip2)-tpres(ip1))
       roatm=rop1+(rop2-rop1)*dpres
       
       return
       end
	
       subroutine readluts(CAMOD,iendarg,tauray,oztransa,wvtransa,
     s                      wvtransb,wvtransc,ogtransa0,ogtransa1,
     s                      ogtransb0,ogtransb1,ogtransc0,ogtransc1,
     s		            tsmax,tsmin,ttv,tts,nbfi,nbfic,indts,
     s                      rolutt,transt,sphalbt,normext,sbandname,err_msg,
     s                      retval,tauraynm,gscoefnm,anglehdf,intrefnm,
     s                      transmnm, spheranm)
C Arguments....     
       parameter (fac = 0.017453293)
       parameter(DFACC_READ = 1)
       character*80 CAMOD
       integer iendarg
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
       character*6 sbandname(16)
       character*80 err_msg
       integer retval
C The following arguments are all names of the LUTs to look up.
       character*256 tauraynm
       character*256 gscoefnm
       character*256 anglehdf
       character*256 intrefnm
       character*256 transmnm
       character*256 spheranm

C local variables....     
       character*6 bandid
       integer sfstart, sfselect, sfrdata, sfendacc, sfend   
       integer sd_id, sds_id, sds_index, status
       integer start(3), edges(3), stride(3)  
       character*10 sbd2
       character*80 fname
       real rolut(7,22,8000)
       real trans(7,22,22),ttsr(22)
       real sphalb(7,22),xx,yy

c Read the gaseous transmission and other constants  
C molecular optical thickness
C       open(1,file="tauray-modis.ASC")
c       write(6,*) "TAURAYNM ",tauraynm
       open(1,file=tauraynm)
       do i=1,8
       read(1,101) bandid,tauray(i)
c       write(6,*) bandid,tauray(i)
       
       if (bandid.ne.sbandname(i)) then
       write(err_msg,*) "Error on rayleigh o.d file,"//
C     s     " check tauray-modis.ASC"
     s     " check ", tauraynm//char(0)
       retval=1
 101   Format(A6,1X,F10.5)       
       return
       endif
       enddo
       close(1)
C gaseous transmission coefficient
C       Open(1,file="gascoef-modis.ASC")
       Open(1,file=gscoefnm)
       Read(1,*) (oztransa(i),i=1,8)
       Read(1,*) (wvtransa(i),i=1,8)
       Read(1,*) (wvtransb(i),i=1,8)
       Read(1,*) (ogtransa1(i),i=1,8)
       Read(1,*) (ogtransb0(i),i=1,8)
       Read(1,*) (ogtransb1(i),i=1,8)
       close(1)
       
C       Write(6,*) "Ozone coeff"
C       Write(6,102) (oztransa(i),i=1,16)
C       Write(6,*) "Water vapor coeff"
C       Write(6,102) (wvtransa(i),i=1,16)
C       Write(6,102) (wvtransb(i),i=1,16)
C       Write(6,102) (wvtransc(i),i=1,16)
C       Write(6,*) "other gases coeff"
C       Write(6,102) (ogtransa0(i),i=1,16)
C       Write(6,102) (ogtransa1(i),i=1,16)
C       Write(6,102) (ogtransb0(i),i=1,16)
C       Write(6,102) (ogtransb1(i),i=1,16)
C       Write(6,102) (ogtransc0(i),i=1,16)
C       Write(6,102) (ogtransc1(i),i=1,16)
 102   FORMAT(16(E12.4,1X))      
       

       xtsmin=0
       xtsstep=4.0
       do i=1,22
       tts(i)=xtsmin+xtsstep*(i-1)
       do j=1,20
       nbfic(i,j)=0
       enddo
       enddo
       xtvmin=2.84090
       xtvstep=(6.52107-2.84090)
C       pi=acos(0.0)*2.0
C
C       sd_id = sfstart("ANGLE_NEW.hdf", DFACC_READ)
C       write(6,*) " here ", anglehdf,DFACC_READ
       sd_id = sfstart(anglehdf, DFACC_READ)
       if (sd_id.eq.-1) then
       write(err_msg,*) "check ANGLE.hdf file ",anglehdf//char(0)
       retval=2
       return
       endif
C
       start(1) = 0
       start(2) = 0
       edges(1) = 22
       edges(2) = 20
       stride(1) = 1
       stride(2) = 1
       sds_index = sfn2index(sd_id, "TSMAX")
       sds_index = 0
       sds_id    = sfselect(sd_id, sds_index)
       status = sfrdata(sds_id, start, stride, edges, tsmax)
       if (status .eq. -1) then
       write(err_msg,*) "check ANGLE.hdf file ",anglehdf//char(0)
       retval=3
       return
       endif
       status = sfendacc(sds_id)
       sds_index = sfn2index(sd_id, "TSMIN")
       sds_index = 1
       sds_id    = sfselect(sd_id, sds_index)
       status = sfrdata(sds_id, start, stride, edges, tsmin)
       status = sfendacc(sds_id)
       sds_index = sfn2index(sd_id, "TTV")
       sds_index = 2
       sds_id    = sfselect(sd_id, sds_index)
       status = sfrdata(sds_id, start, stride, edges, ttv)
       status = sfendacc(sds_id)
       sds_index = sfn2index(sd_id, "NBFI")
       sds_index = 3
       sds_id    = sfselect(sd_id, sds_index)
       status = sfrdata(sds_id, start, stride, edges, nbfi)
       status = sfendacc(sds_id)
       sds_index = sfn2index(sd_id, "NBFIC")
       sds_index = 4
       sds_id    = sfselect(sd_id, sds_index)
       status = sfrdata(sds_id, start, stride, edges, nbfic)
       status = sfendacc(sds_id)
       
C read the one dimension SDS's
              
       sds_index = sfn2index(sd_id, "INDTS")
       sds_index = 5
       sds_id    = sfselect(sd_id, sds_index)
       status = sfrdata(sds_id, start, stride, edges, indts)
       status = sfendacc(sds_id)
       status = sfend(sd_id)  
       sds_index = sfn2index(sd_id, "TTS")
       sds_index = 6
       sds_id    = sfselect(sd_id, sds_index)
       status = sfrdata(sds_id, start, stride, edges, tts)
       status = sfendacc(sds_id)
       
C	write (6,*) 'ANGLES are OK'
C  read LUT from a separate file
C       open(1,file=PAR1,access="direct",recl=100)
C 113   format(4(F10.7,1X),A80)
C       do iaot=1,22   
C       read(1,rec=1+(iaot-1)*322) str2(1:100)
C       read(1,rec=2+(iaot-1)*322) str2(101:124)
C       read(str2,113) lmin,lmax,taul,tau550,amod
C       write(6,*) lmin,lmax,taul,tau550
C       write(6,'(A80)') amod
C       do i=3,322
C       read(1,rec=i+(iaot-1)*322) (rolut(iaot,k),k=(i-2)*25+1,(i-1)*25)
C       enddo
C       enddo
C       close(1)
C
C Begin read all look up table  (intrinsic reflectance)     
       sd_id = sfstart(intrefnm, DFACC_READ)
       do iband=1,8
C       write(sbd2,'(A6)') sbandname(iband)
C       jj=index(sbd2," ")-1
C       write(fname,*) "reslutV2.0-",sbd2(1:jj),
C     &   "-",CAMOD(1:iendarg),".hdf"
C       write(6,*) fname(2:80)
C       sd_id = sfstart(fname(2:80), DFACC_READ)
C       sd_id = sfstart(intrefnm(iband), DFACC_READ)
       start(1) = 0
       start(2) = 0
       start(3) = 0
       edges(1) = 7
       edges(2) = 22
       edges(3) = 8000
       stride(1) = 1
       stride(2) = 1
       stride(3) = 1
C
C   Get sds name...
C       
       if (iband .lt. 10) then
          write(fname,*) "NRLUT_BAND_"//char(iband+ichar('0'))
     &         //char(0)
       else
          write(fname,*) "NRLUT_BAND_1"//char(iband-10+ichar('0'))
     &         //char(0)
       endif
C  But, passing fname through sfn2index doesn't work.      
C       sds_index = sfn2index(sd_id, fname(2:14))
C       write(6,*)fname, sds_index
C     for now just do this, since the SDSs are in order...
       sds_index = iband-1
       sds_id    = sfselect(sd_id, sds_index)
       status = sfrdata(sds_id, start, stride, edges, rolut)
       status = sfendacc(sds_id)
C       status = sfend(sd_id)
       do ipres=1,7  
       do itau=1,22
       do ival=1,8000
       rolutt(iband,ipres,itau,ival)=rolut(ipres,itau,ival)
       enddo
       enddo
       enddo
       enddo
       status = sfend(sd_id)
C       write(6,*) "LOOKUP DEB",rolutt(8,1,1,1),rolutt(8,2,1,1)
C End read all look up table     (intrinsic reflectance)


C Begin read all look up table  (transmission)     
       open(1,file=transmnm)
       do iband=1,8
C       write(sbd2,'(A6)') sbandname(iband)
C       jj=index(sbd2," ")-1
C       write(fname,*) "TRANSLUTV2.0-",sbd2(1:jj),
C     &   "-",CAMOD(1:iendarg),".ASCII"
c       write(6,*) fname
C       open(1,file=fname(2:80))
C
C This first read contains information about the source of the data 
C  ignore for now
       read(1,*)
       do ipres=1,7
C
C This next read contains information about the pressure level of the data
C  ignore for now
       read(1,*)
       do i=1,21
       read(1,*) ttsr(i),(transt(iband,ipres,iaot,i),iaot=1,22)
       if (abs(tts(i)-ttsr(i)).gt.1.E-5) then
          write(err_msg,*) "problem with transmission LUT"//char(0)
          retval=2
          return
	  endif
       enddo
       enddo
C        do ipres=1,7
C       do itau=1,22
C       do ival=1,21
C       transt(iband,ipres,itau,ival)=trans(ipres,itau,ival)
C       enddo
C       enddo
C       enddo
       enddo
      close(1)
C End read all look up table  (transmission)     


C Begin read all look up table  (spherical albedo)     
       open(1,file=spheranm)
       do iband=1,8
C       write(sbd2,'(A6)') sbandname(iband)
C       jj=index(sbd2," ")-1
C       write(fname,*) "AEROLUTV2.0-",sbd2(1:jj),
C     &   "-",CAMOD(1:iendarg),".ASCII"
c       write(6,*) fname
C       open(1,file=fname(2:80))
C       open(1,file=spheranm(iband))
C This first read contains information about the source of the data 
C  ignore for now
       read(1,*)
       do ipres=1,7
C
C This next read contains information about the pressure level of the data
C  ignore for now
       read(1,*)
       do iaot=1,22
       read(1,*) xx,sphalbt(iband,ipres,iaot),normext(iband,ipres,iaot)
       enddo
       enddo
C       close(1)
C       do ipres=1,7
C       do itau=1,22
C       sphalbt(iband,ipres,itau)=sphalb(ipres,itau)
C       enddo
C       enddo
       enddo
       close(1)
C End read all look up table  (spherical albedo) 
    
       retval=0
       return
       end
