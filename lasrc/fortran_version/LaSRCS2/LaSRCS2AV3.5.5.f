        program LaSRCS2A
        
        implicit none
        integer(2), allocatable :: toaband(:,:,:)
        integer(2), allocatable :: ttoaband(:,:)
        integer(2), allocatable :: tempband1(:,:)
        integer(2), allocatable :: tempband2(:,:)
        integer(2), allocatable :: sband(:,:,:)
        integer(2), allocatable :: tband(:,:,:)
        integer(2), allocatable :: oband(:,:)
        integer(2), allocatable :: aotband(:,:)
        integer(2), allocatable :: resband(:,:)
        integer(2), allocatable :: modband(:,:)
        integer(2), allocatable :: wvband(:,:)
        integer(2), allocatable :: wv(:,:)
        real, allocatable :: tlat(:,:)
        real, allocatable :: tlon(:,:)
        real, allocatable :: twv(:,:)
        real, allocatable :: twvi(:,:)
        real, allocatable :: tozi(:,:)
        real, allocatable :: txcmg(:,:)
        real, allocatable :: tycmg(:,:)
        real, allocatable :: tp(:,:)
        real, allocatable :: taero(:,:)
        real, allocatable :: taeros(:,:)
        real, allocatable :: taero3(:,:)
        real, allocatable :: teps(:,:)
        real, allocatable :: tepss(:,:)
        real, allocatable :: tresi(:,:)
        real, allocatable :: slprb1(:,:)
        real, allocatable :: intrb1(:,:)
        real, allocatable :: slprb2(:,:)
        real, allocatable :: intrb2(:,:)
        real, allocatable :: slprb7(:,:)
        real, allocatable :: intrb7(:,:)



        BYTE, allocatable :: oz(:,:)
        BYTE, allocatable :: cloud(:,:)
        BYTE, allocatable :: ipflag(:,:)
        BYTE, allocatable :: aerimp(:,:)
        BYTE, allocatable :: smflag(:,:)
        integer(2), allocatable :: dem(:,:)
        integer(2), allocatable :: ratiob1(:,:)
        integer(2), allocatable :: ratiob2(:,:)
        integer(2), allocatable :: ratiob7(:,:)
        integer (2), allocatable :: tratiob1(:,:)
        integer (2), allocatable :: tratiob2(:,:)
        integer (2), allocatable :: tratiob7(:,:)
        integer (2), allocatable :: tnit(:,:)
        integer(2), allocatable :: intratiob1(:,:)
        integer(2), allocatable :: intratiob2(:,:)
        integer(2), allocatable :: intratiob7(:,:)
        integer(2), allocatable :: slpratiob1(:,:)
        integer(2), allocatable :: slpratiob2(:,:)
        integer(2), allocatable :: slpratiob7(:,:)
        integer(2), allocatable :: andwi(:,:)
        integer(2), allocatable :: sndwi(:,:)
        real th1,th2
C        character*10 ipa1
        integer icor,jdoy
        real om,dsol

        character(400) filename,fname,filenameanc
        character(2) suffix(13)
        data suffix 
     &/"01","02","03","04","05","06","07","08","8a","09","10","11","12"/
        integer ii,nr,nc,ib,als,i,j,k,ierr
        integer nrcmg,nccmg,icmg,jcmg,i20,j20,i60,j60
        real u,v
        real xcmg,ycmg
        real xndwi
        real xts,xfs,cpi,xtv,xfv,xfi
        integer iband
        integer l
        
! BEGIN OF HDF PARAMETER BLOCK        
         integer sfsattr
         integer sfstart, sfselect, sfrdata, sfendacc, sfend,set_qamap
         integer sfwdata,sfcreate
         integer sd_id, sds_id, sds_index, status
         integer start(5), edges(5), stride(5)
         integer DFACC_READ,DFACC_RDWR,DFNT_CHAR8,DFACC_CREATE
         integer DFNT_INT16,DFNT_FLOAT32,DFNT_UINT8
         parameter (DFACC_READ = 1)
         parameter (DFACC_RDWR = 3)
         parameter (DFACC_CREATE = 4)
         parameter (DFNT_INT16 = 22)
         parameter (DFNT_CHAR8 = 4)
         parameter (DFNT_FLOAT32 = 5)
         parameter (DFNT_UINT8 = 21)
         character*80 sds_name
         integer rank
         integer dim_sizes(5)
         integer dim_length(5), comp_type, comp_prm(4)
        real scalefactor,offset
! END OF HDF PARAMETER BLOCK         

!begin block atmospheric correction variables
c   Aerosol model (string) and length of the string
       character*80 CAMOD
       integer iendarg
C Look up table for atmospheric and geometric quantities
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
       real xtsmin,xtsstep,xtvmin,xtvstep,pi
       real aot550nm(22)
       data aot550nm /0.01,0.05,0.10,0.15,0.20,0.30,0.40,0.60,
     &                0.80,1.00,1.20,1.40,1.60,1.80,2.00,
     &                2.30,2.60,3.00,3.50,4.00,4.50,5.00/
       real tpres(7)
       data tpres/1050.,1013.0,900.0,800.0,700.0,600.0,500.0/
c usefull variables       
       character*6 sbandname(16)
       character*80 err_msg
C The following arguments are all names of the LUTs to look up.
       character*256 tauraynm
       character*256 gscoefnm
       character*256 anglehdf
       character*256 intrefnm
       character*256 transmnm
       character*256 spheranm
        real rotoa,raot550nm,uoz,uwv,pres
!        integer ib
        integer retval
        real roslamb
       data (sbandname(i),i=1,16)/"modb1","modb2","modb3","modb4","modb5","modb6",
     s  "modb7","modb8","modb9","modb10","modb11","modb12","modb13","modb14",
     s  "modb15","modb16"/

        real tgo,roatm,ttatmg,satm,xrorayp,next
c        integer cavisind(11)
c        data cavisind /1,2,3,4,4,5,5,5,6,6,7/ 
!end         block atmospheric correction variables
!block sharpening
!end block sharpening
       real x0,y0,gsize,lat,lon
       integer row0,col0,utmzone,row,col
       real pres11,pres12,pres21,pres22
       integer uoz11,uoz21,uoz12,uoz22

!water vapor retrieval band6/band5 ratio coefficient
       real cwv(2)
       data cwv /6.31097,13.5265/
       real fac,xmus
!aerosol retrieval 
       real erelc(16),troatm(16),raot,residual
       real tbtgo(16,22),tbroatm(16,22),tbttatmg(16,22),tbsatm(16,22)
       real corfb
       
       integer iverbose,iband1,iband3,nit   
!calibration
       integer fr,lr,fc,lc
       integer nri,nci,nr20i,nc20i,nr60i,nc60i
        integer xd
        real rb1,rb2
        real slpr11,slpr12,slpr21,slpr22
        real intr11,intr12,intr21,intr22
        real corf,ros4,ros5
        character*400 pfileout,fileout,fileout1,fileout2
        character*26 sstr
        integer iout,iopt
        real eps1,eps2,eps3,residual1,residual2,residual3,epsmin,eps
        real coefc,resepsmin
        real sraot1,sraot2,sraot3
        real rosband1
        integer iaots,ifast
        real xa,xb,xc,xe,xf,coefa,coefb
        real taeroavg,tepsavg
        integer nbaeroavg
        integer nbpixnf,ipass,nbpixtot,wind
        real rsurf




!initialisation for look up table
       iverbose=0
c       call getarg(1,ipa1)
c       read(ipa1,*) icor
       icor=1
       write(6,*) "correction",icor


!initialisation of variable
       fac=atan(1.)*4./180.
!initialisation erecl
!initialisation for look up table
       xtv=0.
       xfi=0.
       xtsmin=0
       xtsstep=4.0
       xtvmin=2.84090
       xtvstep=(6.52107-2.84090)
        CAMOD="URBANCLEAN-V3.0"
        iendarg=10
        tauraynm="MSILUT/tauray-msi.ASC"
        gscoefnm="MSILUT/gascoef-msi.ASC"
        anglehdf="MSILUT/ANGLE_NEW.hdf"
        intrefnm="MSILUT/RES_LUT_V3.0-URBANCLEAN-V3.0.hdf"
        transmnm="MSILUT/TRANS_LUT_V3.0-URBANCLEAN-V3.0.ASCII"
        spheranm="MSILUT/AERO_LUT_V3.0-URBANCLEAN-V3.0.ASCII"
        call readlutsmsi(CAMOD,iendarg,tauray,oztransa,wvtransa,
     s                      wvtransb,wvtransc,ogtransa0,ogtransa1,
     s                      ogtransb0,ogtransb1,ogtransc0,ogtransc1,
     s                            tsmax,tsmin,ttv,tts,nbfi,nbfic,indts,
     s                      rolutt,transt,sphalbt,normext,sbandname,err_msg,
     s                      retval,tauraynm,gscoefnm,anglehdf,intrefnm,
     s                      transmnm, spheranm)
        write(6,*) "the luts for urban clean case v3.0 have been read"
        write(6,*) "we can now perform atmospheric correction"


        pi=atan(1.)*4.
        cpi=pi
        read(5,*) iopt
        read(5,'(A400)') filename
        read(5,'(A400)') filenameanc
        read(5,*) nri,nci,nr20i,nc20i,nr60i,nc60i
        read(5,*) xts,xfs,xtv,xfv
        read(5,*) utmzone,row0,col0,y0,x0,gsize  
        write(6,*) "y0,x0 ",y0,x0
c        read(5,*) (xcal(i),i=1,8)
c        write(6,*) (xcal(i),i=1,8)
c        read(5,*) (ebw(i),i=1,8)
c        write(6,*) (ebw(i),i=1,8)
c        read(5,*) (es(i),i=1,8)
c        write(6,*) (es(i),i=1,8)
        read(5,*) jdoy
        write(6,*) jdoy
        read(5,'(A400)') pfileout
c        read(5,*) fr,lr,fc,lc
        read(5,*,end=19) fr,lr,fc,lc
 19     continue
        if (fr.eq.0) then
        write(fileout,'(A400)') pfileout
        fr=1
        lr=nri
        fc=1
        lc=nci
        iout=0
        ii=index(pfileout," ")-5
        write(fileout1,*) pfileout(1:ii),"-file1.hdf"
        write(fileout2,*) pfileout(1:ii),"-file2.hdf"
        else
        ii=index(pfileout," ")-5
        write(sstr,'(A2,4(A1,I5.5))') "ss","-",fr,"-",lr,"-",fc,"-",lc
        write(fileout1,*) pfileout(1:ii),sstr,"-file1.hdf"
        write(fileout2,*) pfileout(1:ii),sstr,"-file2.hdf"
        iout=2
        endif
c        read(5,*) fpln,lpln,fpcl,lpcl
c        endif
c        ii=index(pfileout," ")-1
c        jj=index(pfileout," ",.true.)
c        write(6,*) fileout(jj:ii),ii,jj,fileout
c        stop
         pres=1013.0
        uoz=0.30
        uwv=0.5
c        gsize=30.
        raot550nm=0.06
        fac=pi/180.
        xmus=cos(xts*fac) 
        om=(.9856*float(jdoy-4))*pi/180.
        dsol=1./((1.-.01673*cos(om))**2)
         
c read the DEM
        nccmg=7200
        nrcmg=3600
       allocate (dem(nccmg,nrcmg),stat=ierr)
       sd_id= sfstart("CMGDEM.hdf",DFACC_READ)
       start(1)=0
       start(2) = 0
       edges(1) = nccmg
       edges(2) = nrcmg
       stride(1) = 1
       stride(2) = 1
       sds_index = 0
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,dem)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
c close HDF file
       status = sfend(sd_id)
       write(6,*) "status sfend ",status
       
       write(6,*) "DEM READ ", dem(2001,1001)
       
c read the RATIOFILE
        nccmg=7200
        nrcmg=3600
       allocate (ratiob1(nccmg,nrcmg),stat=ierr)
       allocate (ratiob2(nccmg,nrcmg),stat=ierr)
       allocate (ratiob7(nccmg,nrcmg),stat=ierr)
       allocate (intratiob1(nccmg,nrcmg),stat=ierr)
       allocate (intratiob2(nccmg,nrcmg),stat=ierr)
       allocate (intratiob7(nccmg,nrcmg),stat=ierr)
       allocate (slpratiob1(nccmg,nrcmg),stat=ierr)
       allocate (slpratiob2(nccmg,nrcmg),stat=ierr)
       allocate (slpratiob7(nccmg,nrcmg),stat=ierr)
       allocate (andwi(nccmg,nrcmg),stat=ierr)
       allocate (sndwi(nccmg,nrcmg),stat=ierr)
      sd_id= sfstart("ratiomapndwiexp.hdf",DFACC_READ)
      start(1)=0
       start(2) = 0
       edges(1) = nccmg
       edges(2) = nrcmg
       stride(1) = 1
       stride(2) = 1
       

       sds_index = 3
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,ratiob2)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       
       sds_index = 2
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,ratiob1)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       
       sds_index = 4
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,ratiob7)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
      
       sds_index = 14
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,sndwi)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       
       sds_index = 21
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,slpratiob1)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       sds_index = 22
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,intratiob1)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status








       sds_index = 24
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,slpratiob2)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       sds_index = 25
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,intratiob2)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       
       
       
       
      
       sds_index = 27
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,slpratiob7)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       sds_index = 28
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,intratiob7)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       
       sds_index = 6
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,andwi)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       
c close HDF file
       status = sfend(sd_id)
       write(6,*) "status sfend ",status
c do not compute ratio based on averaged ndwi (for now!)    
c       do i=1,nccmg
c       do j=1,nrcmg
c       ratiob1(i,j)=int(dble(andwi(i,j)*dble(slpratiob1(i,j)/1000.))+intratiob1(i,j))  
c       ratiob2(i,j)=int(dble(andwi(i,j)*dble(slpratiob2(i,j)/1000.))+intratiob2(i,j))  
c       ratiob7(i,j)=int(dble(andwi(i,j)*dble(slpratiob7(i,j)/1000.))+intratiob7(i,j))
c       enddo
c       enddo  
       
       write(6,*) "RATIO READ ", ratiob1(2001,1001),ratiob2(2001,1001),ratiob7(2001,1001)
       write(6,*) "RATIO READ ", slpratiob1(2001,1001),slpratiob2(2001,1001),slpratiob7(2001,1001)
       write(6,*) "RATIO READ ", intratiob1(2001,1001),intratiob2(2001,1001),intratiob7(2001,1001)
       write(6,*) "RATIO READ ", andwi(2001,1001)
c       stop
        


        
c read ozone and water vapor ancillary        (begin)
        nrcmg=3600
        nccmg=7200
       allocate (oz(nccmg,nrcmg),stat=ierr)
       allocate (wv(nccmg,nrcmg),stat=ierr)
       ii=index(filenameanc," ")-1
        fname=filenameanc(1:ii)
       sd_id= sfstart(fname,DFACC_READ)
      start(1)=0
       start(2) = 0
       edges(1) = nccmg
       edges(2) = nrcmg
       stride(1) = 1
       stride(2) = 1
       sds_index = 0
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,oz)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
c ESPA SDS index is 1 for water vapor in USGS EROS ozone products
cORIG  sds_index = 2
       sds_index = 1
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,wv)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
c close HDF file
       status = sfend(sd_id)
       write(6,*) "status sfend ",status
       
       write(6,*) "Ozone Water vapor read"
c read ozone and water vapor ancillary (end)
        nr=lr-fr+1
        nc=lc-fc+1
        
       allocate(toaband(13,nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(ttoaband(nci,nri),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(tempband1(nc20i,nr20i),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(tempband2(nc60i,nr60i),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(sband(13,nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(tband(13,nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(oband(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(aotband(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(resband(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(tresi(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(taero(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(taeros(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(tepss(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(taero3(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(teps(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(modband(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(wvband(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(tlat(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(tlon(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(txcmg(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(tycmg(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(twv(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(twvi(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(tozi(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(tp(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(tratiob1(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(tratiob2(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(tratiob7(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(tnit(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(cloud(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(ipflag(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(aerimp(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(smflag(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(slprb1(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(intrb1(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(slprb2(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(intrb2(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(slprb7(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"
       allocate(intrb7(nc,nr),STAT=als)
       if (als /= 0) stop "*** Not enough memory ***"

! Getting parameter for atmospheric correction        
!        write(6,*) " aot550nm,pressure [Millibars] ,uoz [cm.atm],uwv [g/cm2]"
         write(6,*) raot550nm,pres,uoz,uwv
c reading the hdf input
      ii=index(filename," ")-1
        fname=filename(1:ii)
       sd_id= sfstart(fname,DFACC_READ)
      do ib=1,13
      start(1)=0
       start(2) = 0
       edges(1) = nc
       edges(2) = nr
       stride(1) = 1
       stride(2) = 1
       sds_index = ib-1
       sds_id    = sfselect(sd_id, sds_index)
C Case 10 meters resolution
       if ((ib.eq.2).or.(ib.eq.3).or.(ib.eq.4).or.(ib.eq.8)) then
       edges(1) = nci
       edges(2) = nri
              write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,ttoaband)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       do i=fc,lc
       do j=fr,lr
       k=i-fc+1
       l=j-fr+1
       toaband(ib,k,l)=ttoaband(i,j)
       enddo
       enddo
       endif
C Case 20 meters resolution
       if ((ib.eq.5).or.(ib.eq.6).or.(ib.eq.7).or.(ib.eq.9).or.(ib.eq.12).or.(ib.eq.13)) then
       edges(1) = nc20i
       edges(2) = nr20i
              write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,tempband1)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       do i=fc,lc
       do j=fr,lr
       k=i-fc+1
       l=j-fr+1
       i20=((i-1)/2)+1
       j20=((j-1)/2)+1
       if (i20.le.0) i20=1
       if (j20.le.0) j20=1
       if (i20.gt.nc20i) i20=nc20i
       if (j20.gt.nr20i) j20=nr20i
       toaband(ib,k,l)=tempband1(i20,j20)
       enddo
       enddo
       endif
C Case 60 meters resolution
       if ((ib.eq.1).or.(ib.eq.10).or.(ib.eq.11)) then
       edges(1) = nc60i
       edges(2) = nr60i
              write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,tempband2)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       do i=fc,lc
       do j=fr,lr
       k=i-fc+1
       l=j-fr+1
       i60=((i-1)/6)+1
       j60=((j-1)/6)+1
       if (i60.le.0) i60=1
       if (j60.le.0) j60=1
       if (i60.gt.nc60i) i60=nc60i
       if (j60.gt.nr60i) j60=nr60i
       toaband(ib,k,l)=tempband2(i60,j60)
       enddo
       enddo
       endif
       
       
       enddo
c close HDF file
       status = sfend(sd_id)
       write(6,*) "status sfend ",status         


c         stop
         
c using scene center to compute atmospheric parameter
         row=fr+nr/2
         col=fc+nc/2
         call utmtodeg(utmzone,row,col,x0,y0,gsize,lat,lon)
c         tlat(i,j)=lat
c         tlon(i,j)=lon
         ycmg=(89.975-lat)/0.05+1.
         xcmg=(179.975+lon)/0.05+1
         icmg=int(ycmg+0.5)
         jcmg=int(xcmg+0.5)
         write(6,*) " icmg jcmg ", icmg,jcmg
         if (wv(jcmg,icmg).ne.0) then
         uwv=wv(jcmg,icmg)/200.
         else
         uwv=0.5
         endif
         
         if (oz(jcmg,icmg).ne.0) then
         uoz=oz(jcmg,icmg)/400.
         else
         uoz=0.3
         endif
         
         if (dem(jcmg,icmg).ne.-9999) then
         pres=1013.*exp(-dem(jcmg,icmg)/8500.)
         else
         pres=1013.
         endif
         raot550nm=0.05
                

        eps=-1.
        xfi=acos(cos((xfs-xfv)*fac))/fac
        do ib=1,13
         iband=ib
         if (iband.ne.10) then
          call atmcorlamb2(xts,xtv,xfi,raot550nm,iband,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       rotoa,roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval,eps,iverbose)
         else
          tgo=1.
          roatm=0.
          ttatmg=1.
          satm=0.
         endif
         if (icor.eq.1) then
          write(6,*) "correction",ib,tgo,roatm,ttatmg,satm,raot550nm,uoz,uwv
         endif
!computing parameter for atmospheric correction        
c call the atmospheric correction
c         if (ib.eq.9) then
c         write(6,*), "ib eq 9", toaband(ib,1,1)
c         stop
c         endif    
         do i=1,nc
         do j=1,nr
          if (toaband(ib,i,j).gt.0.) then
C           roslamb=toaband(ib,i,j)/10000./(1.*dsol)/xmus
           roslamb=toaband(ib,i,j)/10000.
           tband(ib,i,j)=int2(roslamb*10000.)
           if (icor.eq.1) then
            roslamb=roslamb/tgo
            roslamb=roslamb-roatm
            roslamb=roslamb/ttatmg
            roslamb=roslamb/(1.+satm*roslamb)
           endif
           sband(ib,i,j)=int2(roslamb*10000.)
c           sband(ib,i,j)=int((toaband(ib,i,j)*xcal(ib)/100.)*10000.)
c          sband(ib,i,j)=int(toaband(ib,i,j)*10000.)
          else
           sband(ib,i,j)=-100
           tband(ib,i,j)=-100
          endif
         enddo
         enddo
        enddo

c compute water vapor
!note the band 6 and 7 in the tiff file corresponds to 5 and 6 CAVIS
C skip that for now 
!compute the aerosol 
c inverting aerosol
c skip that too
c        goto 997               
        
         do i=1,nc
         if (modulo(i,100).eq.0) write(6,*) "Processing collumn ",i
         do j=1,nr
c         if (i.eq.2376) then
c         write(6,*) "processing collumn ",j
c         endif
         if ((tband(4,i,j).ne.-100).and.(tband(2,i,j).ne.-100).and.(tband(1,i,j).ne.-100)) then
         troatm(1)=0.
         troatm(2)=0.
         troatm(4)=0.
         troatm(13)=0.
         do k=int((i-1)/6)*6+1,int((i-1)/6)*6+6
         do l=int((j-1)/6)*6+1,int((j-1)/6)*6+6
         troatm(1)=troatm(1)+tband(1,k,l)/10000.
         troatm(2)=troatm(2)+tband(2,k,l)/10000.
         troatm(4)=troatm(4)+tband(4,k,l)/10000.
         troatm(13)=troatm(13)+tband(13,k,l)/10000.
         enddo
         enddo
         troatm(1)=troatm(1)/36.
         troatm(2)=troatm(2)/36.
         troatm(4)=troatm(4)/36.
         troatm(13)=troatm(13)/36.
         iband1=4
         iband3=1
         row=fr+(j-1) 
         col=fc+(i-1)
         call utmtodeg(utmzone,row,col,x0,y0,gsize,lat,lon)
         tlat(i,j)=lat
         tlon(i,j)=lon
         ycmg=(89.975-lat)/0.05+1.
         xcmg=(179.975+lon)/0.05+1
         txcmg(i,j)=xcmg
         tycmg(i,j)=ycmg
         icmg=int(ycmg+0.5)
         jcmg=int(xcmg+0.5)
c         twv(i,j)=wv(jcmg,icmg)
         icmg=int(ycmg)
         jcmg=int(xcmg)
         u=(ycmg-icmg)
         v=(xcmg-jcmg)
         twvi(i,j)=wv(jcmg,icmg)*(1.-u)*(1.-v)+wv(jcmg+1,icmg)*(1.-u)*v
     s              +wv(jcmg,icmg+1)*(1.-v)*u+wv(jcmg+1,icmg+1)*u*v
         twvi(i,j)=twvi(i,j)/100.
         
         
         
         if (oz(jcmg,icmg).lt.0) then
         uoz11=256+oz(jcmg,icmg)
         else
         uoz11=oz(jcmg,icmg)
         endif
         
         if (oz(jcmg+1,icmg).lt.0) then
         uoz12=256+oz(jcmg+1,icmg)
         else
         uoz12=oz(jcmg+1,icmg)
         endif
         
         if (oz(jcmg,icmg+1).lt.0) then
         uoz21=256+oz(jcmg,icmg+1)
         else
         uoz21=oz(jcmg,icmg+1)
         endif
         
         if (oz(jcmg+1,icmg+1).lt.0) then
         uoz22=256+oz(jcmg+1,icmg+1)
         else
         uoz22=oz(jcmg+1,icmg+1)
         endif
         
         if (uoz11.eq.0) uoz11=120
         if (uoz12.eq.0) uoz12=120
         if (uoz21.eq.0) uoz21=120
         if (uoz22.eq.0) uoz22=120
         
          
         
         
         tozi(i,j)=uoz11*(1.-u)*(1.-v)+uoz12*(1.-u)*v
     s              +uoz21*(1.-v)*u+uoz22*u*v
     
         tozi(i,j)=tozi(i,j)/400.
         
         if (dem(jcmg,icmg).ne.-9999) then
         pres11=1013.*exp(-dem(jcmg,icmg)/8500.)
         else
         pres11=1013.
         cloud(i,j)=-128
         endif
         if (dem(jcmg+1,icmg).ne.-9999) then
         pres12=1013.*exp(-dem(jcmg+1,icmg)/8500.)
         else
         pres12=1013.
         endif
         if (dem(jcmg,icmg+1).ne.-9999) then
         pres21=1013.*exp(-dem(jcmg,icmg+1)/8500.)
         else
         pres21=1013.
         endif
         if (dem(jcmg+1,icmg+1).ne.-9999) then
         pres22=1013.*exp(-dem(jcmg+1,icmg+1)/8500.)
         else
         pres22=1013.
         endif
         tp(i,j)=pres11*(1.-u)*(1.-v)+pres12*(1.-u)*v
     s              +pres21*(1.-v)*u+pres22*u*v

c inverting aerosol               
        do ib=1,13
        erelc(ib)=-1.
        enddo

        rb1=ratiob1(jcmg,icmg)/1000.
        rb2=ratiob2(jcmg,icmg)/1000.
        if ((rb2.gt.1.0).or.(rb1.gt.1.0).or.(rb2.lt.0.1).or.(rb1.lt.0.1)) then
        slpratiob1(jcmg,icmg)=0
        slpratiob2(jcmg,icmg)=0
        slpratiob7(jcmg,icmg)=0
        intratiob1(jcmg,icmg)=550
        intratiob2(jcmg,icmg)=600
        intratiob7(jcmg,icmg)=2000
        else
        if (sndwi(jcmg,icmg).lt.200) then
        slpratiob1(jcmg,icmg)=0
        slpratiob2(jcmg,icmg)=0
        slpratiob7(jcmg,icmg)=0
        intratiob1(jcmg,icmg)=ratiob1(jcmg,icmg)
        intratiob2(jcmg,icmg)=ratiob2(jcmg,icmg)
        intratiob7(jcmg,icmg)=ratiob7(jcmg,icmg)
        endif
        endif

        

        rb1=ratiob1(jcmg+1,icmg)/1000.
        rb2=ratiob2(jcmg+1,icmg)/1000.
        if ((rb2.gt.1.0).or.(rb1.gt.1.0).or.(rb2.lt.0.1).or.(rb1.lt.0.1)) then
        slpratiob1(jcmg+1,icmg)=0
        slpratiob2(jcmg+1,icmg)=0
        slpratiob7(jcmg+1,icmg)=0
        intratiob1(jcmg+1,icmg)=550
        intratiob2(jcmg+1,icmg)=600
        intratiob7(jcmg+1,icmg)=2000
        else
        if (sndwi(jcmg+1,icmg).lt.200) then
        slpratiob1(jcmg+1,icmg)=0
        slpratiob2(jcmg+1,icmg)=0
        slpratiob7(jcmg+1,icmg)=0
        intratiob1(jcmg+1,icmg)=ratiob1(jcmg+1,icmg)
        intratiob2(jcmg+1,icmg)=ratiob2(jcmg+1,icmg)
        intratiob7(jcmg+1,icmg)=ratiob7(jcmg+1,icmg)
        endif
        endif

        

        rb1=ratiob1(jcmg,icmg+1)/1000.
        rb2=ratiob2(jcmg,icmg+1)/1000.
        if ((rb2.gt.1.0).or.(rb1.gt.1.0).or.(rb2.lt.0.1).or.(rb1.lt.0.1)) then
        slpratiob1(jcmg,icmg+1)=0
        slpratiob2(jcmg,icmg+1)=0
        slpratiob7(jcmg,icmg+1)=0
        intratiob1(jcmg,icmg+1)=550
        intratiob2(jcmg,icmg+1)=600
        intratiob7(jcmg,icmg+1)=2000
        else
        if (sndwi(jcmg,icmg+1).lt.200) then
        slpratiob1(jcmg,icmg+1)=0
        slpratiob2(jcmg,icmg+1)=0
        slpratiob7(jcmg,icmg+1)=0
        intratiob1(jcmg,icmg+1)=ratiob1(jcmg,icmg+1)
        intratiob2(jcmg,icmg+1)=ratiob2(jcmg,icmg+1)
        intratiob7(jcmg,icmg+1)=ratiob7(jcmg,icmg+1)
        endif
        endif



        rb1=ratiob1(jcmg+1,icmg+1)/1000.
        rb2=ratiob2(jcmg+1,icmg+1)/1000.

        if ((rb2.gt.1.0).or.(rb1.gt.1.0).or.(rb2.lt.0.1).or.(rb1.lt.0.1)) then
        slpratiob1(jcmg+1,icmg+1)=0
        slpratiob2(jcmg+1,icmg+1)=0
        slpratiob7(jcmg+1,icmg+1)=0
        intratiob1(jcmg+1,icmg+1)=550
        intratiob2(jcmg+1,icmg+1)=600
        intratiob7(jcmg+1,icmg+1)=2000
        else
        if (sndwi(jcmg+1,icmg+1).lt.200) then
        slpratiob1(jcmg+1,icmg+1)=0
        slpratiob2(jcmg+1,icmg+1)=0
        slpratiob7(jcmg+1,icmg+1)=0
        intratiob1(jcmg+1,icmg+1)=ratiob1(jcmg+1,icmg+1)
        intratiob2(jcmg+1,icmg+1)=ratiob2(jcmg+1,icmg+1)
        intratiob7(jcmg+1,icmg+1)=ratiob7(jcmg+1,icmg+1)
        endif
        endif



         slpr11=slpratiob2(jcmg,icmg)/1000.
         intr11=intratiob2(jcmg,icmg)/1000.
         slpr12=slpratiob2(jcmg+1,icmg)/1000.
         intr12=intratiob2(jcmg+1,icmg)/1000.
         slpr21=slpratiob2(jcmg,icmg+1)/1000.
         intr21=intratiob2(jcmg,icmg+1)/1000.
         slpr22=slpratiob2(jcmg+1,icmg+1)/1000.
         intr22=intratiob2(jcmg+1,icmg+1)/1000.
         slprb2(i,j)=slpr11*(1.-u)*(1.-v)+slpr12*(1.-u)*v
     s              +slpr21*(1.-v)*u+slpr22*u*v
         intrb2(i,j)=intr11*(1.-u)*(1.-v)+intr12*(1.-u)*v
     s              +intr21*(1.-v)*u+intr22*u*v

         slpr11=slpratiob1(jcmg,icmg)/1000.
         intr11=intratiob1(jcmg,icmg)/1000.
         slpr12=slpratiob1(jcmg+1,icmg)/1000.
         intr12=intratiob1(jcmg+1,icmg)/1000.
         slpr21=slpratiob1(jcmg,icmg+1)/1000.
         intr21=intratiob1(jcmg,icmg+1)/1000.
         slpr22=slpratiob1(jcmg+1,icmg+1)/1000.
         intr22=intratiob1(jcmg+1,icmg+1)/1000.
         slprb1(i,j)=slpr11*(1.-u)*(1.-v)+slpr12*(1.-u)*v
     s              +slpr21*(1.-v)*u+slpr22*u*v
         intrb1(i,j)=intr11*(1.-u)*(1.-v)+intr12*(1.-u)*v
     s              +intr21*(1.-v)*u+intr22*u*v

         slpr11=slpratiob7(jcmg,icmg)/1000.
         intr11=intratiob7(jcmg,icmg)/1000.
         slpr12=slpratiob7(jcmg+1,icmg)/1000.
         intr12=intratiob7(jcmg+1,icmg)/1000.
         slpr21=slpratiob7(jcmg,icmg+1)/1000.
         intr21=intratiob7(jcmg,icmg+1)/1000.
         slpr22=slpratiob7(jcmg+1,icmg+1)/1000.
         intr22=intratiob7(jcmg+1,icmg+1)/1000.
         slprb7(i,j)=slpr11*(1.-u)*(1.-v)+slpr12*(1.-u)*v
     s              +slpr21*(1.-v)*u+slpr22*u*v
         intrb7(i,j)=intr11*(1.-u)*(1.-v)+intr12*(1.-u)*v
     s              +intr21*(1.-v)*u+intr22*u*v


        xndwi=(float(sband(9,i,j))-float(sband(13,i,j))/2.)
     s       /(float(sband(9,i,j))+float(sband(13,i,j))/2.)
        th1=((andwi(jcmg,icmg)+2.*sndwi(jcmg,icmg))/1000.)
        th2=((andwi(jcmg,icmg)-2.*sndwi(jcmg,icmg))/1000.)
        if (xndwi.gt.th1) xndwi=th1
        if (xndwi.lt.th2) xndwi=th2
        erelc(1)=(xndwi*slprb1(i,j)+intrb1(i,j))
        erelc(2)=(xndwi*slprb2(i,j)+intrb2(i,j)) 
        erelc(13)=(xndwi*slprb7(i,j)+intrb7(i,j))
         erelc(4)=1.
        tratiob1(i,j)=int2(erelc(1)*1000.)
        tratiob2(i,j)=int2(erelc(2)*1000.)
        tratiob7(i,j)=int2(erelc(13)*1000.)

       if ((modulo(i,6).ne.1).or.(modulo(j,6).ne.1)) goto 1998
       iband1=4
       iband3=1
       pres=tp(i,j)
       uoz=tozi(i,j)
       uwv=twvi(i,j)
c       if ((j.eq.319).and.(i.eq.499)) then
c       iverbose=1
c       else
c       iverbose=0
c       endif
       eps=1
       ifast=0
       iaots=1
       

       call subaeroretv3(iband1,iband3,xts,xtv,xfi,pres,uoz,uwv,erelc,troatm,
     c       tpres,aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,iverbose,raot,residual,nit,next,ierr,
     s       ifast,tbtgo,tbroatm,tbttatmg,tbsatm,iaots,eps)
     
              if (iverbose.eq.1) then
         write(6,*) "aero retrieved i,j,raot,eps,residual" ,i,j,raot,eps,residual
         endif
         
         eps1=eps
         residual1=residual
         sraot1=raot

       eps=1.75
       call subaeroretv3(iband1,iband3,xts,xtv,xfi,pres,uoz,uwv,erelc,troatm,
     c       tpres,aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,iverbose,raot,residual,nit,next,ierr,
     s       ifast,tbtgo,tbroatm,tbttatmg,tbsatm,iaots,eps)
              if (iverbose.eq.1) then
         write(6,*) "aero retrieved i,j,raot,eps,residual" ,i,j,raot,eps,residual
         endif
         eps2=eps
         residual2=residual
         sraot2=raot
         
       eps=2.5
       call subaeroretv3(iband1,iband3,xts,xtv,xfi,pres,uoz,uwv,erelc,troatm,
     c       tpres,aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,iverbose,raot,residual,nit,next,ierr,
     s       ifast,tbtgo,tbroatm,tbttatmg,tbsatm,iaots,eps)
              if (iverbose.eq.1) then
         write(6,*) "aero retrieved i,j,raot,eps,residual" ,i,j,raot,eps,residual
         endif
         eps3=eps
         residual3=residual
         sraot3=raot

c find eps that makes the residual near zero         
        xa=(eps1*eps1)-(eps3*eps3)
        xd=int((eps2*eps2)-(eps3*eps3))
        xb=(eps1-eps3)
        xe=(eps2-eps3)
        xc=residual1-residual3
        xf=residual2-residual3
        coefa=(xc*xe-xb*xf)/(xa*xe-xb*xd)
        coefb=(xa*xf-xc*xd)/(xa*xe-xb*xd)
        coefc=(residual1-coefa*eps1*eps1-coefb*eps1)
c        cdet=(coefb*coefb-4*coefa*coefc)
c      local extremum
       epsmin=-coefb/(2.*coefa)
       resepsmin=xa*epsmin*epsmin+xb*epsmin+xc
       if (iverbose.eq.1) then
       write(6,*) "coef eps",coefa,coefb,coefc,epsmin,resepsmin
       endif
       if ((epsmin.lt.1.0).or.(epsmin.gt.2.5)) then
       if (residual1.lt.residual3) then 
       epsmin=eps1
       else
       epsmin=eps3
       endif
       else
       if ((resepsmin.gt.residual1).or.(resepsmin.gt.residual3)) then
         if (residual1.lt.residual3) then 
         epsmin=eps1
         else
         epsmin=eps3
         endif
       endif
       endif
      
       
       
       
 
       
       if (iverbose.eq.1) then
       write(6,*) "eps retrieved",epsmin
       endif
       eps=epsmin
c       eps=-1.
       call subaeroretv3(iband1,iband3,xts,xtv,xfi,pres,uoz,uwv,erelc,troatm,
     c       tpres,aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,iverbose,raot,residual,nit,next,ierr,
     s       ifast,tbtgo,tbroatm,tbttatmg,tbsatm,iaots,eps)
         
              if (iverbose.eq.1) then
         write(6,*) "aero retrieved i,j,raot,eps,residual" ,i,j,raot,eps,residual
         endif





         tnit(i,j)=int2(nit)
         taero3(i,j)=next*raot
         corf=raot/xmus


c test if band5 makes sense
        corfb=1.
c ESPA the troatm(7) is wrong.  We've changed it.
c        if (residual.lt.((0.015+0.005*corf+0.10*troatm(7))*corfb)) then
        if (residual.lt.((0.015+0.005*corf+0.10*troatm(13))*corfb)) then
        rotoa=0.
        iband=9
         do k=int((i-1)/6)*6+1,int((i-1)/6)*6+6
         do l=int((j-1)/6)*6+1,int((j-1)/6)*6+6
         rotoa=rotoa+tband(9,k,l)/10000.
         enddo
         enddo
         rotoa=rotoa/36.
        raot550nm=raot
        call atmcorlamb2(xts,xtv,xfi,raot550nm,iband,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       rotoa,roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval,eps,iverbose)
        ros5=roslamb
        rotoa=0.
        iband=4
         do k=int((i-1)/6)*6+1,int((i-1)/6)*6+6
         do l=int((j-1)/6)*6+1,int((j-1)/6)*6+6
         rotoa=rotoa+tband(4,k,l)/10000.
         enddo
         enddo
         rotoa=rotoa/36.
        raot550nm=raot
        call atmcorlamb2(xts,xtv,xfi,raot550nm,iband,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       rotoa,roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval,eps,iverbose)
         ros4=roslamb
            if ((ros5.gt.0.1).and.((ros5-ros4)/(ros4+ros5).gt.0)) then
         taero(i,j)=raot
         teps(i,j)=eps
         tresi(i,j)=residual
         ipflag(i,j)=0
         else
c this is water the retrieval needs to be redone
         taero(i,j)=raot
         teps(i,j)=eps
         tresi(i,j)=residual
         ipflag(i,j)=-1
         endif
          else
         taero(i,j)=raot
         teps(i,j)=eps
         tresi(i,j)=residual
         ipflag(i,j)=-1
         endif
c redo the retrieval if water          
        if (ipflag(i,j).eq.-1) then
        do ib=1,16
        erelc(ib)=-1.
        enddo
c        write(6,*) "aeroband1 "
         troatm(1)=0.
         troatm(4)=0.
         troatm(9)=0.
         troatm(13)=0.
         do k=int((i-1)/6)*6+1,int((i-1)/6)*6+6
         do l=int((j-1)/6)*6+1,int((j-1)/6)*6+6
      troatm(1)=troatm(1)+tband(1,k,l)/10000.
c ESPA the troatm(2) is wrong.  We've changed it.
cORIG troatm(4)=troatm(2)+tband(4,k,l)/10000.
      troatm(4)=troatm(4)+tband(4,k,l)/10000.
c ESPA the troatm(4) is wrong.  We've changed it.
cORIG troatm(9)=troatm(4)+tband(9,k,l)/10000.
      troatm(9)=troatm(9)+tband(9,k,l)/10000.
c ESPA the troatm(7) is wrong.  We've changed it.
cORIG troatm(13)=troatm(7)+tband(13,k,l)/10000.
      troatm(13)=troatm(13)+tband(13,k,l)/10000.
         enddo
         enddo
c        write(6,*) "aeroband1 "
         troatm(1)=troatm(1)/36.
         troatm(4)=troatm(4)/36.
         troatm(9)=troatm(9)/36.
         troatm(13)=troatm(13)/36.
        erelc(13)=1.0
        erelc(9)=1.0
        erelc(4)=1.0
        erelc(1)=1.0
        pres=tp(i,j)
        uoz=tozi(i,j)
        uwv=twvi(i,j)

       ifast=0
       eps=1.5
       iaots=1
       call subaeroretwat(iband1,iband3,xts,xtv,xfi,pres,uoz,uwv,erelc,troatm,
     c       tpres,aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,iverbose,raot,residual,nit,next,ierr,
     s       ifast,tbtgo,tbroatm,tbttatmg,tbsatm,iaots,eps)
c         write(6,*) "aeroband1 last "
    
c        if (j.eq.7711 .and. i.eq.3619) then
c          iverbose = 0
c        endif
c              if (iverbose.eq.1) then
c         write(6,*) "aero retrieved i,j,raot,eps,residual" ,i,j,raot,eps,residual
c         endif
         
c              if (iverbose.eq.1) then
c         write(6,*) "aero retrieved i,j,raot,eps,residual" ,i,j,raot,eps,residual
c         endif
c         endif
         tnit(i,j)=int2(nit)
         taero3(i,j)=next*raot
         teps(i,j)=eps
         corf=raot/xmus
         taero(i,j)=raot
         tresi(i,j)=residual

c test band 1 reflectance eliminate negative         
        iband=1
        rotoa=troatm(1)
        raot550nm=raot
        call atmcorlamb2(xts,xtv,xfi,raot550nm,iband,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       rotoa,roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval,eps,iverbose)
        rosband1=roslamb
         
         if ((residual.gt.((0.010+0.005*corf))).or.(rosband1.lt.0)) then
c          write(6,*) "problem in subaeroretwat at i,j ",i,j
         ipflag(i,j)=2
          else         
         tnit(i,j)=int2(nit)
         taero3(i,j)=next*raot
         teps(i,j)=eps
         corf=raot/xmus
         taero(i,j)=raot
         tresi(i,j)=residual
         endif
           endif 
        
        

         endif 
c endif of filled value         
 1998    continue

        if ((modulo(i,6).eq.1).and.(modulo(j,6).eq.1)) then
        do k=i,i+5
        do l=j,j+5
         tnit(k,l)=tnit(i,j)
         taero3(k,l)=taero3(i,j)
         teps(k,l)=teps(i,j)
          taero(k,l)=taero(i,j)
         tresi(k,l)=tresi(i,j)
        enddo
        enddo
        endif
         
          enddo
         enddo

#define WRITE_TAERO
#ifdef WRITE_TAERO
       OPEN(1, FILE="aerosols.img", ACCESS="stream",
     +      FORM="unformatted")
       WRITE(1) taero
       CLOSE(1)
       OPEN(1, FILE="ipflag.img", ACCESS="stream", FORM="unformatted")
       WRITE(1) ipflag
       CLOSE(1)
#endif

c interpolate taero
       do j=1,nr
       do i=1,nc
        do k=i,i+5
        do l=j,j+5
         taero(k,l) = (i+6-k)*(j+6-l)*taero(i,j)
         if(((i+6).le.nc).and.((j+6).le.nr)) then
            taero(k,l)=taero(k,l)+(k-i)*(j+6-l)*taero(i+6,j)+(i+6-k)*(l-j)*taero(i,j+6)+(k-i)*(l-j)*taero(i+6,j+6)
            endif
         if(((i+6).le.nc).and.((j+6).gt.nr)) then
            taero(k,l) = taero(k,l)+(k-i)*(j+6-l)*taero(i+6,j) + (i+6-k)*(l-j)*taero(i,j) + (k-i)*(l-j)*taero(i+6,j)
            endif
         if(((i+6).gt.nc).and.((j+6).le.nr)) then
            taero(k,l)=taero(k,l)+(k-i)*(j+6-l)*taero(i,j)+(i+6-k)*(l-j)*taero(i,j+6)+(k-i)*(l-j)*taero(i,j+6)
            endif
         if(((i+6).gt.nc).and.((j+6).gt.nr)) then
            taero(k,l)=taero(k,l)+(k-i)*(j+6-l)*taero(i,j)+(i+6-k)*(l-j)*taero(i,j)+(k-i)*(l-j)*taero(i,j)
         endif
         taero(k,l) = taero(k,l)/36.
       enddo
       enddo
       enddo
       enddo
         
#ifdef WRITE_TAERO
       OPEN(1, FILE="aerosols2.img", ACCESS="stream",
     +      FORM="unformatted")
       WRITE(1) taero
       CLOSE(1)
       OPEN(1, FILE="ipflag2.img", ACCESS="stream",
     +      FORM="unformatted")
       WRITE(1) ipflag
       CLOSE(1)
#endif

       do j=1,nr
       do i=1,nc
       if ((ipflag(i,j).gt.1).and.(ipflag(i,j).lt.5)) then
       do k=i-12,i+12
       do l=j-12,j+12
       if ((k.ge.1).and.(l.ge.1).and.(k.le.nc).and.(l.le.nr)) then
       if (ipflag(k,l).eq.0) ipflag(k,l)=9
       endif
       enddo
       enddo
       endif
       enddo
       enddo
       
       do j=1,nr
       do i=1,nc
       if (ipflag(i,j).eq.9) ipflag(i,j)=2
       enddo
       enddo

#ifdef WRITE_TAERO
       OPEN(1, FILE="aerosols3.img", ACCESS="stream",
     +      FORM="unformatted")
       WRITE(1) taero
       CLOSE(1)
       OPEN(1, FILE="ipflag3.img", ACCESS="stream",
     +      FORM="unformatted")
       WRITE(1) ipflag
       CLOSE(1)
#endif
       nbpixnf=0
       nbpixtot=0
       do j=1,nr
       do i=1,nc
       smflag(i,j)=0
       if ((tband(4,i,j).ne.-100).and.(tband(2,i,j).ne.-100).and.(tband(1,i,j).ne.-100)) then
       nbpixtot=nbpixtot+1
       taeroavg=0.
       tepsavg=0.
       nbaeroavg=0
       do k=max(j-30,1),min(j+30,nr)
       do l=max(i-30,1),min(i+30,nc)
c       if (j.le.5) write(6,*) "i,j,k ",i,j,k
       if ((ipflag(l,k).eq.0).or.(ipflag(l,k).eq.-1)) then
       nbaeroavg=nbaeroavg+1
       taeroavg=taero(l,k)+taeroavg
       tepsavg=teps(l,k)+tepsavg
       endif
       enddo
       enddo
c       if (j.le.5) write(6,*) "i,j,nbaeroavg",i,j,nbaeroavg
       if (nbaeroavg.gt.20) then
       taeroavg=taeroavg/nbaeroavg
       tepsavg=tepsavg/nbaeroavg
       taeros(i,j)=taeroavg
       tepss(i,j)=tepsavg
       smflag(i,j)=1
       else
c       write(6,*) "Pixels at i,j not filled",i,j
       nbpixnf=nbpixnf+1
       endif
       endif
       enddo
       enddo       
       
       if (nbpixnf.eq.nbpixtot) then
        do j=1,nr
       do i=1,nc
        if ((tband(4,i,j).ne.-100).and.(tband(2,i,j).ne.-100).and.(tband(1,i,j).ne.-100)) then
        taeros(i,j)=0.05
       tepss(i,j)=1.5
       smflag(i,j)=1
       endif
       enddo
       enddo
       goto 121
       endif


       write(6,*) "Second Pass"
c while nbpixnf.ne.0  
       ipass=2
       wind=30
       do while (nbpixnf.ne.0) 
       write(6,*) " Pass number ",ipass, "nbpixnf ",nbpixnf     
       nbpixnf=0
       do j=1,nr
       do i=1,nc
       if((tband(4,i,j).ne.-100).and.(tband(2,i,j).ne.-100).and.(tband(1,i,j).ne.-100).and.(smflag(i,j).eq.0)) then
       tepsavg=0.
       taeroavg=0.
       nbaeroavg=0
       do k=max(j-wind,1),min(j+wind,nr)
       do l=max(i-wind,1),min(i+wind,nc)
       if (smflag(l,k).eq.1) then
       nbaeroavg=nbaeroavg+1
       taeroavg=taeros(l,k)+taeroavg
       tepsavg=tepss(l,k)+tepsavg
      endif
       enddo
       enddo
       if (nbaeroavg.gt.0) then
       taeroavg=taeroavg/nbaeroavg
       tepsavg=tepsavg/nbaeroavg
       taeros(i,j)=taeroavg
       tepss(i,j)=tepsavg
       smflag(i,j)=1
       else
c       write(6,*) "Second pass Pixels at i,j not filled",i,j
       nbpixnf=nbpixnf+1
       endif
       endif
       enddo
       enddo   
       ipass=ipass+1
c       wind=wind+10
       enddo    
       
c end of averaging 

c
 121   continue
       do j=1,nr
       do i=1,nc
       if ((i.eq.408).and.(j.eq.254)) then
       write(6,*) "verbose on in interpolation ",i,j,ipflag(i,j)
       iverbose=1
       else
       iverbose=0
       endif
       
        if ((ipflag(i,j).gt.1).and.(ipflag(i,j).lt.5)) then
       taero(i,j)=taeros(i,j)
       teps(i,j)=tepss(i,j)
       ipflag(i,j)=1
       endif
       enddo
       enddo

#ifdef WRITE_TAERO
       OPEN(1, FILE="aerosols4.img", ACCESS="stream",
     +      FORM="unformatted")
       WRITE(1) taero
       CLOSE(1)
       OPEN(1, FILE="ipflag4.img", ACCESS="stream",
     +      FORM="unformatted")
       WRITE(1) ipflag
       CLOSE(1)
#endif
            
                      
c perform atmospheric correction
         do i=1,nc
         if (modulo(i,100).eq.0) write(6,*) "correcting collumn ",i
         do j=1,nr
c         if (i.eq.2376) then
c         write(6,*) "processing collumn ",j
c         endif
         if ((tband(4,i,j).ne.-100).and.(tband(2,i,j).ne.-100).and.(tband(1,i,j).ne.-100)) then
          raot550nm=taero(i,j)
          eps=teps(i,j)
c          raot550nm=0.9
c          eps=-1.
          pres=tp(i,j)
          uwv=twvi(i,j)
          uoz=tozi(i,j)
        do ib=1,13
         if (ib.ne.11) then
         iband=ib
         rotoa=tband(ib,i,j)/10000.
        call atmcorlamb2(xts,xtv,xfi,raot550nm,iband,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       rotoa,roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval,eps,iverbose)
           if (ib.eq.1) then
           rsurf=sband(ib,i,j)/10000.
           aerimp(i,j)=int(min(int(abs(rsurf-roslamb)*1000),255), 1)
           endif
           sband(ib,i,j)=int2(roslamb*10000.)
           else
           sband(ib,i,j)=tband(ib,i,j)
           endif
           if ((i.eq.1721).and.(j.eq.1905)) then
            write(6,*) "ib ",ib,raot550nm,rotoa,roslamb
           endif
           enddo
           endif
           enddo
           enddo
           
c update cloud mask QA         
         do i=1,nc
         do j=1,nr
c cloud
         if ((sband(11,i,j).gt.30).or.(tresi(i,j).gt.0.05))
     s       cloud(i,j)=int(cloud(i,j)+2, 1)
c shadow
         if ((sband(1,i,j).lt.100).or.(sband(2,i,j).lt.0))
     s       cloud(i,j)=int(cloud(i,j)+8, 1)
c aerosol qa
         if (aerimp(i,j).lt.15) then
         cloud(i,j)=int(cloud(i,j)+16, 1)
         else
         if (aerimp(i,j).lt.30) then
         cloud(i,j)=int(cloud(i,j)+32, 1)
         else
         cloud(i,j)=int(cloud(i,j)+48, 1)
         endif
         endif
         enddo
         enddo
         
                     
         
         
!Saving in an HDF file
 
C 997     write(6,*) "fileout1 ",fileout1
         write(6,*) "fileout1 ",fileout1
         sd_id= sfstart(adjustl(fileout1),DFACC_CREATE)
!writing surface bands
         do ib=1,8
         do i=1,nc
         do j=1,nr
         oband(i,j)=sband(ib,i,j)
           if ((i.eq.171).and.(j.eq.6905)) then
            write(6,*) "ib ",ib,raot550nm,rotoa,roslamb
           endif
           if ((i.eq.1721).and.(j.eq.6905)) then
            write(6,*) "ib ",ib,raot550nm,rotoa,roslamb
           endif
           if ((i.eq.3721).and.(j.eq.6905)) then
            write(6,*) "ib ",ib,raot550nm,rotoa,roslamb
           endif
         enddo
         enddo
         if ((ib.eq.2).or.(ib.eq.3).or.(ib.eq.4)) then
         if (ib.eq.2) write(sds_name,'(A4)') "blue"
         if (ib.eq.3) write(sds_name,'(A5)') "green"
         if (ib.eq.4) write(sds_name,'(A3)') "red"
         else
         write(sds_name,'(A4,A2)') "band",suffix(ib)
         endif
         write(6,*) "writing band ",sds_name
         dim_length(1)=nc
         dim_length(2)=nr
         comp_type=4
         comp_prm(1)=8
         rank=2
         dim_sizes(1)=nc
         dim_sizes(2)=nr
         start(1) = 0
         start(2) = 0
          edges(1) = nc
         edges(2) = nr
         stride(1) = 1
         stride(2) = 1
         sds_id=sfcreate(sd_id,sds_name,DFNT_INT16,rank,dim_sizes)
         status=sfwdata(sds_id,start,stride,edges,oband)
         write(6,*) "status sfwdata ",status
         status = sfendacc(sds_id)
         write(6,*) "status sfendacc ",status
         enddo
         status = sfend(sd_id)         
         
         write(6,*) "fileout2 ",fileout2
         sd_id= sfstart(adjustl(fileout2),DFACC_CREATE)
         do ib=9,13
         do i=1,nc
         do j=1,nr
         oband(i,j)=sband(ib,i,j)
         enddo
         enddo
         if ((ib.eq.2).or.(ib.eq.3).or.(ib.eq.4)) then
         if (ib.eq.2) write(sds_name,'(A4)') "blue"
         if (ib.eq.3) write(sds_name,'(A5)') "green"
         if (ib.eq.4) write(sds_name,'(A3)') "red"
         else
         write(sds_name,'(A4,A2)') "band",suffix(ib)
         endif
         write(6,*) "writing band ",sds_name
         dim_length(1)=nc
         dim_length(2)=nr
         comp_type=4
         comp_prm(1)=8
         rank=2
         dim_sizes(1)=nc
         dim_sizes(2)=nr
         start(1) = 0
         start(2) = 0
          edges(1) = nc
         edges(2) = nr
         stride(1) = 1
         stride(2) = 1
         sds_id=sfcreate(sd_id,sds_name,DFNT_INT16,rank,dim_sizes)
         status=sfwdata(sds_id,start,stride,edges,oband)
         write(6,*) "status sfwdata ",status
         status = sfendacc(sds_id)
         write(6,*) "status sfendacc ",status
         enddo
         write(sds_name,'(A5)') "CLOUD"
         sds_id=sfcreate(sd_id,sds_name,DFNT_UINT8,rank,dim_sizes)
         status=sfwdata(sds_id,start,stride,edges,cloud)
         write(6,*) "status sfwdata ",status
         write(6,*) "writing attribute data for CLOUD band ",sds_id
         status=set_qamap(sds_id)
         write(6,*) "status set qa attribute ",status
         status = sfendacc(sds_id)
         write(6,*) "status sfendacc ",status

         
!writing toa bands
         if (iopt.eq.3) then
c saving toa bands for validation         
          do ib=1,13
          do i=1,nc
          do j=1,nr
          oband(i,j)=tband(ib,i,j)
          enddo
          enddo
         write(sds_name,'(A7,A2)') "toaband",suffix(ib)
          write(6,*) "writing band ",sds_name
          dim_length(1)=nc
          dim_length(2)=nr
          comp_type=4
          comp_prm(1)=8
          rank=2
          dim_sizes(1)=nc
          dim_sizes(2)=nr
          start(1) = 0
          start(2) = 0
           edges(1) = nc
          edges(2) = nr
          stride(1) = 1
          stride(2) = 1
          sds_id=sfcreate(sd_id,sds_name,DFNT_INT16,rank,dim_sizes)
          status=sfwdata(sds_id,start,stride,edges,oband)
          write(6,*) "status sfwdata ",status
          status = sfendacc(sds_id)
          write(6,*) "status sfendacc ",status
          enddo
         endif

! writing water vapor
c         do i=1,nc
c         do j=1,nr
c         oband(i,j)=wvband(i,j)
c         enddo
c         enddo
c        write(sds_name,'(A2)') "wv"
c         write(6,*) "writing band ",sds_name
c         dim_length(1)=nc
c         dim_length(2)=nr
c         comp_type=4
c         comp_prm(1)=8
c         rank=2
c         dim_sizes(1)=nc
c         dim_sizes(2)=nr
c        start(1) = 0
c        start(2) = 0
c         edges(1) = nc
c        edges(2) = nr
c        stride(1) = 1
c        stride(2) = 1
c         sds_id=sfcreate(sd_id,sds_name,DFNT_INT16,rank,dim_sizes)
c        status=sfwdata(sds_id,start,stride,edges,oband)
c         write(6,*) "status sfwdata ",status
c        status = sfendacc(sds_id)
         write(6,*) "status sfendacc ",status    
! writing aot
         if  (iopt.eq.3) then
         write(sds_name,'(A11)') "taero-band1"
        sds_id=sfcreate(sd_id,sds_name,DFNT_FLOAT32,rank,dim_sizes)
         status=sfwdata(sds_id,start,stride,edges,taero3)
        write(6,*) "status sfwdata ",status
         status = sfendacc(sds_id)
         write(sds_name,'(A5)') "taero"
        sds_id=sfcreate(sd_id,sds_name,DFNT_FLOAT32,rank,dim_sizes)
         status=sfwdata(sds_id,start,stride,edges,taero)
        write(6,*) "status sfwdata ",status
         status = sfendacc(sds_id)
         write(sds_name,'(A5)') "tresi"
        sds_id=sfcreate(sd_id,sds_name,DFNT_FLOAT32,rank,dim_sizes)
         status=sfwdata(sds_id,start,stride,edges,tresi)
        write(6,*) "status sfwdata ",status
         status = sfendacc(sds_id)
           write(sds_name,'(A8)') "tratiob1"
           write(6,*) "tratiob1 ",tratiob1(672,709),tratiob1(709,672)
        sds_id=sfcreate(sd_id,sds_name,DFNT_INT16,rank,dim_sizes)
         status=sfwdata(sds_id,start,stride,edges,tratiob1)
         write(6,*) "status sfwdata ",status
c         status = sfendacc(sds_id)
         scalefactor=1000.
         offset=0.0
         status=sfsattr(sds_id, "scale_factor", DFNT_FLOAT32, 1, scalefactor)
         status=sfsattr(sds_id, "offset", DFNT_FLOAT32, 1, offset)
         status = sfendacc(sds_id)
           write(sds_name,'(A8)') "tratiob2"
        sds_id=sfcreate(sd_id,sds_name,DFNT_INT16,rank,dim_sizes)
         status=sfwdata(sds_id,start,stride,edges,tratiob2)
         write(6,*) "status sfwdata ",status
c         status = sfendacc(sds_id)
         scalefactor=1000.
         offset=0.0
         status=sfsattr(sds_id, "scale_factor", DFNT_FLOAT32, 1, scalefactor)
         status=sfsattr(sds_id, "offset", DFNT_FLOAT32, 1, offset)
         status = sfendacc(sds_id)
           write(sds_name,'(A8)') "tratiob7"
        sds_id=sfcreate(sd_id,sds_name,DFNT_INT16,rank,dim_sizes)
         status=sfwdata(sds_id,start,stride,edges,tratiob7)
         write(6,*) "status sfwdata ",status
c         status = sfendacc(sds_id)
         scalefactor=1000.
         offset=0.0
         status=sfsattr(sds_id, "scale_factor", DFNT_FLOAT32, 1, scalefactor)
         status=sfsattr(sds_id, "offset", DFNT_FLOAT32, 1, offset)
         status = sfendacc(sds_id)
         endif
         
         if (iopt.eq.3) then
         write(sds_name,'(A4)') "twvi"
          sds_id=sfcreate(sd_id,sds_name,DFNT_FLOAT32,rank,dim_sizes)
         status=sfwdata(sds_id,start,stride,edges,twvi)
          write(6,*) "status sfwdata ",status
         status = sfendacc(sds_id)
         write(sds_name,'(A4)') "tozi"
          sds_id=sfcreate(sd_id,sds_name,DFNT_FLOAT32,rank,dim_sizes)
         status=sfwdata(sds_id,start,stride,edges,tozi)
          write(6,*) "status sfwdata ",status
         status = sfendacc(sds_id)
         write(sds_name,'(A6)') "tpresi"
          sds_id=sfcreate(sd_id,sds_name,DFNT_FLOAT32,rank,dim_sizes)
         status=sfwdata(sds_id,start,stride,edges,tp)
          write(6,*) "status sfwdata ",status
         status = sfendacc(sds_id)
c writing lat and lon
             write(sds_name,'(A3)') "lat"
          sds_id=sfcreate(sd_id,sds_name,DFNT_FLOAT32,rank,dim_sizes)
             status=sfwdata(sds_id,start,stride,edges,tlat)
          write(6,*) "status sfwdata ",status
             status = sfendacc(sds_id)
             write(sds_name,'(A3)') "lon"
          sds_id=sfcreate(sd_id,sds_name,DFNT_FLOAT32,rank,dim_sizes)
             status=sfwdata(sds_id,start,stride,edges,tlon)
          write(6,*) "status sfwdata ",status
             status = sfendacc(sds_id)
            write(sds_name,'(A8)') "tratiob1"
            write(6,*) "tratiob1 ",tratiob1(672,709),tratiob1(709,672)
         sds_id=sfcreate(sd_id,sds_name,DFNT_INT16,rank,dim_sizes)
          status=sfwdata(sds_id,start,stride,edges,tratiob1)
          write(6,*) "status sfwdata ",status
          status = sfendacc(sds_id)
          scalefactor=1000.
          offset=0.0
          status=sfsattr(sds_id, "scale_factor", DFNT_FLOAT32, 1, scalefactor)
          status=sfsattr(sds_id, "offset", DFNT_FLOAT32, 1, offset)
          status = sfendacc(sds_id)
         write(sds_name,'(A6)') "IPFLAG"
         sds_id=sfcreate(sd_id,sds_name,DFNT_UINT8,rank,dim_sizes)
         status=sfwdata(sds_id,start,stride,edges,ipflag)
         write(6,*) "status sfwdata ",status
          write(sds_name,'(A6)') "AERIMP"
         sds_id=sfcreate(sd_id,sds_name,DFNT_UINT8,rank,dim_sizes)
         status=sfwdata(sds_id,start,stride,edges,aerimp)
         write(6,*) "status sfwdata ",status
        status = sfendacc(sds_id)
         write(sds_name,'(A5)') "tepsa"
        sds_id=sfcreate(sd_id,sds_name,DFNT_FLOAT32,rank,dim_sizes)
         status=sfwdata(sds_id,start,stride,edges,teps)
        write(6,*) "status sfwdata ",status
         status = sfendacc(sds_id)
         endif
         
         
            
!closing the HDF file                 
         status = sfend(sd_id)          
         write(6,*) "status sfend ",status
c         
         
         
        stop
        end
        

