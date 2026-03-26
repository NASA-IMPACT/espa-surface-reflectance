        program PINATUBO
        character*10 ipa1
        real xtau,n(20),k(20)
        real rad(22),dvdlnr(22)
        real cij1,cij2,rm1,rm2,rst1,rst2
	real pi,cte,tmp1,tmp2
	Integer nday
        data rad /0.050000,0.065604,0.086077,0.112939,0.148184,
     &  0.194429,0.255105,0.334716,0.439173,0.576227,0.756052,
     &  0.991996,1.301571,1.707757,2.240702,2.939966,3.857452,
     &  5.061260,6.640745,8.713145,11.432290,15.000000/
C note the optical thickness is at 500nm (but probably close enough to 550nm
C not to worry too much about it.     
        call getarg(1,ipa1)
        read(ipa1,'(I4)') nday
        if (nday.gt.1000) then
	 nday=1000
	 endif
        pi=atan(1.)*4.0
	cte=sqrt(2*pi)
	if (nday.lt.365) then
	   reff=0.17+(2.07E-03-(2.8E-06*nday))*nday
	   else
	   reff=0.17+0.38*exp(-(nday-365.)/1483.)
	   endif
	
	   
	rm1=reff
	rst1=log(1.6)
	rm1=reff/exp(5.*rst1*rst1/2.)
	rm1=exp(log(rm1)+3*rst1*rst1)
	cij1=1.0
	cij2=0.0
	do i=1,20
	n(i)=1.34
	k(i)=0.000
	enddo
	write(6,*) "22"
	do i=1,22
	tmp1=exp(-(log(rad(i))-log(rm1))**2./2./rst1/rst1)
 	tmp2=exp(-(log(rad(i))-log(rm2))**2./2./rst2/rst2)
	dvdlnr(i)=(cij1*tmp1/cte/rst1)+(cij2*tmp2/cte/rst2)
	write(6,*) rad(i),dvdlnr(i)
	enddo
	write(6,101) (n(i),i=1,20)
	write(6,102) (k(i),i=1,20)
 101    format(20(F10.5,1X))
 102    format(20(F10.7,1X))	
	stop
	end
	 
       

