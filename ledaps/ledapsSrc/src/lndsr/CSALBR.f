      subroutine csalbr(xtau,xalb)
!*
!Description: routine used to compute the spherical albedo of a pure
!molecular atmosphere of optical depth xtau

!Input parameters: 

!	float xtau:   the optical depth of the pure molecular
!	atmosphere (unitless).



!Output Parameters: 

!	float xalb:   the spherical albedo (unitless).

!Revision History:
!$LOG: CSLABR.f,v1.0$


!Original version: Wed Dec 13 09:33:38 EST 1995
!E. VERMOTE (eric@kratmos.gsfc.nasa.gov)


!Developer Header:

!DEVELOPPERS:  Eric F Vermote, eric@kratmos.gsfc.nasa.gov, 301-286-6232 
      
      
!References and Credits: 

!Design Notes: This is a 6S routine.
!END****************************************************************************
*/        
      real xtau,xalb,fintexp3
      xalb=(3*xtau-fintexp3(xtau)*(4+2*xtau)+2*exp(-xtau))
      xalb=xalb/(4.+3*xtau)
      return
      end
      real function fintexp3(xtau)
!*
!Description: this function returns the exponential integral of the third
! order for value xtau. It is used by CSALBR routine.

!Input parameters: 

!	float xtau:   argument (unitless).

!Output Parameters: 

!	float fintexp3: The exponential integral of order 3.

!Revision History:
!$LOG: CSLABR.f,v1.0$


!Original version: Wed Dec 13 09:33:38 EST 1995
!E. VERMOTE (eric@kratmos.gsfc.nasa.gov)


!Developer Header:

!DEVELOPPERS:  Eric F Vermote, eric@kratmos.gsfc.nasa.gov, 301-286-6232 
      
      
!References and Credits: 

!Design Notes: This is a 6S routine.
!END****************************************************************************
*/       
      real xx,xtau,fintexp1
      xx=(exp(-xtau)*(1.-xtau)+xtau*xtau*fintexp1(xtau))/2.
      fintexp3=xx
      return
      end
      real function fintexp1(xtau)
!*
!Description: this function returns the exponential integral of the third
! order for value xtau. It is used by CSALBR routine.

!Input parameters: 

!	float xtau:   argument (unitless).



!Output Parameters: 

!	float fintexp1: The exponential integral of order 1 for xtau.

!Revision History:
!$LOG: CSLABR.f,v1.0$


!Original version: Wed Dec 13 09:33:38 EST 1995
!E. VERMOTE (eric@kratmos.gsfc.nasa.gov)


!Developer Header:

!DEVELOPPERS:  Eric F Vermote, eric@kratmos.gsfc.nasa.gov, 301-286-6232 
      
      
!References and Credits: 

!Design Notes: This is a 6S routine.
!END****************************************************************************
*/       
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
