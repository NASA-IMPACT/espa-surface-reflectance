	   subroutine utmtodeg(utmzone,i,j,x0,y0,gsize,lat,lon)
	   IMPLICIT NONE
	   integer i,j
	   real x0,y0,x,y,gsize
	   integer utmzone,zone
	   real lat,lon
           real sa,inv_flattening,sb
	   real e2,e2cuadrada,c
	   real S,v,a,a1,a2,j2,j4,j6
	   real alfa,beta,gama
	   real bm,b,nab,senoheps
	   real epsi,eps,pi
	   real delt,ta0
	   x=x0+(j-1)*gsize
	   y=y0-(i-1)*gsize
           sa = 6378137.
           inv_flattening = 298.257223563;
           sb = sa - (sa / inv_flattening)
           e2=(((sa**2)-(sb**2))**0.5)/sb
           e2cuadrada= (e2**2)
           c=(sa**2)/sb
           x=x-500000
	   if (utmzone.lt.0) then
	   y=y-10000000
	   zone=-utmzone
	   else
	   zone=utmzone
	   endif
           S=((utmzone*6.)-183.)
	   lat=Y/(6366197.724*0.9996);
           v=(c/((1.+(e2cuadrada*(cos(lat))**2)))**0.5)*0.9996;
           a=x/v
	   a1=sin(2*lat)
	   a2=a1*(cos(lat))**2
	   j2=lat+(a1/2)
	   j4=((3.*j2)+a2)/4.
	   j6=((5.*j4)+(a2*(cos(lat))**2))/3.
	   alfa=(3./4.)*e2cuadrada
	   beta=(5./3.)*alfa**2
	   gama=(35./27.)*alfa**3;
           bm=0.9996*c*(lat-alfa*j2+beta*j4-gama*j6)
           b=(y-bm)/v;
           epsi=((e2cuadrada*a**2)/2.)*(cos(lat))**2
           eps=a*(1.-(epsi/3.))
           nab=(b*(1.-epsi))+lat
           senoheps=(exp(eps)-exp(-eps))/2
	   delt=atan(senoheps/(cos(nab)))
	   ta0=atan(cos(delt)*tan(nab))
	   pi=atan(1.)*4.
           lon=(delt*(180./pi))+s
           lat=(lat+(1+e2cuadrada*(cos(lat)**2)-(3./2.)*e2cuadrada*sin(lat)*cos(lat)*(ta0-lat))*(ta0-lat))*(180./pi)
c	   write(6,*) utmzone,x,y,lon,lat
	   return
	   end
	   

