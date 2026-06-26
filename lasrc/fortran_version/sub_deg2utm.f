	   subroutine degtoutm(utmzone,lat,lon,x0,y0,gsize,j,i)
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
	   real deltas,epsilon,latr,lonr,nu,ta
           sa = 6378137.
           sb = 6356752.314245
           e2=(((sa**2)-(sb**2))**0.5)/sb
           e2cuadrada= (e2**2)
           c=(sa**2)/sb
	   pi=atan(1.)*4.
	   latr=lat*pi/180.
	   lonr=lon*pi/180.
	   S=((utmzone*6)-183.)
           deltaS=lonr-(S*(pi/180.))
           a=cos(latr)*sin(deltaS)
           epsilon=0.5*log((1.+a)/(1.-a))
           nu=atan(tan(latr)/cos(deltaS))-latr
           v=(c/((1+(e2cuadrada*(cos(latr))**2.)))**0.5)*0.9996
           ta=(e2cuadrada/2.)*epsilon**2*(cos(latr))**2
           a1=sin(2.*latr)
           a2=a1*(cos(latr))**2.
           j2=latr+(a1/2.)
           j4=((3.*j2)+a2)/4.
           j6=((5.*j4)+(a2*(cos(latr))**2.))/3.
           alfa=(3./4.)*e2cuadrada
           beta=(5./3.)*alfa**2.
           gama=(35./27.)*alfa**3.
           Bm=0.9996*c*(latr-alfa*j2+beta*j4-gama*j6)
           x=epsilon*v*(1+(ta/3.))+500000
           y=nu*v*(1+ta)+Bm
	   i=int((y0-y)/gsize)
	   j=int((x-x0)/gsize)
	   return
	   end
	   

