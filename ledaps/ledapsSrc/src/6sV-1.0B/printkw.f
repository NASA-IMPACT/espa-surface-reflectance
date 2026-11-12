          program printkw
	  real c,wl,r2
	  integer iwl
	  c=0.1
	  do iwl=1,61
	  wl=0.40+.3*(iwl-1)/60
	  call morcasiwat(wl,C,R2)
	  enddo
	  stop
	  end
	  
	  
	  
	        subroutine morcasiwat(wl,C,R2)
C Spectral diffuse attenuation coefficient of Case I Waters as Predicted 
C by MOREL within the spectral range 400-700nm (1988, Journal of Geophysical 
C Research, Vol.93, No C9, pp 10749-10768)
C
C input parameters:	wl wavelength (IN MICROMETERS)
C			C  pigment concentration
C output parameter:	R2  reflectance of water
C
C According Morel,1988, we use:
C
C Kd	spectral value of the attenuation coefficient for 
C	 downwelling irradiance
C	 with: Kd=Kw+Xc*C**e
C Kw	spectral value of the diffuse attenuation coefficient 
C	 for pure oceanic water
C Xc, e	spectral coefficients to compute the diffuse attenuation 
C	 coefficient for pigment
C bb	total backscattering coefficient
C	 with: bb=0.5*bw+bbt*b
C bw	spectral value of the molecular scattering coefficient of water
C bbt,b	parameters to compute the scattering coefficients of pigments
C
C R2	reflectance of water below the surface
C	 with: R2=(0.33/u)*(bb/Kd)	where u is depending of R2
C
      real Kw,Kd
      real tKw(61),tXc(61),te(61),tbw(61)
      real wl,c,r2,xc,e,bw,bb,b,bbt,u1,r1,u2,err
      integer iwl

      data tKw/0.0209,0.0200,0.0196,0.0189,0.0183,
     & 0.0182,0.0171,0.0170,0.0168,0.0166,
     & 0.0168,0.0170,0.0173,0.0174,0.0175,
     & 0.0184,0.0194,0.0203,0.0217,0.0240,
     & 0.0271,0.0320,0.0384,0.0445,0.0490,
     & 0.0505,0.0518,0.0543,0.0568,0.0615,
     & 0.0640,0.0640,0.0717,0.0762,0.0807,
     & 0.0940,0.1070,0.1280,0.1570,0.2000,
     & 0.2530,0.2790,0.2960,0.3030,0.3100,
     & 0.3150,0.3200,0.3250,0.3300,0.3400,
     & 0.3500,0.3700,0.4050,0.4180,0.4300,
     & 0.4400,0.4500,0.4700,0.5000,0.5500,
     & 0.6500/
      data tXc/0.1100,0.1110,0.1125,0.1135,0.1126,
     & 0.1104,0.1078,0.1065,0.1041,0.0996,
     & 0.0971,0.0939,0.0896,0.0859,0.0823,
     & 0.0788,0.0746,0.0726,0.0690,0.0660,
     & 0.0636,0.0600,0.0578,0.0540,0.0498,
     & 0.0475,0.0467,0.0450,0.0440,0.0426,
     & 0.0410,0.0400,0.0390,0.0375,0.0360,
     & 0.0340,0.0330,0.0328,0.0325,0.0330,
     & 0.0340,0.0350,0.0360,0.0375,0.0385,
     & 0.0400,0.0420,0.0430,0.0440,0.0445,
     & 0.0450,0.0460,0.0475,0.0490,0.0515,
     & 0.0520,0.0505,0.0440,0.0390,0.0340,
     & 0.0300/
      data te/0.668,0.672,0.680,0.687,0.693,
     & 0.701,0.707,0.708,0.707,0.704,
     & 0.701,0.699,0.700,0.703,0.703,
     & 0.703,0.703,0.704,0.702,0.700,
     & 0.700,0.695,0.690,0.685,0.680,
     & 0.675,0.670,0.665,0.660,0.655,
     & 0.650,0.645,0.640,0.630,0.623,
     & 0.615,0.610,0.614,0.618,0.622,
     & 0.626,0.630,0.634,0.638,0.642,
     & 0.647,0.653,0.658,0.663,0.667,
     & 0.672,0.677,0.682,0.687,0.695,
     & 0.697,0.693,0.665,0.640,0.620,
     & 0.600/
      data tbw/0.0076,0.0072,0.0068,0.0064,0.0061,
     & 0.0058,0.0055,0.0052,0.0049,0.0047,
     & 0.0045,0.0043,0.0041,0.0039,0.0037,
     & 0.0036,0.0034,0.0033,0.0031,0.0030,
     & 0.0029,0.0027,0.0026,0.0025,0.0024,
     & 0.0023,0.0022,0.0022,0.0021,0.0020,
     & 0.0019,0.0018,0.0018,0.0017,0.0017,
     & 0.0016,0.0016,0.0015,0.0015,0.0014,
     & 0.0014,0.0013,0.0013,0.0012,0.0012,
     & 0.0011,0.0011,0.0010,0.0010,0.0010,
     & 0.0010,0.0009,0.0008,0.0008,0.0008,
     & 0.0007,0.0007,0.0007,0.0007,0.0007,
     & 0.0007/
      if (wl.lt.0.400.or.wl.gt.0.700)then
	R2=0.000
	goto 60
      endif

      iwl=1+nint((wl-0.400)/0.005)
      Kw=tKw(iwl)
      Xc=tXc(iwl)
      e=te(iwl)
      bw=tbw(iwl)
C
      if (abs(C).lt.0.0001)then
         bb=0.5*bw
         Kd=Kw
      else
         b=0.30*C**0.62
         bbt=0.002+0.02*(0.5-0.25*alog10(C))*0.550/wl
         bb=0.5*bw+bbt*b
         Kd=Kw+Xc*C**e
      endif

      u1=0.75
      R1=0.33*bb/u1/Kd

 50   u2=0.90*(1.-R1)/(1.+2.25*R1)
      R2=0.33*bb/u2/Kd
      err=abs((R2-R1)/R2)
      if (err.lt.0.0001)goto 60
      R1=R2
      goto 50
 60   write(6,*) wl,r2,kw,bw*0.5,xc
      return
      end
