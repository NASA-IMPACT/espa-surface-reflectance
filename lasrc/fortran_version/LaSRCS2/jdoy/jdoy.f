      PROGRAM JDOY
      CHARACTER*8 IPA1
      INTEGER JDAY,MONTH,DOY,YEAR
      call getarg(1,ipa1)                       
      IF (IPA1(1:1).eq.' ') Then
      Write(6,*) 'usage:  extracts the year and the julian day number '
      Write(6,*) 'input:  MM/DD/YY (Month,Day in Month,YEAR)'
      Write(6,*) 'output: JJJYY (Julian Day,Year)         '
      Goto 2000
      Endif
      READ(IPA1(1:2),'(I2.2)') MONTH
      READ(IPA1(4:5),'(I3)') JDAY
      READ(IPA1(7:8),'(I3)') YEAR
C     CALCULATION OF THE VARIABILITY OF THE SOLAR CONSTANT DURING THE   
C     YEAR.                                                             
C     JDAY IS THE NUMBER OF THE DAY IN THE MONTH 
      IF (MONTH.LE.2) THEN
                      J=31*(MONTH-1)+JDAY
                      goto 3
                      ENDIF
      IF (MONTH.GT.8) THEN
                      J=31*(MONTH-1)-((MONTH-2)/2)-2+JDAY
                      ELSE
                      J=31*(MONTH-1)-((MONTH-1)/2)-2+JDAY
                      ENDIF
C     IF(year.NE.0 .AND. MOD(year,4).EQ.0) J=J+1 Change to fix 2000 year pbs
      IF(MOD(year,4).EQ.0) J=J+1
    3 DOY=J
      WRITE(6,'(I3.3,I2.2)') DOY,YEAR
2000  STOP                                                            
      END                                                               
