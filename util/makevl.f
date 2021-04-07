      PROGRAM MAKEVL      
*** READ FILE CONTAINING GRID OF VELOCITIES & STANDARD ERRORS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (NUMGRD = 8)
      CHARACTER*10 LONMIN,LONMAX,LATMIN,LATMAX

      COMMON /CDGRID/ GRDLX(NUMGRD), GRDUX(NUMGRD), 
     1          GRDLY(NUMGRD), GRDUY(NUMGRD), 
     1          ICNTX(NUMGRD), ICNTY(NUMGRD), NBASE(NUMGRD)

      I1 = 11
      I2 = 12
      I3 = 13

      OPEN (I1, FILE = 'GRID4.2.5A.txt', STATUS = 'OLD')
      OPEN (I2, FILE = 'initvl.f', STATUS = 'UNKNOWN')
      OPEN (I3, FILE = 'Data4.2.5A.txt', STATUS = 'UNKNOWN',
     1          FORM = 'unformatted')
     
      WRITE(I2, 5)
    5 FORMAT(
     * '      BLOCK DATA INITVL'/
     * '      IMPLICIT DOUBLE PRECISION (A-H, O-Z)'/
     * '      IMPLICIT INTEGER*4 (I-N)'/
     * '      parameter (NUMGRD = 8)'/
     * '      COMMON /CDGRID/ GRDLX(NUMGRD), GRDUX(NUMGRD),'/
     * '     *   GRDLY(NUMGRD), GRDUY(NUMGRD),'/
     * '     *   ICNTX(NUMGRD), ICNTY(NUMGRD), NBASE(NUMGRD)'/
     * '      COMMON /VGRID/ B(800000)'/
     * '      COMMON /SGRID/ C(800000)'/)

      NCD = 0

      DO 100 K = 1, 8
         NBASE(K) = NCD

         READ(I1,10) LONMIN, LONMAX, ICNTX(K),
     &               LATMIN, LATMAX, ICNTY(K)
   10    FORMAT(10X, 2(2A10, I3))  

         CALL RDEG(LONMIN, GRDLX(K), 'W')
         IF(GRDLX(K) .LT. 0.D0) GRDLX(K) = GRDLX(K) + 360.D0

         CALL RDEG(LONMAX, GRDUX(K), 'W')
         IF(GRDUX(K) .LT. 0.D0) GRDUX(K) = GRDUX(K) + 360.D0

         CALL RDEG(LATMIN, GRDLY(K), 'S')
         CALL RDEG(LATMAX, GRDUY(K), 'S')
         
         if (k .gt. 7) NBASE(k) = 0

         WRITE(12, 15) K, GRDLX(K), K, GRDUX(K), K, GRDLY(K),
     *       K, GRDUY(K), K, ICNTX(K), K, ICNTY(K), K, NBASE(K)
   15    FORMAT(
     * '      DATA GRDLX(', I2, ') /', F16.11, ' /,'/
     * '     *     GRDUX(', I2, ') /', F16.11, ' /,'/
     * '     *     GRDLY(', I2, ') /', F16.11, ' /,'/
     * '     *     GRDUY(', I2, ') /', F16.11, ' /,'/
     * '     *     ICNTX(', I2, ') /', I3,     ' /,'/
     * '     *     ICNTY(', I2, ') /', I3,     ' /,'/
     * '     *     NBASE(', I2, ') /', I6,     ' /'  )

         if (k .gt. 7) go to 100
                 
         NUMPNT = (ICNTX(K) + 1) * (ICNTY(K) + 1)
         NCD = NCD + (3 * NUMPNT) 
         IF (NCD .GT. 800000) THEN
           WRITE(6, 39)
   39      FORMAT(' TOO MANY GRID NODES')
           STOP
         ENDIF                  

         DO 50 L = 1, NUMPNT 
            READ(I1,40) VN,SN,VE,SE,VU,SU
   40       FORMAT(32X,6F8.2) 
            write(13) VN,SN,VE,SE,VU,SU
   50    CONTINUE           

  100 CONTINUE

      WRITE(12, 110)
  110 FORMAT('   ' / '      END')

      CLOSE(I1,STATUS = 'KEEP')
      CLOSE(I2, STATUS = 'KEEP')
      CLOSE(I3, STATUS = 'KEEP')

      STOP   
      END
C***************************************************

      SUBROUTINE RDEG (INPUT,VAL,CNEG)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      CHARACTER*10 INPUT
      CHARACTER*1  CNEG
      INTEGER      DEG, MIN, ISEC

      DO 10 I = 1, 9
         IF (INPUT(I:I).EQ.' ') THEN
             INPUT(I:I) = '0'
         ENDIF
   10 CONTINUE

      READ (INPUT,20) DEG, MIN, ISEC

   20 FORMAT(I3,I2,I4)

      SEC = ISEC/100.D0

      VAL = (DEG + (MIN/60.D0) + (SEC/3600.D0) ) 

      IF (INPUT(10:10).EQ.CNEG) THEN
         VAL = -VAL
      ENDIF

      RETURN
      END

C***************************************************

C      INTEGER FUNCTION IUNGRD(IREGN, I, J, IVEC)

C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      IMPLICIT INTEGER*4 (I-N)
C      parameter (NUMGRD = 8)
C      COMMON /CDGRID/ GRDLX(NUMGRD), GRDUX(NUMGRD),
C    1          GRDLY(NUMGRD), GRDUY(NUMGRD),
C    1          ICNTX(NUMGRD), ICNTY(NUMGRD), NBASE(NUMGRD)

C       IUNGRD = NBASE(IREGN) +
C    1      3 * ((J - 1) * (ICNTX(IREGN) + 1) +  (I - 1)) + IVEC

C     RETURN
C     END
     
