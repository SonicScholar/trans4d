      PROGRAM MAKEPS      
*** READ FILE CONTAINING GRID OF POSTSEISMIC AMPLITUDES      

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (NUMPSG = 1)
      CHARACTER*10 LONMIN,LONMAX,LATMIN,LATMAX

      COMMON /PSGRID/ PSGLX(NUMPSG), PSGUX(NUMPSG), 
     1          PSGLY(NUMPSG), PSGUY(NUMPSG), 
     1          ICNTPX(NUMPSG), ICNTPY(NUMPSG), NBASEP(NUMPSG)

      I1 = 11
      I2 = 12

      OPEN (I1, FILE = 'PSGRID3.2.txt', STATUS = 'OLD')
      OPEN (I2, FILE = 'initps.f', STATUS = 'UNKNOWN')

      WRITE(I2, 5)
    5 FORMAT(
     * '      BLOCK DATA INITPS'/
     * '      IMPLICIT DOUBLE PRECISION (A-H, O-Z)'/
     * '      IMPLICIT INTEGER*4 (I-N)'/
     * '      parameter (NUMPSG = 1)'/
     * '      COMMON /PSGRID/ PSGLX(NUMPSG), PSGUX(NUMPSG),'/
     * '     *   PSGLY(NUMPSG), PSGUY(NUMPSG),'/
     * '     *   ICNTPX(NUMPSG), ICNTPY(NUMPSG), NBASEP(NUMPSG)'/
     * '      COMMON /PGRID/ PS(18000)'/)

      NCD = 0

      DO 100 K = 1, NUMPSG
         NBASEP(K) = NCD

         READ(I1,10) LONMIN, LONMAX, ICNTPX(K),
     &               LATMIN, LATMAX, ICNTPY(K)
   10    FORMAT(10X, 2(2A10, I3))  

         CALL RDEG(LONMIN, PSGLX(K), 'W')
         IF(PSGLX(K) .LT. 0.D0) PSGLX(K) = PSGLX(K) + 360.D0

         CALL RDEG(LONMAX, PSGUX(K), 'W')
         IF(PSGUX(K) .LT. 0.D0) PSGUX(K) = PSGUX(K) + 360.D0

         CALL RDEG(LATMIN, PSGLY(K), 'S')
         CALL RDEG(LATMAX, PSGUY(K), 'S')

         WRITE(12, 15) K, PSGLX(K), K, PSGUX(K), K, PSGLY(K),
     *       K, PSGUY(K), K, ICNTPX(K), K, ICNTPY(K), K, NBASEP(K)
   15    FORMAT(
     * '      DATA PSGLX(', I2, ') /', F16.11, ' /,'/
     * '     *     PSGUX(', I2, ') /', F16.11, ' /,'/
     * '     *     PSGLY(', I2, ') /', F16.11, ' /,'/
     * '     *     PSGUY(', I2, ') /', F16.11, ' /,'/
     * '     *     ICNTPX(', I2, ') /', I3,     ' /,'/
     * '     *     ICNTPY(', I2, ') /', I3,     ' /,'/
     * '     *     NBASEP(', I2, ') /', I6,     ' /'  )

         NUMPNT = (ICNTPX(K) + 1) * (ICNTPY(K) + 1)
         NCD = NCD + (3 * NUMPNT)                   

         DO 50 L = 1, NUMPNT 
            READ(I1,40) I,J,VN,SN,VE,SE,VU,SU
   40       FORMAT(2I3,26X,6F8.2)
            INDEX = IPSGRD(K,I,J,1)
            INDEX1 = INDEX + 1
            INDEX2 = INDEX + 2
            WRITE(I2, 45) INDEX, INDEX1, INDEX2, VN, VE, VU
   45       FORMAT(
     *      '      DATA PS(', I6, '), PS(', I6, '), PS(', I6, ')'  / 
     *      5x, '1', 10x, '/ ', F8.2, ',', F8.2, ',', F8.2, ' /')
   50    CONTINUE
  100 CONTINUE

      WRITE(12, 110)
  110 FORMAT('   ' / '      END')

      CLOSE(I1,STATUS = 'KEEP')
      CLOSE(I2, STATUS = 'KEEP')

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

      INTEGER FUNCTION IPSGRD(IGRID, I, J, IVEC)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (NUMPSG = 1)
      COMMON /PSGRID/ PSGLX(NUMPSG), PSGUX(NUMPSG),
     1          PSGLY(NUMPSG), PSGUY(NUMPSG),
     1          ICNTPX(NUMPSG), ICNTPY(NUMPSG), NBASEP(NUMPSG)

      IPSGRD = NBASEP(IGRID) +
     1      3 * ((J - 1) * (ICNTPX(IGRID) + 1) +  (I - 1)) + IVEC

      RETURN
      END
     
