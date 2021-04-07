      PROGRAM MAKEBD

*** Creates BLOCK DATA routine called initbd.f
*** that contains the boundary information for the regions

*** Obtain coordinates for vertices that form the polygons
*** that correspond to the boundaries for the regions.        
*** Region 1 is the San Andreas fault region of central California        
*** Region 2 is western CONUS
*** Region 3 is the southeastern U.S.
*** Region 4 is eastern CONUS & southern Canada
*** Region 5 is mainland Alaska & western Canada
*** Region 6 is Vancouver Island
*** Region 7 is mainland Canada
*** Region 8 is Caribbean & Central America
*** Region 9 is a place holder
*** Region 10 is a placeholder
*** Region 11 is a placeholder
*** Region 12 is a placeholder
*** Region 13 is a placeholder
*** Region 14 is the North American plate 
*** Region 15 is the Caribbean plate
*** Region 16 is the Pacific plate
*** Region 17 is the Juan de Fuca plate
*** Region 18 is the Cocos plate
*** Region 19 is the Mariana plate
*** Region 20 is the Philippine Sea plate
*** Region 21 is the South American Plate
*** Region 22 is the Nazca plate
*** Region 23 is the Panama plate
*** Region 24 is the North Andes plate

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (NMREGN = 24)
      DIMENSION NPOINT(30)

      I1 = 11
      I2 = 12

      OPEN(I1,FILE='BNDRY4.2.4.txt',STATUS='old')
      OPEN(I2, FILE='initbd.f', STATUS = 'UNKNOWN')

      WRITE(I2, 10)
   10 FORMAT(
     * '      BLOCK DATA INITBD'/
     * '      IMPLICIT DOUBLE PRECISION (A-H, O-Z)'/
     * '      IMPLICIT INTEGER*4 (I-N)'/
     * '      parameter (NMREGN = 24)'/
     * '      COMMON /BNDRY/ X(5000), Y(5000), NPOINT(30)'
     * / '      ')

      IREGN = 0
      J = 0
C
C     READ BOUNDARY DATA
C
  100 READ(I1,*,END=110)I, XLAT, XLON       
      IF(I .LE. 0 .OR. I .GT. NMREGN) THEN
         WRITE(LUOUT,101)
  101    FORMAT(' Error in boundary data')
         STOP
      ENDIF
      J = J + 1
      IF (I .ne. IREGN) THEN
         IREGN = I
         NPOINT(IREGN) = J
      ENDIF
      WRITE(I2, 200) J, J, XLAT, XLON
  200 FORMAT(6x, 'DATA  X(', I4, '), Y(',
     1       I4, ') /', F9.5, ',',
     1       F10.5, '/')
      GO TO 100

  110 CONTINUE
      NPOINT(NMREGN + 1) = J + 1
      WRITE(I2, 120)
  120 FORMAT('     '/
     1       '      DATA')

      DO 160 I = 1, NMREGN
         WRITE(I2, 150) I, NPOINT(I)
  150    FORMAT(5x,'* NPOINT(', I2, ') /', I4, '/,')
  160 CONTINUE

      I = NMREGN + 1
      WRITE(12, 165) I, NPOINT(I)
  165 FORMAT(5x,'* NPOINT(', I2, ') /', I4, '/')

      WRITE(I2, 170)
  170 FORMAT( '  ' /
     *        '      END')
      CLOSE(I1,STATUS='KEEP')
      CLOSE(I2, STATUS = 'KEEP')
      STOP  
      END
