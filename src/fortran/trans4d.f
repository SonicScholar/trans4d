      PROGRAM TRANS4D

**************************************************************
*  NAME:       TRANS4D (Transformations in 4 Dimensions)
*
*  WRITTEN BY: Richard Snay & Chris Pearson & Jarir Saleh
*
*  PURPOSE:    Transform coordinates across time
*              and between reference frames
*
****************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      character*10 Trans4D_version   
      CHARACTER*1  OPTION
      character*5  cont

      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
      COMMON /VERSION/ Trans4D_version

C  You must change Trans4D version here if necessary

      Trans4D_version = '0.2.6'

*** Introduce variables for file id's

      LUIN = 5
*       interactive input
      LUOUT = 6
*       interactive output
      I1 = 11
*       input of velocity grid in GETVEL      
*       input of earthquake parameters in GETEQ
*       input of blue-book BFILE in DPLACE, VELOC, AND UPDATE
*       input of blue-book GFILE in TRFBB
      I2 = 12
*       output of predicted displacements in DLACE
*       output of predicted velocities in VELOC
*       output of updated blue-book BFILE in TRFBB
*       output of updated blue-book GFILE in TRFBB
*       output of updated coordinates in TRFBB
      I3 = 13
*       output of point-velocity records in VELOC
      I4 = 14
*       the above file id is no longer used
      I5 = 15
*       storage of earthquake parameters
      I6 = 16
*       the above file id is no longer used      


*** Obtain parameters defining crustal motion model
      CALL MODEL
       
*** Initialize transformation parameters between reference frames
      CALL SETTP 

*** Initialize conversion table between reference frame identifiers
      CALL SETRF

      WRITE(LUOUT,5) Trans4D_version
    5 FORMAT(
     1 ' **************************************************'/
     1 ' *  Trans4D (Transformations in 4 Dimensions)     *'/
     1 ' *  SOFTWARE VERSION ',a10                    /                            
     1 ' *                                                *'/
     1 ' *  AUTHORS:  R. Snay & C. Pearson & J. Saleh     *'/
     1 ' *            Email: rssnay@aol.com               *'/
     1 ' *                                                *'/
     1 ' **************************************************'/)
      WRITE(LUOUT,10)
   10 FORMAT( 
     1 ' This software incorporates numerical models that',/
     3 ' characterize continuous crustal motion as well as ',/
     3 ' the episodic motion associated with earthquakes.'/)
      WRITE(LUOUT,11)
   11 FORMAT(
     5 ' The User Guide contains additional information and a set'/
     5 ' of exercises to familiarize users with the software.'//)
      write (luout, 12)
   12 format(
     1 ' DISCLAIMER' //
     1 ' The Trans4D software and supporting information are '/
     2 ' currently distributed free of charge and are used by' /
     3 ' the recipient with the understanding that the providers'/
     4 ' make no warranties, expressed or implied, concerning'/
     5 ' the accuracy, completeness, reliabilty or suitability'/
     6 ' of this software, of its constituent parts, or of any'/
     7 ' supporting data.' //)
      write (luout, 13)
   13 format(
     1 ' The providers shall be under no liability whatsoever'/
     2 ' resulting from the use of this software. This software'/
     3 ' should not be relied upon as the sole basis for'/
     4 ' solving a problem whose incorrect solution could'/
     5 ' result in injury to person or property.'//
     8 ' Hit ENTER to continue.  ')
      read(luin, '(a5)',err=51,iostat=ios) cont
      if (ios /= 0) goto 51

   25 WRITE(LUOUT,26)
   26 FORMAT(' ***************************************'/
     1 ' MAIN MENU:',/
     6 '    0... Exit software.',/
     7 '    1... Estimate crustal velocities.'/
     8 '    2... Estimate crustal displacements between dates.'/
     9 '    3... Transform positions and/or observations, '/
     & '           entered in Blue Book format, across time'/
     & '           and between reference frames.'/
     & '    4... Transform positions, entered in other formats,'/  
     & '           across time and between reference frames.'/ 
     & '    5... Transform velocities between reference frames.'/)  

   30 READ(LUIN,35,err=52,iostat=ios) OPTION
      if (ios /= 0) goto 52
   35 FORMAT(A1)
      IF(OPTION .EQ. '0') THEN
            GO TO 50                   
      ELSEIF(OPTION .EQ. '2') THEN
	    CALL DPLACE
      ELSEIF(OPTION .EQ. '1') THEN
        call veloc
      ELSEIF(OPTION .EQ. '3') THEN
            call trfbb
      elseif(option .eq. '4') then
            call TRFPOS1
      elseif(option .eq. '5') then
	    call TRFVEL
      ELSE
	    WRITE(LUOUT,40)
   40       FORMAT(' Improper entry--select again  ')
	    GO TO 30
      ENDIF
      GO TO 25
   50 CONTINUE
      stop

   51 write (*,'(/)')
      write (*,*) 'You did not hit enter      : ios = ',ios
      write (*,*) "ABNORMAL TERMINATION"
      STOP

   52 write (*,'(/)')
      write (*,*) 'Problem with reading OPTION: ios = ',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      STOP

      END
*****************************************************************************
      SUBROUTINE MODEL
      
*** Obtain parameters defining crustal motion model

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /TIMREF/ ITREF
      COMMON /GRIDFILES/ IGRID, NeededGrid(8)

      A = 6.378137D06
      F = 1.D0 / 298.257222101D0
      E2 = 0.6694380022903146D-2
      AF = A / (1.D0 -F)
      EPS = F*(2.D0 - F) / ((1.D0 -F)**2)
      PI = 4.D0 * DATAN(1.D0)
      RHOSEC = (180.D0 * 3600.D0) / PI
      TWOPI = PI + PI

C*** Set default reference epoch to Jan. 1, 2010
      IYRREF = 2010
      IMOREF = 1
      IDYREF = 1
      CALL IYMDMJ (IYRREF, IMOREF, IDYREF, MJD)
      ITREF = MJD * 24 * 60

      CALL GETBDY
C***  specify initial storage grid
      IGRID = 0
C***  Sspecify which storage grid is needed for each region.

C***  Storage grid number 1 contains the information for 
C***  all regions except the Caribbean-Central America region
C***  whose information is in storage file number 2

      NeededGrid(1) = 1
      NeededGrid(2) = 1
      NeededGrid(3) = 1
      NeededGrid(4) = 1
      NeededGrid(5) = 1
      NeededGrid(6) = 1
      NeededGrid(7) = 1
      
      NeededGrid(8) = 2
      RETURN
      END

*******************************************************************
      SUBROUTINE GETBDY

*** Obtain coordinates for vertices that form the polygons
*** that define the boundaries for the regions.       
*** Region 1 is the San Andreas fault in central California         
*** Region 2 is western CONUS
*** Region 3 is southeastern U.S.
*** Region 4 is eastern CONUS & southern Canada
*** Region 5 is mainland Alaska
*** Region 6 Vancouver Island
*** Region 7 is mainland Canada
*** Region 8 is the Caribbean & Central America
*** Region 9 is a placeholder region
*** Region 10 is a placeholder region
*** Region 11 is a placeholder region
*** Region 12 is a placeholder region
*** Region 13 is a placeholder region
*** Region 14 is the North American plate 
*** Region 15 is the Caribbean plate
*** Region 16 is the Pacific plate
*** Region 17 is the Juan de Fuca plate
*** Region 18 is the Cocos plate
*** Region 19 is the Mariana plate
*** REGION 20 is the Philippine Sea plate
*** REGION 21 is the South American plate
*** REGION 22 is the Nazca plate
*** REGION 23 is the Panama plate
*** REGION 24 is the North Andes plate

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (NMREGN = 24)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /BNDRY/ X(5000), Y(5000), NPOINT(30)
 
      IEND = NPOINT(NMREGN + 1) - 1  
         DO 10 J = 1, IEND              
           X(J) = (X(J) * 3600.D0)/RHOSEC
           Y(J) = (Y(J) * 3600.D0)/RHOSEC
   10    CONTINUE
      RETURN
      END
*******************************************************************
      SUBROUTINE GETREG(X0,YKEEP,JREGN)

*** Determine the region containing a given point.

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (NMREGN = 24)
      COMMON /BNDRY/ X(5000), Y(5000), NPOINT(30)
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
      COMMON /CONST/ A, F, E2, EPS, AF, PI, TWOPI, RHOSEC

      Y0 = TWOPI - YKEEP
      IF (Y0 .lt. 0.d0) Y0 = Y0 + TWOPI
      IR = 0     
    1 IR = IR + 1
      IF(IR .GT. NMREGN) THEN
         JREGN = 0
         RETURN
      ENDIF
      IBEGIN = NPOINT(IR)
      NUMVER = NPOINT(IR + 1) - IBEGIN
      CALL POLYIN(X0,Y0,X(IBEGIN),Y(IBEGIN), NUMVER, NTEST)
      IF(NTEST .EQ. 0) GO TO 1
      JREGN = IR       

      RETURN
      END
********************************************************
      SUBROUTINE  POLYIN (X0,Y0,X,Y,N,NPC)
C     SUBROUTINE TO DETERMINE IF A POINT AT (X0,Y0) IS INSIDE OR
C     OUTSIDE OF A CLOSED FIGURE DESCRIBED BY A SEQUENCE OF CONNECTED
C     STRAIGHT LINE SEGMENTS WITH VERTICES AT X, Y.
C
C     INPUT -
C         X0, Y0    COORDINATES OF A POINT TO BE TESTED
C                    Y0 corresponds to longitude and must be a number
C                    between 0.0 and 2*PI
C         X, Y      ARRAYS CONTAINING THE VERTICES, IN ORDER, OF A
C                   CLOSED FIGURE DESCRIBED BY STRAIGHT LINE SEGMNENTS.
C                   FOR EACH 'I', THE STRAIGHT LINE FROM (XI),Y(I)) TO
C                   TO (X(I+1),Y(I+1)), IS AN EDGE OF THE FIGURE.
C         N         DIMENSION OF X AND Y, NUMBER OF VERTICES, AND NUMBER
C                   OF STRAIGHT LINE SEGMENTS IN FIGURE.
C     OUTPUT -
C         NPC       NPC=0 WHEN X0,Y0 IS OUTSIDE OF FIGURE DESCRIBED
C                   BY X,Y
C                   NPC=1 WHEN X0,Y0 IS INSIDE FIGURE
C                   NPC=2 WHEN X0,Y0 IS ON BORDER OF FIGURE
C     METHOD -
C     A COUNT IS MADE OF THE NUMBER OF TIMES THE LINE FROM (X0,Y0) TO
C     (X0,+ INFINITY) CROSSES THE BORDER OF THE FIGURE. IF THE COUNT
C     IS ODD, THE POINT IS INSIDE; IF THE COUNT IS EVEN THE POINT
C     IS OUTSIDE.
C     LIMITATIONS -
C     NONE. THE PROGRAM LOGIC IS VALID FOR ALL CLOSED FIGURES,
C     NO MATTER HOW COMPLEX.
C     ACCURACY -
C     MAINTAINS FULL ACCURACY OF INPUT COORDINATES.
C
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),Y(N)
      DATA I6/6/
      IS=0
      NPC=0
C
C     FIND STARTING POINT WHERE X(I).NE.X0
      IP=0
   10 IP=IP+1
      IF(X(IP)-X0) 15,12,16
   12 IF(IP.LE.N) GO TO 10
      WRITE(I6,6001)
 6001 FORMAT('0  POLYGON INPUT ERROR - ALL POINTS ON LINE X = X0')
      STOP
   15 IL=-1
      GO TO 20
   16 IL=1
   20 XL=X(IP)
      YL=Y(IP)
C
C     SET UP SEARCH LOOP
C
      IP1=IP+1
      IPN=IP+N
      DO 100 II=IP1,IPN
      I=II
      IF(I.GT.N) I=I-N
      IF(IL) 30,50,40
   30 IF(X(I)-X0) 90,32,34
   32 IS=-1
      GO TO 60
   34 IL=1
      GO TO 80
   40 IF(X(I)-X0) 42,44,90
   42 IL=-1
      GO TO 80
   44 IS=1
      GO TO 60
   50 IF(X(I)-X0) 52,55,54
   52 IL=-1
      IF(IS) 90,140,80
   54 IL=1
      IF(IS) 80,140,90
   55 IF(Y(I)-Y0) 57,120,58
   57 IF(YL-Y0) 90,120,120
   58 IF(YL-Y0) 120,120,90
C
   60 IL=0
      IF(Y(I)-Y0) 90,120,90
   80 IF(YL-Y0+(Y(I)-YL)*(X0-XL)/(X(I)-XL)) 90,120,85
   85 NPC=NPC+1
   90 XL=X(I)
      YL=Y(I)
  100 CONTINUE
      NPC=MOD(NPC,2)
      RETURN
  120 NPC=2
      RETURN
 140  WRITE(I6,6002)
 6002 FORMAT('0  POLYGON LOGIC ERROR - PROGRAM SHOULD NOT REACH THIS',
     .              ' POINT')
      RETURN
      END
*****************************************************************
      SUBROUTINE RADR8T (YLAT,VN,VE,VNR,VER)

C Convert horizontal velocities from mm/yr to rad/yr

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      CALL RADII (YLAT,RADMER,RADPAR)
      VNR = VN / (1000.D0 * RADMER)
      VER = VE / (1000.D0 * RADPAR)

      RETURN
      END
*****************************************************************
      SUBROUTINE COMVEL(YLAT,YLON,JREGN,VN,VE,VU,SN,SE,SU)
C
C Compute the ITRF2014 velocity at a point in mm/yr
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (NUMGRD = 8)
      parameter (NMREGN = 24)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /FILES/ LUIN, LUOUT, I1,I2,I3,I4,I5,I6
      COMMON /VGRID/ B(800000)
      COMMON /SGRID/ C(800000)
      DIMENSION WEI(2,2), VEL(2,2,3), STDEV(2,2,3)

c     WRITE (6, 1001) JREGN
c1001 FORMAT( 'JREGN = ', I6)

      IF(JREGN .GT. NUMGRD .AND. JREGN .LE. NMREGN) THEN
*** Use tectonic plate model to compute velocity relative
***    to ITRF2014
        IPLATE = JREGN - NUMGRD

*** Subtract the number of placeholder regions
        IPLATE = IPLATE - 5

        ELON = - YLON
        HT = 0.D0
        CALL TOXYZ(YLAT, ELON, HT, X, Y, Z)
        CALL PLATVL(IPLATE, X, Y, Z, VX, VY, VZ)
        VX = VX * 1000.D0
        VY = VY * 1000.D0
        VZ = VZ * 1000.D0

        CALL TOVNEU(YLAT, ELON, VX, VY, VZ, VN, VE, VU)
*** Set standard deviations of VN, VE, and VU to 5.0 mm/yr
        SN = 5.D0
        SE = 5.D0
        SU = 5.D0

      ELSEIF(JREGN .GE. 1 .AND. JREGN .LE. NUMGRD) THEN

C*** Get indices for the lower left hand corner of the grid
C*** and get the weights for the four corners
        CALL GRDWEI (YLON, YLAT, JREGN, I, J, WEI)
        
        CALL GETGRID(JREGN)

C*** Get the velocity vectors at the four corners
        CALL GRDVEC (JREGN, I, J, VEL, B)
C*** Get standard deviations at the four corners
        CALL GRDVEC (JREGN, I, J, STDEV, C)

        VN = WEI(1,1) * VEL(1,1,1) + WEI(1,2) * VEL(1,2,1)
     *     + WEI(2,1) * VEL(2,1,1) + WEI(2,2) * VEL(2,2,1)

        VE = WEI(1,1) * VEL(1,1,2) + WEI(1,2) * VEL(1,2,2)
     *     + WEI(2,1) * VEL(2,1,2) + WEI(2,2) * VEL(2,2,2)
  
        VU = WEI(1,1) * VEL(1,1,3) + WEI(1,2) * VEL(1,2,3)
     *     + WEI(2,1) * VEL(2,1,3) + WEI(2,2) * VEL(2,2,3)

        SN = WEI(1,1)*STDEV(1,1,1) + WEI(1,2)*STDEV(1,2,1)
     *     + WEI(2,1)*STDEV(2,1,1) + WEI(2,2)*STDEV(2,2,1)

        SE = WEI(1,1)*STDEV(1,1,2) + WEI(1,2)*STDEV(1,2,2)
     *     + WEI(2,1)*STDEV(2,1,2) + WEI(2,2)*STDEV(2,2,2)

        SU = WEI(1,1)*STDEV(1,1,3) + WEI(1,2)*STDEV(1,2,3)
     *     + WEI(2,1)*STDEV(2,1,3) + WEI(2,2)*STDEV(2,2,3)

C*** If the point is in one of the first 8 regions, then
c*** the velocity grids contain the ITRF2014 velocity.

      ELSE
        WRITE(LUOUT,100) JREGN
  100   FORMAT(' Improper region identifier ',I4,'in COMVEL.')
        STOP
      ENDIF
      RETURN
      END
***********************************
      subroutine getgrid(jregn)
      
      implicit double precision (a-h, o-z)
      implicit integer (i-n)
      
      parameter (numgrd = 8)
      
      common /FILES/ LUIN,LUOUT,I1,I2,I3,I4,I5,I6
      common /GRIDFILES/ IGRID, Neededgrid(8)
      common /CDGRID/ GRDLX(NUMGRD), GRDUX(NUMGRD),
     *       GRDLY(NUMGRD), GRDUY(NUMGRD),
     *       ICNTX(NUMGRD), ICNTY(NUMGRD), NBASE(NUMGRD)
      common /VGRID/ B(800000)
      common /SGRID/ C(800000)
      
      if(NeededGrid(jregn) .eq. IGRID) return
         
      if(Neededgrid(jregn) .eq. 1) then
          IGRID = 1
          open(i4, file = 'Data4.2.5A.txt',
     1             status = 'old',
     1            form = 'unformatted')
          rewind i4
          do 40 iregn = 1, 7
            do 35 i = 1, icntx(iregn) + 1
              do 30 j = 1, icnty(iregn) + 1
                read(i4) vn,sn,ve,se,vu,SU
                index = iungrd(iregn,i,j,1)
                index1 = index + 1
                index2 = index + 2
                b(index) = vn
                c(index) = SN
                b(index1) = VE
                c(index1) = SE
                b(index2) = vu
                c(index2) = SU
   30         continue
   35        continue
   40       continue
            close(i4, status = 'keep')
          
          
      endif    
          
      if (Neededgrid(jregn) .eq. 2) then 
          IGRID = 2
          open(i4,file = 'Data4.2.5B.txt',
     1             status = 'old',
     1             form = 'unformatted') 
          rewind i4
          k = 0
          do 50 i = 1, 609
             do 45 j = 1, 289
               read (i4) vn,sn,ve,se,vu,su
               index = iungrd(jregn,i,j,1)
               index1 = index + 1
               index2 = index + 2
               b(index) = vn
               c(index) = sn
               b(index1) = ve
               c(index1) = se
               b(index2) = vu
               c(index2) = su
   45         continue
   50     continue
          close(i4, status = 'keep')
      endif
      return
      end
****************************************
      SUBROUTINE PLATVL(IPLATE, X, Y, Z, VX, VY, VZ)
*** Compute the ITRF2014 velocity at point on plate = IPLATE
***    with coordinates X, Y, Z (in meters)
***    The resulting velocities--VX, VY, and VZ--will be in meters/yr
***    References 
***     Altamimi et al. 2017 = JGR (Paper on ITRF2014 plate motion)
***     Kreemer e al. 2014 = Geochem. Geophys. & Geosyst., vol 15
***     DeMets et al. 2010 = Geophysical Journal Int'l, vol 181, 
***     Snay 2003 = SALIS, Vol 63, No 1 (Paper on Frames for Pacific)
***     Bird 2003 = Geochem. Geophys. & Geosyst., (Paper on plate boundaries)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION WX(11), WY(11), WZ(11)
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6

*** IPLATE = 1 --> North America (from Altamimi et al. 2017)
***          2 --> Caribbean (from Kreemer et al. 2014)
***          3 --> Pacific (from Altamimi et al. 2017)
***          4 --> Juan de Fuca (from DeMets et al. 2010)
***          5 --> Cocos (from DeMets et al. 2010)
***          6 --> Mariana (from Snay, 2003)
***          7 --> Philippine Sea (from Kreemer et al. 2014)
***          8 --> South America (from Altamimi et al. 2017)
***          9 --> Nazca (from Altamimi et al. 2017)
***         10 --> Panama (from Kreemer et al. 2014)
***         11 --> North Andes (from Bird 2003)
      DATA WX /0.116D-9,  -0.675D-9,-1.983D-9, 
     1         6.636D-9, -10.380D-9,-0.097D-9,
     2         9.221D-9,  -1.309D-9,-1.614D-9, 
     3         2.088D-9,  -1.872D-9 /
      DATA WY /-3.365D-9, -3.826D-9, 5.076D-9,
     1         11.761D-9,-14.900D-9, 0.509D-9,
     2         -4.963D-9, -1.459D-9,-0.679D-9, 
     3        -23.037D-9, -1.285D-9 /
      DATA WZ /-0.305D-9,  2.910D-9,-10.516D-9, 
     1        -10.630D-9,  9.133D-9, -1.682D-9,
     2        -11.554D-9, -0.679D-9,  7.868D-9, 
     3          6.729D-9, -0.067D-9 /   

      IF (IPLATE .LE. 0 .OR. IPLATE .GT. 11) THEN
          WRITE (LUOUT, 1) IPLATE
    1     FORMAT(' Improper plate ID in PLATVL = ', I6)
          STOP
      ENDIF

      VX = -WZ(IPLATE) * Y + WY(IPLATE) * Z
      VY =  WZ(IPLATE) * X - WX(IPLATE) * Z
      VZ = -WY(IPLATE) * X + WX(IPLATE) * Y

*** The parameters--WX, WY, and WZ--refer to ITRF2000
*** for the Mariana Plate (Snay, 2003). Hence,
*** for this plate, VX, VY, and VZ, correspond to ITRF2000.
*** The following code converts these to ITRF2014 velocities for
*** this plate.
      IF (IPLATE .EQ. 6) THEN
         VX = VX*1000.d0
         VY = VY*1000.d0
         VZ = VZ*1000.d0
         CALL VTRANF(X, Y, Z, VX, VY, VZ, 11, 16)
         VX = VX/1000.d0
         VY = VY/1000.d0
         VZ = VZ/1000.d0
    
      ENDIF

      RETURN
      END
**********************************************************************
      SUBROUTINE PVPRNT(LATD,LATM,SLAT,LOND,LONM,SLON,VN,VE,VU)

***  Print out a point-velocity (PV) record

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /FILES/ LUIN,LUOUT,I1,I2,I3,I4,I5, I6

      LATS = IDINT(SLAT*100.D0 + 0.5D0)
      LONS = IDINT(SLON*100.D0 + 0.5D0)
      IVN  = IDINT(  VN*100.D0 + 0.5D0)
      IVE  = IDINT(  VE*100.D0 + 0.5D0)
      IVU  = IDINT(  VU*100.D0 + 0.5D0)
      JVN  = 300 
      JVE  = 300 
      JVU  = 500

      WRITE(I3,10) LATD,LATM,LATS,LOND,LONM,LONS,
     1             IVN,JVN,IVE,JVE,IVU,JVU
   10 FORMAT('PV',I3,I2.2,I4.4,'N',I3,I2.2,I4.4,'W',6I6)
      RETURN
      END
****************************************************************************
      SUBROUTINE DSDA(ISN, JSN, MIN1, MIN2, DS, DA)

** Compute change in distance and change in azimuth.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (nbbdim = 10000)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
      COMMON /BBINFO/ BLAT(nbbdim),BLON(nbbdim),BEHT(nbbdim),
     *                BVN(nbbdim),BVE(nbbdim),BVU(nbbdim),
     *                K80(nbbdim),K86(nbbdim)

      ylatf = blat(isn)
      ylatt = blat(jsn)
      ylonf = blon(isn)
      ylont = blon(jsn)
      htf = beht(isn)
      htt = beht(jsn)
      vnf = bvn(isn)
      vnt = bvn(jsn)
      vef = bve(isn)
      vet = bve(jsn)
      vuf = bvu(isn)
      vut = bvu(jsn)

      CALL COMPSN(YLATF1,YLONF1,HTF1,YLATF,YLONF,HTF,
     1            MIN1, VNF, VEF, VUF)
      CALL COMPSN(YLATF2,YLONF2,HTF2,YLATF,YLONF,HTF,
     1            MIN2, VNF, VEF, VUF)
      CALL COMPSN(YLATT1,YLONT1,HTT1,YLATT,YLONT,HTT,
     1            MIN1, VNT, VET, VUT)
      CALL COMPSN(YLATT2,YLONT2,HTT2,YLATT,YLONT,HTT,
     1            MIN2, VNT, VET, VUT)

      CALL HELINV(YLATF1,YLONF1,YLATT1,YLONT1,FAZ1,BAZ1,S1)
      CALL HELINV(YLATF2,YLONF2,YLATT2,YLONT2,FAZ2,BAZ2,S2)
      DS = S2 - S1
      DA = FAZ2 - FAZ1
      RETURN   
      END
*********************************************************************
      SUBROUTINE HELINV(GLAT1,GLON1,GLAT2,GLON2,FAZ,BAZ,S)               
C                                                                       
C *** SOLUTION OF THE GEODETIC INVERSE PROBLEM AFTER T.VINCENTY.       
C *** MODIFIED RAINSFORD WITH HELMERT ELLIPTICAL TERMS.                 
C *** EFFECTIVE IN ANY AZIMUTH AND AT ANY DISTANCE SHORT OF ANTIPODAL. 
C *** STANDPOINT/FOREPOINT MUST NOT BE THE GEOGRAPHIC POLE .           
C
C   INPUT
C      GLAT1 = Latitude of from point (radians, positive north)
C      GLON1 = Longitude of from point (radians, positive west)
C      GLAT2 = Latitude of to point (radians, positive north)
C      GLON2 = Longitude of to point (radians, positive west)
C                        
C   OUTPUT
C      FAZ = Foward azimuth (radians, clockwise from north)
C      BAZ = Back azimuth   (radians, clockwise from north)
C      S   = Distance (meters)
C                                               
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      IMPLICIT INTEGER*4 (I-N)
      COMMON/CONST/A,F,EPS2,EPS,AF,PI,TWOPI,RHOSEC                      
      DATA TOL/0.5D-14/                                                 
      R = 1.0D0-F                                                       
      TU1 = R*DSIN(GLAT1)/DCOS(GLAT1)                                   
      TU2 = R*DSIN(GLAT2)/DCOS(GLAT2)                                   
      CU1 = 1.0D0/DSQRT(TU1*TU1+1.0D0)                                  
      SU1 = CU1*TU1                                                     
      CU2 = 1.0D0/DSQRT(TU2*TU2+1.0D0)                                  
      S = CU1*CU2                                                       
      BAZ = S*TU2                                                       
      FAZ = BAZ*TU1                                                     
      X = GLON1-GLON2                                                   
  100 SX = DSIN(X)                                                      
      CX = DCOS(X)                                                      
      TU1 = CU2*SX                                                      
      TU2 = SU1*CU2*CX-BAZ                                              
      SY = DSQRT(TU1*TU1+TU2*TU2)                                       
      CY = S*CX+FAZ                                                     
      Y = DATAN2(SY,CY)                                                 
      SA = S*SX/SY                                                      
      C2A = -SA*SA+1.0D0                                                
      CZ = FAZ+FAZ                                                      
      IF(C2A .GT. 0.0D0)CZ = -CZ/C2A+CY                                 
      E = CZ*CZ*2.0D0-1.0D0                                             
      C = ((-3.0D0*C2A+4.0D0)*F+4.0D0)*C2A*F/16.0D0                     
      D = X                                                             
      X = ((E*CY*C+CZ)*SY*C+Y)*SA                                       
      X = (1.0D0-C)*X*F+GLON1-GLON2                                     
      IF(DABS(D-X) .GT. TOL)GO TO 100                                   
      FAZ = DATAN2(-TU1,TU2)                                            
      BAZ = DATAN2(CU1*SX,BAZ*CX-SU1*CU2)     
      FAZ = FAZ + PI
      BAZ = BAZ + PI                      
      IF(FAZ .LT. 0.0D0)FAZ = FAZ+TWOPI                             
      IF(BAZ .LT. 0.0D0)BAZ = BAZ+TWOPI   
      IF(FAZ .GT. TWOPI)FAZ = FAZ - TWOPI
      IF(BAZ .GT. TWOPI)BAZ = BAZ - TWOPI                            
      X = DSQRT((1.0D0/R/R-1.0D0)*C2A+1.0D0)+1.0D0                      
      X = (X-2.0D0)/X                                                   
      C = 1.0D0-X                                                       
      C = (X*X/4.0D0+1.0D0)/C                                           
      D = (0.375D0*X*X-1.0D0)*X                                         
      X = E*CY                                                          
      S = 1.0D0-E-E                                                     
      S = ((((SY*SY*4.0D0-3.0D0)*S*CZ*D/6.0D0-X)*D/4.0D0+CZ)*SY*D+Y)    
     1    *C*A*R                                                        
      RETURN                                                            
      END
*************************************************************************

      SUBROUTINE TRFDAT(CARD,DATE,IREC12,IYEAR1,IYEAR2,MINS)

C Convert blue-book date to time in minutes

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*80 CARD
      CHARACTER*6 DATE
      COMMON /FILES/ LUIN,LUOUT, I1, I2, I3, I4, I5, I6

      IF (IREC12 .EQ. 0) THEN
         WRITE (LUOUT, 5)
    5    FORMAT(' ABORT: The blue-book needs a valid *12* record.')
         STOP
      ENDIF

      READ (DATE, 10,err=50,iostat=ios) IYEAR, MONTH, IDAY
      if (ios /= 0) goto 50
   10 FORMAT (3I2)

      IF ( IYEAR1 .LE. (1900 + IYEAR) .AND.
     1     (1900 + IYEAR) .LE. IYEAR2) THEN
         IYEAR = 1900 + IYEAR
      ELSEIF ( IYEAR1 .LE. (2000 + IYEAR) .AND.
     1     (2000 + IYEAR) .LE. IYEAR2) THEN
         IYEAR = 2000 + IYEAR
      ELSEIF ( IYEAR1 .LE. (1800 + IYEAR) .AND.
     1     (1800 + IYEAR) .LE. IYEAR2) THEN
         IYEAR = 1800 + IYEAR
      ELSE
         WRITE (LUOUT, 20) CARD
   20    FORMAT(' ABORT: The following record has a date'/
     1          ' which is inconsistent with the *12* record'/
     1          3x, A80)
      ENDIF

      IF (IYEAR .LE. 1906) THEN
         WRITE (LUOUT, 30)
   30    FORMAT(' ***WARNING***'/
     1   ' The blue-book file contains an observation that'/
     1   ' predates 1906.  The TRANS4D model may not be valid'/
     1   ' and the computed corrections may be erroneous.')
      ENDIF
   
      IF (IDAY .EQ. 0) IDAY = 15

      IF (MONTH .EQ. 0) THEN
         MONTH = 7
         IDAY = 1
      ENDIF

C     CALL TOTIME(IYEAR,MONTH,IDAY, MINS)
      CALL IYMDMJ(IYEAR,MONTH,IDAY, MJD)
      MINS = MJD * 24 * 60
      RETURN

  50  write (*,'(/)') 
      write (*,*) 'wrong IYEAR, IMONTH, IDAY:ios =',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop


      END
*********************************************************
      SUBROUTINE TNFDAT(DATE,MINS)

C Convert blue-book date to time in years

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*8 DATE
      COMMON /FILES/ LUIN,LUOUT, I1, I2, I3, I4, I5, I6
      READ(DATE,10) IYEAR,MONTH,IDAY
   10 FORMAT(I4,I2,I2)
      IF(IYEAR .LE. 1906) THEN
	  WRITE(LUOUT,20)
   20     FORMAT(' ***WARNING***'/
     1    ' The blue-book file contains an observation that'/
     2    ' predates 1906.  The TRANS4D model is not valid and the'/
     3    ' computed correction may be erroneous.')
      ENDIF
      IF(IDAY .EQ. 0) IDAY = 15
      IF(MONTH .EQ. 0) THEN
         MONTH = 7
         IDAY = 1
      ENDIF
C     CALL TOTIME(IYEAR,MONTH,IDAY,MINS)
      CALL IYMDMJ(IYEAR,MONTH,IDAY,MJD)
      MINS = MJD * 24 * 60
      RETURN
      END

C************************************************************************
      SUBROUTINE TOXYZ(glat,glon,eht,x,y,z)
 
*** compute x,y,z
*** ref p.17 geometric geodesy notes vol 1, osu, rapp
 
      implicit double precision(a-h,o-z)
      common/CONST/ a,f,e2,ep2,af,pi,twopi,rhosec
 
      slat=dsin(glat)
      clat=dcos(glat)
      w=dsqrt(1.d0-e2*slat*slat)
      en=a/w
 
      x=(en+eht)*clat*dcos(glon)
      y=(en+eht)*clat*dsin(glon)
      z=(en*(1.d0-e2)+eht)*slat
 
      return
      end
C************************************************************************
      logical function FRMXYZ(x,y,z,glat,glon,eht)
 
*** convert x,y,z into geodetic lat, lon, and ellip. ht
*** ref: eq a.4b, p. 132, appendix a, osu #370
*** ref: geom geod notes gs 658, rapp
 
      implicit double precision(a-h,o-z)
      parameter(maxint=10,tol=1.d-13)
      common/CONST/ a,f,e2,ep2,af,pi,twopi,rhosec
 
      ae2=a*e2
 
*** compute initial estimate of reduced latitude  (eht=0)
 
      p=dsqrt(x*x+y*y)
      icount=0
      tgla=z/p/(1.d0-e2)
 
*** iterate to convergence, or to max # iterations
 
    1 if(icount.le.maxint) then
        tglax=tgla
        tgla=z/(p-(ae2/dsqrt(1.d0+(1.d0-e2)*tgla*tgla)))
        icount=icount+1
        if(dabs(tgla-tglax).gt.tol) go to 1
 
*** convergence achieved
 
        frmxyz=.true.
        glat=datan(tgla)
        slat=dsin(glat)
        clat=dcos(glat)
        glon=datan2(y,x)
        w=dsqrt(1.d0-e2*slat*slat)
        en=a/w
        if(dabs(glat).le.0.7854d0) then
          eht=p/clat-en
        else
          eht=z/slat-en+e2*en
        endif
        glon=datan2(y,x)
 
*** too many iterations
 
      else
        frmxyz=.false.
        glat=0.d0
        glon=0.d0
        eht=0.d0
      endif
 
      return
      end
*********************************************************************
      SUBROUTINE RADII(YLAT,RADMER,RADPAR)
C
C  Computes the radius of curvature in the meridian
C  and the radius of curvature in a parallel of latitude
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COSLAT = DCOS(YLAT)
      DENOM = DSQRT(1.D0 + EPS*COSLAT*COSLAT)
      RADMER = AF/(DENOM**3)
      RADPAR = AF*COSLAT/DENOM
      RETURN
      END
*********************************************************************
      SUBROUTINE GETGRD(NAMEG,MINLAT,MAXLAT,ILAT,MINLON,MAXLON,ILON)

*** Interactively obtain specifications for a grid.

      IMPLICIT INTEGER*4 (I-N)
      COMMON /FILES/ LUIN,LUOUT, I1, I2, I3, I4, I5, I6
      CHARACTER*10 NAMEG
	  WRITE(LUOUT,200)
  200     FORMAT(' Enter name for grid (10 character max).  ')
	  READ(LUIN,210,err=50,iostat=ios) NAMEG
          if (ios /= 0) goto 50
  210     FORMAT(A10)
	  WRITE(LUOUT,220)
  220     FORMAT(' Enter minimum latitude for grid in deg-min-sec'/
     1           ' in free format (integer values only).  ')
	  READ(LUIN,*,err=51,iostat=ios) ID1, IM1, IS1
          if (ios /= 0) goto 51
	  WRITE(LUOUT,230)
  230     FORMAT(' Enter maximum latitude in same format.  ')
	  READ(LUIN,*,err=52,iostat=ios) ID2, IM2, IS2
          if (ios /= 0) goto 52
	  WRITE(LUOUT,240)
  240     FORMAT(' Enter latitude increment in seconds ',
     1           ' (integer value).  ')
	  READ(LUIN,*,err=53,iostat=ios) ILAT
          if (ios /= 0) goto 53
	  WRITE(LUOUT,250)
  250     FORMAT(' Enter minimum longitude with positive being west.  ')
	  READ(LUIN,*,err=54,iostat=ios) JD1, JM1, JS1
          if (ios /= 0) goto 54
          WRITE(LUOUT,260)
  260     FORMAT(' Enter maximum longitude.  ')
	  READ(LUIN,*,err=55,iostat=ios) JD2, JM2, JS2
          if (ios /= 0) goto 55
	  WRITE(LUOUT,270)
  270     FORMAT(' Enter longitude increment in seconds ',
     1           ' (integer value).  ')
	  READ(LUIN,*,err=56,iostat=ios) ILON
          if (ios /= 0) goto 56
	  MINLAT = 3600*ID1 + 60*IM1 + IS1
	  MAXLAT = 3600*ID2 + 60*IM2 + IS2
	  MINLON = 3600*JD1 + 60*JM1 + JS1
	  MAXLON = 3600*JD2 + 60*JM2 + JS2
	  RETURN

  50      write (*,'(/)') 
          write (*,*) "Wrong grid name in GETGRD:ios=",ios
          write (*,*) "ABNORMAL TERMINATION"
          write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
          stop

  51      write (*,'(/)') 
          write (*,*) "Wrong min lat in GETGRD:ios=",ios
          write (*,*) "ABNORMAL TERMINATION"
          write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
          stop

  52      write (*,'(/)') 
          write (*,*) "Wrong max lat in GETGRD:ios=",ios
          write (*,*) "ABNORMAL TERMINATION"
          write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
          stop

  53      write (*,'(/)') 
          write (*,*) "Wrong ILAT in GETGRD:ios=",ios
          write (*,*) "ABNORMAL TERMINATION"
          write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
          stop

  54      write (*,'(/)') 
          write (*,*) "Wrong min lon in GETGRD:ios=",ios
          write (*,*) "ABNORMAL TERMINATION"
          write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
          stop

  55      write (*,'(/)') 
          write (*,*) "Wrong max lon in GETGRD:ios=",ios
          write (*,*) "ABNORMAL TERMINATION"
          write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
          stop

  56      write (*,'(/)') 
          write (*,*) "Wrong ILON in GETGRD:ios=",ios
          write (*,*) "ABNORMAL TERMINATION"
          write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
          stop

	  END
*********************************************************************
      SUBROUTINE GETLYN(NAMEG,XLAT,XLON,FAZ,BAZ,XMIN,XMAX,XINC)

*** Interactively obtain the specifications for a line.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /FILES/ LUIN,LUOUT,I1,I2,I3,I4,I5,I6
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      CHARACTER*10 NAMEG

      WRITE(LUOUT,100)
  100 FORMAT(' Enter name for line (10 character max.)  ')
      READ(LUIN,110,err=50,iostat=ios) NAMEG
      if (ios /= 0) goto 50
  110 FORMAT(A10)
      WRITE(LUOUT,120)
  120 Format(' Specify the latitude for the origin of the line'/
     1       ' in deg-min-sec in free format. For example'/
     2       '         35 17 28.3'/
     3       ' Positive is north.  ')
      READ(LUIN,*,err=51,iostat=ios) LATD,LATM,SLAT
      if (ios /= 0) goto 51
      XLAT = (DBLE(3600*LATD + 60*LATM)+SLAT)/RHOSEC
      WRITE(LUOUT,130)
  130 FORMAT(' Specify the longitude for the origin of the line'/
     1       ' in free format with west being positive.  ')
      READ(LUIN,*,err=52,iostat=ios) LOND,LONM,SLON
      if (ios /= 0) goto 52
      XLON = (DBLE(3600*LOND+60*LONM)+SLON)/RHOSEC
      WRITE(LUOUT,140)
  140 FORMAT(' Specify the orientation of the line clockwise from'/
     1       ' north in decimal degrees, eg., 43.7.  ')
      READ(LUIN,*,err=53,iostat=ios) FAZ
      if (ios /= 0) goto 53
      FAZ = 3600.D0*FAZ/RHOSEC
      BAZ = FAZ + PI
      IF(BAZ .GT. TWOPI) BAZ = BAZ - TWOPI
      WRITE(LUOUT,150)
  150 FORMAT(' Specify minimum and maximum distance from origin in'/
     1       ' meters, eg., -40000. 30000. '/
     2       ' NOTE: negative distance corresponds to distance in'/
     3       '       the opposite direction.  ')
      READ(LUIN,*,err=54,iostat=ios) XMIN,XMAX
      if (ios /= 0) goto 54
      WRITE(LUOUT,160)
  160 FORMAT(' Specify distance increment in meters.  ')
      READ(LUIN,*,err=55,iostat=ios) XINC
      if (ios /= 0) goto 55
      RETURN

  50  write (*,'(/)') 
      write (*,*) 'Wrong name of line in GETLYN: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  51  write (*,'(/)') 
      write (*,*) 'Wrong LATD,LATM,SLAT in GETLYN: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  52  write (*,'(/)') 
      write (*,*) 'Wrong LOND,LONM,SLON in GETLYN: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  53  write (*,'(/)') 
      write (*,*) 'Wrong FAZ in GETLYN: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  54  write (*,'(/)') 
      write (*,*) 'Wrong min and max D in GETLYN: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  55  write (*,'(/)') 
      write (*,*) 'Wrong increment in GETLYN: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

      END
*******************************************************************
      SUBROUTINE DIRCT1(GLAT1,GLON1,GLAT2,GLON2,FAZ,BAZ,S)
C
C *** SOLUTION OF THE GEODETIC DIRECT PROBLEM AFTER T.VINCENTY
C *** MODIFIED RAINSFORD'S METHOD WITH HELMERT'S ELLIPTICAL TERMS
C *** EFFECTIVE IN ANY AZIMUTH AND AT ANY DISTANCE SHORT OF ANTIPODAL
C
C *** A IS THE SEMI-MAJOR AXIS OF THE REFERENCE ELLIPSOID
C *** FINV IS THE FLATTENING OF THE REFERENCE ELLIPSOID
C *** LATITUDES AND LONGITUDES IN RADIANS POSITIVE NORTH AND EAST
C *** AZIMUTHS IN RADIANS CLOCKWISE FROM NORTH
C *** GEODESIC DISTANCE S ASSUMED IN UNITS OF SEMI-MAJOR AXIS A
C
C *** PROGRAMMED FOR CDC-6600 BY LCDR L.PFEIFER NGS ROCKVILLE MD 20FEB75
C *** MODIFIED FOR SYSTEM 360 BY JOHN G GERGEN NGS ROCKVILLE MD 750608
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /CONST/ A,FINV,E2,EPI,AF,PI,TWOPI,RHOSEC
      DATA EPS/0.5D-13/
      R=1.D0-FINV
      TU=R*DSIN(GLAT1)/DCOS(GLAT1)
      SF=DSIN(FAZ)
      CF=DCOS(FAZ)
      BAZ=0.D0
      IF(CF.NE.0.D0) BAZ=DATAN2(TU,CF)*2.D0
      CU=1.D0/DSQRT(TU*TU+1.D0)
      SU=TU*CU
      SA=CU*SF
      C2A=-SA*SA+1.D0
      X=DSQRT((1.D0/R/R-1.D0)*C2A+1.D0)+1.D0
      X=(X-2.D0)/X
      C=1.D0-X
      C=(X*X/4.D0+1.D0)/C
      D=(0.375D0*X*X-1.D0)*X
      TU=S/R/A/C
      Y=TU
  100 SY=DSIN(Y)
      CY=DCOS(Y)
      CZ=DCOS(BAZ+Y)
      E=CZ*CZ*2.D0-1.D0
      C=Y
      X=E*CY
      Y=E+E-1.D0
      Y=(((SY*SY*4.D0-3.D0)*Y*CZ*D/6.D0+X)*D/4.D0-CZ)*SY*D+TU
      IF(DABS(Y-C).GT.EPS)GO TO 100
      BAZ=CU*CY*CF-SU*SY
      C=R*DSQRT(SA*SA+BAZ*BAZ)
      D=SU*CY+CU*SY*CF
      GLAT2=DATAN2(D,C)
      C=CU*CY-SU*SY*CF
      X=DATAN2(SY*SF,C)
      C=((-3.D0*C2A+4.D0)*FINV+4.D0)*C2A*FINV/16.D0
      D=((E*CY*C+CZ)*SY*C+Y)*SA
      GLON2=GLON1+X-(1.D0-C)*D*FINV
      IF (GLON2.GE.TWOPI) GLON2=GLON2-TWOPI
      IF(GLON2.LT.0.D0) GLON2=GLON2+TWOPI
      BAZ=DATAN2(SA,BAZ)+PI
      IF (BAZ.GE.TWOPI) BAZ=BAZ-TWOPI
      IF (BAZ.LT.0.D0) BAZ=BAZ+TWOPI
      RETURN
      END
***********************************************************
      SUBROUTINE TOCHAR(ORIG,CHAR14)

*** Convert double precision real number to character*14

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*14 CHAR14

      ORIG = ORIG*10000.D0
      WRITE(CHAR14,10) ORIG
   10 FORMAT(F14.0)
      RETURN
      END
*********************************************************

      SUBROUTINE DDXYZ(ISN, JSN, MIN1, MIN2, 
     1                 DDX, DDY, DDZ)

** Compute change in DX,DY,DZ-vector from time MIN1 to MIN2.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
c     character      PIDs*6
      parameter (nbbdim = 10000)
c     COMMON /ARRAYS/ HT(nbbdim), LOC(nbbdim),PIDs(nbbdim)
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
      COMMON /BBINFO/ BLAT(nbbdim),BLON(nbbdim),BEHT(nbbdim),
     *                BVN(nbbdim),BVE(nbbdim),BVU(nbbdim),
     *                K80(nbbdim),K86(nbbdim)

      glati = blat(isn)
      glatj = blat(jsn)
      gloni = blon(isn)
      glonj = blon(jsn)
      ehti = beht(isn)
      ehtj = beht(jsn)
      vni = bvn(isn)
      vnj = bvn(jsn)
      vei = bve(isn)
      vej = bve(jsn)
      vui = bvu(isn)
      vuj = bvu(jsn)

c     READ(I4,REC=LOC(ISN),err=50,iostat=ios) GLATI,GLONI,VNI,VEI,VUI
c     if (ios /= 0) goto 50
c     READ(I4,REC=LOC(JSN),err=51,iostat=ios) GLATJ,GLONJ,VNJ,VEJ,VUJ
c     if (ios /= 0) goto 51

      CALL COMPSN(YLATI1,YLONI1,HTI1,GLATI,GLONI,ehti,
     1            MIN1,VNI, VEI, VUI)
      CALL COMPSN(YLATI2,YLONI2,HTI2,GLATI,GLONI,ehti,
     1            MIN2, VNI, VEI, VUI)
      CALL COMPSN(YLATJ1,YLONJ1,HTJ1,GLATJ,GLONJ,ehtj,
     1            MIN1, VNJ, VEJ, VUJ)
      CALL COMPSN(YLATJ2,YLONJ2,HTJ2,GLATJ,GLONJ,ehtj,
     1            MIN2, VNJ, VEJ, VUJ)

      XLONI1 = -YLONI1
      XLONI2 = -YLONI2
      XLONJ1 = -YLONJ1
      XLONJ2 = -YLONJ2

      CALL TOXYZ(YLATI1,XLONI1,HTI1,XI1,YI1,ZI1)
      CALL TOXYZ(YLATI2,XLONI2,HTI2,XI2,YI2,ZI2)
      CALL TOXYZ(YLATJ1,XLONJ1,HTJ1,XJ1,YJ1,ZJ1)
      CALL TOXYZ(YLATJ2,XLONJ2,HTJ2,XJ2,YJ2,ZJ2)

      DDX = (XJ2 - XI2) - (XJ1 - XI1)
      DDY = (YJ2 - YI2) - (YJ1 - YI1)
      DDZ = (ZJ2 - ZI2) - (ZJ1 - ZI1)

      RETURN
      END
*************************************************************
      SUBROUTINE DISLOC (YLAT,YLON,STRIKE,HL,EQLAT,EQLON,
     &          SS,DS,DIP,DEPTH,WIDTH,DNORTH,DWEST,DUP)

*** Compute 3-dimensional earthquake displacement at point
*** using dislocation theory
*
*   INPUT:
*         YLAT = Latitude in radians (positive north)
*         YLON = Longitude in radians (positive west)
*         STRIKE = strike in radians clockwise from north such
*                  that the direction of dip is pi/2 radians
*                  counterclockwise from the direction of strike
*         HL   = Half-length in meters
*         EQLAT = Latitude in radians of midpoint of the
*                 rectangle's upper edge (positive north)
*         EQLON = Longitude in radians of midpoint of the
*                 rectangle's upper edge (positive west)
*         SS = strike slip in meters (positive = right lateral)
*         DS = dip slip  in meters (positive = normal faulting)
*         DIP = dip in radians
*         DEPTH = Vertical depth of rectangle's upper edge
*                 in meters
*         WIDTH = width of rectangle in meters
*
*   OUTPUT:
*         DNORTH = northward displacement in radians
*         DWEST = westward displacement in radians
*         DUP = upward displacement in meters
************** 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

*** Compute radii of curvature at fault center
      CALL RADII (EQLAT, RMER, RPAR)

*** Compute planar coordinates in meters
      DLAT = (YLAT - EQLAT) * RMER
      DLON = (YLON - EQLON) * RPAR
      COSSTR = DCOS(STRIKE)
      SINSTR = DSIN(STRIKE)
      X1 = COSSTR*DLAT - SINSTR*DLON
      X2 = SINSTR*DLAT + COSSTR*DLON

*** Compute displacements in fault-oriented coordinates
      CALL OKADA(X1,X2,HL,DEPTH,WIDTH,DIP,U1SS,U2SS,
     &      U3SS,U1DS,U2DS,U3DS)
      U1 = U1SS*SS + U1DS*DS
      U2 = U2SS*SS + U2DS*DS
      DUP = U3SS*SS + U3DS*DS

*** Convert horizontal displacements to radians
*** in north-west coordinate system
      DNORTH = ( COSSTR*U1 + SINSTR*U2) / RMER
      DWEST  = (-SINSTR*U1 + COSSTR*U2) / RPAR

      RETURN
      END
****************************************************
      SUBROUTINE OKADA(X1,X2,XL,DU,W,DIP,
     1     U1SS,U2SS,U3SS,U1DS,U2DS,U3DS)

************************************************************
*  This subroutine computes displacements at the point X1,X2
*  on the Earth's surface due to 1.0 meter of right-lateral 
*  strike slip (SS) and 1.0 meter of normal dip slip (DS) 
*  along a rectangular fault.
*
*  The rectangular fault dips in the direction of the positive
*  X2-axis.  The rectangle's strike parallels the X1-axis.
*  With the X3-axis directed upward out of the Earth, the X1-,
*  X2-, and X3-axes form a right-handed system.
*
*  The equations of dislocation theory are employed whereby
*  Earth is represented an a homogeneous, isotropic half-space
*  with a Poisson ratio of PNU.
*
*  REFERENCE: Okada, Y., Surface deformation due to shear and
*    tensile faults in a half-space, Bulletin of the 
*    Seismological Society of America, vol. 75, pp. 1135-1154 (1985)
*
*  The X3 = 0 plane corresponds to the Earth's surface. The plane's
*  origin is located directly above the midpoint of the rectangle's
*  upper edge.
*
*  INPUT:
*    X1,X2 - Location in meters
*    XL    - Rectangle's half-length in meters
*    DU    - Vertical depth to rectangle's upper edge in meters
*            (always positive or zero)
*    W     - Rectangle's width in meters
*    DIP   - Rectangle's dip in radians (always between 0 and PI/2
*
*  OUTPUT
*    U1SS  - Displacement in X1-direction due to 1.0 meters
*            of right-lateral strike slip
*    U2SS  - Displacement in X2-direction due to 1.0 meters
*            of right-lateral strike slip
*    U3SS  - Displacement in X3-direction due to 1.0 meters
*            of right-lateral strike slip
*    U1DS  - Displacement in X1-direction due to 1.0 meters
*            of normal dip slip
*    U2DS  - Displacement in X2-direction due to 1.0 meters
*            of normal dip slip
*    U3DS  - Displacement in X3-direction due to 1.0 meters
*            of normal dip slip
*******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      LOGICAL VERT   
      PI = 3.141593D0
      TWOPI = PI + PI
      PNU = 0.25D0
      RATIO = 1.D0 - 2.D0*PNU

      IF(DABS(PI/2.D0 - DIP) .LT. .01D0)THEN
               DIPK = -PI/2.D0
               VERT = .TRUE.
      ELSE
               DIPK = -DIP
               VERT = .FALSE.
      ENDIF

      SDIP = DSIN(DIPK)
      CDIP = DCOS(DIPK)
      P = X2*CDIP + DU*SDIP
      Q = X2*SDIP - DU*CDIP

      PSI = X1 + XL
      ETA = P
      CALL OKADAW(PSI,ETA,Q,SDIP,CDIP,RATIO,TWOPI,
     1    VERT,U1SS,U2SS,U3SS,U1DS,U2DS,U3DS)

      PSI = X1 + XL
      ETA = P - W
      CALL OKADAW(PSI,ETA,Q,SDIP,CDIP,RATIO,TWOPI,
     1     VERT,C1SS,C2SS,C3SS,C1DS,C2DS,C3DS)
      U1SS = U1SS - C1SS
      U2SS = U2SS - C2SS
      U3SS = U3SS - C3SS
      U1DS = U1DS - C1DS 
      U2DS = U2DS - C2DS
      U3DS = U3DS - C3DS

      PSI = X1 - XL
      ETA = P
      CALL OKADAW(PSI,ETA,Q,SDIP,CDIP,RATIO,TWOPI,
     1     VERT,C1SS,C2SS,C3SS,C1DS,C2DS,C3DS)
      U1SS = U1SS - C1SS
      U2SS = U2SS - C2SS
      U3SS = U3SS - C3SS
      U1DS = U1DS - C1DS
      U2DS = U2DS - C2DS
      U3DS = U3DS - C3DS

      PSI = X1 - XL
      ETA = P - W
      CALL OKADAW(PSI,ETA,Q,SDIP,CDIP,RATIO,TWOPI,
     1     VERT,C1SS,C2SS,C3SS,C1DS,C2DS,C3DS)
      U1SS = U1SS + C1SS
      U2SS = U2SS + C2SS
      U3SS = U3SS + C3SS
      U1DS = U1DS + C1DS
      U2DS = U2DS + C2DS
      U3DS = U3DS + C3DS
      RETURN
      END
*************************************************************
      SUBROUTINE OKADAW(PSI,ETA,Q,SDIP,CDIP,RATIO,TWOPI,
     1           VERT,U1SS,U2SS,U3SS,U1DS,U2DS,U3DS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      LOGICAL VERT   

      YBAR = ETA*CDIP + Q*SDIP
      DBAR = ETA*SDIP - Q*CDIP
      R = DSQRT(PSI*PSI + ETA*ETA + Q*Q)
      X = DSQRT(PSI*PSI + Q*Q)
      IF(DABS(Q) .LE. 0.1d0) THEN
         TERM = 0.D0
      ELSE
         TERM = DATAN(PSI*ETA/(Q*R))
      ENDIF

      IF(VERT) THEN
         F5 = -RATIO*PSI*SDIP/(R + DBAR)
         F4 = -RATIO*Q/(R + DBAR)
         F3 = 0.5D0*RATIO*(ETA/(R + DBAR)
     1          + YBAR*Q/((R + DBAR)*(R + DBAR))
     2          - DLOG(R + ETA))
         F1 = -0.5D0*RATIO*PSI*Q/
     1        ((R + DBAR)*(R + DBAR))
      ELSE
         IF(DABS(PSI) .LE. 0.1D0) then
            F5 = 0.d0
         ELSE
            F5 = 2.D0*RATIO*
     1      DATAN((ETA*(X+Q*CDIP)+X*(R+X)*SDIP)/(PSI*(R+X)*CDIP))
     2          /CDIP
         ENDIF
         F4 = RATIO*(DLOG(R+DBAR)-SDIP*DLOG(R+ETA))/CDIP
         F3 = RATIO*(YBAR/(CDIP*(R+DBAR)) - DLOG(R+ETA))
     1         + SDIP*F4/CDIP
         F1 = -RATIO*(PSI/(CDIP*(R+DBAR))) - SDIP*F5/CDIP
      ENDIF
         F2 = -RATIO*DLOG(R+ETA) - F3

      U1SS = -(PSI*Q/(R*(R+ETA))
     1         + TERM + F1*SDIP)/TWOPI
      U2SS = -(YBAR*Q/(R*(R+ETA))
     1         + Q*CDIP/(R+ETA)
     2         + F2*SDIP)/TWOPI
      U3SS = -(DBAR*Q/(R*(R+ETA))
     1         + Q*SDIP/(R+ETA)
     2         + F4*SDIP)/TWOPI
      U1DS = -(Q/R - F3*SDIP*CDIP)/TWOPI
      U2DS = -(YBAR*Q/(R*(R+PSI))
     1         + CDIP*TERM - F1*SDIP*CDIP)/TWOPI
      U3DS = -(DBAR*Q/(R*(R+PSI))
     1         + SDIP*TERM - F5*SDIP*CDIP)/TWOPI
      RETURN
      END
*******************************************************************

      SUBROUTINE GRDWEI (YLON, YLAT, JREGN, I, J, WEI)

C
C********1*********2*********3*********4*********5*********6*********7**
C
C NAME:        GRDWEI
C VERSION:     9302.01   (YYMM.DD)
C WRITTEN BY:  MR. C. RANDOLPH PHILIPP
C PURPOSE:     THIS SUBROUTINE RETURNS THE INDICES OF THE LOWER-LEFT
C              HAND CORNER OF THE GRID CELL CONTAINING THE POINT
C              AND COMPUTES NORMALIZED WEIGHTS FOR 
C              BI-LINEAR INTERPOLATION OVER A PLANE
C              
C  INPUT PARAMETERS FROM ARGUMENT LIST:
C  ------------------------------------
C YLON         LONGITUDE OF POINT IN RADIANS, POSITIVE WEST
C YLAT         LATITUDE OF POINT IN RADIANS, POSITIVE NORTH
C JREGN        ID OF GEOGRAPHIC REGION CONTAINING POINT
C
C  OUTPUT PARAMETERS FROM ARGUMENT LIST:
C  -------------------------------------
C I, J         THE COORDINATES OF LOWER LEFT CORNER OF THE GRID
C              CONTAINING THE ABOVE POSITION
C WEI          A TWO BY TWO ARRAY CONTAINING THE NORMALIZED WEIGHTS
C              FOR THE CORNER VECTORS
C
C  GLOBAL VARIABLES AND CONSTANTS:
C  -------------------------------
C NONE
C
C    THIS MODULE CALLED BY:   COMVEL
C
C    THIS MODULE CALLS:       NONE
C
C    INCLUDE FILES USED:      NONE
C
C    COMMON BLOCKS USED:      /CDGRID/, /CONST/
C
C    REFERENCES:  SEE RICHARD SNAY
C
C    COMMENTS:
C
C********1*********2*********3*********4*********5*********6*********7**
C    MOFICATION HISTORY:
C::9302.11, CRP, ORIGINAL CREATION FOR DYNAP
C::9511.09, RAS, MODIFIED FOR HTDP
C::9712.05, RAS, MODIFIED TO ACCOUNT FOR MULTIPLE GRIDS
C********1*********2*********3*********4*********5*********6*********7**
    
C**** COMPUTES THE WEIGHTS FOR AN ELEMENT IN A GRID

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (NUMGRD = 8)
      DIMENSION WEI(2,2)
      COMMON /CDGRID/ GRDLX(NUMGRD), GRDUX(NUMGRD), 
     1          GRDLY(NUMGRD), GRDUY(NUMGRD),
     1          ICNTX(NUMGRD), ICNTY(NUMGRD), NBASE(NUMGRD)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC

C*** Convert input coordinates to degrees
      POSX = (TWOPI - YLON) * 180.D0 / PI
      POSY = YLAT * 180.D0 / PI

C*** Obtain indices for the lower-left corner of the cell
C*** containing the point
      STEPX = (GRDUX(JREGN) - GRDLX(JREGN)) / ICNTX(JREGN)
      STEPY = (GRDUY(JREGN) - GRDLY(JREGN)) / ICNTY(JREGN)
      I = IDINT((POSX - GRDLX(JREGN))/STEPX) + 1
      J = IDINT((POSY - GRDLY(JREGN))/STEPY) + 1
c     write(6,1001) JREGN, I, J
c1001 format(1x, 'jregn = ', I5 /
c    1       1x, ' i = ', I5 /
c    1       1x, ' j = ', I5)

C*** Compute the limits of the grid cell 
      GRLX = GRDLX(JREGN) + (I - 1) * STEPX
      GRUX = GRLX + STEPX                    
      GRLY = GRDLY(JREGN) + (J - 1) * STEPY                
      GRUY = GRLY + STEPY                     

C*** Compute the normalized weights for the point               
      DENOM = (GRUX - GRLX) * (GRUY - GRLY)
      WEI(1,1) = (GRUX - POSX) * (GRUY - POSY) / DENOM
      WEI(2,1) = (POSX - GRLX) * (GRUY - POSY) / DENOM
      WEI(1,2) = (GRUX - POSX) * (POSY - GRLY) / DENOM
      WEI(2,2) = (POSX - GRLX) * (POSY - GRLY) / DENOM

      RETURN
      END

C*********************************************************************
C
      SUBROUTINE GRDVEC (JREGN, I, J, VEL, B)
C
C********1*********2*********3*********4*********5*********6*********7**
C
C NAME:        GRDVEC
C VERSION:     9302.01   (YYMM.DD)
C WRITTEN BY:  MR. C. RANDOLPH PHILIPP
C PURPOSE:     THIS SUBROUTINE RETRIEVES THE APPROXIMATE VALUES OF THE
C              GRID NODE VELOCITIES FOR GRID (I,J) 
C              
C  INPUT PARAMETERS FROM ARGUMENT LIST:
C  ------------------------------------
C JREGN        ID OF GEOGRAPHIC REGION CORRESPONDING TO GRID
C I, J         THE COORDINATES OF LOWER LEFT CORNER OF THE GRID
C              CONTAINING THE ABOVE POSITION
C B            THE ARRAY CONTAINING ALL THE APPROXIMATE VALUES
C              FOR THE ADJUSTMENT
C
C  OUTPUT PARAMETERS FROM ARGUMENT LIST:
C  -------------------------------------
C VEL          A TWO BY TWO ARRAY CONTAINING THE VELOCITY VECTORS
C              FOR THE CORNERS OF THE GRID
C
C  GLOBAL VARIABLES AND CONSTANTS:
C  -------------------------------
C NONE
C
C    THIS MODULE CALLED BY:   COMVEL
C
C    THIS MODULE CALLS:       NONE
C
C    INCLUDE FILES USED:      NONE
C
C    COMMON BLOCKS USED:      NONE     
C
C    REFERENCES:  SEE RICHARD SNAY
C
C    COMMENTS:
C
C********1*********2*********3*********4*********5*********6*********7**
C    MOFICATION HISTORY:
C::9302.11, CRP, ORIGINAL CREATION FOR DYNAP
C::9712.05, RAS, MODIFIED FOR HTDP (version 2.2)
C********1*********2*********3*********4*********5*********6*********7**


      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION VEL(2,2,3), B(*)

      DO 30 II = 0,1
         DO 20 IJ = 0,1
            DO 10 IVEC = 1, 3
               INDEX = IUNGRD(JREGN, I + II, J + IJ, IVEC)
               VEL(II + 1, IJ + 1, IVEC) = B(INDEX)
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE   

      RETURN
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
C*********************************************************
C
      SUBROUTINE TOMNT( IYR, IMON, IDAY, IHR, IMN, MINS )
C
C********1*********2*********3*********4*********5*********6*********7**
C
C NAME:       TOMNT (ORIGINALLY IYMDMJ)
C VERSION:    9004.17
C WRITTEN BY: M. SCHENEWERK
C PURPOSE:    CONVERT DATE TO MODIFIED JULIAN DATE PLUS UT
C
C INPUT PARAMETERS FROM THE ARGUEMENT LIST:
C -----------------------------------------
C IDAY              DAY
C IMON              MONTH
C IYR               YEAR
C
C OUTPUT PARAMETERS FROM ARGUEMENT LIST:
C --------------------------------------
C MINS              MODIFIED JULIAN DATE IN MINUTES
C
C
C LOCAL VARIABLES AND CONSTANTS:
C ------------------------------
C A                 TEMPORARY STORAGE
C B                 TEMPORARY STORAGE
C C                 TEMPORARY STORAGE
C D                 TEMPORARY STORAGE
C IMOP              TEMPORARY STORAGE
C IYRP              TEMPORARY STORAGE
C
C GLOBAL VARIABLES AND CONSTANTS:
C ------------------------------
C
C
C       THIS MODULE CALLED BY: GENERAL USE
C
C       THIS MODULE CALLS:     DINT
C
C       INCLUDE FILES USED:
C
C       COMMON BLOCKS USED:       
C
C       REFERENCES:            DUFFETT-SMITH, PETER  1982, 'PRACTICAL
C                              ASTRONOMY WITH YOUR CALCULATOR', 2ND
C                              EDITION, CAMBRIDGE UNIVERSITY PRESS,
C                              NEW YORK, P.9
C
C       COMMENTS:              THIS SUBROUTINE REQUIRES THE FULL YEAR,
C                              I.E. 1992 RATHER THAN 92.  
C
C********1*********2*********3*********4*********5*********6*********7**
C::LAST MODIFICATION
C::8909.06, MSS, DOC STANDARD IMPLIMENTED
C::9004.17, MSS, CHANGE ORDER YY MM DD
C********1*********2*********3*********4*********5*********6*********7**
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER*4 (I-N)
C
      INTEGER*4     A, B, C, D

      IYRP = IYR
C
C........  0.0  EXPLICIT INITIALIZATION
C
      IF( IMON .LT. 3 ) THEN
        IYRP= IYRP - 1
        IMOP= IMON + 12
      ELSE
        IMOP= IMON
      END IF
C
C........  1.0  CALCULATION
C
      A=  IYRP*0.01D0
      B=  2 - A + DINT( A*0.25D0 )
      C=  365.25D0*IYRP
      D=  30.6001D0*(IMOP + 1)
      MINS =  (B + C + D + IDAY - 679006) * (24 * 60)
     &       + (60 * IHR) + IMN
C      
      RETURN
      END
*****************************************************
      SUBROUTINE TOVNEU(GLAT,GLON,VX,VY,VZ,VN,VE,VU)

*** Convert velocities from vx,vy,vz to vn,ve,vu

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      SLAT = DSIN(GLAT)
      CLAT = DCOS(GLAT)
      SLON = DSIN(GLON)
      CLON = DCOS(GLON)

      VN = -SLAT*CLON*VX - SLAT*SLON*VY + CLAT*VZ
      VE = -SLON*VX + CLON*VY
      VU = CLAT*CLON*VX + CLAT*SLON*VY + SLAT*VZ

      RETURN
      END
***************************************************
      SUBROUTINE TOVXYZ(GLAT,GLON,VN,VE,VU,VX,VY,VZ)
*** Convert velocities from vn,ve,vu to vx,vy,vz
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      SLAT = DSIN(GLAT)
      CLAT = DCOS(GLAT)
      SLON = DSIN(GLON)
      CLON = DCOS(GLON)

      VX = -SLAT*CLON*VN - SLON*VE + CLAT*CLON*VU
      VY = -SLAT*SLON*VN + CLON*VE + CLAT*SLON*VU
      VZ =  CLAT*VN + SLAT*VU

      RETURN
      END
*****************************************************

      subroutine to_std_dev_xyz_velocity(glat,glon,
     &    sn,se,su,sx,sy,sz)

*** Given standard deviations of velocity in north-south-up
*** Compute standard deviations of velocity in x-y-z
*** Assuming the covariances among north-south-up velocity
***     components are zero

      implicit integer*4 (i-n)
      implicit double precision (a-h,o-z)

      slat = dsin(glat)
      clat = dcos(glat)
      slon = dsin(glon)
      clon = dcos(glon)

      sx = dsqrt( (slat*clon*sn)**2
     &           + (slon*se)**2
     &           + (clat*clon*su)**2 )

      sy = dsqrt( (slat*slon*sn)**2
     &           + (clon*se)**2
     &           + (clat*slon*su)**2 )

      sz = dsqrt( (clat*sn)**2 + (slat*su)**2 )

      return
      end
 
**************************************************
      SUBROUTINE DPLACE

*** Predict displacements between two dates.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (numref = 17)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /FILES/ LUIN,LUOUT, I1, I2, I3, I4, I5, I6
      CHARACTER*80 CARD
      CHARACTER    record*120
      CHARACTER*30 NAMEF,NAME,NAMEBB, NAMEIF
      CHARACTER*24 NAME24
      CHARACTER*50 NAME50
      CHARACTER*17 BLAB
      CHARACTER*10 NAMEG
      CHARACTER*4 TYPE
      CHARACTER*1 OPTION,JN,JW, VOPT
      CHARACTER*1 LATDIR, LONDIR
      CHARACTER*1 ANSWER
      character*24 frame1
      LOGICAL TEST

      BLAB = 'OUTSIDE OF REGION'

      WRITE(LUOUT,10)
   10 FORMAT(
     1 ' Displacements will be predicted from time T1 to time T2.'/)
      WRITE(LUOUT,* ) ' Please enter T1 '
   20 CALL GETMDY(MONTH1,IDAY1,IYEAR1,DATE1,MIN1,TEST)

      WRITE(LUOUT,22)
   22 format (/ ' Please enter T2 ')
   35 CALL GETMDY(MONTH2,IDAY2,IYEAR2,DATE2,MIN2,TEST)

      WRITE(LUOUT,45)
   45 FORMAT( /
     1    ' Please enter the name for the file to contain'/
     2    ' the predicted displacements.  ')
      READ(LUIN,50,err=601,iostat=ios) NAMEF
      if (ios /= 0) goto 601
   50 FORMAT(A30)
      OPEN(I2,FILE=NAMEF,STATUS='UNKNOWN')
      CALL HEADER

*** Choosing reference frame for displacements
   56 WRITE(LUOUT,55)
   55 FORMAT(' ************************************************'/
     1   ' Select the reference frame to be used for specifying'/
     2   ' positions and displacements.  '/)
      call MENU1(iopt, frame1)

      IF(IOPT .GE. 1 .AND . IOPT .LE. numref) THEN
          WRITE(I2,57) frame1
   57     FORMAT(' DISPLACEMENTS IN METERS RELATIVE TO ', a24)
      ELSE
          WRITE(LUOUT,70)
   70     FORMAT(' Improper selection--try again.  ')
          GO TO 56
      ENDIF
            
      WRITE(I2,71) MONTH1,IDAY1,IYEAR1,MONTH2,IDAY2,IYEAR2,
     1             DATE1, DATE2
   71 FORMAT (
     1  ' FROM ',I2.2,'-',I2.2,'-',I4,' TO ',
     2           I2.2,'-',I2.2,'-',I4,' (month-day-year)'/
     2  ' FROM ',F8.3, ' TO ',F8.3, ' (decimal years)'//
     3  'NAME OF SITE             LATITUDE          LONGITUDE    ',          
     4        '        NORTH    EAST    UP ')
   75 WRITE(LUOUT,80)
   80 FORMAT(' ********************************'/
     1   ' Displacements will be predicted at each point whose',/
     2   ' horizontal position is specified.',/
     4   ' Please indicate how you wish to supply positions.'/ )
   85 WRITE(LUOUT,86)
   86 FORMAT(
     5   '    0. No more points. Return to main menu.'/
     6   '    1. Individual points entered interactively.'/
     7   '    2. Points on a specified grid.'/
     8   '    3. The *80* records in a specified blue-book file.'/
     9   '    4. Points on a specified line.  ' /
     1   '    5. Batch file of delimited records of form: ' /
     2   '       LAT,LON,TEXT ' /
     4   '       LAT = latitude in degrees (positive north/DBL PREC)' /
     5   '       LON = longitude in degrees (positive west/DBL PREC)' /
     7   '       TEXT = Descriptive text (CHARACTER*24) ' /
     8   '       Example:  ' /
     9   '       40.731671553,112.212671753,SALT AIR '/)
      READ(LUIN,'(A1)',err=602,iostat=ios) OPTION
      if (ios /= 0) goto 602

      IF(OPTION .EQ. '0') THEN
	  GO TO 510               
      ELSEIF(OPTION .EQ. '1') THEN
	  CALL GETPNT(LATD,LATM,SLAT,LATDIR,LOND,LONM,SLON,
     1        LONDIR,NAME24, X,Y,Z,YLAT,YLON,EHT)
	  ELON = - YLON
	  call GETVLY( YLAT, ELON, VX, VY, VZ, VN, VE, VU, VOPT, 210)
	  if (vopt .eq. '0') then
	     call PREDV( ylat, ylon, eht, date1, iopt,
     1                   jregn, vn, ve, vu)
	     if (jregn .eq. 0) then
		write(luout, 140)
		go to 85 
             endif
          endif
          call NEWCOR(ylat, ylon, eht, min1, min2, 
     1       ylatt, ylont, ehtnew, dn, de, du, vn, ve, vu)
  140          FORMAT(' ****************************************'/
     1                ' A displacement can not be predicted because'/
     1                ' the point is outside of the modeled region.'/
     2         ' For additional displacements, please indicate how'/
     3         ' you wish to supply the horizontal coordinates.'/)
	   WRITE(LUOUT,150) DN, DE,DU
  150          FORMAT(' *************************************'/
     1           ' Northward displacement = ',F7.3,' meters.'/
     1           ' Eastward displacement  = ',F7.3,' meters.'/
     1           ' Upward displacement    = ',F7.3,' meters.'/
     1           ' ****************************************'//
     2           ' For additional displacements, please indicate how'/
     3           ' you wish to supply the horizontal coordinates.'/)
	   WRITE(I2,160) NAME24,LATD,LATM,SLAT,LATDIR,
     1	          LOND,LONM,SLON,LONDIR,DN,DE,DU
  160          FORMAT(A24,1X,I2,1X,I2,1X,F8.5,1X,A1,2X,
     1            I3,1X,I2,1X,F8.5,1X,A1,1X,3F8.3)
      ELSEIF(OPTION .EQ. '2') THEN

	  CALL GETGRD(NAMEG,MINLAT,MAXLAT,IDS,MINLON,MAXLON,JDS)
	  I = -1
  280     I = I + 1
	  LAT = MINLAT + I*IDS
	  IF(LAT .GT. MAXLAT) GO TO 296
	  XLAT = DBLE(LAT)/RHOSEC
          CALL TODMSS(XLAT,LATD,LATM,SLAT,ISIGN)
	  LATDIR = 'N'
	  IF (ISIGN .eq. -1) LATDIR = 'S'
	  J = -1
  290     J = J + 1
          YLAT = XLAT
	  LON = MINLON + J*JDS
	  IF(LON .GT. MAXLON) GO TO 280
	  YLON = DBLE(LON)/RHOSEC
          CALL TODMSS(YLON,LOND,LONM,SLON,ISIGN)
	  LONDIR = 'W'
	  IF (ISIGN .eq. -1) LONDIR = 'E'
          EHT = 0.D0
	  call PREDV( ylat,ylon, eht, date1, iopt,
     1       jregn, vn, ve, vu)
          IF(JREGN .eq. 0) THEN     
                WRITE(I2,291)NAMEG,I,J,LATD,LATM,SLAT,LATDIR,LOND,
     1            LOND,LONM,SLON,LONDIR,BLAB
  291           FORMAT(A10,2I4, 7X,I2,1X,I2,1X,F8.5,1X,A1,2X,I3,
     1                 1X,I2,1X,F8.5,1X,A1,1X,A17)
          ELSE
	     call NEWCOR( ylat, ylon, eht, min1, min2,
     1          ylatt, ylont, ehtnew, dn, de, du, vn, ve, vu)
	     WRITE(I2,295)NAMEG,I,J,LATD,LATM,SLAT,LATDIR,
     1                 LOND,LONM,SLON,LONDIR,DN,DE,DU
  295        FORMAT(A10,2I4, 7X,I2,1X,I2,1X,F8.5,1X,A1,2X,I3,1X,I2,
     1                 1X,F8.5,1X,A1,1X,3F8.3)
          ENDIF
	  GO TO 290
  296     WRITE(LUOUT,297)
  297     FORMAT(' ***********************************'/
     1    ' Displacements have been calculated for the specified'/
     2    ' grid.  If you wish to calculate additional displacements,'/
     3    ' please indicate how you will supply the coordinates.'/)
      ELSEIF(OPTION .EQ. '3') THEN
	  VOPT = '0'
	  WRITE(LUOUT,300)
  300     FORMAT(' Enter name of blue-book file  ')
	  READ(LUIN,310,err=603,iostat=ios) NAMEBB
          if (ios /= 0) goto 603
  310     FORMAT(A30)
	  OPEN(I1,FILE=NAMEBB,STATUS='OLD')
  320     READ(I1,330,END=350,err=604,iostat=ios) CARD
          if (ios /= 0) goto 604
  330     FORMAT(A80)
	  TYPE = CARD(7:10)
	  IF(TYPE .EQ. '*80*') THEN
	      READ(CARD,340,err=605,iostat=ios) NAME,LATD,LATM,SLAT,JN,
     1                           LOND,LONM,SLON,JW
              if (ios /= 0) goto 605
  340         FORMAT(BZ,14X,A30,I2,I2,F7.5,A1,I3,I2,F7.5,A1)
              NAME24 = NAME(1:24)
C             IF(JN.EQ.'S' .OR. JW.EQ.'E')GO TO 320
	      YLAT =(DBLE((LATD*60+LATM)*60)+SLAT)/RHOSEC
	      IF (JN .eq. 'S') YLAT = -YLAT
	      YLON =(DBLE((LOND*60+LONM)*60)+SLON)/RHOSEC
	      IF (JW .eq. 'E') YLON = -YLON
              EHT = 0.D0
	      call PREDV( ylat, ylon, eht, date1, iopt,
     1          jregn, vn, ve, vu)
              IF(JREGN .EQ. 0) THEN      
                 WRITE(I2,345)NAME24,LATD,LATM,SLAT,JN,LOND,LONM,SLON,
     1                          JW,BLAB
  345            FORMAT(A24,1X,I2,1X,I2,1X,F8.5,1X,A1,2X,I3,1X,I2,1X,
     1                    F8.5,1X,A1,1X,A17)
              ELSE
		 call NEWCOR( ylat, ylon, eht, min1, min2,
     1             ylatt, ylont, ehtnew, dn, de, du, vn, ve, vu)
	         WRITE(I2,160)NAME24,LATD,LATM,SLAT,JN,LOND,LONM,SLON,
     1             JW,DN,DE,DU
              ENDIF
	  ENDIF
	  GO TO 320
  350     CLOSE(I1,STATUS='KEEP')
	  WRITE(LUOUT,360)
  360     FORMAT(' ************************************'/
     1    ' Displacements have been calculated for the specified'/
     2    ' blue-book file.  If you wish to calculate additional'/
     3    ' displacements, please indicate how you will supply'/
     4    ' the horizontal coordinates.'/)
      ELSEIF(OPTION .EQ. '4') THEN
	  VOPT = '0'
	  CALL GETLYN(NAMEG,XLAT,XLON,FAZ,BAZ,XMIN,XMAX,XINC)
	  XLON = TWOPI - XLON
	  I = -1
  400     I = I + 1
	  S = XMIN + I*XINC
	  IF(S .GT. XMAX) GO TO 430
	  IF(S .LT. 0.0D0) THEN
	      S1 = -S
	      AZ = BAZ
          ELSE
	      S1 = S
	      AZ = FAZ
          ENDIF
	  CALL DIRCT1(XLAT,XLON,YLAT,YLON,AZ,AZ1,S1)
	  YLON = TWOPI - YLON
          CALL TODMSS(YLAT,LATD,LATM,SLAT,ISIGN)
	  LATDIR = 'N'
	  IF (ISIGN .eq. -1) LATDIR = 'S'
          CALL TODMSS(YLON,LOND,LONM,SLON,ISIGN)
	  LONDIR = 'W'
	  IF (ISIGN .eq. -1) LONDIR = 'E'
          EHT = 0.D0
	  call PREDV( ylat, ylon, eht, date1, iopt,
     1       jregn, vn, ve, vu)
          IF(JREGN .eq. 0) THEN
              WRITE(I2,405)NAMEG,I,LATD,LATM,SLAT,LATDIR,
     1               LOND,LONM,SLON,LONDIR, BLAB
  405         FORMAT(A10,I4,11X,I2,1X,I2,1X,F8.5,1X,A1,2X,I3,1X,I2,1X,
     1               F8.5,1X,A1,1X,A17)  
          ELSE
	      call NEWCOR( ylat, ylon, eht, min1, min2,
     1           ylatt, ylont, ehtnew, dn, de, du, vn, ve, vu)
	      WRITE(I2,410)NAMEG,I,LATD,LATM,SLAT,LATDIR,LOND,LONM,
     1                 SLON,LONDIR,DN,DE,DU
  410         FORMAT(A10,I4,11X,I2,1X,I2,1X,F8.5,1X,A1,2X,I3,1X,I2,
     1           1X,F8.5,1X,A1,1X,3F8.3)
          ENDIF
	  GO TO 400
  430     WRITE(LUOUT,440)
  440     FORMAT(' ***************************************'/
     1    ' Displacements have been calculated for the specified'/
     2    ' line.  If you wish to calculate additional displacements,'/
     3    ' please indicate how you will supply the coordinates.'/)
      ELSEIF(OPTION .EQ. '5') THEN
          VOPT = '0'
          EHT = 0.0D0
          write(luout,450)
  450     format(' Enter name of batch file ')
          read(luin, 451,err=606,iostat=ios) NAMEIF
          if (ios /= 0) goto 606
  451     format(a30)
          open(I1,FILE=NAMEIF,STATUS='OLD')
  455     read(I1,'(a)',END=460,err=607,iostat=ios) record
          if (ios /= 0) goto 607
          call interprate_latlon_record (record,XLAT,XLON,name50)
          name24 = name(1:24)
          YLAT = (XLAT*3600.D0)/RHOSEC
          YLON = (XLON*3600.D0)/RHOSEC
          CALL TODMSS(YLAT,LATD,LATM,SLAT,ISIGN)
          IF (ISIGN .EQ. 1) THEN
               JN = 'N'
          ELSE
               JN = 'S'
          ENDIF
          CALL TODMSS(YLON, LOND,LONM,SLON,ISIGN)
          IF (ISIGN .EQ. 1) THEN
               JW = 'W'
          ELSE
               JW = 'E'
          ENDIF
          CALL PREDV(YLAT,YLON,EHT,DATE1,IOPT,
     1               JREGN,VN,VE,VU)
          IF (JREGN .EQ. 0) THEN
             Write(I2,345) NAME24,LATD,LATM,SLAT, JN,
     1                        LOND,LONM,SLON,JW, BLAB
          ELSE
             CALL NEWCOR (YLAT, YLON, EHT, MIN1, MIN2,
     1             YLATT, YLONT, EHTNEW, DN,DE,DU,VN,VE,VU)
             WRITE(I2, 160) NAME24, LATD, LATM, SLAT, JN,
     1               LOND, LONM, SLON, JW, DN, DE , DU
          ENDIF
          GO TO 455
  460     CLOSE(I1, STATUS = 'KEEP')
          write(LUOUT, 470)
  470     format(' ******************************************'/
     1    ' Displacements have been calculated for the specified'/
     1    ' file.  If you wish to calculate additional displacements,'/
     1    ' please indicate how you will supply the horizontal '/
     1    ' coordinates.'/)
      ELSE
	  WRITE(LUOUT,500)
  500     FORMAT(' Improper entry--select again.  ')
      ENDIF
      GO TO 85
  510 CONTINUE
      CLOSE(I2, STATUS = 'KEEP')
      RETURN

  600 write (*,'(/)') 
      write (*,*) 'Wrong answer in DPLACE: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  601 write (*,'(/)') 
      write (*,*) 'Wrong file name in DPLACE: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  602 write (*,'(/)') 
      write (*,*) 'Wrong OPTION in DPLACE: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  603 write (*,'(/)') 
      write (*,*) 'Wrong bbname in DPLACE: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  604 write (*,'(/)') 
      write (*,*) 'Wrong CARD from bbfile in DPLACE:ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  605 write (*,'(/)') 
      write (*,*) 'Wrong CARD80 from bbfile in DPLACE:ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  606 write (*,'(/)') 
      write (*,*) 'Wrong batch file name in DPLACE:ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  607 write (*,'(/)') 
      write (*,*) 'Wrong record from batch in DPLACE:ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

      END
********************************************************************
      SUBROUTINE VELOC

*** Compute velocities at specified locations

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (numref = 17)
      parameter (nrsrch = 0)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /FILES/ LUIN,LUOUT, I1, I2, I3, I4, I5, I6
      CHARACTER*80 CARD
      CHARACTER*30 NAMEF,NAME,NAMEBB,PVFILE
      CHARACTER*24 NAME24
      CHARACTER*4  NAME4
      CHARACTER*17 BLAB
      CHARACTER*10 NAMEG
      CHARACTER*4 TYPE
      CHARACTER*1 OPTION,JN,JW,PVOUT
      CHARACTER*1 LATDIR, LONDIR
      character*24 frame1
      character*120 record
      character*50 name50

      BLAB = 'OUTSIDE OF REGION'

      WRITE(LUOUT,10)
   10 FORMAT(' Please enter name for the file to contain the'/
     1       ' predicted velocities.  ')
      READ(LUIN,20,err=700,iostat=ios) NAMEF
      if (ios /= 0) goto 700
   20 FORMAT(A30)
      OPEN(I2,FILE=NAMEF,STATUS='UNKNOWN')
      CALL HEADER

*** Choosing reference system for velocities
 1000 WRITE(LUOUT,1001)
 1001 FORMAT(' **************************************************'/
     1   ' Select the reference frame to be used for specifying'/
     2   ' positions and velocities. '/ )
    
	   call MENU1(iopt, frame1)
      if (iopt .ge. 1 .and. iopt .le. numref) then
         WRITE(I2,1002) frame1
 1002    FORMAT('Velocities (with standard deviations) in mm/yr',/
     &       '  relative to ',a24 /)

      ELSE
           write(luout, 1080)
 1080      format(' Improper selection -- try again.  ')
           GO TO 1000
      ENDIF

*** Choosing input format for locations where velocities are to be predicted
      WRITE(LUOUT,30)
   30 FORMAT(' ************************************************'/
     1   ' Velocities will be predicted at each point whose'/
     2   ' horizontal position is specified.  Please indicate'/
     4   ' how you wish to supply positions.'/)
   40 WRITE(LUOUT,50)
   50 FORMAT(
     1 '    0... No more points. Return to main menu.'/
     2 '    1... Individual points entered interactively.'/
     3 '    2... Points on a specified grid.'/
     4 '    3... The *80* records in a specified blue-book file.'/
     5 '    4... Points on a specified line.  '/
     6 '    5... Batch file of delimited records of form: '/
     7 '         LAT,LON,TEXT '/
     8 '         LAT = latitude in degrees (positive north/DBL PREC) '/
     9 '         LON = longitude in degrees (positive west/DBL PREC) '/
     1 '         TEXT = Descriptive text (CHARACTER*24) '/
     1 '         Example:  '/
     2 '         40.731671553,112.212671753,SALT AIR '/ )
 
      READ(LUIN,60,err=704,iostat=ios) OPTION
      if (ios /= 0) goto 704
   60 FORMAT(A1)

      IF(OPTION .EQ. '0') THEN
	   GO TO 610               
      ELSEIF(OPTION .EQ. '1') THEN
	       CALL GETPNT(LATD,LATM,SLAT,LATDIR,LOND,LONM,SLON,
     1           LONDIR, NAME24, X,Y,Z,YLAT,YLON,EHT)
           CALL GTOVEL(ylat,ylon,eht,
     1        VN,VE,VU,VX,VY,VZ,JREGN,IOPT,SN,SE,SU,sx,sy,sz)
           IF(JREGN .eq. 0 ) THEN
               WRITE(LUOUT,80)
   80          FORMAT(' *************************************'/
     1           ' A velocity can not be estimated because'/
     1           ' the point is outside of the modeled region.'/
     2           ' For additional velocities, please indicate how'/
     3           ' you wish to supply the horizontal coordinates.'/)
           ELSE
               
	       WRITE(LUOUT,90) VN,SN,VE,SE,VU,SU,
     &                        VX,SX,VY,SY,VZ,SZ
   90          FORMAT(' **************************************'/
     1            ' Northward velocity = ',F6.2,' +/- ',f4.2,' mm/yr'/
     1            ' Eastward velocity  = ',F6.2,' +/- ',f4.2,' mm/yr'/
     1            ' Upward velocity    = ',F6.2,' +/- ',f4.2,' mm/yr'//
     1            ' X-dim. velocity    = ',F6.2,' +/- ',f4.2,' mm/yr'/
     1            ' Y-dim. velocity    = ',F6.2,' +/- ',f4.2,' mm/yr'/
     1            ' Z-dim. velocity    = ',F6.2,' +/- ',f4.2,' mm/yr'/
     1            ' **************************************'//
     2         ' For additional velocities, please indicate how'/
     3         ' you wish to specify the horizontal coordinates.'/)
	       WRITE(I2,1060)NAME24,LATD,LATM,SLAT,LATDIR,VN,SN,LOND,LONM, 
     1                SLON, LONDIR,VE,SE,EHT,VU,SU
               write(i2,1061) X,VX,SX,Y,VY,SY,Z,VZ,SZ
 1060          FORMAT(/10X,A24,/
     1'LATITUDE   = ',2I3,F9.5,1X,A1,'  NORTH VELOCITY =',F7.2,' +/- ',
     & f4.2,/
     2'LONGITUDE  = ',2I3,F9.5,1X,A1,'  EAST VELOCITY  =',F7.2,' +/- ',
     & f4.2,/ 
     3'ELLIPS. HT. = ',F10.3,' m', 6X,'UP VELOCITY    =',F7.2,' +/- ',
     & f4.2)
 1061  Format(
     & 'X =',F13.3,' m',14X,'X VELOCITY     =',F7.2,' +/- ',
     &    f4.2,/
     5 'Y =',F13.3,' m',14X,'Y VELOCITY     =',F7.2,' +/- ',
     &    f4.2,/
     6 'Z =',F13.3,' m',14X,'Z VELOCITY     =',F7.2,' +/- ',
     &    f4.2)
  100          FORMAT(A24,1X,I2,1X,I2,1X,F8.5,1X,A1,2X,I3,1X,I2,1X,F8.5,
     1            1X,A1,1X,3(F8.2,' +/- ',F4.2))
               IF(PVOUT .EQ. 'Y') THEN
                   CALL PVPRNT(LATD,LATM,SLAT,LOND,LONM,SLON,VN,VE,VU)
               ENDIF
           ENDIF
      ELSEIF(OPTION .EQ. '2') THEN
           EHT = 0.D0
           WRITE(I2,105)
  105      FORMAT(
     1     'NAME OF SITE',3X,'LATITUDE          LONGITUDE    ',
     2     '           NORTH           EAST             UP'/)
  106      FORMAT(
     1     'NAME OF SITE', 14X,'LATITUDE          LONGITUDE    ',
     2     '           NORTH            EAST             UP'/)
  107      FORMAT(
     1     'NAME OF LINE', 6X,'LATITUDE          LONGITUDE    ',
     2     '           NORTH            EAST             UP'/)
	   CALL GETGRD(NAMEG,MINLAT,MAXLAT,IDS,MINLON,MAXLON,JDS)
	   I = -1
  110      I = I + 1
	   LAT = MINLAT + I*IDS
	   IF(LAT .GT. MAXLAT) GO TO 150
	   XLAT = DBLE(LAT)/RHOSEC
           CALL TODMSS(XLAT,LATD,LATM,SLAT,ISIGN)
	   LATDIR = 'N'
	   IF (ISIGN .eq. -1) LATDIR = 'S'
	   J = -1
  120  J = J + 1
       YLAT = XLAT
	   LON = MINLON + J*JDS
	   IF(LON .GT. MAXLON) GO TO 110
	   YLON = DBLE(LON)/RHOSEC
       CALL TODMSS(YLON,LOND,LONM,SLON,ISIGN)
	   LONDIR = 'W'
	   IF (ISIGN .eq. -1) LONDIR = 'E'

           CALL GTOVEL(ylat,ylon,eht,
     1           VN,VE,VU,VX,VY,VZ,JREGN,IOPT,SN,SE,SU,sx,sy,sz)
           IF(JREGN .EQ. 0 ) THEN
               WRITE(I2,125)NAMEG,LATD,LATM,SLAT,LATDIR,LOND,LONM,
     1                      SLON,LONDIR,BLAB
  125          FORMAT(A10,4X,I2,1x,i2,1x,F8.5,1X,A1,2X,I3,1X,I2,
     1                1X,F8.5,1X,A1, 1X,A17)
           ELSE
       	       WRITE(I2,130)NAMEG,LATD,LATM,SLAT,LATDIR,
     1                  LOND,LONM,SLON,LONDIR,VN,SN,VE,SE,VU,SU
  130          FORMAT(A10,4X,I2,1X,I2,1X,F8.5,1X,A1,2X,
     1             I3,1X,I2,1X,F8.5,1X,A1,1X,3(F7.2,' +/- ',F4.2))
               IF(PVOUT .EQ. 'Y') THEN
                   CALL PVPRNT(LATD,LATM,SLAT,LOND,LONM,SLON,VN,VE,VU)
               ENDIF
           ENDIF
	  GO TO 120
  150      WRITE(LUOUT,160)
  160      FORMAT(' **********************************************'/
     1     ' Velocities have been calculated for the specified grid.'/
     2     ' If you wish to calculate additional velocities, please'/
     3     ' indicate how you will specify positional coordinates.'/)
      ELSEIF(OPTION .EQ. '3') THEN
           EHT = 0.D0
           WRITE(I2,106)
	   WRITE(LUOUT,200)
  200      FORMAT(' Enter name of blue-book file.  ')
	   READ(LUIN,210,err=705,iostat=ios) NAMEBB
           if (ios /= 0) goto 705
  210      FORMAT(A30)
	   OPEN(I1,FILE=NAMEBB,STATUS='OLD')
  220      READ(I1,230,END=250,err=706,iostat=ios) CARD
           if (ios /= 0) goto 706
  230      FORMAT(A80)
	   TYPE = CARD(7:10)
      IF(TYPE .EQ. '*80*') THEN
		READ(CARD,240,err=707,iostat=ios)NAME,LATD,LATM,SLAT,JN,
     1                             LOND,LONM,SLON,JW
                if (ios /= 0) goto 707
  240           FORMAT(BZ,14X,A30,I2,I2,F7.5,A1,I3,I2,F7.5,A1)
                NAME24 = NAME(1:24)
		YLAT =(DBLE((LATD*60+LATM)*60)+SLAT)/RHOSEC
		YLON =(DBLE((LOND*60+LONM)*60)+SLON)/RHOSEC
		if (jn .eq. 'S') ylat = -ylat
		if (jw .eq. 'E') ylon = -ylon
        CALL GTOVEL(ylat,ylon,eht, 
     1             VN,VE,VU,VX,VY,VZ,JREGN,IOPT,SN,SE,SU,sx,sy,sz)
        IF(JREGN .EQ. 0 ) THEN
              WRITE(I2,245)NAME24,LATD,LATM,SLAT,JN,LOND,LONM,
     1                          SLON,JW,BLAB
  245         FORMAT(A24,1X,I2,1X,I2,1X,F8.5,1X,A1,2X,I3,1X,I2,
     1                    1X,F8.5,1X,A1,1X,A17)
         ELSE
               WRITE(I2,100)NAME24,LATD,LATM,SLAT,JN,LOND,LONM,
     1               SLON,JW,VN,SN,VE,SE,VU,SU  
               IF(PVOUT .EQ. 'Y') THEN
                    CALL PVPRNT(LATD,LATM,SLAT,LOND,LONM,SLON,VN,VE,VU)
               ENDIF
          ENDIF
      ENDIF
      GO TO 220
  250 CLOSE(I1,STATUS='KEEP')
      WRITE(LUOUT,260)
  260      FORMAT(' ****************************************'/
     1     ' Velocities have been calculated for the specified '/
     2     ' blue-book file.  If you wish to calculate additional'/
     3     ' velocities, please indicate how you will supply the'/
     4     ' horizontal coordinates.'/)
      ELSEIF(OPTION .EQ. '4') THEN
           EHT = 0.D0
           WRITE(I2,107)
	   CALL GETLYN(NAMEG,XLAT,XLON,FAZ,BAZ,XMIN,XMAX,XINC)
	   XLON = TWOPI - XLON
	   I = -1
  300      I= I+1
	   S = XMIN + I*XINC
	   IF(S .GT. XMAX) GO TO 330
	   IF(S .LT. 0.0D0) THEN
	       S1 = -S
	       AZ = BAZ
           ELSE
	       S1 = S
	       AZ = FAZ
           ENDIF
	   CALL DIRCT1(XLAT,XLON,YLAT,YLON,AZ,AZ1,S1)
	   YLON = TWOPI - YLON
           CALL TODMSS(YLAT,LATD,LATM,SLAT,ISIGN)
	   LATDIR = 'N'
	   IF (ISIGN .eq. -1) LATDIR = 'S'
           CALL TODMSS(YLON,LOND,LONM,SLON,ISIGN)
	   LONDIR = 'W'
	   IF (ISIGN .eq. -1) LONDIR = 'E'
           CALL GTOVEL(ylat,ylon,eht,
     1           VN,VE,VU,VX,VY,VZ,JREGN,IOPT,SN,SE,SU,sx,sy,sz)
           IF(JREGN .EQ. 0 ) THEN
               WRITE(I2,305)NAMEG,I,LATD,LATM,SLAT,LATDIR,LOND,LONM,
     1                      SLON,LONDIR,BLAB
  305          FORMAT(A10,I4,3X,I2,1X,I2,1X,F8.5,1X,A1,2X,I3,1X,
     1                I2,1X,F8.5,1X,A1,1X,A17)
           ELSE
	       WRITE(I2,310)NAMEG,I,LATD,LATM,SLAT,LATDIR,LOND,LONM,SLON,
     1                  LONDIR,VN,SN,VE,SE,VU,SU
  310          FORMAT(A10,I4,3X,I2,1X,I2,1X,F8.5,1X,A1,2X,I3,1X,I2,1X,
     1                F8.5,1X,A1,1X,3(F8.2,' +/- ',F4.2))
               IF(PVOUT .EQ. 'Y') THEN
                   CALL PVPRNT(LATD,LATM,SLAT,LOND,LONM,SLON,VN,VE,VU)
               ENDIF
           ENDIF
	   GO TO 300
  330      WRITE(LUOUT,340)
  340      FORMAT(' ***********************************************'/
     1     ' Velocities have been calculated for the specified line.'/
     2     ' If you wish to calculate additional velocities, please'/
     3     ' indicate how you will specify positional coordinates.'/)
      ELSEIF(OPTION .EQ. '5') THEN
           EHT = 0.0D0
           WRITE(I2, 106)
           WRITE(LUOUT,400)
  400      FORMAT(' Enter name of input file.  ')
           READ(LUIN,410,err=708,iostat=ios) NAMEBB
           if (ios /= 0) goto 708
  410      FORMAT(A30)
           OPEN(I1,FILE=NAMEBB,STATUS='OLD')
  420      read(i1,'(a)',end=450,err=709,iostat=ios) record           
           if (ios /= 0) goto 709
C          write(i2,'(a)') record
           call interprate_latlon_record(record,xlat,xlon,name50)
C          write(i2,'(a)') name50
           YLAT = (XLAT*3600.D0)/RHOSEC
           YLON = XLON*3600.d0/RHOSEC
           CALL TODMSS(YLAT, LATD, LATM, SLAT, ISIGN)
           IF (ISIGN .eq. 1) then
              JN = 'N'
           Else
              JN = 'S'
           endif
           call TODMSS(YLON, LOND, LONM, SLON, ISIGN)
           IF (ISIGN .eq. 1) then
              JW = 'W'
           else
              JW = 'E'
           endif
           CALL GTOVEL(ylat,ylon,eht, 
     1         VN,VE,VU,VX,VY,VZ,JREGN,IOPT,SN,SE,SU,sx,sy,sz)
           IF (JREGN .eq. 0) then
              write(I2,245)NAME24,LATD,LATM,SLAT,JN,
     1                     LOND,LONM,SLON,JW,BLAB
           ELSE
*** original code
              write(I2,100)NAME24,LATD,LATM,SLAT,JN,
     1                      LOND,LONM,SLON,JW,VN,SN,VE,SE,VU,SU
              IF(PVOUT .EQ. 'Y') then
                call PVPRNT(LATD,LATM,SLAT,LOND,LONM,SLON,VN,VE,VU)
              ENDIF
           ENDIF
           GO TO 420
  450      CLOSE(I1, STATUS = 'KEEP')
           write(LUOUT, 460)
  460      format(' ********************************'/
     1     ' Velocities have been calculated for the specified'/
     2     ' file.  if you wish to calculate additional velocities, '/
     3     ' please indicate how you will supply the coordinates.'/)
     

      ELSE
	   WRITE(LUOUT,600)
  600      FORMAT(' Improper entry--select again.  ')
      ENDIF
      GO TO 40
  610 CONTINUE
      CLOSE(I2, STATUS = 'KEEP')
      IF( PVOUT .EQ. 'Y') CLOSE(I3, STATUS = 'KEEP')
      RETURN

  700 write (*,'(/)') 
      write (*,*) 'Failed to read file name: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  701 write (*,'(/)') 
      write (*,*) 'Failed to read answer: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  702 write (*,'(/)') 
      write (*,*) 'Failed to read file name: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  703 write (*,'(/)') 
      write (*,*) 'Failed to read velocity: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  704 write (*,'(/)') 
      write (*,*) 'Failed to read OPTION: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  705 write (*,'(/)') 
      write (*,*) 'Failed to read name of bbfile: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  706 write (*,'(/)') 
      write (*,*) 'Failed to read card from bbfile: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  707 write (*,'(/)') 
      write (*,*) 'Failed to read card80 from bbfile: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  708 write (*,'(/)') 
      write (*,*) 'Failed to read file name: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  709 write (*,'(/)') 
      write (*,*) 'Failed to read bbfile: ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop


      END
**************************************************************************

      SUBROUTINE GTOVEL(ylat,ylon,eht,
     1   VN,VE,VU,VX,VY,VZ,JREGN,IOPT,SN,SE,SU,sx,sy,sz)

*** Compute velocity in appropriate reference frame for point with
*** latitude YLAT (radians), longitude YLON (radians, positive west)
*** and ellipsoid height EHT (meters) in this reference frame.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

*** Get reference latitude RLAT and reference longitude RLON in ITRF2014
      IF(IOPT .EQ. 16) THEN
          RLAT = YLAT
          RLON = YLON
      ELSE
          ELON = -YLON
          CALL TOXYZ(YLAT,ELON,EHT,X,Y,Z)
          DATE = 2010.0d0
          CALL XTOITRF2014(X,Y,Z,RLAT,RLON,EHT14,DATE,IOPT)
      ENDIF

*** Get velocity in ITRF2014 
      CALL GETREG(RLAT,RLON,JREGN)
      IF (JREGN .EQ. 0) RETURN
      CALL COMVEL(RLAT,RLON,JREGN,VN,VE,VU,SN,SE,SU)

      ELON = -YLON
      CALL TOVXYZ(YLAT,ELON,VN,VE,VU,VX,VY,VZ)
      call to_std_dev_xyz_velocity(ylat,elon,sn,se,su,sx,sy,sz)

*** Convert velocity into another reference frame if needed
      IF(IOPT .NE. 16) THEN
         CALL VTRANF(X,Y,Z,VX,VY,VZ, 16, IOPT)
         CALL TOVNEU(YLAT, ELON, VX, VY, VZ, VN, VE, VU)
      ENDIF
     
      RETURN
      END
****************************************************************************
      SUBROUTINE XTOITRF2014 (X,Y,Z,RLAT,WLON,EHT14,DATE,IOPT)

*** Converts X,Y,Z in specified datum to latitude and
*** longitude (in radians) and height (meters) in ITRF2014
*** datum with longitude positive west.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      LOGICAL FRMXYZ
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC

*** Convert to cartesian coordinates in ITRF2014   
      if (iopt .eq. 16) then
         x1 = x
         y1 = y
         z1 = z
      else   
         call to_itrf2014(x,y,z,x1,y1,z1,date,iopt)
      endif   
         
***Convert to geodetic coordinates
      IF(.NOT.FRMXYZ(X1,Y1,Z1,RLAT,ELON,EHT14))STOP 666

      WLON = -ELON
 100  IF(WLON .LT. 0.D0) THEN
          WLON = WLON + TWOPI
          GO TO 100
      ENDIF
      RETURN
      END
*********************************************************

      subroutine TRFPOS1

*** Transform position across time
*** and between reference frames

*** subroutine TRFPOS1 replaces subroutine TRFPOS

*** This subroutine uses a two-step process:
*** Step 1 transforms the positional coordinates provided
***        in the initial reference frame at the initial time
***        to their corresponding coordinates in this same frame
***        at the ending time
*** Step 2 transforms the positional coordinates in the initial
***        frame at the ending time to the corresponding
***        coordinates in the new reference frame at the
***        ending time
*** Thus, step 1 involves only time transfer and step 2 
*** involves only frame transfer      

      implicit double precision (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (numref = 17)
      
      character*80 card,namebb,nameif,name24
      character    record*120                    
      character*30 namef
      character*24 frame1, frame2
      character*1  jn,jw,LATDIR,LONDIR
      character*1  option, answer, vopt
      character*10 Trans4D_version   
      LOGICAL      FRMXYZ
      LOGICAL      TEST
      COMMON /CONST/ A, F, E2, EPS, AF, PI, TWOPI, RHOSEC
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
      COMMON /VERSION/ Trans4D_version

      write(luout,80)
   80 format(
     1  ' Please enter the name for the file to contain'/
     2  ' the transformed positions.  ')
      read(luin,'(a30)',err=600,iostat=ios) namef
      if (ios /= 0) goto 600
      open(i2, file = namef, status = 'unknown')
      CALL HEADER

      call get_frames(iopt1, iopt2, frame1, frame2)

      write(luout,120)
  120 format (/
     1   ' Enter the reference date of the input positions.')
  121 CALL GETMDY(MONTH1, IDAY1, IYEAR1, DATE1, MIN1, TEST)
      write(luout,130)
  130 format(/
     1  '  Enter the reference date for the output positions.  ')
  131 CALL GETMDY(MONTH2, IDAY2, IYEAR2, DATE2, MIN2, TEST)

      write(i2, 150) frame1, month1, iday1, iyear1, date1,
     1               frame2, month2, iday2, iyear2, date2
  150 format(' TRANSFORMING POSITIONS FROM ',A24,' (EPOCH = ',
     1  I2.2,'-',I2.2,'-',I4,' (',F9.4,'))'/26X,
     1  'TO ', A24, ' (EPOCH = ',I2.2,'-',I2.2,'-',I4,' (',F9.4,'))'/)
 
  160 write(luout, 170)
  170 format (/' ***************************************'/
     1  ' Coordinates will be transformed at each specified point.'/
     1  ' Please indicate how you wish to supply positions.'/
     1  '    0... No more points.   Return to main menu.'/
     1  '    1... Individual points entered interactively.'/
     1  '    2... Transform positions contained in batch file '/
     1  '         of delimited records of the form: '/
     1  '         LAT,LON,EHT,TEXT' /
     1  '         LAT = latitude in degrees (positive north/DBL PREC)'/
     1  '         LON = longitude in degrees (positive west/DBL PREC)'/
     1  '         EHT = ellipsoid height in meters (DBL PREC)'/
     1  '         TEXT = Descriptive text (CHARACTER*24) '/
     1  '         Example:   '/
     1  '         40.731671553,112.212671753,34.241,SALT AIR  '/ 
     1  '    3... Transform positions contained in batch file '/
     1  '         of delimited records of the form: '/
     1  '         X, Y, Z,TEXT' /
     1  '         X, Y, Z are Cartesian coordinates of the station'/
     1  '         TEXT = Descriptive text (CHARACTER*24) '/
     1  '         Example:   '/
     1  '           -86682.104,-5394026.861,3391189.647,SALT AIR  '/)
      read(luin,'(a1)',err=602,iostat=ios) option
      if (ios /= 0) goto 602

      if (option .eq. '0') then
         go to 500

      elseif (option .eq. '1') then
         write(i2,175)
  175    format(14X,'INPUT COORDINATES   OUTPUT COORDINATES'/)
     
	 
  180    call GETPNT(latd, latm, slat, LATDIR, lond, lonm, slon, 
     1        LONDIR, name24, x, y, z, ylat, ylon, eht)

         if (min1 .eq. min2) then
             ylat1 = ylat
             ylon1 = ylon
             eht1 = eht
         else
             call PREDV(ylat,ylon,eht,date1,iopt1,
     1            jregn,vn,ve,vu)        
             if(jregn .eq. 0) then
                write(luout,609)
                go to 220
             endif

             call NEWCOR(ylat,ylon,eht,min1,min2,ylat1,ylon1,eht1,
     1            dn,de,du,vn,ve,vu)
         endif
         elon1 = twopi - ylon1
         call TOXYZ(ylat1,elon1,eht1,x1,y1,z1)
         call to_itrf2014(x1,y1,z1,x2,y2,z2,date2,iopt1)
         call from_itrf2014(x2,y2,z2,x3,y3,z3,date2,iopt2)
         

         if(.not.(FRMXYZ(x3,y3,z3,ylatt,elont,ehtnew))) STOP 666
         ylont = -elont
         if(ylont .lt. 0.0d0) ylont = ylont + twopi
               
         call PRNTTP(x, y, z, x3, y3, z3, ylat, ylatt,
     1      ylon, ylont, eht, ehtnew, name24, 1)
     
  220    write(luout,*) ' Transform more positions? (y/n)  '
         read(luin,'(a1)',err=601,iostat=ios) answer
         if (ios /= 0) goto 601
         if(answer .eq. 'Y' .or. answer .eq. 'y') go to 180
         close (i2, status = 'keep')
         
  330 format(/1x,a,' is outside of the modeled region.')      
         
      elseif (option .eq. '2') then
         
         write (luout, 400)
  400    format (' Enter name of input file: ')     
         read (luin,'(a)',err=608,iostat=ios) nameif
         if (ios /= 0) goto 608
         open (i1, file = nameif, status = 'old')

C  write some comments in the transformed files

         call extract_name (nameif,iii)
c        write (I2,1309) nameif(1:iii)
c        write (I2,1309) trim(nameif)
         write (i2,1310) Trans4D_version
         write (i2,1311) frame2       
         WRITE (I2, 1312) MONTH2, IDAY2, IYEAR2, DATE2
         write (i2, 1313)
         
  410    read (i1,'(a)',end=450,err=607,iostat=ios) record
         if (ios /= 0) goto 607
         call interprate_XYZ_record (record,xlat,xlon,eht,name24)
         ylat = (xlat*3600.d0) / rhosec
         ylon = (xlon*3600.d0) / rhosec
        
         if (min1 .eq. min2) then 
            ylat1 = ylat
            ylon1 = ylon
            eht1 = eht
         else
            call PREDV(ylat,ylon,eht,date1,iopt1,
     1           jregn,vn,ve,vu) 
            if(jregn .eq. 0) then
               write(i2, 330) name24
               go to 410
            endif
            call NEWCOR(ylat,ylon,eht,min1,min2,
     1            ylat1,ylon1,eht1,dn,de,du,vn,ve,vu)
          endif
         
         elon1 = -ylon1
         call TOXYZ(ylat1,elon1,eht1,x1,y1,z1)
         call to_itrf2014(x1,y1,z1,x2,y2,z2,date2,iopt1)
         call from_itrf2014(x2,y2,z2,x3,y3,z3,date2,iopt2)
        

         if(.not.FRMXYZ(x3,y3,z3,ylatt,elont,ehtnew)) STOP 666
         ylont = -elont
         if (ylont .lt. 0.d0) ylont = ylont + twopi
         outlat = ylatt*180.d0/pi
         outlon = ylont*180.d0/pi
         call extract_name (name24,iii)
         write (i2,449) outlat,outlon,ehtnew,name24(1:iii)
 449     format (2f15.10,f10.3,4x,a)
         go to 410
      
  450    close(i1, status = 'keep')
         close(i2, status = 'keep')

      elseif (option .eq. '3') then
         
         write (luout, 4001)
 4001    format (' Enter name of input file: ')
         read (luin, '(a)',err=600,iostat=ios) nameif
         if (ios /= 0) goto 600
         open (i1, file = nameif, status = 'old')

C  write some comments in the transformed files

         call extract_name (nameif,iii)
c        write (i2,1309) nameif(1:iii)
c        write (i2,1309) trim(nameif)
         write (i2,1310) Trans4D_version
         write (i2,1311) frame2       
         WRITE (i2, 1312) MONTH2, IDAY2, IYEAR2, DATE2
         write (i2, 1314)
c        const = 180.d0/PI
  411    read (i1,'(a)',end=451,err=607,iostat=ios) record         
         if (ios /= 0) goto 607
         call interprate_XYZ_record (record,x,y,z,name24)
         
         if (min1 .eq. min2) then
             x1 = x
             y1 = y
             z1 = z
         else
            if(.not.FRMXYZ(x,y,z,ylat,elon,eht)) STOP 666
            ylon = -elon
            if(ylon .lt. 0.d0) ylon = ylon + twopi
            call PREDV(ylat,ylon,eht,date1,iopt1,
     1           jregn,vn,ve,vu)
            if (jregn .eq. 0) then
               write(i2,330) name24
               go to 411
            endif
            call NEWCOR(ylat,ylon,eht,min1,min2,
     1           ylat1,ylon1,eht1,dn,de,du,vn,ve,vu) 
            call TOXYZ(ylat1,-ylon1,eht1,x1,y1,z1)
         endif   
     
          call to_itrf2014(x1,y1,z1,x2,y2,z2,date2,iopt1)
          call from_itrf2014(x2,y2,z2,x3,y3,z3,date2,iopt2)
         

         call extract_name (name24,iii)
         write (i2,1449) x3,y3,z3,name24(1:iii)
 1449    format (3f20.3,4x,a)
         go to 411
      
  451    close(i1, status = 'keep')
         close(i2, status = 'keep')

      else
         write(luout,*) ' Improper selection -- try again.  '
         go to 160
      endif

  500 continue
c     close(i2, status = 'keep')
      return

  600 write (*,'(/)') 
      write (*,*) "Wrong output file name in TRFPOS1: ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  601 write (*,'(/)') 
      write (*,*) "Wrong answer in TRFPOS1: ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  602 write (*,'(/)') 
      write (*,*) "Wrong option in TRFPOS1: ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  603 write (*,'(/)') 
      write (*,*) "Wrong bbname in TRFPOS1: ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  604 write (*,'(/)') 
      write (*,*) "Wrong card in TRFPOS1: ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  607 write (*,*) "Failed reading input file in TRFPOS1: ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop
      
  608 write (*,'(/)')
      write (*,*) "Wrong input file name in TRFPOS1: ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop
      
  609 format(/'********************************************'/
     1 'These coordinates cannot be transformed to a different'/
     1 'date because the point is outside of the modeled region.'/)
      stop
     
 1310    format (' ***CAUTION: This file was processed using Trans4D',
     &           ' version ',a10, '***')
     
 1311    format (' ***CAUTION: Coordinates in this file are in ',    
     &           a24, '***')

 1312    FORMAT(' ***CAUTION: Coordinates in this file have been ',
     *   'updated to ',I2,'-',I2.2,'-',I4, '=(',F8.3,') ***'/)
 1313    format(/'  Latitude_(N)   Longitude_(W)    Ellip_Ht  Name')
 1314    format(/16x,'X',19x,'Y',19x,'Z',7x,'Name')

      end
C*******************************************8*****************************************************
      subroutine PRNTTP(x, y, z, x1, y1, z1, ylat,            
     1  ylatt, ylon, ylont, eht, ehtnew, name24, iprint)
     

** Print updated and/or transformed parameters               

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
      CHARACTER*24 NAME24
      CHARACTER*1 LATDIR, LONDIR, LATDR, LONDR

            CALL TODMSS(YLATT,LATDN,LATMN,SLATN,ISIGN)
            LATDIR = 'N'     
            if (isign .eq. -1) LATDIR = 'S'
            CALL TODMSS(YLONT,LONDN,LONMN,SLONN,ISIGN)
            LONDIR = 'W'
            if (isign .eq. -1) LONDIR = 'E'

      if (iprint .eq. 1) then
            WRITE(LUOUT,1065)LATDN,LATMN,SLATN,LATDIR,LONDN,LONMN,
     1                    SLONN,LONDIR, EHTNEW, X1,Y1,Z1
 1065       FORMAT(' ****************************************'/
     1         ' New latitude   = ',I3,1X,I2,1X,F8.5,1X,A1  /
     2         ' New longitude  = ',I3,1x,I2,1X,F8.5,1X,A1  /
     2         ' New Ellip. Ht. = ', F12.3 ,' meters' /
     3         ' New X          = ',F12.3  ,' meters' /
     4         ' New Y          = ',F12.3  ,' meters' /
     5         ' New Z          = ',F12.3  ,' meters' /
     3         ' ****************************************'/)
      endif

            call TODMSS(ylat,latd,latm,slat,isign)
            latdr = 'N'
            if(isign .eq. -1) latdr = 'S'
            call TODMSS(ylon, lond,lonm,slon,isign)
            londr = 'W'
            if(isign .eq. -1) londr = 'E'
            WRITE(I2,1070)NAME24,
     1          LATD,LATM,SLAT,latdr,LATDN,LATMN,SLATN,latdir,
     1          LOND,LONM,SLON,londr,LONDN,LONMN,SLONN,londir,
     1          EHT, EHTNEW, X, X1,
     1          Y, Y1, Z, Z1
 1070       FORMAT(1X,A24,/
     1   2X, 'LATITUDE  ',2(2X,I3,1X,I2.2,1X,F8.5,1X,A1,2X)/
     1   2X, 'LONGITUDE ',2(2X,I3,1X,I2.2,1X,F8.5,1X,A1,2X)/
     1   2X, 'ELLIP. HT.',2(6X,F14.3),' m  '/ 
     1   2X, 'X         ',2(6X,F14.3),' m  ' /
     1   2X, 'Y         ',2(6X,F14.3),' m  ' /
     1   2X, 'Z         ',2(6X,F14.3),' m  '/ )
    
         RETURN
         END

**********************************************************
      subroutine SETTP

*** Specify transformation parameters from ITRF2014
*** to other reference frames
**************************************



      implicit double precision (a-h, o-z)
      implicit integer*4 (i-n)
      parameter (numref = 17)

      common /const/ a, f, e2, eps, af, pi, twopi, rhosec
      common /tranpa/ tx(numref), ty(numref), tz(numref), 
     &                dtx(numref), dty(numref), dtz(numref),
     &                rx(numref), ry(numref), rz(numref), 
     &                drx(numref), dry(numref), drz(numref),
     &                scale(numref), dscale(numref), refepc(numref)

*** From ITRF2014 to NAD 83(2011) or NAD 83(CORS96)
      tx(1) = 1.00530d0        
      ty(1) = -1.9021d0      
      tz(1) = -.54157d0      
      dtx(1) = 0.00079d0        
      dty(1) = -.00060d0       
      dtz(1) = -.00144d0      
      rx(1) = 0.02678138 / rhosec
      ry(1) = -0.00042027 / rhosec
      rz(1) = 0.01093206 / rhosec
      drx(1) = 0.00006667 / rhosec
      dry(1) = -.00075744 / rhosec
      drz(1) = -.00005133 / rhosec
      scale(1) = 0.36891d-9
      dscale(1) = -0.07201d-9
      refepc(1) = 2010.0d0

*** From ITRF2014 to ITRF88
      tx(2) = 0.0254d0
      ty(2) = -.0005d0
      tz(2) = -.1548d0
      dtx(2) = 0.0001d0    
      dty(2) = -.0005d0
      dtz(2) = -.0033d0
      rx(2) = -.0001d0 / rhosec
      ry(2) = 0.d0                 
      rz(2) = -0.00026d0 / rhosec 
      drx(2) = 0.0d0                  
      dry(2) = 0.0d0              
      drz(2) = -.00002d0 / rhosec                  
      scale(2) = 11.29d-9
      dscale(2) = 0.12d-9
      refepc(2) = 2010.0d0

*** From ITRF2014 to ITRF89
      tx(3) = 0.0304d0
      ty(3) = 0.0355d0
      tz(3) = -.1308d0
      dtx(3) = 0.0001d0    
      dty(3) = -.0005d0
      dtz(3) = -.0033d0
      rx(3) = 0.0d0
      ry(3) = 0.0d0                 
      rz(3) = -.00026d0 / rhosec 
      drx(3) = 0.0d0                  
      dry(3) = 0.0d0              
      drz(3) = -.00002d0 / rhosec                   
      scale(3) = 8.19d-9
      dscale(3) = 0.12d-9
      refepc(3) = 2010.0d0

*** From ITRF2014 to ITRF90
      tx(4) = 0.0254d0
      ty(4) = 0.0115d0
      tz(4) = -.0928d0
      dtx(4) = 0.0001d0    
      dty(4) = -.0005d0
      dtz(4) = -.0033d0
      rx(4) = 0.0d0
      ry(4) = 0.0d0                 
      rz(4) = -.00026d0 / rhosec 
      drx(4) = 0.0d0                  
      dry(4) = 0.0d0              
      drz(4) = -.00002d0 / rhosec                   
      scale(4) = 4.79d-9
      dscale(4) = 0.12d-9
      refepc(4) = 2010.0d0

*** From ITRF2014 to ITRF91
      tx(5) = 0.0274d0
      ty(5) = 0.0155d0
      tz(5) = -.0768d0
      dtx(5) = 0.0001d0    
      dty(5) = -.0005d0
      dtz(5) = -.0033d0
      rx(5) = 0.0d0
      ry(5) = 0.0d0                 
      rz(5) = -.000026d0 /rhosec 
      drx(5) = 0.0d0                  
      dry(5) = 0.0d0              
      drz(5) = -.00002d0 / rhosec                   
      scale(5) = 4.49d-9
      dscale(5) = 0.12d-9
      refepc(5) = 2010.0d0

*** From ITRF2014 to ITRF92
      tx(6) = 0.0154d0
      ty(6) = 0.0015d0
      tz(6) = -.0708d0
      dtx(6) = 0.0001d0    
      dty(6) = -.0005d0
      dtz(6) = -.0033d0
      rx(6) = 0.0d0
      ry(6) = 0.0d0                 
      rz(6) = -.00026d0 / rhosec
      drx(6) = 0.0d0                  
      dry(6) = 0.0d0              
      drz(6) = -.00002d0 / rhosec                   
      scale(6) = 3.09d-9
      dscale(6) = 0.12d-9
      refepc(6) = 2010.0d0

*** From ITRF2014 to ITRF93
      tx(7) = -.0504d0
      ty(7) = 0.0033d0
      tz(7) = -.0602d0
      dtx(7) = -.0028d0
      dty(7) = -.0001d0
      dtz(7) = -.0025d0
      rx(7) = 0.00281d0 / rhosec
      ry(7) = 0.00338d0 / rhosec
      rz(7) = -.00040d0 / rhosec
      drx(7) = .00011d0 / rhosec
      dry(7) = .00019d0 / rhosec
      drz(7) =-.00007d0 / rhosec
      scale(7) = 4.29d-9
      dscale(7) = 0.12d-9
      refepc(7) = 2010.0d0

*** From ITRF2014 to ITRF94 and ITRF96
      tx(8) = 0.0074d0
      ty(8) = -.0005d0
      tz(8) = -.0628d0
      dtx(8) = 0.0001d0
      dty(8) = -.0005d0
      dtz(8) = -.0033d0
      rx(8) = 0.d0
      ry(8) = 0.d0
      rz(8) = -.00026d0 / rhosec
      drx(8) = 0.0d0
      dry(8) = 0.d0
      drz(8) = -.00002d0 / rhosec
      scale(8) = 3.80d-9
      dscale(8) = 0.12d-9
      refepc(8) = 2010.0d0

*** From ITRF2014 to ITRF97 
      tx(9) = 0.0074d0
      ty(9) = -.0005d0
      tz(9) = -.0628d0
      dtx(9) = 0.0001d0
      dty(9) = -.0005d0
      dtz(9) = -.0033d0
      rx(9) = 0.0d0 
      ry(9) = 0.0d0 
      rz(9) = -.00026d0 / rhosec
      drx(9) = 0.0d0 
      dry(9) = 0.0d0 
      drz(9) = -0.00002d0 / rhosec
      scale(9) = 3.80d-9
      dscale(9) = 0.12d-9
      refepc(9) = 2010.0d0

*** From ITRF2014 to ITRF2014-PMM for North America
      tx(10) = 0.0d0
      ty(10) = 0.0d0
      tz(10) = 0.0d0
      dtx(10) = 0.0d0
      dty(10) = 0.0d0
      dtz(10) = 0.0d0
      rx(10) = 0.0d0 
      ry(10) = 0.0d0
      rz(10) = 0.0d0
      drx(10) = +0.000024d0 / rhosec
      dry(10) = -0.000694d0 / rhosec
      drz(10) = -0.000063d0 / rhosec
      scale(10) = 0.0d-9
      dscale(10) = 0.0d-9
      refepc(10) = 2010.0d0

*** From ITRF2014 to ITRF2000
      tx(11) = 0.0007d0
      ty(11) = 0.0012d0
      tz(11) = -.0261d0
      dtx(11) = 0.0001d0
      dty(11) = 0.0001d0
      dtz(11) = -0.0019d0
      rx(11) = 0.0d0 
      ry(11) = 0.0d0 
      rz(11) = 0.0d0 
      drx(11) = 0.0d0 
      dry(11) = 0.0d0 
      drz(11) = 0.0d0 
      scale(11) = 2.12d-9
      dscale(11) = 0.11d-9
      refepc(11) = 2010.0d0

*** From ITRF2014 to PACP00 or PA11
*** Based on the rotation rate of the Pacific plate
***   estimated by Beavan (2002)
      tx(12) = 0.9109d0
      ty(12) = -2.0129d0
      tz(12) = -0.5863d0
      dtx(12) = 0.0001d0
      dty(12) = 0.0001d0
      dtz(12) = -.0019d0
      rx(12) = 0.022749d0 / rhosec
      ry(12) = 0.026560d0 / rhosec
      rz(12) = -.025706d0 / rhosec
      drx(12) = -.000344d0 / rhosec
      dry(12) = 0.001007d0 / rhosec
      drz(12) = -.002186d0 / rhosec
      scale(12) = 2.12d-9
      dscale(12) = 0.11d-9
      refepc(12) = 2010.0d0

*** From ITRF2014 to MARP00 or MA11
*** Based on the velocity of GUAM
      tx(13) = 0.9109d0
      ty(13) = -2.0129d0
      tz(13) = -0.5863d0
      dtx(13) = 0.0001d0
      dty(13) = 0.0001d0
      dtz(13) = -.0019d0
      rx(13) = 0.028711d0 / rhosec
      ry(13) = 0.011785d0 / rhosec
      rz(13) = 0.004417d0 / rhosec
      drx(13) = -.000020d0 / rhosec
      dry(13) = 0.000105d0 / rhosec
      drz(13) = -.000347d0 / rhosec
      scale(13) = 2.12d-9
      dscale(13) = 0.11d-9
      refepc(13) = 2010.0d0

*** From ITRF2014 to ITRF2005
      tx(14) = 0.0026d0
      ty(14) = 0.0010d0
      tz(14) = -.0023d0
      dtx(14) = 0.0003d0
      dty(14) = 0.0000d0
      dtz(14) = -.0001d0
      rx(14) = 0.0d0 
      ry(14) = 0.0d0
      rz(14) = 0.0d0 
      drx(14) = 0.0d0
      dry(14) = 0.0d0
      drz(14) = 0.0d0
      scale(14) = 0.92d-9
      dscale(14) = 0.03d-9
      refepc(14) = 2010.0d0

*** From ITRF2014 to ITRF2008 (also IGS08 and IGB08)
      tx(15) = 0.0016d0
      ty(15) = 0.0019d0
      tz(15) = 0.0024d0
      dtx(15) = 0.0d0
      dty(15) = 0.0d0
      dtz(15) = -.0001d0
      rx(15) = 0.0d0 
      ry(15) = 0.0d0 
      rz(15) = 0.0d0
      drx(15) = 0.0d0 
      dry(15) = 0.0d0 
      drz(15) = 0.0d0 
      scale(15) = -0.02d-9
      dscale(15) = 0.03d-9
      refepc(15) = 2010.0d0

*** From ITRF2014 to ITRF2014
      tx(16)     = 0.0d0
      ty(16)     = 0.0d0
      tz(16)     = 0.0d0
      dtx(16)    = 0.0d0                 
      dty(16)    = 0.0d0                 
      dtz(16)    = 0.0d0                 
      rx(16)     = 0.0d0 
      ry(16)     = 0.0d0
      rz(16)     = 0.0d0 
      drx(16)    = 0.0d0     
      dry(16)    = 0.0d0     
      drz(16)    = 0.0d0      
      scale(16)  = 0.0d0
      dscale(16) =  0.0d0               
      refepc(16) =  2010.0d0                 

*** From ITRF2014 to Pre-CATRF2022 (Caribbean IFVM)

      tx(17) = 0.00d0
      ty(17) = 0.00d0
      tz(17) = 0.00d0
      dtx(17) = 0.000d0
      dty(17) = 0.000d0
      dtz(17) = 0.000d0
      rx(17) = 0.000d0 / rhosec
      ry(17) =  0.000d0 / rhosec
      rz(17) =  0.000d0 / rhosec
      drx(17) = -0.000000000351d0 
      dry(17) = -0.000000004522d0 
      drz(17) = +0.000000002888d0 
      scale(17) = 0.000d0
      dscale(17) = 0.000d0
      refepc(17) = 2010.0d0

      return
      end
*************************************************************

      subroutine from_itrf2014(x1, y1, z1, x2, y2, z2, date, jopt)

*** Converts ITRF2014 cartesian coordinates to cartesian
*** coordinates in the specified reference frame for the
*** given date

*** (x1, y1, z1) --> input ITRF94 coordiates (meters)
*** (x2, y2, z2) --> output coordinates (meters)
*** date --> time (decimal years) to which the input & output
***          coordinates correspond
*** jopt --> input specifier of output reference frame

      implicit double precision (a-h, o-z)
      implicit integer*4 (i-n)
      parameter (numref = 17)

      common /tranpa/ tx(numref), ty(numref), tz(numref), 
     &                dtx(numref), dty(numref), dtz(numref),
     &                rx(numref), ry(numref), rz(numref), 
     &                drx(numref), dry(numref), drz(numref),
     &                scale(numref), dscale(numref), refepc(numref)

      if (jopt .eq. 0) then
         iopt = 1
      else
         iopt = jopt
      endif

      dtime = date - refepc(iopt)
      tranx = tx(iopt) + dtx(iopt)*dtime
      trany = ty(iopt) + dty(iopt)*dtime
      tranz = tz(iopt) + dtz(iopt)*dtime
      rotnx  = rx(iopt) + drx(iopt)*dtime
      rotny  = ry(iopt) + dry(iopt)*dtime
      rotnz  = rz(iopt) + drz(iopt)*dtime
      ds     = 1.d0 + scale(iopt) + dscale(iopt)*dtime

      x2 = tranx + ds*x1 + rotnz*y1 - rotny*z1
      y2 = trany - rotnz*x1 + ds*y1 + rotnx*z1
      z2 = tranz + rotny*x1 - rotnx*y1 + ds*z1

      return
      end

***********************************************************

      subroutine to_itrf2014(x1, y1, z1, x2, y2, z2, date, jopt)

*** Converts cartesian coordinates in a specified reference
*** to ITRF2014 cartesian coordinates for the given date

*** (x1, y1, z1) --> input coordiates (meters)
*** (x2, y2, z2) --> output  ITRF2014 coordinates (meters)
*** date --> time (decimal years) to which the input & output
***          coordinates correspond
*** jopt --> input specifier of input reference frame

      implicit double precision (a-h, o-z)
      implicit integer*4 (i-n)
      parameter (numref = 17)

      common /tranpa/ tx(numref), ty(numref), tz(numref), 
     &                dtx(numref), dty(numref), dtz(numref),
     &                rx(numref), ry(numref), rz(numref), 
     &                drx(numref), dry(numref), drz(numref),
     &                scale(numref), dscale(numref), refepc(numref)

      if (jopt .eq. 0) then
         iopt = 1
      else
         iopt = jopt
      endif

      dtime = date - refepc(iopt)
      tranx = -(tx(iopt) + dtx(iopt)*dtime)
      trany = -(ty(iopt) + dty(iopt)*dtime)
      tranz = -(tz(iopt) + dtz(iopt)*dtime)
      rotnx  = -(rx(iopt) + drx(iopt)*dtime)
      rotny  = -(ry(iopt) + dry(iopt)*dtime)
      rotnz  = -(rz(iopt) + drz(iopt)*dtime)
      ds     = 1.d0 - (scale(iopt) + dscale(iopt)*dtime)

      x2 = tranx + ds*x1 + rotnz*y1 - rotny*z1
      y2 = trany - rotnz*x1 + ds*y1 + rotnx*z1
      z2 = tranz + rotny*x1 - rotnx*y1 + ds*z1

      return
      end
*************************************************

      subroutine MENU1(kopt, mframe)

** Write out options for reference frames

      implicit integer*4 (i-n)
      character*24 mframe, nframe
      dimension nframe(24)
      dimension iframe(24)
      common /files/ luin, luout, i1, i2, i3, i4, i5, i6

      iframe(1) = 1
      nframe(1) = 'NAD_83(2011/CORS96/2007)'       
      iframe(2) = 12
      nframe(2) = 'NAD_83(PA11/PACP00)     '
      iframe(3) = 13
      nframe(3) = 'NAD_83(MA11/MARP00)     '
      iframe(4) = 10
      nframe(4) = 'Stable NA (ITRF2014-PMM)'
      iframe(5) = 1
      nframe(5) = 'WGS_84(transit)         '
c     iframe(6) = 6 (This was incorrect in all versions of HTDP)
      iframe(6) = 5
      nframe(6) = 'WGS_84(G730)            '
      iframe(7) = 8
      nframe(7) = 'WGS_84(G873)            '
      iframe(8) = 11
      nframe(8) = 'WGS_84(G1150)           '
      iframe(9) = 15
      nframe(9) = 'WGS_84(G1674)           '
      iframe(10)= 15
      nframe(10)= 'WGS_84(G1762)           '
      iframe(11)= 17
      nframe(11)= 'Pre-CATRF2022 =Caribbean'
      iframe(12)= 2
      nframe(12)= 'ITRF88                  '
      iframe(13)= 3
      nframe(13)= 'ITRF89                  '
      iframe(14)= 4
      nframe(14)= 'ITRF90/PNEOS_90/NEOS_90 '
      iframe(15)= 5
      nframe(15)= 'ITRF91                  '
      iframe(16)= 6
      nframe(16)= 'ITRF92                  '
      iframe(17)= 7
      nframe(17)= 'ITRF93                  '
      iframe(18)= 8
      nframe(18)= 'ITRF94                  '
      iframe(19)= 8
      nframe(19)= 'ITRF96                  '
      iframe(20)= 9
      nframe(20)= 'ITRF97 or IGS97         '
      iframe(21)= 11
      nframe(21)= 'ITRF2000 or IGS00/IGb00 '
      iframe(22)= 14
      nframe(22)= 'ITRF2005 or IGS05       '
      iframe(23)= 15
      nframe(23)= 'ITRF2008 or IGS08/IGb08 '
      iframe(24)= 16
      nframe(24)= 'ITRF2014 or IGS14       '

      write(luout, 100)
  100 format(
     1'  1...NAD_83(2011/CORS96/2007) (for use near North America) '/
     1'  2...NAD_83(PA11/PACP00)      (for use on Pacific islands) '/
     1'  3...NAD_83(MA11/MARP00)      (for use on the Mariana plate) '/
     1'                                                   '/
     1'  4...Relative to Stable North America according to         '/
     1'      the ITRF2014 plate motion model        '/
     1'                                                      '/
     1'  5...WGS_84(transit) (NAD_83(2011) used)  15...ITRF91 '/
     1'  6...WGS_84(G730) (ITRF91 used)           16...ITRF92 '/
     1'  7...WGS_84(G873) (ITRF94 used)           17...ITRF93 '/
     1'  8...WGS_84(G1150) (ITRF2000 used)        18...ITRF94 '/
     1'  9...WGS_84(G1674) (ITRF2008 used)        19...ITRF96 '/
     1' 10...WGS_84(G1762) (IGb08 used)   20...ITRF97 or IGS97'/
     1' 11...Pre-CATRF2022 CARIBBEAN IFVM 21...ITRF2000 or IGS00/IGb00'/
     1' 12...ITRF88                       22...ITRF2005 or IGS05 '/
     1' 13...ITRF89                       23...ITRF2008 or IGS08/IGb08'/
     1' 14...ITRF90 or (PNEOS_90/NEOS_90) 24...ITRF2014 or IGS14      '/
     1     )

      read (luin, *,err=50,iostat=ios) iopt
      if (ios /= 0) goto 50
      if ( 1 .le. iopt .and. iopt .le. 24) then
	     mframe = nframe(iopt)
	     kopt = iframe(iopt)
      else
	     mframe = '                '
	     kopt = iopt
      endif
      return

 50   write (*,*) 'Failed to read option in MENU1:ios=',ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

      end

***************************************************************
      SUBROUTINE trfbb

** Transform BlueBook positions and/or observations to a
** specified date and reference frame.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (numref = 17)
      parameter (nbbdim = 10000)
      CHARACTER*30 OLDBB,NEWBB, NAMEIF
      CHARACTER*24 NAME24
      CHARACTER*1  OPT, ANSWER, BBTYPE, VOPT
      CHARACTER*1 LATDIR, LONDIR, LATDR, LONDR
      character*24 frame1, frame2
      character*10 Trans4D_version   
      LOGICAL TEST
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
      COMMON /VERSION/ Trans4D_version

      WRITE(LUOUT,20)
   20 FORMAT(' ********************************************'/
     1   ' Please enter the time to which the output'/     
     1   ' positions and/or observations are to correspond.')
   15 CALL GETMDY(MONTH2,IDAY2,IYEAR2,DATE2,MIN2,TEST)

** Choosing reference frame for output positions
   35 WRITE(LUOUT,30)
   30 FORMAT(' **************************************************'/
     1   ' Select the reference frame to be used for the output'/
     2   ' positions and/or observations. '/)
      call MENU1( iopt2, frame2)
      IF(IOPT2 .LT. 1 .OR. IOPT2 .GT. numref) THEN
           WRITE(LUOUT,40)
   40      FORMAT(' Improper selection--try again.  ')
           GO TO 35
      ENDIF

  999 WRITE(LUOUT,1000)
 1000 FORMAT(' ******************************'/
     1   ' Select option:'/
     2   '     0...Return to main menu.'/
     5   '     1...Transform positions for Blue Book stations.'/
     6   '     2...Transform values for Blue Book observations.'/
     7   '     3...Transform both the positions for Blue Book',
     8              ' stations'/
     9   '          and the values for Blue Book',
     9              ' observations.'/ )

      READ(LUIN,'(A1)',err=501,iostat=ios) OPT
      if (ios /= 0) goto 501
      IF(OPT .eq. '0') THEN
         RETURN

*** Transforming a Blue Book 

      ELSEIF(OPT .eq. '1' .or. OPT .eq. '2' .or. OPT .eq. '3') THEN

          WRITE(LUOUT,100)
  100     FORMAT(' Enter name of the input B-FILE.'/)
          READ(LUIN,110,err=502,iostat=ios) OLDBB
          if (ios /= 0) goto 502
  110     FORMAT(A30)
          WRITE(LUOUT,120)
  120     FORMAT(' Enter name for the new B-FILE that is to '/
     1       'contain the transformed information.'/)
          READ(LUIN,110,err=502,iostat=ios) NEWBB
          if (ios /= 0) goto 502
          OPEN(I1,FILE = OLDBB , STATUS = 'OLD')

          
             WRITE(LUOUT,122)
  122        FORMAT(' *********************************************'/
     1           ' Enter the time',
     1           ' to which the input positions correspond.'/)
  123       CALL GETMDY(MONTH1,IDAY1,IYEAR1,DATE1,MIN1,TEST)
  124       write(luout,125)
  125       format (/' Enter the reference frame for ', 
     1               'the input positions')

            call MENU1(iopt1, frame1)
             if (iopt1 .lt. 1 .or. iopt1 .gt. numref) then
               write(luout,*) ' Improper selection -- try again.  '
               go to 124
             endif

*** Retrieve geodetic positions from old blue-book file

         call getpo4(iopt1, date1)

*** Create new blue book file

         OPEN(I2,FILE = NEWBB, STATUS = 'UNKNOWN')

         write (i2,127) Trans4D_version
  127    format (' ***CAUTION: This file was processed using Trans4D',
     &           ' version ',a10, '***')
         if (opt .eq. '2') then
            write (i2,128) frame1
            write (i2,129) month1, iday1,iyear1, date1       
  128       format (' ***CAUTION: Coordinates in this file are in ',    
     &           a24, '***')
         elseif (OPT .EQ. '1' .OR. OPT .EQ. '3') THEN
            write(i2,128) frame2
            WRITE(I2, 129) MONTH2, IDAY2, IYEAR2, DATE2
  129       FORMAT(' ***CAUTION: Coordinates in this file have a ',
     *   'reference date of ',I2,'-',I2.2,'-',I4, ' = (',F8.3,') ***')
         endif

         call upbb4(min1,min2,opt,iopt1,iopt2,date2)

         CLOSE(I1, STATUS = 'KEEP')
         CLOSE(I2, STATUS = 'KEEP')

*** Update G-FILE

         IF(OPT .eq. '2' .or. OPT .eq. '3') THEN
  605      WRITE(LUOUT,610)
  610      FORMAT(/' Is there a G-FILE to be transformed? (y/n)  ')
           READ(LUIN,'(A1)',err=500,iostat=ios) ANSWER
           if (ios /= 0) goto 500
           IF(ANSWER .eq. 'N' .or. ANSWER .eq. 'n') THEN
             CONTINUE
           ELSEIF(ANSWER .eq. 'Y' .or. ANSWER .eq. 'y')THEN
             WRITE(LUOUT,620)
  620        FORMAT(' Enter name of old G-FILE to be updated.  ')
             READ(LUIN,'(A30)',err=502,iostat=ios) OLDBB
             if (ios /= 0) goto 502
             WRITE(LUOUT,630)
  630        FORMAT(' Enter name for the new updated G-FILE.  ')
             READ(LUIN,'(A30)',err=502,iostat=ios) NEWBB
             if (ios /= 0) goto 502
             OPEN(I1,FILE=OLDBB,STATUS='OLD')
             OPEN(I2,FILE=NEWBB,STATUS='UNKNOWN')
             WRITE(I2,130) MONTH2, IDAY2, IYEAR2, DATE2
  130        FORMAT(' ***CAUTION: Observations in this file have been ',
     *       'updated to ',I2,'-',I2.2,'-',I4, ' = (',F8.3,') ***')
        
  634        WRITE(LUOUT, 632)
  632        FORMAT(/' ***************************'/
     *              ' To what reference frame should the GPS'/
     *              ' vectors be transformed?'/
     *              '    -1...Do not transform GPS vectors.')
             CALL MENU1(kopt, frame2)
             IF(KOPT .LT. -1 .OR. KOPT .GT. numref) THEN
                WRITE(LUOUT, 40)
                GO TO 634
             ELSEIF (1 .LE. KOPT .AND. KOPT .LE. numref) then
               write( I2, 640) frame2
  640          format( ' ***CAUTION: All GPS interstation vectors',
     1                 ' have been transformed to ', a24, ' ***')
               write( I2, 641) Trans4D_version
  641          format(' ***CAUTION: Observations were transformed using'
     1                ,' Trans4D version ', a10, ' ***')
             ENDIF
         
             call upgfi4(date2, min2, iopt1,kopt,
     *                   month2, iday2, iyear2)

             CLOSE(I1, STATUS = 'KEEP')
             CLOSE(I2, STATUS = 'KEEP')
           ELSE
             WRITE(LUOUT,700)
             GO TO 605
           ENDIF
         ENDIF
c        CLOSE(I4, STATUS = 'DELETE')
         RETURN
      ELSE
         WRITE(LUOUT,700)
  700    FORMAT(' Improper entry !'/)
         GO TO 999
      ENDIF

  500 write (*,'(/)') 
      write (*,*) "Failed to read answer in UPDATE:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  501 write (*,'(/)') 
      write (*,*) "Failed to read OPTION in UPDATE:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  502 write (*,'(/)') 
      write (*,*) "Failed to read file name in UPDATE:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  504 write (*,'(/)') 
      write (*,*) "Failed to read input in UPDATE:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

      END
*******************************************************************
      SUBROUTINE GETPO4(IOPT, DATE)

*** Retrieve geodetic coordinates from the blue book

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (nbbdim = 10000)
      CHARACTER*1  JN,JW
      CHARACTER*4  TYPE
      CHARACTER*80 CARD
c     CHARACTER*6  PIDs

      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
c     COMMON /ARRAYS/ HT(nbbdim),LOC(nbbdim),PIDs(nbbdim)
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
      COMMON /BBINFO/ BLAT(nbbdim),BLON(nbbdim),BEHT(nbbdim),
     *                BVN(nbbdim),BVE(nbbdim),BVU(nbbdim),
     *                K80(nbbdim),K86(nbbdim)

      DO 100 I = 1, nbbdim
          K80(I) = 0 
          K86(I) = 0
          BEHT(I) = 0.d0
  100 CONTINUE
          
c     JREC = 0
  125 READ(I1,130,END=160,err=200,iostat=ios) CARD
      if (ios /= 0) goto 200
  130 FORMAT(A80)
      TYPE = CARD(7:10)
      IF(TYPE .EQ. '*80*') THEN
c        JREC = JREC + 1
         READ(CARD,140,err=201,iostat=ios) ISN,LATD,LATM,SLAT,
     *        JN,LOND,LONM,SLON,JW
         if (ios /= 0) goto 201
  140    FORMAT(BZ,10X,I4,T45,2I2,F7.5,A1,I3,I2,F7.5,A1)
         if (k80(isn) .eq. 0) then 
             k80(isn) = 1
         else
           write(luout,141) isn
  141      format(' B-FILE contains two *80* records for SSN = ',I4,/
     *            ' Trans4D is aborting',/)
           stop
         endif 

         YLAT = (DBLE((LATD*60+LATM)*60)+SLAT)/RHOSEC
         YLON = (DBLE((LOND*60+LONM)*60)+SLON)/RHOSEC
         IF(JN .EQ. 'S') YLAT = -YLAT
         BLAT(isn) = ylat
         IF(JW .EQ. 'E') YLON = -YLON
         BLON(isn) = YLON

	     call PREDV( ylat, ylon, beht(isn), date, iopt,
     1      IDG, vn, ve, vu)
         IF(IDG .eq. 0) THEN
           WRITE(LUOUT,150) CARD
  150      FORMAT('  **WARNING**  The following point is outside of ',
     1     'the modeled region.'/
     2     ' It will be assumed that this point has not moved.'/A80/)
           vn = 0.d0
           ve = 0.d0
           vu = 0.d0
         ENDIF
         bvn(isn) = vn
         bve(isn) = ve
         bvu(isn) = vu

      ELSEIF(TYPE .eq. '*86*') THEN
         IF(CARD(46:52) .ne. '       ') THEN
            READ(CARD,156,err=202,iostat=ios) ISN, EHT
            if (ios /= 0) goto 202
  156       FORMAT(BZ,10X,I4,T46,F7.3)
         ELSE
            READ(CARD,157,err=202,iostat=ios) ISN, OHT, GHT
            if (ios /= 0) goto 202
  157       FORMAT(BZ,10X,I4,T17,F7.3,T36,F7.3)
            EHT = OHT + GHT
         ENDIF
         if (k86(isn) .eq. 0) then
               k86(isn) = 1
               beht(isn) = eht
         else
           write(luout,158) isn
  158      format(' B-FILE contains two *86* records for SSN = ',I4,/
     *            ' Trans4D is aborting'/)
           stop
         endif
      ENDIF
      GO TO 125
  160 REWIND I1
      
      RETURN

  200 write (*,'(/)') 
      write (*,*) "Failed 1st reading in GETPO4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  201 write (*,'(/)') 
      write (*,*) "Failed reading *80* record in GETPO4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  202 write (*,'(/)') 
      write (*,*) "Failed reading *86* record in GETPO4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

      END

******************************************************************
      SUBROUTINE UPBB4(MIN1,MIN2,OPT,IOPT,IOPT2,DATE2)

*** Update blue book

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      parameter (nbbdim = 10000)

      CHARACTER*80 CARD
      CHARACTER*6 DATE
      CHARACTER*4 TYPE
      CHARACTER*1 OPT,JN,JW
      CHARACTER*1 LATDIR, LONDIR
c     CHARACTER*6 PIDs                  
      LOGICAL TEST, TEST1
c     COMMON /ARRAYS/ HT(nbbdim),LOC(nbbdim) ,PIDs(nbbdim)
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /BBINFO/ BLAT(nbbdim),BLON(nbbdim),BEHT(nbbdim),
     *                BVN(nbbdim),BVE(nbbdim),BVU(nbbdim),
     *                K80(nbbdim),K86(nbbdim)

C***    To update classical observations an *12* record
C***    is needed to identify the correct century, as
C***    classical observation records contain only a 2-digit year
      IREC12 = 0

  170 READ(I1,175,END=600,err=700,iostat=ios) CARD
      if (ios /= 0) goto 700
  175 FORMAT(A80)
      TYPE = CARD(7:10)
      IF(TYPE .EQ. '*A1*' .OR.
     1      TYPE .EQ. '*AA*' .OR.
     1      TYPE .EQ. '*10*' .OR.
     1      TYPE .EQ. '*11*' .OR.
     1      TYPE .EQ. '*13*' .OR.
     1      TYPE .EQ. '*21*' .OR.
     1      TYPE .EQ. '*25*' .OR.
     1      TYPE .EQ. '*26*' .OR.
     1      TYPE .EQ. '*27*' .OR.
     1      TYPE .EQ. '*28*' .OR.
     1      TYPE .EQ. '*29*' ) THEN
              CONTINUE
      ELSEIF(TYPE .EQ. '*31*' .OR.
     1      TYPE .EQ. '*40*' .OR.
     1      TYPE .EQ. '*41*' .OR.
     1      TYPE .EQ. '*42*' .OR.
     1      TYPE .EQ. '*45*' .OR.
     1      TYPE .EQ. '*46*' .OR.
     1      TYPE .EQ. '*47*' .OR.
     1      TYPE .EQ. '*70*' .OR.
     1      TYPE .EQ. '*81*' .OR.
     1      TYPE .EQ. '*82*' .OR.
     1      TYPE .EQ. '*83*' .OR.
     1      TYPE .EQ. '*84*' .OR.
     1      TYPE .EQ. '*85*' .OR.
     1      TYPE .EQ. '*90*') THEN
	      CONTINUE

      ELSEIF (TYPE .EQ. '*12*') THEN
         IF( CARD(11:14) .EQ. '    ' .OR.
     1       CARD(17:20) .EQ. '    ') THEN
             WRITE ( LUOUT, 176 )
  176        FORMAT( ' ERROR: The *12* record does not contain'/
     1               '   appropriate dates.')
         ELSE
             READ ( CARD, 177,err=701,iostat=ios) IYEAR1, IYEAR2
             if (ios /= 0) goto 701
  177        FORMAT ( 10X, I4, 2X, I4)
             IF ( (IYEAR2 - IYEAR1) .LT. 0 .OR.
     1            (IYEAR2 - IYEAR1) .GT. 99 ) THEN
                WRITE (LUOUT, 178)
  178           FORMAT(' ERROR: The dates in the *12* record are in'/
     1            ' error. Either they span more than 99 years or '/
     1            ' the end date preceedes the start date.'/
     1            ' If they span more than 99 years, then the Bluebook'/
     1            ' will need to divided into two or more bluebooks'/
     1            ' with the observations in each spanning no more'/
     1            ' than 99 years.')
             ELSE
                IREC12 = 1
             ENDIF
         ENDIF

*** Regular distance

      ELSEIF(TYPE .EQ. '*50*' .OR.
     1       TYPE .EQ. '*51*' .OR.
     1       TYPE .EQ. '*52*') THEN
         IF(OPT .eq. '2' .or. OPT .eq. '3')THEN
           READ(CARD,180,err=702,iostat=ios) ISN,DATE,JSN,OBS
           if (ios /= 0) goto 702
  180      FORMAT(BZ,10X,I4,T35,A6,T46,I4,T64,F9.4)
           CALL CHECK(ISN,JSN,CARD,TEST)
           IF(TEST) THEN
	      CALL TRFDAT(CARD,DATE,IREC12, IYEAR1, IYEAR2, MINO)
	      CALL DSDA(isn, jsn, MINO, MIN2, DS, DA)     
	      IOBS = IDNINT((OBS+DS)*10000.D0)
	      WRITE(CARD(64:72),190) IOBS
  190         FORMAT(I9)
            ENDIF
         ENDIF

*** Azimuth

      ELSEIF(TYPE .EQ. '*60*' .OR. TYPE .EQ. '*61*') THEN
         IF(OPT .eq. '2' .or. OPT .eq. '3') THEN
         READ(CARD,200,err=703,iostat=ios)ISN,DATE,JSN,IDEG,MIN,SEC
         if (ios /= 0) goto 703
  200    FORMAT(BZ,10X,I4,T40,A6,T51,I4,T64,I3,I2,F3.1)
         CALL CHECK(ISN,JSN,CARD,TEST)
         IF(TEST) THEN
	      CALL TRFDAT(CARD,DATE,IREC12, IYEAR1, IYEAR2, MINO)
	      CALL DSDA(ISN, JSN, MINO, MIN2, DS, DA)    
              OBS = (DBLE((IDEG*60+MIN)*60)+SEC)/RHOSEC
	      OBS = OBS + DA
	      IF(OBS .GE. TWOPI) THEN
	        OBS = OBS - TWOPI
              ELSEIF(OBS .LT. 0.D0) THEN
	        OBS = OBS + TWOPI
              ENDIF
	      CALL TODMSS(OBS,IDEG,MIN,SEC,ISIGN)
	      ISEC = IDNINT(SEC*10.D0)
	      WRITE(CARD(64:71),210) IDEG,MIN,ISEC
  210         FORMAT(I3,I2.2,I3.3)
            ENDIF
         ENDIF

*** Direction observation

      ELSEIF(TYPE .EQ. '*20*') THEN
         IF(OPT .eq. '2' .or. OPT .eq. '3') THEN
            READ(CARD,220,err=704,iostat=ios) ISN0,LIST0,DATE,JSN
            if (ios /= 0) goto 704
  220       FORMAT(BZ,10X,I4,I2,T40,A6,T51,I4)
            CALL CHECK(ISN0,JSN,CARD,TEST) 
            CALL TRFDAT(CARD,DATE,IREC12, IYEAR1, IYEAR2, MINO)
            IF(TEST) THEN
	      CALL DSDA(ISN0, JSN, MINO, MIN2, DS, DA0)     
            ELSE
              DA0 = 0.D0
            ENDIF
         ENDIF
      ELSEIF(TYPE .EQ. '*22*') THEN
         IF(OPT .eq. '2' .or. OPT .eq. '3') THEN
	   READ(CARD,230,err=705,iostat=ios) ISN,LIST,JSN,IDEG,
     &                                       MIN,SEC
           if (ios /= 0) goto 705
  230      FORMAT(BZ,10X,I4,I2,T51,I4,T64,I3,I2,F4.2)
	   IF(LIST.NE.LIST0 .OR. ISN.NE.ISN0) THEN
	        WRITE(LUOUT,240) CARD
  240           FORMAT(' Blue-book file has incorrect structure.'/
     1          ' A *22* record disagrees with its corresponding'/
     2          ' *20* record.  The record reads:'/A80)
	        STOP
           ENDIF
           CALL CHECK(ISN,JSN,CARD,TEST)
           IF(TEST) THEN
	      CALL DSDA(ISN, JSN, MINO, MIN2, DS, DA)    
	      OBS=(DBLE((IDEG*60+MIN)*60)+SEC)/RHOSEC
	      OBS = OBS + DA - DA0
	      IF(OBS .GT. TWOPI) THEN
	        OBS = OBS - TWOPI
              ELSEIF(OBS .LT. 0.D0) THEN
	        OBS = OBS + TWOPI
              ENDIF
	      CALL TODMSS(OBS,IDEG,MIN,SEC,ISIGN)
	      ISEC = IDNINT(SEC*100.D0)
	      WRITE(CARD(64:72),250) IDEG,MIN,ISEC
  250         FORMAT(I3.3,I2.2,I4.4)
            ENDIF
         ENDIF

*** Long distance

      ELSEIF(TYPE .EQ. '*53*' .OR.
     1       TYPE .EQ. '*54*') THEN
         IF(OPT .eq. '2' .or. OPT .eq. '3') THEN
           READ(CARD,260,err=706,iostat=ios) ISN,DATE,JSN,OBS
           if (ios /= 0) goto 706
  260      FORMAT(BZ,10X,I4,T35,A6,T46,I4,T64,F10.3)
           CALL CHECK(ISN,JSN,CARD,TEST)
           IF(TEST) THEN
	      CALL TRFDAT(CARD,DATE,IREC12, IYEAR1, IYEAR2, MINO)
              CALL DSDA(ISN, JSN, MINO, MIN2, DS, DA)    
	      IOBS = IDNINT((OBS + DS)*1000.D0)
	      WRITE(CARD(64:73),270) IOBS
  270         FORMAT(I10)
            ENDIF
         ENDIF
  
*** Horizontal angle

      ELSEIF(TYPE .EQ. '*30*' .OR.
     1       TYPE .EQ. '*32*') THEN
         IF(OPT .eq. '2' .or. OPT .eq. '3') THEN
           READ(CARD,280,err=707,iostat=ios) ISN,JSN,IDEG,MIN,SEC,KSN
           if (ios /= 0) goto 707
  280      FORMAT(BZ,10X,I4,T51,I4,T64,I3,I2,F3.1,I4)
           IF(TYPE.EQ.'*30*') THEN
		  DATE = CARD(40:45)
		  CALL TRFDAT(CARD,DATE,IREC12, IYEAR1, IYEAR2, MINO)
           ENDIF
           CALL CHECK(ISN,JSN,CARD,TEST)
           CALL CHECK(ISN,KSN,CARD,TEST1)
           IF(TEST .and. TEST1) THEN
              CALL DSDA(ISN, JSN, MINO, MIN2, DS, DA0)   
              CALL DSDA(ISN, KSN, MINO, MIN2, DS, DA)     
              OBS=(DBLE((IDEG*60+MIN)*60)+SEC)/RHOSEC
              OBS = OBS + DA - DA0
              IF(OBS .GE. TWOPI) THEN
		   OBS = OBS - TWOPI
              ELSEIF(OBS .LT. 0.D0) THEN
		   OBS = OBS + TWOPI
              ENDIF
              CALL TODMSS(OBS,IDEG,MIN,SEC,ISIGN)
              ISEC = IDNINT(SEC*10.D0)
              WRITE(CARD(64:71),290)IDEG,MIN,SEC
  290         FORMAT(I3.3,I2.2,I3.3)
            ENDIF
         ENDIF

*** position record

      ELSEIF(TYPE .EQ. '*80*') THEN
            IF(OPT .eq. '1' .or. OPT .eq. '3') THEN
              READ(CARD,300,err=708,iostat=ios) ISN,LATD,LATM,
     &                              SLAT,JN,LOND,LONM,SLON,JW
              if (ios /= 0) goto 708
  300         FORMAT(BZ,10X,I4,T45,2I2,F7.5,A1,I3,I2,F7.5,A1)
              YLAT = (DBLE((LATD*60+LATM)*60)+SLAT)/RHOSEC
              YLON = (DBLE((LOND*60+LONM)*60)+SLON)/RHOSEC
              IF(JN .EQ. 'S') YLAT = -YLAT
              IF(JW .EQ. 'E') YLON = -YLON

              eht = beht(isn)
              vn = bvn(isn)
              ve = bve(isn)
              vu = bvu (isn) 
	      call NEWCOR( ylat, ylon, eht, min1, min2,
     1           ylatt, ylont, ehtnew, dn, de, du, vn, ve, vu)
              elont = -ylont
              call tran_frames(ylatt,elont,ehtnew,iopt,
     1             glat2,glon2,eht2,iopt2,date2)
              glon2 = - glon2
              if (glon2 .gt. 360.d0) glon2 = glon2 - 360.d0
              
              CALL TODMSS(glat2,LATDN,LATMN,SLATN,ISIGN)
	      LATDIR = 'N'
	      IF (ISIGN .eq. -1) LATDIR = 'S'
              CALL TODMSS(glon2,LONDN,LONMN,SLONN,ISIGN)
	      LONDIR = 'W'
	      IF (ISIGN .eq. -1) LONDIR = 'E'
              LATS = IDNINT(SLATN*100000.D0)
              LONS = IDNINT(SLONN*100000.D0)
              WRITE(CARD(45:56),310) LATDN,LATMN,LATS,LATDIR
  310         FORMAT(I2,I2.2,I7.7,A1)
              WRITE(CARD(57:69),320) LONDN,LONMN,LONS,LONDIR
  320         FORMAT(I3,I2.2,I7.7,A1)
            ENDIF

*** height record

      ELSEIF(TYPE .EQ. '*86*') THEN
             IF(OPT .eq. '1' .or. OPT .eq. '3') THEN
               READ(CARD,305,err=709,iostat=ios) isn,oht
               if (ios /= 0) go to 709
  305          format(BZ,10x,i4,2x,f7.3)
               
               ylat = blat(isn)
               ylon = blon(isn)
               eht  = beht(isn)
               vn   = bvn(isn)
               ve   = bve(isn)
               vu   = bvu(isn)
               call newcor(ylat,ylon,eht,min1,min2,ylatt,ylont,
     *                ehtnew,dn,de,du,vn,ve,vu)
               elont = -ylont
               call tran_frames(ylatt,elont,ehtnew,iopt,
     *              glat2,glon2,eht2,iopt2,date2)               
               write(card(46:52),'(i7)') nint(eht2*1000)
             ENDIF

*** Unrecognized blue book record

      ELSE
c           WRITE(LUOUT,500) TYPE
  500       FORMAT(' This software does not recognize an ',A4,/
     1       ' record.  The record will be copied to the new'/
     2       ' file without change.')
      ENDIF
      WRITE(I2,175) CARD
      GO TO 170
  600 CONTINUE
      RETURN

  700 write (*,'(/)') 
      write (*,*) "Failed to read card in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  701 write (*,'(/)') 
      write (*,*) "Failed to read card *12* in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  702 write (*,'(/)') 
      write (*,*) "Failed to read *50,51,52* in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  703 write (*,'(/)') 
      write (*,*) "Failed to read *60,61* in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  704 write (*,'(/)') 
      write (*,*) "Failed to read *20* in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  705 write (*,'(/)') 
      write (*,*) "Failed to read *22* in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  706 write (*,'(/)') 
      write (*,*) "Failed to read *53,54* in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  707 write (*,'(/)') 
      write (*,*) "Failed to read *30,32* in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  708 write (*,'(/)') 
      write (*,*) "Failed to read *80* in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  709 write (*,'(/)')
      write (*,*) "FAILED to read *86* in UPBB4:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

      END
***************************************************************
      
      SUBROUTINE CHECK(ISN,JSN,CARD,TEST)

*** Check if stations have *80* records

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (nbbdim = 10000)
      CHARACTER*80 CARD
c     CHARACTER*6  PIDs
      LOGICAL TEST
c     COMMON /ARRAYS/ HT(nbbdim),LOC(nbbdim),PIDs(nbbdim)
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
      COMMON /BBINFO/ BLAT(nbbdim),BLON(nbbdim),BEHT(nbbdim),
     *                BVN(nbbdim),BVE(nbbdim),BVU(nbbdim),
     *                K80(nbbdim),K86(nbbdim)

      TEST = .TRUE.

      IF(k80(isn) .eq. 0) THEN
          WRITE(LUOUT,100) ISN, CARD
  100     FORMAT(' No *80* record for SSN = ',I4/A80/)
          TEST = .FALSE.
      ENDIF

      IF(k86(isn) .eq. 0) THEN
          WRITE(LUOUT,110) isn, CARD
  110     format(' No *86* record for SSN = ',i4/a80/)        
          TEST = .FALSE.
      ENDIF

      if(k80(jsn) .eq. 0) then
         write(luout,100) jsn, card
         test = .false.
      endif

      if(k86(jsn) .eq. 0) then
         write(luout,110) jsn, card
         test = .false.
      endif

      RETURN
      END
C*************************************************
      SUBROUTINE UPGFI4(DATE2, MIN2, IOPT, KOPT,
     *     MONTH, IDAY, IYEAR)

*** Update G-FILE of blue book

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      CHARACTER*120 CARD
      CHARACTER*1 TYPE
      CHARACTER*2 NRF
      CHARACTER*2 ZT
      CHARACTER*14 CHAR14
      LOGICAL TEST
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6

      ZT = 'ZT'

*** Obtain blue-book reference frame identifier
*** corresponding to IOPT
      IF (KOPT .NE. -1) THEN
         CALL RFCON1(KOPT, NBBREF)
         WRITE(NRF,10) NBBREF
   10    FORMAT(I2.2)
      ENDIF

   90 READ(I1,100,END=200,err=300,iostat=ios) CARD
      if (ios /= 0) goto 300
c 100 FORMAT(A80)
  100 FORMAT(A120)
      TYPE = CARD(1:1)
      IF(TYPE .eq. 'A') THEN
         IF (CARD(79:80) .EQ. ZT) THEN
            WRITE(LUOUT, 103)
  103       FORMAT(
     1     ' *****************************************************'/
     1     ' * ERROR: The input GFILE contains the letters, ZT,  *'/
     1     ' * in columns 79-80 of its A-record.  This indicates *'/
     1     ' * that the GPS vectors in this file have already    *'/
     1     ' * been updated to a common date.  This software     *'/
     1     ' * will not further modify the GPS vectors.          *'/
     1     ' *****************************************************'/)
            RETURN
         ENDIF
         WRITE(CARD(79:80),101) ZT
  101    FORMAT(A2)
         WRITE(CARD(4:11), 102) IYEAR, MONTH, IDAY
  102    FORMAT(I4, I2.2, I2.2)
         WRITE(CARD(12:19), 102) IYEAR, MONTH, IDAY
                  
      ELSEIF(TYPE .eq. 'B') THEN
         READ(CARD,110,err=301,iostat=ios)IYEAR1,MONTH1,IDAY1,
     1        IYEAR2,MONTH2,IDAY2,IBBREF
         if (ios /= 0) goto 301
  110    FORMAT(1X,I4,I2,I2,4X,I4,I2,I2,30X,I2)

         CALL IYMDMJ(IYEAR1,MONTH1,IDAY1,MJD1)
         MINO1 = MJD1 * 24 * 60
         CALL IYMDMJ(IYEAR1, 1, 1, MJD0)
         DECYR1 = DBLE(IYEAR1) + DBLE(MJD1 - MJD0)/365.D0
         CALL IYMDMJ(IYEAR2,MONTH2,IDAY2, MJD2)
         MINO2 = MJD2 * 24 * 60
         CALL IYMDMJ(IYEAR2, 1, 1, MJD0)
         DECYR2 = DBLE(IYEAR2) + DBLE(MJD2 - MJD0)/365.D0
         MINO = (MINO1 + MINO2) / 2
	     DECYR = (DECYR1 + DECYR2) / 2.D0
         CALL RFCON(IBBREF, JREF)
         IF (KOPT .NE. -1) THEN
            CARD(52:53) = NRF
         ENDIF
      ELSEIF(TYPE .eq. 'C') THEN
         READ(CARD,120,err=302,iostat=ios)ISN,JSN,DX,DY,DZ
         if (ios /= 0) goto 302
  120    FORMAT(BZ,1X,2I4,F11.4,5X,F11.4,5X,F11.4)
   
         CALL CHECK(ISN, JSN, CARD, TEST)
         IF (TEST) THEN
            CALL DDXYZ(ISN, JSN,                                 
     1              MINO, MIN2, DDX, DDY, DDZ)

            
*** Convert vector from JREF to IOPT frame
                                 
            call TRAVEC( DX, DY, DZ, DECYR, JREF, IOPT)
            DX = DX + DDX
            DY = DY + DDY
            DZ = DZ + DDZ

*** Convert GPS vector from IOPT to KOPT reference frame
            IF (KOPT .NE. -1) THEN
              CALL TRAVEC(DX, DY, DZ, DATE2, IOPT, KOPT)
            ELSE
              CALL TRAVEC(DX, DY, DZ, DATE2, IOPT, JREF)
            ENDIF

*** Rewrite GPS observational record
            CALL TOCHAR(DX,CHAR14)          
            CARD(10:20) = CHAR14(3:13)          

            CALL TOCHAR(DY,CHAR14)
            CARD(26:36) = CHAR14(3:13)    

            CALL TOCHAR(DZ,CHAR14)
            CARD(42:52) = CHAR14(3:13)    
         ENDIF

      ELSEIF(TYPE .eq. 'F') THEN
         READ(CARD,140,err=303,iostat=ios)ISN,JSN,DX,DY,DZ
         if (ios /= 0) goto 303
  140    FORMAT(BZ,1X,2I4,F13.4,5X,F13.4,5X,F13.4)

         CALL CHECK(ISN, JSN, CARD, TEST)
         IF (TEST) THEN
            CALL DDXYZ(ISN, JSN,                               
     1              MINO, MIN2, DDX, DDY, DDZ)

            call TRAVEC( DX, DY, DZ, DECYR, JREF, IOPT)

            DX = DX + DDX
            DY = DY + DDY
            DZ = DZ + DDZ

*** Convert GPS vector to output reference frame
            IF (KOPT .NE. -1) THEN
               
                 CALL TRAVEC(DX, DY, DZ, DATE2, IOPT, KOPT)
               
            ELSE
               
	             CALL TRAVEC(DX, DY, DZ, DATE2, IOPT, JREF)
               
            ENDIF
            
*** Rewrite GPS observational record
            CALL TOCHAR(DX,CHAR14)
            CARD(10:22) = CHAR14(1:13)

            CALL TOCHAR(DY,CHAR14)
            CARD(28:40) = CHAR14(1:13)           

            CALL TOCHAR(DZ,CHAR14)
            CARD(46:58) = CHAR14(1:13)
         ENDIF
      ENDIF
      WRITE(I2,100) CARD
      GO TO 90
  200 CONTINUE
      RETURN

  300 write (*,'(/)') 
      write (*,*) "Failed to read gfile in UPGFIG:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  301 write (*,'(/)') 
      write (*,*) "Failed to read B card in UPGFIG:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  302 write (*,'(/)') 
      write (*,*) "Failed to read C card in UPGFIG:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  303 write (*,'(/)') 
      write (*,*) "Failed to read F card in UPGFIG:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

      END
*********************************************************

      SUBROUTINE TODMSS(val,id,im,s,isign)
 
*** convert position radians to deg,min,sec
*** range is [-twopi to +twopi]
 
      implicit double precision(a-h,o-z)
      IMPLICIT INTEGER*4 (I-N)
      common/CONST/A,F,E2,EP2,AF,PI,TWOPI,RHOSEC
 
    1 if(val.gt.twopi) then
        val=val-twopi
        go to 1
      endif
 
    2 if(val.lt.-twopi) then
        val=val+twopi
        go to 2
      endif
 
      if(val.lt.0.d0) then
        isign=-1
      else
        isign=+1
      endif
 
      s=dabs(val*RHOSEC/3600.D0)
      id=idint(s)
      s=(s-id)*60.d0
      im=idint(s)
      s=(s-im)*60.d0
 
*** account for rounding error
 
      is=idnint(s*1.d5)
      if(is.ge.6000000) then
        s=0.d0
        im=im+1
      endif
      if(im.ge.60) then
        im=0
        id=id+1
      endif
 
      return
      end
*********************************************************
      SUBROUTINE SETRF

*** Specify arrays that may be used to convert
*** a reference frame identifier in the blue book
*** to a reference frame identifier in TRANS4D and back

      IMPLICIT INTEGER*4 (I-N)
      parameter ( numref = 17 )
      COMMON /REFCON/ IRFCON(29), JRFCON(numref)

*** From blue book identifier to HTDP indentifier
*** WGS 72 Precise
      IRFCON(1) = 1

*** WGS 84 (orig) Precise (set  equal to NAD 83)
      IRFCON(2) = 1

*** WGS 72 Broadcast
      IRFCON(3) = 1

*** WGS 84 (orig) Broadcast (set equal to NAD 83)
      IRFCON(4) = 1

*** ITRF89
      IRFCON(5) = 3

*** PNEOS 90 or NEOS 91.25 (set equal to ITRF90)
      IRFCON(6) = 4

*** NEOS 90 (set equal to ITRF90)
      IRFCON(7) = 4

*** ITRF91
      IRFCON(8) = 5

*** SIO/MIT 92.57 (set equal to ITRF91)
      IRFCON(9) = 5

*** ITRF91
      IRFCON(10) = 5

*** ITRF92
      IRFCON(11) = 6

*** ITRF93
      IRFCON(12) = 7

*** WGS 84 (G730) Precise (set equal to ITRF91)
      IRFCON(13) = 5

*** WGS 84 (G730) Broadcast (set equal to ITRF91)
      IRFCON(14) = 5

*** ITRF94
      IRFCON(15) = 8

*** WGS 84 (G873) Precise  (set equal to ITRF94)
      IRFCON(16) = 8

*** WGS 84 (G873) Broadcast (set equal to ITRF94)
      IRFCON(17) = 8

*** ITRF96
      IRFCON(18) = 8

*** ITRF97
      IRFCON(19) = 9

*** IGS97
      IRFCON(20) = 9

*** ITRF00
      IRFCON(21) = 11

*** IGS00
      IRFCON(22) = 11

*** WGS 84 (G1150)
      IRFCON(23) = 11

*** IGb00
      IRFCON(24) = 11

*** ITRF2005
      IRFCON(25) = 14

*** IGS05
      IRFCON(26) = 14

*** ITRF2008 or IGS08
      IRFCON(27) = 15

*** IGb08
      IRFCON(28) = 15

*** ITRF2014
      IRFCON(29) = 16

*** From HTDP identifier to blue book identifier
*** NAD 83 (set equal to WGS 84 (transit))
      JRFCON(1) = 2

*** ITRF88 (set equal to ITRF89)
      JRFCON(2) = 5

*** ITRF89
      JRFCON(3) = 5

*** ITRF90 (set equal to NEOS 90)
      JRFCON(4) = 7

*** ITRF91
      JRFCON(5) = 8

*** ITRF92
      JRFCON(6) = 11

*** ITRF93
      JRFCON(7) = 12

*** ITRF96 (= ITRF94)
      JRFCON(8) = 18

*** ITRF97
      JRFCON(9) = 19

*** NA12
      JRFCON(10) = 0

*** ITRF00
      JRFCON(11) = 21

*** NAD 83(PACP00) or NAD 83(PA11)
      JRFCON(12) = 2

*** NAD 83(MARP00) or NAD 83(MA11)
      JRFCON(13) = 2

*** ITRF2005 or IGS05
      JRFCON(14) = 26

*** ITRF2008 or IGS08/IGb08
      JRFCON(15) = 27

*** ITRF2014
      JRFCON(16) = 29

*** NA_ICE-6G
      JRFCON(17) = 0


      RETURN
      END
***************************************************
      SUBROUTINE RFCON(IBBREF, JREF)

*** Convert reference frame identifier from
*** system used in the blue-book to the
*** system used in HTDP

      IMPLICIT INTEGER*4 (I-N)
      parameter ( numref = 17 )
      COMMON /REFCON/ IRFCON(29), JRFCON(numref)

      IF (1 .LE. IBBREF .AND. IBBREF .LE. 29) THEN
          JREF = IRFCON(IBBREF)
      ELSE
          WRITE(6, 10) IBBREF
   10     FORMAT(' Improper reference frame identifier (=',
     1      I4, ')' /
     1      ' appearing in B-record of the G-FILE')
          STOP
       ENDIF

       RETURN
       END
*******************************************************

      SUBROUTINE RFCON1(JREF, IBBREF)

*** Convert reference frame identifier from
*** system used in HTDP to the  system
*** used in the blue-book

      IMPLICIT INTEGER*4 (I-N)
      parameter ( numref = 17 )
      COMMON /REFCON/ IRFCON(29), JRFCON(numref)

      IF (JREF .EQ. 0) THEN
         I = 1
      ELSE
         I = JREF
      ENDIF

      IF(1. LE. I .AND. I .LE. numref) THEN
          IBBREF = JRFCON(I)
	  IF ( IBBREF .eq. 0) THEN
	     write(6,5) JREF
    5        format(' ERROR: The BlueBook does not recognize '/
     *              'this reference frame, TRANS4D ID = ',I2)
	     stop
          ENDIF

      ELSE
          WRITE(6, 10) JREF
   10     FORMAT(' Improper reference frame identifier (=',
     *       I4, ')' / 'appearing in routine RFCON1')
          STOP
      ENDIF

      RETURN
      END
*****************************************************************

      subroutine TRAVEC(dxi, dyi, dzi, date, jopt1, jopt2)

*** Transform GPS or other vector from one reference frame to another
*** for the given date                                     

*** (dxi, dyi, dzi) --> (input) components of input vector in meters
***                 --> (output) components of transformed vector in meters

*** date --> (input) time (decimal years) to which the input & output
***          vectors correspond

*** jopt1 --> (input) specifier of reference frame for input vector
*** jopt2 --> (input) specifier of reference frame for output vector

      implicit double precision (a-h, o-z)
      implicit integer*4 (i-n)
      parameter (numref = 17)

      common /tranpa/ tx(numref), ty(numref), tz(numref), 
     &                dtx(numref), dty(numref), dtz(numref),
     &                rx(numref), ry(numref), rz(numref), 
     &                drx(numref), dry(numref), drz(numref),
     &                scale(numref), dscale(numref), refepc(numref)
 
*** Transform input vector to ITRF2014 reference frame
      
      iopt = jopt1
           
      dtime = date - refepc(iopt)
      rotnx  = -(rx(iopt) + drx(iopt)*dtime)
      rotny  = -(ry(iopt) + dry(iopt)*dtime)
      rotnz  = -(rz(iopt) + drz(iopt)*dtime)

      ds     = - (scale(iopt) + dscale(iopt)*dtime)
    
      dxt = dxi + ds*dxi + rotnz*dyi - rotny*dzi
      dyt = dyi - rotnz*dxi + ds*dyi + rotnx*dzi
      dzt = dzi + rotny*dxi - rotnx*dyi + ds*dzi
     
*** Transform ITRF2014 vector to new reference frame
      
      iopt = jopt2
      
      dtime = date - refepc(iopt)
      rotnx  = rx(iopt) + drx(iopt)*dtime
      rotny  = ry(iopt) + dry(iopt)*dtime
      rotnz  = rz(iopt) + drz(iopt)*dtime

      ds     = scale(iopt) + dscale(iopt)*dtime
     
      dxi = dxt + ds*dxt + rotnz*dyt - rotny*dzt
      dyi = dyt - rotnz*dxt + ds*dyt + rotnx*dzt
      dzi = dzt + rotny*dxt - rotnx*dyt + ds*dzt
      
      return
      end

***********************************************************
      SUBROUTINE GETPNT(LATD,LATM,SLAT,LATDIR,LOND,LONM,SLON,
     1     LONDIR, NAME, X,Y,Z,XLAT,XLON,EHT)
   
*** Interactively obtain name and coordinates for a point.
*** Output longitude will be positive west.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*24 NAME
      CHARACTER*1 COPT,LATDIR,LONDIR
      LOGICAL FRMXYZ
      COMMON /FILES/ LUIN,LUOUT, I1, I2, I3, I4, I5, I6
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC

      WRITE(LUOUT,100)
  100 FORMAT(' Enter name for point (24 character max).  ')
      READ(LUIN,105,err=200,iostat=ios) NAME
      if (ios /= 0) goto 200
  105 FORMAT(A24)

  110 WRITE(LUOUT,111)
  111 FORMAT(' How do you wish to specify positional coordinates:'/
     1       '     1...geodetic latitude, longitude, ellipsoid height'/
     2       '     2...Cartesian (X,Y,Z) coordinates.  ')
      READ(LUIN,'(A1)',err=201,iostat=ios) COPT
      if (ios /= 0) goto 201

      IF(COPT .EQ. '1') THEN
	  WRITE(LUOUT,115)
  115     FORMAT(
     1    ' Enter latitude degrees-minutes-seconds in free format'/,
     2    ' with north being positive. For example,    35,17,28.3  '/
     2    ' For a point in the southern hemisphere, enter a minus sign'/
     3    ' before each value. For example, -35,-17,-28.3')
	  READ(LUIN,*,err=202,iostat=ios) LATD, LATM, SLAT
          if (ios /= 0) goto 202
	  WRITE(LUOUT,120)
  120     FORMAT(
     1    ' Enter longitude degrees-minutes-seconds in free format'/,
     2    ' with west being positive.  To express a longitude measured'/
     3    ' eastward, enter a minus sign before each value.')
	  READ(LUIN,*,err=203,iostat=ios) LOND, LONM, SLON
          if (ios /= 0) goto 203
          WRITE(LUOUT, 125)
  125     FORMAT(
     1    ' Enter ellipsoid height in meters. (Note that'/,
     1    ' predicted motions are independent of this height.)  ')
          READ(LUIN,*,err=204,iostat=ios) EHT
          if (ios /= 0) goto 204
          XLAT =  (DBLE((LATD*60 + LATM)*60) + SLAT)/RHOSEC
	  LATDIR = 'N'
	  IF (XLAT .lt. 0.0D0) then
	     LATD = - LATD
	     LATM = - LATM
	     SLAT = - SLAT
	     LATDIR = 'S'
          ENDIF
          XLON = (DBLE((LOND*60 + LONM)*60) + SLON)/RHOSEC
	  ELON = -XLON
          CALL TOXYZ(XLAT,ELON,EHT,X,Y,Z)
	  LONDIR = 'W'
	  IF (XLON .lt. 0.0D0) then
	     LOND = -LOND
	     LONM = -LONM
	     SLON = -SLON
	     LONDIR = 'E'
          ENDIF

      ELSEIF(COPT .EQ. '2') THEN
          WRITE(LUOUT,130)
  130     FORMAT(' Enter X coordinate in meters.  ')
          READ(LUIN,*,err=205,iostat=ios) X
          if (ios /= 0) goto 205
          WRITE(LUOUT,140)
  140     FORMAT(' Enter Y coordinate in meters.  ')
          READ(LUIN,*,err=206,iostat=ios) Y
          if (ios /= 0) goto 206
          WRITE(LUOUT,150)
  150     FORMAT(' Enter Z coordinate in meters.  ')
          READ(LUIN,*,err=207,iostat=ios) Z
          if (ios /= 0) goto 207
          IF(.NOT.FRMXYZ(X,Y,Z,XLAT,XLON,EHT)) STOP 666
          XLON = -XLON
          IF(XLON .LT. 0.0D0) XLON = XLON + TWOPI
          CALL TODMSS(XLAT,LATD,LATM,SLAT,ISIGN)
	  LATDIR = 'N'
	  IF (ISIGN .eq. -1) LATDIR = 'S'
          CALL TODMSS(XLON,LOND,LONM,SLON,ISIGN)
	  LONDIR = 'W'
	  IF (ISIGN .eq. -1) LONDIR = 'E'

      ELSE
          WRITE(LUOUT,*) ' Improper response -- try again.  '
          GO TO 110
      ENDIF

      RETURN

  200 write (*,'(/)') 
      write (*,*) "Failed to read point name: ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  201 write (*,'(/)') 
      write (*,*) "Failed to read Coord. form option:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  202 write (*,'(/)') 
      write (*,*) "Failed to read latitude:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  203 write (*,'(/)') 
      write (*,*) "Failed to read longitude:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  204 write (*,'(/)') 
      write (*,*) "Failed to read ellipsoidal height:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  205 write (*,'(/)') 
      write (*,*) "Failed to read the X coordinate:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  206 write (*,'(/)') 
      write (*,*) "Failed to read the Y coordinate:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  207 write (*,'(/)') 
      write (*,*) "Failed to read the Z coordinate:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

      END
*****************************************************************
      SUBROUTINE GETVLY(GLAT, GLON, VX, VY, VZ,
     1      VNORTH, VEAST, VUP, VOPT, IFORM)

*** Interactively optain velocity for a point

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 VOPT
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6

  200 CONTINUE
      IF (IFORM .eq. 210) then
         WRITE(LUOUT, 210)
  210    FORMAT(' What is the velocity of this site relative to '/ 
     1       ' the input reference frame:'/ 
     1       '     0...use velocity provided by this software.'/ 
     1       '     1...user will specify the north-east-up'/ 
     1       '         components of site velocity.'/ 
     1       '     2...user will specify the global x-y-z'/
     1       '         components of site velocity.' )
      else
	     write(luout, 211)
  211    FORMAT(' How do you wish to specify the velocity: '/
     1       '     1...north-east-up components.'/
     1       '     2...global x-y-z components.' )
      endif

      READ(LUIN, '(A1)',err=300,iostat=ios) VOPT
      if (ios /= 0) goto 300

      IF (VOPT .EQ. '0') THEN
	     CONTINUE

      ELSEIF (VOPT .EQ. '1') THEN
	     WRITE(LUOUT, 220)
  220    FORMAT( ' Enter north-south component of velocity in mm/yr'/
     1           ' with north being positive and south being negative')
	     READ(LUIN, *,err=301,iostat=ios) VNORTH
         if (ios /= 0) goto 301
	     WRITE(LUOUT, 221)
  221    FORMAT( ' Enter east-west component of velocity in mm/yr'/
     1           ' with east being positive and west being negative')
	     READ(LUIN, *,err=302,iostat=ios) VEAST
         if (ios /= 0) goto 302
	     WRITE(LUOUT, 222)
  222    FORMAT( ' Enter vertical component of velocity in mm/yr'/
     1           ' with up being positive and down being negative'/
     1           ' or enter 0.0 if unknown.')
	     READ(LUIN, *,err=303,iostat=ios) VUP
         if (ios /= 0) goto 303
	     CALL TOVXYZ( GLAT, GLON, VNORTH, VEAST, VUP, VX, VY, VZ)

      ELSEIF (VOPT .EQ. '2') THEN
	     WRITE(LUOUT, *) ' Enter x-component of velocity in mm/yr.'
	     READ(LUIN, *,err=304,iostat=ios) VX
         if (ios /= 0) goto 304
	     WRITE(LUOUT, *) ' Enter y-component of velocity in mm/yr.'
	     READ(LUIN, *,err=305,iostat=ios) VY
         if (ios /= 0) goto 305
	     WRITE(LUOUT, *) ' Enter z-component of velocity in mm/yr.'
	     READ(LUIN, *,err=306,iostat=ios) VZ
         if (ios /= 0) goto 306
	     CALL TOVNEU(GLAT, GLON, VX, VY, VZ, VNORTH, VEAST, VUP)

      ELSE
	     WRITE(LUOUT, *) 'Improper response -- try again. '
	     GO TO 200

      ENDIF
      RETURN

 300  write (*,'(/)') 
      write (*,*) "Failed to read form of velocity:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

 301  write (*,'(/)') 
      write (*,*) "Failed to read North velocity:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

 302  write (*,'(/)')
      write (*,*) "Failed to read East velocity:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

 303  write (*,'(/)')
      write (*,*) "Failed to read Up velocity:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

 304  write (*,'(/)')
      write (*,*) "Failed to read X velocity:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

 305  write (*,'(/)')
      write (*,*) "Failed to read Y velocity:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

 306  write (*,'(/)')
      write (*,*) "Failed to read Z velocity:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

      END

************************************************************************
      SUBROUTINE COMPSN(YLATT,YLONT,HTT,YLAT,YLON,HT,
     1                  MIN,VN, VE, VU)

*** Compute the position of a point at specified time
*** Upon input of VN, VE, and VU in mm/yr

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (NDLOC = 2195)

      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /TIMREF/ ITREF
      COMMON /QPARM/ STRIKE(NDLOC), HL(NDLOC), EQLAT(NDLOC),
     1          EQLON(NDLOC), SSLIP(NDLOC), DSLIP(NDLOC),
     1          DIP(NDLOC), DEPTH(NDLOC), WIDTH(NDLOC),
     1               EQLATR(50),EQLONR(50),EQRAD(50),
     1               ITEQK(50),NLOC(50),NFP(50),NUMEQ
      COMMON /FILES/ LUIN,LUOUT, I1, I2, I3, I4, I5, I6

** Compute the contribution due to constant velocity
         DTIME = DBLE(MIN - ITREF) / 525960.D0
         CALL RADR8T(YLAT,VN,VE,VNR,VER)
         YLATT = YLAT + VNR*DTIME
         YLONT = YLON - VER*DTIME
         HTT   = HT + ((VU * DTIME) /1000.D0)
       
** Compute the contribution due to earthquakes.
** It is assumed that the components of displacement,
** DNORTH,DWEST,DUP, do not vary from one reference
** frame to another given the accuracy of dislocation
** models.
      DO 10 I = 1, NUMEQ
          IF(ITEQK(I) .GT. ITREF) THEN 
               NTIME = 1
          ELSE
	       NTIME = 0
          ENDIF
	  IF(MIN .LT. ITEQK(I)) NTIME = NTIME - 1
	  IF(NTIME .NE. 0) THEN 
             CALL RADII(EQLATR(I),RADMER,RADPAR)
             DDLAT = (YLAT - EQLATR(I))*RADMER
             DDLON = (YLON - EQLONR(I))*RADPAR
             DIST = DSQRT(DDLAT*DDLAT + DDLON*DDLON)
             IF(DIST .LE. EQRAD(I)) THEN
               ISTART = NLOC(I)
               IEND = NLOC(I) + NFP(I) - 1
               DO 5 JREC = ISTART,IEND
                 CALL DISLOC(YLAT,YLON,STRIKE(JREC),HL(JREC),
     &              EQLAT(JREC),EQLON(JREC),SSLIP(JREC),     
     &              DSLIP(JREC),DIP(JREC),DEPTH(JREC),
     &              WIDTH(JREC),DNORTH,DWEST,DUP)     
                 YLATT = YLATT + NTIME*DNORTH
                 YLONT = YLONT + NTIME*DWEST     
                 HTT   = HTT   + NTIME*DUP
    5          CONTINUE
             ENDIF
          ENDIF  
   10 CONTINUE

*** Compute contribution due to postseismic deformation
      CALL PSDISP(YLAT, YLON, MIN, DNORTH, DEAST, DUP)
c     write(i2,*) 'dnorth, deast, dup = ',dnorth,deast,dup
      CALL RADII (YLAT, RMER, RPAR)
      YLATT = YLATT + DNORTH/RMER
      YLONT = YLONT - DEAST/RPAR
      HTT  = HTT + DUP
      RETURN
      END
***************************************************************
      SUBROUTINE NEWCOR(YLAT,YLON,HTOLD,MIN1,MIN2,
     1         YLAT3,YLON3,HTNEW,DN,DE,DU,VN, VE, VU)

*** Estimate coordinates at time MIN2 given coordinates at time MIN1.
*** Estimate displacements from time MIN1 to time MIN2.
*** Upon input, the velocities--VN, VE, and VU--are in mm/yr.
*** Upon output, the displacements--DN, DE, and DU-- are in meters.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      HT = HTOLD

      CALL COMPSN(YLAT1,YLON1,HT1,YLAT,YLON,HT,
     1            MIN1,VN, VE, VU)

      CALL COMPSN(YLAT2,YLON2,HT2,YLAT,YLON,HT,
     1            MIN2, VN, VE, VU)

      YLAT3 = YLAT + YLAT2 - YLAT1
      YLON3 = YLON + YLON2 - YLON1
      HTNEW   = HT   + HT2   - HT1

      CALL RADII(YLAT,RADMER,RADPAR)

      DN =  RADMER * (YLAT2 - YLAT1)
      DE = -RADPAR * (YLON2 - YLON1)
      DU =  HT2 - HT1

      RETURN
      END
******************************************************************
      subroutine PREDV(ylat, ylon, eht, date, iopt,
     1   jregn, vn, ve, vu) 

** Predict velocity in iopt reference frame       

** ylat       input - north latitude (radians)
** ylon       input - west longitude (radians)
** eht        input - ellipsoid height (meters)
** date       input - date (decimal years)
** iopt       input - reference frame
** jregn      output - deformation region
** vn         output - northward velocity in mm/yr
** ve         output - eastward velocity in mm/yr
** vu         output - upward velocity in mm/yr

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /FILES/ LUIN,LUOUT, I1, I2, I3, I4, I5, I6

** Get reference latitude (RLAT) and reference longitude (RLON)

c*** the following two lines of code were added on July 20, 2015
         elon = -ylon
         call toxyz(ylat,elon,eht,x, y, z)

         IF(IOPT .EQ. 16) THEN
             RLAT = YLAT
             RLON = YLON
         ELSE
             CALL XtoITRF2014(X,Y,Z,RLAT,RLON,EHT2014,DATE,IOPT)
         ENDIF

** Get deformation region

         CALL GETREG(RLAT,RLON,JREGN)
	 IF (JREGN .EQ. 0) THEN
	    VN = 0.D0
	    VE = 0.D0
	    VU = 0.D0
	    RETURN
         ENDIF
	  CALL COMVEL( RLAT, RLON, JREGN, VN, VE, VU,SN,SE,SU)

c*** Convert  velocity to reference of iopt, if iopt != ITRF2014
	  IF (IOPT .NE. 16) THEN
	    CALL TOVXYZ( YLAT, ELON, VN, VE, VU, VX, VY, VZ)
        CALL VTRANF( X, Y, Z, VX, VY, VZ, 16, IOPT)
        CALL TOVNEU( YLAT, ELON, VX, VY, VZ, VN, VE, VU)
      ENDIF
      RETURN
      END

****************************************************************
      subroutine TRFVEL

*** Transform velocities from one reference frame to another

      implicit double precision (a-h, o-z)
      implicit integer*4 (i-n)
      parameter (numref = 17)
      character*80 nameif,name24
      character*30 namef
      character*24 frame1, frame2
      character*1 option
      character*1 vopt, LATDIR, LONDIR
      character record*120
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
      COMMON /CONST/ A, F, E2, EPS, AF, PI, TWOPI, RHOSEC

      write( luout, 100)
  100 format(
     1  ' Please enter the name of the file to contain '/
     1  ' the transformed velocities. ')
      read( luin, '(a30)',err=600,iostat=ios) namef
      if (ios /= 0) goto 600
      
      open( i2, file = namef, status = 'unknown')
      CALL HEADER

  105 write( luout, 110)
  110 format( /'*******************************'/
     1   ' Enter the reference frame of the input velocities.')
      call MENU1(iopt1, frame1)
      if (iopt1 .lt. 1 .or. iopt1 .gt. numref) then
	 write( luout, *) ' Improper selection -- try again.'
	 go to 105
      endif

  115 write( luout, 120)
  120 format( /' Enter the reference frame for the output velocities.')
      call MENU1(iopt2, frame2)
      if (iopt2 .lt. 1 .or. iopt2 .gt. numref) then
	 write( luout, *) 'Improper selection -- try again.'
	 go to 115
      endif

      write( i2, 125) frame1, frame2
  125 format( ' TRANSFORMING VELOCITIES FROM ', A24, ' TO ', a24//)

  130 write( luout, 140)
  140 format( /'**********************************'/
     1  ' Velocities will be transformed at each specified point.'/
     1  ' Please indicate how you wish to input points.'/
     1  '    0...No more points.  Return to main menu.'/
     1  '    1...Individual points entered interactively.'/
     1  '    2...Transform velocities contained in batch file '/
     1  '        of delimited records of the form: '/
     1  '        LAT,LON,VN,VE,VU,TEXT ' /
     1  '        LAT = latitude in degrees (positive north/DBL PREC)'/
     1  '        LON = longitude in degrees (positive west/DBL PREC)'/
     1  '        VN = northward velocity in mm/yr (DBL PREC) '/
     1  '        VE = eastwars velocity in mm/yr (DBL PREC) '/
     1  '        VU = upward velocity in mm/yr (DBL PREC) '/
     1  '        TEXT = descriptive text (CHARACTER*24) '/
     1  '        Example: '/
     1  '        40.731671553,112.212671753, 3.7,3.8,-2.4,SALT AIR '/)
      read( luin,'(a1)',err=601,iostat=ios) option
      if (ios /= 0) goto 601
      if (option .eq. '0') then
	     go to 500

      elseif (option .eq. '1') then
         write(i2,141)
  141    format(16x,' INPUT VELOCITIES       OUTPUT VELOCIITES'/)
	     call GETPNT( latd, latm, slat, LATDIR, lond, lonm, slon,
     1       LONDIR, name24, x, y, z, ylat, ylon, eht)
          elon = - ylon
	     call GETVLY( ylat, elon, vx, vy, vz, vn, ve, vu, vopt, 211)
	     vx1 = vx
	     vy1 = vy
	     vz1 = vz
         
         call VTRANF( x, y, z, vx1, vy1, vz1, iopt1, iopt2)
         
	     call TOVNEU( ylat, elon, vx1, vy1, vz1, vn1, ve1, vu1)
         xlat = (ylat*rhosec)/3600.d0
         xlon = (ylon*rhosec)/3600.d0
	     call PRNTVL( vn, ve, vu, vx, vy, vz, vn1, ve1, vu1,
     1                vx1, vy1, vz1, name24, 1,xlat, xlon)

	     go to 130
      elseif (option .eq. '2') then
         eht = 0.d0
         write (luout, 200)
  200    format(/' Enter name of input file: ')
         read(luin, '(a)',err=602,iostat=ios) nameif
         if (ios /= 0) goto 602
         write(i2,205)
  205    format(/,3x,'INPUT',10x,'INPUT',10x,'NORTHWARD',
     &      ' EASTWARD','    UPWARD  NAME'/
     &      '   LATITUDE       LONGITUDE      VELOCITY',2x,
     &       'VELOCITY  VELOCITY',/)

         open (i1, file = nameif, status = 'old')

  210    read(i1,'(a)',end = 220,err=603,iostat=ios) record     
         if (ios /= 0) goto 603
         call interprate_latlonvel_record(record,glat,glon,
     &                      vn,ve,vu,name24)
         ylat = (glat*3600.d0) / rhosec
         ylon = (glon*3600.d0) / rhosec
         elon = -ylon
         call TOXYZ (ylat, elon, eht, x, y, z)
         call TOVXYZ (ylat, elon, vn,ve,vu,vx,vy,vz)
         vx1 = vx
         vy1 = vy
         vz1 = vz
         
         call VTRANF ( x, y, z, vx1, vy1, vz1, iopt1, iopt2)
         
         call TOVNEU (ylat, elon, vx1, vy1, vz1, vn1, ve1, vu1)

         write(i2,215) glat,glon,vn1,ve1,vu1,name24
  215    format(1x,f14.9,2x,f14.9,2x,f8.2,2x,f8.2,2x,f8.2,2x,a24)
         go to 210
  220    close (i1, status = 'keep')
      endif
  500 continue
      close (i2, status = 'keep')
      write(luout,501)
  501 format(/,' Velocities have been transformed.',/)
      return

  600 write (*,'(/)')
      write (*,*) "Failed to read file name in TRFVEL:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  601 write (*,'(/)')
      write (*,*) "Failed to read option in TRFVEL:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  602 write (*,'(/)')
      write (*,*) "Failed to read input file name in TRFVEL:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

  603 write (*,'(/)')
      write(*,*)"Failed to read input file format in TRFVEL:ios=",ios
      write (*,*) "ABNORMAL TERMINATION"
      write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
      stop

      end
*******************************************************************
      SUBROUTINE VTRANF(X,Y,Z,VX,VY,VZ, IOPT1, IOPT2)

*** Convert velocity from reference frame of IOPT1 to 
*** reference frame of IOPT2.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (numref = 17)
      common /tranpa/ tx(numref), ty(numref), tz(numref),
     &                dtx(numref), dty(numref), dtz(numref),
     &                rx(numref), ry(numref), rz(numref),
     &                drx(numref), dry(numref), drz(numref),
     &                scale(numref), dscale(numref), refepc(numref)
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6

      IF(IOPT1 .le. numref .and. IOPT2. le. numref
     &   .and. IOPT1 .gt. 0 .and. IOPT2 .gt. 0 ) THEN

*** Convert from mm/yr to m/yr
         VX = VX /1000.d0
         VY = VY / 1000.d0
         VZ = VZ / 1000.d0

*** From IOPT1 to ITRF2014 
*** (following equations use approximations assuming
*** that rotations and scale change are small)
         WX = -drx(iopt1)           
         WY = -dry(iopt1)      
         WZ = -drz(iopt1)      
	     DS = -dscale(iopt1)
         VX = VX - dtx(iopt1) + DS*X + WZ*Y - WY*Z
         VY = VY - dty(iopt1) - WZ*X  +DS*Y + WX*Z
         VZ = VZ - dtz(iopt1) + WY*X - WX*Y + DS*Z

*** From ITRF2014 to IOPT2 reference frame
*** (following equations use approximations assuming
***  that rotations and scale change are small)
         WX = drx(iopt2)
         WY = dry(iopt2)
         WZ = drz(iopt2)
	     DS = dscale(iopt2)
         VX = VX + dtx(iopt2) + DS*X + WZ*Y - WY*Z
         VY = VY + dty(iopt2) - WZ*X + DS*Y + WX*Z 
         VZ = VZ + dtz(iopt2) + WY*X - WX*Y + DS*Z

*** FROM m/yr to mm/yr
         VX = VX * 1000.d0
         VY = VY * 1000.d0
         VZ = VZ * 1000.d0

      ELSE
         write(luout,*) ' Improper reference frame in routine vtranf'
         stop
      ENDIF

      RETURN
      END

******************************************************
      subroutine PRNTVL(VN, VE, VU, VX, VY, VZ, VN1, VE1, VU1,
     1                  VX1, VY1, VZ1, NAME24, IPRINT,XLAT,XLON)

*** Print transformed velocities

      implicit double precision (a-h, o-z)
      implicit integer*4 (i-n)
      character*80 name24
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6

      if (iprint .eq. 1) then
      write( luout, 100) vn1, ve1, vu1, vx1, vy1, vz1
  100 format( ' ****************************************'/
     1  ' New northward velocity = ', f8.2, ' mm/yr' /
     1  ' New eastward velocity  = ', f8.2, ' mm/yr'/
     1  ' New upward velocity    = ', f8.2, ' mm/yr'/
     1  ' New x velocity         = ', f8.2, ' mm/yr'/
     1  ' New y velocity         = ', f8.2, ' mm/yr'/
     1  ' New z velocity         = ', f8.2, ' mm/yr'/)
      endif 

      write( i2, 200) name24, xlat, xlon,
     1    vn, vn1, ve, ve1, vu, vu1,
     1    vx, vx1, vy, vy1, vz, vz1
  200 format(1x, a24 /
     1   1x, 'latitude = ',F14.9,2x,'longitude = ',F14.9 /
     1   5x, 'northward velocity ', f8.2, 6x, f8.2, ' mm/yr' /
     1   5x, 'eastward velocity  ', f8.2, 6x, f8.2, ' mm/yr' /
     1   5x, 'upward velocity    ', f8.2, 6x, f8.2, ' mm/yr' /
     1   5x, 'x velocity         ', f8.2, 6x, f8.2, ' mm/yr' /
     1   5x, 'y velocity         ', f8.2, 6x, f8.2, ' mm/yr' /
     1   5x, 'z velocity         ', f8.2, 6x, f8.2, ' mm/yr' /)

      return 
      end
*******************************************************************
       SUBROUTINE HEADER

       IMPLICIT DOUBLE PRECISION (A-H, O-Z)
       IMPLICIT INTEGER*4 (I-N)
       character  Trans4D_version*10
       COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6
       COMMON /VERSION/ Trans4D_version                        

       WRITE(I2, 10) Trans4D_version
   10  FORMAT(' Trans4D (VERSION ',a,') OUTPUT' / )
       RETURN
       END
*********************************************
      SUBROUTINE GETMDY(MONTH, IDAY, IYEAR, DATE, MINS, TEST)

*** Read month-day-year and convert to decimal years
*** and Julian time in minutes      
***    MONTH      output - number from 1 to 12
***    IDAY       output - number from 1 to 31
***    IYEAR      output - must be after 1906
***    DATE       output - corresponding time in decimal years
***    MINS       output - corresponding julian time in minutes
***    TEST       output - if (true) then there is an error

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION M(12)
      LOGICAL TEST
      CHARACTER*1 TOPT, answer
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6

      M(1) = 31
      M(2) = 28
      M(3) = 31
      M(4) = 30
      M(5) = 31
      M(6) = 30
      M(7) = 31
      M(8) = 31
      M(9) = 30
      M(10) = 31
      M(11) = 30
      M(12) = 31

   3  WRITE (LUOUT,1)
   1  FORMAT (' How do you wish to enter the time?'/
     1        '     1. month-day-year in free format'/
     1        '        for example, 5,12,1979'/
     1        '        represents May 12, 1979'/
     1        '     2. decimal year'/
     1        '        for example the entry 1979.359'/
     1        '        represents UTC midnight at the begining'/
     1        '        of May 12, 1979.'/) 

      READ (LUIN,5,err=100,iostat=ios) TOPT
      if (ios /= 0) goto 100
   5  format (A1)

      IF (TOPT .EQ. '1') Then
  121     write (luout,*) ' Enter month-day-year '
          READ(LUIN,*,err=101,iostat=ios) MONTH,IDAY,IYEAR
          if (ios /= 0)  goto 101
    
          IF(IYEAR .le. 1906) THEN
            WRITE(LUOUT,10)
   10       FORMAT(' The model is not valid for dates prior ',
     1           'to 1906.'/)
            write(luout,*) ' Do you wish to re-enter the date? (y/n)'
            read(luin, '(A1)',err=100,iostat=ios) ANSWER
            if (ios /= 0) goto 100
            if(ANSWER .eq. 'y' .or. ANSWER .eq. 'Y')GO TO 121
            go to 100
          ENDIF

          IF(MONTH .le. 0 .or. MONTH .gt. 12) THEN
            WRITE(LUOUT,20)
   20       FORMAT(' Improper month specified.'/)
            write(luout,*) ' Do you wish to re-enter the date? (y/n)'
            read(luin, '(A1)',err=100,iostat=ios) ANSWER
            if (ios /= 0) goto 100
            if(ANSWER .eq. 'y' .or. ANSWER .eq. 'Y')GO TO 121
            go to 100
          ENDIF

          IF(IDAY .le. 0 .or. IDAY .gt. 31) THEN
            WRITE(LUOUT,30)
   30       FORMAT(' Improper day specified.'/)
            write(luout,*) ' Do you wish to re-enter the date? (y/n)'
            read(luin, '(A1)',err=100,iostat=ios) ANSWER
            if (ios /= 0) goto 100
            if(ANSWER .eq. 'y' .or. ANSWER .eq. 'Y')GO TO 121
            go to 100
          ENDIF

          CALL IYMDMJ(IYEAR, MONTH, IDAY, MJD)
          CALL IYMDMJ(IYEAR, 1, 1, MJD0)
          IYEAR1 = IYEAR + 1
          CALL IYMDMJ(IYEAR1, 1, 1, MJD1)
          DAY = DBLE(MJD - MJD0)
          DENOM = DBLE(MJD1 - MJD0)
          DATE = DBLE(IYEAR) + (DAY / DENOM)
          MINS = MJD * 24 * 60
          TEST = .FALSE.
          RETURN

      ELSEIF (TOPT .EQ. '2') then
  122     write(luout,*) ' Enter decimal year '
          READ (LUIN, *,err=102,iostat=ios) DATE
          if (ios /= 0)  goto 102
          
          IF (DATE .lt. 1906.0d0) then
             write (luout, 10)
             write(luout,*) ' Do you wish to re-enter the date? (y/n)'
             read(luin, '(A1)',err=100,iostat=ios) ANSWER
             if (ios /= 0) goto 100
             if(ANSWER .eq. 'y' .or. ANSWER .eq. 'Y')GO TO 122
             go to 100
          ENDIF
**** add small increment to circumvent round-off error
*         DATE = DATE + 0.0004D0
****
          IYEAR = DATE
          CALL IYMDMJ(IYEAR, 1, 1, MJD0)
          IYEAR1 = IYEAR + 1
          CALL IYMDMJ(IYEAR1, 1, 1, MJD1)
          LEAP = 0
          IF ((MJD1 - MJD0) .eq. 366) LEAP = 1
          REMDAY = (DATE - IYEAR)* (MJD1 - MJD0)
          IBEGIN = 0
          ITOTAL = 31
          IF (REMDAY .LT. ITOTAL) then
              MONTH = 1
              IDAY = REMDAY - IBEGIN + 1
              CALL IYMDMJ(IYEAR,MONTH,IDAY,MJD)
              MINS = MJD * 24 * 60
              TEST = .FALSE.
              RETURN
          ENDIF
          IBEGIN = ITOTAL
          ITOTAL = ITOTAL + LEAP
          DO I = 2, 12
             ITOTAL = ITOTAL + M(I)
             IF (REMDAY .LT. ITOTAL) then
                  MONTH = I
                  IDAY = REMDAY - IBEGIN + 1
                  TEST = .FALSE.
                  CALL IYMDMJ(IYEAR,MONTH,IDAY,MJD)
                  MINS = MJD * 24 * 60
                  RETURN
             ENDIF
             IBEGIN = ITOTAL
           ENDDO
           Write(LUOUT, 60)
   60      Format (' Error could not convert Decimal years to '
     1             ,'month-day-year')
           go to 100
       ELSE
         write (luout, 70)
   70    format (' Improper entry')
         Go TO 3
       ENDIF

 100   write (*,'(/)')
       write (*,*) "Failed to read option in GETDMY:ios=",ios
       write (*,*) "ABNORMAL TERMINATION"
       write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
       stop

 101   write (*,'(/)')
       write (*,*) "Wrong MDY input in GETDMY:ios=",ios
       write (*,*) "ABNORMAL TERMINATION"
       write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
       stop

 102   write (*,'(/)')
       write (*,*) "Wrong decimal year in GETDMY:ios=",ios
       write (*,*) "ABNORMAL TERMINATION"
       write (*,*) "PLEASE CHECK YOUR INPUT FILE AND TRY AGAIN"
       stop

       END
************************************
      SUBROUTINE IYMDMJ( IYR, IMON, IDAY, MJD )
C
C********1*********2*********3*********4*********5*********6*********7**
C
C NAME:       IYMDMJ
C VERSION:    Sep. 17, 2010
C WRITTEN BY: R. SNAY (after M. SCHENEWERK)
C PURPOSE:    CONVERT DATE TO MODIFIED JULIAN DATE 
C
C INPUT PARAMETERS FROM THE ARGUEMENT LIST:
C -----------------------------------------
C IDAY              DAY
C IMON              MONTH
C IYR               YEAR
C
C OUTPUT PARAMETERS FROM ARGUEMENT LIST:
C --------------------------------------
C MJD               MODIFIED JULIAN DATE 
C
C
C LOCAL VARIABLES AND CONSTANTS:
C ------------------------------
C A                 TEMPORARY STORAGE
C B                 TEMPORARY STORAGE
C C                 TEMPORARY STORAGE
C D                 TEMPORARY STORAGE
C IMOP              TEMPORARY STORAGE
C IYRP              TEMPORARY STORAGE
C
C GLOBAL VARIABLES AND CONSTANTS:
C ------------------------------
C
C
C       THIS MODULE CALLED BY: GENERAL USE
C
C       THIS MODULE CALLS:     DINT
C
C       INCLUDE FILES USED:
C
C       COMMON BLOCKS USED:       
C
C       REFERENCES:            DUFFETT-SMITH, PETER  1982, 'PRACTICAL
C                              ASTRONOMY WITH YOUR CALCULATOR', 2ND
C                              EDITION, CAMBRIDGE UNIVERSITY PRESS,
C                              NEW YORK, P.9
C
C       COMMENTS:              THIS SUBROUTINE REQUIRES THE FULL YEAR,
C                              I.E. 1992 RATHER THAN 92.  
C
C********1*********2*********3*********4*********5*********6*********7**
C::LAST MODIFICATION
C::8909.06, MSS, DOC STANDARD IMPLIMENTED
C::9004.17, MSS, CHANGE ORDER YY MM DD
C********1*********2*********3*********4*********5*********6*********7**
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER*4 (I-N)
C
      INTEGER*4     A, B, C, D

      IYRP = IYR
C
C........  0.0  EXPLICIT INITIALIZATION
C
      IF( IMON .LT. 3 ) THEN
        IYRP= IYRP - 1
        IMOP= IMON + 12
      ELSE
        IMOP= IMON
      END IF
C
C........  1.0  CALCULATION
C
      A=  IYRP*0.01D0
      B=  2 - A + DINT( A*0.25D0 )
      C=  365.25D0*IYRP
      D=  30.6001D0*(IMOP + 1)
      MJD =  (B + C + D + IDAY - 679006) 
C      
      RETURN
      END
*****************************************************
      SUBROUTINE PSDISP(YLAT, YLON, MIN, DNORTH, DEAST, DUP)
********
*   Compute total postseismic displacement for all earthquakes
*
* INPUT
*   YLAT       latitude of point in radians, positive north
*   YLON       longitude of point in radians, positive west
*   MIN        modified julian date of reference epoch for new coordinates
*               in minutes
*
*   DNORTH     Total northward postseismic displacement at point during
*              period from ITREF to MIN in meters
*   DEAST      TOTAL eastward postseismic displacement
*   DUP        Total upward postseismic displacement 
*******

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER*4 (I-N)
      LOGICAL INSIDE

      parameter (NUMPSG = 1)
      COMMON /CONST/ A, F,E2,EPS,AF,PI,TWOPI,RHOSEC
      COMMON /TIMREF/ ITREF
      COMMON /PSGRID/ PSGLX(NUMPSG), PSGUX(NUMPSG),
     1          PSGLY(NUMPSG), PSGUY(NUMPSG),
     1          ICNTPX(NUMPSG), ICNTPY(NUMPSG), NBASEP(NUMPSG)
      COMMON /PGRID/ PS(18000)
      DIMENSION ITEQ(NUMPSG)      
      DIMENSION TAU(NUMPSG)   
      DIMENSION WEI(2,2)
      DIMENSION AMP(2,2,3)

*** Relaxation constant (in years) for 2002 Denali earthquake
      TAU(1) = 5.0D0

*** Modofied Julian Date (in minutes) for the 2002 Denali earthquake
      IYEAR = 2002
      IMO = 11
      IDAY = 3
      CALL IYMDMJ(IYEAR,IMO,IDAY, MJD)
      ITEQ(1) = MJD*60*24

      DNORTH = 0.0D0
      DEAST = 0.0D0
      DUP = 0.0D0

      DO K = 1, NUMPSG
*** Check if the point is inside the grid
         POSX = YLON*180.d0/PI
         POSX = 360.d0 - POSX
         IF (POSX .GT. 360.D0) POSX = POSX - 360.D0
         POSY = YLAT*180.D0/PI
         CALL GRDCHK(POSX, POSY, PSGLX(K), PSGUX(K),
     1            PSGLY(K), PSGUY(K), INSIDE)

         IF (INSIDE ) THEN
*** Get the indices for the lower left-hand corner of the grid
         CALL PSGWEI(POSX,POSY,K,I,J,WEI)
*** Get the displacement amplitude at the four corners
         CALL GRDAMP(K,I,J,AMP,PS)

         ANORTH = WEI(1,1)*AMP(1,1,1) + WEI(1,2)*AMP(1,2,1)
     1          + WEI(2,1)*AMP(2,1,1) + WEI(2,2)*AMP(2,2,1)     
         AEAST  = WEI(1,1)*AMP(1,1,2) + WEI(1,2)*AMP(1,2,2)
     1          + WEI(2,1)*AMP(2,1,2) + WEI(2,2)*AMP(2,2,2)
         AUP    = WEI(1,1)*AMP(1,1,3) + WEI(1,2)*AMP(1,2,3)
     1          + WEI(2,1)*AMP(2,1,3) + WEI(2,2)*AMP(2,2,3)

c        write (6, 30) anorth, aeast, aup
c  30    format (1x, 3f15.5)

*** Convert amplitudes from mm to meters
         ANORTH = ANORTH / 1000.D0
         AEAST =  AEAST  / 1000.D0
         AUP =    AUP    / 1000.D0

         IF (MIN .GT. ITEQ(K)) THEN
            DTIME = DBLE(MIN - ITEQ(K))/(60.D0*24.D0*365.D0)
            FACTOR = 1.D0 - DEXP(-DTIME/TAU(K))
            DNORTH = DNORTH + ANORTH*FACTOR
            DEAST = DEAST + AEAST*FACTOR
            DUP = DUP + AUP*FACTOR
         ENDIF
         IF (ITREF .GT. ITEQ(K)) THEN
            DTIME = DBLE(ITREF - ITEQ(K))/(60.D0*24.D0*365.D0)
            FACTOR = 1.D0 - DEXP(-DTIME/TAU(K))
            DNORTH = DNORTH - ANORTH*FACTOR
            DEAST = DEAST - AEAST*FACTOR
            DUP = DUP - AUP*FACTOR
         ENDIF
         ENDIF
      ENDDO
      RETURN
      END

****************************************************
      SUBROUTINE GRDCHK (POSX, POSY, GRDLX, GRDUX,
     1             GRDLY, GRDUY, INSIDE)

C
C ROUTINE CHECKS IF THE POINT HAVING COORDINATES (POSX, POSY)
C IS WITHIN THE REGION SPANNED BY THE GRID
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      LOGICAL INSIDE

      INSIDE = .TRUE.

      IF (POSX .LT. GRDLX .OR. POSX .GT. GRDUX) THEN
         INSIDE = .FALSE.
      ENDIF
      IF (POSY .LT. GRDLY .OR. POSY .GT. GRDUY) THEN
         INSIDE = .FALSE.
      ENDIF

      RETURN
      END
*******************************************************************
      SUBROUTINE PSGWEI (POSX, POSY, K, I, J, WEI)

C
C********1*********2*********3*********4*********5*********6*********7**
C
C PURPOSE:     THIS SUBROUTINE RETURNS THE INDICES OF THE LOWER-LEFT
C              HAND CORNER OF THE GRID CELL CONTAINING THE POINT
C              AND COMPUTES NORMALIZED WEIGHTS FOR 
C              BI-LINEAR INTERPOLATION OVER A PLANE
C              
C  INPUT PARAMETERS FROM ARGUMENT LIST:
C  ------------------------------------
C POSX         LONGITUDE OF POINT IN DEGREES, POSITIVE EAST
C POSY         LATITUDE OF POINT IN DEGREES, POSITIVE NORTH
C K            ID OF EARTHQUAKE GRID                   
C
C  OUTPUT PARAMETERS FROM ARGUMENT LIST:
C  -------------------------------------
C I, J         THE COORDINATES OF LOWER LEFT CORNER OF THE GRID
C              CONTAINING THE ABOVE POSITION
C WEI          A TWO BY TWO ARRAY CONTAINING THE NORMALIZED WEIGHTS
C              FOR THE CORNER VECTORS
C
C  GLOBAL VARIABLES AND CONSTANTS:
C  -------------------------------
C NONE
C
C    THIS MODULE CALLED BY:   PSDISP
C
C    THIS MODULE CALLS:       NONE
C
C    INCLUDE FILES USED:      NONE
C
C    COMMON BLOCKS USED:      /PSGRID/, /CONST/
C
C    REFERENCES:  SEE RICHARD SNAY
C
C    COMMENTS:
C
C********1*********2*********3*********4*********5*********6*********7**
C    MOFICATION HISTORY:
C::9302.11, CRP, ORIGINAL CREATION FOR DYNAP
C::9511.09, RAS, MODIFIED FOR HTDP
C::9712.05, RAS, MODIFIED TO ACCOUNT FOR MULTIPLE GRIDS
C********1*********2*********3*********4*********5*********6*********7**
    
C**** COMPUTES THE WEIGHTS FOR AN ELEMENT IN A GRID

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (NUMPSG = 1)
      DIMENSION WEI(2,2)
      COMMON /PSGRID/ PSGLX(NUMPSG), PSGUX(NUMPSG), 
     1          PSGLY(NUMPSG), PSGUY(NUMPSG),
     1          ICNTPX(NUMPSG), ICNTPY(NUMPSG), NBASEP(NUMPSG)
      COMMON /CONST/ A,F,E2,EPS,AF,PI,TWOPI,RHOSEC

C*** Obtain indices for the lower-left corner of the cell
C*** containing the point
      STEPX = (PSGUX(K) - PSGLX(K)) / ICNTPX(K)
      STEPY = (PSGUY(K) - PSGLY(K)) / ICNTPY(K)
      I = IDINT((POSX - PSGLX(K))/STEPX) + 1
      J = IDINT((POSY - PSGLY(K))/STEPY) + 1
c     write(6,1001) K, I, J
c1001 format(1x, 'quake = ', I5 /
c    1       1x, ' i = ', I5 /
c    1       1x, ' j = ', I5)

C*** Compute the limits of the grid cell 
      GRLX = PSGLX(K) + (I - 1) * STEPX
      GRUX = GRLX + STEPX                    
      GRLY = PSGLY(K) + (J - 1) * STEPY                
      GRUY = GRLY + STEPY                     

C*** Compute the normalized weights for the point               
      DENOM = (GRUX - GRLX) * (GRUY - GRLY)
      WEI(1,1) = (GRUX - POSX) * (GRUY - POSY) / DENOM
      WEI(2,1) = (POSX - GRLX) * (GRUY - POSY) / DENOM
      WEI(1,2) = (GRUX - POSX) * (POSY - GRLY) / DENOM
      WEI(2,2) = (POSX - GRLX) * (POSY - GRLY) / DENOM

      RETURN
      END

C*********************************************************************
      SUBROUTINE GRDAMP (K, I, J, AMP, PS)
C********1*********2*********3*********4*********5*********6*********7**
C
C PURPOSE:     THIS SUBROUTINE RETRIEVES THE AMPLITUDES OF THE FOUR
C              GRID NODES OFGRID K WHERE I,J ARE THE INDICES OF
C              THE LOWER LEFT HAND CORNER
C              
C  INPUT PARAMETERS FROM ARGUMENT LIST:
C  ------------------------------------
C
C K            ID OF EARTHQUAKE CORRESPONDING TO GRID
C I, J         THE COORDINATES OF LOWER LEFT CORNER OF THE GRID
C              CONTAINING THE ABOVE POSITION
C PS           THE ARRAY CONTAINING ALL THE GRIDDED AMPLITUDES
C
C  OUTPUT PARAMETERS FROM ARGUMENT LIST:
C  -------------------------------------
C AMP          A TWO BY TWO ARRAY CONTAINING THE 3D AMPLITUDES
C              FOR THE CORNERS OF THE GRID
C
C  GLOBAL VARIABLES AND CONSTANTS:
C  -------------------------------
C NONE
C
C    THIS MODULE CALLED BY:   PSDISP
C
C    THIS MODULE CALLS:       NONE
C
C    INCLUDE FILES USED:      NONE
C
C    COMMON BLOCKS USED:      NONE     
C
C    REFERENCES:  SEE RICHARD SNAY
C
C    COMMENTS:
C
C********1*********2*********3*********4*********5*********6*********7**
C    MOFICATION HISTORY:
C::2011.08.17, RAS, ORIGINAL CREATION FOR TRANS4D
C********1*********2*********3*********4*********5*********6*********7**


      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION AMP(2,2,3), PS(*)

      DO 30 II = 0,1
         DO 20 IJ = 0,1
            DO 10 IVEC = 1, 3
               INDEX = IPSGRD(K, I + II, J + IJ, IVEC)
               AMP(II + 1, IJ + 1, IVEC) = PS(INDEX)
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE   

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
C---------------------------------------------------------------------
      subroutine extract_name (name,i)

C  In f90, use trim(name) to truncate all spaces after the name.
C  But "trim" were not available in f77. This is what this subroutine is for

      implicit none
      integer*4 i
      character name*80,scratch*80

      do i=80,1,-1
        if (name(i:i) /= ' ') exit 
      enddo
      scratch = name(1:i)
      name    = scratch

      return
      end
C-----------------------------------------------------------------------------------
         subroutine interprate_XYZ_record (record,x,y,z,name)

         implicit none

         integer*4   i,j,length
         real*8      x,y,z
         character   name*80,record*120,record1*120,chars*120
         character   xxxx*80,yyyy*80,zzzz*80

         record1 = trim(adjustl(record))
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == '	') then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 

C  Get x or phi

         xxxx = '0'
         do i=1,length
           if (chars(i:i) == '1') then
             xxxx(i:i) = record1(i:i)
             record1(i:i) = ' '
           elseif (chars(i:i) == '0') then
             exit
           endif
         enddo
         record = trim(adjustl(record1))
         read (xxxx,*) x
         if (record(1:1) == ',') then
           record(1:1) = ' '
         endif
         record1 = adjustl(record)
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == "	") then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 
           
C  Get y or lamda

         yyyy = '0'
         do i=1,length
           if (chars(i:i) == '1') then
             yyyy(i:i) = record1(i:i)
             record1(i:i) = ' '
           elseif (chars(i:i) == '0') then
             exit
           endif
         enddo
         record = trim(adjustl(record1))
         read (yyyy,*) y
         if (record(1:1) == ',') then
           record(1:1) = ' '
         endif
         record1 = adjustl(record)
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == "	") then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 

C  Get z or h     

         zzzz = '0'
         do i=1,length
           if (chars(i:i) == '1') then
             zzzz(i:i) = record1(i:i)
             record1(i:i) = ' '
           elseif (chars(i:i) == '0') then
             exit
           endif
         enddo
         record = trim(adjustl(record1))
         read (zzzz,*) z
         if (record(1:1) == ',') then
           record(1:1) = ' '
         endif
         record1 = adjustl(record)
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == "	") then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 

C  Get name

         name = trim(record1)

C  Done

         return 
         end
C-----------------------------------------------------------------------------------
      subroutine interprate_latlon_record (record,x,y,name)

         implicit none

         integer*4   i,j,length
         real*8      x,y,z
         character   record*120,record1*120,chars*120
         character   xxxx*80,yyyy*80,zzzz*80
         character   name*50
C        character   name*24

         record1 = trim(adjustl(record))
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == '	') then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 

C  Get x or phi

         xxxx = '0'
         do i=1,length
           if (chars(i:i) == '1') then
             xxxx(i:i) = record1(i:i)
             record1(i:i) = ' '
           elseif (chars(i:i) == '0') then
             exit
           endif
         enddo
         record = trim(adjustl(record1))
         read (xxxx,*) x
         if (record(1:1) == ',') then
           record(1:1) = ' '
         endif
         record1 = adjustl(record)
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == "	") then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 
           
C  Get y or lamda

         yyyy = '0'
         do i=1,length
           if (chars(i:i) == '1') then
             yyyy(i:i) = record1(i:i)
             record1(i:i) = ' '
           elseif (chars(i:i) == '0') then
             exit
           endif
         enddo
         record = trim(adjustl(record1))
         read (yyyy,*) y
         if (record(1:1) == ',') then
           record(1:1) = ' '
         endif
         record1 = adjustl(record)
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == "	") then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 

C  Get name

         name = trim(record1)

C  Done

         return 
         end

C********************

      subroutine interprate_latlonvel_record (record,glat,glon,
     *                        vn,ve,vu,name)

         implicit none

         integer*4   i,j,length
         real*8      glat,glon,vn,ve,vu
         character   name*80,record*120,record1*120,chars*120
         character   clat*80,clon*80,cvn*80,cve*80,cvu*80

         record1 = trim(adjustl(record))
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == '	') then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 

C  Get glat (latitude)

         clat = '0'
         do i=1,length
           if (chars(i:i) == '1') then
             clat(i:i) = record1(i:i)
             record1(i:i) = ' '
           elseif (chars(i:i) == '0') then
             exit
           endif
         enddo
         record = trim(adjustl(record1))
         read (clat,*) glat
         if (record(1:1) == ',') then
           record(1:1) = ' '
         endif
         record1 = adjustl(record)
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == "	") then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 
           
C  Get glon (longitude)

         clon = '0'
         do i=1,length
           if (chars(i:i) == '1') then
             clon(i:i) = record1(i:i)
             record1(i:i) = ' '
           elseif (chars(i:i) == '0') then
             exit
           endif
         enddo
         record = trim(adjustl(record1))
         read (clon,*) glon
         if (record(1:1) == ',') then
           record(1:1) = ' '
         endif
         record1 = adjustl(record)
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == "	") then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 

C  Get vn (northward velocity)     

         cvn = '0'
         do i=1,length
           if (chars(i:i) == '1') then
             cvn(i:i) = record1(i:i)
             record1(i:i) = ' '
           elseif (chars(i:i) == '0') then
             exit
           endif
         enddo
         record = trim(adjustl(record1))
         read (cvn,*) vn
         if (record(1:1) == ',') then
           record(1:1) = ' '
         endif
         record1 = adjustl(record)
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == "	") then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 

C  Get ve (eastward velocity)     

         cve = '0'
         do i=1,length
           if (chars(i:i) == '1') then
             cve(i:i) = record1(i:i)
             record1(i:i) = ' '
           elseif (chars(i:i) == '0') then
             exit
           endif
         enddo
         record = trim(adjustl(record1))
         read (cve,*) ve
         if (record(1:1) == ',') then
           record(1:1) = ' '
         endif
         record1 = adjustl(record)
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == "	") then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 

C  Get vu (upward velocity)     

         cvu = '0'
         do i=1,length
           if (chars(i:i) == '1') then
             cvu(i:i) = record1(i:i)
             record1(i:i) = ' '
           elseif (chars(i:i) == '0') then
             exit
           endif
         enddo
         record = trim(adjustl(record1))
         read (cvu,*) vu
         if (record(1:1) == ',') then
           record(1:1) = ' '
         endif
         record1 = adjustl(record)
         length = len_trim(record1)
         chars = ' '
         do i=1,length
           if (record1(i:i) == ' ' .or. record1(i:i) == ',' .or. 
     &         record1(i:i) == "	") then
             chars(i:i) = '0'
           else
             chars(i:i) = '1'
           endif
         enddo 

C  Get name

         name = trim(record1)

C  Done

         return 
         end

C***********

      subroutine get_frames(iopt1, iopt2, frame1, frame2)

*** get reference frames for input positions and 
*** output positions

      implicit double precision (a-h, o-z)
      implicit integer*4 (i-n)
      character*24 frame1, frame2
      parameter (numref = 17)
      COMMON /FILES/ LUIN, LUOUT, I1, I2, I3, I4, I5, I6

   95 write(luout,100)
  100 format (' *******************************************'/
     1  ' Enter the reference frame of the input positions')
      call MENU1(iopt1, frame1)
      if (iopt1 .lt. 1 .or. iopt1 .gt. numref) then
         write(luout,*) ' Improper selection -- try again.  '
         go to 95
      endif

  105 write(luout,110)
  110 format (/' Enter the reference frame for the output positions')
      call MENU1(iopt2, frame2)
      if (iopt2 .lt. 1 .or. iopt2 .gt. numref) then
         write(luout,*) ' Improper selection -- try again.  '
         go to 105
      endif
  
      return 
      end

****************************************

      subroutine tran_frames(xlatin,xlonin,ehtin,iframein,
     *     xlatout,xlonout,ehtout,iframeout,date)

*** transform geodetic coordinates between frames 
*** at a given time

*** xlatin     = input latitude in radians (positive north)
*** xlonin     = input longitude in radians (positive east)
*** ehtin      = input ellipsoid height in meters
*** iframein   = input reference frame
*** xlatout    = output latitude in radians (positive north)
*** xlonout    = output longitude in radians (positive east)
*** ehtout     = output ellipsoid height in meters
*** iframeout  = output reference frame
*** date       = reference date in decimal years

      implicit double precision (a-h, o-z)
      implicit integer*4 (i-n)
      logical frmxyz
     
      call toxyz(xlatin,xlonin,ehtin,xin,yin,zin)
      
      call to_itrf2014(xin,yin,zin,xtemp,ytemp,ztemp,date,iframein)
      call from_itrf2014(xtemp,ytemp,ztemp,
     1         xout,yout,zout,date,iframeout)
  
      if(.not.frmxyz(xout,yout,zout,xlatout,xlonout,ehtout)) then
        write(6,1)
    1   format(' Error in subroutine tran_frames'/
     *         ' cannot convert x,y,z to lat,lon,eht'/
     *         ' program is aborting')
        stop
      endif

      return
      end
         
******************      
      INTEGER FUNCTION IUNGRD(IREGN, I, J, IVEC)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (NUMGRD = 8)
      COMMON /CDGRID/ GRDLX(NUMGRD), GRDUX(NUMGRD),
     1          GRDLY(NUMGRD), GRDUY(NUMGRD),
     1          ICNTX(NUMGRD), ICNTY(NUMGRD), NBASE(NUMGRD)

      IUNGRD = NBASE(IREGN) +
     1      3 * ((J - 1) * (ICNTX(IREGN) + 1) +  (I - 1)) + IVEC

      RETURN
      END