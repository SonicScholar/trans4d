      PROGRAM MAKEEQ

*** Create BLOCK DATA routine called initeq.f
*** that contains the earthquake information

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*20 QNAME
      CHARACTER*81 CARD

      I1 = 11
      I2 = 12
      PI = 4.D0 * DATAN(1.0D0)
      RHOSEC = (180.D0 * 3600.D0) /PI

      OPEN(I1, FILE = 'QUAKE4.1.txt', STATUS = 'OLD')
      OPEN(I2, FILE = 'initeq.f', STATUS = 'UNKNOWN')

      WRITE(I2, 10)
   10 FORMAT(
     * '      BLOCK DATA INITEQ'/
     * '      IMPLICIT DOUBLE PRECISION (A-H, O-Z)'/
     * '      IMPLICIT INTEGER*4 (I-N)'/
     * '      parameter (NDLOC = 2195)'/
     * '      COMMON /QPARM/ STRIKE(NDLOC), HL(NDLOC), EQLAT(NDLOC),'/
     * '     *   EQLON(NDLOC), SSLIP(NDLOC), DSLIP(NDLOC), DIP(NDLOC),'/
     * '     *   DEPTH(NDLOC), WIDTH(NDLOC), EQLATR(50), EQLONR(50),'/
     * '     *   EQRAD(50), ITEQK(50), NLOC(50), NFP(50), NUMEQ '/)

      NUMEQ = 0
      JREC  = 0

   15 READ(I1, 20, END = 500) CARD
   20 FORMAT(A81)

      IF (CARD(1:1) .EQ. 'H') THEN
         NUMEQ = NUMEQ + 1
         READ (CARD, 30) QNAME, IYEAR, MONTH, IDAY, DEGLAT, DEGLON,
     *      RADIUS, NFP, NREF
   30    FORMAT(5X, A20, 4X, I4, 1X, I2, 1X, I2, 4X, F8.5, 1X,
     *          F9.5, 1X, F5.0, 1X, I3, 1X, I2)
C        CALL TOTIME(IYEAR, MONTH, IDAY, ITEQK)
         CALL IYMDMJ(IYEAR, MONTH, IDAY, MJD)
         ITEQK = MJD * 24 * 60
         EQLATR = DEGLAT * 3600.D0/ RHOSEC
         EQLONR = DEGLON * 3600.D0/ RHOSEC
         EQRAD  = RADIUS * 1000.D0
         DO 40 I = 1, NREF
            READ(I1, 20) CARD
   40    CONTINUE
         NLOC = JREC + 1
         WRITE (I2, 41) NUMEQ, EQLATR, NUMEQ, EQLONR, NUMEQ,
     *      EQRAD, NUMEQ, ITEQK, NUMEQ, NLOC, NUMEQ, NFP
   41    FORMAT(6X,'DATA EQLATR(',I3,')  /', E16.10, '/'/
     *          5X,'*   EQLONR(',I3,')  /', E16.10, '/'/
     *          5X,'*   EQRAD(',I3,')   /', E16.10, '/'/
     *          5X,'*   ITEQK(',I3,')   /', I12, '/'/
     *          5X,'*   NLOC(',I3,')    /', I4, '/'/
     *          5X,'*   NFP(',I3,')     /', I4, '/')
         DO 50 I = 1, NFP
            JREC = JREC + 1
            READ(I1, *) KEQ, KREC, STRIKE, HL, DEGLAT, DEGLON,
     *         SSLIP, DSLIP, DIP, DEPTH, WIDTH
            IF ( KEQ .NE. NUMEQ .OR. 
     *           KREC .NE. (JREC - NLOC + 1)) THEN
               WRITE(6, 45) QNAME, KREC
   45          FORMAT ( ' ERROR in QUAKE file for ',A20,/
     *                  ' record = ', I5)
               STOP
            ENDIF
            STRIKE = STRIKE * 3600.D0 / RHOSEC
            EQLAT  = DEGLAT * 3600.D0 / RHOSEC
            EQLON  = DEGLON * 3600.D0 / RHOSEC
            DIP    = DIP    * 3600.D0 / RHOSEC
            WRITE(I2, 46) JREC, STRIKE, JREC, HL, JREC, EQLAT,
     *        JREC, EQLON, JREC, SSLIP, JREC, DSLIP,
     *        JREC, DIP, JREC, DEPTH, JREC, WIDTH
   46       FORMAT(6X,'DATA STRIKE(',I4,') /', E16.10, '/'/
     *             5X,'*    HL(',I4,') /', E16.10, '/'/
     *             5X,'*    EQLAT(',I4,') /', E16.10, '/'/
     *             5X,'*    EQLON(',I4,') /', E16.10, '/'/
     *             5X,'*    SSLIP(',I4,') /', E16.10, '/'/
     *             5X,'*    DSLIP(',I4,') /', E16.10, '/'/
     *             5X,'*    DIP(',I4,') /', E16.10, '/'/
     *             5X,'*    DEPTH(',I4,') /', E16.10, '/'/
     *             5X,'*    WIDTH(',I4,') /', E16.10, '/')
   50    CONTINUE
      ELSE
         WRITE(6, 60)
   60    FORMAT ('QUAKE file has improper format')
         STOP
      ENDIF
      GO TO 15

  500 CONTINUE
      CLOSE (I1, STATUS = 'KEEP')
      WRITE (I2, 510) NUMEQ
  510 FORMAT (6X,'DATA NUMEQ  /', I4, '/'/
     *        6X,'   '/
     *        6X,'END')
      CLOSE (I2, STATUS = 'KEEP')
      STOP
      END

C********************************************************* 
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
