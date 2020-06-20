      BLOCK DATA INITVL
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter (NUMGRD = 8)
      COMMON /CDGRID/ GRDLX(NUMGRD), GRDUX(NUMGRD),
     *   GRDLY(NUMGRD), GRDUY(NUMGRD),
     *   ICNTX(NUMGRD), ICNTY(NUMGRD), NBASE(NUMGRD)
      COMMON /VGRID/ B(800000)
      COMMON /SGRID/ C(800000)

      DATA GRDLX( 1) / 238.20000000000 /,
     *     GRDUX( 1) / 239.49000000000 /,
     *     GRDLY( 1) /  35.80000000000 /,
     *     GRDUY( 1) /  36.79000000000 /,
     *     ICNTX( 1) /129 /,
     *     ICNTY( 1) / 99 /,
     *     NBASE( 1) /     0 /
      DATA GRDLX( 2) / 235.00000000000 /,
     *     GRDUX( 2) / 253.00000000000 /,
     *     GRDLY( 2) /  31.00000000000 /,
     *     GRDUY( 2) /  49.00000000000 /,
     *     ICNTX( 2) /288 /,
     *     ICNTY( 2) /288 /,
     *     NBASE( 2) / 39000 /
      DATA GRDLX( 3) / 253.00000000000 /,
     *     GRDUX( 3) / 286.50000000000 /,
     *     GRDLY( 3) /  24.00000000000 /,
     *     GRDUY( 3) /  40.00000000000 /,
     *     ICNTX( 3) /536 /,
     *     ICNTY( 3) /256 /,
     *     NBASE( 3) /289563 /
      DATA GRDLX( 4) / 253.00000000000 /,
     *     GRDUX( 4) / 294.00000000000 /,
     *     GRDLY( 4) /  24.00000000000 /,
     *     GRDUY( 4) /  50.00000000000 /,
     *     ICNTX( 4) / 82 /,
     *     ICNTY( 4) / 52 /,
     *     NBASE( 4) /703590 /
      DATA GRDLX( 5) / 190.00000000000 /,
     *     GRDUX( 5) / 230.00000000000 /,
     *     GRDLY( 5) /  53.00000000000 /,
     *     GRDUY( 5) /  73.00000000000 /,
     *     ICNTX( 5) /160 /,
     *     ICNTY( 5) / 80 /,
     *     NBASE( 5) /716787 /
      DATA GRDLX( 6) / 231.00000000000 /,
     *     GRDUX( 6) / 240.00000000000 /,
     *     GRDLY( 6) /  49.00000000000 /,
     *     GRDUY( 6) /  51.00000000000 /,
     *     ICNTX( 6) / 36 /,
     *     ICNTY( 6) /  8 /,
     *     NBASE( 6) /755910 /
      DATA GRDLX( 7) / 230.00000000000 /,
     *     GRDUX( 7) / 308.00000000000 /,
     *     GRDLY( 7) /  42.00000000000 /,
     *     GRDUY( 7) /  78.00000000000 /,
     *     ICNTX( 7) / 52 /,
     *     ICNTY( 7) / 36 /,
     *     NBASE( 7) /756909 /
      DATA GRDLX( 8) / 265.00000000000 /,
     *     GRDUX( 8) / 303.00000000000 /,
     *     GRDLY( 8) /   6.00000000000 /,
     *     GRDUY( 8) /  24.00000000000 /,
     *     ICNTX( 8) /608 /,
     *     ICNTY( 8) /288 /,
     *     NBASE( 8) /     0 /
   
      END
