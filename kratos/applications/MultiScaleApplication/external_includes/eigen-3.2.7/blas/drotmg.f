      SUBROUTINE DROTMG(DD1,DD2,DX1,DY1,DPARAM)
*     .. Scalar Arguments ..
      DOUBLE PRECISION DD1,DD2,DX1,DY1
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DPARAM(5)
*     ..
*
*  Purpose
*  =======
*
*     CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
*     THE SECOND COMPONENT OF THE 2-VECTOR  (DSQRT(DD1)*DX1,DSQRT(DD2)*
*     DY2)**T.
*     WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
*
*     DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
*
*       (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
*     H=(          )    (          )    (          )    (          )
*       (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
*     LOCATIONS 2-4 OF DPARAM CONTAIN DH11, DH21, DH12, AND DH22
*     RESPECTIVELY. (VALUES OF 1.D0, -1.D0, OR 0.D0 IMPLIED BY THE
*     VALUE OF DPARAM(1) ARE NOT STORED IN DPARAM.)
*
*     THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
*     INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
*     OF DD1 AND DD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.
*
*
*  Arguments
*  =========
*
*  DD1    (input/output) DOUBLE PRECISION
*
*  DD2    (input/output) DOUBLE PRECISION 
*
*  DX1    (input/output) DOUBLE PRECISION 
*
*  DY1    (input) DOUBLE PRECISION
*
*  DPARAM (input/output)  DOUBLE PRECISION array, dimension 5
*     DPARAM(1)=DFLAG
*     DPARAM(2)=DH11
*     DPARAM(3)=DH21
*     DPARAM(4)=DH12
*     DPARAM(5)=DH22
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DFLAG,DH11,DH12,DH21,DH22,DP1,DP2,DQ1,DQ2,DTEMP,
     +                 DU,GAM,GAMSQ,ONE,RGAMSQ,TWO,ZERO
      INTEGER IGO
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DABS
*     ..
*     .. Data statements ..
*
      DATA ZERO,ONE,TWO/0.D0,1.D0,2.D0/
      DATA GAM,GAMSQ,RGAMSQ/4096.D0,16777216.D0,5.9604645D-8/
*     ..

      IF (.NOT.DD1.LT.ZERO) GO TO 10
*       GO ZERO-H-D-AND-DX1..
      GO TO 60
   10 CONTINUE
*     CASE-DD1-NONNEGATIVE
      DP2 = DD2*DY1
      IF (.NOT.DP2.EQ.ZERO) GO TO 20
      DFLAG = -TWO
      GO TO 260
*     REGULAR-CASE..
   20 CONTINUE
      DP1 = DD1*DX1
      DQ2 = DP2*DY1
      DQ1 = DP1*DX1
*
      IF (.NOT.DABS(DQ1).GT.DABS(DQ2)) GO TO 40
      DH21 = -DY1/DX1
      DH12 = DP2/DP1
*
      DU = ONE - DH12*DH21
*
      IF (.NOT.DU.LE.ZERO) GO TO 30
*         GO ZERO-H-D-AND-DX1..
      GO TO 60
   30 CONTINUE
      DFLAG = ZERO
      DD1 = DD1/DU
      DD2 = DD2/DU
      DX1 = DX1*DU
*         GO SCALE-CHECK..
      GO TO 100
   40 CONTINUE
      IF (.NOT.DQ2.LT.ZERO) GO TO 50
*         GO ZERO-H-D-AND-DX1..
      GO TO 60
   50 CONTINUE
      DFLAG = ONE
      DH11 = DP1/DP2
      DH22 = DX1/DY1
      DU = ONE + DH11*DH22
      DTEMP = DD2/DU
      DD2 = DD1/DU
      DD1 = DTEMP
      DX1 = DY1*DU
*         GO SCALE-CHECK
      GO TO 100
*     PROCEDURE..ZERO-H-D-AND-DX1..
   60 CONTINUE
      DFLAG = -ONE
      DH11 = ZERO
      DH12 = ZERO
      DH21 = ZERO
      DH22 = ZERO
*
      DD1 = ZERO
      DD2 = ZERO
      DX1 = ZERO
*         RETURN..
      GO TO 220
*     PROCEDURE..FIX-H..
   70 CONTINUE
      IF (.NOT.DFLAG.GE.ZERO) GO TO 90
*
      IF (.NOT.DFLAG.EQ.ZERO) GO TO 80
      DH11 = ONE
      DH22 = ONE
      DFLAG = -ONE
      GO TO 90
   80 CONTINUE
      DH21 = -ONE
      DH12 = ONE
      DFLAG = -ONE
   90 CONTINUE
      GO TO IGO(120,150,180,210)
*     PROCEDURE..SCALE-CHECK
  100 CONTINUE
  110 CONTINUE
      IF (.NOT.DD1.LE.RGAMSQ) GO TO 130
      IF (DD1.EQ.ZERO) GO TO 160
      ASSIGN 120 TO IGO
*              FIX-H..
      GO TO 70
  120 CONTINUE
      DD1 = DD1*GAM**2
      DX1 = DX1/GAM
      DH11 = DH11/GAM
      DH12 = DH12/GAM
      GO TO 110
  130 CONTINUE
  140 CONTINUE
      IF (.NOT.DD1.GE.GAMSQ) GO TO 160
      ASSIGN 150 TO IGO
*              FIX-H..
      GO TO 70
  150 CONTINUE
      DD1 = DD1/GAM**2
      DX1 = DX1*GAM
      DH11 = DH11*GAM
      DH12 = DH12*GAM
      GO TO 140
  160 CONTINUE
  170 CONTINUE
      IF (.NOT.DABS(DD2).LE.RGAMSQ) GO TO 190
      IF (DD2.EQ.ZERO) GO TO 220
      ASSIGN 180 TO IGO
*              FIX-H..
      GO TO 70
  180 CONTINUE
      DD2 = DD2*GAM**2
      DH21 = DH21/GAM
      DH22 = DH22/GAM
      GO TO 170
  190 CONTINUE
  200 CONTINUE
      IF (.NOT.DABS(DD2).GE.GAMSQ) GO TO 220
      ASSIGN 210 TO IGO
*              FIX-H..
      GO TO 70
  210 CONTINUE
      DD2 = DD2/GAM**2
      DH21 = DH21*GAM
      DH22 = DH22*GAM
      GO TO 200
  220 CONTINUE
      IF (DFLAG) 250,230,240
  230 CONTINUE
      DPARAM(3) = DH21
      DPARAM(4) = DH12
      GO TO 260
  240 CONTINUE
      DPARAM(2) = DH11
      DPARAM(5) = DH22
      GO TO 260
  250 CONTINUE
      DPARAM(2) = DH11
      DPARAM(3) = DH21
      DPARAM(4) = DH12
      DPARAM(5) = DH22
  260 CONTINUE
      DPARAM(1) = DFLAG
      RETURN
      END
