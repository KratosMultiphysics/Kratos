      SUBROUTINE SROTMG(SD1,SD2,SX1,SY1,SPARAM)
*     .. Scalar Arguments ..
      REAL SD1,SD2,SX1,SY1
*     ..
*     .. Array Arguments ..
      REAL SPARAM(5)
*     ..
*
*  Purpose
*  =======
*
*     CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
*     THE SECOND COMPONENT OF THE 2-VECTOR  (SQRT(SD1)*SX1,SQRT(SD2)*
*     SY2)**T.
*     WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS..
*
*     SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0
*
*       (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)
*     H=(          )    (          )    (          )    (          )
*       (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).
*     LOCATIONS 2-4 OF SPARAM CONTAIN SH11,SH21,SH12, AND SH22
*     RESPECTIVELY. (VALUES OF 1.E0, -1.E0, OR 0.E0 IMPLIED BY THE
*     VALUE OF SPARAM(1) ARE NOT STORED IN SPARAM.)
*
*     THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
*     INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
*     OF SD1 AND SD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.
*
*
*  Arguments
*  =========
*
*
*  SD1    (input/output) REAL
*
*  SD2    (input/output) REAL
*
*  SX1    (input/output) REAL
*
*  SY1    (input) REAL
*
*
*  SPARAM (input/output)  REAL array, dimension 5
*     SPARAM(1)=SFLAG
*     SPARAM(2)=SH11
*     SPARAM(3)=SH21
*     SPARAM(4)=SH12
*     SPARAM(5)=SH22
*
*  =====================================================================
*
*     .. Local Scalars ..
      REAL GAM,GAMSQ,ONE,RGAMSQ,SFLAG,SH11,SH12,SH21,SH22,SP1,SP2,SQ1,
     +     SQ2,STEMP,SU,TWO,ZERO
      INTEGER IGO
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS
*     ..
*     .. Data statements ..
*
      DATA ZERO,ONE,TWO/0.E0,1.E0,2.E0/
      DATA GAM,GAMSQ,RGAMSQ/4096.E0,1.67772E7,5.96046E-8/
*     ..

      IF (.NOT.SD1.LT.ZERO) GO TO 10
*       GO ZERO-H-D-AND-SX1..
      GO TO 60
   10 CONTINUE
*     CASE-SD1-NONNEGATIVE
      SP2 = SD2*SY1
      IF (.NOT.SP2.EQ.ZERO) GO TO 20
      SFLAG = -TWO
      GO TO 260
*     REGULAR-CASE..
   20 CONTINUE
      SP1 = SD1*SX1
      SQ2 = SP2*SY1
      SQ1 = SP1*SX1
*
      IF (.NOT.ABS(SQ1).GT.ABS(SQ2)) GO TO 40
      SH21 = -SY1/SX1
      SH12 = SP2/SP1
*
      SU = ONE - SH12*SH21
*
      IF (.NOT.SU.LE.ZERO) GO TO 30
*         GO ZERO-H-D-AND-SX1..
      GO TO 60
   30 CONTINUE
      SFLAG = ZERO
      SD1 = SD1/SU
      SD2 = SD2/SU
      SX1 = SX1*SU
*         GO SCALE-CHECK..
      GO TO 100
   40 CONTINUE
      IF (.NOT.SQ2.LT.ZERO) GO TO 50
*         GO ZERO-H-D-AND-SX1..
      GO TO 60
   50 CONTINUE
      SFLAG = ONE
      SH11 = SP1/SP2
      SH22 = SX1/SY1
      SU = ONE + SH11*SH22
      STEMP = SD2/SU
      SD2 = SD1/SU
      SD1 = STEMP
      SX1 = SY1*SU
*         GO SCALE-CHECK
      GO TO 100
*     PROCEDURE..ZERO-H-D-AND-SX1..
   60 CONTINUE
      SFLAG = -ONE
      SH11 = ZERO
      SH12 = ZERO
      SH21 = ZERO
      SH22 = ZERO
*
      SD1 = ZERO
      SD2 = ZERO
      SX1 = ZERO
*         RETURN..
      GO TO 220
*     PROCEDURE..FIX-H..
   70 CONTINUE
      IF (.NOT.SFLAG.GE.ZERO) GO TO 90
*
      IF (.NOT.SFLAG.EQ.ZERO) GO TO 80
      SH11 = ONE
      SH22 = ONE
      SFLAG = -ONE
      GO TO 90
   80 CONTINUE
      SH21 = -ONE
      SH12 = ONE
      SFLAG = -ONE
   90 CONTINUE
      GO TO IGO(120,150,180,210)
*     PROCEDURE..SCALE-CHECK
  100 CONTINUE
  110 CONTINUE
      IF (.NOT.SD1.LE.RGAMSQ) GO TO 130
      IF (SD1.EQ.ZERO) GO TO 160
      ASSIGN 120 TO IGO
*              FIX-H..
      GO TO 70
  120 CONTINUE
      SD1 = SD1*GAM**2
      SX1 = SX1/GAM
      SH11 = SH11/GAM
      SH12 = SH12/GAM
      GO TO 110
  130 CONTINUE
  140 CONTINUE
      IF (.NOT.SD1.GE.GAMSQ) GO TO 160
      ASSIGN 150 TO IGO
*              FIX-H..
      GO TO 70
  150 CONTINUE
      SD1 = SD1/GAM**2
      SX1 = SX1*GAM
      SH11 = SH11*GAM
      SH12 = SH12*GAM
      GO TO 140
  160 CONTINUE
  170 CONTINUE
      IF (.NOT.ABS(SD2).LE.RGAMSQ) GO TO 190
      IF (SD2.EQ.ZERO) GO TO 220
      ASSIGN 180 TO IGO
*              FIX-H..
      GO TO 70
  180 CONTINUE
      SD2 = SD2*GAM**2
      SH21 = SH21/GAM
      SH22 = SH22/GAM
      GO TO 170
  190 CONTINUE
  200 CONTINUE
      IF (.NOT.ABS(SD2).GE.GAMSQ) GO TO 220
      ASSIGN 210 TO IGO
*              FIX-H..
      GO TO 70
  210 CONTINUE
      SD2 = SD2/GAM**2
      SH21 = SH21*GAM
      SH22 = SH22*GAM
      GO TO 200
  220 CONTINUE
      IF (SFLAG) 250,230,240
  230 CONTINUE
      SPARAM(3) = SH21
      SPARAM(4) = SH12
      GO TO 260
  240 CONTINUE
      SPARAM(2) = SH11
      SPARAM(5) = SH22
      GO TO 260
  250 CONTINUE
      SPARAM(2) = SH11
      SPARAM(3) = SH21
      SPARAM(4) = SH12
      SPARAM(5) = SH22
  260 CONTINUE
      SPARAM(1) = SFLAG
      RETURN
      END
