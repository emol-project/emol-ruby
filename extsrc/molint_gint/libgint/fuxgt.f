  !     SUBROUTINE FUXGT
  !     IMPLICIT REAL*8 (A-H,O-Z)
  !     COMMON/FFDG/T,G(24),DFAC(24),GICA(81),M
  !     COMMON/FPR/ISW
  !     DIMENSION TABG(24,101)
  !     DIMENSION GZ(24)
  !     EQUIVALENCE (TABG,GZ)
  !     SAVE TABG
  !     IF(ISW.EQ.0) CALL FPRSET(TABG)
  !     IF(T.LT.1.0D-17)GO TO 8
  !     IF(T.LT.20.0999D+00) GO TO 6
  !     X=0.5D+00/T
  !     IF(T.GE.35.0D0) GO TO 33
  !                       AT1=  DEXP(-T)
  !     G(1)=.88622692545275801D0/DSQRT(T) -
  !    &  AT1*X*(1.0D0-X*(1.0D0-X*(3.0D0-X*(5.0D0-X*7.0D0))))
  !  34 IF(M.EQ.1) RETURN
  !     DO 2 I=2,M
  !   2 G(I)=G(I-1)*X*DFLOAT(2*I-3)-AT1   *X
  !     RETURN
  !  33 AT1=0.0D+00
  !     G(1)=.88622692545275801D0/DSQRT(T)
  !     GO TO 34
  !   6 AT1=T*5.0D+00
  !     L=AT1+1.5D+00
  !     X=T-0.2D+00*DFLOAT(L-1)
  !     G(M)=(((((((2.480158730D-5*TABG(M+8,L)*X-1.984126984D-4*TABG(
  !    &  M+7,L))*X  +
  !    &      .13888888888888889D-2*TABG(M+6,L))*X
  !    &  -.83333333333333333D-2*
  !    &  TABG(M+5,L))*X+.41666666666666667D-1*TABG(M+4,L))*X-
  !    &   .16666666666666667D0
  !    &  *TABG(M+3,L))*X+0.5D+00*TABG(M+2,L))*X-TABG(M+1,L))*X
  !    &  +TABG(M,L)
  !     IF(M.EQ.1) RETURN
  !     X=T*2.0D+00
  !     CEX=DEXP(-T)
  !     DO 100 I=2,M
  !     NMI=M-I+1
  ! 100 G(NMI)=(X*G(NMI+1)+CEX)/DFLOAT(2*NMI-1)
  !     RETURN
  !   8 DO 9 L=1,M
  !   9 G(L)=GZ(L)
  !     RETURN
  !     END

C ------ libgint/fuxgt.f ------
      SUBROUTINE FUXGT
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FFDG/T,G(24),DFAC(24),GICA(81),M
      COMMON/FPR/ISW
      REAL*8 F_TEMP(0:24)
      INTEGER I
      INTEGER MAX_L_VAL, M_VAL
      REAL*8 COEF

C     Create Table
      IF(ISW.EQ.0) THEN
         MAX_L_VAL = 4
         CALL FMT_INITIALIZE(MAX_L_VAL)
         ISW = 1
      END IF

C     T ~= 0 (F_m(0) = 1 / (2m+1))
      IF (T .LT. 1.0D-17) THEN
         DO 10 I = 1, M
            G(I) = 1.0D0 / DFLOAT(2*(I-1) + 1)
   10    CONTINUE
         RETURN
      END IF

C     call fast fmt
      M_VAL = M - 1
      COEF = 1.0D0
      CALL FMT(F_TEMP, M_VAL, T, COEF)

C
      DO 20 I = 1, M
         G(I) = F_TEMP(I-1)
   20 CONTINUE

      RETURN
      END