      SUBROUTINE FIAD(U,XA,G,NI,NF)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(7)
      DO 1 I=NI,NF
      XA(I+1)=XA(I)*G
      T=0.0D+00
      IF(I.EQ.1) GO TO 3
      DO 2 JX=2,I
      J=I-JX+2
      B=XA(J)*DFLOAT(J-1)
      XA(J)=XA(J-1)*G+XA(J)*U+T
    2 T=B
    3 XA(1)=XA(1)*U+T
    1 CONTINUE
      RETURN
      END
