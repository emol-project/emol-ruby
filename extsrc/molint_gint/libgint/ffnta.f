      SUBROUTINE FFNTA(XF,NT,SX,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FFDG/XX(49),GICA(82)
      DIMENSION XF(2),SX(2),W(2)
      J=0
      DO 2 L=1,NT
      W(L)=XF(L)*SX(L)
      J=J+1
      M=L-2
      IF(M) 2,2,4
    4 DO 6 K=3,L,2
      IT=L-K/2
      J=J+1
      W(IT)=W(IT)+XF(L)*SX(M)*GICA(J)
    6 M=M-2
    2 CONTINUE
      RETURN
      END
