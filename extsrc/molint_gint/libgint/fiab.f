      SUBROUTINE FIAB(X,Y,NX,NY,Z,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(2),Y(2),Z(2)
      IF(NY) 1,1,2
    1 DO 3 K=1,N
    3 Z(K)=Y(1)*X(K)
      RETURN
    2 DO 4 K=1,NY+1
    4 Z(K)=X(1)*Y(K)
      IF(NX) 5,5,6
    5 RETURN
    6 DO 7 K=1,NX
      M=1
      KMAX=K+NY-1
      DO 8 L=K,KMAX
      Z(L+1)=Z(L+1)+X(K+1)*Y(M)
    8 M=M+1
    7 Z(KMAX+2)=X(K+1)*Y(M)
      RETURN
      END