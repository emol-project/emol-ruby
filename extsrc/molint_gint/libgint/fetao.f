C     ================================================================
C     Original code contained:
C       SUBROUTINE FETAO
C       ENTRY FETREQ(KQ)
C       ENTRY FVALUE
C     Split into separate subroutines for gfortran compatibility.
C     ================================================================

      SUBROUTINE FETAO
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FFLA/XA(3),XB(3),XC(3),XD(3),AOE(20),ACE(20),BOE(20),
     &  BCE(20),COE(20),CCE(20),DOE(20),DCE(20),NA(3),NB(3),NC(3),
     &  ND(3),ITA,ITB,ITC,ITD
      COMMON/FAREA/OE(20,102),CE(20,102),XSD(3,100),NQND(4,102),NAO,
     &  NPS,LA,LPA,LB,LPB
      COMMON/FNTV/W,WS,WT,ONA(100),IA,IPA,IB,IPB,IC,IPC,ID,IPD,
     &  NAP(100),ITIME
      DIMENSION NQW(3),CET(20)
      COMMON/FNTVS/WM(3),WN(6)
      COMMON/FNTVT/ONAF(3,100),ONAFG(6,100)
      DOUBLE PRECISION VDX
      COMMON /FOCI/ VDX,WOCI,IUX
C     DATA LA/-1/
      NMOD=0
      IU=0
      KFORMT=0
C     READ(5,*) NAO,NPS,NMOD,IU,KFORMT
      READ(5,100) NAO,NPS,NMOD,IU,KFORMT
      IUX=IU-1
      IF(NMOD.EQ.0) GO TO 5
      IF(KFORMT.EQ.0) READ(5,111) (ONA(I),I=1,NMOD)
      IF(KFORMT.EQ.1) READ(5,101) (ONA(I),I=1,NMOD)
      IF(KFORMT.EQ.2) READ(5,*) (ONA(I),I=1,NMOD)
      WRITE(6,200) (I,ONA(I),I=1,NMOD)
  200 FORMAT(///5X,53('*'),' ORBITAL MODIFIERS '
     &  /(/5X,6(I3,F10.5)))
  100 FORMAT(10I5)
  101 FORMAT(6F12.6)
  111 FORMAT(5F16.8)
    5 WRITE(6,202)
  202 FORMAT(///5X,55('*'),' ATOMIC ORBITALS  ')
      DO 3 I=1,NAO
      IF(KFORMT.NE.2) READ(5,100) (NQND(J,I),J=1,4),NMFI
      IF(KFORMT.EQ.2) READ(5,*) (NQND(J,I),J=1,4),NMFI
      NTET=NQND(4,I)
      IF(NTET.GT.20) WRITE(6,999) NTET
  999 FORMAT(20('*'),'ERROR ***** NO OF TERMS',I5)
      IF(KFORMT.EQ.0) READ(5,111) (OE(J,I),J=1,NTET)
      IF(KFORMT.EQ.0) READ(5,111) (CE(J,I),J=1,NTET)
      IF(KFORMT.EQ.1) READ(5,101) (OE(J,I),J=1,NTET)
      IF(KFORMT.EQ.1) READ(5,101) (CE(J,I),J=1,NTET)
      IF(KFORMT.EQ.2) READ(5,*) (OE(J,I),J=1,NTET)
      IF(KFORMT.EQ.2) READ(5,*) (CE(J,I),J=1,NTET)
      WRITE(6,201) I,(NQND(J,I),J=1,3),NMFI,(OE(J,I),J=1,NTET)
  201 FORMAT(///I8,I13,2I3,I4/(/5X,6F16.10))
      WRITE(6,203) (CE(J,I),J=1,NTET)
  203 FORMAT(/5X,6F16.10)
      NTERM=NQND(4,I)
      DO 4 J=1,NTERM
      IF(NMFI.EQ.0) GO TO 4
      OE(J,I)=OE(J,I)*ONA(NMFI)**2
    4 CE(J,I)=CE(J,I)*FNFB(NQND(1,I),OE(J,I))
    3 CONTINUE
      DO 1 J=1,NPS
      IF(KFORMT.EQ.0) READ(5,111) (XSD(I,J),I=1,3)
      IF(KFORMT.EQ.1) READ(5,101) (XSD(I,J),I=1,3)
    1 IF(KFORMT.EQ.2) READ(5,*) (XSD(I,J),I=1,3)
      WRITE(6,204) (J,(XSD(I,J),I=1,3),J=1,NPS)
  204 FORMAT(///5X, 53('*'),' NUCLEAR POSITIONS '
     &  /(/I8,3F16.10))
      WRITE(6,205)
  205 FORMAT(///5X,61('*'),' INTEGRALS ')
      RETURN
C
C     --- Original ENTRY FETREQ(KQ) was here ---
C     ENTRY FETREQ(KQ)
C     READ(5,100) IA,IPA,IB,IPB,IC,IPC,ID,IPD
C     KQ=IA
C     IF(IA.EQ.0) RETURN
C     IF(IPC.NE.0.OR.IC.EQ.0) GO TO 6
C     READ(5,100) (NAP(I),I=1,IC)
C
C     --- Original ENTRY FVALUE was here ---
C     ENTRY FVALUE
C     (falls through to shared code at label 6)
C
      END

C     ================================================================
C     FETREQ: split from ENTRY FETREQ(KQ) in original FETAO
C     ================================================================
      SUBROUTINE FETREQ(KQ)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FFLA/XA(3),XB(3),XC(3),XD(3),AOE(20),ACE(20),BOE(20),
     &  BCE(20),COE(20),CCE(20),DOE(20),DCE(20),NA(3),NB(3),NC(3),
     &  ND(3),ITA,ITB,ITC,ITD
      COMMON/FAREA/OE(20,102),CE(20,102),XSD(3,100),NQND(4,102),NAO,
     &  NPS,LA,LPA,LB,LPB
      COMMON/FNTV/W,WS,WT,ONA(100),IA,IPA,IB,IPB,IC,IPC,ID,IPD,
     &  NAP(100),ITIME
      COMMON/FNTVS/WM(3),WN(6)
      COMMON/FNTVT/ONAF(3,100),ONAFG(6,100)
      DOUBLE PRECISION VDX
      COMMON /FOCI/ VDX,WOCI,IUX

      READ(5,100) IA,IPA,IB,IPB,IC,IPC,ID,IPD
  100 FORMAT(10I5)
      KQ=IA
      IF(IA.EQ.0) RETURN
      IF(IPC.NE.0.OR.IC.EQ.0) GO TO 6
      READ(5,100) (NAP(I),I=1,IC)
C     --- Fall through to shared code (originally ENTRY FVALUE) ---
    6 CONTINUE
      CALL FVALUE
      RETURN
      END

C     ================================================================
C     FVALUE: split from ENTRY FVALUE in original FETAO
C     This contains the shared computation code used by both
C     FETREQ and FVALUE entry points.
C     ================================================================
      SUBROUTINE FVALUE
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FFLA/XA(3),XB(3),XC(3),XD(3),AOE(20),ACE(20),BOE(20),
     &  BCE(20),COE(20),CCE(20),DOE(20),DCE(20),NA(3),NB(3),NC(3),
     &  ND(3),ITA,ITB,ITC,ITD
      COMMON/FAREA/OE(20,102),CE(20,102),XSD(3,100),NQND(4,102),NAO,
     &  NPS,LA,LPA,LB,LPB
      COMMON/FNTV/W,WS,WT,ONA(100),IA,IPA,IB,IPB,IC,IPC,ID,IPD,
     &  NAP(100),ITIME
      DIMENSION CET(20)
      COMMON/FNTVS/WM(3),WN(6)
      COMMON/FNTVT/ONAF(3,100),ONAFG(6,100)
      DOUBLE PRECISION VDX
      COMMON /FOCI/ VDX,WOCI,IUX

C   6 CALL CLOCKM(ISTIM)
      IF((IA-NAO-1)*IA) 7,998,998
    7 IF((IB-NAO-1)*IB) 8,998,998
    8 IF((IPA-NPS-1)*IPA) 9,998,998
    9 IF((IPB-NPS-1)*IPB) 10,998,998
   10 IF(IA.EQ.LA.AND.IB.EQ.LB.AND.IPA.EQ.LPA.AND.IPB.EQ.LPB)
     &  GO TO 24
      CALL FTRA(IA,IPA,1)
      CALL FTRA(IB,IPB,2)
      LA=IA
      LB=IB
      LPA=IPA
      LPB=IPB
      KST=1
   25 IF(IPC.EQ.0) GO TO 11
      IF((IC-NAO-1)*IC.GE.0.OR.(IPC-NPS-1)*IPC.GE.0.OR.
     &   (ID-NAO-1)*ID.GE.0.OR.(IPD-NPS-1)*IPD.GE.0) GO TO 998
      CALL FTRA(IC,IPC,3)
      CALL FTRA(ID,IPD,4)
      GO TO (21,22),KST
   21 CALL FINTA
      GO TO 50
   22 CALL FINTB
      GO TO 50
   24 KST=2
      GO TO 25
   11 IF(ID.EQ.3) GO TO 60
      DO 12 I=1,ITA
   12 CET(I)=ACE(I)
      IUX=-1
      WT=0.0D+00
      DO 14 IT=1,3
      NA(IT)=NA(IT)-2
      IF( NA(IT).LT.0) GO TO 15
      DO 16 I=1,ITA
   16 ACE(I)=-CET(I)*FLOAT((NA(IT)+1)*(NA(IT)+2))
      CALL FINTC
      WT=WT+WS
   15 NA(IT)=NA(IT)+4
      DO 17 I=1,ITA
   17 ACE(I)=-CET(I)*4.0D+00*AOE(I)**2
      CALL FINTD
      WT=WT+WS
   14 NA(IT)=NA(IT)-2
      DO 19 I=1,ITA
   19 ACE(I)=CET(I)
      CALL FINTC
      WT=   W*FLOAT(2*(NA(1)+NA(2)+NA(3))+3)+WT*0.5D+00
      IF(IPD.EQ.0) GO TO 30
      SSS=WS
      DO 31 IT=1,3
      NA(IT)=NA(IT)+1
      CALL FINTD
      WM(IT)=WS
   31 NA(IT)=NA(IT)-1
      IW=0
      DO 32 IT=1,3
      NA(IT)=NA(IT)+1
      XYZT=XSD(IT,IPA)-XSD(IT,IPD)
      DO 33 IU=IT,3
      NA(IU)=NA(IU)+1
      XYZU=XSD(IU,IPA)-XSD(IU,IPD)
      IW=IW+1
      CALL FINTD
      WN(IW)=WS+XYZU*WM(IT)+XYZT*WM(IU)+XYZT*XYZU*SSS
   33 NA(IU)=NA(IU)-1
   32 NA(IT)=NA(IT)-1
      DO 34 IT=1,3
      XYZ=XSD(IT,IPA)-XSD(IT,IPD)
   34 WM(IT)=WM(IT)+XYZ*SSS
      WS=SSS
      IF(IC.NE.0) CALL FINTC
   30 CONTINUE
      IF(IC.EQ.0) GO TO 50
      DO 20 IT=1,IC
      NAPIT=NAP(IT)
      DO 35 I=1,3
      XC(I)=XSD(I,NAPIT)
   35 NC(I)=0
      CALL FINTE
      ONA(IT)=W
      IF(ID.EQ.0) GO TO 20
      DO 36 I=1,3
      NC(I)=1
      CALL FINTE
      ONAF(I,IT)=W
   36 NC(I)=0
      IF(ID.EQ.1) GO TO 20
      IU=0
      DO 37 I=1,3
      NC(I)=1
      DO 38 J=I,3
      NC(J)=NC(J)+1
      IU=IU+1
      CALL FINTE
      ONAFG(IU,IT)=W
   38 NC(J)=NC(J)-1
   37 NC(I)=0
   20 CONTINUE
      GO TO 50
   60 CALL FINTC
      IUX=-1
      IF(IC.LE.0) GO TO 998
      DO 62 IT=1,IC
      NAPIT=NAP(IT)
      DO 64 I=1,3
      XC(I)=XSD(I,NAPIT)
   64 NC(I)=0
      DO 63 I=1,3
      NC(I)=1
      CALL FINTE
      ONAF(I,IT)=W
   63 NC(I)=0
      IU=0
      DO 66 I=1,3
      NC(I)=1
      DO 68 J=I,3
      NC(J)=NC(J)+1
      IU=IU+1
      CALL FINTE
      ONAFG(IU,IT)=W
   68 NC(J)=NC(J)-1
   66 NC(I)=0
   62 CONTINUE
      GO TO 50
  998 WRITE(6,997) IA,IPA,IB,IPB,IC,IPC,ID,IPD
  997 FORMAT(20('*'),'ERROR ***** REQUESTED INTEGRAL',8I5)
      W=1.0D+40
      WS=1.0D+40
      WT=1.0D+40
C  50 CALL CLOCKM(IFTIM)
   50 CONTINUE
C     ITIME=IFTIM-ISTIM
      RETURN
      END
