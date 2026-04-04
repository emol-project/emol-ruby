C     ================================================================
C     Original code contained:
C       SUBROUTINE FINTA
C       ENTRY FINTB
C       ENTRY FINTC
C       ENTRY FINTD
C       ENTRY FINTE
C     Split into separate subroutines for gfortran compatibility.
C     A new COMMON /FIKI/ is introduced to share the local variables
C     IKI and W that were shared across all entry points.
C     ================================================================

C     ----------------------------------------------------------------
C     FINTA: main entry - computes AB pair loop, then falls through
C            to FINTB logic.
C     ----------------------------------------------------------------
      SUBROUTINE FINTA
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FFLA/XA(3),XB(3),XC(3),XD(3),AOE(20),ACE(20),BOE(20),
     &  BCE(20),COE(20),CCE(20),DOE(20),DCE(20),NA(3),NB(3),NC(3),
     &  ND(3),ITA,ITB,ITC,ITD
      COMMON/FFDG/XX,FNU(24),DFAC(24),GICA(81),NTSUM1
      COMMON/FNTV/E,ES,ET,ONA(100),IA,IPA,IB,IPB,IC,IPC,ID,IPD,
     &  NAP(100),ITIME
      COMMON/FPAR/V,XEA(5200),XFA(17,3),XF(17,3),XPA(3,400),
     &  SP(17,3),XQ(3),G1A(400),G2,AMF(400),MXDT,NAB(3),NCD(3),
     &  NT(3),NKT(3),NTXYM1
      DOUBLE PRECISION VDX,WDX
      COMMON /FOCI/ VDX,WOCI,IUX
      COMMON /FIKI/ W,IKI

      IKI=0
      CALL FINTA_CORE
C     --- Original code fell through to ENTRY FINTB here ---
      CALL FINTB
      RETURN
C
C     --- Original ENTRY FINTB was here ---
C     ENTRY FINTB
C
C     --- Original ENTRY FINTC was here ---
C     ENTRY FINTC
C     IKI=1
C     W=0.0D+00
C     V=0.0D+00
C     GO TO 41
C  43 E =W*0.66566768190019485D0
C     ES=V*0.66566768190019485D0
C     RETURN
C
C     --- Original ENTRY FINTD was here ---
C     ENTRY FINTD
C     IKI=2
C     W=0.0D+00
C     V=0.0D+00
C     GO TO 41
C
C     --- Original ENTRY FINTE was here ---
C     ENTRY FINTE
C     (separate subroutine below)
C
      END

C     ----------------------------------------------------------------
C     FINTA_CORE: the AB-pair loop (originally label 41 onward in
C     FINTA). Used by FINTA, FINTC, and FINTD with different IKI.
C     W is shared via COMMON /FIKI/ since it was a local variable
C     in the original code shared across all entry points.
C     ----------------------------------------------------------------
      SUBROUTINE FINTA_CORE
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FFLA/XA(3),XB(3),XC(3),XD(3),AOE(20),ACE(20),BOE(20),
     &  BCE(20),COE(20),CCE(20),DOE(20),DCE(20),NA(3),NB(3),NC(3),
     &  ND(3),ITA,ITB,ITC,ITD
      COMMON/FFDG/XX,FNU(24),DFAC(24),GICA(81),NTSUM1
      COMMON/FNTV/E,ES,ET,ONA(100),IA,IPA,IB,IPB,IC,IPC,ID,IPD,
     &  NAP(100),ITIME
      COMMON/FPAR/V,XEA(5200),XFA(17,3),XF(17,3),XPA(3,400),
     &  SP(17,3),XQ(3),G1A(400),G2,AMF(400),MXDT,NAB(3),NCD(3),
     &  NT(3),NKT(3),NTXYM1
      COMMON /FIKI/ W,IKI
      EQUIVALENCE (NKT,NTM1),(NKT(2),NTM2),(NKT(3),NTM3)

C     --- label 41 in original ---
      SAB=0.0D+00
      DO 30 I=1,3
      NAB(I)=NA(I)+NB(I)
   30 SAB=SAB+(XA(I)-XB(I))**2
      MXD=0
      KT=1
      DO 3 JA=1,ITA
      IF(IKI.EQ.2) KT=1
      EPA=AOE(JA)
      DO 3 JB=1,ITB
      EPB=BOE(JB)
      G1=EPA+EPB
      AFX=SAB*EPA*EPB/G1
      IF(AFX.GT.100.0D+0) GO TO 3
      MXD=MXD+1
      AMF(MXD)=DEXP(-AFX)*ACE(JA)*BCE(JB)
      Q1=0.5D+00/G1
      G1A(MXD)=G1
      DO 4 I=1,3
      IF(XA(I).EQ.XB(I)) GO TO 17
      XPA(I,MXD)=(EPA*XA(I)+EPB*XB(I))/G1
      GO TO 18
   17 XPA(I,MXD)=XA(I)
   18 XPI=XPA(I,MXD)
      XEA(KT)=1.0D+00
      IF(NA(I).NE.0) CALL FIAD(XPI-XA(I),XEA(KT),Q1,1,NA(I))
      IF(NB(I).NE.0)
     &  CALL FIAD(XPI-XB(I),XEA(KT),Q1,NA(I)+1,NAB(I))
      NKT(I)=KT
      KT=KT+1+NAB(I)
    4 CONTINUE
      IF(IKI) 3,3,40
   40 U=AMF(MXD)*XEA(NTM1)*XEA(NTM2)*XEA(NTM3)/(DSQRT(G1)*G1)
      W=W+U*EPA
      V=V+U
    3 CONTINUE
      IF(KT.GT.5201) WRITE(6,999) KT
  999 FORMAT(20('*'),' ERROR  (KT) EXCEEDS 5200 ',I10,' *****')
      MXDT=MXD
      IF(IKI.GT.0) THEN
C       --- label 43 in original ---
        E =W*0.66566768190019485D0
        ES=V*0.66566768190019485D0
      ENDIF
      RETURN
      END

C     ----------------------------------------------------------------
C     FINTB: split from ENTRY FINTB in original FINTA.
C     CD-pair loop with FYCLE integration.
C     ----------------------------------------------------------------
      SUBROUTINE FINTB
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FFLA/XA(3),XB(3),XC(3),XD(3),AOE(20),ACE(20),BOE(20),
     &  BCE(20),COE(20),CCE(20),DOE(20),DCE(20),NA(3),NB(3),NC(3),
     &  ND(3),ITA,ITB,ITC,ITD
      COMMON/FFDG/XX,FNU(24),DFAC(24),GICA(81),NTSUM1
      COMMON/FNTV/E,ES,ET,ONA(100),IA,IPA,IB,IPB,IC,IPC,ID,IPD,
     &  NAP(100),ITIME
      COMMON/FPAR/V,XEA(5200),XFA(17,3),XF(17,3),XPA(3,400),
     &  SP(17,3),XQ(3),G1A(400),G2,AMF(400),MXDT,NAB(3),NCD(3),
     &  NT(3),NKT(3),NTXYM1
      DOUBLE PRECISION VDX,WDX
      COMMON /FOCI/ VDX,WOCI,IUX
      COMMON /FIKI/ W,IKI
      EQUIVALENCE (NKT(1),NTM1),(NKT(2),NTM2),(NKT(3),NTM3)

      SCD=0.0D+00
      DO 1 I=1,3
      SP(1,I)=1.0D+00
      NCD(I)=NC(I)+ND(I)
      NT(I)=NAB(I)+NCD(I)+1
    1 SCD=SCD+(XC(I)-XD(I))**2
      NTSUM1=NT(1)+NT(2)+NT(3)-2
      NTM1=NT(1)-1
      NTM2=NT(2)-1
      NTM3=NT(3)-1
      NTXYM1=NTM1+NTM2
      W=0.0D+00
      WDX=0.0D+00
      DO 6 JC=1,ITC
      EPC=COE(JC)
      DO 6 JD=1,ITD
      EPD=DOE(JD)
      G2=EPC+EPD
      AW=SCD*EPC*EPD/G2
      IF(AW.GT.100.0D+0) GO TO 6
      Q2=0.5D+00/G2
      DO 10 I=1,3
      IF(XC(I).EQ.XD(I)) GO TO 15
      XQ(I)=(EPC*XC(I)+EPD*XD(I))/G2
      GO TO 16
   15 XQ(I)=XC(I)
   16 XFA(1,I)=1.0D+00
      IF(NC(I).NE.0)
     &  CALL FIAD(XQ(I)-XC(I),XFA(1,I),Q2,1,NC(I))
      IF(ND(I).NE.0)
     &  CALL FIAD(XQ(I)-XD(I),XFA(1,I),Q2,NC(I)+1,NCD(I))
      NDCT=NCD(I)
      IF(NDCT.EQ.0) GO TO 10
      DO 22 JT=1,NDCT,2
   22 XFA(JT+1,I)=-XFA(JT+1,I)
   10 CONTINUE
      CALL FYCLE
      ATUV=DEXP(-AW)*CCE(JC)*DCE(JD)*Q2
      IF(IUX) 300,301,302
  300 W=W+V*ATUV
      GO TO 6
  302 W=W+V*ATUV
  301 WDX=WDX+VDX*ATUV
    6 CONTINUE
      E=W
      WOCI=WDX
      RETURN
      END

C     ----------------------------------------------------------------
C     FINTC: split from ENTRY FINTC in original FINTA.
C     Sets IKI=1, then calls FINTA_CORE.
C     ----------------------------------------------------------------
      SUBROUTINE FINTC
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FPAR/V,XEA(5200),XFA(17,3),XF(17,3),XPA(3,400),
     &  SP(17,3),XQ(3),G1A(400),G2,AMF(400),MXDT,NAB(3),NCD(3),
     &  NT(3),NKT(3),NTXYM1
      COMMON /FIKI/ W,IKI

      IKI=1
      W=0.0D+00
      V=0.0D+00
      CALL FINTA_CORE
      RETURN
      END

C     ----------------------------------------------------------------
C     FINTD: split from ENTRY FINTD in original FINTA.
C     Sets IKI=2, then calls FINTA_CORE.
C     ----------------------------------------------------------------
      SUBROUTINE FINTD
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FPAR/V,XEA(5200),XFA(17,3),XF(17,3),XPA(3,400),
     &  SP(17,3),XQ(3),G1A(400),G2,AMF(400),MXDT,NAB(3),NCD(3),
     &  NT(3),NKT(3),NTXYM1
      COMMON /FIKI/ W,IKI

      IKI=2
      W=0.0D+00
      V=0.0D+00
      CALL FINTA_CORE
      RETURN
      END

C     ----------------------------------------------------------------
C     FINTE: split from ENTRY FINTE in original FINTA.
C     ----------------------------------------------------------------
      SUBROUTINE FINTE
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FFLA/XA(3),XB(3),XC(3),XD(3),AOE(20),ACE(20),BOE(20),
     &  BCE(20),COE(20),CCE(20),DOE(20),DCE(20),NA(3),NB(3),NC(3),
     &  ND(3),ITA,ITB,ITC,ITD
      COMMON/FFDG/XX,FNU(24),DFAC(24),GICA(81),NTSUM1
      COMMON/FNTV/E,ES,ET,ONA(100),IA,IPA,IB,IPB,IC,IPC,ID,IPD,
     &  NAP(100),ITIME
      COMMON/FPAR/V,XEA(5200),XFA(17,3),XF(17,3),XPA(3,400),
     &  SP(17,3),XQ(3),G1A(400),G2,AMF(400),MXDT,NAB(3),NCD(3),
     &  NT(3),NKT(3),NTXYM1
      EQUIVALENCE (NKT,NTM1),(NKT(2),NTM2),(NKT(3),NTM3)

      G2=2.0D+40
      DO 50 I=1,3
      SP(1,I)=1.0D+00
      NKT(I)=NAB(I)+NC(I)
      NT(I)=NKT(I)+1
      XQ(I)=XC(I)
      NCIT=NC(I)
      NCD(I)=NCIT
      DO 51 JT=1,NCIT
   51 XFA(JT,I)=0.0D+00
   50 XFA(NCIT+1,I)=1.0D+00
      NTXYM1=NTM1+NTM2
      NTSUM1=NTXYM1+NT(3)
      CALL FYCLE
      E=V*1.0306560687383531D+10**2
      RETURN
      END
