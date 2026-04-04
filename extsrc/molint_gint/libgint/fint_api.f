C     ================================================================
C     GINT_API: libcint-compatible Fortran API for the GINT library
C
C     This module provides a shell-based interface similar to libcint's
C     Cartesian GTO integral routines, wrapping the GINT integral
C     engine (F. Sasaki).
C
C     Data layout follows libcint conventions:
C
C       atm(1:6, 1:natm)  -- atom info
C         atm(1,i) = nuclear charge Z
C         atm(2,i) = offset in env for coordinates (0-based)
C         atm(3,i) = nuclear model (unused, set 0)
C         atm(4,i) = offset in env for zeta        (unused, set 0)
C         atm(5,i) = 0  (reserved)
C         atm(6,i) = 0  (reserved)
C
C       bas(1:8, 1:nbas)  -- shell (basis) info
C         bas(1,i) = atom index (1-based, unlike libcint's 0-based)
C         bas(2,i) = angular momentum l
C         bas(3,i) = number of primitive GTOs
C         bas(4,i) = number of contractions (must be 1 for now)
C         bas(5,i) = 0  (kappa, unused)
C         bas(6,i) = offset in env for exponents   (0-based)
C         bas(7,i) = offset in env for coefficients (0-based)
C         bas(8,i) = 0  (reserved)
C
C       env(*)  -- flat array holding coordinates, exponents, coeffs
C         env(atm(2,i)+1 : atm(2,i)+3) = x, y, z of atom i
C         env(bas(6,i)+1 : bas(6,i)+nprim) = exponents
C         env(bas(7,i)+1 : bas(7,i)+nprim) = contraction coefficients
C                                             (for normalized primitives)
C
C     NOTE: offsets in atm/bas are 0-based (libcint convention) but
C           Fortran arrays are 1-based, so env(offset+1) is used.
C           Atom indices in bas(1,i) are 1-based for Fortran naturalness.
C
C     ================================================================

C     ----------------------------------------------------------------
C     GINT_SHELL_MAP: internal mapping from shell to GINT AO indices
C
C       ISHAO(i)  = starting GINT AO index for shell i
C       NCART(i)  = number of Cartesian components in shell i
C       NSHELL    = total number of shells registered
C       NAOALL    = total number of AOs registered
C       NPSALL    = total number of centers registered
C       IATOM(i)  = atom index for center i (used for charge lookup)
C
C     Cartesian component ordering (same as libcint _cart):
C       l=0: s                                       (1 component)
C       l=1: px, py, pz                              (3 components)
C       l=2: dxx, dxy, dxz, dyy, dyz, dzz            (6 components)
C       l=3: fxxx,fxxy,fxxz,fxyy,fxyz,fxzz,
C            fyyy,fyyz,fyzz,fzzz                     (10 components)
C     General: (l+1)*(l+2)/2 components, ordered by
C       descending nx, then descending ny for each nx.
C     ----------------------------------------------------------------
      BLOCK DATA GINT_SHELL_INIT
      IMPLICIT NONE
      INTEGER ISHAO, NCART, NSHELL, NAOALL, NPSALL, IATOM
      COMMON/GSHMAP/ISHAO(100),NCART(100),NSHELL,NAOALL,NPSALL,
     &              IATOM(100)
      DATA NSHELL/0/, NAOALL/0/, NPSALL/0/
      END

C     ================================================================
C     GINT_INIT:  Initialize GINT from libcint-style arrays
C
C     This fills COMMON/FAREA/ (OE, CE, XSD, NQND, NAO, NPS) and
C     the shell-mapping table COMMON/GSHMAP/.
C
C     Arguments:
C       ATM(6, NATM)  -- atom table
C       NATM          -- number of atoms
C       BAS(8, NBAS)  -- shell table
C       NBAS          -- number of shells
C       ENV(*)        -- environment array (coords, exponents, coeffs)
C     ================================================================
      SUBROUTINE GINT_INIT(ATM, NATM, BAS, NBAS, ENV)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER ATM(6,*), NATM, BAS(8,*), NBAS
      REAL*8  ENV(*)

C     --- GINT internal common blocks ---
      COMMON/FAREA/OE(20,102),CE(20,102),XSD(3,100),NQND(4,102),NAO,
     &  NPS,LA,LPA,LB,LPB
      COMMON/GSHMAP/ISHAO(100),NCART(100),NSHELL,NAOALL,NPSALL,
     &              IATOM(100)

C     --- Local variables ---
      INTEGER NQ(3)
      INTEGER ISH, IAT, LANG, NPRIM, NCTR, IOFF_E, IOFF_C, IOFF_X
      INTEGER IAO, IPS, I, J, NX, NY, NZ
      INTEGER JP, JQ, KL, JIAO, IAO_L00
      REAL*8  CN, CNORM, APQ, GAUSL, DFAC_L
      REAL*8  PIPI32, K_CONST
      PARAMETER (PIPI32 = 5.568327996831708D0)
C     K_CONST = C3/(16*(pi/2)^(3/2)) normalizes FNFB^2 × GAUSL to true WS
C     K_CONST = C3 * 2^(3/2) / (16 * pi^(3/2)) = 0.119545343284187
      PARAMETER (K_CONST = 0.119545343284187D0)
      DOUBLE PRECISION FNFB
      EXTERNAL FNFB

C     --- 1. Set up centers (nuclear positions) ---
C     Copy atom coordinates into XSD.
C     Each atom becomes one center. Additional centers (for potentials
C     at non-nuclear positions) can be added later.
      NPS = NATM
      NPSALL = NATM
      DO 100 IAT = 1, NATM
        IOFF_X = ATM(2, IAT)
        XSD(1, IAT) = ENV(IOFF_X + 1)
        XSD(2, IAT) = ENV(IOFF_X + 2)
        XSD(3, IAT) = ENV(IOFF_X + 3)
        IATOM(IAT) = IAT
  100 CONTINUE

C     --- 2. Expand shells into individual Cartesian AOs ---
      IAO = 0
      NSHELL = NBAS

      DO 200 ISH = 1, NBAS
        IAT   = BAS(1, ISH)
        LANG  = BAS(2, ISH)
        NPRIM = BAS(3, ISH)
        NCTR  = BAS(4, ISH)
        IOFF_E = BAS(6, ISH)
        IOFF_C = BAS(7, ISH)

C       Sanity checks
        IF (NPRIM .GT. 20) THEN
          WRITE(6,*) 'GINT_INIT ERROR: NPRIM > 20 for shell', ISH
          STOP
        END IF
        IF (NCTR .NE. 1) THEN
          WRITE(6,*) 'GINT_INIT ERROR: NCTR must be 1, shell', ISH
          STOP
        END IF

C       Record shell->AO mapping
        ISHAO(ISH) = IAO + 1
        NCART(ISH) = (LANG + 1) * (LANG + 2) / 2

C       Loop over Cartesian components: nx+ny+nz = LANG
C       Ordering: descending nx, for each nx descending ny
        DO 210 NX = LANG, 0, -1
          DO 220 NY = LANG - NX, 0, -1
            NZ = LANG - NX - NY
            IAO = IAO + 1

            IF (IAO .GT. 102) THEN
              WRITE(6,*) 'GINT_INIT ERROR: too many AOs (>102)'
              STOP
            END IF

C           Set quantum numbers
            NQND(1, IAO) = NX
            NQND(2, IAO) = NY
            NQND(3, IAO) = NZ
            NQND(4, IAO) = NPRIM

C           Copy exponents and compute unnormalized coefficients
            NQ(1) = NX
            NQ(2) = NY
            NQ(3) = NZ
            DO 230 J = 1, NPRIM
              OE(J, IAO) = ENV(IOFF_E + J)
              CE(J, IAO) = ENV(IOFF_C + J) * FNFB(NQ, OE(J, IAO))
  230       CONTINUE

C           Zero out remaining slots
            DO 240 J = NPRIM + 1, 20
              OE(J, IAO) = 0.0D0
              CE(J, IAO) = 0.0D0
  240       CONTINUE
  220     CONTINUE
  210   CONTINUE

C       --- Contracted shell normalization ---
C       See comments in GINT_INIT preamble for the formula.
      IAO_L00 = ISHAO(ISH)
      CN = 0.0D0
      DO 250 JP = 1, NPRIM
        DO 251 JQ = 1, NPRIM
          APQ = OE(JP, IAO_L00) + OE(JQ, IAO_L00)
C         GAUSL = (2*LANG-1)!! * pi^(3/2) / ((2*APQ)^LANG * APQ^(3/2))
C         (2*LANG-1)!! : 1,1,3,15,105 for LANG=0,1,2,3,4
          IF (LANG .EQ. 0) DFAC_L = 1.0D0
          IF (LANG .EQ. 1) DFAC_L = 1.0D0
          IF (LANG .EQ. 2) DFAC_L = 3.0D0
          IF (LANG .EQ. 3) DFAC_L = 15.0D0
          IF (LANG .EQ. 4) DFAC_L = 105.0D0
          GAUSL = PIPI32 * DFAC_L
          GAUSL = GAUSL / (APQ * DSQRT(APQ))
          DO 252 KL = 1, LANG
            GAUSL = GAUSL / (2.0D0 * APQ)
  252     CONTINUE
          CN = CN + CE(JP, IAO_L00) * CE(JQ, IAO_L00) * GAUSL
  251   CONTINUE
  250 CONTINUE
      CN = CN * K_CONST
      CNORM = 1.0D0 / DSQRT(CN)
C     Apply normalization to all Cartesian AOs in this shell
      DO 260 JIAO = ISHAO(ISH), IAO
        DO 261 J = 1, NPRIM
          CE(J, JIAO) = CE(J, JIAO) * CNORM
  261   CONTINUE
  260 CONTINUE

  200 CONTINUE

C     --- 3. Set global counts ---
      NAO = IAO
      NAOALL = IAO

C     --- 4. Invalidate GINT's internal cache ---
C     LA, LPA, LB, LPB are used by FVALUE to skip redundant FTRA
C     calls. Set them to impossible values to force recomputation.
      LA  = -9
      LPA = -9
      LB  = -9
      LPB = -9

C     (Removed: print statement suppressed for repeated internal calls)
C     To check init results, read NAOALL, NPSALL from COMMON/GSHMAP/

      RETURN
      END

C     ================================================================
C     GINT_CGTO_CART:  Number of Cartesian GTOs in a shell
C
C     Arguments:
C       SHELL_ID  -- shell index (1-based)
C       BAS(8,*)  -- shell table
C     Returns:
C       (l+1)*(l+2)/2
C     ================================================================
      INTEGER FUNCTION GINT_CGTO_CART(SHELL_ID, BAS)
      IMPLICIT NONE
      INTEGER SHELL_ID, BAS(8,*)
      INTEGER L
      L = BAS(2, SHELL_ID)
      GINT_CGTO_CART = (L + 1) * (L + 2) / 2
      RETURN
      END

C     ================================================================
C     GINT_SHELL_TO_AO: Get the starting GINT AO index for a shell
C
C     Utility for users who need direct access to GINT internals.
C     ================================================================
      INTEGER FUNCTION GINT_SHELL_TO_AO(SHELL_ID)
      IMPLICIT NONE
      INTEGER SHELL_ID
      INTEGER ISHAO, NCART, NSHELL, NAOALL, NPSALL, IATOM
      COMMON/GSHMAP/ISHAO(100),NCART(100),NSHELL,NAOALL,NPSALL,
     &              IATOM(100)
      GINT_SHELL_TO_AO = ISHAO(SHELL_ID)
      RETURN
      END

C     ================================================================
C     GINT_1E_STV: Compute S, T, V together for a shell pair
C
C     FVALUE computes overlap (WS), kinetic (WT), and nuclear
C     attraction (ONA) in a single call. This routine calls FVALUE
C     exactly once per Cartesian component pair, avoiding the
C     state-corruption issue that occurs with repeated FVALUE calls.
C
C     Arguments:
C       SBUF(di*dj)  -- overlap integrals (output)
C       TBUF(di*dj)  -- kinetic energy integrals (output)
C       VBUF(di*dj)  -- nuclear attraction integrals (output)
C                        V = sum_A (-Z_A) <i|1/r_A|j>
C       SHLS(2)      -- shell indices (1-based)
C       ATM,NATM,BAS,NBAS,ENV -- standard libcint-style arrays
C
C     Any of SBUF, TBUF, VBUF may point to the same array if only
C     one component is needed (the others will be overwritten).
C     ================================================================
      SUBROUTINE GINT_1E_STV(SBUF, TBUF, VBUF,
     &                        SHLS, ATM, NATM, BAS, NBAS, ENV)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  SBUF(*), TBUF(*), VBUF(*)
      INTEGER SHLS(2), ATM(6,*), NATM, BAS(8,*), NBAS
      REAL*8  ENV(*)

      COMMON/FNTV/W,WS,WT,ONA(100),IA,IPA,IB,IPB,IC,IPC,ID,IPD,
     &  NAP(100),ITIME
      INTEGER NQND, NAO, NPS
      INTEGER LA, LPA, LB, LPB
      COMMON/FAREA/OE(20,102),CE(20,102),XSD(3,100),
     &  NQND(4,102),NAO,NPS,
     &  LA,LPA,LB,LPB
      COMMON/GSHMAP/ISHAO(100),NCART(100),NSHELL,NAOALL,NPSALL,
     &              IATOM(100)

      INTEGER ISH, JSH, DI, DJ, II, JJ, IAO0, JAO0, ICNTR, JCNTR
      INTEGER KK, IDX
      REAL*8  VTOT
      INTEGER LI, LJ
      REAL*8  SCALE, GINT_SCALE
      EXTERNAL GINT_SCALE

      ISH = SHLS(1)
      JSH = SHLS(2)
      DI  = NCART(ISH)
      DJ  = NCART(JSH)
      IAO0  = ISHAO(ISH) - 1
      JAO0  = ISHAO(JSH) - 1
      ICNTR = BAS(1, ISH)
      JCNTR = BAS(1, JSH)

C     --- 角運動量からスケール係数を計算 ---
      LI = BAS(2, ISH)
      LJ = BAS(2, JSH)
      SCALE = GINT_SCALE(LI) * GINT_SCALE(LJ)

      LA  = -9
      LPA = -9
      LB  = -9
      LPB = -9
      IC  = NATM
      IPC = 0
      ID  = 0
      IPD = 0
      DO 5 KK = 1, NATM
        NAP(KK) = KK
    5 CONTINUE

      DO 20 JJ = 1, DJ
        IB  = JAO0 + JJ
        IPB = JCNTR
        DO 10 II = 1, DI
          IA  = IAO0 + II
          IPA = ICNTR

          CALL FVALUE

          IDX = (JJ - 1) * DI + II
C         --- スケール係数を適用 ---
          SBUF(IDX) = WS * SCALE
          TBUF(IDX) = WT * SCALE

          VTOT = 0.0D0
          DO 15 KK = 1, NATM
            VTOT = VTOT - DBLE(ATM(1, KK)) * ONA(KK)
   15     CONTINUE
          VBUF(IDX) = VTOT * SCALE
   10   CONTINUE
   20 CONTINUE

      RETURN
      END

C     ================================================================
C     GINT_1E_OVL:  Overlap integrals <i|j> for a shell pair
C     Thin wrapper around GINT_1E_STV.
C     ================================================================
      SUBROUTINE GINT_1E_OVL(BUF, SHLS, ATM, NATM, BAS, NBAS, ENV)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  BUF(*)
      INTEGER SHLS(2), ATM(6,*), NATM, BAS(8,*), NBAS
      REAL*8  ENV(*)
C     Scratch buffers (max shell size = 10f x 10f = 100)
      REAL*8  TDUM(100), VDUM(100)
      CALL GINT_1E_STV(BUF, TDUM, VDUM,
     &                  SHLS, ATM, NATM, BAS, NBAS, ENV)
      RETURN
      END

C     ================================================================
C     GINT_1E_KIN:  Kinetic energy integrals <i|-1/2 nabla^2|j>
C     Thin wrapper around GINT_1E_STV.
C     ================================================================
      SUBROUTINE GINT_1E_KIN(BUF, SHLS, ATM, NATM, BAS, NBAS, ENV)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  BUF(*)
      INTEGER SHLS(2), ATM(6,*), NATM, BAS(8,*), NBAS
      REAL*8  ENV(*)
      REAL*8  SDUM(100), VDUM(100)
      CALL GINT_1E_STV(SDUM, BUF, VDUM,
     &                  SHLS, ATM, NATM, BAS, NBAS, ENV)
      RETURN
      END

C     ================================================================
C     GINT_1E_NUC:  Nuclear attraction integrals <i|V_nuc|j>
C     Thin wrapper around GINT_1E_STV.
C     ================================================================
      SUBROUTINE GINT_1E_NUC(BUF, SHLS, ATM, NATM, BAS, NBAS, ENV)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  BUF(*)
      INTEGER SHLS(2), ATM(6,*), NATM, BAS(8,*), NBAS
      REAL*8  ENV(*)
      REAL*8  SDUM(100), TDUM(100)
      CALL GINT_1E_STV(SDUM, TDUM, BUF,
     &                  SHLS, ATM, NATM, BAS, NBAS, ENV)
      RETURN
      END

C     ================================================================
C     GINT_1E_VNUC_RAW:  Raw nuclear attraction <i|1/|r-R_C||j>
C
C     Unlike GINT_1E_NUC, this returns the bare potential at each
C     center separately, without multiplying by -Z.
C
C     BUF(di, dj, natm) -- output: BUF(a, b, C) = <a|1/|r-R_C||b>
C
C     This is analogous to libcint's cint1e_rinv_cart when the
C     rinv origin is set to each nucleus in turn.
C     ================================================================
      SUBROUTINE GINT_1E_VNUC_RAW(BUF, SHLS, ATM, NATM,
     &                             BAS, NBAS, ENV)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  BUF(*)
      INTEGER SHLS(2), ATM(6,*), NATM, BAS(8,*), NBAS
      REAL*8  ENV(*)

      COMMON/FNTV/W,WS,WT,ONA(100),IA,IPA,IB,IPB,IC,IPC,ID,IPD,
     &  NAP(100),ITIME
      COMMON/GSHMAP/ISHAO(100),NCART(100),NSHELL,NAOALL,NPSALL,
     &              IATOM(100)

      INTEGER ISH, JSH, DI, DJ, II, JJ, IAO0, JAO0, ICNTR, JCNTR
      INTEGER KK, DIDJ
      INTEGER LI, LJ
      REAL*8  SCALE, GINT_SCALE
      EXTERNAL GINT_SCALE

      ISH = SHLS(1)
      JSH = SHLS(2)
      DI  = NCART(ISH)
      DJ  = NCART(JSH)
      IAO0  = ISHAO(ISH) - 1
      JAO0  = ISHAO(JSH) - 1
      ICNTR = BAS(1, ISH)
      JCNTR = BAS(1, JSH)

C     --- 角運動量からスケール係数を計算 ---
      LI = BAS(2, ISH)
      LJ = BAS(2, JSH)
      SCALE = GINT_SCALE(LI) * GINT_SCALE(LJ)

      CALL GINT_FLUSH_CACHE
      DIDJ  = DI * DJ

      IC  = NATM
      IPC = 0
      ID  = 0
      IPD = 0
      DO 5 KK = 1, NATM
        NAP(KK) = KK
    5 CONTINUE

      DO 20 JJ = 1, DJ
        IB  = JAO0 + JJ
        IPB = JCNTR
        DO 10 II = 1, DI
          IA  = IAO0 + II
          IPA = ICNTR

          CALL FVALUE

          DO 15 KK = 1, NATM
C           --- スケール係数を適用 ---
            BUF((KK - 1) * DIDJ + (JJ - 1) * DI + II) = ONA(KK) * SCALE
   15     CONTINUE
   10   CONTINUE
   20 CONTINUE

      RETURN
      END

C     ================================================================
C     GINT_1E_MOMENT:  Dipole/quadrupole moment integrals
C
C     Computes moment integrals around a specified center:
C       WM(1:3)  = <i| (x-Xc), (y-Yc), (z-Zc) |j>
C       WN(1:6)  = <i| (x-Xc)^2, (x-Xc)(y-Yc), ... |j>
C
C     Arguments:
C       BUFM(di, dj, 3)  -- dipole moment integrals (output)
C       BUFQ(di, dj, 6)  -- quadrupole moment integrals (output)
C                            Set to NULL/unused if IORD=1
C       SHLS(2)           -- shell indices (1-based)
C       ICENTER           -- center index for moments (1-based)
C       IORD              -- 1 = dipole only, 2 = dipole + quadrupole
C       ATM,NATM,BAS,NBAS,ENV -- standard arrays
C     ================================================================
      SUBROUTINE GINT_1E_MOMENT(BUFM, BUFQ, SHLS, ICENTER, IORD,
     &                           ATM, NATM, BAS, NBAS, ENV)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  BUFM(*), BUFQ(*)
      INTEGER SHLS(2), ICENTER, IORD
      INTEGER ATM(6,*), NATM, BAS(8,*), NBAS
      REAL*8  ENV(*)

      COMMON/FNTV/W,WS,WT,ONA(100),IA,IPA,IB,IPB,IC,IPC,ID,IPD,
     &  NAP(100),ITIME
      COMMON/FNTVS/WM(3),WN(6)
      COMMON/GSHMAP/ISHAO(100),NCART(100),NSHELL,NAOALL,NPSALL,
     &              IATOM(100)

      INTEGER ISH, JSH, DI, DJ, II, JJ, IAO0, JAO0, ICNTR, JCNTR
      INTEGER KK, DIDJ, IDX

      ISH = SHLS(1)
      JSH = SHLS(2)
      DI  = NCART(ISH)
      DJ  = NCART(JSH)
      IAO0  = ISHAO(ISH) - 1
      JAO0  = ISHAO(JSH) - 1
      ICNTR = BAS(1, ISH)
      JCNTR = BAS(1, JSH)
      CALL GINT_FLUSH_CACHE
      DIDJ  = DI * DJ

C     Set 1e mode: compute V (for state consistency) + moments
      IC  = NATM
      IPC = 0
      ID  = 0
      IPD = ICENTER
      DO 5 KK = 1, NATM
        NAP(KK) = KK
    5 CONTINUE

      DO 20 JJ = 1, DJ
        IB  = JAO0 + JJ
        IPB = JCNTR
        DO 10 II = 1, DI
          IA  = IAO0 + II
          IPA = ICNTR

          CALL FVALUE

          IDX = (JJ - 1) * DI + II
C         Dipole: x, y, z
          DO 12 KK = 1, 3
            BUFM((KK - 1) * DIDJ + IDX) = WM(KK)
   12     CONTINUE
C         Quadrupole: xx, xy, xz, yy, yz, zz
          IF (IORD .GE. 2) THEN
            DO 14 KK = 1, 6
              BUFQ((KK - 1) * DIDJ + IDX) = WN(KK)
   14       CONTINUE
          END IF
   10   CONTINUE
   20 CONTINUE

      RETURN
      END

C     ================================================================
C     GINT_1E_FIELD:  Electric field integrals  dV/dx, dV/dy, dV/dz
C
C     BUF(di, dj, 3, natm) -- output
C       BUF(:,:,k,A) = d/d{x_k} <i|1/|r-R_A||j>
C
C     Same shell-pair interface.
C     ================================================================
      SUBROUTINE GINT_1E_FIELD(BUF, SHLS, ATM, NATM,
     &                          BAS, NBAS, ENV)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  BUF(*)
      INTEGER SHLS(2), ATM(6,*), NATM, BAS(8,*), NBAS
      REAL*8  ENV(*)

      COMMON/FNTV/W,WS,WT,ONA(100),IA,IPA,IB,IPB,IC,IPC,ID_COM,IPD_COM,
     &  NAP(100),ITIME
      COMMON/FNTVT/ONAF(3,100),ONAFG(6,100)
      COMMON/GSHMAP/ISHAO(100),NCART(100),NSHELL,NAOALL,NPSALL,
     &              IATOM(100)

      INTEGER ISH, JSH, DI, DJ, II, JJ, IAO0, JAO0, ICNTR, JCNTR
      INTEGER KK, LL, DIDJ, DIDJ3
      INTEGER LI, LJ
      REAL*8  SCALE, GINT_SCALE
      EXTERNAL GINT_SCALE

      ISH = SHLS(1)
      JSH = SHLS(2)
      DI  = NCART(ISH)
      DJ  = NCART(JSH)
      IAO0  = ISHAO(ISH) - 1
      JAO0  = ISHAO(JSH) - 1
      ICNTR = BAS(1, ISH)
      JCNTR = BAS(1, JSH)

      LI = BAS(2, ISH)
      LJ = BAS(2, JSH)
      SCALE = GINT_SCALE(LI) * GINT_SCALE(LJ)

      CALL GINT_FLUSH_CACHE
      DIDJ  = DI * DJ
      DIDJ3 = DIDJ * 3

      IC  = NATM
      IPC = 0
C     --- 修正点: 共通ブロックのモード設定を ID_COM = 1 にする ---
      ID_COM  = 1
      IPD_COM = 0
      DO 5 KK = 1, NATM
        NAP(KK) = KK
    5 CONTINUE

      DO 20 JJ = 1, DJ
        IB  = JAO0 + JJ
        IPB = JCNTR
        DO 10 II = 1, DI
          IA  = IAO0 + II
          IPA = ICNTR
          CALL FVALUE
          DO 15 KK = 1, NATM
C           --- 核電荷 -Z を掛けてスケールを適用 ---
            DO 16 LL = 1, 3
              BUF((KK-1)*DIDJ3 + (LL-1)*DIDJ + (JJ-1)*DI + II)
     &          = -DBLE(ATM(1, KK)) * ONAF(LL, KK) * SCALE
   16       CONTINUE
   15     CONTINUE
   10   CONTINUE
   20 CONTINUE
      RETURN
      END

C     ================================================================
C     GINT_1E_FLDGRD:  Electric field gradient integrals
C
C     BUF(di, dj, 6, natm) -- output
C       BUF(:,:,k,A) = d2V/d{comp_k} at R_A
C       comp order: xx, xy, xz, yy, yz, zz
C     ================================================================
      SUBROUTINE GINT_1E_FLDGRD(BUF, SHLS, ATM, NATM,
     &                         BAS, NBAS, ENV)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  BUF(*)
      INTEGER SHLS(2), ATM(6,*), NATM, BAS(8,*), NBAS
      REAL*8  ENV(*)

      COMMON/FNTV/W,WS,WT,ONA(100),IA,IPA,IB,IPB,IC,IPC,ID_COM,IPD_COM,
     &  NAP(100),ITIME
      COMMON/FNTVT/ONAF(3,100),ONAFG(6,100)
      COMMON/GSHMAP/ISHAO(100),NCART(100),NSHELL,NAOALL,NPSALL,
     &              IATOM(100)

      INTEGER ISH, JSH, DI, DJ, II, JJ, IAO0, JAO0, ICNTR, JCNTR
      INTEGER KK, LL, DIDJ, DIDJ6
      INTEGER LI, LJ
      REAL*8  SCALE, GINT_SCALE
      EXTERNAL GINT_SCALE

      ISH = SHLS(1)
      JSH = SHLS(2)
      DI  = NCART(ISH)
      DJ  = NCART(JSH)
      IAO0  = ISHAO(ISH) - 1
      JAO0  = ISHAO(JSH) - 1
      ICNTR = BAS(1, ISH)
      JCNTR = BAS(1, JSH)

      LI = BAS(2, ISH)
      LJ = BAS(2, JSH)
      SCALE = GINT_SCALE(LI) * GINT_SCALE(LJ)

      CALL GINT_FLUSH_CACHE
      DIDJ  = DI * DJ
      DIDJ6 = DIDJ * 6

      IC  = NATM
      IPC = 0
C     --- 修正点: 共通ブロックのモード設定を ID_COM = 2 にする ---
      ID_COM  = 2
      IPD_COM = 0
      DO 5 KK = 1, NATM
        NAP(KK) = KK
    5 CONTINUE

      DO 20 JJ = 1, DJ
        IB  = JAO0 + JJ
        IPB = JCNTR
        DO 10 II = 1, DI
          IA  = IAO0 + II
          IPA = ICNTR
          CALL FVALUE
          DO 15 KK = 1, NATM
C           --- 核電荷 -Z を掛けてスケールを適用 ---
            DO 16 LL = 1, 6
              BUF((KK-1)*DIDJ6 + (LL-1)*DIDJ + (JJ-1)*DI + II)
     &          = -DBLE(ATM(1, KK)) * ONAFG(LL, KK) * SCALE
   16       CONTINUE
   15     CONTINUE
   10   CONTINUE
   20 CONTINUE
      RETURN
      END

C     ================================================================
C     GINT_FLUSH_CACHE: Force FVALUE to rerun FTRA on next call
C     ================================================================
      SUBROUTINE GINT_FLUSH_CACHE
      IMPLICIT NONE
      REAL*8  OE, CE, XSD
      INTEGER NQND, NAO, NPS, LA, LPA, LB, LPB
      COMMON/FAREA/OE(20,102),CE(20,102),XSD(3,100),NQND(4,102),NAO,
     &  NPS,LA,LPA,LB,LPB
      LA  = -9
      LPA = -9
      LB  = -9
      LPB = -9
      RETURN
      END

C     ================================================================
C     GINT_SCALE: PySCF-compatible Cartesian GTO normalization scale
C     Returns sqrt(4*pi / (2l+1)) for l >= 2, and 1.0 for l <= 1.
C     ================================================================
      DOUBLE PRECISION FUNCTION GINT_SCALE(L)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER L
      REAL*8 FOURPI
      PARAMETER (FOURPI=12.56637061435917295D0)
      
      IF (L .LE. 1) THEN
        GINT_SCALE = 1.0D0
      ELSE
        GINT_SCALE = DSQRT(FOURPI / DFLOAT(2*L + 1))
      END IF
      RETURN
      END

C     ================================================================
C     GINT_2E: Two-electron repulsion integrals (ij|kl) for a shell
C              quartet. (PySCF normalization scaling applied)
C     ================================================================
      SUBROUTINE GINT_2E(BUF, SHLS, ATM, NATM, BAS, NBAS, ENV)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  BUF(*)
      INTEGER SHLS(4), ATM(6,*), NATM, BAS(8,*), NBAS
      REAL*8  ENV(*)

      COMMON/FNTV/W,WS,WT,ONA(100),IA,IPA,IB,IPB,IC,IPC,ID,IPD,
     &  NAP(100),ITIME
      COMMON/GSHMAP/ISHAO(100),NCART(100),NSHELL,NAOALL,NPSALL,
     &              IATOM(100)
      DOUBLE PRECISION VDX
      COMMON/FOCI/VDX,WOCI,IUX

      INTEGER ISH,JSH,KSH,LSH
      INTEGER DI,DJ,DK,DL
      INTEGER IAO0,JAO0,KAO0,LAO0
      INTEGER ICNTR,JCNTR,KCNTR,LCNTR
      INTEGER II,JJ,KK,LL,IDX
      INTEGER LI,LJ,LK,LL_ANG
      REAL*8  SCALE, GINT_SCALE
      EXTERNAL GINT_SCALE

      ISH = SHLS(1)
      JSH = SHLS(2)
      KSH = SHLS(3)
      LSH = SHLS(4)

      DI = NCART(ISH)
      DJ = NCART(JSH)
      DK = NCART(KSH)
      DL = NCART(LSH)

      IAO0  = ISHAO(ISH) - 1
      JAO0  = ISHAO(JSH) - 1
      KAO0  = ISHAO(KSH) - 1
      LAO0  = ISHAO(LSH) - 1

      ICNTR = BAS(1,ISH)
      JCNTR = BAS(1,JSH)
      KCNTR = BAS(1,KSH)
      LCNTR = BAS(1,LSH)

C     --- 各軌道の角運動量からスケール係数を計算 ---
      LI = BAS(2,ISH)
      LJ = BAS(2,JSH)
      LK = BAS(2,KSH)
      LL_ANG = BAS(2,LSH)

      SCALE = GINT_SCALE(LI) * GINT_SCALE(LJ) *
     &        GINT_SCALE(LK) * GINT_SCALE(LL_ANG)

      IUX = -1

      DO 40 LL = 1, DL
        ID  = LAO0 + LL
        IPD = LCNTR
        DO 30 KK = 1, DK
          IC  = KAO0 + KK
          IPC = KCNTR
          CALL GINT_FLUSH_CACHE
          DO 20 JJ = 1, DJ
            IB  = JAO0 + JJ
            IPB = JCNTR
            CALL GINT_FLUSH_CACHE
            DO 10 II = 1, DI
              IA  = IAO0 + II
              IPA = ICNTR

              CALL FVALUE

              IDX = (LL-1)*DI*DJ*DK + (KK-1)*DI*DJ + (JJ-1)*DI + II
C             --- スケール係数を適用 ---
              BUF(IDX) = W * SCALE

   10       CONTINUE
   20     CONTINUE
   30   CONTINUE
   40 CONTINUE

      RETURN
      END

C     ================================================================
C     GINT_2E_SYMM: Two-electron repulsion integrals with intra-shell
C                   permutation symmetry. (PySCF scaling applied)
C     ================================================================
      SUBROUTINE GINT_2E_SYMM(BUF, SHLS, ATM, NATM, BAS, NBAS, ENV)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  BUF(*)
      INTEGER SHLS(4), ATM(6,*), NATM, BAS(8,*), NBAS
      REAL*8  ENV(*)

      COMMON/FNTV/W,WS,WT,ONA(100),IA,IPA,IB,IPB,IC,IPC,ID,IPD,
     &  NAP(100),ITIME
      COMMON/FAREA/OE(20,102),CE(20,102),XSD(3,100),
     &  NQND(4,102),NAO,NPS,
     &  LA,LPA,LB,LPB
      COMMON/GSHMAP/ISHAO(100),NCART(100),NSHELL,NAOALL,NPSALL,
     &              IATOM(100)
      DOUBLE PRECISION VDX
      COMMON/FOCI/VDX,WOCI,IUX

      INTEGER ISH,JSH,KSH,LSH
      INTEGER DI,DJ,DK,DL,DIDJ,DKDL
      INTEGER IAO0,JAO0,KAO0,LAO0
      INTEGER ICNTR,JCNTR,KCNTR,LCNTR
      INTEGER II,JJ,KK,LL
      INTEGER IDX,IDX_SYM
      LOGICAL SYM_IJ,SYM_KL,SYM_BRAKET
      REAL*8  WVAL
      INTEGER LI,LJ,LK,LL_ANG
      REAL*8  SCALE, GINT_SCALE
      EXTERNAL GINT_SCALE

      ISH = SHLS(1)
      JSH = SHLS(2)
      KSH = SHLS(3)
      LSH = SHLS(4)

      DI = NCART(ISH)
      DJ = NCART(JSH)
      DK = NCART(KSH)
      DL = NCART(LSH)
      DIDJ = DI * DJ
      DKDL = DK * DL

      IAO0  = ISHAO(ISH) - 1
      JAO0  = ISHAO(JSH) - 1
      KAO0  = ISHAO(KSH) - 1
      LAO0  = ISHAO(LSH) - 1

      ICNTR = BAS(1,ISH)
      JCNTR = BAS(1,JSH)
      KCNTR = BAS(1,KSH)
      LCNTR = BAS(1,LSH)

C     --- 各軌道の角運動量からスケール係数を計算 ---
      LI = BAS(2,ISH)
      LJ = BAS(2,JSH)
      LK = BAS(2,KSH)
      LL_ANG = BAS(2,LSH)

      SCALE = GINT_SCALE(LI) * GINT_SCALE(LJ) *
     &        GINT_SCALE(LK) * GINT_SCALE(LL_ANG)

      IUX = -1

      SYM_IJ     = (ISH .EQ. JSH)
      SYM_KL     = (KSH .EQ. LSH)
      SYM_BRAKET = (ISH .EQ. KSH) .AND. (JSH .EQ. LSH)

      IF (.NOT.SYM_IJ .AND. .NOT.SYM_KL .AND. .NOT.SYM_BRAKET) THEN
        CALL GINT_2E(BUF, SHLS, ATM, NATM, BAS, NBAS, ENV)
        RETURN
      END IF

      DO 5 IDX = 1, DI*DJ*DK*DL
    5   BUF(IDX) = 0.0D0

      DO 40 LL = 1, DL
        ID  = LAO0 + LL
        IPD = LCNTR
        DO 30 KK = 1, DK
          IC  = KAO0 + KK
          IPC = KCNTR

          IF (SYM_KL .AND. KK.GT.LL) GO TO 30

          CALL GINT_FLUSH_CACHE
          DO 20 JJ = 1, DJ
            IB  = JAO0 + JJ
            IPB = JCNTR
            CALL GINT_FLUSH_CACHE
            DO 10 II = 1, DI
              IA  = IAO0 + II
              IPA = ICNTR

              IF (SYM_IJ .AND. II.GT.JJ) GO TO 10

              IF (SYM_BRAKET) THEN
                IF ((KK-1)*DL+LL .GT. (II-1)*DJ+JJ) GO TO 10
              END IF

              CALL FVALUE
C             --- スケール係数を適用 ---
              WVAL = W * SCALE

              IDX = (LL-1)*DIDJ*DK + (KK-1)*DIDJ + (JJ-1)*DI + II
              BUF(IDX) = WVAL

              IF (SYM_IJ .AND. II.NE.JJ) THEN
                IDX_SYM = (LL-1)*DIDJ*DK+(KK-1)*DIDJ+(II-1)*DI+JJ
                BUF(IDX_SYM) = WVAL
              END IF

              IF (SYM_KL .AND. KK.NE.LL) THEN
                IDX_SYM = (KK-1)*DIDJ*DL+(LL-1)*DIDJ+(JJ-1)*DI+II
                BUF(IDX_SYM) = WVAL
              END IF

              IF (SYM_IJ .AND. SYM_KL .AND.
     &            II.NE.JJ .AND. KK.NE.LL) THEN
                IDX_SYM = (KK-1)*DIDJ*DL+(LL-1)*DIDJ+(II-1)*DI+JJ
                BUF(IDX_SYM) = WVAL
              END IF

              IF (SYM_BRAKET) THEN
                IF ((KK-1)*DL+LL .NE. (II-1)*DJ+JJ) THEN
                  IDX_SYM = (JJ-1)*DKDL*DK+(II-1)*DKDL+(LL-1)*DK+KK
                  BUF(IDX_SYM) = WVAL
                  IF (SYM_IJ .AND. II.NE.JJ) THEN
                    IDX_SYM=(II-1)*DKDL*DK+(JJ-1)*DKDL+(LL-1)*DK+KK
                    BUF(IDX_SYM) = WVAL
                  END IF
                  IF (SYM_KL .AND. KK.NE.LL) THEN
                    IDX_SYM=(JJ-1)*DKDL*DL+(II-1)*DKDL+(KK-1)*DK+LL
                    BUF(IDX_SYM) = WVAL
                  END IF
                  IF (SYM_IJ .AND. SYM_KL .AND.
     &                II.NE.JJ .AND. KK.NE.LL) THEN
                    IDX_SYM=(II-1)*DKDL*DL+(JJ-1)*DKDL+(KK-1)*DK+LL
                    BUF(IDX_SYM) = WVAL
                  END IF
                END IF
              END IF

   10       CONTINUE
   20     CONTINUE
   30   CONTINUE
   40 CONTINUE

      RETURN
      END
C     ================================================================
C     GINT_CGTO_SPH:  Number of Spherical GTOs in a shell  (= 2l+1)
C     ================================================================
      INTEGER FUNCTION GINT_CGTO_SPH(SHELL_ID, BAS)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER SHELL_ID, BAS(8,*)
      GINT_CGTO_SPH = 2 * BAS(2, SHELL_ID) + 1
      RETURN
      END

C     ================================================================
C     GINT_C2S_1E:  Apply Cart->Sph transformation to one 1e block
C
C     Computes  SBUF(m,j) = sum_c T_l[m,c] * CBUF(c,j)
C
C     Arguments:
C       L     -- angular momentum of the shell being transformed
C       CBUF  -- Cartesian integrals (DC x DJ), column-major
C       DC    -- (L+1)(L+2)/2
C       DJ    -- number of ket components
C       SBUF  -- Spherical output   (DS x DJ), column-major
C       DS    -- 2L+1
C
C     Cart order (GINT): descending nx, for each nx descending ny
C       l=2: xx,xy,xz,yy,yz,zz
C       l=3: xxx,xxy,xxz,xyy,xyz,xzz,yyy,yyz,yzz,zzz
C       l=4: xxxx,xxxy,xxxz,xxyy,xxyz,xxzz,xyyy,xyyz,xyzz,xzzz,
C            yyyy,yyyz,yyzz,yzzz,zzzz
C     Sph order: m = -l, -l+1, ..., +l
C
C     Matrices satisfy T @ S_gram @ T^T = I (Gram-orthonormal).
C     Verified: max off-diagonal < 3e-16 for all l = 0..4.
C     ================================================================
      SUBROUTINE GINT_C2S_1E(L, CBUF, DC, DJ, SBUF, DS)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER L, DC, DJ, DS
      REAL*8  CBUF(DC, DJ), SBUF(DS, DJ)

C     ---- l=2: T2(5,6) column-major ----
      REAL*8 T2(5,6)
      DATA T2/
     & +0.000000000000000D0,+0.000000000000000D0,-0.500000000000000D0,
     & +0.000000000000000D0,+0.866025403784438D0,+1.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+1.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,-0.500000000000000D0,
     & +0.000000000000000D0,-0.866025403784438D0,+0.000000000000000D0,
     & +1.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +1.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0
     & /

C     ---- l=3: T3(7,10) column-major ----
      REAL*8 T3(7,10)
      DATA T3/
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,-0.612372435695794D0,+0.000000000000000D0,
     & +0.790569415042095D0,+1.060660171779821D0,+0.000000000000000D0,
     & -0.273861278752583D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,-0.670820393249937D0,
     & +0.000000000000000D0,+0.866025403784438D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,-0.273861278752583D0,+0.000000000000000D0,
     & -1.060660171779821D0,+0.000000000000000D0,+1.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +1.095445115010332D0,+0.000000000000000D0,+0.000000000000000D0,
     & -0.790569415042095D0,+0.000000000000000D0,-0.612372435695794D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,-0.670820393249937D0,+0.000000000000000D0,
     & -0.866025403784438D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+1.095445115010332D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +1.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0
     & /

C     ---- l=4: T4(9,15) column-major ----
      REAL*8 T4(9,15)
      DATA T4/
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.375000000000000D0,+0.000000000000000D0,
     & -0.559016994374948D0,+0.000000000000000D0,+0.739509972887452D0,
     & +1.118033988749895D0,+0.000000000000000D0,-0.422577127364258D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,-0.896421457000795D0,
     & +0.000000000000000D0,+0.790569415042095D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.219577516413420D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,-1.299038105676658D0,
     & +0.000000000000000D0,+1.060660171779821D0,+0.000000000000000D0,
     & -0.400891862868637D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,-0.878310065653680D0,+0.000000000000000D0,
     & +0.981980506061966D0,+0.000000000000000D0,+0.000000000000000D0,
     & -1.118033988749895D0,+0.000000000000000D0,-0.422577127364258D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,-0.400891862868637D0,
     & +0.000000000000000D0,-1.060660171779821D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+1.133893419027682D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+1.195228609334394D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.375000000000000D0,+0.000000000000000D0,
     & +0.559016994374948D0,+0.000000000000000D0,+0.739509972887452D0,
     & +0.000000000000000D0,-0.790569415042095D0,+0.000000000000000D0,
     & -0.896421457000795D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,-0.878310065653680D0,+0.000000000000000D0,
     & -0.981980506061966D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +1.195228609334394D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+1.000000000000000D0,+0.000000000000000D0,
     & +0.000000000000000D0,+0.000000000000000D0,+0.000000000000000D0
     & /

      INTEGER MI, CI, JJ

C     --- L=0: trivial ---
      IF (L .EQ. 0) THEN
        DO 10 JJ = 1, DJ
          SBUF(1, JJ) = CBUF(1, JJ)
   10   CONTINUE
        RETURN
      END IF

C     --- L=1: permutation px,py,pz -> m=-1(py), m=0(pz), m=+1(px) ---
      IF (L .EQ. 1) THEN
        DO 20 JJ = 1, DJ
          SBUF(1, JJ) = CBUF(2, JJ)
          SBUF(2, JJ) = CBUF(3, JJ)
          SBUF(3, JJ) = CBUF(1, JJ)
   20   CONTINUE
        RETURN
      END IF

C     --- L=2 ---
      IF (L .EQ. 2) THEN
        DO 30 JJ = 1, DJ
          DO 31 MI = 1, 5
            SBUF(MI, JJ) = 0.0D0
            DO 32 CI = 1, 6
              SBUF(MI, JJ) = SBUF(MI, JJ) + T2(MI, CI) * CBUF(CI, JJ)
   32       CONTINUE
   31     CONTINUE
   30   CONTINUE
        RETURN
      END IF

C     --- L=3 ---
      IF (L .EQ. 3) THEN
        DO 40 JJ = 1, DJ
          DO 41 MI = 1, 7
            SBUF(MI, JJ) = 0.0D0
            DO 42 CI = 1, 10
              SBUF(MI,JJ) = SBUF(MI,JJ) + T3(MI,CI)*CBUF(CI,JJ)
   42       CONTINUE
   41     CONTINUE
   40   CONTINUE
        RETURN
      END IF

C     --- L=4 ---
      IF (L .EQ. 4) THEN
        DO 50 JJ = 1, DJ
          DO 51 MI = 1, 9
            SBUF(MI, JJ) = 0.0D0
            DO 52 CI = 1, 15
              SBUF(MI,JJ) = SBUF(MI,JJ) + T4(MI,CI)*CBUF(CI,JJ)
   52       CONTINUE
   51     CONTINUE
   50   CONTINUE
        RETURN
      END IF

      RETURN
      END

C     ================================================================
C     GINT_1E_STV_SPH:  S, T, V integrals in Spherical GTO basis
C
C     Same interface as GINT_1E_STV.
C     Output buffers: (dsi x dsj), dsN = 2*lN+1.
C
C     All intermediate buffers are 1-dimensional to avoid Fortran
C     leading-dimension mismatches when passing sub-arrays.
C     ================================================================
      SUBROUTINE GINT_1E_STV_SPH(SBUF, TBUF, VBUF, SHLS,
     &                            ATM, NATM, BAS, NBAS, ENV)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  SBUF(*), TBUF(*), VBUF(*)
      INTEGER SHLS(2), ATM(6,*), NATM, BAS(8,*), NBAS
      REAL*8  ENV(*)
      INTEGER ISH, JSH, LI, LJ, DCI, DCJ, DSI, DSJ, II, JJ, IDX
C     1D work buffers (max size: 15*15 = 225 for g-g block)
      REAL*8  SCART(225), TCART(225), VCART(225)
      REAL*8  SI(225), TI(225), VI(225)
      REAL*8  SJ(225), TJ(225), VJ(225)
      REAL*8  ST(225), TT(225), VT(225)

      ISH = SHLS(1); JSH = SHLS(2)
      LI  = BAS(2, ISH); LJ = BAS(2, JSH)
      DCI = (LI+1)*(LI+2)/2; DCJ = (LJ+1)*(LJ+2)/2
      DSI = 2*LI+1;           DSJ = 2*LJ+1

C     Step 1: Cartesian integrals -> SCART(DCI, DCJ) col-major
      CALL GINT_1E_STV(SCART, TCART, VCART, SHLS, ATM, NATM,
     &                  BAS, NBAS, ENV)

C     Step 2: Transform bra i: (DCI, DCJ) -> SI(DSI, DCJ)
      CALL GINT_C2S_1E(LI, SCART, DCI, DCJ, SI, DSI)
      CALL GINT_C2S_1E(LI, TCART, DCI, DCJ, TI, DSI)
      CALL GINT_C2S_1E(LI, VCART, DCI, DCJ, VI, DSI)

C     Step 3: Transpose SI(DSI,DCJ) -> SJ(DCJ,DSI), transform ket j
      DO 10 JJ = 1, DCJ
        DO 11 II = 1, DSI
          SJ((II-1)*DCJ + JJ) = SI((JJ-1)*DSI + II)
          TJ((II-1)*DCJ + JJ) = TI((JJ-1)*DSI + II)
          VJ((II-1)*DCJ + JJ) = VI((JJ-1)*DSI + II)
   11   CONTINUE
   10 CONTINUE
      CALL GINT_C2S_1E(LJ, SJ, DCJ, DSI, ST, DSJ)
      CALL GINT_C2S_1E(LJ, TJ, DCJ, DSI, TT, DSJ)
      CALL GINT_C2S_1E(LJ, VJ, DCJ, DSI, VT, DSJ)

C     Step 4: Transpose ST(DSJ,DSI) -> output SBUF(DSI,DSJ)
      DO 20 JJ = 1, DSJ
        DO 21 II = 1, DSI
          IDX = (JJ-1)*DSI + II
          SBUF(IDX) = ST((II-1)*DSJ + JJ)
          TBUF(IDX) = TT((II-1)*DSJ + JJ)
          VBUF(IDX) = VT((II-1)*DSJ + JJ)
   21   CONTINUE
   20 CONTINUE
      RETURN
      END

C     GINT_2E_SPH:  Two-electron ERI in Spherical GTO basis
C
C     Same interface as GINT_2E.
C     Output buffer: (dsi x dsj x dsk x dsl), dsN = 2*lN+1.
C
C     Algorithm: Cartesian ERI computed, then GINT_C2S_1E applied
C     sequentially to each index.  All work arrays are 1-dimensional
C     to avoid Fortran rank-mismatch errors.
C     ================================================================
      SUBROUTINE GINT_2E_SPH(BUF, SHLS, ATM, NATM, BAS, NBAS, ENV)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  BUF(*)
      INTEGER SHLS(4), ATM(6,*), NATM, BAS(8,*), NBAS
      REAL*8  ENV(*)
      INTEGER ISH,JSH,KSH,LSH,LI,LJ,LK,LL
      INTEGER DCI,DCJ,DCK,DCL,DSI,DSJ,DSK,DSL
      INTEGER II,JJ,KK,MM,IDX,ICOL
      PARAMETER (MXBUF=50625)
      REAL*8  CBUF(MXBUF), W1(MXBUF), W2(MXBUF)

      ISH=SHLS(1); JSH=SHLS(2); KSH=SHLS(3); LSH=SHLS(4)
      LI=BAS(2,ISH); LJ=BAS(2,JSH); LK=BAS(2,KSH); LL=BAS(2,LSH)
      DCI=(LI+1)*(LI+2)/2; DCJ=(LJ+1)*(LJ+2)/2
      DCK=(LK+1)*(LK+2)/2; DCL=(LL+1)*(LL+2)/2
      DSI=2*LI+1; DSJ=2*LJ+1; DSK=2*LK+1; DSL=2*LL+1

C     Step 1: Cartesian ERI -> CBUF(DCI, DCJ*DCK*DCL) column-major
      CALL GINT_2E(CBUF, SHLS, ATM, NATM, BAS, NBAS, ENV)

C     Step 2: Transform index I (fastest): DCI -> DSI
C     Treat CBUF as (DCI, DCJ*DCK*DCL) -> W1 as (DSI, DCJ*DCK*DCL)
      CALL GINT_C2S_1E(LI, CBUF, DCI, DCJ*DCK*DCL, W1, DSI)
C     W1 layout: W1[ II + (JJ-1)*DSI + (KK-1)*DSI*DCJ + (MM-1)*DSI*DCJ*DCK ]

C     Step 3: Transform index J: DCJ -> DSJ
C     Rearrange W1(DSI,DCJ,DCK,DCL) -> W2(DCJ, DSI*DCK*DCL)
      DO 30 MM=1,DCL
        DO 31 KK=1,DCK
          DO 32 JJ=1,DCJ
            DO 33 II=1,DSI
              ICOL = (II-1) + (KK-1)*DSI + (MM-1)*DSI*DCK
              W2(JJ + ICOL*DCJ) =
     &          W1(II + (JJ-1)*DSI + (KK-1)*DSI*DCJ
     &             + (MM-1)*DSI*DCJ*DCK)
   33       CONTINUE
   32     CONTINUE
   31   CONTINUE
   30 CONTINUE
C     Apply: (DCJ, DSI*DCK*DCL) -> (DSJ, DSI*DCK*DCL)
      CALL GINT_C2S_1E(LJ, W2, DCJ, DSI*DCK*DCL, W1, DSJ)
C     Restore W2(DSI,DSJ,DCK,DCL) from W1(DSJ, DSI*DCK*DCL)
      DO 40 MM=1,DCL
        DO 41 KK=1,DCK
          DO 42 JJ=1,DSJ
            DO 43 II=1,DSI
              ICOL = (II-1) + (KK-1)*DSI + (MM-1)*DSI*DCK
              W2(II + (JJ-1)*DSI + (KK-1)*DSI*DSJ
     &           + (MM-1)*DSI*DSJ*DCK) = W1(JJ + ICOL*DSJ)
   43       CONTINUE
   42     CONTINUE
   41   CONTINUE
   40 CONTINUE

C     Step 4: Transform index K: DCK -> DSK
C     Rearrange W2(DSI,DSJ,DCK,DCL) -> W1(DCK, DSI*DSJ*DCL)
      DO 50 MM=1,DCL
        DO 51 KK=1,DCK
          DO 52 JJ=1,DSJ
            DO 53 II=1,DSI
              ICOL = (II-1) + (JJ-1)*DSI + (MM-1)*DSI*DSJ
              W1(KK + ICOL*DCK) =
     &          W2(II + (JJ-1)*DSI + (KK-1)*DSI*DSJ
     &             + (MM-1)*DSI*DSJ*DCK)
   53       CONTINUE
   52     CONTINUE
   51   CONTINUE
   50 CONTINUE
C     Apply: (DCK, DSI*DSJ*DCL) -> (DSK, DSI*DSJ*DCL)
      CALL GINT_C2S_1E(LK, W1, DCK, DSI*DSJ*DCL, W2, DSK)
C     Restore W1(DSI,DSJ,DSK,DCL) from W2(DSK, DSI*DSJ*DCL)
      DO 60 MM=1,DCL
        DO 61 KK=1,DSK
          DO 62 JJ=1,DSJ
            DO 63 II=1,DSI
              ICOL = (II-1) + (JJ-1)*DSI + (MM-1)*DSI*DSJ
              W1(II + (JJ-1)*DSI + (KK-1)*DSI*DSJ
     &           + (MM-1)*DSI*DSJ*DSK) = W2(KK + ICOL*DSK)
   63       CONTINUE
   62     CONTINUE
   61   CONTINUE
   60 CONTINUE

C     Step 5: Transform index L: DCL -> DSL
C     Rearrange W1(DSI,DSJ,DSK,DCL) -> W2(DCL, DSI*DSJ*DSK)
      DO 70 MM=1,DCL
        DO 71 KK=1,DSK
          DO 72 JJ=1,DSJ
            DO 73 II=1,DSI
              ICOL = (II-1) + (JJ-1)*DSI + (KK-1)*DSI*DSJ
              W2(MM + ICOL*DCL) =
     &          W1(II + (JJ-1)*DSI + (KK-1)*DSI*DSJ
     &             + (MM-1)*DSI*DSJ*DSK)
   73       CONTINUE
   72     CONTINUE
   71   CONTINUE
   70 CONTINUE
C     Apply: (DCL, DSI*DSJ*DSK) -> (DSL, DSI*DSJ*DSK)
      CALL GINT_C2S_1E(LL, W2, DCL, DSI*DSJ*DSK, W1, DSL)
C     Write BUF(DSI,DSJ,DSK,DSL) from W1(DSL, DSI*DSJ*DSK)
      DO 80 MM=1,DSL
        DO 81 KK=1,DSK
          DO 82 JJ=1,DSJ
            DO 83 II=1,DSI
              ICOL = (II-1) + (JJ-1)*DSI + (KK-1)*DSI*DSJ
              IDX  = II + (JJ-1)*DSI + (KK-1)*DSI*DSJ
     &             + (MM-1)*DSI*DSJ*DSK
              BUF(IDX) = W1(MM + ICOL*DSL)
   83       CONTINUE
   82     CONTINUE
   81   CONTINUE
   80 CONTINUE
      RETURN
      END
