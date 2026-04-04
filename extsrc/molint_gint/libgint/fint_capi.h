/*
 * ==================================================================
 * gint_mol.h
 * ==================================================================
 * Public C API for libgint: converts emol-style molecular data into
 * GINT internal format and computes 1e/2e integrals.
 *
 * This header is copied to the parent directory (molint/) by
 * libgint/Makefile so that molint_gint.c can include it as:
 *   #include "gint_mol.h"
 *
 * Input convention (emol):
 *   prim_coe[i] = raw_coef * gto_norm(l, alpha) * sqrt((2l+1)/(4*pi))
 *
 * Output buffer layout (all flat 1D, column-major / Fortran order):
 *   1E:  buf[jp * nci + ip]                          (ip fastest)
 *   2E:  buf[((lp*nck+kp)*ncj+jp)*nci+ip]            (ip fastest)
 *   -> caller transposes 2E to emol order if needed:
 *        em[((ip*ncj+jp)*nck+kp)*ncl+lp]
 *
 * Shell indices passed to calc_* are 0-based (C convention).
 * Internally converted to 1-based for Fortran.
 *
 * angular_type:
 *   GINT_CART (0) -> Cartesian GTOs, (l+1)(l+2)/2 components
 *   GINT_SPH  (1) -> Spherical GTOs, 2l+1 components
 * ==================================================================
 */
#ifndef GINT_MOL_H
#define GINT_MOL_H

#ifdef __cplusplus
extern "C" {
#endif

/* Angular type flags */
#define GINT_CART 0
#define GINT_SPH  1

/* Maximum Cartesian component buffer size: 15*15 for l=4 (g-shell) */
#define GINT_MAX_CART 225

/* ================================================================
 * gint_mol_init
 *
 * Initialize libgint from emol-style molecular data.
 * Must be called before any calc_* function.
 * Can be called multiple times (re-initializes state).
 *
 * Arguments:
 *   numb_atom   : number of atoms
 *   numb_shell  : number of shells
 *   numb_prim   : total number of primitives (sum over shells)
 *   shel_lqn    : [numb_shell] angular momentum per shell
 *   shel_tem    : [numb_shell] number of primitives per shell
 *   shel_atm    : [numb_shell] atom index per shell (0-based)
 *   shel_add    : [numb_shell] offset into prim_exp/prim_coe
 *   atom_charg  : [numb_atom]  nuclear charge Z
 *   atom_xyz    : [numb_atom*3] coordinates in Bohr (x,y,z interleaved)
 *   prim_exp    : [numb_prim]  primitive exponents
 *   prim_coe    : [numb_prim]  emol-normalized contraction coefficients
 *   angular_type: GINT_CART or GINT_SPH
 * ================================================================ */
void gint_mol_init(
    int     numb_atom,
    int     numb_shell,
    int     numb_prim,
    int    *shel_lqn,
    int    *shel_tem,
    int    *shel_atm,
    int    *shel_add,
    double *atom_charg,
    double *atom_xyz,
    double *prim_exp,
    double *prim_coe,
    int     angular_type);

/* ================================================================
 * gint_mol_finalize
 *
 * Release all internally allocated memory.
 * Optional but recommended when done.
 * ================================================================ */
void gint_mol_finalize(void);

/* ================================================================
 * gint_mol_shell_comp
 *
 * Return number of components for shell ish in the current mode
 * (Cartesian: (l+1)(l+2)/2,  Spherical: 2l+1).
 * ================================================================ */
int gint_mol_shell_comp(int ish);

/* ================================================================
 * gint_mol_shell_lqn
 *
 * Return angular momentum of shell ish.
 * ================================================================ */
int gint_mol_shell_lqn(int ish);

/* ================================================================
 * 1E integral routines
 *
 * All three integrals (S, T, V) are computed in a single call.
 * Pass NULL for any buffer you don't need.
 *
 * Output size: nci * ncj  where nci = gint_mol_shell_comp(ish)
 *                                ncj = gint_mol_shell_comp(jsh)
 * Layout: buf[jp * nci + ip]  (ip = bra AO, fastest index)
 *
 * V = sum_A (-Z_A) <i | 1/|r-R_A| | j>
 * ================================================================ */
void gint_mol_1e_stv(
    double *sbuf,   /* overlap       (may be NULL) */
    double *tbuf,   /* kinetic       (may be NULL) */
    double *vbuf,   /* nuc.attract.  (may be NULL) */
    int     ish,
    int     jsh);

/* Return number of atoms (set by gint_mol_init) */
int gint_mol_natm(void);

/* Raw nuclear attraction: buf[k*nci*ncj + jp*nci + ip]
 * k = nucleus index (0-based), size = nci*ncj*natm            */
void gint_mol_1e_vnuc_raw(double *buf, int ish, int jsh);

/* ================================================================
 * 2E integral routine
 *
 * Computes ERI (ish jsh | ksh lsh) in the current angular mode.
 *
 * Output size: nca * ncb * ncc * ncd
 * Layout (GINT/Fortran col-major, ip fastest):
 *   buf[((lp*ncc+kp)*ncb+jp)*nca+ip]
 *
 * To convert to emol order ((ip*ncb+jp)*ncc+kp)*ncd+lp,
 * use gint_mol_2e_transpose() or do it in the caller.
 * ================================================================ */
void gint_mol_2e(
    double *buf,
    int     ish,
    int     jsh,
    int     ksh,
    int     lsh);

/* ================================================================
 * gint_mol_2e_transpose
 *
 * In-place transpose of a 2E buffer from GINT col-major order to
 * emol row-major order.
 * buf must have nca*ncb*ncc*ncd elements.
 * ================================================================ */
void gint_mol_2e_transpose(
    double *buf,
    int nca, int ncb, int ncc, int ncd);

#ifdef __cplusplus
}
#endif

#endif /* GINT_MOL_H */
