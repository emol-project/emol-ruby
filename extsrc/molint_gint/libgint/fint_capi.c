/*
 * ==================================================================
 * gint_mol.c
 * ==================================================================
 * Implementation of the gint_mol C API.
 * Converts emol-style molecular data to GINT internal format,
 * then calls Fortran GINT subroutines for integral evaluation.
 *
 * Coefficient conversion (emol -> GINT raw):
 *   emol_coe = raw * gto_norm(l, alpha) * sqrt((2l+1)/(4*pi))
 *   raw      = emol_coe / gto_norm(l, alpha) / sqrt((2l+1)/(4*pi))
 *   GINT_INIT then internally applies FNFB normalization.
 * ==================================================================
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fint_capi.h"

/* ================================================================
 * Fortran GINT entry points
 * (gfortran default: lowercase name + trailing underscore,
 *  all arguments passed by pointer)
 * ================================================================ */
extern void fmt_initialize_(const int *max_lqn);
extern void gint_init_(int *atm, int *natm, int *bas, int *nbas, double *env);
extern void gint_1e_stv_(double *s, double *t, double *v,
                          int *shls, int *atm, int *natm,
                          int *bas,  int *nbas, double *env);
extern void gint_1e_stv_sph_(double *s, double *t, double *v,
                               int *shls, int *atm, int *natm,
                               int *bas,  int *nbas, double *env);
extern void gint_2e_(double *buf, int *shls,
                      int *atm, int *natm, int *bas, int *nbas, double *env);
extern void gint_2e_sph_(double *buf, int *shls,
                          int *atm, int *natm, int *bas, int *nbas, double *env);
extern void gint_1e_vnuc_raw_(double *buf, int *shls,
                               int *atm, int *natm,
                               int *bas, int *nbas, double *env);


/* ================================================================
 * Internal state
 * ================================================================ */
static int     g_angular_type = GINT_CART;
static int     g_natm         = 0;
static int     g_nbas         = 0;
static int    *g_atm          = NULL;   /* [natm * 6] */
static int    *g_bas          = NULL;   /* [nbas * 8], atom index is 1-based */
static double *g_env          = NULL;
static int    *g_shel_lqn     = NULL;   /* [nbas] angular momentum */
static int     g_fmt_ready    = 0;

/* ================================================================
 * Helper: radial normalization factor  (PySCF / libcint convention)
 *   N(l, alpha) = 1 / sqrt( Gamma((2l+3)/2) / (2*(2*alpha)^((2l+3)/2)) )
 * ================================================================ */
static double gto_norm(int l, double alpha) {
    int    n  = 2 * l + 2;
    double n1 = (n + 1) * 0.5;
    double gi = tgamma(n1) / (2.0 * pow(2.0 * alpha, n1));
    return 1.0 / sqrt(gi);
}

/* ================================================================
 * Helper: GINT spherical scale factor
 *
 * GINT's Cartesian-to-Spherical transformation matrix (GINT_C2S_1E)
 * includes a factor of SCALE(l) = sqrt(4*pi/(2l+1)).
 * This means GINT_SPH output = libcint output * SCALE(li) * SCALE(lj).
 * We divide it out here so that the output matches libcint convention.
 *
 * l=0,1: SCALE=1 (no correction needed)
 * l>=2:  SCALE = sqrt(4*pi/(2l+1))
 * ================================================================ */
static double gint_scale(int l) {
    if (l <= 1) return 1.0;
    return sqrt(4.0 * M_PI / (2 * l + 1));
}

/* ================================================================
 * gint_mol_init
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
    int     angular_type)
{
    gint_mol_finalize();

    g_angular_type = angular_type;
    g_natm         = numb_atom;
    g_nbas         = numb_shell;

    /* ---- allocate env ----
     * coords:       3 * natm
     * per shell:    nprim exponents + nprim coefficients
     * slack:        64
     */
    int env_size = 3 * numb_atom + 2 * numb_prim + 64;
    g_env = (double *)calloc(env_size, sizeof(double));

    /* ---- atm array: 6 ints per atom ----
     * [0] = nuclear charge Z
     * [1] = 0-based env offset for coordinates
     * [2..5] = 0
     */
    g_atm = (int *)calloc(numb_atom * 6, sizeof(int));
    int env_ptr = 0;

    for (int i = 0; i < numb_atom; i++) {
        g_atm[i * 6 + 0] = (int)atom_charg[i];
        g_atm[i * 6 + 1] = env_ptr;
        g_env[env_ptr + 0] = atom_xyz[i * 3 + 0];
        g_env[env_ptr + 1] = atom_xyz[i * 3 + 1];
        g_env[env_ptr + 2] = atom_xyz[i * 3 + 2];
        env_ptr += 3;
    }

    /* ---- bas array: 8 ints per shell ----
     * [0] = atom index (*** 1-based *** for GINT Fortran)
     * [1] = angular momentum l
     * [2] = nprim
     * [3] = nctr (must be 1)
     * [4] = kappa (0)
     * [5] = 0-based env offset for exponents
     * [6] = 0-based env offset for coefficients
     * [7] = 0
     */
    g_bas      = (int *)calloc(numb_shell * 8, sizeof(int));
    g_shel_lqn = (int *)malloc(numb_shell * sizeof(int));

    for (int ish = 0; ish < numb_shell; ish++) {
        int    l     = shel_lqn[ish];
        int    nprim = shel_tem[ish];
        int    iatom = shel_atm[ish];   /* 0-based from emol */
        int    poff  = shel_add[ish];

        g_shel_lqn[ish] = l;

        g_bas[ish * 8 + 0] = iatom + 1; /* 1-based for GINT */
        g_bas[ish * 8 + 1] = l;
        g_bas[ish * 8 + 2] = nprim;
        g_bas[ish * 8 + 3] = 1;
        g_bas[ish * 8 + 4] = 0;
        g_bas[ish * 8 + 5] = env_ptr;   /* exponent offset */

        for (int ip = 0; ip < nprim; ip++)
            g_env[env_ptr + ip] = prim_exp[poff + ip];
        env_ptr += nprim;

        g_bas[ish * 8 + 6] = env_ptr;   /* coefficient offset */

        /* Convert emol_coe -> raw for GINT_INIT.
         * emol_coe = raw * gto_norm(l,alpha) * sqrt((2l+1)/(4*pi))
         * GINT_INIT will multiply by FNFB internally.
         */
        double ang_factor = sqrt((2.0 * l + 1.0) / (4.0 * M_PI));
        for (int ip = 0; ip < nprim; ip++) {
            double alpha = prim_exp[poff + ip];
            g_env[env_ptr + ip] =
                prim_coe[poff + ip] / gto_norm(l, alpha) / ang_factor;
        }
        env_ptr += nprim;
    }

    /* Initialize Boys function table (once is enough) */
    if (!g_fmt_ready) {
        int max_lqn = 4;
        fmt_initialize_(&max_lqn);
        g_fmt_ready = 1;
    }

    /* Initialize GINT internal state */
    gint_init_(g_atm, &g_natm, g_bas, &g_nbas, g_env);
}

/* ================================================================
 * gint_mol_finalize
 * ================================================================ */
void gint_mol_finalize(void) {
    if (g_atm)      { free(g_atm);      g_atm      = NULL; }
    if (g_bas)      { free(g_bas);      g_bas      = NULL; }
    if (g_env)      { free(g_env);      g_env      = NULL; }
    if (g_shel_lqn) { free(g_shel_lqn);g_shel_lqn = NULL; }
    g_natm = 0;
    g_nbas = 0;
}

/* ================================================================
 * gint_mol_shell_comp
 * ================================================================ */
int gint_mol_shell_comp(int ish) {
    int l = g_shel_lqn[ish];
    return (g_angular_type == GINT_SPH) ? (2 * l + 1) : (l + 1) * (l + 2) / 2;
}

/* ================================================================
 * gint_mol_shell_lqn
 * ================================================================ */
int gint_mol_shell_lqn(int ish) {
    return g_shel_lqn[ish];
}

int gint_mol_natm(void) {
    return g_natm;
}

/* ================================================================
 * gint_mol_1e_stv
 * ================================================================ */
void gint_mol_1e_stv(double *sbuf, double *tbuf, double *vbuf,
                      int ish, int jsh)
{
    int f_shls[2] = { ish + 1, jsh + 1 };   /* 1-based for Fortran */

    /* Scratch buffers for unwanted outputs */
    double scratch1[GINT_MAX_CART];
    double scratch2[GINT_MAX_CART];

    double *s = sbuf ? sbuf : scratch1;
    double *t = tbuf ? tbuf : scratch2;
    double *v = vbuf ? vbuf : scratch1;

    /* If two of the three pointers are NULL and would alias scratch,
     * we need a second scratch.  Simplest: always use two scratches. */
    double tmp_s[GINT_MAX_CART];
    double tmp_t[GINT_MAX_CART];
    double tmp_v[GINT_MAX_CART];

    double *ps = sbuf ? sbuf : tmp_s;
    double *pt = tbuf ? tbuf : tmp_t;
    double *pv = vbuf ? vbuf : tmp_v;

    if (g_angular_type == GINT_SPH) {
        gint_1e_stv_sph_(ps, pt, pv,
                          f_shls, g_atm, &g_natm, g_bas, &g_nbas, g_env);

        /* Remove GINT's built-in SCALE factor from the C2S transform.
         * GINT_SPH output = libcint output * SCALE(li) * SCALE(lj)
         * Divide out so the result matches libcint/PySCF convention.  */
        int li   = g_shel_lqn[ish];
        int lj   = g_shel_lqn[jsh];
        int nci  = 2 * li + 1;
        int ncj  = 2 * lj + 1;
        double inv = 1.0 / (gint_scale(li) * gint_scale(lj));
        if (inv != 1.0) {
            for (int k = 0; k < nci * ncj; k++) {
                ps[k] *= inv;
                pt[k] *= inv;
                pv[k] *= inv;
            }
        }
    } else {
        gint_1e_stv_(ps, pt, pv,
                      f_shls, g_atm, &g_natm, g_bas, &g_nbas, g_env);
    }
    (void)s; (void)t; (void)v;
}
                               
void gint_mol_1e_vnuc_raw(double *buf, int ish, int jsh)
{
    int f_shls[2] = { ish + 1, jsh + 1 };
    gint_1e_vnuc_raw_(buf, f_shls, g_atm, &g_natm, g_bas, &g_nbas, g_env);
}

/* ================================================================
 * gint_mol_2e
 * ================================================================ */
void gint_mol_2e(double *buf, int ish, int jsh, int ksh, int lsh)
{
    int f_shls[4] = { ish + 1, jsh + 1, ksh + 1, lsh + 1 };

    if (g_angular_type == GINT_SPH) {
        gint_2e_sph_(buf, f_shls, g_atm, &g_natm, g_bas, &g_nbas, g_env);

        /* Remove GINT's built-in SCALE factors from the C2S transforms.
         * GINT_SPH output = libcint output * SCALE(la)*SCALE(lb)*SCALE(lc)*SCALE(ld) */
        int la = g_shel_lqn[ish], lb = g_shel_lqn[jsh];
        int lc = g_shel_lqn[ksh], ld = g_shel_lqn[lsh];
        int nca = 2*la+1, ncb = 2*lb+1, ncc = 2*lc+1, ncd = 2*ld+1;
        double inv = 1.0 / (gint_scale(la) * gint_scale(lb)
                          * gint_scale(lc) * gint_scale(ld));
        if (inv != 1.0) {
            int ntot = nca * ncb * ncc * ncd;
            for (int k = 0; k < ntot; k++)
                buf[k] *= inv;
        }
    } else {
        gint_2e_(buf, f_shls, g_atm, &g_natm, g_bas, &g_nbas, g_env);
    }
}

/* ================================================================
 * gint_mol_2e_transpose
 *
 * Convert GINT output layout:
 *   src[((lp*ncc+kp)*ncb+jp)*nca+ip]   (ip fastest)
 * to emol layout:
 *   dst[((ip*ncb+jp)*ncc+kp)*ncd+lp]   (lp fastest... no: ip fastest)
 *
 * Done in-place via a temporary copy.
 * ================================================================ */
void gint_mol_2e_transpose(double *buf,
                             int nca, int ncb, int ncc, int ncd)
{
    int ntot = nca * ncb * ncc * ncd;
    double *tmp = (double *)malloc(ntot * sizeof(double));
    memcpy(tmp, buf, ntot * sizeof(double));

    for (int ip = 0; ip < nca; ip++)
    for (int jp = 0; jp < ncb; jp++)
    for (int kp = 0; kp < ncc; kp++)
    for (int lp = 0; lp < ncd; lp++) {
        int gint_idx = ((lp * ncc + kp) * ncb + jp) * nca + ip;
        int em_idx   = ((ip * ncb + jp) * ncc + kp) * ncd + lp;
        buf[em_idx]  = tmp[gint_idx];
    }
    free(tmp);
}
