/*
 * molint_libcint.c
 *
 * Ruby C extension that replaces the original OpenFMO-based Molint FFI
 * with libcint as the integral engine.
 *
 * Provides the same Ruby classes and interface:
 *   Sint0  : init(14 args), calc(ish, jsh)    -> [inttype, nsize, array]
 *   Tint0  : init(14 args), calc(ish, jsh)    -> [inttype, nsize, array]
 *   Vint0  : init(14 args), calc(ish, jsh)    -> [inttype, nsize, array]
 *   Eri0   : init(14 args), calc(ish,jsh,ksh,lsh) -> [inttype, nsize, array]
 *            exchange_order()                  -> [a, b, c, d]
 *
 * Integral convention is selectable at init() time via the 14th argument:
 *   angular_type = 0  ->  Cartesian (6d, 10f), GAMESS ordering and normalization
 *   angular_type = 1  ->  Spherical (5d, 7f),  libcint real-spherical-harmonic ordering
 *
 * Author: Generated for emol project
 * Integral engine: libcint (Qiming Sun)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ruby.h"
#include "cint.h"

/* Try to include cint_funcs.h (present in v6+, absent in v3) */
#ifdef HAVE_CINT_FUNCS_H
#include "cint_funcs.h"
#endif

/* ================================================================
 * Angular type flag (0 = Cartesian, 1 = Spherical)
 * Set once per init() call.
 * ================================================================ */
static int g_use_spherical = 0;

/* ================================================================
 * Global data: libcint's atm, bas, env arrays
 * ================================================================ */
static FINT *g_atm = NULL;
static FINT  g_natm = 0;
static FINT *g_bas = NULL;
static FINT  g_nbas = 0;
static double *g_env = NULL;
static int    g_env_size = 0;

/* Shell info for emol compatibility */
static int *g_shel_lqn = NULL;   /* angular momentum per shell */
static int *g_shel_ini = NULL;   /* starting basis function index per shell */
static int  g_nbf = 0;          /* total number of basis functions */

/* ================================================================
 * Component count helpers
 * ================================================================ */

/* Cartesian component count: (l+1)*(l+2)/2 */
static int cart_comp(int l) {
    return (l + 1) * (l + 2) / 2;
}

/* Spherical component count: 2*l + 1 */
static int sph_comp(int l) {
    return 2 * l + 1;
}

/* Dispatch: return component count based on current angular type */
static int shell_comp(int l) {
    return g_use_spherical ? sph_comp(l) : cart_comp(l);
}

/* ================================================================
 * GAMESS normalization factor for Cartesian component x^a y^b z^c
 *   factor = sqrt(l! / (a! b! c!))  where l = a+b+c
 *
 * Only used in Cartesian mode.
 * ================================================================ */
static double gamess_cart_norm(int a, int b, int c) {
    static const double fact[] = {1,1,2,6,24,120,720,5040};
    int l = a + b + c;
    return sqrt(fact[l] / (fact[a] * fact[b] * fact[c]));
}

/* ================================================================
 * Cartesian component ordering tables (Cartesian mode only)
 *
 * libcint order: lexicographic (x highest first)
 *   d: xx,xy,xz,yy,yz,zz
 *   f: xxx,xxy,xxz,xyy,xyz,xzz,yyy,yyz,yzz,zzz
 *
 * GAMESS/emol order:
 *   d: xx,yy,zz,xy,xz,yz
 *   f: xxx,yyy,zzz,xxy,xxz,xyy,yyz,xzz,yzz,xyz
 * ================================================================ */

typedef struct { int a, b, c; } CartPow;

static void libcint_cart_powers(int l, CartPow *out) {
    int idx = 0;
    for (int ax = l; ax >= 0; ax--) {
        for (int ay = l - ax; ay >= 0; ay--) {
            int az = l - ax - ay;
            out[idx].a = ax;
            out[idx].b = ay;
            out[idx].c = az;
            idx++;
        }
    }
}

static const int gamess_order_d[][3] = {
    {2,0,0},{0,2,0},{0,0,2},{1,1,0},{1,0,1},{0,1,1}
};
static const int gamess_order_f[][3] = {
    {3,0,0},{0,3,0},{0,0,3},{2,1,0},{2,0,1},{1,2,0},
    {0,2,1},{1,0,2},{0,1,2},{1,1,1}
};
static const int gamess_order_g[][3] = {
    {4,0,0},{0,4,0},{0,0,4},{3,1,0},{3,0,1},{1,3,0},
    {0,3,1},{1,0,3},{0,1,3},{2,2,0},{2,0,2},{0,2,2},
    {2,1,1},{1,2,1},{1,1,2}
};

static void build_cart_mapping(int l, int *perm, double *norm_fac) {
    int ncomp = cart_comp(l);
    CartPow *lc = (CartPow *)malloc(ncomp * sizeof(CartPow));
    libcint_cart_powers(l, lc);

    const int (*gamess_ord)[3];
    if (l == 0 || l == 1) {
        for (int i = 0; i < ncomp; i++) {
            perm[i] = i;
            norm_fac[i] = 1.0;
        }
        free(lc);
        return;
    } else if (l == 2) {
        gamess_ord = gamess_order_d;
    } else if (l == 3) {
        gamess_ord = gamess_order_f;
    } else {
        gamess_ord = gamess_order_g;
    }

    for (int i = 0; i < ncomp; i++) {
        int ga = gamess_ord[i][0];
        int gb = gamess_ord[i][1];
        int gc = gamess_ord[i][2];
        norm_fac[i] = gamess_cart_norm(ga, gb, gc);
        for (int j = 0; j < ncomp; j++) {
            if (lc[j].a == ga && lc[j].b == gb && lc[j].c == gc) {
                perm[i] = j;
                break;
            }
        }
    }
    free(lc);
}

/* Pre-computed Cartesian mapping tables (up to l=4) */
#define MAX_L 4
static int    g_perm[MAX_L+1][15+1];
static double g_norm[MAX_L+1][15+1];
static int g_cart_mapping_ready = 0;

static void ensure_cart_mapping_tables(void) {
    if (g_cart_mapping_ready) return;
    for (int l = 0; l <= MAX_L; l++) {
        build_cart_mapping(l, g_perm[l], g_norm[l]);
    }
    g_cart_mapping_ready = 1;
}

/* ================================================================
 * gto_norm: PySCF/libcint normalization for radial part
 * ================================================================ */
static double gto_norm(int l, double alpha) {
    int n = 2 * l + 2;
    double n1 = (n + 1) * 0.5;
    double gi = tgamma(n1) / (2.0 * pow(2.0 * alpha, n1));
    return 1.0 / sqrt(gi);
}

/* ================================================================
 * Build libcint data structures from emol's shell data
 * ================================================================ */
static void build_libcint_data(
    int numb_atom, int numb_shell, int numb_prim,
    int *shel_lqn, int *shel_tem, int *shel_atm, int *shel_add,
    double *atom_charg, double *atom_xyz,
    double *prim_exp, double *prim_coe)
{
    /* Free previous data */
    if (g_atm) { free(g_atm); g_atm = NULL; }
    if (g_bas) { free(g_bas); g_bas = NULL; }
    if (g_env) { free(g_env); g_env = NULL; }
    if (g_shel_lqn) { free(g_shel_lqn); g_shel_lqn = NULL; }
    if (g_shel_ini) { free(g_shel_ini); g_shel_ini = NULL; }

    g_natm = numb_atom;
    g_nbas = numb_shell;

    /* --- Build env array --- */
    int env_ptr = PTR_ENV_START;

    int total_env = env_ptr;
    total_env += numb_atom * 4;
    total_env += numb_prim * 2;
    total_env += 100;

    g_env = (double *)calloc(total_env, sizeof(double));
    g_env_size = total_env;

    /* --- Build atm array --- */
    g_atm = (FINT *)calloc(numb_atom * ATM_SLOTS, sizeof(FINT));
    for (int i = 0; i < numb_atom; i++) {
        g_atm[i * ATM_SLOTS + CHARGE_OF] = (FINT)atom_charg[i];
        g_atm[i * ATM_SLOTS + PTR_COORD] = env_ptr;
        g_env[env_ptr + 0] = atom_xyz[i * 3 + 0];
        g_env[env_ptr + 1] = atom_xyz[i * 3 + 1];
        g_env[env_ptr + 2] = atom_xyz[i * 3 + 2];
        env_ptr += 3;
        g_atm[i * ATM_SLOTS + NUC_MOD_OF] = 1;
        g_atm[i * ATM_SLOTS + PTR_ZETA] = env_ptr;
        g_env[env_ptr] = 0.0;
        env_ptr += 1;
    }

    /* --- Build bas array --- */
    g_bas = (FINT *)calloc(numb_shell * BAS_SLOTS, sizeof(FINT));

    g_shel_lqn = (int *)malloc(numb_shell * sizeof(int));
    g_shel_ini = (int *)malloc(numb_shell * sizeof(int));
    int bf_count = 0;

    for (int ish = 0; ish < numb_shell; ish++) {
        int l = shel_lqn[ish];
        int nprim = shel_tem[ish];
        int iatom = shel_atm[ish];
        int prim_offset = shel_add[ish];

        g_shel_lqn[ish] = l;
        g_shel_ini[ish] = bf_count;
        bf_count += shell_comp(l);   /* cart or sph depending on g_use_spherical */

        g_bas[ish * BAS_SLOTS + ATOM_OF]   = iatom;
        g_bas[ish * BAS_SLOTS + ANG_OF]    = l;
        g_bas[ish * BAS_SLOTS + NPRIM_OF]  = nprim;
        g_bas[ish * BAS_SLOTS + NCTR_OF]   = 1;
        g_bas[ish * BAS_SLOTS + KAPPA_OF]  = 0;
        g_bas[ish * BAS_SLOTS + PTR_EXP]   = env_ptr;

        /* Store exponents */
        for (int ip = 0; ip < nprim; ip++) {
            g_env[env_ptr + ip] = prim_exp[prim_offset + ip];
        }
        env_ptr += nprim;

        g_bas[ish * BAS_SLOTS + PTR_COEFF] = env_ptr;

        /* Store coefficients: convert from emol to libcint convention.
         *
         * emol_coe = original_coef * gto_norm(l,alpha) * sqrt((2l+1)/(4pi))
         * libcint wants: original_coef * gto_norm(l,alpha)
         * So: coef_libcint = emol_coe / sqrt((2l+1)/(4pi))
         */
        double ang_factor = sqrt((2.0 * l + 1.0) / (4.0 * M_PI));
        for (int ip = 0; ip < nprim; ip++) {
            g_env[env_ptr + ip] = prim_coe[prim_offset + ip] / ang_factor;
        }

        /* Contraction renormalization */
        {
            double *es = &prim_exp[prim_offset];
            double *cs = &g_env[env_ptr];
            double ovlp = 0.0;
            for (int i = 0; i < nprim; i++) {
                for (int j = 0; j < nprim; j++) {
                    int n = 2 * l + 2;
                    double n1 = (n + 1) * 0.5;
                    double gi = tgamma(n1) / (2.0 * pow(es[i] + es[j], n1));
                    ovlp += cs[i] * cs[j] * gi;
                }
            }
            if (ovlp > 0.0) {
                double renorm = 1.0 / sqrt(ovlp);
                for (int ip = 0; ip < nprim; ip++) {
                    cs[ip] *= renorm;
                }
            }
        }

        env_ptr += nprim;
    }

    g_nbf = bf_count;

    /* Build Cartesian mapping tables if needed */
    if (!g_use_spherical) {
        ensure_cart_mapping_tables();
    }
}

/* ================================================================
 * 1-electron integral: compute shell pair (ish, jsh)
 *
 * Cartesian mode: reorder libcint output to GAMESS ordering.
 * Spherical mode: use libcint _sph output directly.
 *
 * Output layout: dint[jp * ncomp_i + ip]  (ip fast, jp slow)
 * ================================================================ */
static int calc_1e(
    CINTIntegralFunction *int1e_func_cart,
    CINTIntegralFunction *int1e_func_sph,
    int ish, int jsh,
    int *out_inttype, int *out_nsize, double *out_dint)
{
    int li = g_shel_lqn[ish];
    int lj = g_shel_lqn[jsh];
    int nci = shell_comp(li);
    int ncj = shell_comp(lj);

    *out_inttype = (li + 1) * 10 + (lj + 1);
    *out_nsize = nci * ncj;

    double *buf = (double *)calloc(nci * ncj, sizeof(double));
    FINT shls[2] = { ish, jsh };

    CACHE_SIZE_T not0;
    if (g_use_spherical) {
        not0 = int1e_func_sph(buf, NULL, shls,
                              g_atm, g_natm, g_bas, g_nbas, g_env, NULL, NULL);
    } else {
        not0 = int1e_func_cart(buf, NULL, shls,
                               g_atm, g_natm, g_bas, g_nbas, g_env, NULL, NULL);
    }

    if (!not0) {
        memset(out_dint, 0, nci * ncj * sizeof(double));
        free(buf);
        return 0;
    }

    if (g_use_spherical) {
        /* Spherical: no reordering needed */
        memcpy(out_dint, buf, nci * ncj * sizeof(double));
    } else {
        /* Cartesian: reorder from libcint order to GAMESS order */
        int *perm_i = g_perm[li];
        int *perm_j = g_perm[lj];
        for (int jp = 0; jp < ncj; jp++) {
            for (int ip = 0; ip < nci; ip++) {
                int lc_idx = perm_j[jp] * nci + perm_i[ip];
                out_dint[jp * nci + ip] = buf[lc_idx];
            }
        }
    }

    free(buf);
    return 0;
}

/* ================================================================
 * 2-electron integral: compute shell quartet (ish, jsh, ksh, lsh)
 *
 * Shell reordering for exchange_order() compatibility is preserved.
 *
 * Cartesian mode: reorder components from libcint to GAMESS ordering.
 * Spherical mode: transpose index convention only (no component reorder).
 *
 * Output: flat array in order ((ip*ncj+jp)*nck+kp)*ncl+lp
 * ================================================================ */

static int g_last_a = 0, g_last_b = 0, g_last_c = 0, g_last_d = 0;

static int calc_2e(int ish, int jsh, int ksh, int lsh,
                   int *out_inttype, int *out_nsize, double *out_eris)
{
    int a = ish, b = jsh, c = ksh, d = lsh;

    /* Ensure La >= Lb */
    if (g_shel_lqn[a] < g_shel_lqn[b]) { int t = a; a = b; b = t; }
    /* Ensure Lc >= Ld */
    if (g_shel_lqn[c] < g_shel_lqn[d]) { int t = c; c = d; d = t; }
    /* Ensure La+Lb >= Lc+Ld */
    if (g_shel_lqn[a] + g_shel_lqn[b] < g_shel_lqn[c] + g_shel_lqn[d]) {
        int t;
        t = a; a = c; c = t;
        t = b; b = d; d = t;
    }

    g_last_a = a; g_last_b = b; g_last_c = c; g_last_d = d;

    int la = g_shel_lqn[a], lb = g_shel_lqn[b];
    int lc = g_shel_lqn[c], ld = g_shel_lqn[d];
    int nca = shell_comp(la), ncb = shell_comp(lb);
    int ncc = shell_comp(lc), ncd = shell_comp(ld);
    int ntot = nca * ncb * ncc * ncd;

    *out_inttype = ((la+1)*10 + (lb+1))*100 + (lc+1)*10 + (ld+1);
    *out_nsize = ntot;

    double *buf = (double *)calloc(ntot, sizeof(double));
    FINT shls[4] = { a, b, c, d };

    CACHE_SIZE_T not0;
    if (g_use_spherical) {
        not0 = int2e_sph(buf, NULL, shls,
                         g_atm, g_natm, g_bas, g_nbas, g_env, NULL, NULL);
    } else {
        not0 = int2e_cart(buf, NULL, shls,
                          g_atm, g_natm, g_bas, g_nbas, g_env, NULL, NULL);
    }

    if (!not0) {
        memset(out_eris, 0, ntot * sizeof(double));
        free(buf);
        return 0;
    }

    /*
     * libcint 2e output: column-major [i,j,k,l].
     * buf[ ((l_idx * ncc + k_idx) * ncb + j_idx) * nca + i_idx ]
     *   where i is fastest.
     *
     * emol output: ((ip*ncb+jp)*ncc+kp)*ncd+lp
     *
     * In Cartesian mode, additionally apply the GAMESS permutation.
     * In Spherical mode, only the index transpose is needed.
     */

    if (g_use_spherical) {
        /* Spherical: no component reordering, just index transpose */
        for (int ip = 0; ip < nca; ip++) {
            for (int jp = 0; jp < ncb; jp++) {
                for (int kp = 0; kp < ncc; kp++) {
                    for (int lp = 0; lp < ncd; lp++) {
                        int lc_idx = ((lp * ncc + kp) * ncb + jp) * nca + ip;
                        int em_idx = ((ip * ncb + jp) * ncc + kp) * ncd + lp;
                        out_eris[em_idx] = buf[lc_idx];
                    }
                }
            }
        }
    } else {
        /* Cartesian: apply GAMESS permutation + index transpose */
        int *perm_a = g_perm[la], *perm_b = g_perm[lb];
        int *perm_c = g_perm[lc], *perm_d = g_perm[ld];

        for (int ip = 0; ip < nca; ip++) {
            for (int jp = 0; jp < ncb; jp++) {
                for (int kp = 0; kp < ncc; kp++) {
                    for (int lp = 0; lp < ncd; lp++) {
                        int lc_i = perm_a[ip], lc_j = perm_b[jp];
                        int lc_k = perm_c[kp], lc_l = perm_d[lp];
                        int lc_idx = ((lc_l * ncc + lc_k) * ncb + lc_j) * nca + lc_i;
                        int em_idx = ((ip * ncb + jp) * ncc + kp) * ncd + lp;
                        out_eris[em_idx] = buf[lc_idx];
                    }
                }
            }
        }
    }

    free(buf);
    return 0;
}

/* ================================================================
 * Ruby wrapper functions
 * ================================================================ */

/* --- Generic init wrapper (same for all 4 classes) ---
 * 14 arguments: 13 original + angular_type (0=cart, 1=sph)
 */
static VALUE generic_wrapper_init(VALUE self,
    VALUE _numb_atom, VALUE _numb_shell, VALUE _numb_prim,
    VALUE _shel_lqn, VALUE _shel_tem, VALUE _shel_atm, VALUE _shel_add,
    VALUE _atom_charg, VALUE _atom_xyz, VALUE _prim_exp, VALUE _prim_coe,
    VALUE _tol, VALUE _level_dbg_print, VALUE _angular_type)
{
    int i, numb_atom, numb_shell, numb_prim;

    numb_atom  = NUM2INT(_numb_atom);
    numb_shell = NUM2INT(_numb_shell);
    numb_prim  = NUM2INT(_numb_prim);

    /* Set angular type: 0 = Cartesian, 1 = Spherical */
    g_use_spherical = NUM2INT(_angular_type);

    int *shel_lqn = (int *)malloc(numb_shell * sizeof(int));
    int *shel_tem = (int *)malloc(numb_shell * sizeof(int));
    int *shel_atm = (int *)malloc(numb_shell * sizeof(int));
    int *shel_add = (int *)malloc(numb_shell * sizeof(int));
    double *atom_charg = (double *)malloc(numb_atom * sizeof(double));
    double *atom_xyz   = (double *)malloc(numb_atom * 3 * sizeof(double));
    double *prim_exp   = (double *)malloc(numb_prim * sizeof(double));
    double *prim_coe   = (double *)malloc(numb_prim * sizeof(double));

    for (i = 0; i < numb_shell; i++) {
        shel_lqn[i] = NUM2INT(*(RARRAY_PTR(_shel_lqn) + i));
        shel_tem[i] = NUM2INT(*(RARRAY_PTR(_shel_tem) + i));
        shel_atm[i] = NUM2INT(*(RARRAY_PTR(_shel_atm) + i));
        shel_add[i] = NUM2INT(*(RARRAY_PTR(_shel_add) + i));
    }
    for (i = 0; i < numb_atom; i++) {
        atom_charg[i] = NUM2DBL(*(RARRAY_PTR(_atom_charg) + i));
        for (int j = 0; j < 3; j++) {
            atom_xyz[i*3+j] = NUM2DBL(*(RARRAY_PTR(_atom_xyz) + i*3+j));
        }
    }
    for (i = 0; i < numb_prim; i++) {
        prim_exp[i] = NUM2DBL(*(RARRAY_PTR(_prim_exp) + i));
        prim_coe[i] = NUM2DBL(*(RARRAY_PTR(_prim_coe) + i));
    }

    build_libcint_data(numb_atom, numb_shell, numb_prim,
                       shel_lqn, shel_tem, shel_atm, shel_add,
                       atom_charg, atom_xyz, prim_exp, prim_coe);

    free(shel_lqn); free(shel_tem); free(shel_atm); free(shel_add);
    free(atom_charg); free(atom_xyz); free(prim_exp); free(prim_coe);

    return Qnil;
}

/* --- Overlap calc --- */
static VALUE sint_wrapper_calc(VALUE self, VALUE _ish, VALUE _jsh) {
    int ish = NUM2INT(_ish), jsh = NUM2INT(_jsh);
    int inttype, nsize;
    double dint[225];

    calc_1e(int1e_ovlp_cart, int1e_ovlp_sph, ish, jsh, &inttype, &nsize, dint);

    VALUE eris = rb_ary_new2(nsize);
    for (int i = 0; i < nsize; i++)
        rb_ary_store(eris, i, rb_float_new(dint[i]));

    VALUE results = rb_ary_new2(3);
    rb_ary_store(results, 0, INT2NUM(inttype));
    rb_ary_store(results, 1, INT2NUM(nsize));
    rb_ary_store(results, 2, eris);
    return results;
}

/* --- Kinetic calc --- */
static VALUE tint_wrapper_calc(VALUE self, VALUE _ish, VALUE _jsh) {
    int ish = NUM2INT(_ish), jsh = NUM2INT(_jsh);
    int inttype, nsize;
    double dint[225];

    calc_1e(int1e_kin_cart, int1e_kin_sph, ish, jsh, &inttype, &nsize, dint);

    VALUE eris = rb_ary_new2(nsize);
    for (int i = 0; i < nsize; i++)
        rb_ary_store(eris, i, rb_float_new(dint[i]));

    VALUE results = rb_ary_new2(3);
    rb_ary_store(results, 0, INT2NUM(inttype));
    rb_ary_store(results, 1, INT2NUM(nsize));
    rb_ary_store(results, 2, eris);
    return results;
}

/* --- Nuclear attraction calc --- */
static VALUE vint_wrapper_calc(VALUE self, VALUE _ish, VALUE _jsh) {
    int ish = NUM2INT(_ish), jsh = NUM2INT(_jsh);
    int inttype, nsize;
    double dint[225];

    calc_1e(int1e_nuc_cart, int1e_nuc_sph, ish, jsh, &inttype, &nsize, dint);

    VALUE eris = rb_ary_new2(nsize);
    for (int i = 0; i < nsize; i++)
        rb_ary_store(eris, i, rb_float_new(dint[i]));

    VALUE results = rb_ary_new2(3);
    rb_ary_store(results, 0, INT2NUM(inttype));
    rb_ary_store(results, 1, INT2NUM(nsize));
    rb_ary_store(results, 2, eris);
    return results;
}

/* --- ERI calc --- */
static VALUE eri0_wrapper_calc_lc(VALUE self,
    VALUE _ish, VALUE _jsh, VALUE _ksh, VALUE _lsh)
{
    int ish = NUM2INT(_ish), jsh = NUM2INT(_jsh);
    int ksh = NUM2INT(_ksh), lsh = NUM2INT(_lsh);
    int inttype, nsize;
    double eris_buf[10000];

    calc_2e(ish, jsh, ksh, lsh, &inttype, &nsize, eris_buf);

    VALUE eris = rb_ary_new2(nsize);
    for (int i = 0; i < nsize; i++)
        rb_ary_store(eris, i, rb_float_new(eris_buf[i]));

    VALUE results = rb_ary_new2(3);
    rb_ary_store(results, 0, INT2NUM(inttype));
    rb_ary_store(results, 1, INT2NUM(nsize));
    rb_ary_store(results, 2, eris);
    return results;
}

/* --- exchange_order --- */
static VALUE eri0_exchange_order_lc(VALUE self) {
    VALUE results = rb_ary_new2(4);
    rb_ary_store(results, 0, INT2NUM(g_last_a));
    rb_ary_store(results, 1, INT2NUM(g_last_b));
    rb_ary_store(results, 2, INT2NUM(g_last_c));
    rb_ary_store(results, 3, INT2NUM(g_last_d));
    return results;
}

/* ================================================================
 * Init_Molint: Ruby extension entry point
 *
 * init now takes 14 arguments (was 13):
 *   13 original args + angular_type (0=Cartesian, 1=Spherical)
 * ================================================================ */
void Init_Molint(void)
{
    VALUE cl_sint = rb_define_class("Sint0", rb_cObject);
    VALUE cl_tint = rb_define_class("Tint0", rb_cObject);
    VALUE cl_vint = rb_define_class("Vint0", rb_cObject);
    VALUE cl_eri0 = rb_define_class("Eri0",  rb_cObject);

    rb_define_method(cl_sint, "init", generic_wrapper_init, 14);
    rb_define_method(cl_sint, "calc", sint_wrapper_calc, 2);

    rb_define_method(cl_tint, "init", generic_wrapper_init, 14);
    rb_define_method(cl_tint, "calc", tint_wrapper_calc, 2);

    rb_define_method(cl_vint, "init", generic_wrapper_init, 14);
    rb_define_method(cl_vint, "calc", vint_wrapper_calc, 2);

    rb_define_method(cl_eri0, "init", generic_wrapper_init, 14);
    rb_define_method(cl_eri0, "calc", eri0_wrapper_calc_lc, 4);
    rb_define_method(cl_eri0, "exchange_order", eri0_exchange_order_lc, 0);
}
