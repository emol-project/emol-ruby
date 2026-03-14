/*
 * molint_libcint.c
 *
 * Ruby C extension that replaces the original OpenFMO-based Molint FFI
 * with libcint as the integral engine.
 *
 * Provides the same Ruby classes and interface:
 *   Sint0  : init(13 args), calc(ish, jsh)    -> [inttype, nsize, array]
 *   Tint0  : init(13 args), calc(ish, jsh)    -> [inttype, nsize, array]
 *   Vint0  : init(13 args), calc(ish, jsh)    -> [inttype, nsize, array]
 *   Eri0   : init(13 args), calc(ish,jsh,ksh,lsh) -> [inttype, nsize, array]
 *            exchange_order()                  -> [a, b, c, d]
 *
 * Integral convention: Cartesian (6d, 10f), GAMESS ordering and normalization.
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
#include "cint_funcs.h"

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

/* Cartesian component count: (l+1)*(l+2)/2 */
static int cart_comp(int l) {
    return (l + 1) * (l + 2) / 2;
}

/* ================================================================
 * GAMESS normalization factor for Cartesian component x^a y^b z^c
 *   factor = sqrt(l! / (a! b! c!))  where l = a+b+c
 *
 * For s,p this is always 1.0.
 * For d: xx,yy,zz -> 1.0;  xy,xz,yz -> sqrt(2)
 * For f: xxx,yyy,zzz -> 1.0; xxy,xxz,xyy,yyz,xzz,yzz -> sqrt(3); xyz -> sqrt(6)
 * ================================================================ */
static double gamess_cart_norm(int a, int b, int c) {
    /* Simple factorial for small numbers */
    static const double fact[] = {1,1,2,6,24,120,720,5040};
    int l = a + b + c;
    return sqrt(fact[l] / (fact[a] * fact[b] * fact[c]));
}

/* ================================================================
 * Cartesian component ordering tables
 *
 * libcint order: lexicographic (x highest first)
 *   d: xx,xy,xz,yy,yz,zz
 *   f: xxx,xxy,xxz,xyy,xyz,xzz,yyy,yyz,yzz,zzz
 *
 * GAMESS/emol order:
 *   d: xx,yy,zz,xy,xz,yz
 *   f: xxx,yyy,zzz,xxy,xxz,xyy,yyz,xzz,yzz,xyz
 *
 * We store: for each l, a permutation array and normalization factors.
 * ================================================================ */

/* (a,b,c) powers for libcint Cartesian ordering */
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

/* GAMESS/emol Cartesian ordering (hardcoded for l=0..4) */
static const int gamess_order_d[][3] = {
    {2,0,0},{0,2,0},{0,0,2},{1,1,0},{1,0,1},{0,1,1}
};
static const int gamess_order_f[][3] = {
    {3,0,0},{0,3,0},{0,0,3},{2,1,0},{2,0,1},{1,2,0},
    {0,2,1},{1,0,2},{0,1,2},{1,1,1}
};
static const int gamess_order_g[][3] = {
    {4,0,0},{0,4,0},{0,0,4},{3,1,0},{3,0,1},{0,3,1},
    {0,3,0},{1,0,3},{0,1,3},{2,2,0},{2,0,2},{0,2,2},
    {2,1,1},{1,2,1},{1,1,2}
};

/*
 * Build permutation: gamess_component[i] = libcint_component[perm[i]]
 * Also build GAMESS normalization factor array.
 */
static void build_cart_mapping(int l, int *perm, double *norm_fac) {
    int ncomp = cart_comp(l);
    CartPow *lc = (CartPow *)malloc(ncomp * sizeof(CartPow));
    libcint_cart_powers(l, lc);

    const int (*gamess_ord)[3];
    if (l == 0 || l == 1) {
        /* s, p: same ordering, norm = 1.0 */
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
        /* find this component in libcint ordering */
        for (int j = 0; j < ncomp; j++) {
            if (lc[j].a == ga && lc[j].b == gb && lc[j].c == gc) {
                perm[i] = j;
                break;
            }
        }
    }
    free(lc);
}

/* Pre-computed mapping tables (up to l=4) */
#define MAX_L 4
static int  g_perm[MAX_L+1][15+1];      /* max 15 components for g */
static double g_norm[MAX_L+1][15+1];
static int g_mapping_ready = 0;

static void ensure_mapping_tables(void) {
    if (g_mapping_ready) return;
    for (int l = 0; l <= MAX_L; l++) {
        build_cart_mapping(l, g_perm[l], g_norm[l]);
    }
    g_mapping_ready = 1;
}

/* ================================================================
 * gto_norm: PySCF/libcint normalization for radial part
 *   N = 1/sqrt( Gamma((2l+3)/2) / (2 * (2*alpha)^((2l+3)/2)) )
 * ================================================================ */
static double gto_norm(int l, double alpha) {
    /* gaussian_int(n, a) = Gamma((n+1)/2) / (2 * a^((n+1)/2)) */
    /* gto_norm = 1/sqrt(gaussian_int(2l+2, 2*alpha)) */
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
    /* Layout: [PTR_ENV_START ... coords ... exps ... coefs ...] */
    int env_ptr = PTR_ENV_START;

    /* Calculate total size needed */
    int total_env = env_ptr;
    total_env += numb_atom * 4;   /* coords(3) + zeta(1) per atom */
    total_env += numb_prim * 2;   /* exponents + coefficients */
    total_env += 100;             /* safety margin */

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
        g_atm[i * ATM_SLOTS + NUC_MOD_OF] = 1; /* point nucleus */
        g_atm[i * ATM_SLOTS + PTR_ZETA] = env_ptr;
        g_env[env_ptr] = 0.0;
        env_ptr += 1;
    }

    /* --- Build bas array --- */
    g_bas = (FINT *)calloc(numb_shell * BAS_SLOTS, sizeof(FINT));

    /*
     * emol stores prim_coe as already-normalized coefficients.
     * We need to convert to libcint convention:
     *   coef_libcint = original_coef * gto_norm(l, alpha)
     *
     * emol's prim_coe = original_coef * emol_full_norm(l, alpha)
     * emol_full_norm = gto_norm * sqrt((2l+1)/(4pi))
     *
     * So: coef_libcint = prim_coe / sqrt((2l+1)/(4pi))
     *
     * But actually, since we apply GAMESS normalization in the
     * post-processing step (per-component factor sqrt(l!/(a!b!c!))),
     * we need to make sure the combination is correct.
     *
     * Cleanest approach: use original XML coefficients.
     * But we don't have them here — we only have prim_coe.
     *
     * Conversion: coef_for_libcint = prim_coe / sqrt((2l+1)/(4pi))
     */

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
        bf_count += cart_comp(l);

        g_bas[ish * BAS_SLOTS + ATOM_OF]   = iatom;
        g_bas[ish * BAS_SLOTS + ANG_OF]    = l;
        g_bas[ish * BAS_SLOTS + NPRIM_OF]  = nprim;
        g_bas[ish * BAS_SLOTS + NCTR_OF]   = 1;  /* 1 contraction */
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

        /* Contraction renormalization (same as PySCF's _normalize_contracted_ao).
         * Compute the self-overlap of this contracted shell using libcint's
         * normalization, then divide all coefficients by sqrt(overlap).
         * This ensures <chi|chi> = 1 for the basic Cartesian component.
         */
        {
            double *es = &prim_exp[prim_offset];
            double *cs = &g_env[env_ptr];
            double ovlp = 0.0;
            for (int i = 0; i < nprim; i++) {
                for (int j = 0; j < nprim; j++) {
                    /* 3D overlap of x^l*exp(-ai*r^2) with x^l*exp(-aj*r^2):
                     * Uses the fact that libcint's cart integrals include
                     * per-component normalization internally.
                     * The overlap kernel is: (pi/(ai+aj))^1.5 * dfact(2l-1) / (2(ai+aj))^l
                     * But since coefficients already contain gto_norm, we compute:
                     * S = sum_ij ci*cj * gaussian_int(2l+2, ai+aj) */
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
    ensure_mapping_tables();
}

/* ================================================================
 * 1-electron integral: compute shell pair (ish, jsh)
 *
 * Returns integral values in GAMESS order with GAMESS normalization.
 * The output array is in column-major order matching the original FFI:
 *   dint[jp * ncomp_i + ip]  (ip fast, jp slow)
 * ================================================================ */
static int calc_1e(
    CINTIntegralFunction *int1e_func,
    int ish, int jsh,
    int *out_inttype, int *out_nsize, double *out_dint)
{
    int li = g_shel_lqn[ish];
    int lj = g_shel_lqn[jsh];
    int nci = cart_comp(li);
    int ncj = cart_comp(lj);

    *out_inttype = (li + 1) * 10 + (lj + 1);
    *out_nsize = nci * ncj;

    /* libcint buffer */
    double *buf = (double *)calloc(nci * ncj, sizeof(double));
    FINT shls[2] = { ish, jsh };

    /* Call libcint (CINTIntegralFunction signature: 10 args with dims and cache) */
    CACHE_SIZE_T not0 = int1e_func(buf, NULL, shls, g_atm, g_natm, g_bas, g_nbas, g_env, NULL, NULL);

    if (!not0) {
        /* All zero */
        memset(out_dint, 0, nci * ncj * sizeof(double));
        free(buf);
        return 0;
    }

    /*
     * libcint _cart output: column-major [i_fast, j_slow] in libcint order.
     * buf[j_lc * nci + i_lc]
     *
     * We need to reorder to GAMESS order.
     * libcint _cart already normalizes each Cartesian component (GAMESS style),
     * so we only need the permutation, not the normalization factors.
     */
    int *perm_i = g_perm[li];
    int *perm_j = g_perm[lj];

    for (int jp = 0; jp < ncj; jp++) {
        for (int ip = 0; ip < nci; ip++) {
            int lc_idx = perm_j[jp] * nci + perm_i[ip];
            out_dint[jp * nci + ip] = buf[lc_idx];
        }
    }

    free(buf);
    return 0;
}

/* ================================================================
 * 2-electron integral: compute shell quartet (ish, jsh, ksh, lsh)
 *
 * The original FFI does shell reordering internally via
 * eri0_exchange_to_atomic_order and returns the reordered shell indices
 * via exchange_order(). We mimic this.
 *
 * Output: flat array in order ((ip*ncj+jp)*nck+kp)*ncl+lp
 * with GAMESS normalization and ordering.
 * ================================================================ */

/* Store last exchange order for exchange_order() method */
static int g_last_a = 0, g_last_b = 0, g_last_c = 0, g_last_d = 0;

static int calc_2e(int ish, int jsh, int ksh, int lsh,
                   int *out_inttype, int *out_nsize, double *out_eris)
{
    /*
     * Shell reordering: libcint requires no specific ordering,
     * but the original FFI reorders shells so that La >= Lb, Lc >= Ld,
     * and the (ij) pair has La+Lb >= Lc+Ld.
     * We keep the same convention for exchange_order() compatibility.
     */
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
    int nca = cart_comp(la), ncb = cart_comp(lb);
    int ncc = cart_comp(lc), ncd = cart_comp(ld);
    int ntot = nca * ncb * ncc * ncd;

    *out_inttype = ((la+1)*10 + (lb+1))*100 + (lc+1)*10 + (ld+1);
    *out_nsize = ntot;

    /* libcint buffer */
    double *buf = (double *)calloc(ntot, sizeof(double));
    FINT shls[4] = { a, b, c, d };

    CACHE_SIZE_T not0 = int2e_cart(buf, NULL, shls, g_atm, g_natm, g_bas, g_nbas, g_env, NULL, NULL);

    if (!not0) {
        memset(out_eris, 0, ntot * sizeof(double));
        free(buf);
        return 0;
    }

    /*
     * libcint 2e output: column-major [i,j,k,l] in libcint Cartesian order.
     * buf[ ((l_lc * nck + k_lc) * ncb + j_lc) * nca + i_lc ]
     *   where i is fastest.
     *
     * Original FFI output: ((ip*ncj+jp)*nck+kp)*ncl+lp
     *   where ip corresponds to shell a (fastest in emol convention).
     *
     * Actually looking at emol_int.rb mk_eris:
     *   nn = ((ip * ncomp_j + jp) * ncomp_k + kp) * ncomp_l + lp
     * So the order is ip fastest... but wait, in the original FFI calc,
     * the flat array just stores nsize values sequentially.
     * Let's match: the original stores as
     *   eris[((ip*ncb+jp)*ncc+kp)*ncd+lp]
     * ip = component of shell a, jp = shell b, kp = shell c, lp = shell d.
     */

    int *perm_a = g_perm[la], *perm_b = g_perm[lb];
    int *perm_c = g_perm[lc], *perm_d = g_perm[ld];

    for (int ip = 0; ip < nca; ip++) {
        for (int jp = 0; jp < ncb; jp++) {
            for (int kp = 0; kp < ncc; kp++) {
                for (int lp = 0; lp < ncd; lp++) {
                    /* libcint index: i fastest */
                    int lc_i = perm_a[ip], lc_j = perm_b[jp];
                    int lc_k = perm_c[kp], lc_l = perm_d[lp];
                    int lc_idx = ((lc_l * ncc + lc_k) * ncb + lc_j) * nca + lc_i;

                    /* emol output index */
                    int em_idx = ((ip * ncb + jp) * ncc + kp) * ncd + lp;

                    out_eris[em_idx] = buf[lc_idx];
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

/* --- Generic init wrapper (same for all 4 classes) --- */
static VALUE generic_wrapper_init(VALUE self,
    VALUE _numb_atom, VALUE _numb_shell, VALUE _numb_prim,
    VALUE _shel_lqn, VALUE _shel_tem, VALUE _shel_atm, VALUE _shel_add,
    VALUE _atom_charg, VALUE _atom_xyz, VALUE _prim_exp, VALUE _prim_coe,
    VALUE _tol, VALUE _level_dbg_print)
{
    int i, numb_atom, numb_shell, numb_prim;

    numb_atom  = NUM2INT(_numb_atom);
    numb_shell = NUM2INT(_numb_shell);
    numb_prim  = NUM2INT(_numb_prim);

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
    double dint[225]; /* max for l=4: 15*15 = 225 */

    calc_1e(int1e_ovlp_cart, ish, jsh, &inttype, &nsize, dint);

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

    calc_1e(int1e_kin_cart, ish, jsh, &inttype, &nsize, dint);

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

    calc_1e(int1e_nuc_cart, ish, jsh, &inttype, &nsize, dint);

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
    double eris_buf[10000]; /* _N_ERI_ from original */

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
 * Defines exactly the same classes and methods as the original FFI.
 * ================================================================ */
void Init_Molint(void)
{
    VALUE cl_sint = rb_define_class("Sint0", rb_cObject);
    VALUE cl_tint = rb_define_class("Tint0", rb_cObject);
    VALUE cl_vint = rb_define_class("Vint0", rb_cObject);
    VALUE cl_eri0 = rb_define_class("Eri0",  rb_cObject);

    /* All classes share the same init (builds libcint data once) */
    rb_define_method(cl_sint, "init", generic_wrapper_init, 13);
    rb_define_method(cl_sint, "calc", sint_wrapper_calc, 2);

    rb_define_method(cl_tint, "init", generic_wrapper_init, 13);
    rb_define_method(cl_tint, "calc", tint_wrapper_calc, 2);

    rb_define_method(cl_vint, "init", generic_wrapper_init, 13);
    rb_define_method(cl_vint, "calc", vint_wrapper_calc, 2);

    rb_define_method(cl_eri0, "init", generic_wrapper_init, 13);
    rb_define_method(cl_eri0, "calc", eri0_wrapper_calc_lc, 4);
    rb_define_method(cl_eri0, "exchange_order", eri0_exchange_order_lc, 0);
}
