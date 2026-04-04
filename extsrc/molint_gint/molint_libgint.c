/*
 * ==================================================================
 * molint_gint.c
 * ==================================================================
 * Ruby C extension providing the same interface as molint_libcint.c,
 * but using libgint (GINT Fortran engine) as the integral backend.
 *
 * This file contains ONLY Ruby glue code.
 * All integral logic lives in gint_mol.c (compiled into libgint.a).
 *
 * Ruby interface (identical to molint_libcint.c):
 *   Sint0 : init(14 args), calc(ish, jsh)         -> [inttype, nsize, array]
 *   Tint0 : init(14 args), calc(ish, jsh)         -> [inttype, nsize, array]
 *   Vint0 : init(14 args), calc(ish, jsh)         -> [inttype, nsize, array]
 *   Eri0  : init(14 args), calc(ish,jsh,ksh,lsh)  -> [inttype, nsize, array]
 *           exchange_order()                       -> [a, b, c, d]
 *
 * angular_type (14th init arg):
 *   0 -> Cartesian (GAMESS ordering and normalization)
 *   1 -> Spherical (real spherical harmonics, libcint/PySCF ordering)
 *
 * Cartesian output is reordered from GINT order to GAMESS order and
 * GAMESS normalization factors are applied (same as molint_libcint.c).
 * ==================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ruby.h"
#include "fint_capi.h"

/* ================================================================
 * Cartesian component ordering and GAMESS normalization
 *
 * GINT/libcint order: descending nx, for each nx descending ny
 *   d: xx,xy,xz,yy,yz,zz
 *   f: xxx,xxy,xxz,xyy,xyz,xzz,yyy,yyz,yzz,zzz
 *   g: xxxx,xxxy,...,zzzz
 *
 * GAMESS/emol order:
 *   d: xx,yy,zz,xy,xz,yz
 *   f: xxx,yyy,zzz,xxy,xxz,xyy,yyz,xzz,yzz,xyz
 *   g: xxxx,yyyy,zzzz,xxxy,xxxz,xyyy,yyyz,xzzz,yzzz,
 *      xxyy,xxzz,yyzz,xxyz,xyyz,xyzz
 *
 * GAMESS norm factor for x^a y^b z^c: sqrt(l! / (a! b! c!))
 * ================================================================ */

typedef struct { int a, b, c; } CartPow;

static void gint_cart_powers(int l, CartPow *out) {
    int idx = 0;
    for (int ax = l; ax >= 0; ax--)
        for (int ay = l - ax; ay >= 0; ay--) {
            out[idx].a = ax;
            out[idx].b = ay;
            out[idx].c = l - ax - ay;
            idx++;
        }
}

static double gamess_cart_norm(int a, int b, int c) {
    static const double fact[] = {1,1,2,6,24,120,720,5040};
    int l = a + b + c;
    return sqrt(fact[l] / (fact[a] * fact[b] * fact[c]));
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

#define MAX_L 4
static int    g_perm[MAX_L+1][15+1];
static double g_norm[MAX_L+1][15+1];
static int    g_cart_tables_ready = 0;

static void build_cart_tables(void) {
    if (g_cart_tables_ready) return;
    for (int l = 0; l <= MAX_L; l++) {
        int ncomp = (l+1)*(l+2)/2;
        CartPow lc[15];
        gint_cart_powers(l, lc);
        if (l <= 1) {
            for (int i = 0; i < ncomp; i++) {
                g_perm[l][i] = i;
                g_norm[l][i] = 1.0;
            }
            continue;
        }
        const int (*go)[3] = (l==2) ? gamess_order_d :
                             (l==3) ? gamess_order_f : gamess_order_g;
        for (int i = 0; i < ncomp; i++) {
            int ga = go[i][0], gb = go[i][1], gc = go[i][2];
            g_norm[l][i] = gamess_cart_norm(ga, gb, gc);
            for (int j = 0; j < ncomp; j++)
                if (lc[j].a == ga && lc[j].b == gb && lc[j].c == gc) {
                    g_perm[l][i] = j; break;
                }
        }
    }
    g_cart_tables_ready = 1;
}

static void apply_cart_reorder_1e(
    double *dst, const double *src,
    int li, int lj, int nci, int ncj)
{
    int    *pi = g_perm[li], *pj = g_perm[lj];
    double *ni = g_norm[li], *nj = g_norm[lj];
    for (int jp = 0; jp < ncj; jp++)
        for (int ip = 0; ip < nci; ip++)
            dst[jp * nci + ip] =
                src[pj[jp] * nci + pi[ip]] * ni[ip] * nj[jp];
}

static void apply_cart_reorder_2e(
    double *out, const double *src,
    int la, int lb, int lc, int ld,
    int nca, int ncb, int ncc, int ncd)
{
    int    *pa = g_perm[la], *pb = g_perm[lb];
    int    *pc = g_perm[lc], *pd = g_perm[ld];
    double *na = g_norm[la], *nb = g_norm[lb];
    double *nc_= g_norm[lc], *nd = g_norm[ld];

    for (int ip=0; ip<nca; ip++)
    for (int jp=0; jp<ncb; jp++)
    for (int kp=0; kp<ncc; kp++)
    for (int lp=0; lp<ncd; lp++) {
        int gi = pa[ip], gj = pb[jp], gk = pc[kp], gl = pd[lp];
        int gint_idx = ((gl*ncc+gk)*ncb+gj)*nca+gi;
        int em_idx   = ((ip*ncb+jp)*ncc+kp)*ncd+lp;
        out[em_idx]  = src[gint_idx] * na[ip]*nb[jp]*nc_[kp]*nd[lp];
    }
}

/* ================================================================
 * Global state
 * ================================================================ */
static int g_use_spherical = 0;
static int g_last_a = 0, g_last_b = 0, g_last_c = 0, g_last_d = 0;

/* ================================================================
 * Ruby: generic init
 * ================================================================ */
static VALUE generic_wrapper_init(VALUE self,
    VALUE _numb_atom,  VALUE _numb_shell, VALUE _numb_prim,
    VALUE _shel_lqn,   VALUE _shel_tem,   VALUE _shel_atm, VALUE _shel_add,
    VALUE _atom_charg, VALUE _atom_xyz,
    VALUE _prim_exp,   VALUE _prim_coe,
    VALUE _tol,        VALUE _level_dbg_print, VALUE _angular_type)
{
    int numb_atom  = NUM2INT(_numb_atom);
    int numb_shell = NUM2INT(_numb_shell);
    int numb_prim  = NUM2INT(_numb_prim);
    int angular    = NUM2INT(_angular_type);

    g_use_spherical = angular;

    int    *shel_lqn   = (int    *)malloc(numb_shell * sizeof(int));
    int    *shel_tem   = (int    *)malloc(numb_shell * sizeof(int));
    int    *shel_atm   = (int    *)malloc(numb_shell * sizeof(int));
    int    *shel_add   = (int    *)malloc(numb_shell * sizeof(int));
    double *atom_charg = (double *)malloc(numb_atom  * sizeof(double));
    double *atom_xyz   = (double *)malloc(numb_atom  * 3 * sizeof(double));
    double *prim_exp   = (double *)malloc(numb_prim  * sizeof(double));
    double *prim_coe   = (double *)malloc(numb_prim  * sizeof(double));

    for (int i = 0; i < numb_shell; i++) {
        shel_lqn[i] = NUM2INT(*(RARRAY_PTR(_shel_lqn) + i));
        shel_tem[i] = NUM2INT(*(RARRAY_PTR(_shel_tem) + i));
        shel_atm[i] = NUM2INT(*(RARRAY_PTR(_shel_atm) + i));
        shel_add[i] = NUM2INT(*(RARRAY_PTR(_shel_add) + i));
    }
    for (int i = 0; i < numb_atom; i++) {
        atom_charg[i] = NUM2DBL(*(RARRAY_PTR(_atom_charg) + i));
        for (int j = 0; j < 3; j++)
            atom_xyz[i*3+j] = NUM2DBL(*(RARRAY_PTR(_atom_xyz) + i*3+j));
    }
    for (int i = 0; i < numb_prim; i++) {
        prim_exp[i] = NUM2DBL(*(RARRAY_PTR(_prim_exp) + i));
        prim_coe[i] = NUM2DBL(*(RARRAY_PTR(_prim_coe) + i));
    }

    gint_mol_init(numb_atom, numb_shell, numb_prim,
                  shel_lqn, shel_tem, shel_atm, shel_add,
                  atom_charg, atom_xyz, prim_exp, prim_coe,
                  angular);

    if (!g_use_spherical) build_cart_tables();

    free(shel_lqn); free(shel_tem); free(shel_atm); free(shel_add);
    free(atom_charg); free(atom_xyz); free(prim_exp); free(prim_coe);
    return Qnil;
}

/* ================================================================
 * make_result
 * ================================================================ */
static VALUE make_result(int inttype, int nsize, double *dint) {
    VALUE arr = rb_ary_new2(nsize);
    for (int i = 0; i < nsize; i++)
        rb_ary_store(arr, i, rb_float_new(dint[i]));
    VALUE res = rb_ary_new2(3);
    rb_ary_store(res, 0, INT2NUM(inttype));
    rb_ary_store(res, 1, INT2NUM(nsize));
    rb_ary_store(res, 2, arr);
    return res;
}

/* ================================================================
 * 1E calc
 * ================================================================ */
typedef enum { WHICH_S, WHICH_T, WHICH_V } Which1E;

static VALUE calc_1e(VALUE _ish, VALUE _jsh, Which1E which) {
    int ish = NUM2INT(_ish);
    int jsh = NUM2INT(_jsh);
    int li  = gint_mol_shell_lqn(ish);
    int lj  = gint_mol_shell_lqn(jsh);
    int nci = gint_mol_shell_comp(ish);
    int ncj = gint_mol_shell_comp(jsh);
    int n   = nci * ncj;
    int inttype = (li+1)*10 + (lj+1);

    double s_raw[GINT_MAX_CART], t_raw[GINT_MAX_CART], v_raw[GINT_MAX_CART];
    gint_mol_1e_stv(s_raw, t_raw, v_raw, ish, jsh);

    double *raw = (which == WHICH_S) ? s_raw :
                  (which == WHICH_T) ? t_raw : v_raw;
    double out[GINT_MAX_CART];

    if (g_use_spherical) {
        memcpy(out, raw, n * sizeof(double));
    } else {
        apply_cart_reorder_1e(out, raw, li, lj, nci, ncj);
    }
    return make_result(inttype, n, out);
}

static VALUE sint_calc(VALUE self, VALUE ish, VALUE jsh) {
    return calc_1e(ish, jsh, WHICH_S);
}
static VALUE tint_calc(VALUE self, VALUE ish, VALUE jsh) {
    return calc_1e(ish, jsh, WHICH_T);
}
static VALUE vint_calc(VALUE self, VALUE ish, VALUE jsh) {
    return calc_1e(ish, jsh, WHICH_V);
}

/* Vint0#calc_raw(ish, jsh) -> [inttype, nsize, array]
 * array: buf[k*nci*ncj + jp*nci + ip], k=nucleus(0-based)
 * nsize = nci * ncj * natm
 */
static VALUE vint_calc_raw(VALUE self, VALUE _ish, VALUE _jsh) {
    int ish = NUM2INT(_ish);
    int jsh = NUM2INT(_jsh);
    int li  = gint_mol_shell_lqn(ish);
    int lj  = gint_mol_shell_lqn(jsh);
    int nci = gint_mol_shell_comp(ish);
    int ncj = gint_mol_shell_comp(jsh);
    int nsize = nci * ncj * gint_mol_natm();  /* g_natm -> gint_mol_natm() */
    int inttype = (li+1)*10 + (lj+1);

    double *buf = (double *)calloc(nsize, sizeof(double));
    gint_mol_1e_vnuc_raw(buf, ish, jsh);

    VALUE result = make_result(inttype, nsize, buf);
    free(buf);
    return result;
}


/* ================================================================
 * 2E calc
 * ================================================================ */
static VALUE eri0_calc(VALUE self,
    VALUE _ish, VALUE _jsh, VALUE _ksh, VALUE _lsh)
{
    int a = NUM2INT(_ish), b = NUM2INT(_jsh);
    int c = NUM2INT(_ksh), d = NUM2INT(_lsh);

    if (gint_mol_shell_lqn(a) < gint_mol_shell_lqn(b)) { int t=a; a=b; b=t; }
    if (gint_mol_shell_lqn(c) < gint_mol_shell_lqn(d)) { int t=c; c=d; d=t; }
    if (gint_mol_shell_lqn(a)+gint_mol_shell_lqn(b) <
        gint_mol_shell_lqn(c)+gint_mol_shell_lqn(d)) {
        int t; t=a;a=c;c=t; t=b;b=d;d=t;
    }
    g_last_a=a; g_last_b=b; g_last_c=c; g_last_d=d;

    int la=gint_mol_shell_lqn(a), lb=gint_mol_shell_lqn(b);
    int lc=gint_mol_shell_lqn(c), ld=gint_mol_shell_lqn(d);
    int nca=gint_mol_shell_comp(a), ncb=gint_mol_shell_comp(b);
    int ncc=gint_mol_shell_comp(c), ncd=gint_mol_shell_comp(d);
    int ntot = nca*ncb*ncc*ncd;
    int inttype = ((la+1)*10+(lb+1))*100+(lc+1)*10+(ld+1);

    double *buf = (double *)malloc(ntot * sizeof(double));
    double *out = (double *)malloc(ntot * sizeof(double));

    gint_mol_2e(buf, a, b, c, d);

    if (g_use_spherical) {
        for (int ip=0;ip<nca;ip++)
        for (int jp=0;jp<ncb;jp++)
        for (int kp=0;kp<ncc;kp++)
        for (int lp=0;lp<ncd;lp++) {
            out[((ip*ncb+jp)*ncc+kp)*ncd+lp] =
                buf[((lp*ncc+kp)*ncb+jp)*nca+ip];
        }
    } else {
        apply_cart_reorder_2e(out, buf, la, lb, lc, ld, nca, ncb, ncc, ncd);
    }

    free(buf);
    VALUE result = make_result(inttype, ntot, out);
    free(out);
    return result;
}

static VALUE eri0_exchange_order(VALUE self) {
    VALUE res = rb_ary_new2(4);
    rb_ary_store(res, 0, INT2NUM(g_last_a));
    rb_ary_store(res, 1, INT2NUM(g_last_b));
    rb_ary_store(res, 2, INT2NUM(g_last_c));
    rb_ary_store(res, 3, INT2NUM(g_last_d));
    return res;
}

/* ================================================================
 * Init_Molint
 * ================================================================ */
void Init_Molint(void) {
    VALUE cl_sint = rb_define_class("Sint0", rb_cObject);
    VALUE cl_tint = rb_define_class("Tint0", rb_cObject);
    VALUE cl_vint = rb_define_class("Vint0", rb_cObject);
    VALUE cl_eri0 = rb_define_class("Eri0",  rb_cObject);

    rb_define_method(cl_sint, "init", generic_wrapper_init, 14);
    rb_define_method(cl_sint, "calc", sint_calc,            2);

    rb_define_method(cl_tint, "init", generic_wrapper_init, 14);
    rb_define_method(cl_tint, "calc", tint_calc,            2);

    rb_define_method(cl_vint, "init", generic_wrapper_init, 14);
    rb_define_method(cl_vint, "calc", vint_calc,            2);
    rb_define_method(cl_vint, "calc_raw", vint_calc_raw,    2);

    rb_define_method(cl_eri0, "init",           generic_wrapper_init, 14);
    rb_define_method(cl_eri0, "calc",           eri0_calc,            4);
    rb_define_method(cl_eri0, "exchange_order", eri0_exchange_order,  0);
}
