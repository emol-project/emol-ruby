#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double MOE_LD;
typedef long   MOE_LI;

//==============================================================================
static MOE_LD FMT_DZERO = 0.0e+0;
static MOE_LD FMT_DONE  = 1.0e+0;
static MOE_LD FMT_DTWO  = 2.0e+0;
static MOE_LD DHALF     = 0.5e+0;
static MOE_LD THRZERO   = 1.0e-17;
static MOE_LD DTHIRD, TWO_PI;
static MOE_LD DINV31, DINV29;
static MOE_LD DINV27, DINV25, DINV23, DINV21, DINV19;
static MOE_LD DINV17, DINV15, DINV13, DINV11, DINV09;
static MOE_LD DINV07, DINV05, DINV03;
static MOE_LI IONE = 1L;
static MOE_LI ITWO = 2L;

// 誤差関数表作成に必要な定数
static MOE_LI FMT_LAST_M_VALUE    = 17L;
static MOE_LI FMT_N_AREA          = 3L;
static MOE_LI FMT_INIT_T_VALUES[] = { 0L, 22L, 45L };
static MOE_LI FMT_LAST_T_VALUES[] = { 22L, 45L, 90L };
static MOE_LI FMT_N_STEPS[]       = { 1024L, 512L, 256L };

//==============================================================================
static MOE_LD *FMT_TABLE = NULL;
static MOE_LD *FMT_EXP_TABLE = NULL;
static MOE_LD *FMT_STORE = NULL;
static MOE_LD FMT_1_SIZE, FMT_2_SIZE, FMT_3_SIZE;
static MOE_LD FMT_1_LASTT, FMT_2_LASTT, FMT_3_LASTT;
static MOE_LI FMT_TABLE_SIZE, FMT_EXP_TABLE_SIZE, FMT_M_SIZE;
static MOE_LI FMT_1_STEPS, FMT_2_STEPS, FMT_3_STEPS;
static MOE_LI FMT_1_GRIDS, FMT_12_GRIDS;
static FILE   *FMT_OUTFILE;

//==============================================================================
//==============================================================================
//==============================================================================
MOE_LI FMT_free_table() {
    if (FMT_TABLE != NULL) { free(FMT_TABLE); FMT_TABLE = NULL; }
    if (FMT_EXP_TABLE != NULL) { free(FMT_EXP_TABLE); FMT_EXP_TABLE = NULL; }
    return 0;
}

MOE_LI FMT_set_store_area(MOE_LD *Mem) {
    FMT_STORE = Mem;
    return 0;
}

MOE_LD FMT_value_slow(MOE_LI Mvalue, MOE_LD Tvalue) {
    MOE_LD tvalue2, expminust, thrfmt, oddnumber, term;
    MOE_LD fmt = FMT_DZERO;
    MOE_LI nowmvalue;
    if ((Tvalue < FMT_DZERO) || (Mvalue < 0)) return FMT_DZERO;
    if (Tvalue == FMT_DZERO) {
        fmt = FMT_DONE / (ITWO * Mvalue + IONE);
    } else {
        tvalue2   = Tvalue + Tvalue;
        expminust = exp(-Tvalue);
        thrfmt    = (expminust / tvalue2) * THRZERO;
        nowmvalue = Tvalue + 0.4;
        if (Mvalue > nowmvalue) nowmvalue = Mvalue;
        nowmvalue += IONE;
        oddnumber = nowmvalue + nowmvalue + FMT_DONE;
        term      = tvalue2 / oddnumber;
        fmt       = FMT_DZERO;
        do {
            fmt       += term;
            oddnumber += FMT_DTWO;
            term      *= tvalue2 / oddnumber;
        } while ((term / fmt) >= thrfmt);
        fmt = (fmt + term) * (expminust / tvalue2);
        while (nowmvalue-- > Mvalue) {
            fmt = (fmt * tvalue2 + expminust) / (ITWO * nowmvalue + IONE);
        }
    }
    return fmt;
}

//==============================================================================
MOE_LI FMT_make_table(FILE *Outfile) {
    MOE_LD *FmtTable, *ExpTable;
    MOE_LD Steplength, Tvalue;
    MOE_LI nowcase, NowSteps, nowTvalue, i, j;

    FMT_M_SIZE    = FMT_LAST_M_VALUE + IONE; 
    FMT_1_STEPS   = FMT_N_STEPS[0];
    FMT_2_STEPS   = FMT_N_STEPS[1];
    FMT_3_STEPS   = FMT_N_STEPS[2];
    FMT_1_SIZE    = FMT_DONE / FMT_1_STEPS;
    FMT_2_SIZE    = FMT_DONE / FMT_2_STEPS;
    FMT_3_SIZE    = FMT_DONE / FMT_3_STEPS;
    FMT_1_LASTT   = FMT_LAST_T_VALUES[0];
    FMT_2_LASTT   = FMT_LAST_T_VALUES[1];
    FMT_3_LASTT   = FMT_LAST_T_VALUES[2];
    FMT_1_GRIDS   = (FMT_1_LASTT - FMT_INIT_T_VALUES[0]) * FMT_1_STEPS;
    FMT_12_GRIDS  = (FMT_2_LASTT - FMT_INIT_T_VALUES[1]) * FMT_2_STEPS + FMT_1_GRIDS;
    TWO_PI        = 8.0 * atan(FMT_DONE);
    
    //==============================================================================
    DTHIRD = FMT_DONE / 3.0;
    DINV31 = FMT_DONE / 31.0;
    DINV29 = FMT_DONE / 29.0;
    DINV27 = FMT_DONE / 27.0;
    DINV25 = FMT_DONE / 25.0;  DINV23 = FMT_DONE / 23.0;
    DINV21 = FMT_DONE / 21.0;  DINV19 = FMT_DONE / 19.0;
    DINV17 = FMT_DONE / 17.0;  DINV15 = FMT_DONE / 15.0;
    DINV13 = FMT_DONE / 13.0;  DINV11 = FMT_DONE / 11.0;
    DINV09 = FMT_DONE /  9.0;  DINV07 = FMT_DONE /  7.0;
    DINV05 = FMT_DONE /  5.0;  DINV03 = FMT_DONE /  3.0;
    FMT_OUTFILE = Outfile;

    FMT_EXP_TABLE_SIZE = 1 
        + FMT_N_STEPS[0] * (FMT_LAST_T_VALUES[0] - FMT_INIT_T_VALUES[0])
        + FMT_N_STEPS[1] * (FMT_LAST_T_VALUES[1] - FMT_INIT_T_VALUES[1])
        + FMT_N_STEPS[2] * (FMT_LAST_T_VALUES[2] - FMT_INIT_T_VALUES[2]);
    FMT_TABLE_SIZE = (FMT_LAST_M_VALUE + IONE) * FMT_EXP_TABLE_SIZE;

    FMT_EXP_TABLE = (MOE_LD *)malloc(FMT_EXP_TABLE_SIZE * sizeof(MOE_LD));
    if (!FMT_EXP_TABLE) { if(Outfile) fprintf(Outfile, "Mem-alloc-error\n"); exit(10); }
    FMT_TABLE = (MOE_LD *)malloc(FMT_TABLE_SIZE * sizeof(MOE_LD));
    if (!FMT_TABLE) { if(Outfile) fprintf(Outfile, "Mem-alloc-error\n"); exit(20); }

    FmtTable = FMT_TABLE;
    ExpTable = FMT_EXP_TABLE;

    for (i = 0; i <= FMT_LAST_M_VALUE; i++) *FmtTable++ = FMT_DONE / (ITWO * i + IONE);
    *ExpTable++ = FMT_DONE;

    for (nowcase = 0; nowcase < FMT_N_AREA; nowcase++) {
        NowSteps = FMT_N_STEPS[nowcase];
        Steplength = FMT_DONE / NowSteps;
        for (nowTvalue = FMT_INIT_T_VALUES[nowcase]; nowTvalue < FMT_LAST_T_VALUES[nowcase]; nowTvalue++) {
            for (i = 1; i <= NowSteps; i++) {
                Tvalue = nowTvalue + i * Steplength;
                *ExpTable++ = exp(-Tvalue);
                for (j = 0; j <= FMT_LAST_M_VALUE; j++) {
                    *FmtTable++ = FMT_value_slow(j, Tvalue);
                }
            }
        }
    }
    return 0;
}

//==============================================================================
MOE_LI FMT_value_fast(MOE_LI Mvalue, MOE_LD Tvalue, MOE_LD Const) {
    MOE_LD *fmt, *tbl, *exptbl, evl, tv2, del, t_inv, tm;
    MOE_LI i, k;

    //==============================================================================
    if (Mvalue > 16) return 1;

    if (Tvalue <= FMT_1_LASTT) {
        k = Tvalue * FMT_1_STEPS + DHALF;
        del = k * FMT_1_SIZE - Tvalue;
        tbl = FMT_TABLE + k * FMT_M_SIZE;
        if (Mvalue > 1) exptbl = FMT_EXP_TABLE + k;
    } else if (Tvalue <= FMT_2_LASTT) {
        k = (Tvalue - FMT_1_LASTT) * FMT_2_STEPS + DHALF;
        del = k * FMT_2_SIZE + FMT_1_LASTT - Tvalue;
        tbl = FMT_TABLE + (k + FMT_1_GRIDS) * FMT_M_SIZE;
        if (Mvalue > 1) exptbl = FMT_EXP_TABLE + k + FMT_1_GRIDS;
    } else if (Tvalue <= FMT_3_LASTT) {
        k = (Tvalue - FMT_2_LASTT) * FMT_3_STEPS + DHALF;
        del = k * FMT_3_SIZE + FMT_2_LASTT - Tvalue;
        tbl = FMT_TABLE + (k + FMT_12_GRIDS) * FMT_M_SIZE;
        if (Mvalue > 1) exptbl = FMT_EXP_TABLE + k + FMT_12_GRIDS;
    } else {
        fmt = FMT_STORE;
        t_inv = DHALF / Tvalue;
        del = Const * DHALF * sqrt(TWO_PI * t_inv);
        *fmt++ = del;
        k = IONE;
        for (i = 1; i <= Mvalue; i++) {
            del *= (t_inv * k);
            *fmt++ = del;
            k += ITWO;
        }
        return 0;
    }

    if (Mvalue < 1) {
        *FMT_STORE = Const * (*tbl + del * (*(tbl+1) + del * DHALF * (*(tbl+2) + del * DTHIRD * (*(tbl+3)))));
    } else if (Mvalue < 2) {
        FMT_STORE[0] = Const * (*tbl + del * (*(tbl+1) + del * DHALF * (*(tbl+2) + del * DTHIRD * (*(tbl+3)))));
        FMT_STORE[1] = Const * (*(tbl+1) + del * (*(tbl+2) + del * DHALF * (*(tbl+3) + del * DTHIRD * (*(tbl+4)))));
    } else {
        tbl += Mvalue;
        fmt = FMT_STORE + Mvalue;
        *fmt = Const * (*tbl + del * (*(tbl+1) + del * DHALF * (*(tbl+2) + del * DTHIRD * (*(tbl+3)))));
        evl = Const * (*exptbl) * (FMT_DONE + del * (FMT_DONE + del * DHALF * (FMT_DONE + del * DTHIRD)));
        tv2 = Tvalue + Tvalue;
        
        switch (Mvalue) {
            case 16 : { tm = *fmt--; *fmt = DINV31 * (evl + tv2 * tm); }
            case 15 : { tm = *fmt--; *fmt = DINV29 * (evl + tv2 * tm); }
            case 14 : { tm = *fmt--; *fmt = DINV27 * (evl + tv2 * tm); }
            case 13 : { tm = *fmt--; *fmt = DINV25 * (evl + tv2 * tm); }
            case 12 : { tm = *fmt--; *fmt = DINV23 * (evl + tv2 * tm); }
            case 11 : { tm = *fmt--; *fmt = DINV21 * (evl + tv2 * tm); }
            case 10 : { tm = *fmt--; *fmt = DINV19 * (evl + tv2 * tm); }
            case  9 : { tm = *fmt--; *fmt = DINV17 * (evl + tv2 * tm); }
            case  8 : { tm = *fmt--; *fmt = DINV15 * (evl + tv2 * tm); }
            case  7 : { tm = *fmt--; *fmt = DINV13 * (evl + tv2 * tm); }
            case  6 : { tm = *fmt--; *fmt = DINV11 * (evl + tv2 * tm); }
            case  5 : { tm = *fmt--; *fmt = DINV09 * (evl + tv2 * tm); }
            case  4 : { tm = *fmt--; *fmt = DINV07 * (evl + tv2 * tm); }
            case  3 : { tm = *fmt--; *fmt = DINV05 * (evl + tv2 * tm); }
            case  2 : { tm = *fmt--; *fmt = DINV03 * (evl + tv2 * tm); }
            case  1 : { tm = *fmt--; *fmt =          (evl + tv2 * tm); }
        }
    }
    return 0;
}

//==============================================================================
//==============================================================================
//==============================================================================

void fmt_initialize_(const int *max_lqn) {
    static int initialized = 0;
    // max_lqn は無視されます（自動的に17までのテーブルが確保されます）
    if (!initialized) {
        FMT_make_table(stderr); // エラーメッセージは標準エラー出力へ
        atexit((void (*)(void))FMT_free_table); // プログラム終了時にメモリ解放
        initialized = 1;
    }
}

void fmt_(double f[], const int *pm, const double *pt, const double *pcoef) {
    FMT_set_store_area(f);
    FMT_value_fast((MOE_LI)*pm, (MOE_LD)*pt, (MOE_LD)*pcoef);
}