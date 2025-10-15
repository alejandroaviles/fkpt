#ifndef CLPT_TYPES_H
#define CLPT_TYPES_H

#include "globaldefs.h"

/* ---------- r-space tables (filled in compute_clpt) ---------- */
typedef struct r_arrays {
    real *rTab;
    real *xi_00T, *xi_10T, *xi_20T, *xi_01T, *xi_02T, *xi_11T, *xi_eftT;
    real *xi_001T, *xi_002T, *xi_101T, *xi_011T;
    real *xi_LT, *xi_zaT, *xi_AT, *xi_WT;
    real *v_00T, *v_10T, *v_20T, *v_01T, *v_11T, *v_001T, *v_101T, *v_eftT;
    real *spar_00T, *spar_10T, *spar_20T, *spar_01T, *spar_001T;
    real *sperp_00T, *sperp_10T, *sperp_20T, *sperp_01T, *sperp_001T;
} r_arrays;

/* ---------- q-space tables (used by X/Y/U/iBessel/nabla2xi helpers) ---------- */
typedef struct q_arrays {
    real *qTab;
    real *XLT, *YLT, *XloopT, *YloopT, *preVT, *TT, *X10T, *Y10T;
    real *ULT, *UloopT, *U11T, *U20T;
    real *dotVaT, *ddotVaT, *dotVbT;
    real *V10T;
    real *iBessel0T, *iBessel2T, *iBessel4T;
    real *nabla2xiT;
} q_arrays;

/* Derivative/aux arrays in q-space (field-for-field mirror of subset above) */
typedef struct q_arraysd {
    real *XLT, *YLT, *XloopT, *YloopT, *preVT, *TT, *X10T, *Y10T;
    real *ULT, *UloopT, *U11T, *U20T;
    real *dotVaT, *ddotVaT, *dotVbT;
    real *V10T;
    real *iBessel0T, *iBessel2T, *iBessel4T;
    real *nabla2xiT;
} q_arraysd;

/* ---------- ZA correlation container ---------- */
typedef struct {
    real r;   /* separation */
    real xi;  /* xi_ZA */
} global_zacorrfunctions;
typedef global_zacorrfunctions* global_zacorrfunctions_ptr;

#ifndef rzacorrfun
#define rzacorrfun(p)        ((p)->r)
#endif
#ifndef xizacorrfun
#define xizacorrfun(p)       ((p)->xi)
#endif

/* ---------- CLPT correlation container ---------- */
typedef struct {
    /* independent variable */
    real r;

    /* xi-type terms */
    real xiA, xiW, xi10, xi20, xi01, xi02, xi11, nabla2xi;
    real xi001, xi002, xi101, xi011;

    /* pairwise velocity moments v12_* */
    real v12_00, v12_10, v12_20, v12_01, v12_11, v12_001, v12_101, v12_eft;

    /* parallel dispersion s12par_* */
    real s12par_00, s12par_10, s12par_20, s12par_01, s12par_001;

    /* perpendicular dispersion s12perp_* */
    real s12perp_00, s12perp_10, s12perp_20, s12perp_01, s12perp_001;
} global_clptcorrfunctions;
typedef global_clptcorrfunctions* global_clptcorrfunctions_ptr;

/* Convenient access macros used in clpt.c */
#ifndef rclptcorrfun
#define rclptcorrfun(p)         ((p)->r)
#endif

#ifndef xiAclptcorrfun
#define xiAclptcorrfun(p)       ((p)->xiA)
#endif
#ifndef xiWclptcorrfun
#define xiWclptcorrfun(p)       ((p)->xiW)
#endif
#ifndef xi10clptcorrfun
#define xi10clptcorrfun(p)      ((p)->xi10)
#endif
#ifndef xi20clptcorrfun
#define xi20clptcorrfun(p)      ((p)->xi20)
#endif
#ifndef xi01clptcorrfun
#define xi01clptcorrfun(p)      ((p)->xi01)
#endif
#ifndef xi02clptcorrfun
#define xi02clptcorrfun(p)      ((p)->xi02)
#endif
#ifndef xi11clptcorrfun
#define xi11clptcorrfun(p)      ((p)->xi11)
#endif
#ifndef nabla2xiclptcorrfun
#define nabla2xiclptcorrfun(p)  ((p)->nabla2xi)
#endif
#ifndef xi001clptcorrfun
#define xi001clptcorrfun(p)     ((p)->xi001)
#endif
#ifndef xi002clptcorrfun
#define xi002clptcorrfun(p)     ((p)->xi002)
#endif
#ifndef xi101clptcorrfun
#define xi101clptcorrfun(p)     ((p)->xi101)
#endif
#ifndef xi011clptcorrfun
#define xi011clptcorrfun(p)     ((p)->xi011)
#endif

#ifndef v12_00_vfun
#define v12_00_vfun(p)          ((p)->v12_00)
#endif
#ifndef v12_10_vfun
#define v12_10_vfun(p)          ((p)->v12_10)
#endif
#ifndef v12_20_vfun
#define v12_20_vfun(p)          ((p)->v12_20)
#endif
#ifndef v12_01_vfun
#define v12_01_vfun(p)          ((p)->v12_01)
#endif
#ifndef v12_11_vfun
#define v12_11_vfun(p)          ((p)->v12_11)
#endif
#ifndef v12_001_vfun
#define v12_001_vfun(p)         ((p)->v12_001)
#endif
#ifndef v12_101_vfun
#define v12_101_vfun(p)         ((p)->v12_101)
#endif
#ifndef v12_eft_vfun
#define v12_eft_vfun(p)         ((p)->v12_eft)
#endif

#ifndef s12par_00_sfun
#define s12par_00_sfun(p)       ((p)->s12par_00)
#endif
#ifndef s12par_10_sfun
#define s12par_10_sfun(p)       ((p)->s12par_10)
#endif
#ifndef s12par_20_sfun
#define s12par_20_sfun(p)       ((p)->s12par_20)
#endif
#ifndef s12par_01_sfun
#define s12par_01_sfun(p)       ((p)->s12par_01)
#endif
#ifndef s12par_001_sfun
#define s12par_001_sfun(p)      ((p)->s12par_001)
#endif

#ifndef s12perp_00_sfun
#define s12perp_00_sfun(p)      ((p)->s12perp_00)
#endif
#ifndef s12perp_10_sfun
#define s12perp_10_sfun(p)      ((p)->s12perp_10)
#endif
#ifndef s12perp_20_sfun
#define s12perp_20_sfun(p)      ((p)->s12perp_20)
#endif
#ifndef s12perp_01_sfun
#define s12perp_01_sfun(p)      ((p)->s12perp_01)
#endif
#ifndef s12perp_001_sfun
#define s12perp_001_sfun(p)     ((p)->s12perp_001)
#endif

#endif /* CLPT_TYPES_H */