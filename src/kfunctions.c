/*==============================================================================
 NAME: kfunctions.c				[code for fk - Perturbation Theory]
 Alejandro Aviles (avilescervantes@gmail.com)
 ================================================================================ 
*/

#include "globaldefs.h"
#include "protodefs.h"
#include "models.h"
#include <stdlib.h>   /* calloc, free, abort */
#include <stdio.h>    /* fprintf */
#include <math.h>     /* fabs */

/* keep track of current capacity so we can reuse between runs */
static int kF_cap = 0;

#define ALLOC1(arr) do{ \
  (arr) = (real*)calloc(cmd.Nk + 1, sizeof(real)); \
  if(!(arr)){ fprintf(stderr,"[alloc_kFarrays] OOM for %s (Nk=%d)\n", #arr, cmd.Nk); abort(); } \
}while(0)

#define FREE1(arr) do{ if(arr){ free(arr); (arr)=NULL; } }while(0)

global void fkpt_free_kFarrays(void){
  /* wiggle */
  FREE1(kFArrays.kT);
  FREE1(kFArrays.pklT);
  FREE1(kFArrays.P22ddT);
  FREE1(kFArrays.P22duT);
  FREE1(kFArrays.P22uuT);
  FREE1(kFArrays.I1udd1AT);
  FREE1(kFArrays.I2uud1AT);
  FREE1(kFArrays.I2uud2AT);
  FREE1(kFArrays.I3uuu2AT);
  FREE1(kFArrays.I3uuu3AT);
  FREE1(kFArrays.I2uudd1BpCT);
  FREE1(kFArrays.I2uudd2BpCT);
  FREE1(kFArrays.I3uuud2BpCT);
  FREE1(kFArrays.I3uuud3BpCT);
  FREE1(kFArrays.I4uuuu2BpCT);
  FREE1(kFArrays.I4uuuu3BpCT);
  FREE1(kFArrays.I4uuuu4BpCT);
  FREE1(kFArrays.Pb1b2T);
  FREE1(kFArrays.Pb1bs2T);
  FREE1(kFArrays.Pb22T);
  FREE1(kFArrays.Pb2s2T);
  FREE1(kFArrays.Ps22T);
  FREE1(kFArrays.Pb2thetaT);
  FREE1(kFArrays.Pbs2thetaT);
  FREE1(kFArrays.P13ddT);
  FREE1(kFArrays.P13duT);
  FREE1(kFArrays.P13uuT);
  FREE1(kFArrays.sigma32PSLT);
  FREE1(kFArrays.fkT);

  /* no-wiggle */
  FREE1(kFArrays_nw.kT);
  FREE1(kFArrays_nw.pklT);
  FREE1(kFArrays_nw.P22ddT);
  FREE1(kFArrays_nw.P22duT);
  FREE1(kFArrays_nw.P22uuT);
  FREE1(kFArrays_nw.I1udd1AT);
  FREE1(kFArrays_nw.I2uud1AT);
  FREE1(kFArrays_nw.I2uud2AT);
  FREE1(kFArrays_nw.I3uuu2AT);
  FREE1(kFArrays_nw.I3uuu3AT);
  FREE1(kFArrays_nw.I2uudd1BpCT);
  FREE1(kFArrays_nw.I2uudd2BpCT);
  FREE1(kFArrays_nw.I3uuud2BpCT);
  FREE1(kFArrays_nw.I3uuud3BpCT);
  FREE1(kFArrays_nw.I4uuuu2BpCT);
  FREE1(kFArrays_nw.I4uuuu3BpCT);
  FREE1(kFArrays_nw.I4uuuu4BpCT);
  FREE1(kFArrays_nw.Pb1b2T);
  FREE1(kFArrays_nw.Pb1bs2T);
  FREE1(kFArrays_nw.Pb22T);
  FREE1(kFArrays_nw.Pb2s2T);
  FREE1(kFArrays_nw.Ps22T);
  FREE1(kFArrays_nw.Pb2thetaT);
  FREE1(kFArrays_nw.Pbs2thetaT);
  FREE1(kFArrays_nw.P13ddT);
  FREE1(kFArrays_nw.P13duT);
  FREE1(kFArrays_nw.P13uuT);
  FREE1(kFArrays_nw.sigma32PSLT);
  FREE1(kFArrays_nw.fkT);

  kF_cap = 0;
}

static void alloc_kFarrays(void){
  if (kF_cap == cmd.Nk && kFArrays.kT && kFArrays_nw.kT) return; /* already ok */
  if (kF_cap) fkpt_free_kFarrays();

  /* wiggle */
  ALLOC1(kFArrays.kT);
  ALLOC1(kFArrays.pklT);
  ALLOC1(kFArrays.P22ddT);
  ALLOC1(kFArrays.P22duT);
  ALLOC1(kFArrays.P22uuT);
  ALLOC1(kFArrays.I1udd1AT);
  ALLOC1(kFArrays.I2uud1AT);
  ALLOC1(kFArrays.I2uud2AT);
  ALLOC1(kFArrays.I3uuu2AT);
  ALLOC1(kFArrays.I3uuu3AT);
  ALLOC1(kFArrays.I2uudd1BpCT);
  ALLOC1(kFArrays.I2uudd2BpCT);
  ALLOC1(kFArrays.I3uuud2BpCT);
  ALLOC1(kFArrays.I3uuud3BpCT);
  ALLOC1(kFArrays.I4uuuu2BpCT);
  ALLOC1(kFArrays.I4uuuu3BpCT);
  ALLOC1(kFArrays.I4uuuu4BpCT);
  ALLOC1(kFArrays.Pb1b2T);
  ALLOC1(kFArrays.Pb1bs2T);
  ALLOC1(kFArrays.Pb22T);
  ALLOC1(kFArrays.Pb2s2T);
  ALLOC1(kFArrays.Ps22T);
  ALLOC1(kFArrays.Pb2thetaT);
  ALLOC1(kFArrays.Pbs2thetaT);
  ALLOC1(kFArrays.P13ddT);
  ALLOC1(kFArrays.P13duT);
  ALLOC1(kFArrays.P13uuT);
  ALLOC1(kFArrays.sigma32PSLT);
  ALLOC1(kFArrays.fkT);

  /* no-wiggle */
  ALLOC1(kFArrays_nw.kT);
  ALLOC1(kFArrays_nw.pklT);
  ALLOC1(kFArrays_nw.P22ddT);
  ALLOC1(kFArrays_nw.P22duT);
  ALLOC1(kFArrays_nw.P22uuT);
  ALLOC1(kFArrays_nw.I1udd1AT);
  ALLOC1(kFArrays_nw.I2uud1AT);
  ALLOC1(kFArrays_nw.I2uud2AT);
  ALLOC1(kFArrays_nw.I3uuu2AT);
  ALLOC1(kFArrays_nw.I3uuu3AT);
  ALLOC1(kFArrays_nw.I2uudd1BpCT);
  ALLOC1(kFArrays_nw.I2uudd2BpCT);
  ALLOC1(kFArrays_nw.I3uuud2BpCT);
  ALLOC1(kFArrays_nw.I3uuud3BpCT);
  ALLOC1(kFArrays_nw.I4uuuu2BpCT);
  ALLOC1(kFArrays_nw.I4uuuu3BpCT);
  ALLOC1(kFArrays_nw.I4uuuu4BpCT);
  ALLOC1(kFArrays_nw.Pb1b2T);
  ALLOC1(kFArrays_nw.Pb1bs2T);
  ALLOC1(kFArrays_nw.Pb22T);
  ALLOC1(kFArrays_nw.Pb2s2T);
  ALLOC1(kFArrays_nw.Ps22T);
  ALLOC1(kFArrays_nw.Pb2thetaT);
  ALLOC1(kFArrays_nw.Pbs2thetaT);
  ALLOC1(kFArrays_nw.P13ddT);
  ALLOC1(kFArrays_nw.P13duT);
  ALLOC1(kFArrays_nw.P13uuT);
  ALLOC1(kFArrays_nw.sigma32PSLT);
  ALLOC1(kFArrays_nw.fkT);

  kF_cap = cmd.Nk;

  /* debug
  fprintf(stderr,"[alloc_kFarrays] cap=%d k=%p k_nw=%p pklin=%p pklin_nw=%p\n",
          cmd.Nk,(void*)kFArrays.kT,(void*)kFArrays_nw.kT,
          (void*)kFArrays.pklT,(void*)kFArrays_nw.pklT);
  */
}

#undef ALLOC1
#undef FREE1

local void quadrature(real ki);
local void k_functions(void);
//~ local void pk_non_wiggle(void);

local real Interpolation_nr(real k, double kPS[], double pPS[], int nPS, double pPS2[]);

local real get_sigma8(void);
local real sigma28_function_int(real y);

local real sigma2L_function_int(real y);
local real Sigma2_int(real y);
local real deltaSigma2_int(real y);
local real sigma2v_function_int(real y);
local real sigma_constants(void);

#define KMIN    1.0e-20

/* ===== debug helpers ===== */
#ifndef KF_DEBUG
#define KF_DEBUG 0   /* set to 0 to silence without touching calls */
#endif

#if KF_DEBUG
  #define KDBG(fmt, ...) do { \
      fprintf(stderr, "[fkpt:set_external_fk] " fmt "\n", ##__VA_ARGS__); \
      fflush(stderr); \
  } while (0)
#else
  #define KDBG(...) do{}while(0)
#endif

// CGQ MOD: get growth from ISiTGR
global void set_external_fk(int n, const double *k_in, const double *fk_in, double f0)
{
    KDBG("ENTER n=%d  k_in=%p  fk_in=%p  f0=%g", n, (void*)k_in, (void*)fk_in, f0);
    KDBG("state before: nPSLT=%d  kPS=%p  fkT=%p  fkT2=%p", nPSLT, (void*)kPS, (void*)fkT, (void*)fkT2);

    if (nPSLT <= 0 || kPS == NULL) {
        KDBG("ABORT: internal kPS not ready (nPSLT=%d, kPS=%p)", nPSLT, (void*)kPS);
        fprintf(stderr, "[set_external_fk] internal kPS not ready (nPSLT=%d)\n", nPSLT);
        abort();
    }
    if (!fkT || !fkT2) {
        KDBG("ABORT: fkT/fkT2 not allocated yet (fkT=%p fkT2=%p) — call PSLTable first", (void*)fkT, (void*)fkT2);
        fprintf(stderr, "[set_external_fk] fkT/fkT2 not allocated yet (call PSLTable first)\n");
        abort();
    }
    if (n < 2 || !k_in || !fk_in) {
        KDBG("ABORT: bad inputs (n=%d, k_in=%p, fk_in=%p)", n, (void*)k_in, (void*)fk_in);
        fprintf(stderr, "[set_external_fk] bad inputs (n=%d, k_in=%p, fk_in=%p)\n", n, (void*)k_in, (void*)fk_in);
        abort();
    }

    /* quick peek at ends to catch garbage pointers */
    KDBG("k_in[0]=%.6e fk_in[0]=%.6e   k_in[n-1]=%.6e fk_in[n-1]=%.6e",
         k_in[0], fk_in[0], k_in[n-1], fk_in[n-1]);

    /* verify monotonic k_in */
    for (int i = 1; i < n; ++i) {
        if (!(k_in[i] > k_in[i-1])) {
            KDBG("ABORT: k_in not strictly increasing at i=%d: %.16e <= %.16e",
                 i, k_in[i], k_in[i-1]);
            fprintf(stderr, "[set_external_fk] k_in must be strictly increasing at i=%d (%.16e <= %.16e)\n",
                    i, k_in[i], k_in[i-1]);
            abort();
        }
    }
    KDBG("k_in monotonic OK");

    /* try direct grid match */
    int direct_ok = (n == nPSLT);
    if (direct_ok) {
        const double tol0 = 1e-10;
        for (int i = 1; i <= nPSLT; ++i) {
            const double kin  = k_in[i-1];
            const double kref = kPS[i];
            const double tol  = tol0 * (kref > 0 ? kref : 1.0);
            if (fabs(kin - kref) > tol) { direct_ok = 0; break; }
        }
    }
    KDBG("grid match? %s", direct_ok ? "YES (direct copy)" : "NO (will resample)");

    if (direct_ok) {
        for (int i = 1; i <= nPSLT; ++i) fkT[i] = fk_in[i-1];
        KDBG("direct copy done; fkT[1]=%.6e fkT[%d]=%.6e", fkT[1], nPSLT, fkT[nPSLT]);
    } else {
        /* resample via temporary spline */
        KDBG("alloc temp x/y/y2 for spline (n=%d)", n);
        double *x  = dvector(1, n);
        double *y  = dvector(1, n);
        double *y2 = dvector(1, n);
        if (!x || !y || !y2) { KDBG("ABORT: temp alloc failed"); abort(); }

        for (int i = 1; i <= n; ++i) { x[i] = k_in[i-1]; y[i] = fk_in[i-1]; }
        KDBG("build spline on input grid");
        spline(x, y, n, 1.0e30, 1.0e30, y2);

        KDBG("evaluate spline on kPS (nPSLT=%d)", nPSLT);
        for (int i = 1; i <= nPSLT; ++i) {
            double fki, kc = kPS[i];
            if (kc < x[1]) kc = x[1];
            if (kc > x[n]) kc = x[n];
            splint(x, y, y2, n, kc, &fki);
            fkT[i] = fki;
        }
        KDBG("resample done; fkT[1]=%.6e fkT[%d]=%.6e", fkT[1], nPSLT, fkT[nPSLT]);

        KDBG("free temp x/y/y2");
        free_dvector(y2, 1, n);
        free_dvector(y,  1, n);
        free_dvector(x,  1, n);
    }

    KDBG("rebuild fk spline on internal grid");
    spline(kPS, fkT, nPSLT, 1.0e30, 1.0e30, fkT2);

    /* set f0 and flag */
    gd.f0 = (isfinite(f0) && f0 != 0.0) ? f0 : fkT[1];
    KDBG("set gd.f0=%.6e  (fkT[1]=%.6e)", gd.f0, fkT[1]);
    if (gd.f0 == 0.0) {
        KDBG("ABORT: f0=0");
        fprintf(stderr, "[set_external_fk] fatal: f0=0 after fallback\n");
        abort();
    }

    gd.use_external_fk = 1;
    KDBG("set gd.use_external_fk=1; EXIT OK");
}
// CGQ MOD: get growth from ISiTGR

global void compute_kfunctions(void)
{
    stream outstrQsRs, outtables;
    real dk;
    real bTime;
    real kBAOmin=0.005, kBAOmax=1.0, epsquadsave;
    int iBAOmin, iBAOmax;
    global_D2v2_ptr ptmp;
    global_D3v2_ptr ptmpR1;
    real fR0save;
    //~ real sigma2psi;

    bTime = second();
 
    // before: gd.f0 = f_growth_LCDM();
    if (!gd.use_external_fk) {
        gd.f0 = f_growth_LCDM();
        /* also build fkT internally if that's your current default */
    }
    /* else: fkT and gd.f0 were set by your set_external_fk() wrapper */
    // else: gd.f0 already set by set_external_fk()
    //gd.f0=f_growth_LCDM();  //Modificacion-LCDMfk
    
    ptmp = DsSecondOrder_func(KMIN, KMIN, KMIN);
    KA_LCDM = DA2D2(ptmp) / ( (3.0/7.0) * Dpk1D2(ptmp) * Dpk2D2(ptmp) );
    KAp_LCDM = DA2primeD2(ptmp) 
                / ( (3.0/7.0) * Dpk1D2(ptmp) * Dpk2D2(ptmp) ) -
               2.0 * DA2D2(ptmp) 
                / ( (3.0/7.0) * Dpk1D2(ptmp) * Dpk2D2(ptmp) )* gd.f0; //Modificacion-LCDMfk
    KB_LCDM = KA_LCDM;

    ptmpR1 = DsThirdOrder_func(0.0000001, KMIN, KMIN);   
    KR1_LCDM = (21.0/5.0)*D3symmD3(ptmpR1)
        /( DpkD3(ptmpR1)*DppD3(ptmpR1)*DppD3(ptmpR1) );     
    KR1p_LCDM = (21.0/5.0)*D3symmprimeD3(ptmpR1)
        /( DpkD3(ptmpR1)*DppD3(ptmpR1)*DppD3(ptmpR1) )/(3.*gd.f0);   //Modificacion-LCDMfk     

	get_sigma8();
	if(cmd.chatty==3) fprintf(stdout,"%g\n",gd.sigma8);
    
    /* debug
    fprintf(stderr,"[compute_kfunctions] about to enter sigma_constants()\n"); fflush(stderr);
    fprintf(stderr,"[compute_kfunctions] ptrs: kPS=%p pPS=%p pPS2=%p  fkT=%p fkT2=%p  pPS_nw=%p pPS2_nw=%p  nPSLT=%d\n",
            (void*)kPS,(void*)pPS,(void*)pPS2,(void*)fkT,(void*)fkT2,(void*)pPS_nw,(void*)pPS2_nw,nPSLT);
    fflush(stderr);
    */

    sigma_constants();
    
    /* debug
    fprintf(stderr,"[compute_kfunctions] sigma_constants() done: sigma2L=%g sigma2v=%g Sigma2=%g deltaSigma2=%g\n",
            gd.sigma2L, gd.sigma2v, gd.Sigma2, gd.deltaSigma2);
    fflush(stderr);
    */
    
    if (cmd.chatty==1){
        fprintf(stdout,"\nA_LCDM=%g, Ap_LCDM=%g, KR1_LCDM = %g, KR1p_LCDM = %g",
                KA_LCDM, KAp_LCDM, KR1_LCDM, KR1p_LCDM);
        fprintf(stdout,"\nsigma quadratures from kmin = %g to kmax = %g", kPS[1],kPS[nPSLT]);
        fprintf(stdout,"\ns2psi = %g,   s2v = %g,   Sigma2 = %g,   deltaSigma2 = %g",
                gd.sigma2L, gd.sigma2v, gd.Sigma2, gd.deltaSigma2);
        fprintf(stdout,"\nk-functions:");
        fprintf(stdout," Nk=%d values from kmin=%g to kmax=%g ",
                cmd.Nk, cmd.kmin, cmd.kmax);
    }
    
    /* debug
    fprintf(stderr,"[compute_kfunctions] about to enter k_functions(): Nk=%d kmin=%g kmax=%g nquadSteps=%d\n",
            cmd.Nk, cmd.kmin, cmd.kmax, cmd.nquadSteps);
    fflush(stderr);
    */
    
    alloc_kFarrays();
    k_functions();
    
    /* debug
    fprintf(stderr,"[compute_kfunctions] returned from k_functions()\n"); fflush(stderr);
    */
    
    
    if(cmd.chatty==1) fprintf(stdout,"...time = %g seconds",second()-bTime);    

}
#undef KMIN



#define _K_LOGSPACED_  1
local void k_functions(void)
{
    global_kFs qrs, qrs_nw;
    real kval, ki, pkl, pkl_nw, fk;
    int i;
    real dk;

    /* Handy printf with flush (disabled) */
    #define PRF(...) /* debug */

    /* debug
    PRF("[k_functions] ENTER: Nk=%d kmin=%g kmax=%g nPSLT=%d nquadSteps=%d\n",
        cmd.Nk, cmd.kmin, cmd.kmax, nPSLT, cmd.nquadSteps);

    PRF("[k_functions] OUT ptrs (wiggle): ...\n", ...);
    PRF("[k_functions] OUT ptrs (no-wiggle): ...\n", ...);
    */

    if (_K_LOGSPACED_ == 1){
        dk = (rlog10(cmd.kmax/cmd.kmin))/((real)(cmd.Nk - 1));
    } else {
        dk = (cmd.kmax-cmd.kmin)/((real)(cmd.Nk - 1));
    }
    /* debug PRF("[k_functions] dk=%g (logspaced=%d)\n", dk, (int)_K_LOGSPACED_); */

    /* --- verify output arrays are allocated --- */
    #define REQ(ptr,name) do{ \
      if(!(ptr)){ fprintf(stderr,"[k_functions] NULL %s\n", name); fflush(stderr); abort(); } \
    }while(0)

    /* wiggle */
    REQ(kFArrays.kT, "kFArrays.kT");
    REQ(kFArrays.pklT, "kFArrays.pklT");
    REQ(kFArrays.P22ddT, "kFArrays.P22ddT");
    REQ(kFArrays.P22duT, "kFArrays.P22duT");
    REQ(kFArrays.P22uuT, "kFArrays.P22uuT");
    REQ(kFArrays.I1udd1AT, "kFArrays.I1udd1AT");
    REQ(kFArrays.I2uud1AT, "kFArrays.I2uud1AT");
    REQ(kFArrays.I2uud2AT, "kFArrays.I2uud2AT");
    REQ(kFArrays.I3uuu2AT, "kFArrays.I3uuu2AT");
    REQ(kFArrays.I3uuu3AT, "kFArrays.I3uuu3AT");
    REQ(kFArrays.I2uudd1BpCT, "kFArrays.I2uudd1BpCT");
    REQ(kFArrays.I2uudd2BpCT, "kFArrays.I2uudd2BpCT");
    REQ(kFArrays.I3uuud2BpCT, "kFArrays.I3uuud2BpCT");
    REQ(kFArrays.I3uuud3BpCT, "kFArrays.I3uuud3BpCT");
    REQ(kFArrays.I4uuuu2BpCT, "kFArrays.I4uuuu2BpCT");
    REQ(kFArrays.I4uuuu3BpCT, "kFArrays.I4uuuu3BpCT");
    REQ(kFArrays.I4uuuu4BpCT, "kFArrays.I4uuuu4BpCT");
    REQ(kFArrays.Pb1b2T, "kFArrays.Pb1b2T");
    REQ(kFArrays.Pb1bs2T, "kFArrays.Pb1bs2T");
    REQ(kFArrays.Pb22T, "kFArrays.Pb22T");
    REQ(kFArrays.Pb2s2T, "kFArrays.Pb2s2T");
    REQ(kFArrays.Ps22T, "kFArrays.Ps22T");
    REQ(kFArrays.Pb2thetaT, "kFArrays.Pb2thetaT");
    REQ(kFArrays.Pbs2thetaT, "kFArrays.Pbs2thetaT");
    REQ(kFArrays.P13ddT, "kFArrays.P13ddT");
    REQ(kFArrays.P13duT, "kFArrays.P13duT");
    REQ(kFArrays.P13uuT, "kFArrays.P13uuT");
    REQ(kFArrays.sigma32PSLT, "kFArrays.sigma32PSLT");
    REQ(kFArrays.fkT, "kFArrays.fkT");

    /* no-wiggle */
    REQ(kFArrays_nw.kT, "kFArrays_nw.kT");
    REQ(kFArrays_nw.pklT, "kFArrays_nw.pklT");
    REQ(kFArrays_nw.P22ddT, "kFArrays_nw.P22ddT");
    REQ(kFArrays_nw.P22duT, "kFArrays_nw.P22duT");
    REQ(kFArrays_nw.P22uuT, "kFArrays_nw.P22uuT");
    REQ(kFArrays_nw.I1udd1AT, "kFArrays_nw.I1udd1AT");
    REQ(kFArrays_nw.I2uud1AT, "kFArrays_nw.I2uud1AT");
    REQ(kFArrays_nw.I2uud2AT, "kFArrays_nw.I2uud2AT");
    REQ(kFArrays_nw.I3uuu2AT, "kFArrays_nw.I3uuu2AT");
    REQ(kFArrays_nw.I3uuu3AT, "kFArrays_nw.I3uuu3AT");
    REQ(kFArrays_nw.I2uudd1BpCT, "kFArrays_nw.I2uudd1BpCT");
    REQ(kFArrays_nw.I2uudd2BpCT, "kFArrays_nw.I2uudd2BpCT");
    REQ(kFArrays_nw.I3uuud2BpCT, "kFArrays_nw.I3uuud2BpCT");
    REQ(kFArrays_nw.I3uuud3BpCT, "kFArrays_nw.I3uuud3BpCT");
    REQ(kFArrays_nw.I4uuuu2BpCT, "kFArrays_nw.I4uuuu2BpCT");
    REQ(kFArrays_nw.I4uuuu3BpCT, "kFArrays_nw.I4uuuu3BpCT");
    REQ(kFArrays_nw.I4uuuu4BpCT, "kFArrays_nw.I4uuuu4BpCT");
    REQ(kFArrays_nw.Pb1b2T, "kFArrays_nw.Pb1b2T");
    REQ(kFArrays_nw.Pb1bs2T, "kFArrays_nw.Pb1bs2T");
    REQ(kFArrays_nw.Pb22T, "kFArrays_nw.Pb22T");
    REQ(kFArrays_nw.Pb2s2T, "kFArrays_nw.Pb2s2T");
    REQ(kFArrays_nw.Ps22T, "kFArrays_nw.Ps22T");
    REQ(kFArrays_nw.Pb2thetaT, "kFArrays_nw.Pb2thetaT");
    REQ(kFArrays_nw.Pbs2thetaT, "kFArrays_nw.Pbs2thetaT");
    REQ(kFArrays_nw.P13ddT, "kFArrays_nw.P13ddT");
    REQ(kFArrays_nw.P13duT, "kFArrays_nw.P13duT");
    REQ(kFArrays_nw.P13uuT, "kFArrays_nw.P13uuT");
    REQ(kFArrays_nw.sigma32PSLT, "kFArrays_nw.sigma32PSLT");
    REQ(kFArrays_nw.fkT, "kFArrays_nw.fkT");

    #undef REQ

    for (i = 1; i <= cmd.Nk; i++) {
        if (_K_LOGSPACED_ == 1){
            kval = rlog10(cmd.kmin) + dk*((real)(i - 1));
            ki   = rpow(10.0,kval);
        } else {
            ki   = cmd.kmin + dk*((real)(i - 1));
        }

        /* debug PRF("[k_functions] i=%d  ki=%e  (pre qrs)\n", i, ki); */

        /* compute both wiggle and no-wiggle blocks */
        qrs    = ki_functions_driver(ki, kPS, pPS,    nPSLT, pPS2);
        qrs_nw = ki_functions_driver(ki, kPS, pPS_nw, nPSLT, pPS2_nw);

        /* debug
        PRF("[k_functions] i=%d  (post qrs)  P22dd=%g  P13dd=%g\n", i, qrs.P22dd, qrs.P13dd);
        */

        /* scalar helpers */
        pkl    = psInterpolation_nr(ki, kPS, pPS,    nPSLT);
        fk     = Interpolation_nr(ki, kPS, fkT,   nPSLT, fkT2);
        pkl_nw = Interpolation_nr(ki, kPS, pPS_nw, nPSLT, pPS2_nw);

        /* debug
        PRF("[k_functions] i=%d  RHS base:  pkl=%g  fk=%g  pkl_nw=%g\n", i, pkl, fk, pkl_nw);
        */

        /* ---- writes: use [i], not [i-1] ---- */
		kFArrays.kT[i-1]           =  ki ;
		kFArrays.P22ddT[i-1]       =  qrs.P22dd;
		kFArrays.P22duT[i-1]       =  qrs.P22du;
		kFArrays.P22uuT[i-1]       =  qrs.P22uu;
		// A 
		kFArrays.I1udd1AT[i-1]     =  qrs.I1udd1A;
		kFArrays.I2uud1AT[i-1]     =  qrs.I2uud1A;
		kFArrays.I2uud2AT[i-1]     =  qrs.I2uud2A;
		kFArrays.I3uuu2AT[i-1]     =  qrs.I3uuu2A;
		kFArrays.I3uuu3AT[i-1]     =  qrs.I3uuu3A;
		//  B plus C   
		kFArrays.I2uudd1BpCT[i-1]  =  qrs.I2uudd1BpC;
		kFArrays.I2uudd2BpCT[i-1]  =  qrs.I2uudd2BpC;
		kFArrays.I3uuud2BpCT[i-1]  =  qrs.I3uuud2BpC;
		kFArrays.I3uuud3BpCT[i-1]  =  qrs.I3uuud3BpC;
		kFArrays.I4uuuu2BpCT[i-1]  =  qrs.I4uuuu2BpC;
		kFArrays.I4uuuu3BpCT[i-1]  =  qrs.I4uuuu3BpC;
		kFArrays.I4uuuu4BpCT[i-1]  =  qrs.I4uuuu4BpC;
		//  Bias
		kFArrays.Pb1b2T[i-1]       =  qrs.Pb1b2;
		kFArrays.Pb1bs2T[i-1]      =  qrs.Pb1bs2;
		kFArrays.Pb22T[i-1]        =  qrs.Pb22;
		kFArrays.Pb2s2T[i-1]       =  qrs.Pb2s2;
		kFArrays.Ps22T[i-1]        =  qrs.Ps22;
		kFArrays.Pb2thetaT[i-1]    =  qrs.Pb2theta;
		kFArrays.Pbs2thetaT[i-1]   =  qrs.Pbs2theta;
		//
		kFArrays.P13ddT[i-1]       =  qrs.P13dd;
		kFArrays.P13duT[i-1]       =  qrs.P13du;
		kFArrays.P13uuT[i-1]       =  qrs.P13uu;
		kFArrays.sigma32PSLT[i-1]  =  qrs.sigma32PSL;
		kFArrays.pklT[i-1]         =  pkl;
		kFArrays.fkT[i-1]          =  fk;
		
		kFArrays_nw.kT[i-1]           =  ki ;
		kFArrays_nw.P22ddT[i-1]       =  qrs_nw.P22dd;
		kFArrays_nw.P22duT[i-1]       =  qrs_nw.P22du;
		kFArrays_nw.P22uuT[i-1]       =  qrs_nw.P22uu;
		// A 
		kFArrays_nw.I1udd1AT[i-1]     =  qrs_nw.I1udd1A;
		kFArrays_nw.I2uud1AT[i-1]     =  qrs_nw.I2uud1A;
		kFArrays_nw.I2uud2AT[i-1]     =  qrs_nw.I2uud2A;
		kFArrays_nw.I3uuu2AT[i-1]     =  qrs_nw.I3uuu2A;
		kFArrays_nw.I3uuu3AT[i-1]     =  qrs_nw.I3uuu3A;
		//  B plus C   
		kFArrays_nw.I2uudd1BpCT[i-1]  =  qrs_nw.I2uudd1BpC;
		kFArrays_nw.I2uudd2BpCT[i-1]  =  qrs_nw.I2uudd2BpC;
		kFArrays_nw.I3uuud2BpCT[i-1]  =  qrs_nw.I3uuud2BpC;
		kFArrays_nw.I3uuud3BpCT[i-1]  =  qrs_nw.I3uuud3BpC;
		kFArrays_nw.I4uuuu2BpCT[i-1]  =  qrs_nw.I4uuuu2BpC;
		kFArrays_nw.I4uuuu3BpCT[i-1]  =  qrs_nw.I4uuuu3BpC;
		kFArrays_nw.I4uuuu4BpCT[i-1]  =  qrs_nw.I4uuuu4BpC;
		//  Bias
		kFArrays_nw.Pb1b2T[i-1]       =  qrs_nw.Pb1b2;
		kFArrays_nw.Pb1bs2T[i-1]      =  qrs_nw.Pb1bs2;
		kFArrays_nw.Pb22T[i-1]        =  qrs_nw.Pb22;
		kFArrays_nw.Pb2s2T[i-1]       =  qrs_nw.Pb2s2;
		kFArrays_nw.Ps22T[i-1]        =  qrs_nw.Ps22;
		kFArrays_nw.Pb2thetaT[i-1]    =  qrs_nw.Pb2theta;
		kFArrays_nw.Pbs2thetaT[i-1]   =  qrs_nw.Pbs2theta;
		//
		kFArrays_nw.P13ddT[i-1]       =  qrs_nw.P13dd;
		kFArrays_nw.P13duT[i-1]       =  qrs_nw.P13du;
		kFArrays_nw.P13uuT[i-1]       =  qrs_nw.P13uu;
		kFArrays_nw.sigma32PSLT[i-1]  =  qrs_nw.sigma32PSL;
		kFArrays_nw.pklT[i-1]         =  pkl_nw;
		kFArrays_nw.fkT[i-1]          =  fk;	
    }

    /* debug PRF("[k_functions] EXIT OK\n"); */

    #undef PRF
}
#undef _K_LOGSPACED_


#define QROMBERG     qromo
#define KK  5

//Modificacion-LCDMfk
// Q and R functions quadrature
// kk is the inner integration moment: 
// kk = k * r, so usual notation kk = p
//~ global_kFs ki_functions(real eta, real ki)
global_kFs ki_functions(real ki, double kPKL[], double pPKL[], int nPKLT, double pPKL2[])
{
    int i, j;
    real pkl_k;
     
    real PSLA, PSLB, psl;
	real fk, fp, fkmp, pklp, pklkmp;
    real rmin, rmax;
    real r, deltar, r2, y, y2;
    real mumin, mumax;
    real x, w, x2;
    real psl1;
    int Nx, nGL;
    real ypi, dk;
    real *xxGL, *wwGL, *xGL, *wGL;
    real kmin, kmax, pmin, pmax;
	real AngleEvQ, S2evQ, G2evQ, F2evQ, G2evR, F2evR; 
	real Gamma2evR, Gamma2fevR, C3Gamma3, C3Gamma3f, G3K, F3K, AngleEvR; 
    real A, ApOverf0, CFD3, CFD3p;
    
    
    if (gd.kernels_beyond_eds == 1) {
		A=KA_LCDM; ApOverf0 = KAp_LCDM/gd.f0;
		CFD3 = KR1_LCDM;  CFD3p = KR1p_LCDM;
	} else {		
		A=1; ApOverf0 = 0;
		CFD3 = 1;  CFD3p = 1;
	}
    
	real P22dd_p =0.0, P22dd_A = 0.0, P22dd_B = 0.0;
    real P22du_p =0.0, P22du_A = 0.0, P22du_B = 0.0;
    real P22uu_p =0.0, P22uu_A = 0.0, P22uu_B = 0.0;
    
    real I1udd1tA_p =0.0, I1udd1tA_A = 0.0, I1udd1tA_B = 0.0;
    real I2uud1tA_p =0.0, I2uud1tA_A = 0.0, I2uud1tA_B = 0.0;
    real I2uud2tA_p =0.0, I2uud2tA_A = 0.0, I2uud2tA_B = 0.0;
    real I3uuu2tA_p =0.0, I3uuu2tA_A = 0.0, I3uuu2tA_B = 0.0;
    real I3uuu3tA_p =0.0, I3uuu3tA_A = 0.0, I3uuu3tA_B = 0.0;
    
    real I2uudd1BpC_p =0.0, I2uudd1BpC_A =0.0, I2uudd1BpC_B =0.0;
    real I2uudd2BpC_p =0.0, I2uudd2BpC_A =0.0, I2uudd2BpC_B =0.0 ;
    real I3uuud2BpC_p =0.0, I3uuud2BpC_A =0.0, I3uuud2BpC_B =0.0 ;
    real I3uuud3BpC_p =0.0, I3uuud3BpC_A =0.0, I3uuud3BpC_B =0.0 ;
    real I4uuuu2BpC_p =0.0, I4uuuu2BpC_A =0.0, I4uuuu2BpC_B =0.0 ;
    real I4uuuu3BpC_p =0.0, I4uuuu3BpC_A =0.0, I4uuuu3BpC_B =0.0 ;
    real I4uuuu4BpC_p =0.0, I4uuuu4BpC_A =0.0, I4uuuu4BpC_B =0.0 ;
      
    real  Pb1b2_p=0.0,     Pb1b2_A =0.0,     Pb1b2_B =0.0;
    real  Pb1bs2_p=0.0,    Pb1bs2_A =0.0,    Pb1bs2_B =0.0; 
    real  Pb22_p=0.0,      Pb22_A =0.0,      Pb22_B =0.0; 
    real  Pb2s2_p=0.0,     Pb2s2_A =0.0,     Pb2s2_B =0.0; 
    real  Ps22_p=0.0,      Ps22_A =0.0,      Ps22_B =0.0; 
    real  Pb2theta_p=0.0,  Pb2theta_A =0.0,  Pb2theta_B =0.0; 
    real  Pbs2theta_p=0.0, Pbs2theta_A =0.0, Pbs2theta_B =0.0; 

    real KP22dd, KP22du, KP22uu;
    real KI1udd1tA, KI2uud1tA, KI2uud2tA, KI3uuu2tA, KI3uuu3tA;
    real KI2uudd1BpC, KI2uudd2BpC, KI3uuud2BpC, KI3uuud3BpC; 
    real KI4uuuu2BpC, KI4uuuu3BpC, KI4uuuu4BpC;
    real KPb1b2, KPb1bs2, KPb22, KPb2s2, KPs22, KPb2theta, KPbs2theta;
    
	real P13dd_p =0.0, P13dd_A = 0.0, P13dd_B = 0.0;
    real P13du_p =0.0, P13du_A = 0.0, P13du_B = 0.0;
    real P13uu_p =0.0, P13uu_A = 0.0, P13uu_B = 0.0;

    real sigma32PSL_p =0.0, sigma32PSL_A = 0.0, sigma32PSL_B = 0.0;
    
    real I1udd1a_p =0.0, I1udd1a_A = 0.0, I1udd1a_B = 0.0;
    real I2uud1a_p =0.0, I2uud1a_A = 0.0, I2uud1a_B = 0.0;
    real I2uud2a_p =0.0, I2uud2a_A = 0.0, I2uud2a_B = 0.0;
    real I3uuu2a_p =0.0, I3uuu2a_A = 0.0, I3uuu2a_B = 0.0;
    real I3uuu3a_p =0.0, I3uuu3a_A = 0.0, I3uuu3a_B = 0.0;
    
    real KP13dd, KP13du, KP13uu, Ksigma32PSL;
    real KI1udd1a, KI2uud1a, KI2uud2a, KI3uuu2a, KI3uuu3a;
    
    real *kk, *dkk;
    pointPSTableptr p;
    //
    /* store result on the stack to avoid a heap leak */
    global_kFs QRstmp_stack;
    global_kFs_ptr QRstmp = &QRstmp_stack;
 
    kmin = kPS[1];
    kmax = kPS[nPSLT];
    pmin = MAX(kmin,0.01*cmd.kmin);
    pmax = MIN(kmax,16.0*cmd.kmax);
   
    
    dk = (rlog10(pmax) - rlog10(pmin))/((real)(cmd.nquadSteps - 1));
    kk=dvector(1,cmd.nquadSteps);
    dkk=dvector(1,cmd.nquadSteps);
    kk[1] = rpow(10.0,rlog10(pmin));
    for (i=2; i<cmd.nquadSteps; i++) {
        ypi = rlog10(pmin) + dk*((real)(i - 1));
        kk[i] = rpow(10.0,ypi);
        dkk[i] = (kk[i]-kk[i-1]);
    }

// Q functions

    PSLA = 0.0;
    rmax = kmax/ki;
    rmin = kmin/ki;
    Nx=10;
    xxGL=dvector(1,Nx);
    wwGL=dvector(1,Nx);
    
	for (i=2; i<cmd.nquadSteps; i++) {
		r = kk[i]/ki;
        r2= r*r;
		PSLB = Interpolation_nr(kk[i], kPKL, pPKL, nPKLT, pPKL2);
		
		pklp=PSLB;
		fp = Interpolation_nr(kk[i], kPS, fkT, nPSLT, fkT2);
		fp /= gd.f0;
		mumin = MAX( -1.0, (1.0 + rsqr(r) - rsqr(rmax)) / (2.0*r)  );
		mumax = MIN(   1.0, (1.0  + rsqr(r) - rsqr(rmin)) / (2.0*r)  );
        
		if (r>=0.5)
			mumax = 0.5/r;
			gauleg(mumin,mumax,xxGL,wwGL,Nx);
        for (j=1; j<=Nx; j++) {
            x = xxGL[j];
            w = wwGL[j];
            x2= x*x;
            
            y2=1.0 + r2 - 2.0 * r * x;
            y = rsqrt(y2);
            psl = Interpolation_nr(ki * y, kPKL, pPKL, nPKLT, pPKL2); 
            pklkmp=psl;
			fkmp = Interpolation_nr(ki * y, kPS, fkT, nPSLT, fkT2);
			fkmp /= gd.f0;
			 
			AngleEvQ = (x - r)/y;
			S2evQ = AngleEvQ *AngleEvQ - 1./3.;
			F2evQ = 1./2. + 3./14. * A + (1./2. - 3./14. * A) * AngleEvQ*AngleEvQ + 
					AngleEvQ / 2. * (y/r + r/y);
			G2evQ = 3./14. * A * (fp + fkmp) + 3./14. * ApOverf0 
					+ (1./2. * (fp + fkmp) - 3./14. *  A *  (fp + fkmp) - 
					3./14. * ApOverf0)*AngleEvQ*AngleEvQ  
					+ AngleEvQ/2. *  (fkmp * y/r + fp * r/y);
			

			KP22dd = 2*r2*F2evQ*F2evQ;
			KP22du = 2*r2*F2evQ*G2evQ;
			KP22uu = 2*r2*G2evQ*G2evQ;
// A


			KI1udd1tA = 2.* (fp * r * x + fkmp * r2 * (1. - r * x)/y2 ) * F2evQ ;

			KI2uud1tA = - fp*fkmp * r2 *(1. - x2)/y2 * F2evQ;

			KI2uud2tA = 2.* (fp* r * x + fkmp *  r2*(1. - r * x)/y2) * G2evQ + 
						fp * fkmp * ( r2*(1. - 3.* x2) + 2.* r * x )/y2 * F2evQ;
						
			KI3uuu2tA = fp * fkmp *  r2 * (x2 - 1.)/y2 * G2evQ ;

			KI3uuu3tA = fp * fkmp * ( r2*(1. - 3.* x2) + 2.* r * x )/y2 * G2evQ;


     
// B+C     
			//~ KI2uudd1BpC =  fp*( fp * (1.-x2) + fkmp * r2 * (-1. + x2) / y2 ) / 2.;
			KI2uudd1BpC = 1/4. * (1.-x2)*(fp*fp + fkmp*fkmp*r2*r2/y2/y2) 
						+ fkmp*fp *r2 *(-1.+x2)/y2/2.; //Modificacion_Julio1
			
			KI2uudd2BpC = ( fp*fp*(-1. + 3.*x2) + 
				 2 * fkmp * fp*r * (r + 2*x - 3*r*x2) / y2 + 
				 fkmp * fkmp * r2*
				  ( 2 - 4*r*x + r2*(-1 + 3*x2) )/(y2*y2)    )/4. ;
  
			KI3uuud2BpC = -(   fkmp*fp*(fkmp*(-2 + 3*r*x)*r2 - 
				 fp*(-1 + 3*r*x)*(1 - 2*r*x + r2))*(-1 + x2 )    )/ (2.* y2*y2)  ; 
  
  			KI3uuud3BpC =   (  fkmp*fp*( -(fp*(1 - 2*r*x + r2)*(1 - 3*x2 + r*x*(-3 + 5*x2))) + 
				   fkmp*r*(2*x + r*(2 - 6*x2 + r*x*(-3 + 5*x2))) ) ) / (2.* y2*y2)  ;
				    
			KI4uuuu2BpC =  (3*rpow(fkmp,2)*rpow(fp,2)*r2*rpow(-1 + x2,2)) / (16. *y2*y2 )  ;         
				   
			KI4uuuu3BpC = -(rpow(fkmp,2)*rpow(fp,2)*(-1 + x2)*
				  (2 + 3*r*(-4*x + r*(-1 + 5*x2))))  / (8. *y2*y2 );
				   
			KI4uuuu4BpC =  (rpow(fkmp,2)*rpow(fp,2)*(-4 + 8*r*x*(3 - 5*x2) + 
				   12*x2 + r2*(3 - 30*x2 + 35*rpow(x,4))))  / (16. *y2*y2 );			
			
			       
  // Bias    
       
			KPb1b2  = r2 * F2evQ;
			KPb1bs2 = r2 * F2evQ*S2evQ;
			KPb22   = 1./2. * r2 * (1./2. * (1. - pklp/pklkmp) 
							+ 1./2. * (1. - pklkmp/pklp));
			KPb2s2  = 1./2. * r2 * (  1./2.* (S2evQ - 2./3. * pklp/pklkmp) 
							        + 1./2.* (S2evQ - 2./3. * pklkmp/pklp));
			KPs22   = 1./2. * r2 * (1./2. * (S2evQ * S2evQ - 4./9. * pklp/pklkmp) 
							+ 1./2.* (S2evQ * S2evQ - 4./9. * pklkmp/pklp));
			KPb2theta  = r2 * G2evQ;
			KPbs2theta = r2 * S2evQ * G2evQ;      
       
            //
    
            P22dd_B +=   w*KP22dd*psl;
            P22du_B +=   w*KP22du*psl;
            P22uu_B +=   w*KP22uu*psl;
            
            I1udd1tA_B +=   w*KI1udd1tA*psl;
            I2uud1tA_B +=   w*KI2uud1tA*psl;
            I2uud2tA_B +=   w*KI2uud2tA*psl;
            I3uuu2tA_B +=   w*KI3uuu2tA*psl;
            I3uuu3tA_B +=   w*KI3uuu3tA*psl;
          
            I2uudd1BpC_B +=   w*KI2uudd1BpC*psl;
            I2uudd2BpC_B +=   w*KI2uudd2BpC*psl;
            I3uuud2BpC_B +=   w*KI3uuud2BpC*psl;
            I3uuud3BpC_B +=   w*KI3uuud3BpC*psl;
            I4uuuu2BpC_B +=   w*KI4uuuu2BpC*psl;
            I4uuuu3BpC_B +=   w*KI4uuuu3BpC*psl;
            I4uuuu4BpC_B +=   w*KI4uuuu4BpC*psl;

            
            Pb1b2_B     +=   w*KPb1b2*psl;
            Pb1bs2_B    +=   w*KPb1bs2*psl;
            Pb22_B      +=   w*KPb22*psl;
            Pb2s2_B     +=   w*KPb2s2*psl;
            Ps22_B      +=   w*KPs22*psl;
            Pb2theta_B  +=   w*KPb2theta*psl;
            Pbs2theta_B +=   w*KPbs2theta*psl;

        }
        
        P22dd_p   += dkk[i]*(P22dd_A*PSLA + P22dd_B*PSLB)/2.0;
        P22du_p   += dkk[i]*(P22du_A*PSLA + P22du_B*PSLB)/2.0;
        P22uu_p   += dkk[i]*(P22uu_A*PSLA + P22uu_B*PSLB)/2.0;
        
        I1udd1tA_p   += dkk[i]*(I1udd1tA_A*PSLA + I1udd1tA_B*PSLB)/2.0;
        I2uud1tA_p   += dkk[i]*(I2uud1tA_A*PSLA + I2uud1tA_B*PSLB)/2.0;
        I2uud2tA_p   += dkk[i]*(I2uud2tA_A*PSLA + I2uud2tA_B*PSLB)/2.0;
        I3uuu2tA_p   += dkk[i]*(I3uuu2tA_A*PSLA + I3uuu2tA_B*PSLB)/2.0;
        I3uuu3tA_p   += dkk[i]*(I3uuu3tA_A*PSLA + I3uuu3tA_B*PSLB)/2.0;
        
        I2uudd1BpC_p   += dkk[i]*(I2uudd1BpC_A*PSLA + I2uudd1BpC_B*PSLB)/2.0;
        I2uudd2BpC_p   += dkk[i]*(I2uudd2BpC_A*PSLA + I2uudd2BpC_B*PSLB)/2.0;
        I3uuud2BpC_p   += dkk[i]*(I3uuud2BpC_A*PSLA + I3uuud2BpC_B*PSLB)/2.0;
        I3uuud3BpC_p   += dkk[i]*(I3uuud3BpC_A*PSLA + I3uuud3BpC_B*PSLB)/2.0;
        I4uuuu2BpC_p   += dkk[i]*(I4uuuu2BpC_A*PSLA + I4uuuu2BpC_B*PSLB)/2.0;
        I4uuuu3BpC_p   += dkk[i]*(I4uuuu3BpC_A*PSLA + I4uuuu3BpC_B*PSLB)/2.0;
        I4uuuu4BpC_p   += dkk[i]*(I4uuuu4BpC_A*PSLA + I4uuuu4BpC_B*PSLB)/2.0;

        Pb1b2_p       += dkk[i]*(Pb1b2_A    *PSLA + Pb1b2_B    *PSLB)/2.0;
        Pb1bs2_p      += dkk[i]*(Pb1bs2_A   *PSLA + Pb1bs2_B   *PSLB)/2.0;
        Pb22_p        += dkk[i]*(Pb22_A     *PSLA + Pb22_B     *PSLB)/2.0;
        Pb2s2_p       += dkk[i]*(Pb2s2_A    *PSLA + Pb2s2_B    *PSLB)/2.0;
        Ps22_p        += dkk[i]*(Ps22_A     *PSLA + Ps22_B     *PSLB)/2.0;
        Pb2theta_p    += dkk[i]*(Pb2theta_A *PSLA + Pb2theta_B *PSLB)/2.0;
        Pbs2theta_p   += dkk[i]*(Pbs2theta_A*PSLA + Pbs2theta_B*PSLB)/2.0;

     
        P22dd_A =   P22dd_B;      P22dd_B = 0.0;
        P22du_A =   P22du_B;      P22du_B = 0.0;
        P22uu_A =   P22uu_B;      P22uu_B = 0.0;
        
        I1udd1tA_A   =   I1udd1tA_B;        I1udd1tA_B   = 0.0;
        I2uud1tA_A   =   I2uud1tA_B;        I2uud1tA_B   = 0.0;
        I2uud2tA_A   =   I2uud2tA_B;        I2uud2tA_B   = 0.0;
        I3uuu2tA_A   =   I3uuu2tA_B;        I3uuu2tA_B   = 0.0;
        I3uuu3tA_A   =   I3uuu3tA_B;        I3uuu3tA_B   = 0.0;
         
        I2uudd1BpC_A =   I2uudd1BpC_B;      I2uudd1BpC_B = 0.0;
        I2uudd2BpC_A =   I2uudd2BpC_B;      I2uudd2BpC_B = 0.0;
        I3uuud2BpC_A =   I3uuud2BpC_B;      I3uuud2BpC_B = 0.0;
        I3uuud3BpC_A =   I3uuud3BpC_B;      I3uuud3BpC_B = 0.0;
        I4uuuu2BpC_A =   I4uuuu2BpC_B;      I4uuuu2BpC_B = 0.0;
        I4uuuu3BpC_A =   I4uuuu3BpC_B;      I4uuuu3BpC_B = 0.0;
        I4uuuu4BpC_A =   I4uuuu4BpC_B;      I4uuuu4BpC_B = 0.0;
        
        Pb1b2_A     =   Pb1b2_B;        Pb1b2_B     = 0.0;
        Pb1bs2_A    =   Pb1bs2_B;       Pb1bs2_B    = 0.0;
        Pb22_A      =   Pb22_B;         Pb22_B      = 0.0;
        Pb2s2_A     =   Pb2s2_B;        Pb2s2_B     = 0.0;
        Ps22_A      =   Ps22_B;         Ps22_B      = 0.0;
        Pb2theta_A  =   Pb2theta_B;     Pb2theta_B  = 0.0;
        Pbs2theta_A =   Pbs2theta_B;    Pbs2theta_B = 0.0;
       
        PSLA = PSLB;
    }
		
	P22dd_p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	P22du_p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	P22uu_p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
    
	I1udd1tA_p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I2uud1tA_p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I2uud2tA_p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I3uuu2tA_p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I3uuu3tA_p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
        
    I2uudd1BpC_p    *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I2uudd2BpC_p    *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I3uuud2BpC_p    *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I3uuud3BpC_p    *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I4uuuu2BpC_p    *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I4uuuu3BpC_p    *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I4uuuu4BpC_p    *= 2.0*(rpow(ki,3)/FOURPI2)/ki;

	Pb1b2_p        *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	Pb1bs2_p       *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	Pb22_p         *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	Pb2s2_p        *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	Ps22_p         *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	Pb2theta_p     *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	Pbs2theta_p    *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
    
    free_dvector(wwGL,1,Nx);
    free_dvector(xxGL,1,Nx);

//  R functions

    nGL=10;
    xGL=dvector(1,nGL);
    wGL=dvector(1,nGL);
    gauleg(-1.0,1.0,xGL,wGL,nGL);
    
    //~ fk = 1.0; // This is f(k)/f0 we need to interpolate f(k) 
	fk = Interpolation_nr(ki, kPS, fkT, nPSLT, fkT2);
	fk /= gd.f0;    
    for (i=2; i<cmd.nquadSteps; i++) {
        r = kk[i]/ki;
        r2= r*r;
        //~ psl = psInterpolation_nr(kk[i], kPS, pPS, nPSLT);
        psl = Interpolation_nr(kk[i], kPKL, pPKL, nPKLT, pPKL2); ;
        fp = Interpolation_nr(kk[i], kPS, fkT, nPSLT, fkT2);
		fp /= gd.f0; 
        for (j=1; j<=nGL; j++) {
            x = xGL[j];
            w = wGL[j];
            x2 =x*x;     
            y2=1.0 + r2 - 2.0 * r * x;
            
            
			Gamma2evR  = A *(1. - x2);
			Gamma2fevR = A *(1. - x2)*(fk + fp)/2. + 1./2. * ApOverf0 *(1 - x2);

			C3Gamma3  = 2.*5./21. * CFD3  *(1 - x2)*(1 - x2)/y2;
			C3Gamma3f = 2.*5./21. * CFD3p *(1 - x2)*(1 - x2)/y2 *(fk + 2 * fp)/3.;

			G3K = C3Gamma3f/ 2. + (2 * Gamma2fevR * x)/(7. * r) - (fk  * x2)/(6 * r2) 
				+ fp * Gamma2evR*(1 - r * x)/(7 * y2) 
				- 1./7.*(fp * Gamma2evR + 2 *Gamma2fevR) * (1. - x2)/y2;
			F3K = C3Gamma3/6. - x2/(6 * r2) + (Gamma2evR * x *(1 - r * x))/(7. *r *y2);
			   
			AngleEvR = -x;
			F2evR = 1./2. + 3./14. * A + (1./2. - 3./14. * A) * AngleEvR*AngleEvR 
				+ AngleEvR/2. *(1./r + r);

			G2evR = 3./14.* A *(fp + fk) + 3./14.* ApOverf0 
				+ ((fp + fk)/2. - 3./14.* A *(fp + fk) - 3./14.* ApOverf0)*AngleEvR*AngleEvR 
				+ AngleEvR/2. * (fk/r + fp * r);             
            
			KP13dd = 6.* r2 * F3K;
			KP13du = 3.* r2 * G3K + 3.* r2 * F3K * fk;
			KP13uu = 6.* r2 * G3K * fk;

			Ksigma32PSL = ( 5.0* r2 * (7. - 2*r2 + 4*r*x + 6*(-2 + r2)*x2 - 
				12*r*x2*x + 9*x2*x2)) / (24.0 * y2 ) ;
				
			KI1udd1a = 2.*r2*(1 - r * x)/y2 *G2evR + 2*fp *r* x* F2evR ;

			KI2uud1a = -fp * r2 * (1 - x2)/y2*G2evR ;

			KI2uud2a = ( (r2 *(1 - 3.*x2) + 2.* r* x) /y2*fp + 
				  fk*2.*r2*(1 - r * x)/y2)* G2evR + 2*x*r*fp*fk*F2evR;

			KI3uuu2a = fk * KI2uud1a;

			KI3uuu3a = (r2 *(1 - 3.*x2) + 2.* r* x)/y2 * fp*fk * G2evR;              

    
            P13dd_B +=   w*KP13dd*psl;
            P13du_B +=   w*KP13du*psl;
            P13uu_B +=   w*KP13uu*psl;
            
            sigma32PSL_B +=   w*Ksigma32PSL*psl;
            
            I1udd1a_B +=   w*KI1udd1a*psl;
            I2uud1a_B +=   w*KI2uud1a*psl;
            I2uud2a_B +=   w*KI2uud2a*psl;
            I3uuu2a_B +=   w*KI3uuu2a*psl;
            I3uuu3a_B +=   w*KI3uuu3a*psl;           
            
        }
        
        P13dd_p   += dkk[i]*(P13dd_A + P13dd_B) /  (2.0*ki);
        P13du_p   += dkk[i]*(P13du_A + P13du_B) /  (2.0*ki);
        P13uu_p   += dkk[i]*(P13uu_A + P13uu_B) /  (2.0*ki);
        
        sigma32PSL_p   += dkk[i]*(sigma32PSL_A + sigma32PSL_B) /  (2.0*ki);
        
        I1udd1a_p   += dkk[i]*(I1udd1a_A + I1udd1a_B) /  (2.0*ki);
        I2uud1a_p   += dkk[i]*(I2uud1a_A + I2uud1a_B) /  (2.0*ki);
        I2uud2a_p   += dkk[i]*(I2uud2a_A + I2uud2a_B) /  (2.0*ki);
        I3uuu2a_p   += dkk[i]*(I3uuu2a_A + I3uuu2a_B) /  (2.0*ki);
        I3uuu3a_p   += dkk[i]*(I3uuu3a_A + I3uuu3a_B) /  (2.0*ki);
        
              
        P13dd_A =   P13dd_B;      P13dd_B = 0.0;
        P13du_A =   P13du_B;      P13du_B = 0.0;
        P13uu_A =   P13uu_B;      P13uu_B = 0.0;
        sigma32PSL_A =   sigma32PSL_B;      sigma32PSL_B = 0.0;
        
        I1udd1a_A =   I1udd1a_B;      I1udd1a_B = 0.0;
        I2uud1a_A =   I2uud1a_B;      I2uud1a_B = 0.0;
        I2uud2a_A =   I2uud2a_B;      I2uud2a_B = 0.0;
        I3uuu2a_A =   I3uuu2a_B;      I3uuu2a_B = 0.0;
        I3uuu3a_A =   I3uuu3a_B;      I3uuu3a_B = 0.0;
      
    }

    pkl_k = Interpolation_nr(ki, kPKL, pPKL, nPKLT, pPKL2);
    P13dd_p      *= (rpow(ki,3.0)/FOURPI2)*pkl_k;
    P13du_p      *= (rpow(ki,3.0)/FOURPI2)*pkl_k;
    P13uu_p      *= (rpow(ki,3.0)/FOURPI2)*pkl_k;
    sigma32PSL_p *= (rpow(ki,3.0)/FOURPI2)*pkl_k;
    I1udd1a_p    *= (rpow(ki,3.0)/FOURPI2)*pkl_k;
    I2uud1a_p    *= (rpow(ki,3.0)/FOURPI2)*pkl_k;
    I2uud2a_p    *= (rpow(ki,3.0)/FOURPI2)*pkl_k;
    I3uuu2a_p    *= (rpow(ki,3.0)/FOURPI2)*pkl_k;
    I3uuu3a_p    *= (rpow(ki,3.0)/FOURPI2)*pkl_k;

	kFs(QRstmp)    = ki;
		
	P22dd(  QRstmp)      = P22dd_p;
	P22du(  QRstmp)      = P22du_p;
	P22uu(  QRstmp)      = P22uu_p;
	// A TNS
	I1udd1A(  QRstmp)      = I1udd1tA_p + 2.0*I1udd1a_p;
	I2uud1A(  QRstmp)      = I2uud1tA_p + 2.0*I2uud1a_p;
	I2uud2A(  QRstmp)      = I2uud2tA_p + 2.0*I2uud2a_p;
	I3uuu2A(  QRstmp)      = I3uuu2tA_p + 2.0*I3uuu2a_p;
	I3uuu3A(  QRstmp)      = I3uuu3tA_p + 2.0*I3uuu3a_p;
	// D function: B + C - G  
	I2uudd1BpC(  QRstmp)   = I2uudd1BpC_p 
								- ki*ki*gd.sigma2v*pkl_k;
	//~ I2uudd1BpC(  QRstmp)   = I2uudd1BpC_p;
	I2uudd2BpC(  QRstmp)   = I2uudd2BpC_p;
	I3uuud2BpC(  QRstmp)   = I3uuud2BpC_p
								- 2.0*ki*ki*gd.sigma2v*fk*pkl_k;
	//~ I3uuud2BpC(  QRstmp)   = I3uuud2BpC_p;
	I3uuud3BpC(  QRstmp)   = I3uuud3BpC_p;
	I4uuuu2BpC(  QRstmp)   = I4uuuu2BpC_p;
	I4uuuu3BpC(  QRstmp)   = I4uuuu3BpC_p
								- ki*ki*gd.sigma2v*fk*fk*pkl_k;
	//~ I4uuuu3BpC(  QRstmp)   = I4uuuu3BpC_p;
	I4uuuu4BpC(  QRstmp)   = I4uuuu4BpC_p;
	
	//  Bias
	Pb1b2(    QRstmp) = Pb1b2_p;
	Pb1bs2(   QRstmp) = Pb1bs2_p;
	Pb22(     QRstmp) = Pb22_p;
	Pb2s2(    QRstmp) = Pb2s2_p;
	Ps22(     QRstmp) = Ps22_p;
	Pb2theta( QRstmp) = Pb2theta_p;
	Pbs2theta(QRstmp) = Pbs2theta_p;
	//
	P13dd(  QRstmp)      = P13dd_p;
	P13du(  QRstmp)      = P13du_p;
	P13uu(  QRstmp)      = P13uu_p;
	
	sigma32PSL(QRstmp)   = sigma32PSL_p;
    
    /* free Gauss–Legendre buffers for R-part */
    free_dvector(wGL,1,nGL);
    free_dvector(xGL,1,nGL);

    /* free k-grid buffers */
    free_dvector(dkk,1,cmd.nquadSteps);
    free_dvector(kk,1,cmd.nquadSteps);
    
    /* return the stack object by value (no heap leak) */
    return QRstmp_stack;
}

// END Qs and Rs

//~ global_kFs qrs;
global_kFs kfunctions;
//~ global_kFs qrs_nw;

//~ global global_kFs ki_functions_driver(real eta, real ki)
//~ global global_kFs ki_functions_driver(real ki)
global global_kFs ki_functions_driver(real ki, double kPKL[], double pPKL[], int nPKLT, double pPKL2[])
{
    //~ quadrature(ki);
    kfunctions = ki_functions(ki,kPKL,pPKL,nPKLT,pPKL2);
    //~ return qrs;
    return kfunctions;
}


#define ROMO 1
#define NULLMETHOD 0
#define TRAPEZOID 2
#define TRAPEZOID3 5

void quadraturemethod_string_to_int(string method_str,int *method_int)
{
    *method_int=-1;
    if (strcmp(method_str,"romberg") == 0) {
        *method_int = ROMO;
        strcpy(gd.quadraturemethod_comment, "romberg open quadrature method");
    }
//
    if (strcmp(method_str,"trapezoid") == 0) {
        *method_int = TRAPEZOID;
        strcpy(gd.quadraturemethod_comment, "trapezoid quadrature method");
    }
//
    if (strcmp(method_str,"trapezoid3") == 0) {
        *method_int = TRAPEZOID3;
        strcpy(gd.quadraturemethod_comment, "trapezoid3 quadrature method");
    }
//
    if (strnull(method_str)) {
        *method_int = NULLMETHOD;
        strcpy(gd.quadraturemethod_comment,
               "given null quadrature method ... running deafult (trapezoid)");
        fprintf(stdout,"\n\tintegration: default integration method (trapezoid)...\n");
    }
//
    if (*method_int == -1) {
        *method_int = TRAPEZOID;
        strcpy(gd.quadraturemethod_comment,
               "Unknown quadrature method ... running deafult (trapezoid)");
        fprintf(stdout,"\n\tquadrature: Unknown method... %s ",cmd.quadratureMethod);
        fprintf(stdout,
                "\n\trunning default quadrature method (trapezoid)...\n");
    }
}

#undef ROMO
#undef TRAPEZOID
#undef TRAPEZOID3
#undef NULLMETHOD

#undef KK
#undef QROMBERG

local  real Interpolation_nr(real k, double kPS[], double pPS[], int nPS, double pPS2[])
{
    real psftmp;   
    splint(kPS,pPS,pPS2,nPS,k,&psftmp);  
    return (psftmp);
}

local real sigma28_function_int(real y)
{
    real p;
    real PSL,fk,j1, W;

    p = rpow(10.0,y);
    PSL = psInterpolation_nr(p, kPS, pPS, nPSLT);
    j1 = rj1Bessel(p*8.0);
    W = 3.0*j1/(p*8.0);
    
    return p*p*p*PSL*W*W;
}

local real sigma2L_function_int(real y)
{
    real p;
    real PSL;

    p = rpow(10.0,y);
    PSL = psInterpolation_nr(p, kPS, pPS, nPSLT);

    return p*PSL;
}

local real sigma2v_function_int(real y)
{
    real p;
    real PSL,fk;

    p = rpow(10.0,y);
    PSL = psInterpolation_nr(p, kPS, pPS, nPSLT);
    fk = Interpolation_nr(p, kPS, fkT, nPSLT, fkT2);
	fk /= gd.f0;
    return p*PSL*fk*fk;
}

local real Sigma2_int(real y)
{
    real p, PSL_nw, kosc;
    
    kosc=1.0/104.;
   
    p = rpow(10.0,y);
    PSL_nw = Interpolation_nr(p, kPS, pPS_nw, nPSLT,pPS2_nw);
    return p * PSL_nw * (1- rj0Bessel(p/kosc) + 2.*rj2Bessel(p/kosc) );
}

local real deltaSigma2_int(real y)
{
    real p, PSL_nw, kosc;
    
    kosc=1.0/104.;
   
    p = rpow(10.0,y);
    PSL_nw = Interpolation_nr(p, kPS, pPS_nw, nPSLT,pPS2_nw);
    return 3.0 * p * PSL_nw * rj2Bessel(p/kosc) ;
}

local real sigma_constants(void)
{
    real sigma2v, sigma2L, Sigma2, deltaSigma2;
    real kmin, kmax;
    real ymin, ymax, ymaxSigma;
    real EPSQ = 0.000001;
    int KK = 5;
    real ks;

    /* debug
    fprintf(stderr,"[sigma_constants] ENTER\n"); fflush(stderr);
    */

    kmin = kPS[1];
    kmax = kPS[nPSLT];
    ymin = rlog10(kmin);
    ymax = rlog10(kmax);

    /* debug
    fprintf(stderr,"[sigma_constants] ranges: k=[%g,%g] y=[%g,%g]\n", kmin,kmax,ymin,ymax);
    fprintf(stderr,"[sigma_constants] ptrs: kPS=%p pPS=%p pPS2=%p  fkT=%p fkT2=%p  pPS_nw=%p pPS2_nw=%p\n",
            (void*)kPS,(void*)pPS,(void*)pPS2,(void*)fkT,(void*)fkT2,(void*)pPS_nw,(void*)pPS2_nw);
    fflush(stderr);
    fprintf(stderr,"[sigma_constants] qromo sigma2v...\n"); fflush(stderr);
    */
    sigma2v = (1.0/SIXPI2)*rlog(10.0)*qromo(sigma2v_function_int,ymin,ymax,midpnt,EPSQ,KK);
    gd.sigma2v = sigma2v;
    /* debug fprintf(stderr,"[sigma_constants] sigma2v=%g\n", sigma2v); fflush(stderr); */

    /* debug fprintf(stderr,"[sigma_constants] qromo sigma2L...\n"); fflush(stderr); */
    sigma2L = (1.0/SIXPI2)*rlog(10.0)*qromo(sigma2L_function_int,ymin,ymax,midpnt,EPSQ,KK);
    gd.sigma2L = sigma2L;
    /* debug fprintf(stderr,"[sigma_constants] sigma2L=%g\n", sigma2L); fflush(stderr); */

    ks = 0.4;
    ymaxSigma = rlog10(ks);

    /* debug fprintf(stderr,"[sigma_constants] qromo Sigma2 (ymaxSigma=%g)...\n", ymaxSigma); fflush(stderr); */
    Sigma2 = (1.0/SIXPI2)*rlog(10.0)*qromo(Sigma2_int,ymin,ymaxSigma,midpnt,EPSQ,KK);
    gd.Sigma2 = Sigma2;
    /* debug fprintf(stderr,"[sigma_constants] Sigma2=%g\n", Sigma2); fflush(stderr); */

    /* debug fprintf(stderr,"[sigma_constants] qromo deltaSigma2 (ymaxSigma=%g)...\n", ymaxSigma); fflush(stderr); */
    deltaSigma2 = (1.0/SIXPI2)*rlog(10.0)*qromo(deltaSigma2_int,ymin,ymaxSigma,midpnt,EPSQ,KK);
    gd.deltaSigma2 = deltaSigma2;
    /* debug fprintf(stderr,"[sigma_constants] deltaSigma2=%g\n", deltaSigma2); fflush(stderr);
    fprintf(stderr,"[sigma_constants] EXIT OK\n"); fflush(stderr);
    */
}

local real get_sigma8(void)
{

    real sigma28, sigma8;
    real kmin, kmax;
    real ymin, ymax, ymaxSigma;
    real EPSQ = 0.000001;
    int KK = 5;
	real ks; 

    kmin = kPS[1];
    //~ kmax = kPS[nPSLT];
    kmax = 1.0;
    ymin = rlog10(kmin);
    ymax = rlog10(kmax);
	
	sigma28 = 3*(1.0/SIXPI2)*rlog(10.0)*
				qromo(sigma28_function_int,ymin,ymaxSigma,midpnt,EPSQ,KK);	
	gd.sigma8 = rsqrt(sigma28);  
};