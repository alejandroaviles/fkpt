/*==============================================================================
 MODULE: startrun.c
==============================================================================*/
/* expose POSIX prototypes like strdup() on glibc when compiling as C99 */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "globaldefs.h"
#include "protodefs.h"
#include <stdlib.h>   /* malloc/free, etc. */
#include <stdio.h>
#include <string.h>
#include <math.h>     /* isfinite, etc. */

/* debug macro (disabled) */
#ifndef STARTRUN_DEBUG
#define STARTRUN_DEBUG 0
#endif
#if STARTRUN_DEBUG
  #define SDBG(fmt, ...) do { fprintf(stderr, "[startrun] " fmt "\n", ##__VA_ARGS__); fflush(stderr); } while(0)
#else
  #define SDBG(...) do{}while(0)
#endif

/* at the top of startrun.c, near other file-scope statics */
static double *fkT_alloc  = NULL;  /* we free only this, never fkT directly */
static double *fkT2_alloc = NULL;

/* Use only the last component of a path (avoids slashes in output filenames) */
static const char *basename_noslash(const char *p){
    if (!p) return "ps";
    const char *b = strrchr(p, '/');
    return b ? (b + 1) : p;
}

/* Extract value for key=value from any of the 4 headlines */
static const char *extract_from_headlines(const char *key) {
    static char buf[256];
    const char *hs[4] = { gd.headline0, gd.headline1, gd.headline2, gd.headline3 };
    size_t klen = strlen(key);
    for (int i = 0; i < 4; ++i) {
        const char *h = hs[i];
        if (!h) continue;
        const char *p = h;
        while ((p = strstr(p, key)) != NULL) {
            if ((p == h || p[-1] == ' ') && p[klen] == '=') {
                p += klen + 1;
                size_t j = 0;
                while (p[j] && p[j] != ' ' && j < sizeof(buf)-1) { buf[j] = p[j]; j++; }
                buf[j] = '\0';
                return buf;
            }
            ++p;
        }
    }
    return NULL;
}

local void ReadParameterFile(char *);
local void PrintParameterFile(char *);

local void startrun_parameterfile(void);
local void startrun_cmdline(void);
void ReadParametersCmdline(void);
local void startrun_Common(void);
local void startrun_ParamStat(void);
local void CheckParameters(void);

/* Safe 2-column ASCII reader (1-based NR-style arrays) */
static int read_two_col_table(const char *path, double **x_out, double **y_out, int *n_out){
    FILE *fp = fopen(path, "r");
    #if STARTRUN_DEBUG
        fprintf(stderr, "[read_two_col_table] opening '%s'\n", path);
    #endif
    if (!fp) return -1;

    int cap = 1024, n = 0;
    double *tx = (double*)malloc(cap * sizeof(double));
    double *ty = (double*)malloc(cap * sizeof(double));
    if (!tx || !ty) { fclose(fp); free(tx); free(ty); return -2; }

    char line[4096];
    while (fgets(line, sizeof(line), fp)) {
        char *s = line;
        while (*s==' ' || *s=='\t') ++s;
        if (*s=='#' || *s=='\0' || *s=='\n') continue;
        double a,b;
        if (sscanf(s, "%lf %lf", &a, &b) == 2) {
            if (!isfinite(a) || !isfinite(b) || a <= 0.0 || b <= 0.0) continue;
            if (n == cap) {
                cap *= 2;
                double *nx = (double*)realloc(tx, cap * sizeof(double));
                double *ny = (double*)realloc(ty, cap * sizeof(double));
                if (!nx || !ny) { fclose(fp); free(tx); free(ty); return -3; }
                tx = nx; ty = ny;
            }
            tx[n] = a; ty[n] = b; n++;
        }
    }
    fclose(fp);
    if (n < 2) { free(tx); free(ty); return -4; }

    double *x = dvector(1, n);
    double *y = dvector(1, n);
    for (int i=1; i<=n; ++i) { x[i] = tx[i-1]; y[i] = ty[i-1]; }
    free(tx); free(ty);

    *x_out = x; *y_out = y; *n_out = n;
    return 0;
}

local void InputPSTable(void);
local void PSLTable(void);
local void PSLTableNW(void);
local void GaussLegendrePoints(void);
local void ClearRunState(void);

extern void model_reset_all(void);
extern void model_bind_globals(void);

/* Free any global allocations from a previous run in the same process */
local void ClearRunState(void) {

    /* 1) Nuke model/kernels cached pointers first */
    model_reset_all();

    /* Close log if open */
    if (gd.outlog) { fclose(gd.outlog); gd.outlog = NULL; }
    gd.logfilePath[0] = '\0';

    /* PS tables (log/linear, splines, no-wiggle) */
    if (PSLCDMtab)       { free(PSLCDMtab);       PSLCDMtab = NULL; }
    if (PSLCDMLogtab)    { free(PSLCDMLogtab);    PSLCDMLogtab = NULL; }
    if (PSLT)            { free(PSLT);            PSLT = NULL; }
    if (kPS)             { free_dvector(kPS,  1, nPSLT);  kPS  = NULL; }
    if (pPS)             { free_dvector(pPS,  1, nPSLT);  pPS  = NULL; }
    if (pPS2)            { free_dvector(pPS2, 1, nPSLT);  pPS2 = NULL; }
    if (pPS_nw)          { free_dvector(pPS_nw,  1, nPSLT); pPS_nw  = NULL; }
    if (pPS2_nw)         { free_dvector(pPS2_nw, 1, nPSLT); pPS2_nw = NULL; }

    /* growth/f(k) tables — free only what we actually allocated here */
    if (fkT_alloc)  { free_dvector(fkT_alloc,  1, nPSLT); fkT_alloc  = NULL; }
    if (fkT2_alloc) { free_dvector(fkT2_alloc, 1, nPSLT); fkT2_alloc = NULL; }
    /* always null the working pointers (they may have been overwritten to external memory) */
    fkT  = NULL;
    fkT2 = NULL;
    if (DplusT)          { free_dvector(DplusT,  1, nPSLT); DplusT = NULL; }
    if (DplusT2)         { free_dvector(DplusT2, 1, nPSLT); DplusT2 = NULL; }

    /* Gauss–Legendre */
    if (pGL) {
        if (xGL(pGL)) free_dvector(xGL(pGL), 1, nGL(pGL));
        if (wGL(pGL)) free_dvector(wGL(pGL), 1, nGL(pGL));
        free(pGL);
        pGL = NULL;
    }

    /* hard reset the public outputs (avoid stale pointers) */
    memset(&kFArrays,     0, sizeof(kFArrays));
    memset(&kFArraysd,    0, sizeof(kFArraysd));
    memset(&kFArrays_nw,  0, sizeof(kFArrays_nw));
    memset(&kFArraysd_nw, 0, sizeof(kFArraysd_nw));

    /* Reset counters so next allocation sizes are clean */
    nPSTable = 0; nPSLogT = 0; nPSLT = 0;

    /* Reset flags that might persist across runs */
    gd.use_external_fk = 0;
    gd.f0 = 0.0;
    /* release paramfile if we strdup'ed it last run */
    if (cmd.paramfile && *cmd.paramfile) {
        free((void*)cmd.paramfile);
    }
    cmd.paramfile = NULL;
    gd.headline0 = gd.headline1 = gd.headline2 = gd.headline3 = NULL;
}

void StartRun(string head0, string head1, string head2, string head3)
{
    real aTime = second();

    ClearRunState();
    gd.headline0 = head0; gd.headline1 = head1;
    gd.headline2 = head2; gd.headline3 = head3;

    SDBG("ENTER");
    SDBG("h0='%s'", head0 ? head0 : "(null)");
    SDBG("h1='%s'", head1 ? head1 : "(null)");
    SDBG("h2='%s'", head2 ? head2 : "(null)");
    SDBG("h3='%s'", head3 ? head3 : "(null)");

    /* DO NOT CALL GetParam() here — we haven't InitParam()'d. */
    const char *pf = extract_from_headlines("paramfile");
    if (pf && *pf && strcmp(pf, "/dev/null") != 0 && strcmp(pf, "NONE") != 0) {
        SDBG("paramfile from headlines -> '%s' (parameter-file branch)", pf);
        cmd.paramfile = strdup(pf);
        startrun_parameterfile();   /* reads file, sets cmd/gd */
    } else {
        SDBG("no real paramfile; using cmdline/headlines branch");
        cmd.paramfile = "";          /* mark as “no file” */
        startrun_cmdline();          /* must parse from headlines, not GetParam */
    }

    StartOutput();

    fprintf(gd.outlog, "\nStartRun elapsed time: %g sec.\n\n", second() - aTime);
    fflush(gd.outlog);
}

local void startrun_parameterfile(void)
{
    ReadParameterFile(cmd.paramfile);
    startrun_ParamStat();
    startrun_Common();
    PrintParameterFile(cmd.paramfile);
}

#define parameter_null "parameters_null-mgpt"

/* --- startrun_cmdline --- */
local void startrun_cmdline(void)
{
    SDBG("startrun_cmdline ENTER");
    ReadParametersCmdline();
    SDBG("startrun_cmdline after ReadParametersCmdline: "
         "fnamePS=%s Nk=%d kmin=%g kmax=%g model=%s",
         cmd.fnamePS ? cmd.fnamePS : "(null)", cmd.Nk, cmd.kmin, cmd.kmax,
         cmd.mgmodel ? cmd.mgmodel : "(null)");
    startrun_Common();
    PrintParameterFile(parameter_null);
    SDBG("startrun_cmdline EXIT");
}

/* --- ReadParametersCmdline --- */
void ReadParametersCmdline(void)
{
    SDBG("ReadParametersCmdline ENTER");

    /* helper macros so every read is bracketed by logs */
    #define STAT(key) SDBG("  GetParamStat(%s)=0x%x", key, GetParamStat(key))
    #define GETI(field, key) do{ \
        SDBG("  -> GetiParam(%s)", key); \
        cmd.field = GetiParam(key); \
        SDBG("     %s = %d", key, cmd.field); \
    }while(0)
    #define GETD(field, key) do{ \
        SDBG("  -> GetdParam(%s)", key); \
        cmd.field = GetdParam(key); \
        SDBG("     %s = %g", key, cmd.field); \
    }while(0)
    #define GETS(field, key) do{ \
        SDBG("  -> GetParam(%s)", key); \
        cmd.field = GetParam(key); \
        SDBG("     %s = '%s'", key, cmd.field ? cmd.field : "(null)"); \
    }while(0)

    const char *expect[] = {
        "chatty","b1","b2","bs2","b3nl","alpha0","alpha2","alpha4",
        "ctilde","pshotp","alpha0shot","alpha2shot",
        "model","suffixModel","modelParamfile","fR0","is_PS_input_LCDM",
        "fnamePS","kmin","kmax","Nk",
        "Om","h","zout","nquadSteps",
        /* MG params (optional) */
        "mg1","mg2","mg3","mg4","mg5",
        "mg_variant","mgv",
        "mu0","c1","c2","Lambda",
        "beta_1","lambda_1","exp_s","beta_2","lambda_2",
        "mu1","mu2","mu3","mu4",
        "z_div","z_TGR","z_tw","k_tw","k_c",
        NULL
    };
    for (int i=0; expect[i]; ++i) STAT((char*)expect[i]);

    GETI(chatty, "chatty");

    /* bias & cterms */
    GETD(b1,"b1"); GETD(b2,"b2"); GETD(bs2,"bs2"); GETD(b3nl,"b3nl");
    GETD(alpha0,"alpha0"); GETD(alpha2,"alpha2"); GETD(alpha4,"alpha4");
    GETD(ctilde,"ctilde"); GETD(PshotP,"pshotp");
    GETD(alpha0shot,"alpha0shot"); GETD(alpha2shot,"alpha2shot");

    /* fixed defaults */
    cmd.c1eft = 0.0; cmd.c2eft = 0.0; cmd.s2eft = 0.0;

    /* output arrays */
    cmd.smin = 1; cmd.smax = 130; cmd.Ns = 100;

    /* q-functions */
    cmd.NqperLogDecade = 100;
    cmd.Nk_qFunctionsQuad = 1200;

    /* GSM */
    cmd.gsm_width = 100.0;
    cmd.gsm_sizeyT = 160;
    cmd.gsm_NGL = 16;

    /* model */
    GETS(mgmodel,"model");
    GETS(suffixModel,"suffixModel");
    GETS(model_paramfile,"modelParamfile");
    GETD(fR0,"fR0");
    cmd.screening = 1.0;
    cmd.eps_DGP = -1.0; cmd.rc_DGP = 1.0;

    /* ---- HDKI sub-variant selector ---- */
    cmd.mg_variant = "baseline";
    if (GetParamStat("mg_variant") & ARGPARAM) cmd.mg_variant = GetParam("mg_variant");
    else if (GetParamStat("mgv") & ARGPARAM)   cmd.mg_variant = GetParam("mgv");
    SDBG("     mg_variant = '%s'", cmd.mg_variant ? cmd.mg_variant : "(null)");

    /* ---- MG general parameters: defaults + optional overrides ---- */
    /* mg1..mg5 */
    cmd.mg1 = 1.0; cmd.mg2 = 0.0; cmd.mg3 = 0.0; cmd.mg4 = 0.0; cmd.mg5 = 0.0;
    if (GetParamStat("mg1") & ARGPARAM) cmd.mg1 = GetdParam("mg1");
    if (GetParamStat("mg2") & ARGPARAM) cmd.mg2 = GetdParam("mg2");
    if (GetParamStat("mg3") & ARGPARAM) cmd.mg3 = GetdParam("mg3");
    if (GetParamStat("mg4") & ARGPARAM) cmd.mg4 = GetdParam("mg4");
    if (GetParamStat("mg5") & ARGPARAM) cmd.mg5 = GetdParam("mg5");

    /* mu0, c1, c2, Lambda */
    cmd.mu0    = 1.0;
    cmd.c1     = 1.0;
    cmd.c2     = 1.0;
    cmd.Lambda = 0.0;
    if (GetParamStat("mu0") & ARGPARAM)     cmd.mu0    = GetdParam("mu0");
    if (GetParamStat("c1") & ARGPARAM)      cmd.c1     = GetdParam("c1");
    if (GetParamStat("c2") & ARGPARAM)      cmd.c2     = GetdParam("c2");
    if (GetParamStat("Lambda") & ARGPARAM)  cmd.Lambda = GetdParam("Lambda");

    /* BZ-like time/scale dependence */
    cmd.beta_1   = 1.0;
    cmd.lambda_1 = 0.0;
    cmd.exp_s    = 1.0;
    cmd.beta_2   = 1.0;
    cmd.lambda_2 = 0.0;
    if (GetParamStat("beta_1") & ARGPARAM)   cmd.beta_1   = GetdParam("beta_1");
    if (GetParamStat("lambda_1") & ARGPARAM) cmd.lambda_1 = GetdParam("lambda_1");
    if (GetParamStat("exp_s") & ARGPARAM)    cmd.exp_s    = GetdParam("exp_s");
    if (GetParamStat("beta_2") & ARGPARAM)   cmd.beta_2   = GetdParam("beta_2");
    if (GetParamStat("lambda_2") & ARGPARAM) cmd.lambda_2 = GetdParam("lambda_2");

    /* binned mu/eta (or mu/Sigma) */
    cmd.mu1 = 1.0; cmd.mu2 = 1.0; cmd.mu3 = 1.0; cmd.mu4 = 1.0;
    cmd.z_div=1.0; cmd.z_TGR=2.0; cmd.z_tw=0.05; cmd.k_tw=0.001; cmd.k_c=0.01;
    if (GetParamStat("mu1") & ARGPARAM) cmd.mu1 = GetdParam("mu1");
    if (GetParamStat("mu2") & ARGPARAM) cmd.mu2 = GetdParam("mu2");
    if (GetParamStat("mu3") & ARGPARAM) cmd.mu3 = GetdParam("mu3");
    if (GetParamStat("mu4") & ARGPARAM) cmd.mu4 = GetdParam("mu4");
    if (GetParamStat("z_div") & ARGPARAM) cmd.z_div = GetdParam("z_div");
    if (GetParamStat("z_TGR") & ARGPARAM) cmd.z_TGR = GetdParam("z_TGR");
    if (GetParamStat("z_tw") & ARGPARAM) cmd.z_tw = GetdParam("z_tw");
    if (GetParamStat("k_tw") & ARGPARAM) cmd.k_tw = GetdParam("k_tw");
    if (GetParamStat("k_c") & ARGPARAM) cmd.k_c = GetdParam("k_c");

    /* PS table */
    cmd.is_PS_input_LCDM = 1;  /* sensible default */
    if (GetParamStat("is_PS_input_LCDM") & ARGPARAM)
        cmd.is_PS_input_LCDM = GetiParam("is_PS_input_LCDM");
    GETS(fnamePS,"fnamePS");
    GETD(kmin,"kmin"); GETD(kmax,"kmax"); GETI(Nk,"Nk");

    /* CLPT (fixed here) */
    cmd.rmin = 1.; cmd.rmax = 210.; cmd.Nr = 210;

    /* background */
    GETD(om,"Om");
    cmd.olstr = "1 - Om";
    GETD(h,"h");

    /* ODE / redshift */
    cmd.x = -4.0;
    cmd.dxstr = "2/5";
    cmd.dxmin = 0.0;
    cmd.eps   = 1.0e-4;
    GETD(xstop,"zout");
    cmd.maxnsteps = 10000; cmd.integration_method = "rkqs";

    /* quadrature */
    cmd.quadratureMethod = "trapezoid3";
    GETI(nquadSteps,"nquadSteps");
    cmd.ngausslegpoints = 16;
    cmd.epsquad = 1.0e-6;

    /* post */
    cmd.postprocessing = FALSE;
    cmd.options = "";

    if (cmd.chatty == 1) {
        printf("\n \t\t fkpt - compute power spectrum in modified gravity models \n\n");
        fprintf(stdout,"Reading pkl from file %s in the form column1: k[h/Mpc], 2: pkl[Mpc/h]^3 \n", cmd.fnamePS);
        fprintf(stdout,"redshift: z=%g\n",    cmd.xstop);
        fprintf(stdout,"OmegaM=%g\n",         cmd.om);
        fprintf(stdout,"h=%g\n",              cmd.h);
        fprintf(stdout,"b1=%g\n",             cmd.b1);
        fprintf(stdout,"b2=%g\n",             cmd.b2);
        fprintf(stdout,"bs2=%g\n",            cmd.bs2);
        fprintf(stdout,"b3nl=%g\n",           cmd.b3nl);
        fprintf(stdout,"alpha0=%g\n",         cmd.alpha0);
        fprintf(stdout,"alpha2=%g\n",         cmd.alpha2);
        fprintf(stdout,"alpha4=%g\n",         cmd.alpha4);
        fprintf(stdout,"ctilde=%g\n",         cmd.ctilde);
        fprintf(stdout,"Pshotp=%g\n",         cmd.PshotP);
        fprintf(stdout,"alpha0shot=%g\n",     cmd.alpha0shot);
        fprintf(stdout,"alpha2shot=%g\n",     cmd.alpha2shot);
    }

    SDBG("ReadParametersCmdline DONE");
    #undef STAT
    #undef GETI
    #undef GETD
    #undef GETS
}

#undef parameter_null

/* --- common setup --- */
local void startrun_Common(void)
{
    real dx1=0., dx2=0.;
    char *ep = NULL;

    SDBG("startrun_Common ENTER");

    /* zero path buffers if empty */
    gd.logfilePath[0] = gd.logfilePath[0] ? gd.logfilePath[0] : '\0';
    gd.tmpDir[0]      = gd.tmpDir[0]      ? gd.tmpDir[0]      : '\0';
    gd.clptDir[0]     = gd.clptDir[0]     ? gd.clptDir[0]     : '\0';

    SDBG("pre: tmpDir='%s' clptDir='%s' logfilePath='%s'",
         gd.tmpDir, gd.clptDir, gd.logfilePath);

    /* log paths */
    SDBG("call setFilesDirs_log()");
    setFilesDirs_log();
    SDBG("after setFilesDirs_log: tmpDir='%s' clptDir='%s' logfilePath='%s'",
         gd.tmpDir, gd.clptDir, gd.logfilePath);

    /* open log */
    strcpy(gd.mode,"w");
    SDBG("opening outlog='%s' mode='%s'", gd.logfilePath, gd.mode);
    if(!(gd.outlog=fopen(gd.logfilePath, gd.mode))) {
        strcpy(gd.logfilePath, "rk_fallback.log");
        gd.outlog = fopen(gd.logfilePath, gd.mode);
    }
    if (!gd.outlog) error("start_Common: error opening logfile (even fallback)");

    /* cosmology + dx */
    SDBG("cosmo: olstr='%s' om=%g h=%g", cmd.olstr, cmd.om, cmd.h);

    ep = strchr(cmd.olstr ? cmd.olstr : "", '-');
    if (ep == NULL) {
        gd.ol = GetdParam("OL");
        fprintf(gd.outlog,"\nOL string input without '-' :: %s\n",cmd.olstr ? cmd.olstr : "(null)");
        fprintf(gd.outlog,"\nOLambda, Om and sum  : %g %g %g\n",gd.ol, cmd.om, gd.ol+cmd.om);
    } else {
        ep = strchr(cmd.olstr, '1');
        if (ep == NULL)
            error("\nstart_Common: OL not in the format '1 - Om'\n");
        ep = strchr(cmd.olstr, 'O');
        if (ep == NULL || *(ep+1)!='m')
            error("\nstart_Common: OL not in the format '1 - Om'\n");
        gd.ol = 1. - cmd.om;
        fprintf(gd.outlog,"\nFound Om; OLambda and Om : %g %g\n",gd.ol, cmd.om);
    }

    SDBG("dxstr='%s'", cmd.dxstr);
    gd.dx = (sscanf(cmd.dxstr, "%lf/%lf", &dx1, &dx2) == 2 ? dx1/dx2 : atof(cmd.dxstr));
    if ( dx2 == 0. ) error("startrun_Common: dx : dx2 must be finite");

    /* checks + tables */
    SDBG("CheckParameters()");
    CheckParameters();

    SDBG("GaussLegendrePoints()");
    GaussLegendrePoints();

    SDBG("quadraturemethod_string_to_int()");
    quadraturemethod_string_to_int(cmd.quadratureMethod, &gd.quadmethod_int);

    SDBG("integration_method_string_to_int()");
    integration_method_string_to_int(cmd.integration_method, &gd.method_int);

    gd.xnow = cmd.x;
    gd.xout = gd.xnow;
    gd.xoutinfo = gd.xnow;

    gd.xstop = rlog(1.0/(1.0+cmd.xstop));  /* convert redshift -> ln a here */
    SDBG("xnow=%g xstop=%g (ln a)", gd.xnow, gd.xstop);

    SDBG("set_model()");
    set_model();

    SDBG("call setFilesDirs()");
    setFilesDirs();
    SDBG("after setFilesDirs: tmpDir='%s' clptDir='%s'", gd.tmpDir, gd.clptDir);
    /* NOTE: do not forcibly reset gd.use_external_fk here;
       if external fk will be used, set_external_fk() will flip it later. */

    /* PS input + derived tables */
    if (!strnull(cmd.fnamePS)) {
        /* If cmd.fnamePS is an absolute path or a proc/dev-fd handle, use it verbatim.
           Otherwise, treat it as a basename under Input/. */
        if (cmd.fnamePS[0] == '/' ||
            strncmp(cmd.fnamePS, "/proc/", 6) == 0 ||
            strncmp(cmd.fnamePS, "/dev/fd/", 8) == 0) {
            snprintf(gd.fnamePSPath, sizeof(gd.fnamePSPath), "%s", cmd.fnamePS);
        } else {
            snprintf(gd.fnamePSPath, sizeof(gd.fnamePSPath), "Input/%s", cmd.fnamePS);
        }
        SDBG("fnamePSPath='%s'", gd.fnamePSPath);

        SDBG("InputPSTable()");
        fprintf(gd.outlog, "[safe reader] fopen('%s')...\n", gd.fnamePSPath); fflush(gd.outlog);
        fprintf(gd.outlog, "\n\nReading power spectrum from file %s...\n", gd.fnamePSPath);
        InputPSTable();
        SDBG("InputPSTable done: nPSTable=%d nPSLogT=%d", nPSTable, nPSLogT);

        SDBG("PSLTable()");
        PSLTable();
        SDBG("PSLTable done: nPSLT=%d", nPSLT);

        SDBG("PSLTableNW()");
        PSLTableNW();
        SDBG("PSLTableNW done");
    } else {
        SDBG("fnamePS is empty -> skipping InputPSTable/PSL*");
    }

    SDBG("startrun_Common EXIT");
}

local void startrun_ParamStat(void)
{
    real dx1, dx2;
    /* output array params: */
    if (GetParamStat("smin") & ARGPARAM)
        cmd.smin = GetdParam("smin");
    if (GetParamStat("smax") & ARGPARAM)
        cmd.smax = GetdParam("smax");
    if (GetParamStat("Nk") & ARGPARAM)
        cmd.Nk = GetiParam("Nk");
    /* bias and counterterms */
    if (GetParamStat("b1") & ARGPARAM)
        cmd.b1 = GetdParam("b1");
    if (GetParamStat("b2") & ARGPARAM)
        cmd.b2 = GetdParam("b2");
    if (GetParamStat("bs2") & ARGPARAM)
        cmd.bs2 = GetdParam("bs2");
    if (GetParamStat("c1eft") & ARGPARAM)
        cmd.c1eft = GetdParam("c1eft");
    if (GetParamStat("c2eft") & ARGPARAM)
        cmd.c2eft = GetdParam("c2eft");
    if (GetParamStat("sigma2eft") & ARGPARAM)
        cmd.s2eft = GetdParam("sigma2eft");
    /* q functions: */
    if (GetParamStat("NqperLogDecade") & ARGPARAM)
        cmd.NqperLogDecade = GetiParam("NqperLogDecade");
    if (GetParamStat("Nk_qFunctionsQuad") & ARGPARAM)
        cmd.Nk_qFunctionsQuad = GetiParam("Nk_qFunctionsQuad");
    /* GSM: */
    if (GetParamStat("gsm_width") & ARGPARAM)
        cmd.gsm_width = GetdParam("gsm_width");
    if (GetParamStat("gsm_sizeyT") & ARGPARAM)
        cmd.gsm_sizeyT = GetiParam("gsm_sizeyT");
    if (GetParamStat("gsm_NGL") & ARGPARAM)
        cmd.gsm_NGL = GetiParam("gsm_NGL");
    /* PS table: */
    if (GetParamStat("fnamePS") & ARGPARAM)
        cmd.fnamePS = GetParam("fnamePS");
    if (GetParamStat("is_PS_input_LCDM") & ARGPARAM)
        cmd.is_PS_input_LCDM = GetiParam("is_PS_input_LCDM");
    if (GetParamStat("kmin") & ARGPARAM)
        cmd.kmin = GetdParam("kmin");
    if (GetParamStat("kmax") & ARGPARAM)
        cmd.kmax = GetdParam("kmax");
    if (GetParamStat("Nk") & ARGPARAM)
        cmd.Nk = GetiParam("Nk");
    /* CLPT correlation functions table: */
    cmd.rmin = 1.;
    cmd.rmax = 210.;
    cmd.Nr = 210;
    /* MG model: */
    cmd.mgmodel = "LCDM";
    if (GetParamStat("suffixModel") & ARGPARAM)
        cmd.suffixModel = GetParam("suffixModel");
    cmd.fR0 = 1.0e-10;
    cmd.screening = 1.0;
    /* DGP: */
    cmd.eps_DGP = -1.0;
    cmd.rc_DGP = 1.0;
    /* HDKI sub-variant (default + optional override) */
    cmd.mg_variant = "baseline";
    if (GetParamStat("mg_variant") & ARGPARAM)
        cmd.mg_variant = GetParam("mg_variant");
    if (GetParamStat("mgv") & ARGPARAM)
        cmd.mg_variant = GetParam("mgv");
    /* background: */
    if (GetParamStat("Om") & ARGPARAM)
        cmd.om = GetdParam("Om");
    cmd.olstr = "1 - Om";
    if (GetParamStat("h") & ARGPARAM)
        cmd.h = GetdParam("h");
    /* DE evolution: */
    cmd.x = -4.0;
    cmd.dxstr = "2/5";
    gd.dx = (sscanf(cmd.dxstr, "%lf/%lf", &dx1, &dx2) == 2 ?
                    dx1/dx2 : atof(cmd.dxstr));
    if ( dx2 == 0. )
        error("\n\nstartrun_ParamStat: deta : deta2 must be finite\n");

    if (GetParamStat("zout") & ARGPARAM) {
        double z = GetdParam("zout");
        cmd.xstop = z;  /* store redshift; startrun_Common() converts to ln a */
    }
    cmd.dxmin = (cmd.dxmin > 0 ? cmd.dxmin : 1e-4);

    cmd.eps = 0.0001;
    cmd.maxnsteps = 10000;

    cmd.integration_method = "rkqs";
    /* Quadrature parameters: */
    cmd.quadratureMethod = "trapezoid3";
    if (GetParamStat("nquadSteps") & ARGPARAM)
        cmd.nquadSteps = GetiParam("nquadSteps");
    cmd.ngausslegpoints = 16;
    cmd.epsquad = 1.0e-6;

    /* Post processing parameters: */
    cmd.postprocessing = FALSE;
    cmd.options = "";
}

local void CheckParameters(void)
{
    /* Power spectrum table: */
    if (strnull(cmd.fnamePS))
        error("CheckParameters: You should give a power spectrum filename\n");
    if (cmd.kmin < 0.0)
        error("CheckParameters: absurd value for kmin\n");
    if (cmd.kmax < 0.0)
        error("CheckParameters: absurd value for kmax\n");
    if (cmd.kmin > cmd.kmax)
        error("CheckParameters: kmin can not be greater than kmax\n");
    if (cmd.Nk < 0)
        error("CheckParameters: absurd value for Nk\n");
    /* CLPT correlation functions table: */
    if (cmd.rmin < 0.0)
        error("CheckParameters: absurd value for rmin\n");
    if (cmd.rmax < 0.0)
        error("CheckParameters: absurd value for rmax\n");
    if (cmd.rmin > cmd.rmax)
        error("CheckParameters: rmin can not be greater than rmax\n");
    if (cmd.Nr < 0)
        error("CheckParameters: absurd value for Nr\n");
    /* Background cosmology: */
    if (cmd.om > 1.0 || cmd.om < 0.0)
        error("CheckParameters: absurd value for om\n");
    if ( gd.ol < 0. )
        error("\n\nstartrun_ParamStat: OL (=%g) : must be positive\n",gd.ol);
    if (cmd.h < 0.0)
        error("CheckParameters: absurd value for h\n");
    /* Differential equations evolution parameters: */
    if (gd.dx == 0)
        error("CheckParameters: absurd value for deta\n");
    if (cmd.x == rlog(1.0/(1.0+cmd.xstop)))
        error("\n\nstartrun_Common: etaini and etaout=exp(-zout)-1 must be different\n");

    if (cmd.eps > 1.0e-4 || cmd.eps <= 0)
        error("CheckParameters: inapropriate or absurd value for epsSolver\n");
    if (cmd.maxnsteps < 1)
        error("CheckParameters: absurd value for maxnsteps\n");
    /* Quadrature parameters: */
    if (cmd.nquadSteps <= 1)
        error("CheckParameters: absurd value for nquadSteps\n");
    if (cmd.ngausslegpoints <= 1)
        error("CheckParameters: absurd value for ngausslegpoints\n");
    if (cmd.epsquad >= 1.0e-1 || cmd.epsquad <= 0)
        error("CheckParameters: absurd value for epsquad\n");
}

local void ReadParameterFile(char *fname)
{
#define DOUBLE 1
#define STRING 2
#define INT 3
#define BOOLEAN 4
#define MAXTAGS 300

  FILE *fd,*fdout;

  char buf[200],buf1[200],buf2[200],buf3[200];
  int  i,j,nt;
  int  id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int  errorFlag=0;

  nt=0;
    SPName(cmd.fnamePS,"fnamePS",100);
    IPName(cmd.is_PS_input_LCDM,"is_PS_input_LCDM");
    RPName(cmd.xstop,"zout");
    RPName(cmd.om,"Om");
    RPName(cmd.h,"h");
    /* bias and counterterms */
    RPName(cmd.b1,"b1");
    RPName(cmd.b2,"b2");
    RPName(cmd.bs2,"bs2");
    RPName(cmd.c1eft,"c1eft");
    RPName(cmd.c2eft,"c2eft");
    RPName(cmd.s2eft,"sigma2eft");
    /* output array params: */
    RPName(cmd.smin,"smin");
    RPName(cmd.smax,"smax");
    IPName(cmd.Ns,"Ns");
    /* q functions: */
    IPName(cmd.NqperLogDecade,"NqperLogDecade");
    IPName(cmd.Nk_qFunctionsQuad,"Nk_qFunctionsQuad");
    /* GSM: */
    RPName(cmd.gsm_width,"gsm_width");
    IPName(cmd.gsm_sizeyT,"gsm_sizeyT");
    IPName(cmd.gsm_NGL,"gsm_NGL");
    /* PS table: */
    RPName(cmd.kmin,"kmin");
    RPName(cmd.kmax,"kmax");
    IPName(cmd.Nk,"Nk");
    /* CLPT correlation functions table: */
    cmd.rmin = 1.;
    cmd.rmax = 210.;
    cmd.Nr = 210;
    /* MG model parameters: */
    SPName(cmd.suffixModel,"suffixModel",100);
    SPName(cmd.model_paramfile,"modelParamfile",100);
    SPName(cmd.mg_variant,"mg_variant",100);
    SPName(cmd.mg_variant,"mgv",100);
    cmd.mgmodel = "LCDM";
    cmd.fR0 = 1.0e-10;
    cmd.screening =1.0;
    /* DGP: */
    cmd.eps_DGP = -1.0;
    cmd.rc_DGP = 1.0;
    cmd.olstr ="1 - Om";

    /* ---- MG general parameters (defaults) ---- */
    cmd.mg1 = 1.0; cmd.mg2 = 0.0; cmd.mg3 = 0.0; cmd.mg4 = 0.0; cmd.mg5 = 0.0;
    RPName(cmd.mg1,"mg1");
    RPName(cmd.mg2,"mg2");
    RPName(cmd.mg3,"mg3");
    RPName(cmd.mg4,"mg4");
    RPName(cmd.mg5,"mg5");

    /* mu0, c1, c2, Lambda */
    cmd.mu0    = 1.0;
    cmd.c1     = 1.0;
    cmd.c2     = 1.0;
    cmd.Lambda = 0.0;
    RPName(cmd.mu0,"mu0");
    RPName(cmd.c1,"c1");
    RPName(cmd.c2,"c2");
    RPName(cmd.Lambda,"Lambda");

    /* BZ-like */
    cmd.beta_1   = 1.0;
    cmd.lambda_1 = 0.0;
    cmd.exp_s    = 1.0;
    cmd.beta_2   = 1.0;
    cmd.lambda_2 = 0.0;
    RPName(cmd.beta_1,"beta_1");
    RPName(cmd.lambda_1,"lambda_1");
    RPName(cmd.exp_s,"exp_s");
    RPName(cmd.beta_2,"beta_2");
    RPName(cmd.lambda_2,"lambda_2");

    /* binned */
    cmd.mu1 = 1.0; cmd.mu2 = 1.0; cmd.mu3 = 1.0; cmd.mu4 = 1.0;
    cmd.z_div=1.0; cmd.z_TGR=2.0; cmd.z_tw=0.05; cmd.k_tw=0.001; cmd.k_c=0.01;
    RPName(cmd.mu1,"mu1");
    RPName(cmd.mu2,"mu2");
    RPName(cmd.mu3,"mu3");
    RPName(cmd.mu4,"mu4");
    RPName(cmd.z_div,"z_div");
    RPName(cmd.z_TGR,"z_TGR");
    RPName(cmd.z_tw,"z_tw");
    RPName(cmd.k_tw,"k_tw");
    RPName(cmd.k_c,"k_c");

    /* DE evolution parameters: */
    cmd.dxstr = "2/5";
    cmd.eps = 0.0001;
    cmd.maxnsteps = 10000;
    /* Quadrature parameters: */
    IPName(cmd.nquadSteps,"nquadSteps");
    cmd.integration_method = "rkqs";
    cmd.quadratureMethod = "trapezoid3";
    cmd.ngausslegpoints = 16;
    cmd.epsquad = 1.0e-6;
    /* Post processing parameters: */
    cmd.postprocessing = FALSE;
    SPName(cmd.options,"options",100);

    if((fd=fopen(fname,"r"))) {
        while(!feof(fd)) {
            fgets(buf,200,fd);
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<1)
                continue;
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<2)
                *buf2='\0';
            if(buf1[0]=='%')
                continue;
            for(i=0,j=-1;i<nt;i++)
                if(strcmp(buf1,tag[i])==0) {
                    j=i;
                    tag[i][0]=0;
                    break;
                }
            if(j>=0) {
                switch(id[j]) {
                    case DOUBLE:
                        *((double*)addr[j])=atof(buf2);
                        break;
                    case STRING:
                        strcpy(addr[j],buf2);
                        break;
                    case INT:
                        *((int*)addr[j])=atoi(buf2);
                        break;
                    case BOOLEAN:
                        if (strchr("tTyY1", *buf2) != NULL) {
                            *((bool*)addr[j])=TRUE;
                        } else
                            if (strchr("fFnN0", *buf2) != NULL)  {
                                *((bool*)addr[j])=FALSE;
                            } else {
                                error("getbparam: %s=%s not bool\n",buf1,buf2);
                            }
                        break;
                }
            } else {
                fprintf(stdout, "Error in file %s: Tag '%s' %s.\n",
                    fname, buf1, "not allowed or multiple defined");
                errorFlag=1;
            }
        }
        fclose(fd);
    } else {
        fprintf(stdout,"Parameter file %s not found.\n", fname);
        errorFlag=1;
        exit(1);
    }

    for(i=0;i<nt;i++) {
        if(*tag[i]) {
            fprintf(stdout,
                "Error. I miss a value for tag '%s' in parameter file '%s'.\n",
                tag[i],fname);
            exit(0);
        }
    }
#undef DOUBLE
#undef STRING
#undef INT
#undef BOOLEAN
#undef MAXTAGS
}

#define FMTT "%-35s%s\n"
#define FMTI "%-35s%d\n"
#define FMTR "%-35s%g\n"

local void PrintParameterFile(char *fname)
{
    /* Disk writes disabled in pyfkpt wrapper */
    (void)fname;
    return;
}

#undef FMTT
#undef FMTI
#undef FMTR

#define NPT 10
#define SPREAD 1.0
local void InputPSTable(void)
{
    stream outstr;
    pointPSTableptr p, plog, pn;
    pointPSTableptr PSLCDMtabtmp;
    int nPSTabletmp;
    int i;
    real dk, kval, PSval, kmin, kmax;
    int mwt;
    double al,bl,chi2,q,siga,sigb,*x,*y,*sig;
    double au, bu;
    real *kPStmp = NULL;
    real *pPStmp = NULL;
    real *pPS2tmp = NULL;
    char namebuf[256];
    real kminext, kmaxext, dktmp, kmn, kmx;
    int Nkext=800, NkL=50, NkU=50;
    real kminT=1.0e-5, kmaxT=400.0;

    SDBG(" -> InputPSTable: about to call reader on '%s'", gd.fnamePSPath);
    fprintf(gd.outlog, "[safe reader] fopen('%s')...\n", gd.fnamePSPath); fflush(gd.outlog);
    fprintf(gd.outlog, "\n\nReading power spectrum from file %s...\n", gd.fnamePSPath);

    double *lx = NULL, *ly = NULL;
    int rc = read_two_col_table(gd.fnamePSPath, &lx, &ly, &nPSTabletmp);
    if (rc != 0) {
        SDBG("FATAL: read_two_col_table('%s') failed rc=%d", gd.fnamePSPath, rc);
        error("InputPSTable: failed to read '%s' (rc=%d)\n", gd.fnamePSPath, rc);
    }
    SDBG(" -> InputPSTable: read_two_col_table returned n=%d", nPSTabletmp);

    if (nPSTabletmp < 1)
        error("\n\nInputPSTable: nPSTable = %d is absurd\n\n", nPSTabletmp);

    PSLCDMtabtmp = (pointPSTableptr) allocate(nPSTabletmp * sizeof(pointPSTable));
    fprintf(gd.outlog,"nPSTable : %d\n", nPSTabletmp);

    i = 1;
    for (p=PSLCDMtabtmp; p<PSLCDMtabtmp+nPSTabletmp; ++p, ++i) {
        kPos(p) = lx[i];
        PS(p)   = ly[i];
    }
    free_dvector(lx, 1, nPSTabletmp);
    free_dvector(ly, 1, nPSTabletmp);

    fprintf(gd.outlog,"\n\nCreating log power spectrum...\n");

    PSLCDMLogtab = (pointPSTableptr) allocate(nPSTabletmp * sizeof(pointPSTable));
    plog = PSLCDMLogtab;
    nPSLogT=0;
    for (p=PSLCDMtabtmp; p<PSLCDMtabtmp+nPSTabletmp; p++) {
        kPos(plog) = rlog10(kPos(p));
        PS(plog) = rlog10(PS(p));
        plog++;
        nPSLogT++;
    }

    fprintf(gd.outlog,"\n\nTotal numbers in Log PS: %d %ld\n",nPSLogT,plog-PSLCDMLogtab);
    fprintf(gd.outlog,"Total numbers in Normal PS: %d %ld\n\n",nPSTabletmp,p-PSLCDMtabtmp);

    fprintf(gd.outlog,"\n\nLinear fit (a + b x) to log-log power spectrum at minset and maxset...\n");
    fprintf(gd.outlog,"\n\nLinear fit to log-log power spectrum at minset and maxset...\n");

    x=dvector(1,NPT);
    y=dvector(1,NPT);
    sig=dvector(1,NPT);

    /* Low-k */
    fprintf(gd.outlog,"\nAt low-ks of the spectrum...\n");

    plog = PSLCDMLogtab;
    for (i=1;i<=NPT;i++) {
        x[i]=kPos(plog);
        y[i]=PS(plog);
        sig[i]=SPREAD;
        plog++;
    }
    for (mwt=0;mwt<=1;mwt++) {
        fit(x,y,NPT,sig,mwt,&al,&bl,&siga,&sigb,&chi2,&q);
        if (mwt == 0)
            fprintf(gd.outlog,"\nIgnoring standard deviations\n");
        else
            fprintf(gd.outlog,"\nIncluding standard deviations\n");
        fprintf(gd.outlog,"%12s %9.6f %18s %9.6f \n",
               "a  =  ",al,"uncertainty:",siga);
        fprintf(gd.outlog,"%12s %9.6f %18s %9.6f \n",
               "b  =  ",bl,"uncertainty:",sigb);
        fprintf(gd.outlog,"%19s %14.6f \n","chi-squared: ",chi2);
        fprintf(gd.outlog,"%23s %10.6f \n","goodness-of-fit: ",q);
    }

    /* High-k */
    fprintf(gd.outlog,"\nAt high-ks of the spectrum...\n");

    plog = PSLCDMLogtab+nPSLogT-1;
    for (i=1;i<=NPT;i++) {
        x[i]=kPos(plog);
        y[i]=PS(plog);
        sig[i]=SPREAD;
        plog--;
    }
    for (mwt=0;mwt<=1;mwt++) {
        fit(x,y,NPT,sig,mwt,&au,&bu,&siga,&sigb,&chi2,&q);
        if (mwt == 0)
            fprintf(gd.outlog,"\nIgnoring standard deviations\n");
        else
            fprintf(gd.outlog,"\nIncluding standard deviations\n");
        fprintf(gd.outlog,"%12s %9.6f %18s %9.6f \n",
               "a  =  ",au,"uncertainty:",siga);
        fprintf(gd.outlog,"%12s %9.6f %18s %9.6f \n",
               "b  =  ",bu,"uncertainty:",sigb);
        fprintf(gd.outlog,"%19s %14.6f \n","chi-squared: ",chi2);
        fprintf(gd.outlog,"%23s %10.6f \n","goodness-of-fit: ",q);
    }

    /* Extend PS */
    kPStmp = dvector(1,nPSLogT);
    pPStmp = dvector(1,nPSLogT);
    pPS2 = dvector(1,nPSLogT);
    pn = PSLCDMLogtab;
    i=1;
    for (pn = PSLCDMLogtab; pn<PSLCDMLogtab+nPSLogT; pn++) {
        kPStmp[i] = kPos(pn);
        pPStmp[i] = PS(pn);
        i++;
    }
    spline(kPStmp,pPStmp,nPSLogT,1.0e30,1.0e30,pPS2);

    kmin = kPos(PSLCDMtabtmp);
    kmax = kPos(PSLCDMtabtmp+nPSTabletmp-1);
    fprintf(gd.outlog,
        "\nkmin, kmax of the given power spectrum (with %d values): %g %g",
        nPSTabletmp, kmin, kmax);
    dktmp = (rlog10(kmax) - rlog10(kmin))/((real)(nPSTabletmp - 1));
    kminext = rpow(10.0, rlog10(kmin)-((real)NkL)*dktmp);
    kmaxext = rpow(10.0, rlog10(kmax)+((real)NkU)*dktmp);

    fprintf(gd.outlog,"\n\nNkL, NkU: %d %d\n",NkL, NkU);
    NkL = rlog10(kmin/kminT)/dktmp;
    NkU = rlog10(kmaxT/kmax)/dktmp;
    fprintf(gd.outlog,"\n\nNkL, NkU targets: %d %d\n",NkL,NkU);
    kminext = rpow(10.0, rlog10(kmin)-((real)NkL)*dktmp);
    kmaxext = rpow(10.0, rlog10(kmax)+((real)NkU)*dktmp);

    fprintf(gd.outlog,"\nkmin, kmax of the extended power spectrum (first try): %g %g",
            kminext, kmaxext);

    kmn = MIN(kminext,cmd.kmin);
    kmx = MAX(kmaxext,cmd.kmax);
    fprintf(gd.outlog,"\nkmin, kmax of the extended power spectrum (second try): %g %g\n",
            kmn, kmx);

    nPSTable = Nkext;
    fprintf(gd.outlog,"\n\nCreating new PSTable with %d values\n",nPSTable);
    PSLCDMtab = (pointPSTableptr) allocate(nPSTable * sizeof(pointPSTable));
    dk = (rlog10(kmx) - rlog10(kmn))/((real)(nPSTable - 1));
    p = PSLCDMtab;

    for (i=1; i<=nPSTable; i++) {
        kval = rlog10(kmn) + dk*((real)(i - 1));
        if (rpow(10.0,kval) >= kmin && rpow(10.0,kval) <= kmax)
            PSval = psInterpolation_nr(kval, kPStmp, pPStmp, nPSLogT);
        else
            if (rpow(10.0,kval) < kmin)
                PSval = al + bl*kval;
            else
                if (rpow(10.0,kval) > kmax)
                    PSval = au + bu*kval;
                else
                    error("\n\nError: InputPSTable :: kmin, kmax, kval: %g %g %g", kmin, kmax, kval);
        kPos(p) = rpow(10.0,kval);
        PS(p) = rpow(10.0,PSval);
        p++;
    }

    free_dvector(pPS2,1,nPSLogT);
    free_dvector(pPStmp,1,nPSLogT);
    free_dvector(kPStmp,1,nPSLogT);

    free_dvector(sig,1,NPT);
    free_dvector(y,1,NPT);
    free_dvector(x,1,NPT);
    free(PSLCDMtabtmp);
}
#undef NPT
#undef SPREAD

local void PSLTable(void)
{
    char namebuf[256];
    stream outstr;
    real kmin, Dpkmin, Dpk;
    pointPSTableptr p, pn;
    int i;

    real xstoptmp, Dp0, Dpzout_LCDM, Dpzout, fac;

    xstoptmp = gd.xstop;
    gd.xstop = 0.;
    Dp0 = DpFunction_LCDM(0.);
    fprintf(gd.outlog,"\n\n Dp(0) = %g",Dp0);
    gd.xstop = xstoptmp;
    Dpzout_LCDM = DpFunction_LCDM(0.);
    Dpzout = DpFunction(0.);
    fprintf(gd.outlog,"\n Dp(%g) = %g\n",cmd.xstop,Dpzout);

    gd.Dplus=Dpzout/Dpzout_LCDM;

    kmin = kPos(PSLCDMtab);
    Dpkmin = DpFunction(kmin);
    fprintf(gd.outlog,"\n\n Dpkmin = %g\n",Dpkmin);

    PSLT = (pointPSTableptr) allocate(nPSTable * sizeof(pointPSTable));
    nPSLT = 0;
    pn = PSLT;
    if (cmd.is_PS_input_LCDM == 1) {
        if (strcasecmp(cmd.mgmodel, "LCDM") == 0) {
            // LCDM case: do not rescale, pass through
            for (p = PSLCDMtab; p < PSLCDMtab + nPSTable; p++) {
                kPos(pn) = kPos(p);
                PS(pn)   = PS(p);
                pn++; nPSLT++;
            }
        } else {
            const real A = rsqr(gd.Dplus);  // (Dp_MG(z_out)/Dp_LCDM(z_out))^2
            for (p = PSLCDMtab; p < PSLCDMtab + nPSTable; p++) {
                kPos(pn) = kPos(p);
                const real Dpk = DpFunction(kPos(p));   // MG D+(k)
                PS(pn) = A * rsqr(Dpk / Dpkmin) * PS(p);
                pn++; nPSLT++;
            }
        }
    } else {
        // No rescaling: pass input PS through
        for (p = PSLCDMtab; p < PSLCDMtab + nPSTable; p++) {
            kPos(pn) = kPos(p);
            PS(pn)   = PS(p);
            pn++; nPSLT++;
        }
    }

    kPS = dvector(1,nPSLT);
    pPS = dvector(1,nPSLT);
    pPS2 = dvector(1,nPSLT);
    pn = PSLT;
    i=1;
    for (pn = PSLT; pn<PSLT+nPSLT; pn++) {
        kPS[i] = kPos(pn);
        pPS[i] = PS(pn);
        i++;
    }

    fkT  = dvector(1,nPSLT);
    fkT2 = dvector(1,nPSLT);
    /* record the blocks we own so ClearRunState can free safely */
    fkT_alloc  = fkT;
    fkT2_alloc = fkT2;
    DplusT  = dvector(1,nPSLT);
    DplusT2 = dvector(1,nPSLT);

    if (!gd.use_external_fk) {
        for (int i2 = 1; i2 <= nPSLT; i2++) {
            fkT[i2]    = f_growth(kPS[i2]);
            DplusT[i2] = DpFunction(kPS[i2]);
        }
        /* Build spline for fkT only in internal mode */
        spline(kPS, fkT, nPSLT, 1.0e30, 1.0e30, fkT2);
        gd.f0 = fkT[1];                 /* internal default */
    } else {
        /* external path: set_external_fk() will fill fkT and fkT2 and set gd.f0 */
        /* if you also want DplusT here, compute it: */
        for (int i2 = 1; i2 <= nPSLT; i2++) {
            DplusT[i2] = DpFunction(kPS[i2]);
        }
    }

    /* bind current globals to the model layer for this run */
    model_bind_globals();

    /* always spline the PS itself */
    spline(kPS, pPS, nPSLT, 1.0e30, 1.0e30, pPS2);
}

real Interpolation_local(real k, double kPS[], double pPS[], int nPS, double pPS2[]);

local void PSLTableNW(void)
{
    char namebuf[256];
    stream outstr;
    double ksmin;
    double ksmax;
    double dk;
    int Nks;
    int cutmin;
    int cutmax;
    int i;
    double *kPSNW;
    double *pPSNW;
    double *pPSNW2;

    double *kPSNWtmp;
    double *pPSNWtmp;
    double *pPSNWtmp2;

    double *logkpkT;
    double *logkpkTeven;
    double *logkpkTodd;

    double *mlogkpkTevencutted;
    double *logkpkTevencutted;
    double *logkpkTevencutted2;

    double *mlogkpkToddcutted;
    double *logkpkToddcutted;
    double *logkpkToddcutted2;

    double *preT;

    double *logkpkTsave;
    double *logkpkTsave2;
    double PSL;
    int FSTtype;
    int FSTsizestep;

    int j;
    int Nksevencutted;
    int Nksoddcutted;

    double outkmin, outkmax, doutk;
    int outsize;
    double kval;
    double ki;
    double kPSNW8;
    double PNWval;
    double PSLkval;
    double DeltaAppf;
    double fac1;
    real aTime;
    aTime = second();

    outkmin = 0.00001;
    outkmax = 10;
    outsize = nPSLT;
    (void)outkmin; (void)outkmax; (void)doutk; /* currently unused */

    ksmin = 7.0e-5 / cmd.h;
    ksmax = 7.0 / cmd.h;
    Nks = 65536;            /* 2^16 */
    cutmin = 120;
    cutmax = 240;

    if (Nks==1)
        dk = 0.;
    else
        dk = (ksmax - ksmin)/((double)(Nks - 1));

    fprintf(gd.outlog," -> PSLTableNW: dk value: %g\n",dk);

    logkpkT = dvector(1,Nks);

    logkpkTsave = dvector(1,Nks);
    logkpkTsave2 = dvector(1,Nks);
    logkpkTeven = dvector(1,Nks/2);
    logkpkTodd = dvector(1,Nks/2);

    Nksevencutted = Nks/2 - (cutmax-cutmin)-4;
    Nksoddcutted = Nks/2 - (cutmax-cutmin)-1;

    fprintf(gd.outlog," -> PSTableNW : numbers discarded %d %d\n", Nksevencutted, Nksoddcutted);

    mlogkpkTevencutted = dvector(1,Nksevencutted);
    logkpkTevencutted = dvector(1,Nksevencutted);
    logkpkTevencutted2 = dvector(1,Nksevencutted);

    mlogkpkToddcutted = dvector(1,Nksoddcutted);
    logkpkToddcutted = dvector(1,Nksoddcutted);
    logkpkToddcutted2 = dvector(1,Nksoddcutted);

    preT = dvector(1,Nks);

    kPSNWtmp = dvector(1,Nks);
    pPSNWtmp = dvector(1,Nks);
    pPSNWtmp2 = dvector(1,Nks);

    kPSNW = dvector(1,outsize);
    pPSNW = dvector(1,outsize);
    pPSNW2 = dvector(1,outsize);

    for (i=1; i<=Nks; i++) {
        kval = ksmin + (double)(i-1)*dk;
        PSL = Interpolation_local(kval, kPS, pPS, nPSLT, pPS2);
        logkpkT[i] = rlog10(kval * PSL);
        logkpkTsave[i] = logkpkT[i];
    }

    dsinft(logkpkT, Nks);
    for (j=1; j<=Nks/2; j++) {
        logkpkTeven[j] = logkpkT[2*j];
        logkpkTodd[j] = logkpkT[2*j-1];
    }

    i=1;
    for (j=2; j<=cutmin-2; j++) {
        mlogkpkTevencutted[i] = (double)j;
        logkpkTevencutted[i] = logkpkTeven[j];
        ++i;
    }
    for (j=cutmax+2; j<=Nks/2; j++) {
        mlogkpkTevencutted[i] = (double)j;
        logkpkTevencutted[i] = logkpkTeven[j];
        ++i;
    }
    fprintf(gd.outlog," -> PSTableNW : even cutted %d\n", i-1);

    /* odd */
    i=1;
    for (j=1; j<=cutmin-1; j++) {
        mlogkpkToddcutted[j] = (double)j;
        logkpkToddcutted[j] = logkpkTodd[j];
        ++i;
    }
    for (j=cutmax+1; j<=Nks/2; j++) {
        mlogkpkToddcutted[i] = (double)j;
        logkpkToddcutted[i] = logkpkTodd[j];
        ++i;
    }
    fprintf(gd.outlog," -> PSTableNW : odd cutted %d\n", i-1);

    spline(mlogkpkTevencutted,logkpkTevencutted,Nksevencutted,1.0e30,1.0e30,logkpkTevencutted2);
    spline(mlogkpkToddcutted,logkpkToddcutted,Nksoddcutted,1.0e30,1.0e30,logkpkToddcutted2);

    double yy, dyy;
    unsigned long jj;

    for (i=1; i<=Nks/2; i++) {
        if(cutmin < i && i < cutmax) {
            locate(mlogkpkTevencutted, Nksevencutted, (double)(i+1), &jj);
            polint(&mlogkpkTevencutted[jj-2],&logkpkTevencutted[jj-2],5,(double)(i+1),&yy,&dyy);
            preT[2*i] = yy;
            locate(mlogkpkToddcutted, Nksoddcutted, (double)(i), &jj);
            polint(&mlogkpkToddcutted[jj-2],&logkpkToddcutted[jj-2],5,(double)(i),&yy,&dyy);
            preT[2*i-1] = yy;
        } else{
            preT[2*i] = logkpkT[2*i];
            preT[2*i-1] = logkpkT[2*i-1];
        }
    }

    dsinft(preT, Nks);

    fac1 = (2.0/(double) Nks);
    for (i=1; i<=Nks; i++) {
        preT[i] *= fac1;
    }

    for (i=1; i<=Nks; i++) {
        kPSNWtmp[i] = ksmin + (double)(i-1)*dk;
        pPSNWtmp[i] = rpow(10.0,preT[i])/kPSNWtmp[i];
    }
    spline(kPSNWtmp,pPSNWtmp,Nks,1.0e30,1.0e30,pPSNWtmp2);

    kPSNW8 = kPSNWtmp[8];

    pPS_nw   = dvector(1,nPSLT);
    pPS2_nw  = dvector(1,nPSLT);

    for (i=1; i<=outsize; i++) {
        ki = kPS[i];

        if (ki <= 0.001) {
            PSL       = Interpolation_local(kPSNW8, kPS, pPS, nPSLT, pPS2);
            PNWval    = Interpolation_local(kPSNW8, kPSNWtmp, pPSNWtmp, Nks, pPSNWtmp2);
            PSLkval   = Interpolation_local(ki, kPS, pPS, nPSLT, pPS2);
            DeltaAppf = (kPS[i]*(PSL - PNWval) / PNWval) /kPSNW8;
            pPS_nw[i] = PSLkval/(DeltaAppf + 1);
        } else if (ki > 0.001 && ki < ksmax){
            pPS_nw[i] = Interpolation_local(ki, kPSNWtmp, pPSNWtmp, Nks, pPSNWtmp2);
        } else {
            pPS_nw[i] = Interpolation_local(ki, kPS, pPS, nPSLT, pPS2);
        }
    }

    /* Prepare second-derivative table for no-wiggle interpolation */
    spline(kPS, pPS_nw, nPSLT, 1.0e30, 1.0e30, pPS2_nw);

    free_dvector(pPSNW2,1,outsize);
    free_dvector(pPSNW ,1,outsize);
    free_dvector(kPSNW ,1,outsize);

    free_dvector(pPSNWtmp2,1,Nks);
    free_dvector(pPSNWtmp,1,Nks);
    free_dvector(kPSNWtmp,1,Nks);

    free_dvector(preT,1,Nks);

    free_dvector(logkpkToddcutted2,1,Nksoddcutted);
    free_dvector(logkpkToddcutted,1,Nksoddcutted);
    free_dvector(mlogkpkToddcutted,1,Nksoddcutted);

    free_dvector(logkpkTevencutted2,1,Nksevencutted);
    free_dvector(logkpkTevencutted,1,Nksevencutted);
    free_dvector(mlogkpkTevencutted,1,Nksevencutted);

    free_dvector(logkpkTodd,1,Nks/2);
    free_dvector(logkpkTeven,1,Nks/2);
    free_dvector(logkpkTsave2,1,Nks);
    free_dvector(logkpkTsave,1,Nks);

    free_dvector(logkpkT,1,Nks);

    fprintf(gd.outlog," -> PSLTableNW CPU time %g...\n",second()-aTime);
}

real Interpolation_local(real k, double kPS[], double pPS[], int nPS, double pPS2[])
{
    real psftmp;
    splint(kPS,pPS,pPS2,nPS,k,&psftmp);
    return (psftmp);
}

local void GaussLegendrePoints(void)
{
    real x1=-1.0, x2=1.0;

    pGL = (global_GL_ptr) allocate(sizeof(global_GL));

    nGL(pGL) = cmd.ngausslegpoints;
    xGL(pGL)=dvector(1,cmd.ngausslegpoints);
    wGL(pGL)=dvector(1,cmd.ngausslegpoints);
    x1GL(pGL) = -1.0;
    x2GL(pGL) = 1.0;

    fprintf(gd.outlog,"\nComputation of GLs...\n");
    gauleg(x1GL(pGL),x2GL(pGL),xGL(pGL),wGL(pGL),nGL(pGL));
    fprintf(gd.outlog,"\nEnd of Gauss-Legendre computing (StartRun)\n");
}

global real psInterpolation_nr(real k, double kPS[], double pPS[], int nPS)
{
    real psftmp;
    splint(kPS,pPS,pPS2,nPS,k,&psftmp);
    return (psftmp);
}