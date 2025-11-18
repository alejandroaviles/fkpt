// pyfkpt/_fkptmodule.c  (DIAGNOSTIC VERSION WITH PRINTFs)
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#if defined(__linux__)
  #define _GNU_SOURCE
  #include <fcntl.h>
  #include <sys/mman.h>
#endif
#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef I
#undef I
#endif
#define write fkpt_write

#include "globaldefs.h"
#include "protodefs.h"
#include "models.h"
#include "getparam.h"
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

/* === FKPT internals we reference (1-based NR arrays!) === */
extern int    nPSLT;   /* length of FKPT's PS table */
extern double *kPS;    /* PS k-grid (dvector 1..nPSLT) */
/* Optional: access growth buffers for printing after set_external_fk */
extern double *fkT;
extern double *fkT2;

/* Quick, consistent prints */
#define PFX "[pyfkpt] "

#ifndef PYFKPT_DEBUG
#define PYFKPT_DEBUG 0   /* set to 1 (or compile with -DPYFKPT_DEBUG=1) to re-enable prints */
#endif

#if PYFKPT_DEBUG
  #define P(...) do{ fprintf(stderr, PFX __VA_ARGS__); fflush(stderr); }while(0)
#else
  #define P(...) do{}while(0)
#endif

#ifndef FKPT_EXPOSES_KERNELS
#define FKPT_EXPOSES_KERNELS 1
#endif
#ifndef FKPT_EXPOSES_NW
#define FKPT_EXPOSES_NW 1
#endif

/* needed prototypes from src/startrun.c */
void ReadParametersCmdline(void);
void StartRun(char *h0, char *h1, char *h2, char *h3);

/* --- tiny helpers --- */
static int is_1d_double(PyArrayObject *arr){
    return PyArray_NDIM(arr)==1 && PyArray_TYPE(arr)==NPY_DOUBLE;
}

/* Python API: compute_tables(...) */
static PyObject* py_compute_tables(PyObject* self, PyObject* args, PyObject* kwargs){
    P("ENTER compute_tables()\n");

    static char *kwlist[] = {
        "k","pk","z","Om","h",
        "b1","b2","bs2","b3nl",
        "alpha0","alpha2","alpha4",
        "ctilde","PshotP","alpha0shot","alpha2shot",
        "kmin","kmax","Nk","model","chatty","nquadSteps",
        /* keyword-only (all optional) */
        "fR0",
        "is_PS_input_LCDM",
        "mg1","mg2","mg3","mg4","mg5",
        "mg_variant",
        "mu0","c1","c2","Lambda",
        "beta_1","lambda_1","exp_s","beta_2","lambda_2",
        "mu1","mu2","mu3","mu4",
        "z_div","z_TGR","z_tw","k_tw","k_c",
        "fk","f0",
        "use_beyond_eds_kernels",
        NULL
    };

    /* required */
    PyObject *ok=NULL, *opk=NULL;
    double z, Om, h;

    /* galaxy bias / counterterms (defaults) */
    double b1=2.0, b2=0.0, bs2=0.0, b3nl=0.0;
    double alpha0=0.0, alpha2=0.0, alpha4=0.0;
    double ctilde=0.0, PshotP=0.0, alpha0shot=0.0, alpha2shot=0.0;

    /* k-range and integration controls */
    double kmin=0.0, kmax=0.0;
    int    Nk=12, chatty=0, nquadSteps=80;
    const char *model="LCDM";

    /* MG & knobs */
    double fR0 = 1e-10;
    int    is_PS_input_LCDM = 1;

    double mg1 = 1.0, mg2 = 0.0, mg3 = 0.0, mg4 = 0.0, mg5 = 0.0;
    const char *mg_variant = "baseline";

    double mu0 = 1.0;
    double c1 = 1.0, c2 = 1.0, Lambda = 0.0;
    double beta_1 = 1.0, lambda_1 = 0.0, exp_s = 1.0;
    double beta_2 = 1.0, lambda_2 = 0.0;

    double mu1 = 1.0, mu2 = 1.0, mu3 = 1.0, mu4 = 1.0;
    double z_div = 1.0, z_TGR = 2.0, z_tw = 0.05, k_tw = 0.001, k_c = 0.01;

    /* optional external f(k) and f0 */
    PyObject *ofk = NULL;
    double    f0_in = 0.0;

    /* kernels toggle */
    int use_beyond_eds_kernels = 0;

    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs,
            "OOddd"
            "|dddd"
            "ddd"
            "dddd"
            "dd"
            "i"
            "s"
            "i"
            "i"
            "$"
            "d"
            "i"
            "ddddd"
            "s"
            "dddd"
            "ddddd"
            "dddd"
            "ddddd"
            "O"
            "d"
            "p",
            kwlist,
            &ok,&opk,&z,&Om,&h,
            &b1,&b2,&bs2,&b3nl,
            &alpha0,&alpha2,&alpha4,
            &ctilde,&PshotP,&alpha0shot,&alpha2shot,
            &kmin,&kmax,&Nk,&model,&chatty,&nquadSteps,
            &fR0,
            &is_PS_input_LCDM, &mg1, &mg2, &mg3, &mg4, &mg5,
            &mg_variant,
            &mu0,&c1,&c2,&Lambda,
            &beta_1,&lambda_1,&exp_s,&beta_2,&lambda_2,
            &mu1,&mu2,&mu3,&mu4,
            &z_div,&z_TGR,&z_tw,&k_tw,&k_c,
            &ofk,&f0_in,
            &use_beyond_eds_kernels))
    {
        P("PyArg_ParseTupleAndKeywords FAILED\n");
        return NULL;
    }
    P("parsed args: z=%.6g Om=%.6g h=%.6g Nk=%d chatty=%d nquadSteps=%d model=%s\n",
      z,Om,h,Nk,chatty,nquadSteps, model);

    PyArrayObject *ak = (PyArrayObject*)PyArray_FROM_OTF(ok, NPY_DOUBLE, NPY_ARRAY_C_CONTIGUOUS);
    PyArrayObject *apk= (PyArrayObject*)PyArray_FROM_OTF(opk,NPY_DOUBLE, NPY_ARRAY_C_CONTIGUOUS);
    if(!ak || !apk){ P("FROM_OTF failed\n"); Py_XDECREF(ak); Py_XDECREF(apk); return NULL; }
    if(!is_1d_double(ak) || !is_1d_double(apk)){
        Py_DECREF(ak); Py_DECREF(apk);
        PyErr_SetString(PyExc_TypeError, "k and pk must be 1D double arrays");
        P("k/pk not 1D double\n");
        return NULL;
    }
    if(PyArray_SIZE(ak) != PyArray_SIZE(apk)){
        Py_DECREF(ak); Py_DECREF(apk);
        PyErr_SetString(PyExc_ValueError, "k and pk must have same length");
        P("k/pk size mismatch\n");
        return NULL;
    }

    int n = (int)PyArray_SIZE(ak);
    const double *k  = (const double*)PyArray_DATA(ak);
    const double *pk = (const double*)PyArray_DATA(apk);
    if (kmin <= 0.0) kmin = k[0];
    if (kmax <= 0.0) kmax = k[n-1];
    P("input arrays ok: n=%d  k[0]=%.6e  k[n-1]=%.6e  pk[0]=%.6e  pk[n-1]=%.6e\n",
      n, k[0], k[n-1], pk[0], pk[n-1]);

    char pathPS[256];
    int fd = -1;
    FILE *fp = NULL;

    #if defined(__linux__)
        fd = memfd_create("pyfkpt_ps", MFD_CLOEXEC);
        if (fd == -1) {
            Py_DECREF(ak); Py_DECREF(apk);
            PyErr_SetFromErrno(PyExc_RuntimeError);
            P("memfd_create failed: %s\n", strerror(errno));
            return NULL;
        }
        fp = fdopen(fd, "w+");
        if (!fp) {
            Py_DECREF(ak); Py_DECREF(apk);
            close(fd);
            PyErr_SetFromErrno(PyExc_RuntimeError);
            P("fdopen failed: %s\n", strerror(errno));
            return NULL;
        }
        for (int i = 0; i < n; i++) fprintf(fp, "%.17g %.17g\n", k[i], pk[i]);
        fflush(fp);
        fseek(fp, 0, SEEK_SET);
        snprintf(pathPS, sizeof(pathPS), "/proc/self/fd/%d", fd);
    #else
        fp = tmpfile();
        if (!fp) {
            Py_DECREF(ak); Py_DECREF(apk);
            PyErr_SetFromErrno(PyExc_RuntimeError);
            P("tmpfile failed: %s\n", strerror(errno));
            return NULL;
        }
        fd = fileno(fp);
        for (int i = 0; i < n; i++) fprintf(fp, "%.17g %.17g\n", k[i], pk[i]);
        fflush(fp);
        fseek(fp, 0, SEEK_SET);
        snprintf(pathPS, sizeof(pathPS), "/dev/fd/%d", fd);
    #endif
    P("PS table staged at %s (fd=%d)\n", pathPS, fd);

    char tokens[8192];
    #if defined(__linux__)
        const char *tmpdir = "/dev/shm";
    #else
        const char *tmpdir = "/tmp";
    #endif

    snprintf(tokens, sizeof(tokens),
        "paramfile=/dev/null "
        "tmpDir=%s clptDir=%s suffixModel=pyfkpt "
        "Om=%g h=%g zout=%g "
        "fnamePS=%s Nk=%d kmin=%g kmax=%g "
        "nquadSteps=%d model=%s chatty=%d "
        "b1=%g b2=%g bs2=%g b3nl=%g "
        "alpha0=%g alpha2=%g alpha4=%g "
        "ctilde=%g pshotp=%g alpha0shot=%g alpha2shot=%g "
        "fR0=%g is_PS_input_LCDM=%d "
        "mg1=%g mg2=%g mg3=%g mg4=%g mg5=%g mg_variant=%s "
        "mu0=%g c1=%g c2=%g Lambda=%g "
        "beta_1=%g lambda_1=%g exp_s=%g beta_2=%g lambda_2=%g "
        "mu1=%g mu2=%g mu3=%g mu4=%g "
        "z_div=%g z_TGR=%g z_tw=%g k_tw=%g k_c=%g",
        tmpdir, tmpdir,
        Om, h, z, pathPS, Nk, kmin, kmax,
        nquadSteps, model, chatty,
        b1,b2,bs2,b3nl,
        alpha0,alpha2,alpha4,
        ctilde,PshotP,alpha0shot,alpha2shot,
        fR0, is_PS_input_LCDM,
        mg1, mg2, mg3, mg4, mg5, mg_variant,
        mu0, c1, c2, Lambda,
        beta_1, lambda_1, exp_s, beta_2, lambda_2,
        mu1, mu2, mu3, mu4,
        z_div, z_TGR, z_tw, k_tw, k_c
    );
    P("tokens built (len=%zu)\n", strlen(tokens));

    /* Split tokens into 4 headlines */
    char h0[200] = "", h1[200] = "", h2[200] = "", h3[200] = "";
    {
        char *outs[4] = { h0, h1, h2, h3 };
        const size_t cap = sizeof(h0) - 1;
        int cur = 0;
        char *dup = strdup(tokens);
        if (!dup) { Py_DECREF(ak); Py_DECREF(apk); PyErr_NoMemory(); P("strdup tokens failed\n"); return NULL; }
        char *save = NULL;
        for (char *tok = strtok_r(dup, " ", &save); tok; tok = strtok_r(NULL, " ", &save)) {
            size_t need = (outs[cur][0] ? 1 : 0) + strlen(tok);
            if (strlen(outs[cur]) + need >= cap) {
                ++cur;
                if (cur >= 4) {
                    free(dup);
                    Py_DECREF(ak); Py_DECREF(apk);
                    PyErr_SetString(PyExc_RuntimeError, "argument string too long for fkpt headlines");
                    P("headlines overflow\n");
                    return NULL;
                }
            }
            if (outs[cur][0]) strncat(outs[cur], " ", cap - strlen(outs[cur]) - 1);
            strncat(outs[cur], tok, cap - strlen(outs[cur]) - 1);
        }
        free(dup);
    }
    P("h0='%s'\n", h0);
    P("h1='%s'\n", h1);
    P("h2='%s'\n", h2);
    P("h3='%s'\n", h3);

    char *hh0 = strdup(h0), *hh1 = strdup(h1), *hh2 = strdup(h2), *hh3 = strdup(h3);
    if (!hh0 || !hh1 || !hh2 || !hh3) {
        free(hh0); free(hh1); free(hh2); free(hh3);
        Py_DECREF(ak); Py_DECREF(apk);
        PyErr_NoMemory();
        P("strdup(h[0..3]) failed\n");
        return NULL;
    }

    /* argv-style InitParam */
    char *argv_dup = strdup(tokens);
    if(!argv_dup){ Py_DECREF(ak); Py_DECREF(apk); PyErr_NoMemory(); P("strdup tokens(2) failed\n"); return NULL; }
    size_t cap_tok = 64, argc = 1;
    char **argv_vec = (char**)calloc(cap_tok, sizeof(char*));
    if(!argv_vec){ free(argv_dup); Py_DECREF(ak); Py_DECREF(apk); PyErr_NoMemory(); P("calloc argv_vec failed\n"); return NULL; }
    argv_vec[0] = (char*)"pyfkpt";

    char *save2 = NULL;
    for(char *t = strtok_r(argv_dup, " ", &save2); t; t = strtok_r(NULL, " ", &save2)) {
        if(argc + 2 > cap_tok){
            size_t newcap = cap_tok * 2;
            char **tmp = (char**)realloc(argv_vec, newcap * sizeof(*argv_vec));
            if(!tmp){
                for(size_t i=1;i<argc;i++) free(argv_vec[i]);
                free(argv_vec); free(argv_dup);
                Py_DECREF(ak); Py_DECREF(apk);
                PyErr_NoMemory(); P("realloc argv_vec failed\n"); return NULL;
            }
            argv_vec = tmp; cap_tok = newcap;
        }
        argv_vec[argc] = strdup(t);
        if(!argv_vec[argc]){
            for(size_t i=1;i<argc;i++) free(argv_vec[i]);
            free(argv_vec); free(argv_dup);
            Py_DECREF(ak); Py_DECREF(apk);
            PyErr_NoMemory(); P("strdup argv token failed\n"); return NULL;
        }
        argc++;
    }
    argv_vec[argc] = NULL;
    P("InitParam(argc=%zu)\n", argc);

    static char *defv_vec[] = {
        "; defaults injected by pyfkpt wrapper",
        "paramfile=???", "tmpDir=/tmp", "clptDir=/tmp", "suffixModel=pyfkpt", "chatty=0",
        "b1=2","b2=0","bs2=0","b3nl=0",
        "alpha0=0","alpha2=0","alpha4=0",
        "ctilde=0","pshotp=0","alpha0shot=0","alpha2shot=0",
        "model=LCDM","modelParamfile=",
        "fR0=1e-10",
        "is_PS_input_LCDM=1",
        "mg1=1","mg2=0","mg3=0","mg4=0","mg5=0",
        "mg_variant=baseline",
        "mu0=1",
        "c1=1","c2=1","Lambda=0",
        "beta_1=1","lambda_1=0","exp_s=1","beta_2=1","lambda_2=0",
        "mu1=1","mu2=1","mu3=1","mu4=1",
        "z_div=1","z_TGR=2","z_tw=0.05","k_tw=0.001","k_c=0.01",
        "fnamePS=???","kmin=???","kmax=???","Nk=???",
        "Om=???","h=???","zout=???",
        "nquadSteps=80",
        NULL
    };

    InitParam(argv_vec, defv_vec);
    P("calling StartRun()\n");
    StartRun(hh0, hh1, hh2, hh3);
    P("StartRun() returned: nPSLT=%d kPS=%p fkT=%p fkT2=%p use_ext=%d f0=%.6e\n",
      nPSLT, (void*)kPS, (void*)fkT, (void*)fkT2, (int)gd.use_external_fk, gd.f0);
    if (nPSLT > 0 && kPS){
        /* NR arrays 1..nPSLT */
        P("kPS[1]=%.6e  kPS[nPSLT]=%.6e\n", kPS[1], kPS[nPSLT]);
    }

    /* optional external fk */
    PyArrayObject *afk = NULL;
    if (ofk && ofk != Py_None) {
        afk = (PyArrayObject*)PyArray_FROM_OTF(ofk, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
        if(!afk){
            for(size_t i=1;i<argc;i++) free(argv_vec[i]);
            free(argv_vec); free(argv_dup);
            free(hh0); free(hh1); free(hh2); free(hh3);
            Py_DECREF(ak); Py_DECREF(apk);
            P("FROM_OTF(fk) failed\n");
            return PyErr_NoMemory();
        }
        if (PyArray_NDIM(afk) != 1) {
            Py_DECREF(afk);
            for(size_t i=1;i<argc;i++) free(argv_vec[i]);
            free(argv_vec); free(argv_dup);
            free(hh0); free(hh1); free(hh2); free(hh3);
            Py_DECREF(ak); Py_DECREF(apk);
            PyErr_SetString(PyExc_TypeError, "fk must be a 1D double array");
            P("fk not 1D\n");
            return NULL;
        }
        if (PyArray_SIZE(afk) < 2) {
            Py_DECREF(afk);
            for(size_t i=1;i<argc;i++) free(argv_vec[i]);
            free(argv_vec); free(argv_dup);
            free(hh0); free(hh1); free(hh2); free(hh3);
            Py_DECREF(ak); Py_DECREF(apk);
            PyErr_SetString(PyExc_ValueError, "fk must have at least 2 points");
            P("fk too short\n");
            return NULL;
        }

        const double *k_in  = (const double*)PyArray_DATA(ak);
        const double *fk_in = (const double*)PyArray_DATA(afk);
        const int nin = (int)PyArray_SIZE(afk);
        P("about to set_external_fk(n=%d)  k_in=%p fk_in=%p  f0_in=%.6e\n",
          nin, (void*)k_in, (void*)fk_in, f0_in);
        set_external_fk(nin, k_in, fk_in, f0_in);
        P("returned from set_external_fk: nPSLT=%d kPS=%p fkT=%p fkT2=%p f0=%.6e use_ext=%d\n",
          nPSLT, (void*)kPS, (void*)fkT, (void*)fkT2, gd.f0, (int)gd.use_external_fk);
        if (nPSLT>0 && fkT){
            P("fkT[1]=%.6e fkT[nPSLT]=%.6e\n", fkT[1], fkT[nPSLT]);
        }

        /* (we also do a local resample + free; kept here for parity) */
        double *x1 = dvector(1, nin);
        double *y1 = dvector(1, nin);
        double *y2 = dvector(1, nin);
        for (int i=1; i<=nin; ++i) { x1[i] = k_in[i-1]; y1[i] = fk_in[i-1]; }
        spline(x1, y1, nin, 1.0e30, 1.0e30, y2);
        double *fk_res = dvector(1, nPSLT);
        const double slope_lo = (y1[2]   - y1[1])     / (x1[2]   - x1[1]);
        const double slope_hi = (y1[nin] - y1[nin-1]) / (x1[nin] - x1[nin-1]);
        for (int i=1; i<=nPSLT; ++i) {
            const double kk = kPS[i];
            double yy;
            if (kk < x1[1])       yy = y1[1]   + slope_lo * (kk - x1[1]);
            else if (kk > x1[nin])yy = y1[nin] + slope_hi * (kk - x1[nin]);
            else                  splint(x1, y1, y2, nin, kk, &yy);
            fk_res[i] = yy;
        }
        P("local resample done (fk_res[1]=%.6e, fk_res[n]=%.6e) — freeing temps\n", fk_res[1], fk_res[nPSLT]);
        free_dvector(fk_res, 1, nPSLT);
        free_dvector(y2, 1, nin);
        free_dvector(y1, 1, nin);
        free_dvector(x1, 1, nin);
        Py_DECREF(afk);
    } else {
        P("no external fk provided\n");
    }

    P("kernels_beyond_eds=%d -> gd.kernels_beyond_eds=%d\n",
      use_beyond_eds_kernels, use_beyond_eds_kernels ? 1:0);
    gd.kernels_beyond_eds = use_beyond_eds_kernels ? 1 : 0;

    P("about to compute_kfunctions()\n");
    compute_kfunctions();
    P("compute_kfunctions() returned\n");

    if (!kFArrays.kT) {
        Py_DECREF(ak); Py_DECREF(apk);
        PyErr_SetString(PyExc_RuntimeError, "fkpt did not produce k-grid (kFArrays.kT==NULL)");
        P("ERROR: kFArrays.kT is NULL\n");
        return NULL;
    }
    P("kFArrays.kT=%p; first=%.6e last=%.6e (Nk=%d)\n",
      (void*)kFArrays.kT, kFArrays.kT[0], kFArrays.kT[Nk-1], Nk);

    for(size_t i=1;i<argc;i++) free(argv_vec[i]);
    free(argv_vec);
    free(argv_dup);

    free(hh0); free(hh1); free(hh2); free(hh3);

    npy_intp dims[1] = { cmd.Nk };
    PyObject *out = PyDict_New();
    if(!out){
        Py_DECREF(ak); Py_DECREF(apk);
        P("PyDict_New failed\n");
        return NULL;
    }

    #define SAFE_PUT_VEC(name, src) do{                                         \
        if ((src)) {                                                            \
            PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);             \
            if(!arr){ Py_DECREF(out); Py_DECREF(ak); Py_DECREF(apk); P("alloc arr fail: %s\n", name); return NULL; } \
            double *dst = (double*)PyArray_DATA((PyArrayObject*)arr);           \
            for(int i=0;i<cmd.Nk;i++) dst[i] = (src)[i];                        \
            PyDict_SetItemString(out, (name), arr);                             \
            Py_DECREF(arr);                                                     \
        } else {                                                                \
            P("note: %s source is NULL\n", name);                               \
        }                                                                       \
    } while(0)

    #define SAFE_PUT_VEC_SCALE(name, src, scale) do{                            \
        if ((src)) {                                                            \
            PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);             \
            if(!arr){ Py_DECREF(out); Py_DECREF(ak); Py_DECREF(apk); P("alloc arr fail: %s\n", name); return NULL; } \
            double *dst = (double*)PyArray_DATA((PyArrayObject*)arr);           \
            for(int i=0;i<cmd.Nk;i++) dst[i] = (src)[i] * (scale);              \
            PyDict_SetItemString(out, (name), arr);                             \
            Py_DECREF(arr);                                                     \
        } else {                                                                \
            P("note: %s source is NULL\n", name);                               \
        }                                                                       \
    } while(0)

    #define SAFE_PUT_CONST(name, value) do{                                     \
        PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);                  \
        if(!arr){ Py_DECREF(out); Py_DECREF(ak); Py_DECREF(apk); P("alloc const fail: %s\n", name); return NULL; }  \
        double *dst = (double*)PyArray_DATA((PyArrayObject*)arr);                \
        for(int i=0;i<cmd.Nk;i++) dst[i] = (value);                              \
        PyDict_SetItemString(out, (name), arr);                                  \
        Py_DECREF(arr);                                                          \
    } while(0)

    #define SAFE_PUT_VEC_OR_ZERO(name, src) do{                                  \
        PyObject *arr__ = PyArray_SimpleNew(1, dims, NPY_DOUBLE);                \
        if(!arr__){ Py_DECREF(out); Py_DECREF(ak); Py_DECREF(apk); P("alloc zero fail: %s\n", name); return NULL; } \
        double *dst__ = (double*)PyArray_DATA((PyArrayObject*)arr__);            \
        if (src){ for(int i__=0;i__< cmd.Nk;i__++) dst__[i__] = (src)[i__]; }    \
        else     { for(int i__=0;i__< cmd.Nk;i__++) dst__[i__] = 0.0; }          \
        PyDict_SetItemString(out, (name), arr__);                                \
        Py_DECREF(arr__);                                                        \
    } while(0)

    /* ================= WIGGLE (W) ================= */
    SAFE_PUT_VEC_OR_ZERO("k",          kFArrays.kT);
    SAFE_PUT_VEC_OR_ZERO("pklin",      kFArrays.pklT);

    /* 1-loop SPT (W) */
    SAFE_PUT_VEC_OR_ZERO("P22dd",      kFArrays.P22ddT);
    SAFE_PUT_VEC_OR_ZERO("P22du",      kFArrays.P22duT);
    SAFE_PUT_VEC_OR_ZERO("P22uu",      kFArrays.P22uuT);
    SAFE_PUT_VEC_OR_ZERO("P13dd",      kFArrays.P13ddT);
    SAFE_PUT_VEC_OR_ZERO("P13du",      kFArrays.P13duT);
    SAFE_PUT_VEC_OR_ZERO("P13uu",      kFArrays.P13uuT);

    /* Growth (W) */
    SAFE_PUT_VEC_SCALE("f_over_f0",  kFArrays.fkT, 1.0/gd.f0);
    SAFE_PUT_VEC_OR_ZERO("Fk",       kFArrays.fkT);
    SAFE_PUT_CONST("f0_row",    gd.f0);

    /* === Export velocity-dispersion IR damping (scalar) === */
    SAFE_PUT_CONST("Sigma2",      gd.Sigma2);
    SAFE_PUT_CONST("deltaSigma2", gd.deltaSigma2);
    SAFE_PUT_CONST("sigma2v",     gd.sigma2v);
    SAFE_PUT_CONST("sigma2v_nw",  gd.sigma2v);     // same in C

    #if FKPT_EXPOSES_KERNELS
    SAFE_PUT_VEC_OR_ZERO("I1udd_1",    kFArrays.I1udd1AT);
    SAFE_PUT_VEC_OR_ZERO("I2uud_1",    kFArrays.I2uud1AT);
    SAFE_PUT_VEC_OR_ZERO("I2uud_2",    kFArrays.I2uud2AT);
    SAFE_PUT_VEC_OR_ZERO("I3uuu_2",    kFArrays.I3uuu2AT);
    SAFE_PUT_VEC_OR_ZERO("I3uuu_3",    kFArrays.I3uuu3AT);

    SAFE_PUT_VEC_OR_ZERO("I2uudd_1D",  kFArrays.I2uudd1BpCT);
    SAFE_PUT_VEC_OR_ZERO("I2uudd_2D",  kFArrays.I2uudd2BpCT);
    SAFE_PUT_VEC_OR_ZERO("I3uuud_2D",  kFArrays.I3uuud2BpCT);
    SAFE_PUT_VEC_OR_ZERO("I3uuud_3D",  kFArrays.I3uuud3BpCT);
    SAFE_PUT_VEC_OR_ZERO("I4uuuu_2D",  kFArrays.I4uuuu2BpCT);
    SAFE_PUT_VEC_OR_ZERO("I4uuuu_3D",  kFArrays.I4uuuu3BpCT);
    SAFE_PUT_VEC_OR_ZERO("I4uuuu_4D",  kFArrays.I4uuuu4BpCT);

    SAFE_PUT_VEC_OR_ZERO("Pb1b2",      kFArrays.Pb1b2T);
    SAFE_PUT_VEC_OR_ZERO("Pb1bs2",     kFArrays.Pb1bs2T);
    SAFE_PUT_VEC_OR_ZERO("Pb22",       kFArrays.Pb22T);
    SAFE_PUT_VEC_OR_ZERO("Pb2s2",      kFArrays.Pb2s2T);
    SAFE_PUT_VEC_OR_ZERO("Ps22",       kFArrays.Ps22T);
    SAFE_PUT_VEC_OR_ZERO("Pb2theta",   kFArrays.Pb2thetaT);
    SAFE_PUT_VEC_OR_ZERO("Pbs2theta",  kFArrays.Pbs2thetaT);
    SAFE_PUT_VEC_OR_ZERO("sigma32PSL", kFArrays.sigma32PSLT);
    #endif

    #if FKPT_EXPOSES_NW && FKPT_EXPOSES_KERNELS
    SAFE_PUT_VEC_OR_ZERO("I1udd_1_nw",    kFArrays_nw.I1udd1AT);
    SAFE_PUT_VEC_OR_ZERO("I2uud_1_nw",    kFArrays_nw.I2uud1AT);
    SAFE_PUT_VEC_OR_ZERO("I2uud_2_nw",    kFArrays_nw.I2uud2AT);
    SAFE_PUT_VEC_OR_ZERO("I3uuu_2_nw",    kFArrays_nw.I3uuu2AT);
    SAFE_PUT_VEC_OR_ZERO("I3uuu_3_nw",    kFArrays_nw.I3uuu3AT);

    SAFE_PUT_VEC_OR_ZERO("I2uudd_1D_nw",  kFArrays_nw.I2uudd1BpCT);
    SAFE_PUT_VEC_OR_ZERO("I2uudd_2D_nw",  kFArrays_nw.I2uudd2BpCT);
    SAFE_PUT_VEC_OR_ZERO("I3uuud_2D_nw",  kFArrays_nw.I3uuud2BpCT);
    SAFE_PUT_VEC_OR_ZERO("I3uuud_3D_nw",  kFArrays_nw.I3uuud3BpCT);
    SAFE_PUT_VEC_OR_ZERO("I4uuuu_2D_nw",  kFArrays_nw.I4uuuu2BpCT);
    SAFE_PUT_VEC_OR_ZERO("I4uuuu_3D_nw",  kFArrays_nw.I4uuuu3BpCT);
    SAFE_PUT_VEC_OR_ZERO("I4uuuu_4D_nw",  kFArrays_nw.I4uuuu4BpCT);

    SAFE_PUT_VEC_OR_ZERO("Pb1b2_nw",      kFArrays_nw.Pb1b2T);
    SAFE_PUT_VEC_OR_ZERO("Pb1bs2_nw",     kFArrays_nw.Pb1bs2T);
    SAFE_PUT_VEC_OR_ZERO("Pb22_nw",       kFArrays_nw.Pb22T);
    SAFE_PUT_VEC_OR_ZERO("Pb2s2_nw",      kFArrays_nw.Pb2s2T);
    SAFE_PUT_VEC_OR_ZERO("Ps22_nw",       kFArrays_nw.Ps22T);
    SAFE_PUT_VEC_OR_ZERO("Pb2theta_nw",   kFArrays_nw.Pb2thetaT);
    SAFE_PUT_VEC_OR_ZERO("Pbs2theta_nw",  kFArrays_nw.Pbs2thetaT);

    SAFE_PUT_VEC_OR_ZERO("sigma32PSL_nw", kFArrays_nw.sigma32PSLT);
    #endif

    #if FKPT_EXPOSES_NW
    SAFE_PUT_VEC_OR_ZERO("pklin_nw",     kFArrays_nw.pklT);
    SAFE_PUT_VEC_OR_ZERO("P22dd_nw",     kFArrays_nw.P22ddT);
    SAFE_PUT_VEC_OR_ZERO("P22du_nw",     kFArrays_nw.P22duT);
    SAFE_PUT_VEC_OR_ZERO("P22uu_nw",     kFArrays_nw.P22uuT);
    SAFE_PUT_VEC_OR_ZERO("P13dd_nw",     kFArrays_nw.P13ddT);
    SAFE_PUT_VEC_OR_ZERO("P13du_nw",     kFArrays_nw.P13duT);
    SAFE_PUT_VEC_OR_ZERO("P13uu_nw",     kFArrays_nw.P13uuT);
    SAFE_PUT_VEC_SCALE("f_over_f0_nw",   kFArrays_nw.fkT, 1.0/gd.f0);
    #endif

    Py_DECREF(ak);
    Py_DECREF(apk);
    if (fp) fclose(fp);
    P("EXIT compute_tables() OK\n");
    return out;
}

static PyMethodDef Methods[] = {
    {"compute_tables", (PyCFunction)py_compute_tables, METH_VARARGS | METH_KEYWORDS,
     "Run fkpt on a (k, Pk) table and return k-grid and k-function blocks."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_fkpt", 0, -1, Methods,
};

PyMODINIT_FUNC PyInit__fkpt(void){
    import_array();
    P("PyInit__fkpt()\n");
    return PyModule_Create(&moduledef);
}