#include "../src/globaldefs.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "stdinc.h"
#include "mathfns.h"
#include "diffeqs.h"
#include "numrec.h"

/* ------------------------- debugging helpers ------------------------- */

#ifndef DIFFEQS_DEBUG_ENV
#define DIFFEQS_DEBUG_ENV "FKPT_DIFFEQS_DEBUG"
#endif

#ifndef DIFFEQS_PRINTV_ENV
#define DIFFEQS_PRINTV_ENV "FKPT_DIFFEQS_PRINTV"
#endif

#ifndef RKCK_EVERY_ENV
#define RKCK_EVERY_ENV "FKPT_RKCK_EVERY"
#endif

static int diffeqs_dbg_enabled(void) {
    static int init=0, on=0;
    if (!init) {
        const char *e = getenv(DIFFEQS_DEBUG_ENV);
        on = (e && *e && *e != '0');
        init = 1;
    }
    return on;
}

static int diffeqs_printv(void) {
    static int init=0, n=3;
    if (!init) {
        const char *e = getenv(DIFFEQS_PRINTV_ENV);
        if (e && *e) {
            long t = strtol(e, NULL, 10);
            if (t > 0 && t < 10000) n = (int)t;
        }
        init = 1;
    }
    return n;
}

static long rkck_every(void) {
    static int init=0;
    static long v=0; /* 0 => off (no prints from rkck) */
    if (!init) {
        const char *e = getenv(RKCK_EVERY_ENV);
        v = (e && *e) ? strtol(e, NULL, 10) : 0;
        if (v < 0) v = 0;
        init = 1;
    }
    return v;
}

#define DDBG(fmt, ...) \
    do { if (diffeqs_dbg_enabled()) fprintf(stderr, "[diffeqs:%s] " fmt "\n", __func__, ##__VA_ARGS__); } while(0)

static void print_vec(const char *tag, const double v[], int n) {
    if (!diffeqs_dbg_enabled()) return;
    int m = diffeqs_printv();
    if (m > n) m = n;
    fprintf(stderr, "[diffeqs:%s] %s[1:%d]:", __func__, tag, m);
    for (int i=1; i<=m; ++i) fprintf(stderr, " % .3e", v[i]);
    if (n > m) fprintf(stderr, " ...");
    fputc('\n', stderr);
}

static void stats_vec(const char *tag, const double v[], int n) {
    if (!diffeqs_dbg_enabled()) return;
    double vmin = +HUGE_VAL, vmax = -HUGE_VAL, vavg = 0.0;
    for (int i=1; i<=n; ++i) {
        if (v[i] < vmin) vmin = v[i];
        if (v[i] > vmax) vmax = v[i];
        vavg += v[i];
    }
    vavg /= (double)n;
    fprintf(stderr, "[diffeqs:%s] %s{min=%.3e max=%.3e mean=%.3e}\n", __func__, tag, vmin, vmax, vavg);
}

/* ------------------------------ rk4 ---------------------------------- */

void rk4(double y[], double dydx[], int n, double x, double h, double yout[],
               void (*derivs)(double, double [], double []))
{
    int i;
    double xh,hh,h6,*dym,*dyt,*yt;

    // if (diffeqs_dbg_enabled()) {
    //     DDBG("ENTER x=% .15e h=% .3e n=%d", x, h, n);
    //     print_vec("y", y, n);
    //     print_vec("dydx", dydx, n);
    // }

    dym=dvector(1,n);
    dyt=dvector(1,n);
    yt=dvector(1,n);
    hh=h*0.5;
    h6=h/6.0;
    xh=x+hh;
    for (i=1;i<=n;i++) yt[i]=y[i]+hh*dydx[i];
    (*derivs)(xh,yt,dyt);
    for (i=1;i<=n;i++) yt[i]=y[i]+hh*dyt[i];
    (*derivs)(xh,yt,dym);
    for (i=1;i<=n;i++) {
        yt[i]=y[i]+h*dym[i];
        dym[i] += dyt[i];
    }
    (*derivs)(x+h,yt,dyt);
    for (i=1;i<=n;i++)
        yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);

    // if (diffeqs_dbg_enabled()) {
    //     print_vec("yout", yout, n);
    //     DDBG("EXIT");
    // }
    free_dvector(yt,1,n);
    free_dvector(dyt,1,n);
    free_dvector(dym,1,n);
}


/* ------------------------------ odeint ------------------------------- */

#define MAXSTP 10000
#define TINY 1.0e-30

void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
            double hmin, int *nok, int *nbad, int maxnsteps,
            void (*derivsin)(double, double [], double []),
            void (*rkqsin)(double [], double [], int, double *, double, double, double [],
                         double *, double *, void (*)(double, double [], double [])))
{
    int nstp,i;
    double xsav,x,hnext,hdid,h;
    double *yscal,*y,*dydx;

    // if (diffeqs_dbg_enabled()) {
    //     DDBG("ENTER nvar=%d x1=% .15e x2=% .15e eps=%.3e h1=%.3e hmin=%.3e maxnsteps=%d",
    //          nvar, x1, x2, eps, h1, hmin, maxnsteps);
    //     print_vec("ystart", ystart, nvar);
    // }

    yscal=dvector(1,nvar);
    y=dvector(1,nvar);
    dydx=dvector(1,nvar);
    x=x1;
    h=SIGN(h1,x2-x1);
    *nok = (*nbad) = kount = 0;
    for (i=1;i<=nvar;i++) y[i]=ystart[i];
    if (kmax > 0) xsav=x-dxsav*2.0;
    for (nstp=1;nstp<=maxnsteps;nstp++) {
        (*derivsin)(x,y,dydx);
        for (i=1;i<=nvar;i++)
            yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;

        // if (diffeqs_dbg_enabled()) {
        //     DDBG("STEP %d: x=% .15e h=% .3e (toward % .15e)", nstp, x, h, x2);
        //     stats_vec("yscal", yscal, nvar);
        //     stats_vec("dydx ", dydx , nvar);
        //     print_vec("y    ", y, nvar);
        // }

        if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
            xp[++kount]=x;
            for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
            xsav=x;
        }
        if ((x+h-x2)*(x+h-x1) > 0.0) {
            // if (diffeqs_dbg_enabled()) DDBG("trim last step: x=%.15e h=%.3e -> %.3e", x, h, x2-x);
            h=x2-x;
        }

        (*rkqsin)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivsin);

        if (hdid == h) ++(*nok); else ++(*nbad);

        // if (diffeqs_dbg_enabled()) {
        //     DDBG("AFTER rkqsin: x=% .15e hdid=% .3e hnext=% .3e nok=%d nbad=%d",
        //          x, hdid, hnext, *nok, *nbad);
        //     print_vec("y    ", y, nvar);
        // }

        if ((x-x2)*(x2-x1) >= 0.0) {
            for (i=1;i<=nvar;i++) ystart[i]=y[i];
            if (kmax) {
                xp[++kount]=x;
                for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
            }
            // if (diffeqs_dbg_enabled()) {
            //     DDBG("DONE at x=%.15e, total steps=%d, nok=%d, nbad=%d", x, nstp, *nok, *nbad);
            // }
            free_dvector(dydx,1,nvar);
            free_dvector(y,1,nvar);
            free_dvector(yscal,1,nvar);
            return;
        }
        if (fabs(hnext) <= hmin) {
            // DDBG("ERROR: hnext=%.3e <= hmin=%.3e at x=%.15e", fabs(hnext), hmin, x);
            nrerror("Step size too small in odeint");
        }
        h=hnext;
    }
    // DDBG("ERROR: exceeded maxnsteps=%d", maxnsteps);
    nrerror("Too many steps in routine odeint");
}
#undef MAXSTP
#undef TINY


/* ------------------------------ rkqs --------------------------------- */

#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

/* Restored to the original (BOSS) numerical logic */
void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
          double yscal[], double *hdid, double *hnext,
          void (*derivs)(double, double [], double []))
{
    void rkck(double y[], double dydx[], int n, double x, double h,
              double yout[], double yerr[], void (*derivs)(double, double [], double []));
    int i;
    double errmax,h,htemp,xnew,*yerr,*ytemp;

    yerr=dvector(1,n);
    ytemp=dvector(1,n);
    h=htry;

    // if (diffeqs_dbg_enabled()) {
    //     DDBG("ENTER x=% .15e htry=% .3e eps=%.3e n=%d", *x, htry, eps, n);
    //     print_vec("y    ", y, n);
    //     print_vec("dydx ", dydx, n);
    //     stats_vec("yscal", yscal, n);
    // }

    for (;;) {
        rkck(y,dydx,n,*x,h,ytemp,yerr,derivs);
        errmax=0.0;
        for (i=1;i<=n;i++) errmax=DMAX(errmax,fabs(yerr[i]/yscal[i]));
        errmax /= eps;
        if (errmax <= 1.0) break;

        // if (diffeqs_dbg_enabled()) {
        //     DDBG("REJECT: x=% .15e h=% .3e errmax/eps=%.6g", *x, h, errmax);
        // }

        htemp=SAFETY*h*rpow(errmax,PSHRNK);
        h=(h >= 0.0 ? DMAX(htemp,0.1*h) : FMIN(htemp,0.1*h));
        xnew=(*x)+h;
        if (xnew == *x) nrerror("stepsize underflow in rkqs");
    }
    if (errmax > ERRCON) *hnext=SAFETY*h*rpow(errmax,PGROW);
    else *hnext=5.0*h;
    *x += (*hdid=h);
    for (i=1;i<=n;i++) y[i]=ytemp[i];

    // if (diffeqs_dbg_enabled()) {
    //     DDBG("ACCEPT: x=% .15e hdid=% .3e hnext=% .3e errmax/eps=%.3e", *x, *hdid, *hnext, errmax);
    //     print_vec("y    ", y, n);
    // }

    free_dvector(ytemp,1,n);
    free_dvector(yerr,1,n);
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON


/* ------------------------------ rkck --------------------------------- */

void rkck(double y[], double dydx[], int n, double x, double h, double yout[],
          double yerr[], void (*derivs)(double, double [], double []))
{
    int i;
    static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
    b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
    b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
    b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
    b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
    c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
    dc5 = -277.00/14336.0;
    double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
    dc4=c4-13525.0/55296.0,dc6=c6-0.25;
    double *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;

    // static long calls = 0;
    // long every = rkck_every();
    // int do_print = diffeqs_dbg_enabled() && (every > 0) && (((++calls) % every) == 0);
    // if (do_print) {
    //     DDBG("ENTER x=% .15e h=% .3e n=%d (sampled call %ld)", x, h, n, calls);
    //     print_vec("y    ", y, n);
    //     print_vec("dydx ", dydx, n);
    // }

    ak2=dvector(1,n);
    ak3=dvector(1,n);
    ak4=dvector(1,n);
    ak5=dvector(1,n);
    ak6=dvector(1,n);
    ytemp=dvector(1,n);
    for (i=1;i<=n;i++)
        ytemp[i]=y[i]+b21*h*dydx[i];
    (*derivs)(x+a2*h,ytemp,ak2);
    for (i=1;i<=n;i++)
        ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
    (*derivs)(x+a3*h,ytemp,ak3);
    for (i=1;i<=n;i++)
        ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
    (*derivs)(x+a4*h,ytemp,ak4);
    for (i=1;i<=n;i++)
        ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
    (*derivs)(x+a5*h,ytemp,ak5);
    for (i=1;i<=n;i++)
        ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
    (*derivs)(x+a6*h,ytemp,ak6);
    for (i=1;i<=n;i++)
        yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
    for (i=1;i<=n;i++)
        yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);

    // if (do_print) {
    //     stats_vec("yout", yout, n);
    //     stats_vec("yerr", yerr, n);
    //     DDBG("EXIT");
    // }

    free_dvector(ytemp,1,n);
    free_dvector(ak6,1,n);
    free_dvector(ak5,1,n);
    free_dvector(ak4,1,n);
    free_dvector(ak3,1,n);
    free_dvector(ak2,1,n);
}


/* ------------------------------ bsstep ------------------------------- */

#define KMAXX 8
#define IMAXX (KMAXX+1)
#define SAFE1 0.25
#define SAFE2 0.7
#define REDMAX 1.0e-5
#define REDMIN 0.7
#define TINY 1.0e-30
#define SCALMX 0.1

double **d,*x;

void bsstep(double y[], double dydx[], int nv, double *xx, double htry, double eps,
            double yscal[], double *hdid, double *hnext,
            void (*derivs)(double, double [], double []))
{
    void mmid(double y[], double dydx[], int nvar, double xs, double htot,
              int nstep, double yout[], void (*derivs)(double, double[], double[]));
    void pzextr(int iest, double xest, double yest[], double yz[], double dy[],
                int nv);
    int i,iq,k,kk,km;
    static int first=1,kmax,kopt;
    static double epsold = -1.0,xnew;
    double eps1,errmax,fact,h,red,scale,work,wrkmin,xest;
    double *err,*yerr,*ysav,*yseq;
    static double a[IMAXX+1];
    static double alf[KMAXX+1][KMAXX+1];
    static int nseq[IMAXX+1]={0,2,4,6,8,10,12,14,16,18};
    int reduct,exitflag=0;

    // if (diffeqs_dbg_enabled()) {
    //     DDBG("ENTER Bulirsch–Stoer: nv=%d x=% .15e htry=%.3e eps=%.3e", nv, *xx, htry, eps);
    //     stats_vec("yscal", yscal, nv);
    // }

    d=dmatrix(1,nv,1,KMAXX);
    err=dvector(1,KMAXX);
    x=dvector(1,KMAXX);
    yerr=dvector(1,nv);
    ysav=dvector(1,nv);
    yseq=dvector(1,nv);
    if (eps != epsold) {
        *hnext = xnew = -1.0e29;
        eps1=SAFE1*eps;
        a[1]=nseq[1]+1;
        for (k=1;k<=KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
        for (iq=2;iq<=KMAXX;iq++) {
            for (k=1;k<iq;k++)
                alf[k][iq]=rpow(eps1,(a[k+1]-a[iq+1])/
                               ((a[iq+1]-a[1]+1.0)*(2*k+1)));
        }
        epsold=eps;
        for (kopt=2;kopt<KMAXX;kopt++)
            if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
        kmax=kopt;
    }
    h=htry;
    for (i=1;i<=nv;i++) ysav[i]=y[i];
    if (*xx != xnew || h != (*hnext)) {
        first=1;
        kopt=kmax;
    }
    reduct=0;
    for (;;) {
        for (k=1;k<=kmax;k++) {
            xnew=(*xx)+h;
            if (xnew == (*xx)) nrerror("step size underflow in bsstep");
            mmid(ysav,dydx,nv,*xx,h,nseq[k],yseq,derivs);
            xest=DSQR(h/nseq[k]);
            pzextr(k,xest,yseq,y,yerr,nv);
            if (k != 1) {
                errmax=TINY;
                for (i=1;i<=nv;i++) errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
                errmax /= eps;
                km=k-1;
                err[km]=rpow(errmax/SAFE1,1.0/(2*km+1));
            }
            // if (diffeqs_dbg_enabled()) {
            //     DDBG("k=%d errmax/eps=%.3e h=%.3e", k, (k!=1)?(errmax):0.0, h);
            // }
            if (k != 1 && (k >= kopt-1 || first)) {
                if (errmax < 1.0) {
                    exitflag=1;
                    break;
                }
                if (k == kmax || k == kopt+1) {
                    red=SAFE2/err[km];
                    break;
                }
                else if (k == kopt && alf[kopt-1][kopt] < err[km]) {
                    red=1.0/err[km];
                    break;
                }
                else if (kopt == kmax && alf[km][kmax-1] < err[km]) {
                    red=alf[km][kmax-1]*SAFE2/err[km];
                    break;
                }
                else if (alf[km][kopt] < err[km]) {
                    red=alf[km][kopt-1]/err[km];
                    break;
                }
            }
        }
        if (exitflag) break;
        red=FMIN(red,REDMIN);
        red=FMAX(red,REDMAX);
        // if (diffeqs_dbg_enabled()) DDBG("shrink h *= %.3e", red);
        h *= red;
        reduct=1;
    }
    *xx=xnew;
    *hdid=h;
    first=0;
    wrkmin=1.0e35;
    for (kk=1;kk<=km;kk++) {
        fact=FMAX(err[kk],SCALMX);
        work=fact*a[kk+1];
        if (work < wrkmin) {
            scale=fact;
            wrkmin=work;
            kopt=kk+1;
        }
    }
    *hnext=h/scale;
    if (kopt >= k && kopt != kmax && !reduct) {
        fact=FMAX(scale/alf[kopt-1][kopt],SCALMX);
        if (a[kopt+1]*fact <= wrkmin) {
            *hnext=h/fact;
            kopt++;
        }
    }
    // if (diffeqs_dbg_enabled()) {
    //     DDBG("EXIT: x=%.15e hdid=%.3e hnext=%.3e", *xx, *hdid, *hnext);
    //     print_vec("y    ", y, nv);
    // }
    free_dvector(yseq,1,nv);
    free_dvector(ysav,1,nv);
    free_dvector(yerr,1,nv);
    free_dvector(x,1,KMAXX);
    free_dvector(err,1,KMAXX);
    free_dmatrix(d,1,nv,1,KMAXX);
}
#undef KMAXX
#undef IMAXX
#undef SAFE1
#undef SAFE2
#undef REDMAX
#undef REDMIN
#undef TINY
#undef SCALMX


/* ------------------------------- mmid -------------------------------- */

void mmid(double y[], double dydx[], int nvar, double xs, double htot, int nstep,
          double yout[], void (*derivs)(double, double[], double[]))
{
    int n,i;
    double x,swap,h2,h,*ym,*yn;

    // if (diffeqs_dbg_enabled()) {
    //     DDBG("ENTER xs=% .15e htot=%.3e nstep=%d", xs, htot, nstep);
    //     print_vec("y    ", y, nvar);
    //     print_vec("dydx ", dydx, nvar);
    // }

    ym=dvector(1,nvar);
    yn=dvector(1,nvar);
    h=htot/nstep;
    for (i=1;i<=nvar;i++) {
        ym[i]=y[i];
        yn[i]=y[i]+h*dydx[i];
    }
    x=xs+h;
    (*derivs)(x,yn,yout);
    h2=2.0*h;
    for (n=2;n<=nstep;n++) {
        for (i=1;i<=nvar;i++) {
            swap=ym[i]+h2*yout[i];
            ym[i]=yn[i];
            yn[i]=swap;
        }
        x += h;
        (*derivs)(x,yn,yout);
    }
    for (i=1;i<=nvar;i++)
        yout[i]=0.5*(ym[i]+yn[i]+h*yout[i]);

    // if (diffeqs_dbg_enabled()) {
    //     print_vec("yout ", yout, nvar);
    //     DDBG("EXIT x=%.15e", x);
    // }
    free_dvector(yn,1,nvar);
    free_dvector(ym,1,nvar);
}


/* ------------------------------- pzextr ------------------------------ */

void pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv)
{
    int k1,j;
    double q,f2,f1,delta,*c;

    // if (diffeqs_dbg_enabled()) {
    //     DDBG("ENTER iest=%d xest=%.3e", iest, xest);
    //     /* small vector, printing can be helpful here */
    //     print_vec("yest ", yest, nv);
    // }

    c=dvector(1,nv);
    x[iest]=xest;
    for (j=1;j<=nv;j++) dy[j]=yz[j]=yest[j];
    if (iest == 1) {
        for (j=1;j<=nv;j++) d[j][1]=yest[j];
    } else {
        for (j=1;j<=nv;j++) c[j]=yest[j];
        for (k1=1;k1<iest;k1++) {
            delta=1.0/(x[iest-k1]-xest);
            f1=xest*delta;
            f2=x[iest-k1]*delta;
            for (j=1;j<=nv;j++) {
                q=d[j][k1];
                d[j][k1]=dy[j];
                delta=c[j]-q;
                dy[j]=f1*delta;
                c[j]=f2*delta;
                yz[j] += dy[j];
            }
        }
        for (j=1;j<=nv;j++) d[j][iest]=dy[j];
    }

    // if (diffeqs_dbg_enabled()) {
    //     print_vec("yz   ", yz, nv);
    //     print_vec("dy   ", dy, nv);
    //     DDBG("EXIT");
    // }
    free_dvector(c,1,nv);
}