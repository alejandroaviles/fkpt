
#ifndef _mathfns_h
#define _mathfns_h

#include <math.h>
#include "stdinc.h"
#include "vectdefs.h"


#if defined(MIXEDPREC) || defined(DOUBLEPREC)
#define rsqrt    sqrt
#define rcbrt    cbrt
#define rsin     sin
#define rcos     cos
#define rtan     tan
#define rasin    asin
#define racos    acos
#define ratan    atan
#define ratan2   atan2
#define rlog     log
#define rexp     exp
#define rlog10   log10
#define rsinh    sinh
#define rcosh    cosh
#define rtanh    tanh
#define rpow     pow
#define rabs     fabs
#define rfloor   floor
#define rceil    ceil
#endif

#if defined(SINGLEPREC)
#define rsqrt    fsqrt
#define rsin     fsin
#define rcos     fcos
#define rtan     ftan
#define rasin    fasin
#define racos    facos
#define ratan    fatan
#define ratan2   fatan2
#define rlog     flog
#define rexp     fexp
#define rlog10   log10f
#define rsinh    fsinh
#define rcosh    fcosh
#define rtanh    ftanh
#define rpow     powf
#define rabs     fabsf
#define rfloor   ffloor
#define rceil    fceil
#endif


#if defined(MIXEDPREC) || defined(DOUBLEPREC)
#define rsqr     sqr
#define rqbe     qbe
#define rlog2    log2
#define rexp2    exp2
#define rdex     dex
#endif

#if defined(SINGLEPREC)
#define rsqr     fsqr
#define rqbe     fqbe
#define rlog2    flog2
#define rexp2    fexp2
#define rdex     fdex
#define rcbrt    fcbrt
#endif

real rsqr(real);
real rqbe(real);
/* 
Las siguientes dos funciones causan el warning de definici'on
multiple. Parece que ya estan incluidas en las librerias del
sistema.
*/

/* No estan incluidas en VEGA	*/

/*
real rlog2(real);
real rexp2(real);
*/
real rdex(real);

#if defined(SINGLEPREC)
float rcbrt(float);
#endif

void xsrandom(long);
double xrandom(double, double);
double grandom(double, double);


#if defined(MIXEDPREC)
#define pickshell mpickshell
#define pickball  mpickball
#define pickbox   mpickbox
#endif

#if defined(DOUBLEPREC)
#define pickshell dpickshell
#define pickball  dpickball
#define pickbox   dpickbox
#endif

#if defined(SINGLEPREC)
#define pickshell fpickshell
#define pickball  fpickball
#define pickbox   fpickbox
#endif


void pickshell(real *, int, real);
void pickball(real *, int, real);
void pickbox(real *, int, real);

/*
real Heaveside(real);
*/

int nint(real);

void bessel(real, realptr, realptr, realptr, realptr);
real erfcc(real);
real expint(real);

// Bessel functions from Numerical Recipes ...
double bessi0(double);
double bessi1(double);
double bessj0(double);
double bessj1(double);
double bessj(int, double);
double bessk0(double);
double bessk1(double);
// End Bessel functions from Numerical Recipes


double gammln(double xx);  // Ln(Gamma(x))
double beta(double z, double w);
double erffnr(double x);            // Error function NR 2nd E
double erffcnr(double x);           // Complementary Error function NR 2nd E
//double erfcc(double x);               // Complementary Error function NR 2nd E
                                    // using Chebyshev
double gammp(double a, double x);
double gammq(double a, double x);
void gser(double *gamser, double a, double x, double *gln);
void gcf(double *gammcf, double a, double x, double *gln);

void VRand(vector p);

real rj0Bessel(real);
real rj1Bessel(real);
real rj2Bessel(real);
real rj3Bessel(real);
real rj4Bessel(real);

#endif	/* ! _mathfns_h	*/

