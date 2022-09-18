
#ifndef _quad_h
#define _quad_h

#include "stdinc.h"


void gauleg(double x1, double x2, double x[], double w[], int n);
double qgauss(double (*func)(double), double a, double b,
              double x[], double w[], int n);

//double qromo(double (*func)(double), double a, double b,
//             double (*choose)(double (*)(double), double, double, int));
//double qromo(double (*func)(double), double a, double b,
//             double (*choose)(double (*)(double), double, double, int),
//             double epsq);
double qromo(double (*func)(double), double a, double b,
             double (*choose)(double (*)(double), double, double, int),
             double epsq, int KK);

void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
double midexp(double (*funk)(double), double aa, double bb, int n);
double midinf(double (*funk)(double), double aa, double bb, int n);
double midpnt(double (*func)(double), double a, double b, int n);
//double midpnt(double (*func)(double), double a, double b, int n);
double midsql(double (*funk)(double), double aa, double bb, int n);
double midsqu(double (*funk)(double), double aa, double bb, int n);

void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);

void locate(double xx[], unsigned long n, double x, unsigned long *j);

#endif // ! _quad_h


