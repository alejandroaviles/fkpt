
#ifndef _stdinc_h
#define _stdinc_h

#include <stdio.h>
#include <stdlib.h>

#define copyright	"Copyright (c) .... "

//long idum;				// seed for random generators

#if !defined(NULL)
#define NULL 0L
#endif

#define local     static

typedef short int bool;

#if !defined(TRUE)
#define TRUE  ((bool) 1)
#define FALSE ((bool) 0)
#endif

typedef unsigned char byte;

typedef char *string;

typedef FILE *stream;							
												
#include "precision.h"

#if !defined(MIXEDPREC) && !defined(SINGLEPREC) && !defined(DOUBLEPREC)
#define SINGLEPREC
#endif

#if defined(DOUBLEPREC)
#undef SINGLEPREC
#undef MIXEDPREC
typedef double real, *realptr;
#define Precision "DOUBLEPREC"
#endif

#if defined(MIXEDPREC)
#undef DOUBLEPREC
#undef SINGLEPREC
typedef float *realptr, real;
#define Precision "MIXEDPREC"
#endif

#if defined(SINGLEPREC)
#undef DOUBLEPREC
#undef MIXEDPREC
typedef float real, *realptr;
#define Precision "SINGLEPREC"
#endif

// Complex definitions

typedef struct {
  real R, I;
} Cmplx;

#define CSet(a, x, y)                                       \
   a.R = x,                                                 \
   a.I = y
#define CAdd(a, b, c)                                       \
   a.R = b.R + c.R,                                         \
   a.I = b.I + c.I
#define CSub(a, b, c)                                       \
   a.R = b.R - c.R,                                         \
   a.I = b.I - c.I
#define CMul(a, b, c)                                       \
  a.R = b.R * c.R - b.I * c.I,                              \
  a.I = b.R * c.I + b.I * c.R

#ifndef PI
#define PI		   3.141592653589793238462643383279502884197
#endif

#define TWO_PI     6.28318530717958647693
#define FOUR_PI   12.56637061435917295385
#define HALF_PI    1.57079632679489661923
#define FRTHRD_PI  4.18879020478639098462

#if !defined(M_LN2)
#define M_LN2		0.69314718055994530942
#endif
#if !defined(M_LN10)
#define M_LN10		2.30258509299404568402
#endif


#if !defined (HZ)
#define HZ        100
#endif


#define streq(x,y) (strcmp((x), (y)) == 0)
#define strnull(x) (strcmp((x), "") == 0)


#define ABS(x)   (((x)<0)?-(x):(x))

#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif
#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif

#define SignR(x,y)  (((y) >= 0) ? (x) : (- (x)))


//void *allocate(int);			// Definicion original ...
void *allocate(long int);		// Correccion para trabajar con mas de 16x10^6 particulas...			
real *AllocVecR(int);			// Define un arreglo de reales como en fortran
int *AllocVecI(int);			// Define un arreglo de enteros como en fortran
int *AllocVecINormal(int);		// Define un arreglo de enteros como en C (comenzando en cero)
void FreeVecR(real *);
void FreeVecI(int *);
void FreeVecINormal(int *);


#define AllocMem(a, n, t)  a = (t *) malloc ((n) * sizeof (t))
#define AllocMem2(a, n1, n2, t)                             \
   AllocMem (a, n1, t *);                                   \
   AllocMem (a[0], (n1) * (n2), t);                         \
   for (k = 1; k < n1; k ++) a[k] = a[k - 1] + n2;

#define Cube(x)    ((x) * (x) * (x))

#define Nint(x)                                             \
   (((x) < 0.) ? (- (int) (0.5 - (x))): ((int) (0.5 + (x))))

double cputime(void);

void error(string, ...);
void   endrun(int);
												
void eprintf(string, ...);						

bool scanopt(string, string);					

stream stropen(string, string);					

double second(void);
double timediff(double t0,double t1);

												
#endif  /* ! _stdinc_h	*/

