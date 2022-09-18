
#ifndef _mathutil_h
#define _mathutil_h

#include "stdinc.h"
#include "numrec.h"


void Rotation3D(real *, real, real, real);
real radian(real);
real degree(real);
real angle(real, real, real, real);

/* 
dmax es usado como variable en datanaly.c. 
Usar el formato de funcion y subroutina o pasar 
al formato de Macro MAX, MIN como ya definido en
stdinc.h. Cual sea mas eficiente o conveniente.
*/


double dmax(double,double);								
double dmin(double,double);								
int    imax(int,int);									
int    imin(int,int);

void int_piksrt(int n, int arr[]);	// Ordena un arreglo de enteros
									// de menor a mayor.
void HeapSort (real *, int *, int); // Ordena un arreglo de reales (Ascending)
void HeapSortInt (int *, int *, int); // Ordena un arreglo de enteros (Ascending)
void HeapSortDescend (real *, int *, int); // Ordena un arreglo de reales (Descending)
										// Esta rutina no es necesaria ya que con
										// un loop descendiente podemos operar.
										// Pero se deja por completitud.
void HeapSortIntDescend (int *, int *, int); // Ordena un arreglo de enteros (Descending)
										// Esta rutina no es necesaria ya que con
										// un loop descendiente podemos operar.
										// Pero se deja por completitud.

void FftComplex (Cmplx *, int);		// Fast Fourier Transform of complex array data
real Integrate (real *, int);

void sort3arrays(unsigned long n, double arr[], double brr[], double crr[]);


// Fitting
void fit(double x[], double y[], int ndata, double sig[], int mwt,
         double *a, double *b, double *siga, double *sigb, double *chi2, double *q);

#endif  


