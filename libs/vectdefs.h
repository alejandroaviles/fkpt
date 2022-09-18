
//	It is very important to clean the object and library 
//	files, and remake all from scratch, when switching dimensions. 
//	Not doing so cause memory problems. 


#ifndef _vectdefs_h
#define _vectdefs_h

//#include "../general/stdinc.h"

#if !defined(NDIM) && !defined(TWODIM) && !defined(THREEDIM) && !defined(ONEDIM)
#define THREEDIM
//#define TWODIM
//#define ONEDIM
#endif

#if defined(THREEDIM) || (NDIM==3)
#undef  ONEDIM
#undef  TWODIM
#define THREEDIM
#define NDIM 3
#endif

#if defined(TWODIM) || (NDIM==2)
#undef  ONEDIM
#undef  THREEDIM
#define TWODIM
#define NDIM 2
#endif

#if defined(ONEDIM) || (NDIM==1)
#undef  THREEDIM
#undef  TWODIM
#define ONEDIM
#define NDIM 1
#endif

typedef real vector[NDIM];
typedef real matrix[NDIM][NDIM];

typedef int vectorI[NDIM];

// To be compatible with R...
typedef struct {real x, y;} VecR2;
typedef struct {real x, y, z;} VecR3;
typedef struct {int x, y;} VecI2;
typedef struct {int x, y, z;} VecI3;

#if NDIM == 2
typedef VecR2 VecR;
typedef VecI2 VecI;
#endif

#if NDIM == 3
typedef VecR3 VecR;
typedef VecI3 VecI;
#endif

#endif  // ! _vectdefs_h

