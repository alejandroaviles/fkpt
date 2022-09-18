
#ifndef _vectmath_h
#define _vectmath_h

#include "vectdefs.h"

#ifndef NDIM
#define NDIM  3
#endif


#define CLRV(v)									                        \
{                                                                       \
    int _i;                                                             \
    for (_i = 0; _i < NDIM; _i++)                                       \
        (v)[_i] = 0.0;                                                  \
}

#define UNITV(v,j)                                                      \
{                                                                       \
    int _i;                                                             \
    for (_i = 0; _i < NDIM; _i++)                                       \
        (v)[_i] = (_i == (j) ? 1.0 : 0.0);                              \
}

#define SETV(v,u)                                                       \
{                                                                       \
    int _i;                                                             \
    for (_i = 0; _i < NDIM; _i++)                                       \
        (v)[_i] = (u)[_i];                                              \
}

#if defined(THREEDIM)

#define ADDV(v,u,w)                                                     \
{                                                                       \
    (v)[0] = (u)[0] + (w)[0];                                           \
    (v)[1] = (u)[1] + (w)[1];                                           \
    (v)[2] = (u)[2] + (w)[2];                                           \
}

#define ADD2VS(v,u,w,s)                                                 \
{                                                                       \
    (v)[0] = (u)[0] + (w)[0] * (s);                                     \
    (v)[1] = (u)[1] + (w)[1] * (s);                                     \
    (v)[2] = (u)[2] + (w)[2] * (s);                                     \
}

#define SUBV(v,u,w)                                                     \
{                                                                       \
    (v)[0] = (u)[0] - (w)[0];                                           \
    (v)[1] = (u)[1] - (w)[1];                                           \
    (v)[2] = (u)[2] - (w)[2];                                           \
}

#define MULVS(v,u,s)                                                    \
{                                                                       \
    (v)[0] = (u)[0] * s;                                                \
    (v)[1] = (u)[1] * s;                                                \
    (v)[2] = (u)[2] * s;                                                \
}

#else

#define ADDV(v,u,w)                                                     \
{                                                                       \
    int _i;                                                             \
    for (_i = 0; _i < NDIM; _i++)                                       \
        (v)[_i] = (u)[_i] + (w)[_i];                                    \
}

#define ADD2VS(v,u,w,s)                                                 \
{                                                                       \
    int _i;                                                             \
    for (_i = 0; _i < NDIM; _i++)                                       \
        (v)[_i] = (u)[_i] + (w)[_i] * (s);                              \
}

#define SUBV(v,u,w)                                                     \
{                                                                       \
    int _i;                                                             \
    for (_i = 0; _i < NDIM; _i++)                                       \
        (v)[_i] = (u)[_i] - (w)[_i];                                    \
}

#define MULVS(v,u,s)                                                    \
{                                                                       \
    int _i;                                                             \
    for (_i = 0; _i < NDIM; _i++)                                       \
        (v)[_i] = (u)[_i] * (s);                                        \
}

#endif

#define DIVVS(v,u,s)                                                    \
{                                                                       \
    int _i;                                                             \
    for (_i = 0; _i < NDIM; _i++)                                       \
        (v)[_i] = (u)[_i] / (s);                                        \
}

#if defined(THREEDIM)

#define DOTVP(s,v,u)                                                    \
{                                                                       \
    (s) = (v)[0]*(u)[0] + (v)[1]*(u)[1] + (v)[2]*(u)[2];                \
}

#else

#define DOTVP(s,v,u)                                                    \
{                                                                       \
    int _i;                                                             \
    (s) = 0.0;                                                          \
    for (_i = 0; _i < NDIM; _i++)                                       \
        (s) += (v)[_i] * (u)[_i];                                       \
}

#endif

#define ABSV(s,v)                                                       \
{                                                                       \
    real _tmp;                                                          \
    int _i;                                                             \
    _tmp = 0.0;                                                         \
    for (_i = 0; _i < NDIM; _i++)                                       \
        _tmp += (v)[_i] * (v)[_i];                                      \
    (s) = rsqrt(_tmp);                                                  \
}

#define DISTV(s,u,v)                                                    \
{                                                                       \
    real _tmp;                                                          \
    int _i;                                                             \
    _tmp = 0.0;                                                         \
    for (_i = 0; _i < NDIM; _i++)                                       \
        _tmp += ((u)[_i]-(v)[_i]) * ((u)[_i]-(v)[_i]);                  \
    (s) = rsqrt(_tmp);                                                  \
}

#if defined(ONEDIM)

#define CROSSVP(s,v,u)                                                  \
{                                                                       \
    (s) = 0.0;															\
/*    (s)[0] = 0.0;	*/														\
}

#endif

#if defined(TWODIM)

#define CROSSVP(s,v,u)                                                  \
{																		\
    (s) = (v)[0]*(u)[1] - (v)[1]*(u)[0];                                \
}

#endif

#if defined(THREEDIM)

#define CROSSVP(v,u,w)                                                  \
{                                                                       \
    (v)[0] = (u)[1]*(w)[2] - (u)[2]*(w)[1];                             \
    (v)[1] = (u)[2]*(w)[0] - (u)[0]*(w)[2];                             \
    (v)[2] = (u)[0]*(w)[1] - (u)[1]*(w)[0];                             \
}

#endif


#define CLRM(p)                                                         \
{                                                                       \
    int _i, _j;                                                         \
    for (_i = 0; _i < NDIM; _i++)                                       \
        for (_j = 0; _j < NDIM; _j++)                                   \
            (p)[_i][_j] = 0.0;                                          \
}

#define SETMI(p)                                                        \
{                                                                       \
    int _i, _j;                                                         \
    for (_i = 0; _i < NDIM; _i++)                                       \
        for (_j = 0; _j < NDIM; _j++)                                   \
            (p)[_i][_j] = (_i == _j ? 1.0 : 0.0);                       \
}

#define SETM(p,q)                                                       \
{                                                                       \
    int _i, _j;                                                         \
    for (_i = 0; _i < NDIM; _i++)                                       \
        for (_j = 0; _j < NDIM; _j++)                                   \
            (p)[_i][_j] = (q)[_i][_j];                                  \
}

#define TRANM(p,q)                                                      \
{                                                                       \
    int _i, _j;                                                         \
    for (_i = 0; _i < NDIM; _i++)                                       \
        for (_j = 0; _j < NDIM; _j++)                                   \
            (p)[_i][_j] = (q)[_j][_i];                                  \
}

#define ADDM(p,q,r)                                                     \
{                                                                       \
    int _i, _j;                                                         \
    for (_i = 0; _i < NDIM; _i++)                                       \
        for (_j = 0; _j < NDIM; _j++)                                   \
            (p)[_i][_j] = (q)[_i][_j] + (r)[_i][_j];                    \
}

#define SUBM(p,q,r)                                                     \
{                                                                       \
    int _i, _j;                                                         \
    for (_i = 0; _i < NDIM; _i++)                                       \
        for (_j = 0; _j < NDIM; _j++)                                   \
            (p)[_i][_j] = (q)[_i][_j] - (r)[_i][_j];                    \
}

#define MULM(p,q,r)                                                     \
{                                                                       \
    int _i, _j, _k;                                                     \
    for (_i = 0; _i < NDIM; _i++)                                       \
        for (_j = 0; _j < NDIM; _j++) {                                 \
            (p)[_i][_j] = 0.0;                                          \
            for (_k = 0; _k < NDIM; _k++)                               \
                (p)[_i][_j] += (q)[_i][_k] * (r)[_k][_j];               \
        }                                                               \
}

#define MULMS(p,q,s)                                                    \
{                                                                       \
    int _i, _j;                                                         \
    for (_i = 0; _i < NDIM; _i++)                                       \
        for (_j = 0; _j < NDIM; _j++)                                   \
            (p)[_i][_j] = (q)[_i][_j] * (s);                            \
}

#define DIVMS(p,q,s)                                                    \
{                                                                       \
    int _i, _j;                                                         \
    for (_i = 0; _i < NDIM; _i++)                                       \
        for (_j = 0; _j < NDIM; _j++)                                   \
            (p)[_i][_j] = (q)[_i][_j] / (s);                            \
}

#define MULMV(v,p,u)                                                    \
{                                                                       \
    int _i, _j;                                                         \
    for (_i = 0; _i < NDIM; _i++) {                                     \
        (v)[_i] = 0.0;                                                  \
        for (_j = 0; _j < NDIM; _j++)                                   \
            (v)[_i] += (p)[_i][_j] * (u)[_j];                           \
    }                                                                   \
}

#define OUTVP(p,v,u)                                                    \
{                                                                       \
    int _i, _j;                                                         \
    for (_i = 0; _i < NDIM; _i++)                                       \
        for (_j = 0; _j < NDIM; _j++)                                   \
            (p)[_i][_j] = (v)[_i] * (u)[_j];                            \
}

#define TRACEM(s,p)                                                     \
{                                                                       \
    int _i;                                                             \
    (s) = 0.0;                                                          \
    for (_i = 0.0; _i < NDIM; _i++)                                     \
        (s) += (p)[_i][_i];                                             \
}


#if defined(THREEDIM)

#define DOTPSUBV(s,v,u,w)										        \
{                                                                       \
    (v)[0] = (u)[0] - (w)[0];    (s)  = (v)[0] * (v)[0];                \
    (v)[1] = (u)[1] - (w)[1];    (s) += (v)[1] * (v)[1];                \
    (v)[2] = (u)[2] - (w)[2];    (s) += (v)[2] * (v)[2];                \
}

#define DOTPMULMV(s,v,p,u)											    \
{                                                                       \
    DOTVP(v[0], p[0], u);    (s)  = (v)[0] * (u)[0];                    \
    DOTVP(v[1], p[1], u);    (s) += (v)[1] * (u)[1];                    \
    DOTVP(v[2], p[2], u);    (s) += (v)[2] * (u)[2];                    \
}

#define ADDMULVS(v,u,s)												    \
{                                                                       \
    (v)[0] += (u)[0] * (s);                                             \
    (v)[1] += (u)[1] * (s);                                             \
    (v)[2] += (u)[2] * (s);                                             \
}

#define ADDMULVS2(v,u,s,w,r)										    \
{                                                                       \
    (v)[0] += (u)[0] * (s) + (w)[0] * (r);                              \
    (v)[1] += (u)[1] * (s) + (w)[1] * (r);                              \
    (v)[2] += (u)[2] * (s) + (w)[2] * (r);                              \
}

#define ADDVMULVS(v,u,w,s)      /* MUL Vect by Scalar, ADD to Vects */ \
{                                                                      \
(v)[0] = (u)[0] + (w)[0] * (s);                                    \
(v)[1] = (u)[1] + (w)[1] * (s);                                    \
(v)[2] = (u)[2] + (w)[2] * (s);                                    \
}

#endif

#if defined(TWODIM)

#define DOTPSUBV(s,v,u,w)										        \
{                                                                       \
    (v)[0] = (u)[0] - (w)[0];    (s)  = (v)[0] * (v)[0];                \
    (v)[1] = (u)[1] - (w)[1];    (s) += (v)[1] * (v)[1];                \
}

#define DOTPMULMV(s,v,p,u)											    \
{                                                                       \
    DOTVP(v[0], p[0], u);    (s)  = (v)[0] * (u)[0];                    \
    DOTVP(v[1], p[1], u);    (s) += (v)[1] * (u)[1];                    \
}

#define ADDMULVS(v,u,s)										   		    \
{                                                                       \
    (v)[0] += (u)[0] * (s);                                             \
    (v)[1] += (u)[1] * (s);                                             \
}

#define ADDMULVS2(v,u,s,w,r)										    \
{                                                                       \
    (v)[0] += (u)[0] * (s) + (w)[0] * (r);                              \
    (v)[1] += (u)[1] * (s) + (w)[1] * (r);                              \
}

#define ADDVMULVS(v,u,w,s)      /* MUL Vect by Scalar, ADD to Vects */ \
{                                                                      \
(v)[0] = (u)[0] + (w)[0] * (s);                                    \
(v)[1] = (u)[1] + (w)[1] * (s);                                    \
}

#endif

#if defined(ONEDIM)

#define DOTPSUBV(s,v,u,w)										        \
{                                                                       \
    (v)[0] = (u)[0] - (w)[0];    (s)  = (v)[0] * (v)[0];                \
}

#define DOTPMULMV(s,v,p,u)											    \
{                                                                       \
    DOTVP(v[0], p[0], u);    (s)  = (v)[0] * (u)[0];                    \
}

#define ADDMULVS(v,u,s)												    \
{                                                                       \
    (v)[0] += (u)[0] * (s);                                             \
}

#define ADDMULVS2(v,u,s,w,r)										    \
{                                                                       \
    (v)[0] += (u)[0] * (s) + (w)[0] * (r);                              \
}

#define ADDVMULVS(v,u,w,s)      /* MUL Vect by Scalar, ADD to Vects */ \
{                                                                      \
(v)[0] = (u)[0] + (w)[0] * (s);                                    \
}

#endif

#define DOTPSUBV2(s,u,v,w,z)                                            \
{                                                                       \
    real _tmp;                                                          \
    int _i;                                                             \
    _tmp = 0.0;                                                         \
    for (_i = 0; _i < NDIM; _i++)                                       \
        _tmp += ((u)[_i]-(v)[_i]) * ((w)[_i]-(z)[_i]);                  \
    (s) = _tmp;															\
}
 

#define SETVS(v,s)                                                      \
{                                                                       \
    int _i;                                                             \
    for (_i = 0; _i < NDIM; _i++)                                       \
        (v)[_i] = (s);                                                  \
}

#define ADDVS(v,u,s)                                                    \
{                                                                       \
    int _i;                                                             \
    for (_i = 0; _i < NDIM; _i++)                                       \
        (v)[_i] = (u)[_i] + (s);                                        \
}

#define SETMS(p,s)                                                      \
{                                                                       \
    int _i, _j;                                                         \
    for (_i = 0; _i < NDIM; _i++)                                       \
        for (_j = 0; _j < NDIM; _j++)                                   \
            (p)[_i][_j] = (s);                                          \
}

#endif  

// Check equivalences above...

#define VDiv(v1, v2, v3)												\
{                                                                       \
    int _i;                                                             \
    for (_i = 0; _i < NDIM; _i++)                                       \
		(v1)[_i] = (v2)[_i] / (v3)[_i];									\
}

#if NDIM == 3

#define VProd(v)														\
   ((v)[0] * (v)[1] * (v)[2])
#define VLinear(p, s)													\
   (((p)[2] * (s)[1] + (p)[1]) * (s)[0] + (p)[0])
#define VSet(v, sx, sy, sz)												\
   (v)[0] = sx,															\
   (v)[1] = sy,															\
   (v)[2] = sz
#define VScale(v, s)													\
   (v)[0] *= s,															\
   (v)[1] *= s,															\
   (v)[2] *= s
#define VSCopy(v2, s1, v1)												\
   (v2)[0] = (s1) * (v1)[0],											\
   (v2)[1] = (s1) * (v1)[1],											\
   (v2)[2] = (s1) * (v1)[2]
#define VDot(v1, v2)													\
   ((v1)[0] * (v2)[0] + (v1)[1] * (v2)[1] + (v1)[2] * (v2)[2])
#define VAdd(v1, v2, v3)												\
   (v1)[0] = (v2)[0] + (v3)[0],											\
   (v1)[1] = (v2)[1] + (v3)[1],											\
   (v1)[2] = (v2)[2] + (v3)[2]
#define VSetAll(v, s)													\
   VSet (v, s, s, s)
#define VComp(v, k)														\
   *((k == 0) ? &(v)[0] : ((k == 1) ? &(v)[1] : &(v)[2]))
#define VSub(v1, v2, v3)												\
   (v1)[0] = (v2)[0] - (v3)[0],											\
   (v1)[1] = (v2)[1] - (v3)[1],											\
   (v1)[2] = (v2)[2] - (v3)[2]

#endif

#if NDIM == 2

#define VProd(v)														\
   ((v)[0] * (v)[1])
#define VLinear(p, s)													\
   ((p)[1] * (s)[0] + (p)[0])
#define VSet(v, sx, sy)													\
   (v)[0] = sx,															\
   (v)[1] = sy
#define VScale(v, s)													\
   (v)[0] *= s,															\
   (v)[1] *= s
#define VSCopy(v2, s1, v1)												\
   (v2)[0] = (s1) * (v1)[0],											\
   (v2)[1] = (s1) * (v1)[1]
#define VDot(v1, v2)													\
   ((v1)[0] * (v2)[0] + (v1)[1] * (v2)[1])
#define VAdd(v1, v2, v3)												\
   (v1)[0] = (v2)[0] + (v3)[0],											\
   (v1)[1] = (v2)[1] + (v3)[1]
#define VSetAll(v, s)													\
   VSet (v, s, s)
#define VComp(v, k)														\
   *((k == 0) ? &(v)[0] : &(v)[1])

#define VSub(v1, v2, v3)												\
   (v1)[0] = (v2)[0] - (v3)[0],											\
   (v1)[1] = (v2)[1] - (v3)[1]

#endif

#define VZero(v)  VSetAll (v, 0)
#define VVAdd(v1, v2)  VAdd (v1, v1, v2)

#define VMul(v1, v2, v3)												\
{                                                                       \
    int _i;                                                             \
    for (_i = 0; _i < NDIM; _i++)                                       \
		(v1)[_i] = (v2)[_i] * (v3)[_i];									\
}

#define VVSub(v1, v2)  SUBV(v1, v1, v2)
#define VLenSq(v)  VDot (v, v)
#define VVSAdd(v1, s2, v2) ADD2VS(v1, v1, v2, s2)
