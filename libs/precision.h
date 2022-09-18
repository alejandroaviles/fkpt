

#define DOUBLEPREC
/* #undef DOUBLEPREC */

#define SINGLEPREC
#undef SINGLEPREC

#define MIXEDPREC
#undef MIXEDPREC

#if defined(MIXEDPREC)
#undef MIXEDPREC
#undef SINGLEPREC
#endif

#if defined(SINGLEPREC)
#undef SINGLEPREC
#undef MIXEDPREC
#endif

