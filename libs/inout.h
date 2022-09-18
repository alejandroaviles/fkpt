/*==============================================================================
	HEADER: inout.h				[General_libs]
==============================================================================*/

#ifndef _inout_h
#define _inout_h

#include "vectdefs.h"

// ------------[	inout normal definitions	 	]------------

void in_int(stream, int *);
void in_short(stream, short *);
void in_real(stream, real *);
void in_vector(stream, vector);
void out_int(stream, int);
void out_short(stream, short);
void out_real(stream, real);
void out_vector(stream, vector);


// ------------[	inout mar definitions	 		]------------

void out_int_mar(stream, int);
void out_short_mar(stream, short);
void out_bool_mar(stream, bool);
void out_real_mar(stream, real);          
void out_vector_mar(stream, vector);      


// ------------[	inout binary definitions	 	]------------

void in_int_bin(stream, int *);               
void in_short_bin(stream, short *);
void in_real_bin(stream, real *);             
void in_vector_bin(stream, vector);           
void out_int_bin(stream, int);
void out_short_bin(stream, short);
void out_real_bin(stream, real);
void out_vector_bin(stream, vector);
void out_bool_mar_bin(stream, bool);


// ------------[	inout other definitions	 		]------------

void in_vector_ruben(stream, vector);           

void out_vector_ndim(stream, double *, int);
void in_vector_ndim(stream, real *, int);

// Macros for binary in/out

#define safewrite(ptr,len,str)                  \
    if (fwrite((void *) ptr, len, 1, str) != 1) \
        error("safewrite: fwrite failed\n")

#define saferead(ptr,len,str)                  \
    if (fread((void *) ptr, len, 1, str) != 1) \
        error("saferead: fread failed\n")


// Macros and routines for binary gdgt in/out


size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream);
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream);

size_t gdgt_fread(void *ptr, size_t size, size_t nmemb, FILE *stream);
size_t gdgt_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream);

#define SKIP gdgt_fread(&blklen,sizeof(int4byte),1,fd);
#define BLKLEN gdgt_fwrite(&blklen, sizeof(blklen), 1, fd);

// IN/OUT ROUTINES TO HANDLE STRINGS ...

void ReadInString(stream, char *);
void ReadInLineString(stream, char *);


// InputData (como el que está en nplot2d)
// PARA IMPLEMENTAR LECTURA GENERAL DE ARCHIVOS DE DATOS CON FORMATO DE COLUMNAS

void inout_InputData(string, int, int, int *);
void inout_InputData_1c(string, int, int *);
void inout_InputData_3c(string filename, int col1, int col2, int col3,
                  int *npts);
void inout_InputData_4c(string filename, int col1, int col2, int col3, int col4,
                        int *npts);
//real *inout_xval;
//real *inout_yval;
//real *inout_zval;
//real *inout_wval;

// END: PARA IMPLEMENTAR LECTURA GENERAL DE ARCHIVOS DE DATOS CON FORMATO DE COLUMNAS

#endif	// ! _inout_h
