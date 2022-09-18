/*==============================================================================
 HEADER: globaldefs.h		[gsm]
 ==============================================================================*/

#ifndef _globaldefs_h
#define _globaldefs_h


#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
//

#include "../libs/stdinc.h"
#include "../libs/numrec.h"
#include "../libs/diffeqs.h"
#include "../libs/quads.h"
#include "../libs/mathfns.h"
#include "../libs/mathutil.h"
#include "../libs/inout.h"
#include "../libs/vectmath.h"
#include "../libs/getparam.h"
#include "../libs/machines.h"
#include "../libs/strings.h"

#if !defined(global)
#  define global extern
#endif

#define IPName(param,paramtext)    \
{strcpy(tag[nt],paramtext);    \
addr[nt]=&(param);    \
id[nt++]=INT;}

#define RPName(param,paramtext)    \
{strcpy(tag[nt],paramtext);    \
addr[nt]=&param;    \
id[nt++]=DOUBLE;}

#define BPName(param,paramtext)    \
{strcpy(tag[nt],paramtext);    \
addr[nt]=&param;    \
id[nt++]=BOOLEAN;}

#define SPName(param,paramtext,n)    \
{strcpy(tag[nt],paramtext);    \
param=(string) malloc(n);    \
addr[nt]=param;    \
id[nt++]=STRING;}
//

#include "models.h"


#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#define invH0     2997.92458   //This is H_0^{-1} in Mpc/h units
#define FOURPI2   39.4784176043574    //but see 0903.5321
#define PI2     9.8696044010893586188
#define TWOPI2     19.739208802178716
#define SIXPI2  59.21762640653615
#define INVSQRTDTWOPI 0.39894228040143267794

typedef struct {
// Background cosmology:
    real om;
    string olstr;
    real h;
// bias parameters:
    real b1;
    real b2;
//    real bs;
    real bs2;
    real c1eft;
    real c2eft;
    real s2eft;
//
    real b3nl;
    real alpha0;
    real alpha2;
    real alpha4;
    real ctilde;
    real PshotP;
    real alpha0shot;
    real alpha2shot;
// Differential equations evolution parameters:
	real x;
	string dxstr;
	real xstop; 
    int maxnsteps;
	string integration_method;
    real dxmin;
    real eps;

// k table
    string fnamePS;
    real kmin;
    real kmax;
    int Nk;
//
// 2pcf rsd multipoles array:
    real smin;
    real smax;
    int Ns;
    real gsm_width;
    int gsm_sizeyT;
    int gsm_NGL;
// q functions:
	int NqperLogDecade;
	int Nk_qFunctionsQuad;
// CLPT correlation functions output array:
    real rmin;
    real rmax;
    int Nr;
// Post processing parameters:
    bool postprocessing;
    string options;
//
    string paramfile;
// Modified gravity model parameters:
    string mgmodel;
    string suffixModel;
    string model_paramfile;
    int nHS;
    real fR0;
    real omegaBD;
    real screening;
// DGP:
    real eps_DGP;
    real rc_DGP;
// Quadrature parameters:
    string quadratureMethod;
    int nquadSteps;
    int ngausslegpoints;
    real epsquad;
//
} cmdline_data, *cmdline_data_ptr;




typedef struct {
	real cpuinit;
	real dx;
    int method_int;
    int quadmethod_int;

// Modified gravity model parameters:
    real beta2;
//
    char integration_method_comment[100];
    char quadraturemethod_comment[100];

	string headline0;
	string headline1;
	string headline2;
	string headline3;

    char model_comment[100];

	FILE *outlog;
    
    real ol;

    real xnow;
    real xout;
    real xoutinfo;
    real xstop;

	char mode[2];
// 
    char fnamePSPath[100];
    char logfilePath[100];
    char clptDir[100];
    char tmpDir[100];
    char fpfnamekfun[100];
    char fpfnamekfunnw[100];
    char fpfnameSPTPowerSpectrum[100];
    char fpfnameqfunctions[100];
    char fpfnameclptfunctions[100];
    char fpfnamersd[100];
    char fpfnamev12[100];
    char fpfnamesigma12_parallel[100];
    char fpfnamesigma12_perp[100];
    //~ char fpfnameTables[100];
    char fpfnameParams[100];


    //~ char fpfnamekfun2[100];
    //~ char fpfnameclptfunctions2[100];    

    real kf;
    real k1;
    real k2;

    real x;
    real k;
    real p;
    
    real f0;
    real Dplus;
    real particles_meanPath;
    real sigma2L;
    real sigma8;
    real sigma2v;
    
    real Sigma2;
    real deltaSigma2;
    
} global_data, *global_data_ptr;


global global_data gd;
global cmdline_data cmd;

global real *yout;
#define NEQS3Order    10
#define NEQS2Order    8
#define NEQS1Order      2

typedef struct _pointPSTable {
    real k;
    real ps;
} pointPSTable, *pointPSTableptr;

global int nPSTable;
global pointPSTableptr PSLCDMtab;
global int nPSLogT;
global pointPSTableptr PSLCDMLogtab;

global int nPSLT;
global pointPSTableptr PSLT;
global int nPSLTLog;
global pointPSTableptr PSLTLog;

global real *kPS;
global real *pPS;
global real *pPS2;

global real *pPS_nw;
global real *pPS2_nw;




#define kPos(x)    (((pointPSTableptr) (x))->k)
#define PS(x)    (((pointPSTableptr) (x))->ps)


global real *fkT;
global real *fkT2;
global real *DplusT;
global real *DplusT2;




typedef struct {
    real eta;
    real y1;
    real y2;
    real y3;
    real y4;
    real y5;
    real y6;
    real y7;
    real y8;
} global_D2v2, *global_D2v2_ptr;

#define etaD2(x)    (((global_D2v2_ptr) (x))->eta)
#define Dpk1D2(x)    (((global_D2v2_ptr) (x))->y1)
#define Dpk2D2(x)    (((global_D2v2_ptr) (x))->y3)
#define DA2D2(x)    (((global_D2v2_ptr) (x))->y5)
#define DA2primeD2(x)    (((global_D2v2_ptr) (x))->y6)   //Modificacion-LCDMfk
#define DB2D2(x)    (((global_D2v2_ptr) (x))->y7)
//

typedef struct {
    real eta;
    real y1;
    real y2;
    real y3;
    real y4;
    real y5;
    real y6;
    real y7;
    real y8;
    real y9;
    real y10;
} global_D3v2, *global_D3v2_ptr;


#define etaD3(x)    (((global_D3v2_ptr) (x))->eta)
#define DpkD3(x)    (((global_D3v2_ptr) (x))->y1)
#define DppD3(x)    (((global_D3v2_ptr) (x))->y3)
#define D2fD3(x)    (((global_D3v2_ptr) (x))->y5)
#define D2mfD3(x)    (((global_D3v2_ptr) (x))->y7)
#define D3symmD3(x)    (((global_D3v2_ptr) (x))->y9)
#define D3symmprimeD3(x)    (((global_D3v2_ptr) (x))->y10) //Modificacion-LCDMfk
//

// GL structure
typedef struct {
    int npts;
    real x1;
    real x2;
    real *xgl;
    real *wgl;
} global_GL, *global_GL_ptr;

// STATIC problem: gcc version 11
//global_GL_ptr pGL;
global global_GL_ptr pGL;
//~ global_GL_ptr pGL;

#define nGL(x)    (((global_GL_ptr) (x))->npts)
#define x1GL(x)    (((global_GL_ptr) (x))->x1)
#define x2GL(x)    (((global_GL_ptr) (x))->x2)
#define xGL(x)    (((global_GL_ptr) (x))->xgl)
#define wGL(x)    (((global_GL_ptr) (x))->wgl)




typedef struct {
    real k;
    real P0;
    real P2;
    real P4;
} global_rsdmultipoles, *global_rsdmultipoles_ptr;

#define k_rsdmultipoles(x)    (((global_rsdmultipoles_ptr) (x))->k)
#define P0_rsdmultipoles(x)    (((global_rsdmultipoles_ptr) (x))->P0)
#define P2_rsdmultipoles(x)    (((global_rsdmultipoles_ptr) (x))->P2)
#define P4_rsdmultipoles(x)    (((global_rsdmultipoles_ptr) (x))->P4)


// END :: CLPT correlation auxiliary functions and structures



// NEW for v2

typedef struct {
real kmin;
real kmax;
int Nk;
real qmin;
real qmax;
int Nq;		
real rmin;
real rmax;
int Nr;
real smin;
real smax;
int Ns;
} global_output_lists, *global_output_lists_ptr;

global global_output_lists golists;





global real *inout_xval;
global real *inout_yval;
global real *inout_zval;
global real *inout_wval;

global double dxsav,*xp,**yp;
global int kmax,kount;
global int nrhs;

global long idum;                // seed for random generators







typedef struct {
int eta;
real k;
//
real P22dd;
real P22du;
real P22uu;
// A 
real I1udd1A;
real I2uud1A;
real I2uud2A;
real I3uuu2A;
real I3uuu3A;
//  B plus C   
real I2uudd1BpC;
real I2uudd2BpC;
real I3uuud2BpC;
real I3uuud3BpC;
real I4uuuu2BpC;
real I4uuuu3BpC;
real I4uuuu4BpC;
//  Bias
real Pb1b2;
real Pb1bs2;
real Pb22;
real Pb2s2;
real Ps22;
real Pb2theta;
real Pbs2theta;
//
real P13dd;
real P13du;
real P13uu;
real sigma32PSL;
} global_kFs, *global_kFs_ptr;

#define kFs(x)    (((global_kFs_ptr) (x))->k)
#define P22dd(x)    (((global_kFs_ptr) (x))->P22dd)
#define P22du(x)    (((global_kFs_ptr) (x))->P22du)
#define P22uu(x)    (((global_kFs_ptr) (x))->P22uu)
// A 
#define I1udd1A(x)    (((global_kFs_ptr) (x))->I1udd1A)
#define I2uud1A(x)    (((global_kFs_ptr) (x))->I2uud1A)
#define I2uud2A(x)    (((global_kFs_ptr) (x))->I2uud2A)
#define I3uuu2A(x)    (((global_kFs_ptr) (x))->I3uuu2A)
#define I3uuu3A(x)    (((global_kFs_ptr) (x))->I3uuu3A)
//  B plus C   
#define I2uudd1BpC(x)    (((global_kFs_ptr) (x))->I2uudd1BpC)
#define I2uudd2BpC(x)    (((global_kFs_ptr) (x))->I2uudd2BpC)
#define I3uuud2BpC(x)    (((global_kFs_ptr) (x))->I3uuud2BpC)
#define I3uuud3BpC(x)    (((global_kFs_ptr) (x))->I3uuud3BpC)
#define I4uuuu2BpC(x)    (((global_kFs_ptr) (x))->I4uuuu2BpC)
#define I4uuuu3BpC(x)    (((global_kFs_ptr) (x))->I4uuuu3BpC)
#define I4uuuu4BpC(x)    (((global_kFs_ptr) (x))->I4uuuu4BpC)
//  Bias
#define Pb1b2(x)    (((global_kFs_ptr) (x))->Pb1b2)
#define Pb1bs2(x)    (((global_kFs_ptr) (x))->Pb1bs2)
#define Pb22(x)    (((global_kFs_ptr) (x))->Pb22)
#define Pb2s2(x)    (((global_kFs_ptr) (x))->Pb2s2)
#define Ps22(x)    (((global_kFs_ptr) (x))->Ps22)
#define Pb2theta(x)    (((global_kFs_ptr) (x))->Pb2theta)
#define Pbs2theta(x)    (((global_kFs_ptr) (x))->Pbs2theta)
//
#define P13dd(x)    (((global_kFs_ptr) (x))->P13dd)
#define P13du(x)    (((global_kFs_ptr) (x))->P13du)
#define P13uu(x)    (((global_kFs_ptr) (x))->P13uu)
#define sigma32PSL(x)    (((global_kFs_ptr) (x))->sigma32PSL)



typedef struct {
int eta;
real *kT;
//
real *P22ddT;
real *P22duT;
real *P22uuT;
// A 
real *I1udd1AT;
real *I2uud1AT;
real *I2uud2AT;
real *I3uuu2AT;
real *I3uuu3AT;
//  B plus C   
real *I2uudd1BpCT;
real *I2uudd2BpCT;
real *I3uuud2BpCT;
real *I3uuud3BpCT;
real *I4uuuu2BpCT;
real *I4uuuu3BpCT;
real *I4uuuu4BpCT;
//  Bias
real *Pb1b2T;
real *Pb1bs2T;
real *Pb22T;
real *Pb2s2T;
real *Ps22T;
real *Pb2thetaT;
real *Pbs2thetaT;
//
real *P13ddT;
real *P13duT;
real *P13uuT;
real *sigma32PSLT;
real *pklT;
real *fkT;
} global_kFArrays, *global_kFArrays_ptr;


global global_kFArrays kFArrays;
global global_kFArrays kFArraysd;

global global_kFArrays kFArrays_nw;
global global_kFArrays kFArraysd_nw;

#endif // ! _globaldefs_h

