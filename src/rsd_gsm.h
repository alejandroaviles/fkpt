#include "globaldefs.h"
#include "protodefs.h"



global_rsdmultipoles CF_rsd_multipoles(real s, real b1, real b2, real bs, real alpha1eft, real alpha2eft, real sigma2eft);

local void InputrfunctionsTable(void);
local  real Interpolation_nr(real k, double kPS[], double pPS[], int nPS, double pPS2[]);


local real xi12F(real r, real b1, real b2, real bs, real alpha1eft, real alpha2eft, real sigma2eft);
local real v12F(real r, real b1, real b2, real bs, real alpha1eft, real alpha2eft, real sigma2eft);
local real sigma2par12F(real r, real b1, real b2, real bs, real alpha1eft, real alpha2eft, real sigma2eft);
local real sigma2perp12F(real r, real b1, real b2, real bs, real alpha1eft, real alpha2eft, real sigma2eft);

local real xi_zaF(real r);

local real xi_00F(real r);
local real xi_10F(real r);
local real xi_20F(real r);
local real xi_01F(real r);
local real xi_02F(real r);
local real xi_11F(real r);
local real xi_eftF(real r);
local real xi_001F(real r);
local real xi_002F(real r);
local real xi_101F(real r);
local real xi_011F(real r);

local real v_00F(real r);
local real v_10F(real r);
local real v_20F(real r);
local real v_01F(real r);
local real v_11F(real r);
local real v_001F(real r);
local real v_101F(real r);
local real v_eftF(real r);


local real spar_00F(real r);
local real spar_10F(real r);
local real spar_20F(real r);
local real spar_01F(real r);
local real spar_001F(real r);

local real sperp_00F(real r);
local real sperp_10F(real r);
local real sperp_20F(real r);
local real sperp_01F(real r);
local real sperp_001F(real r);


// rfunctions table structure
typedef struct _pointrfunctionsTable {
    real r;         //1
    real xi_00;  //2   (xi_clpt: matter - clpt)
    real xi_10;  //3
    real xi_20;  //4
    real xi_01;   //5
    real xi_02;   //6
    real xi_11;    //7
    real xi_eft;   //8
    real xi_001;    //9
    real xi_002;   //10
    real xi_101;    //11
    real xi_011;     //12
	//
	real xi_za; //14
    //
    real v_00;   //17
    real v_10;    //18
    real v_20;     //19
    real v_01;     //20
    real v_11;     //21
    real v_001;    //22
    real v_101;     //23
    real v_eft;      //24
	// 
	real spar_00;      //25
	real spar_10;      //26
	real spar_20;      //27
	real spar_01;      //28
	real spar_001;        //29
	//
	real sperp_00;    //30
	real sperp_10;   //31
	real sperp_20;    //32
	real sperp_01;     //33
	real sperp_001;     //34
} pointrfunctionsTable, *pointrfunctionsTableptr;

//~ local int nrfunctionsTable;
local pointrfunctionsTableptr rfunctionstab;


