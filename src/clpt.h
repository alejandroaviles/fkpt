#include "globaldefs.h"
#include "protodefs.h"


//~ #define KK  5
//~ #define EPSQ 1.0e-6

local  real Interpolation_nr(real k, double kPS[], double pPS[], int nPS, double pPS2[]);

local real sigma2L_function_int(real y);
local real sigma2L_function(void);

// BEGIN :: CLPT correlation auxiliary functions and structures
local real XLF(real q);
local real YLF(real q);
local real XloopF(real q);
local real YloopF(real q);
local real preVF(real q);
local real TF(real q);
local real X10F(real q);
local real Y10F(real q);
local real ULF(real q);
local real UloopF(real q);
local real U11F(real q);
local real U20F(real q);
local real dotVaF(real q);
local real ddotVaF(real q);
local real dotVbF(real q);
local real V10F(real q);
local real iBessel0F(real q);
local real iBessel2F(real q);
local real iBessel4F(real q);
local real nabla2xiF(real q);

//~ DERIVED
local real UF(real q);
local real fXF(real q);
local real hYF(real q);



local real MZAF(real y, real qminusmur, real fx, real hy);



//~ // qfunctions table structure
//~ typedef struct _pointqfunctionsTable {
    //~ real q;         //1
    //~ real XL;        //2
    //~ real YL;        //3
    //~ real Xloop;     //4
    //~ real Yloop;     //5
    //~ real preV;         //6
    //~ real T;         //7
    //~ real X10;       //8
    //~ real Y10;       //9
    //~ real UL;        //10
    //~ real Uloop;     //11
    //~ real U11;       //12
    //~ real U20;       //13
    //~ real dotVa;       //14
    //~ real ddotVa;     //15
    //~ real dotVb;  //16
    //~ real V10;  //17
    //~ real iBessel0;  //18
    //~ real iBessel2;  //19
    //~ real iBessel4;  //20
    //~ real nabla2xi;  // 21   
//~ } pointqfunctionsTable, *pointqfunctionsTableptr;

//~ local int nqfunctionsTable;
//~ local pointqfunctionsTableptr qfunctionsstab;

//~ #define qqfuncs(x)      (((pointqfunctionsTableptr) (x))->q)
//~ #define XLqfuncs(x)     (((pointqfunctionsTableptr) (x))->XL)
//~ #define YLqfuncs(x)     (((pointqfunctionsTableptr) (x))->YL)
//~ #define Xloopqfuncs(x)     (((pointqfunctionsTableptr) (x))->Xloop)
//~ #define Yloopqfuncs(x)     (((pointqfunctionsTableptr) (x))->Yloop)
//~ #define preVqfuncs(x)     (((pointqfunctionsTableptr) (x))->preV)
//~ #define Tqfuncs(x)     (((pointqfunctionsTableptr) (x))->T)
//~ #define X10qfuncs(x)     (((pointqfunctionsTableptr) (x))->X10)
//~ #define Y10qfuncs(x)    (((pointqfunctionsTableptr) (x))->Y10)
//~ #define ULqfuncs(x)    (((pointqfunctionsTableptr) (x))->UL)
//~ #define Uloopqfuncs(x)    (((pointqfunctionsTableptr) (x))->Uloop)
//~ #define U11qfuncs(x)     (((pointqfunctionsTableptr) (x))->U11)
//~ #define U20qfuncs(x)     (((pointqfunctionsTableptr) (x))->U20)
//~ #define dotVaqfuncs(x)     (((pointqfunctionsTableptr) (x))->dotVa)
//~ #define ddotVaqfuncs(x)     (((pointqfunctionsTableptr) (x))->ddotVa)
//~ #define dotVbqfuncs(x)     (((pointqfunctionsTableptr) (x))->dotVb)
//~ #define V10qfuncs(x)     (((pointqfunctionsTableptr) (x))->V10)
//~ #define iBessel0qfuncs(x)     (((pointqfunctionsTableptr) (x))->iBessel0)
//~ #define iBessel2qfuncs(x)     (((pointqfunctionsTableptr) (x))->iBessel2)
//~ #define iBessel4qfuncs(x)     (((pointqfunctionsTableptr) (x))->iBessel4)
//~ #define nabla2xiqfuncs(x)     (((pointqfunctionsTableptr) (x))->nabla2xi)



//~ local pointqfunctionsTableptr Pqfunctab;

//~ global_zacorrfunctions zacorrelation_functions(real ri, real *muGL, real *wGL, int Nmu);
global_zacorrfunctions zacorrelation_functions(real ri, real sigmav);

global_clptcorrfunctions clptcorrelation_functions(real ri, real sigmav);

local void InputqfunctionsTable(void);



