


#include <math.h>
#include "globaldefs.h"
#include "protodefs.h"

local  real Interpolation_nr(real k, double kPS[], double pPS[], int nPS, double pPS2[]);

local real PSLF(real k);
local real Q1F(real k);
local real Q2F(real k);
local real Q3F(real k);
local real Q5F(real k);
local real Q8F(real k);
local real Qs2F(real k);
local real R1F(real k);
local real R2F(real k);

local void InputQsRsTable(void);

global_qfunctions qfunctions(real qi);
global_corrfunctions correlation_functions(real qi);

local void postprocess_string_to_int(string, int *);

// Qs and Rs table structure
typedef struct _pointQsRsTable {
    real k;               //1
    real Q1;            //2
    real Q2;            //3
    real Q3;            //4
    real Q5;            //5
    real Q8;            //6
    real Qs2;          //7
    real R1;            //8
    real R2;            //9
    real Dpk;          //10
    real PSMGL;     //11
} pointQsRsTable, *pointQsRsTableptr;

local int nQsRsTable;
local pointQsRsTableptr QsRstab;

#define kQsRs(x)    (((pointQsRsTableptr) (x))->k)
#define Q1QsRs(x)    (((pointQsRsTableptr) (x))->Q1)
#define Q2QsRs(x)    (((pointQsRsTableptr) (x))->Q2)
#define Q3QsRs(x)    (((pointQsRsTableptr) (x))->Q3)
#define Q5QsRs(x)    (((pointQsRsTableptr) (x))->Q5)
#define Q8QsRs(x)    (((pointQsRsTableptr) (x))->Q8)
#define Qs2QsRs(x)    (((pointQsRsTableptr) (x))->Qs2)
#define R1QsRs(x)    (((pointQsRsTableptr) (x))->R1)
#define R2QsRs(x)    (((pointQsRsTableptr) (x))->R2)
#define DpkQsRs(x)    (((pointQsRsTableptr) (x))->Dpk)
#define PSMGLQsRs(x)    (((pointQsRsTableptr) (x))->PSMGL)

local pointQsRsTableptr PQsRstab;

local real *kTab;
local real *Q1T;
local real *Q1T2;
local real *Q2T;
local real *Q2T2;
local real *Q3T;
local real *Q3T2;
local real *Q5T;
local real *Q5T2;
local real *Q8T;
local real *Q8T2;
local real *Qs2T;
local real *Qs2T2;
local real *R1T;
local real *R1T2;
local real *R2T;
local real *R2T2;
local real *PSLMGT;
local real *PSLMGT2;
