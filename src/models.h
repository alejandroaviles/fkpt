/*==============================================================================
	HEADER: models.h			[model]
==============================================================================*/

#ifndef _models_h
#define _models_h

global int model_int_flag;
#define LCDM 2

global real KA_LCDM;
global real KB_LCDM;
global real KR1_LCDM;
global real KAp_LCDM; //Modificacion-LCDMfk
global real KR1p_LCDM; //Modificacion-LCDMfk

global void set_model(void);

// functions (mgpt_fns)
global real OmM(real eta);
global real H(real eta);
global real f1(real eta);
global real f2(real eta);
global real A0(real eta);
//

// ==========================================
// Model global HEADERS
global real mu(real eta, real k);
global real sourceA(real eta, real kf, real k1, real k2);
global real sourceb(real eta, real kf, real k1, real k2);
global real kpp(real x, real k, real p);
global real SD2(real eta, real x, real k, real p);
global real S3I(real eta, real x, real k, real p, real Dpk, real Dpp,
                real D2f, real D2mf);
global real S3II(real eta, real x, real k, real p, real Dpk, real Dpp,
                 real D2f, real D2mf);
global real S3FL(real eta, real x, real k, real p, real Dpk, real Dpp, real D2f, real D2mf);
global real S3dI(real eta, real x, real k, real p, real Dpk, real Dpp,
                 real D2f, real D2mf);
// ==========================================

// User model :: reading/printing parameters
global void ReadMGModelParameterFile(char *);
global void PrintMGModelParameterFile(char *);

#endif // !_models_h
