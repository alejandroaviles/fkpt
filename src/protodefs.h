/*==============================================================================
 HEADER: protodefs.h				[fkpt]
 ==============================================================================*/

#ifndef _protodefs_h
#define _protodefs_h

void integration_method_string_to_int(string,int *);
void quadraturemethod_string_to_int(string,int *);

//void output(void);

void MainLoop(void);
void StartRun(string, string, string, string);
void StartOutput(void);
void EndRun(void);

global void global_variables(void);
global void compute_kfunctions(void);
global void compute_rsdmultipoles(void);
global void compute_qfunctions(void);
global void compute_clpt(void);
global void compute_gsm(void);
global void write(void);
global void free_variables(void);


//~ global void PostProcessing(void); //BORRAR;

// MGLPT DIFFEQS and QUAD
global real DpFunction(real k);
global real DpFunction_LCDM(real k);
global real f_growth(real k);
global real f_growth_LCDM(void);
global global_D2v2_ptr DsSecondOrder_func(real kf, real k1, real k2);
global global_D3v2_ptr DsThirdOrder_func(real x, real k, real p);

// MGLPT QUADS
//~ global_kFs ki_functions_driver(real eta, real ki);
//~ global_kFs ki_functions_driver(real ki);
global_kFs ki_functions_driver(real ki, double kPKL[], double pPKL[], int nPKLT, double pPKL2[]);

// I/O directories:
global void setFilesDirs_log(void);
global void setFilesDirs(void);

global  real psInterpolation_nr(real k, double kPS[], double pPS[], int nPS);
//

#endif // ! _protodefs_h
