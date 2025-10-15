/*==============================================================================
 NAME: global.c				[code for redshift space correlation function - GSM]
 Alejandro Aviles (avilescervantes@gmail.com), ...
 ================================================================================ 
*/

/* src/global.c */
#define global
#include "clpt_types.h"   // brings full definitions
#include "globaldefs.h"
#include "models.h"
#include "protodefs.h"

q_arrays  qArrays  = {0};
q_arraysd qArraysd = {0};
r_arrays  rArrays  = {0};

void global_variables(void){


	
	golists.rmin=1.0;
	golists.rmax= cmd.smax + 40.;
	golists.Nr =(int)( golists.rmax - golists.rmin + 1.0);
	
	golists.qmin = 0.001;     
    golists.qmax = golists.rmax + 50.;    
    golists.Nq = (int)floor(log10(golists.qmax/golists.qmin)*(real)cmd.NqperLogDecade); 
    
    golists.kmin = cmd.kmin;
    golists.kmax = cmd.kmax;
    golists.Nk = cmd.Nk;
    
    
    golists.smin = cmd.smin;
    golists.smax = cmd.smax;
    golists.Ns = cmd.Ns;
	
	//~ fprintf(stdout,"\nglists.qmin=%g, glists.qmax=%g, glists.Nq=%d",golists.rmin, golists.rmax, golists.Nr);
	//~ fprintf(stdout,"\nglists.rmin=%g, glists.rmax=%g, glists.Nr=%d",golists.qmin, golists.qmax, golists.Nq);

	kFArrays.kT = malloc( golists.Nk * sizeof(real));
	kFArrays.P22ddT = malloc( golists.Nk * sizeof(real));
	kFArrays.P22duT = malloc( golists.Nk * sizeof(real));
	kFArrays.P22uuT = malloc( golists.Nk * sizeof(real));
		// A 
	kFArrays.I1udd1AT = malloc( golists.Nk * sizeof(real));
	kFArrays.I2uud1AT = malloc( golists.Nk * sizeof(real));
	kFArrays.I2uud2AT = malloc( golists.Nk * sizeof(real));
	kFArrays.I3uuu2AT = malloc( golists.Nk * sizeof(real));
	kFArrays.I3uuu3AT = malloc( golists.Nk * sizeof(real));
	//  B plus C   
	kFArrays.I2uudd1BpCT = malloc( golists.Nk * sizeof(real));
	kFArrays.I2uudd2BpCT = malloc( golists.Nk * sizeof(real));
	kFArrays.I3uuud2BpCT = malloc( golists.Nk * sizeof(real));
	kFArrays.I3uuud3BpCT = malloc( golists.Nk * sizeof(real));
	kFArrays.I4uuuu2BpCT = malloc( golists.Nk * sizeof(real));
	kFArrays.I4uuuu3BpCT = malloc( golists.Nk * sizeof(real));
	kFArrays.I4uuuu4BpCT = malloc( golists.Nk * sizeof(real));
	//  Bias
	kFArrays.Pb1b2T = malloc( golists.Nk * sizeof(real));
	kFArrays.Pb1bs2T = malloc( golists.Nk * sizeof(real));
	kFArrays.Pb22T = malloc( golists.Nk * sizeof(real));
	kFArrays.Pb2s2T = malloc( golists.Nk * sizeof(real));
	kFArrays.Ps22T = malloc( golists.Nk * sizeof(real));
	kFArrays.Pb2thetaT = malloc( golists.Nk * sizeof(real));
	kFArrays.Pbs2thetaT = malloc( golists.Nk * sizeof(real));
	//
	kFArrays.P13ddT = malloc( golists.Nk * sizeof(real));
	kFArrays.P13duT = malloc( golists.Nk * sizeof(real));
	kFArrays.P13uuT = malloc( golists.Nk * sizeof(real));
	kFArrays.sigma32PSLT = malloc( golists.Nk * sizeof(real));
	kFArrays.pklT = malloc( golists.Nk * sizeof(real));
	kFArrays.fkT = malloc( golists.Nk * sizeof(real));
	
	
	kFArraysd.kT = malloc( golists.Nk * sizeof(real));
	kFArraysd.P22ddT = malloc( golists.Nk * sizeof(real));
	kFArraysd.P22duT = malloc( golists.Nk * sizeof(real));
	kFArraysd.P22uuT = malloc( golists.Nk * sizeof(real));
		// A 
	kFArraysd.I1udd1AT = malloc( golists.Nk * sizeof(real));
	kFArraysd.I2uud1AT = malloc( golists.Nk * sizeof(real));
	kFArraysd.I2uud2AT = malloc( golists.Nk * sizeof(real));
	kFArraysd.I3uuu2AT = malloc( golists.Nk * sizeof(real));
	kFArraysd.I3uuu3AT = malloc( golists.Nk * sizeof(real));
	//  B plus C   
	kFArraysd.I2uudd1BpCT = malloc( golists.Nk * sizeof(real));
	kFArraysd.I2uudd2BpCT = malloc( golists.Nk * sizeof(real));
	kFArraysd.I3uuud2BpCT = malloc( golists.Nk * sizeof(real));
	kFArraysd.I3uuud3BpCT = malloc( golists.Nk * sizeof(real));
	kFArraysd.I4uuuu2BpCT = malloc( golists.Nk * sizeof(real));
	kFArraysd.I4uuuu3BpCT = malloc( golists.Nk * sizeof(real));
	kFArraysd.I4uuuu4BpCT = malloc( golists.Nk * sizeof(real));
	//  Bias
	kFArraysd.Pb1b2T = malloc( golists.Nk * sizeof(real));
	kFArraysd.Pb1bs2T = malloc( golists.Nk * sizeof(real));
	kFArraysd.Pb22T = malloc( golists.Nk * sizeof(real));
	kFArraysd.Pb2s2T = malloc( golists.Nk * sizeof(real));
	kFArraysd.Ps22T = malloc( golists.Nk * sizeof(real));
	kFArraysd.Pb2thetaT = malloc( golists.Nk * sizeof(real));
	kFArraysd.Pbs2thetaT = malloc( golists.Nk * sizeof(real));
	//
	kFArraysd.P13ddT = malloc( golists.Nk * sizeof(real));
	kFArraysd.P13duT = malloc( golists.Nk * sizeof(real));
	kFArraysd.P13uuT = malloc( golists.Nk * sizeof(real));
	kFArraysd.sigma32PSLT = malloc( golists.Nk * sizeof(real));
	kFArraysd.pklT = malloc( golists.Nk * sizeof(real));
	kFArraysd.fkT = malloc( golists.Nk * sizeof(real));






	kFArrays_nw.kT = malloc( golists.Nk * sizeof(real));
	kFArrays_nw.P22ddT = malloc( golists.Nk * sizeof(real));
	kFArrays_nw.P22duT = malloc( golists.Nk * sizeof(real));
	kFArrays_nw.P22uuT = malloc( golists.Nk * sizeof(real));
		// A 
	kFArrays_nw.I1udd1AT = malloc( golists.Nk * sizeof(real));
	kFArrays_nw.I2uud1AT = malloc( golists.Nk * sizeof(real));
	kFArrays_nw.I2uud2AT = malloc( golists.Nk * sizeof(real));
	kFArrays_nw.I3uuu2AT = malloc( golists.Nk * sizeof(real));
	kFArrays_nw.I3uuu3AT = malloc( golists.Nk * sizeof(real));
	//  B plus C   
	kFArrays_nw.I2uudd1BpCT = malloc( golists.Nk * sizeof(real));
	kFArrays_nw.I2uudd2BpCT = malloc( golists.Nk * sizeof(real));
	kFArrays_nw.I3uuud2BpCT = malloc( golists.Nk * sizeof(real));
	kFArrays_nw.I3uuud3BpCT = malloc( golists.Nk * sizeof(real));
	kFArrays_nw.I4uuuu2BpCT = malloc( golists.Nk * sizeof(real));
	kFArrays_nw.I4uuuu3BpCT = malloc( golists.Nk * sizeof(real));
	kFArrays_nw.I4uuuu4BpCT = malloc( golists.Nk * sizeof(real));
	//  Bias
	kFArrays_nw.Pb1b2T = malloc( golists.Nk * sizeof(real));
	kFArrays_nw.Pb1bs2T = malloc( golists.Nk * sizeof(real));
	kFArrays_nw.Pb22T = malloc( golists.Nk * sizeof(real));
	kFArrays_nw.Pb2s2T = malloc( golists.Nk * sizeof(real));
	kFArrays_nw.Ps22T = malloc( golists.Nk * sizeof(real));
	kFArrays_nw.Pb2thetaT = malloc( golists.Nk * sizeof(real));
	kFArrays_nw.Pbs2thetaT = malloc( golists.Nk * sizeof(real));
	//
	kFArrays_nw.P13ddT = malloc( golists.Nk * sizeof(real));
	kFArrays_nw.P13duT = malloc( golists.Nk * sizeof(real));
	kFArrays_nw.P13uuT = malloc( golists.Nk * sizeof(real));
	kFArrays_nw.sigma32PSLT = malloc( golists.Nk * sizeof(real));
	kFArrays_nw.pklT = malloc( golists.Nk * sizeof(real));
	kFArrays_nw.fkT = malloc( golists.Nk * sizeof(real));
	
	
	kFArraysd_nw.kT = malloc( golists.Nk * sizeof(real));
	kFArraysd_nw.P22ddT = malloc( golists.Nk * sizeof(real));
	kFArraysd_nw.P22duT = malloc( golists.Nk * sizeof(real));
	kFArraysd_nw.P22uuT = malloc( golists.Nk * sizeof(real));
		// A 
	kFArraysd_nw.I1udd1AT = malloc( golists.Nk * sizeof(real));
	kFArraysd_nw.I2uud1AT = malloc( golists.Nk * sizeof(real));
	kFArraysd_nw.I2uud2AT = malloc( golists.Nk * sizeof(real));
	kFArraysd_nw.I3uuu2AT = malloc( golists.Nk * sizeof(real));
	kFArraysd_nw.I3uuu3AT = malloc( golists.Nk * sizeof(real));
	//  B plus C   
	kFArraysd_nw.I2uudd1BpCT = malloc( golists.Nk * sizeof(real));
	kFArraysd_nw.I2uudd2BpCT = malloc( golists.Nk * sizeof(real));
	kFArraysd_nw.I3uuud2BpCT = malloc( golists.Nk * sizeof(real));
	kFArraysd_nw.I3uuud3BpCT = malloc( golists.Nk * sizeof(real));
	kFArraysd_nw.I4uuuu2BpCT = malloc( golists.Nk * sizeof(real));
	kFArraysd_nw.I4uuuu3BpCT = malloc( golists.Nk * sizeof(real));
	kFArraysd_nw.I4uuuu4BpCT = malloc( golists.Nk * sizeof(real));
	//  Bias
	kFArraysd_nw.Pb1b2T = malloc( golists.Nk * sizeof(real));
	kFArraysd_nw.Pb1bs2T = malloc( golists.Nk * sizeof(real));
	kFArraysd_nw.Pb22T = malloc( golists.Nk * sizeof(real));
	kFArraysd_nw.Pb2s2T = malloc( golists.Nk * sizeof(real));
	kFArraysd_nw.Ps22T = malloc( golists.Nk * sizeof(real));
	kFArraysd_nw.Pb2thetaT = malloc( golists.Nk * sizeof(real));
	kFArraysd_nw.Pbs2thetaT = malloc( golists.Nk * sizeof(real));
	//
	kFArraysd_nw.P13ddT = malloc( golists.Nk * sizeof(real));
	kFArraysd_nw.P13duT = malloc( golists.Nk * sizeof(real));
	kFArraysd_nw.P13uuT = malloc( golists.Nk * sizeof(real));
	kFArraysd_nw.sigma32PSLT = malloc( golists.Nk * sizeof(real));
	kFArraysd_nw.pklT = malloc( golists.Nk * sizeof(real));
	kFArraysd_nw.fkT = malloc( golists.Nk * sizeof(real));










/*

    qArrays.qTab = malloc( golists.Nq * sizeof(real));
    qArrays.XLT = malloc( golists.Nq * sizeof(real));
    qArrays.YLT = malloc( golists.Nq * sizeof(real));
    qArrays.XloopT = malloc( golists.Nq * sizeof(real));
    qArrays.YloopT = malloc( golists.Nq * sizeof(real));
    qArrays.preVT = malloc( golists.Nq * sizeof(real));
    qArrays.TT = malloc( golists.Nq * sizeof(real));
    qArrays.X10T = malloc( golists.Nq * sizeof(real));
    qArrays.Y10T = malloc( golists.Nq * sizeof(real));
    qArrays.ULT = malloc( golists.Nq * sizeof(real));
    qArrays.UloopT = malloc( golists.Nq * sizeof(real));
    qArrays.U11T = malloc( golists.Nq * sizeof(real));
    qArrays.U20T = malloc( golists.Nq * sizeof(real));
    qArrays.dotVaT = malloc( golists.Nq * sizeof(real));
    qArrays.ddotVaT = malloc( golists.Nq * sizeof(real));
    qArrays.dotVbT = malloc( golists.Nq * sizeof(real)); 
    qArrays.V10T = malloc( golists.Nq * sizeof(real));
    qArrays.iBessel0T = malloc( golists.Nq * sizeof(real));
    qArrays.iBessel2T = malloc( golists.Nq * sizeof(real));
    qArrays.iBessel4T = malloc( golists.Nq * sizeof(real));
    qArrays.nabla2xiT = malloc( golists.Nq * sizeof(real));

    qArraysd.qTab = malloc( golists.Nq * sizeof(real));
    qArraysd.XLT = malloc( golists.Nq * sizeof(real));
    qArraysd.YLT = malloc( golists.Nq * sizeof(real));
    qArraysd.XloopT = malloc( golists.Nq * sizeof(real));
    qArraysd.YloopT = malloc( golists.Nq * sizeof(real));
    qArraysd.preVT = malloc( golists.Nq * sizeof(real));
    qArraysd.TT = malloc( golists.Nq * sizeof(real));
    qArraysd.X10T = malloc( golists.Nq * sizeof(real));
    qArraysd.Y10T = malloc( golists.Nq * sizeof(real));
    qArraysd.ULT = malloc( golists.Nq * sizeof(real));
    qArraysd.UloopT = malloc( golists.Nq * sizeof(real));
    qArraysd.U11T = malloc( golists.Nq * sizeof(real));
    qArraysd.U20T = malloc( golists.Nq * sizeof(real));
    qArraysd.dotVaT = malloc( golists.Nq * sizeof(real));
    qArraysd.ddotVaT = malloc( golists.Nq * sizeof(real));
    qArraysd.dotVbT = malloc( golists.Nq * sizeof(real)); 
    qArraysd.V10T = malloc( golists.Nq * sizeof(real));
    qArraysd.iBessel0T = malloc( golists.Nq * sizeof(real));
    qArraysd.iBessel2T = malloc( golists.Nq * sizeof(real));
    qArraysd.iBessel4T = malloc( golists.Nq * sizeof(real));
    qArraysd.nabla2xiT = malloc( golists.Nq * sizeof(real));
	
	    rArrays.rTab= malloc(golists.Nr * sizeof(real) );
		rArrays.xi_00T= malloc(golists.Nr * sizeof(real) );
		rArrays.xi_10T= malloc(golists.Nr * sizeof(real) );
		rArrays.xi_20T= malloc(golists.Nr * sizeof(real) );
		rArrays.xi_01T= malloc(golists.Nr * sizeof(real) );
		rArrays.xi_02T= malloc(golists.Nr * sizeof(real) );
		rArrays.xi_11T= malloc(golists.Nr * sizeof(real) );
		rArrays.xi_eftT= malloc(golists.Nr * sizeof(real) );
		rArrays.xi_001T= malloc(golists.Nr * sizeof(real) );
		rArrays.xi_002T= malloc(golists.Nr * sizeof(real) );
		rArrays.xi_101T= malloc(golists.Nr * sizeof(real) );
		rArrays.xi_011T= malloc(golists.Nr * sizeof(real) );         
		rArrays.xi_LT= malloc(golists.Nr * sizeof(real) );
		rArrays.xi_zaT= malloc(golists.Nr * sizeof(real) );
		rArrays.xi_AT= malloc(golists.Nr * sizeof(real) ); 
		rArrays.xi_WT= malloc(golists.Nr * sizeof(real) ); 
		rArrays.v_00T= malloc(golists.Nr * sizeof(real) );
		rArrays.v_10T= malloc(golists.Nr * sizeof(real) );
		rArrays.v_20T= malloc(golists.Nr * sizeof(real) );
		rArrays.v_01T= malloc(golists.Nr * sizeof(real) );
		rArrays.v_11T= malloc(golists.Nr * sizeof(real) );
		rArrays.v_001T= malloc(golists.Nr * sizeof(real) );
		rArrays.v_101T= malloc(golists.Nr * sizeof(real) );
		rArrays.v_eftT= malloc(golists.Nr * sizeof(real) );         
		rArrays.spar_00T= malloc(golists.Nr * sizeof(real) );
		rArrays.spar_10T= malloc(golists.Nr * sizeof(real) );
		rArrays.spar_20T= malloc(golists.Nr * sizeof(real) );
		rArrays.spar_01T= malloc(golists.Nr * sizeof(real) );
		rArrays.spar_001T= malloc(golists.Nr * sizeof(real) );            
		rArrays.sperp_00T= malloc(golists.Nr * sizeof(real) );
		rArrays.sperp_10T= malloc(golists.Nr * sizeof(real) );
		rArrays.sperp_20T= malloc(golists.Nr * sizeof(real) );
		rArrays.sperp_01T= malloc(golists.Nr * sizeof(real) );
		rArrays.sperp_001T= malloc(golists.Nr * sizeof(real) );
		
	    rArraysd.rTab= malloc(golists.Nr * sizeof(real) );
		rArraysd.xi_00T= malloc(golists.Nr * sizeof(real) );
		rArraysd.xi_10T= malloc(golists.Nr * sizeof(real) );
		rArraysd.xi_20T= malloc(golists.Nr * sizeof(real) );
		rArraysd.xi_01T= malloc(golists.Nr * sizeof(real) );
		rArraysd.xi_02T= malloc(golists.Nr * sizeof(real) );
		rArraysd.xi_11T= malloc(golists.Nr * sizeof(real) );
		rArraysd.xi_eftT= malloc(golists.Nr * sizeof(real) );
		rArraysd.xi_001T= malloc(golists.Nr * sizeof(real) );
		rArraysd.xi_002T= malloc(golists.Nr * sizeof(real) );
		rArraysd.xi_101T= malloc(golists.Nr * sizeof(real) );
		rArraysd.xi_011T= malloc(golists.Nr * sizeof(real) );         
		rArraysd.xi_LT= malloc(golists.Nr * sizeof(real) );
		rArraysd.xi_zaT= malloc(golists.Nr * sizeof(real) );
		rArraysd.xi_AT= malloc(golists.Nr * sizeof(real) ); 
		rArraysd.xi_WT= malloc(golists.Nr * sizeof(real) ); 
		rArraysd.v_00T= malloc(golists.Nr * sizeof(real) );
		rArraysd.v_10T= malloc(golists.Nr * sizeof(real) );
		rArraysd.v_20T= malloc(golists.Nr * sizeof(real) );
		rArraysd.v_01T= malloc(golists.Nr * sizeof(real) );
		rArraysd.v_11T= malloc(golists.Nr * sizeof(real) );
		rArraysd.v_001T= malloc(golists.Nr * sizeof(real) );
		rArraysd.v_101T= malloc(golists.Nr * sizeof(real) );
		rArraysd.v_eftT= malloc(golists.Nr * sizeof(real) );         
		rArraysd.spar_00T= malloc(golists.Nr * sizeof(real) );
		rArraysd.spar_10T= malloc(golists.Nr * sizeof(real) );
		rArraysd.spar_20T= malloc(golists.Nr * sizeof(real) );
		rArraysd.spar_01T= malloc(golists.Nr * sizeof(real) );
		rArraysd.spar_001T= malloc(golists.Nr * sizeof(real) );            
		rArraysd.sperp_00T= malloc(golists.Nr * sizeof(real) );
		rArraysd.sperp_10T= malloc(golists.Nr * sizeof(real) );
		rArraysd.sperp_20T= malloc(golists.Nr * sizeof(real) );
		rArraysd.sperp_01T= malloc(golists.Nr * sizeof(real) );
		rArraysd.sperp_001T= malloc(golists.Nr * sizeof(real) );		
		
	
	*/
	
	
}
