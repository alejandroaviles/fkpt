/*==============================================================================
 NAME: kfunctions.c				[code for fk - Perturbation Theory]
 Alejandro Aviles (avilescervantes@gmail.com)
 ================================================================================ 
*/



#include "globaldefs.h"
#include "protodefs.h"
#include "models.h"



local void quadrature(real ki);
local void k_functions(void);
//~ local void pk_non_wiggle(void);

local real Interpolation_nr(real k, double kPS[], double pPS[], int nPS, double pPS2[]);

local real sigma2L_function_int(real y);
local real Sigma2_int(real y);
local real deltaSigma2_int(real y);
local real sigma2v_function_int(real y);
local real sigma_constants(void);

#define KMIN    1.0e-20
#define _KERNELS_LCDMfk_ 0

global void compute_kfunctions(void)
{
    stream outstrQsRs, outtables;
    real dk;
    real bTime;
    real kBAOmin=0.005, kBAOmax=1.0, epsquadsave;
    int iBAOmin, iBAOmax;
    global_D2v2_ptr ptmp;
    global_D3v2_ptr ptmpR1;
    real fR0save;
    //~ real sigma2psi;


    bTime = second();
 
    gd.f0=f_growth_LCDM();  //Modificacion-LCDMfk
    
    ptmp = DsSecondOrder_func(KMIN, KMIN, KMIN);
    KA_LCDM = DA2D2(ptmp) / ( (3.0/7.0) * Dpk1D2(ptmp) * Dpk2D2(ptmp) );
    KAp_LCDM = DA2primeD2(ptmp) 
                / ( (3.0/7.0) * Dpk1D2(ptmp) * Dpk2D2(ptmp) ) -
               2.0 * DA2D2(ptmp) 
                / ( (3.0/7.0) * Dpk1D2(ptmp) * Dpk2D2(ptmp) )* gd.f0; //Modificacion-LCDMfk
    KB_LCDM = KA_LCDM;

    ptmpR1 = DsThirdOrder_func(0.0000001, KMIN, KMIN);   
    KR1_LCDM = (21.0/5.0)*D3symmD3(ptmpR1)
        /( DpkD3(ptmpR1)*DppD3(ptmpR1)*DppD3(ptmpR1) );     
    KR1p_LCDM = (21.0/5.0)*D3symmprimeD3(ptmpR1)
        /( DpkD3(ptmpR1)*DppD3(ptmpR1)*DppD3(ptmpR1) )/(3.*gd.f0);   //Modificacion-LCDMfk     
        
   fprintf(stdout,"\nA_LCDM=%g, Ap_LCDM=%g, KR1_LCDM = %g, KR1p_LCDM = %g"
   ,KA_LCDM, KAp_LCDM, KR1_LCDM, KR1p_LCDM);

  
    sigma_constants();
        fprintf(stdout,"\nsigma quadratures from kmin = %g to kmax = %g", kPS[1],kPS[nPSLT]);
    
    //~ gd.sigma2L = 21.0466;  
    //~ gd.sigma2v = 21.0466;
    
    //~ fprintf(stdout,"\nA_LCDM=%g,  KR1_LCDM = %g,  f0 = %g",
			//~ KA_LCDM, KR1_LCDM, gd.f0);
    fprintf(stdout,"\ns2psi = %g,   s2v = %g,   Sigma2 = %g,   deltaSigma2 = %g",
			gd.sigma2L, gd.sigma2v, gd.Sigma2, gd.deltaSigma2);
    
    
    
			//~ gd.Sigma2=33.4765970542241;
			//~ gd.deltaSigma2 = 8.084891440805126;
    fprintf(stdout,"\ns2psi = %g,   s2v = %g,   Sigma2 = %g,   deltaSigma2 = %g",
			gd.sigma2L, gd.sigma2v, gd.Sigma2, gd.deltaSigma2);
			
    fprintf(stdout,"\nk-functions:");
    fprintf(stdout," Nk=%d values from kmin=%g to kmax=%g ",
            cmd.Nk, cmd.kmin, cmd.kmax);


    k_functions();
    
    fprintf(stdout,"...time = %g seconds",second()-bTime);    

}
#undef KMIN



#define _K_LOGSPACED_  1
local void k_functions(void)
{
    global_kFs qrs, qrs_nw;
    real aTime;
    real kval, ki, pkl, pkl_nw, fk;
    //~ int counter;
    int i;
    real dk;
    if (_K_LOGSPACED_ ==1){
		dk = (rlog10(cmd.kmax/cmd.kmin))/((real)(cmd.Nk - 1));
	} else {
		dk = (cmd.kmax-cmd.kmin)/((real)(cmd.Nk - 1));
	}

    for (i=1; i<=cmd.Nk; i++) {
		if (_K_LOGSPACED_ ==1){
			kval = rlog10(cmd.kmin) + dk*((real)(i - 1));
			ki = rpow(10.0,kval);
		} else {
			ki = cmd.kmin + dk*((real)(i - 1));
		}
        qrs = ki_functions_driver(ki,kPS, pPS, nPSLT, pPS2);
        qrs_nw = ki_functions_driver(ki,kPS, pPS_nw, nPSLT, pPS2_nw);

        pkl = psInterpolation_nr(ki, kPS, pPS, nPSLT);
        fk = Interpolation_nr(ki, kPS, fkT, nPSLT, fkT2);   
        pkl_nw = Interpolation_nr(ki, kPS, pPS_nw, nPSLT, pPS2_nw);      
		

		kFArrays.kT[i-1]           =  ki ;
		kFArrays.P22ddT[i-1]       =  qrs.P22dd;
		kFArrays.P22duT[i-1]       =  qrs.P22du;
		kFArrays.P22uuT[i-1]       =  qrs.P22uu;
		// A 
		kFArrays.I1udd1AT[i-1]     =  qrs.I1udd1A;
		kFArrays.I2uud1AT[i-1]     =  qrs.I2uud1A;
		kFArrays.I2uud2AT[i-1]     =  qrs.I2uud2A;
		kFArrays.I3uuu2AT[i-1]     =  qrs.I3uuu2A;
		kFArrays.I3uuu3AT[i-1]     =  qrs.I3uuu3A;
		//  B plus C   
		kFArrays.I2uudd1BpCT[i-1]  =  qrs.I2uudd1BpC;
		kFArrays.I2uudd2BpCT[i-1]  =  qrs.I2uudd2BpC;
		kFArrays.I3uuud2BpCT[i-1]  =  qrs.I3uuud2BpC;
		kFArrays.I3uuud3BpCT[i-1]  =  qrs.I3uuud3BpC;
		kFArrays.I4uuuu2BpCT[i-1]  =  qrs.I4uuuu2BpC;
		kFArrays.I4uuuu3BpCT[i-1]  =  qrs.I4uuuu3BpC;
		kFArrays.I4uuuu4BpCT[i-1]  =  qrs.I4uuuu4BpC;
		//  Bias
		kFArrays.Pb1b2T[i-1]       =  qrs.Pb1b2;
		kFArrays.Pb1bs2T[i-1]      =  qrs.Pb1bs2;
		kFArrays.Pb22T[i-1]        =  qrs.Pb22;
		kFArrays.Pb2s2T[i-1]       =  qrs.Pb2s2;
		kFArrays.Ps22T[i-1]        =  qrs.Ps22;
		kFArrays.Pb2thetaT[i-1]    =  qrs.Pb2theta;
		kFArrays.Pbs2thetaT[i-1]   =  qrs.Pbs2theta;
		//
		kFArrays.P13ddT[i-1]       =  qrs.P13dd;
		kFArrays.P13duT[i-1]       =  qrs.P13du;
		kFArrays.P13uuT[i-1]       =  qrs.P13uu;
		kFArrays.sigma32PSLT[i-1]  =  qrs.sigma32PSL;
		kFArrays.pklT[i-1]         =  pkl;
		kFArrays.fkT[i-1]          =  fk;
		
	
		
		
		kFArrays_nw.kT[i-1]           =  ki ;
		kFArrays_nw.P22ddT[i-1]       =  qrs_nw.P22dd;
		kFArrays_nw.P22duT[i-1]       =  qrs_nw.P22du;
		kFArrays_nw.P22uuT[i-1]       =  qrs_nw.P22uu;
		// A 
		kFArrays_nw.I1udd1AT[i-1]     =  qrs_nw.I1udd1A;
		kFArrays_nw.I2uud1AT[i-1]     =  qrs_nw.I2uud1A;
		kFArrays_nw.I2uud2AT[i-1]     =  qrs_nw.I2uud2A;
		kFArrays_nw.I3uuu2AT[i-1]     =  qrs_nw.I3uuu2A;
		kFArrays_nw.I3uuu3AT[i-1]     =  qrs_nw.I3uuu3A;
		//  B plus C   
		kFArrays_nw.I2uudd1BpCT[i-1]  =  qrs_nw.I2uudd1BpC;
		kFArrays_nw.I2uudd2BpCT[i-1]  =  qrs_nw.I2uudd2BpC;
		kFArrays_nw.I3uuud2BpCT[i-1]  =  qrs_nw.I3uuud2BpC;
		kFArrays_nw.I3uuud3BpCT[i-1]  =  qrs_nw.I3uuud3BpC;
		kFArrays_nw.I4uuuu2BpCT[i-1]  =  qrs_nw.I4uuuu2BpC;
		kFArrays_nw.I4uuuu3BpCT[i-1]  =  qrs_nw.I4uuuu3BpC;
		kFArrays_nw.I4uuuu4BpCT[i-1]  =  qrs_nw.I4uuuu4BpC;
		//  Bias
		kFArrays_nw.Pb1b2T[i-1]       =  qrs_nw.Pb1b2;
		kFArrays_nw.Pb1bs2T[i-1]      =  qrs_nw.Pb1bs2;
		kFArrays_nw.Pb22T[i-1]        =  qrs_nw.Pb22;
		kFArrays_nw.Pb2s2T[i-1]       =  qrs_nw.Pb2s2;
		kFArrays_nw.Ps22T[i-1]        =  qrs_nw.Ps22;
		kFArrays_nw.Pb2thetaT[i-1]    =  qrs_nw.Pb2theta;
		kFArrays_nw.Pbs2thetaT[i-1]   =  qrs_nw.Pbs2theta;
		//
		kFArrays_nw.P13ddT[i-1]       =  qrs_nw.P13dd;
		kFArrays_nw.P13duT[i-1]       =  qrs_nw.P13du;
		kFArrays_nw.P13uuT[i-1]       =  qrs_nw.P13uu;
		kFArrays_nw.sigma32PSLT[i-1]  =  qrs_nw.sigma32PSL;
		kFArrays_nw.pklT[i-1]         =  pkl_nw;
		kFArrays_nw.fkT[i-1]          =  fk;		
		
		//~ fprintf(stdout,"(k=%e, f=%e), ",ki,kFArrays.fkT[i-1]);
    }
}

#undef _K_LOGSPACED_


#define QROMBERG     qromo
#define KK  5


/*
// BEGIN Qs and Rs
// kk is the inner integration moment p. 
// kk = ki * r, so usually kk is called p
//~ global_kFs ki_functions(real eta, real ki)
global_kFs ki_functions(real ki, double kPKL[], double pPKL[], int nPKLT, double pPKL2[])
{
    int i, j;
    real pkl_k;
     
    real PSLA, PSLB, psl;
	real fk, fp, fkmp, pklp, pklkmp;
    real rmin, rmax;
    real r, deltar, r2, kmp2;
    real mumin, mumax;
    real x, w, x2;
    real psl1;
    int Nx, nGL;
    real ypi, dk;
    real *xxGL, *wwGL, *xGL, *wGL;
    real kmin, kmax, pmin, pmax;
    
	real P22dd_p =0.0, P22dd_A = 0.0, P22dd_B = 0.0;
    real P22du_p =0.0, P22du_A = 0.0, P22du_B = 0.0;
    real P22uu_p =0.0, P22uu_A = 0.0, P22uu_B = 0.0;
    
    real I1udd1tA_p =0.0, I1udd1tA_A = 0.0, I1udd1tA_B = 0.0;
    real I2uud1tA_p =0.0, I2uud1tA_A = 0.0, I2uud1tA_B = 0.0;
    real I2uud2tA_p =0.0, I2uud2tA_A = 0.0, I2uud2tA_B = 0.0;
    real I3uuu2tA_p =0.0, I3uuu2tA_A = 0.0, I3uuu2tA_B = 0.0;
    real I3uuu3tA_p =0.0, I3uuu3tA_A = 0.0, I3uuu3tA_B = 0.0;
    

    real I2uudd1BpC_p =0.0, I2uudd1BpC_A =0.0, I2uudd1BpC_B =0.0;
    real I2uudd2BpC_p =0.0, I2uudd2BpC_A =0.0, I2uudd2BpC_B =0.0 ;
    real I3uuud2BpC_p =0.0, I3uuud2BpC_A =0.0, I3uuud2BpC_B =0.0 ;
    real I3uuud3BpC_p =0.0, I3uuud3BpC_A =0.0, I3uuud3BpC_B =0.0 ;
    real I4uuuu2BpC_p =0.0, I4uuuu2BpC_A =0.0, I4uuuu2BpC_B =0.0 ;
    real I4uuuu3BpC_p =0.0, I4uuuu3BpC_A =0.0, I4uuuu3BpC_B =0.0 ;
    real I4uuuu4BpC_p =0.0, I4uuuu4BpC_A =0.0, I4uuuu4BpC_B =0.0 ;
      
    real  Pb1b2_p=0.0,     Pb1b2_A =0.0,     Pb1b2_B =0.0;
    real  Pb1bs2_p=0.0,    Pb1bs2_A =0.0,    Pb1bs2_B =0.0; 
    real  Pb22_p=0.0,      Pb22_A =0.0,      Pb22_B =0.0; 
    real  Pb2s2_p=0.0,     Pb2s2_A =0.0,     Pb2s2_B =0.0; 
    real  Ps22_p=0.0,      Ps22_A =0.0,      Ps22_B =0.0; 
    real  Pb2theta_p=0.0,  Pb2theta_A =0.0,  Pb2theta_B =0.0; 
    real  Pbs2theta_p=0.0, Pbs2theta_A =0.0, Pbs2theta_B =0.0; 

    real KP22dd, KP22du, KP22uu;
    real KI1udd1tA, KI2uud1tA, KI2uud2tA, KI3uuu2tA, KI3uuu3tA;
    real KI2uudd1BpC, KI2uudd2BpC, KI3uuud2BpC, KI3uuud3BpC; 
    real KI4uuuu2BpC, KI4uuuu3BpC, KI4uuuu4BpC;
    real KPb1b2, KPb1bs2, KPb22, KPb2s2, KPs22, KPb2theta, KPbs2theta;
    
	real P13dd_p =0.0, P13dd_A = 0.0, P13dd_B = 0.0;
    real P13du_p =0.0, P13du_A = 0.0, P13du_B = 0.0;
    real P13uu_p =0.0, P13uu_A = 0.0, P13uu_B = 0.0;

    real sigma32PSL_p =0.0, sigma32PSL_A = 0.0, sigma32PSL_B = 0.0;
    
    real I1udd1a_p =0.0, I1udd1a_A = 0.0, I1udd1a_B = 0.0;
    real I2uud1a_p =0.0, I2uud1a_A = 0.0, I2uud1a_B = 0.0;
    real I2uud2a_p =0.0, I2uud2a_A = 0.0, I2uud2a_B = 0.0;
    real I3uuu2a_p =0.0, I3uuu2a_A = 0.0, I3uuu2a_B = 0.0;
    real I3uuu3a_p =0.0, I3uuu3a_A = 0.0, I3uuu3a_B = 0.0;
    
    real KP13dd, KP13du, KP13uu, Ksigma32PSL;
    real KI1udd1a, KI2uud1a, KI2uud2a, KI3uuu2a, KI3uuu3a;
    
    real *kk, *dkk;
    //
    pointPSTableptr p;
    //
    global_kFs_ptr QRstmp;
    
    QRstmp = (global_kFs_ptr) allocate(1 * sizeof(global_kFs));
    
    
    //Modificacion   
    kmin = kPS[1];
    kmax = kPS[nPSLT];
    pmin = MAX(kmin,0.01*cmd.kmin);
    pmax = MIN(kmax,16.0*cmd.kmax);
    //~ fprintf(stdout," %f, %f, ",pmin,pmax);
   
    
    dk = (rlog10(pmax) - rlog10(pmin))/((real)(cmd.nquadSteps - 1));
    kk=dvector(1,cmd.nquadSteps);
    dkk=dvector(1,cmd.nquadSteps);
    kk[1] = rpow(10.0,rlog10(pmin));
    for (i=2; i<cmd.nquadSteps; i++) {
        ypi = rlog10(pmin) + dk*((real)(i - 1));
        kk[i] = rpow(10.0,ypi);
        dkk[i] = (kk[i]-kk[i-1]);
    }
    
    
    
//
// Q functions


//
    PSLA = 0.0;
    rmax = kmax/ki;
    rmin = rmin/ki;
    Nx=10;
    xxGL=dvector(1,Nx);
    wwGL=dvector(1,Nx);
    
	for (i=2; i<cmd.nquadSteps; i++) {
		r = kk[i]/ki;
        r2= r*r;
		//~ PSLB = psInterpolation_nr(kk[i], kPS, pPKL, nPSLT);
		PSLB = Interpolation_nr(kk[i], kPKL, pPKL, nPKLT, pPKL2);
		
		pklp=PSLB;
		//~ fp = 1.0;  // This is f(kk[i])/f0  We need to interpolate f
		fp = Interpolation_nr(kk[i], kPS, fkT, nPSLT, fkT2);
		fp /= gd.f0;
		mumin = MAX( -1.0, (1.0 + rsqr(r) - rsqr(rmax)) / (2.0*r)  );
		mumax = MIN(   1.0, (1.0  + rsqr(r) - rsqr(rmin)) / (2.0*r)  );
        
		if (r>=0.5)
			mumax = 0.5/r;
			gauleg(mumin,mumax,xxGL,wwGL,Nx);
        for (j=1; j<=Nx; j++) {
            x = xxGL[j];
            w = wwGL[j];
            x2= x*x;
            
            kmp2=1.0 + r2 - 2.0 * r * x;
            //~ psl = Interpolation_nr(ki * rsqrt(kmp2), kPKL, pPKL, nPKLT, pPKLT2); 
            psl = Interpolation_nr(ki * rsqrt(kmp2), kPKL, pPKL, nPKLT, pPKL2); 
            pklkmp=psl;
			//~ fkmp = 1.0;  // This is f(ki * rsqrt(kmp2))/f0  We need to interpolate f
			fkmp = Interpolation_nr(ki * rsqrt(kmp2), kPS, fkT, nPSLT, fkT2);
			fkmp /= gd.f0;

			KP22dd = rpow( 7 * x + r * (3 - 10 * x2), 2) / (98. * kmp2 * kmp2);
			KP22du =( -7*x + r*(-3 + 10*x2) ) * (fp*r*(-3 + 7*r*x - 4*x2) + 
				fkmp*(4*r - 7*x*(1 + r2) + 10*r*x2)) / (98. * kmp2 * kmp2);
			KP22uu = rpow( fkmp * ( 7*x + r*(-4 + 7*r*x - 10*x2) ) 
			    + fp*r*(3 - 7*r*x + 4*x2), 2) / (98. * kmp2 * kmp2);

// A
			KI1udd1tA =  (fkmp*r*(-1.+r*x) + fp * x * (-1. -r2 + 2.*r*x) ) *  (-7*x + r*(-3. + 10*x2) )    
				/ (7.* kmp2*kmp2);

			KI2uud1tA = -(fkmp*fp*r*(-1 + x2)* (-7*x + r*(-3 + 10*x2))) / (14. * kmp2 * kmp2);

			KI2uud2tA = (-2*r*(-1 + r*x)*rpow(fkmp,2)*(7*x + r*(-4 + 7*r*x - 10*x2)) - 
			   2*r*x*rpow(fp,2)*(1 - 2*r*x + r2)*(-3 + 7*r*x - 4*x2) + 
			   fkmp*fp*(r*x*(5 - 89*x2) + 28*x2 + 
			   28*rpow(r,4)*x2 - 28*rpow(r,3)*(x + 2*rpow(x,3)) + 
			   r2*(9 + 33*x2 + 70*rpow(x,4))) ) /  (14. * kmp2 * kmp2);
				   
			KI3uuu2tA =  ( fkmp*fp*r*(-1 + x2) * (fkmp* (7*x + r*(-4 + 7*r*x - 10*x2)) + 
				 fp*r*(3 - 7*r*x + 4*x2))  ) /  (14. * kmp2 * kmp2);    
				 
			KI3uuu3tA = -(   fkmp*fp*(-2*x + r*(-1 + 3*x2))*
			   (fkmp*(7*x + r*(-4 + 7*r*x - 10*x2)) + 
				 fp*r*(3 - 7*r*x + 4*x2))  ) /  (14. * kmp2*kmp2);  
     
// B+C     
			// KI2uudd1BpC =  fp*( fp * (1.-x2) + fkmp * r2 * (-1. + x2) / y2 ) / 2.;
			KI2uudd1BpC = 1/4. * (1.-x2)*(fp*fp + fkmp*fkmp*r2*r2/y2/y2) 
						+ fkmp*fp *r2 *(-1.+x2)/y2/2.; //Modificacion_Julio1
			
			KI2uudd2BpC = ( fp*fp*(-1 + 3*x2) + 
				 2 * fkmp * fp*r * (r + 2*x - 3*r*x2) / kmp2 + 
				 fkmp * fkmp * r2*
				  ( 2 - 4*r*x + r2*(-1 + 3*x2) )/(kmp2*kmp2)    )/4. ;
  
			KI3uuud2BpC = -(   fkmp*fp*(fkmp*(-2 + 3*r*x)*r2 - 
				 fp*(-1 + 3*r*x)*(1 - 2*r*x + r2))*(-1 + x2 )    )/ (2.* kmp2*kmp2)  ; 
  
  			KI3uuud3BpC =   (  fkmp*fp*( -(fp*(1 - 2*r*x + r2)*(1 - 3*x2 + r*x*(-3 + 5*x2))) + 
				   fkmp*r*(2*x + r*(2 - 6*x2 + r*x*(-3 + 5*x2))) ) ) / (2.* kmp2*kmp2)  ;
				    
			KI4uuuu2BpC =  (3*rpow(fkmp,2)*rpow(fp,2)*r2*rpow(-1 + x2,2)) / (16. *kmp2*kmp2 )  ;         
				   
			KI4uuuu3BpC = -(rpow(fkmp,2)*rpow(fp,2)*(-1 + x2)*
				  (2 + 3*r*(-4*x + r*(-1 + 5*x2))))  / (8. *kmp2*kmp2 )     ;
				   
			KI4uuuu4BpC =  (rpow(fkmp,2)*rpow(fp,2)*(-4 + 8*r*x*(3 - 5*x2) + 
				   12*x2 + r2*(3 - 30*x2 + 35*rpow(x,4))))  / (16. *kmp2*kmp2 )  ;
       
  // Bias    
			KPb1b2 = (r*(7*x + r*(3 - 10*x2))) /14. / kmp2;
		   
			KPb1bs2 = -(r*(-1 - 4*r*x + 2*r2 + 3*x2)*
				(-7*x + r*(-3 + 10*x2)))/42. / (kmp2*kmp2);
		   
			KPb22 =  -(rpow(pklkmp - pklp,2)*r2)/4. / (pklkmp * pklp ) ;
		   
			KPb2s2 = -(r2*(kmp2*(rpow(pklkmp,2) + rpow(pklp,2)) + 
				pklkmp*pklp*(1 + 4*r*x - 2*r2 - 3*x2)))/ (6. * pklkmp * pklp * kmp2);
		   
			KPs22 = r2/4. * (  -4*pklp/ (9. * pklkmp ) - 4*pklkmp / (9. * pklp) 
				+ 2. * rpow(-1./3. + (r-x)*(r-x)/ kmp2 ,2) ) ;
		   
		   
			KPb2theta = (r*(fkmp*(7*x + r*(-4 + 7*r*x - 10*x2)) + 
				fp*r*(3 - 7*r*x + 4*x2)))/14. / kmp2  ;
		   
		   
			KPbs2theta = (r*(-1 - 4*r*x + 2*r2 + 3*x2)* 
				( fkmp*(7*x + r*(-4 + 7*r*x - 10*x2)) + 
				fp*r*(3 - 7*r*x + 4*x2)) )  /42. / (kmp2*kmp2)  ;
       
            //
    
            P22dd_B +=   w*KP22dd*psl;
            P22du_B +=   w*KP22du*psl;
            P22uu_B +=   w*KP22uu*psl;
            
            I1udd1tA_B +=   w*KI1udd1tA*psl;
            I2uud1tA_B +=   w*KI2uud1tA*psl;
            I2uud2tA_B +=   w*KI2uud2tA*psl;
            I3uuu2tA_B +=   w*KI3uuu2tA*psl;
            I3uuu3tA_B +=   w*KI3uuu3tA*psl;
          
            I2uudd1BpC_B +=   w*KI2uudd1BpC*psl;
            I2uudd2BpC_B +=   w*KI2uudd2BpC*psl;
            I3uuud2BpC_B +=   w*KI3uuud2BpC*psl;
            I3uuud3BpC_B +=   w*KI3uuud3BpC*psl;
            I4uuuu2BpC_B +=   w*KI4uuuu2BpC*psl;
            I4uuuu3BpC_B +=   w*KI4uuuu3BpC*psl;
            I4uuuu4BpC_B +=   w*KI4uuuu4BpC*psl;

            
            Pb1b2_B     +=   w*KPb1b2*psl;
            Pb1bs2_B    +=   w*KPb1bs2*psl;
            Pb22_B      +=   w*KPb22*psl;
            Pb2s2_B     +=   w*KPb2s2*psl;
            Ps22_B      +=   w*KPs22*psl;
            Pb2theta_B  +=   w*KPb2theta*psl;
            Pbs2theta_B +=   w*KPbs2theta*psl;

        }
        
        P22dd_p   += dkk[i]*(P22dd_A*PSLA + P22dd_B*PSLB)/2.0;
        P22du_p   += dkk[i]*(P22du_A*PSLA + P22du_B*PSLB)/2.0;
        P22uu_p   += dkk[i]*(P22uu_A*PSLA + P22uu_B*PSLB)/2.0;
        
        I1udd1tA_p   += dkk[i]*(I1udd1tA_A*PSLA + I1udd1tA_B*PSLB)/2.0;
        I2uud1tA_p   += dkk[i]*(I2uud1tA_A*PSLA + I2uud1tA_B*PSLB)/2.0;
        I2uud2tA_p   += dkk[i]*(I2uud2tA_A*PSLA + I2uud2tA_B*PSLB)/2.0;
        I3uuu2tA_p   += dkk[i]*(I3uuu2tA_A*PSLA + I3uuu2tA_B*PSLB)/2.0;
        I3uuu3tA_p   += dkk[i]*(I3uuu3tA_A*PSLA + I3uuu3tA_B*PSLB)/2.0;
        
        I2uudd1BpC_p   += dkk[i]*(I2uudd1BpC_A*PSLA + I2uudd1BpC_B*PSLB)/2.0;
        I2uudd2BpC_p   += dkk[i]*(I2uudd2BpC_A*PSLA + I2uudd2BpC_B*PSLB)/2.0;
        I3uuud2BpC_p   += dkk[i]*(I3uuud2BpC_A*PSLA + I3uuud2BpC_B*PSLB)/2.0;
        I3uuud3BpC_p   += dkk[i]*(I3uuud3BpC_A*PSLA + I3uuud3BpC_B*PSLB)/2.0;
        I4uuuu2BpC_p   += dkk[i]*(I4uuuu2BpC_A*PSLA + I4uuuu2BpC_B*PSLB)/2.0;
        I4uuuu3BpC_p   += dkk[i]*(I4uuuu3BpC_A*PSLA + I4uuuu3BpC_B*PSLB)/2.0;
        I4uuuu4BpC_p   += dkk[i]*(I4uuuu4BpC_A*PSLA + I4uuuu4BpC_B*PSLB)/2.0;

        Pb1b2_p       += dkk[i]*(Pb1b2_A*PSLA + Pb1b2_B*PSLB)/2.0;
        Pb1bs2_p      += dkk[i]*(Pb1bs2_A*PSLA + Pb1bs2_B*PSLB)/2.0;
        Pb22_p        += dkk[i]*(Pb22_A*PSLA + Pb22_B*PSLB)/2.0;
        Pb2s2_p       += dkk[i]*(Pb2s2_A*PSLA + Pb2s2_B*PSLB)/2.0;
        Ps22_p        += dkk[i]*(Ps22_A*PSLA + Ps22_B*PSLB)/2.0;
        Pb2theta_p    += dkk[i]*(Pb2theta_A*PSLA + Pb2theta_B*PSLB)/2.0;
        Pbs2theta_p   += dkk[i]*(Pbs2theta_A*PSLA + Pbs2theta_B*PSLB)/2.0;

     
        P22dd_A =   P22dd_B;      P22dd_B = 0.0;
        P22du_A =   P22du_B;      P22du_B = 0.0;
        P22uu_A =   P22uu_B;      P22uu_B = 0.0;
        
        I1udd1tA_A =   I1udd1tA_B;      I1udd1tA_B = 0.0;
        I2uud1tA_A =   I2uud1tA_B;      I2uud1tA_B = 0.0;
        I2uud2tA_A =   I2uud2tA_B;      I2uud2tA_B = 0.0;
        I3uuu2tA_A =   I3uuu2tA_B;      I3uuu2tA_B = 0.0;
        I3uuu3tA_A =   I3uuu3tA_B;      I3uuu3tA_B = 0.0;
         
        I2uudd1BpC_A =   I2uudd1BpC_B;      I2uudd1BpC_B = 0.0;
        I2uudd2BpC_A =   I2uudd2BpC_B;      I2uudd2BpC_B = 0.0;
        I3uuud2BpC_A =   I3uuud2BpC_B;      I3uuud2BpC_B = 0.0;
        I3uuud3BpC_A =   I3uuud3BpC_B;      I3uuud3BpC_B = 0.0;
        I4uuuu2BpC_A =   I4uuuu2BpC_B;      I4uuuu2BpC_B = 0.0;
        I4uuuu3BpC_A =   I4uuuu3BpC_B;      I4uuuu3BpC_B = 0.0;
        I4uuuu4BpC_A =   I4uuuu4BpC_B;      I4uuuu4BpC_B = 0.0;
        
        Pb1b2_A     =   Pb1b2_B;        Pb1b2_B     = 0.0;
        Pb1bs2_A    =   Pb1bs2_B;       Pb1bs2_B    = 0.0;
        Pb22_A      =   Pb22_B;         Pb22_B      = 0.0;
        Pb2s2_A     =   Pb2s2_B;        Pb2s2_B     = 0.0;
        Ps22_A      =   Ps22_B;         Ps22_B      = 0.0;
        Pb2theta_A  =   Pb2theta_B;     Pb2theta_B  = 0.0;
        Pbs2theta_A =   Pbs2theta_B;    Pbs2theta_B = 0.0;
       
 
        
        
        PSLA = PSLB;
    }
		
	P22dd_p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	P22du_p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	P22uu_p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
    
	I1udd1tA_p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I2uud1tA_p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I2uud2tA_p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I3uuu2tA_p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I3uuu3tA_p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
        
    I2uudd1BpC_p    *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I2uudd2BpC_p    *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I3uuud2BpC_p    *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I3uuud3BpC_p    *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I4uuuu2BpC_p    *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I4uuuu3BpC_p    *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I4uuuu4BpC_p    *= 2.0*(rpow(ki,3)/FOURPI2)/ki;

	Pb1b2_p        *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	Pb1bs2_p       *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	Pb22_p         *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	Pb2s2_p        *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	Ps22_p         *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	Pb2theta_p     *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	Pbs2theta_p    *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
    
    free_dvector(wwGL,1,Nx);
    free_dvector(xxGL,1,Nx);

//  R functions



    nGL=10;
    xGL=dvector(1,nGL);
    wGL=dvector(1,nGL);
    gauleg(-1.0,1.0,xGL,wGL,nGL);
    
    //~ fk = 1.0; // This is f(k)/f0 we need to interpolate f(k) 
	fk = Interpolation_nr(ki, kPS, fkT, nPSLT, fkT2);
	fk /= gd.f0;    
    for (i=2; i<cmd.nquadSteps; i++) {
        r = kk[i]/ki;
        r2= r*r;
        //~ psl = psInterpolation_nr(kk[i], kPS, pPS, nPSLT);
        psl = Interpolation_nr(kk[i], kPKL, pPKL, nPKLT, pPKL2); ;
        fp = Interpolation_nr(kk[i], kPS, fkT, nPSLT, fkT2);
		fp /= gd.f0; 
        //~ for (j=1; j<=nGL(pGL); j++) { //Modificacion
        for (j=1; j<=nGL; j++) {
            x = xGL[j];
            w = wGL[j];
            x2 =x*x;     
            kmp2=1.0 + r2 - 2.0 * r * x;
              
			KP13dd = (-21*x2 + 6*r*x*(3 + 4*x2) + 
				r2*(10 - 59*x2 + 28*x2*x2))/21. / kmp2;

			KP13du = (fp*r*(-1 + x2)*(-9*x + r*(-1 + 10*x2)) + 
				fk*(-21*x2 - 9*r2*r*x*(-1 + x2) + 
				6*r*x*(3 + 4*x2) + 
				r2*(1 - 50*x2 + 28*x2*x2)))/21. / kmp2;

			KP13uu = -(fk*(-2*fp*r*(-1 + x2)*(-9*x + r*(-1 + 10*x2)) + 
				fk*(21*x2 + 18*r2*r*x*(-1 + x2) - 
				6*r*x*(3 + 4*x2) + 
				r2*(8 + 41*x2 - 28*x2*x2))))/21. / kmp2;
                 
			Ksigma32PSL = ( 5.0* r2 * (7. - 2*r2 + 4*r*x + 6*(-2 + r2)*x2 - 
				12*r*x2*x + 9*x2*x2)) / (24.0 * kmp2 ) ;
            
			KI1udd1a = 2*r*(fp*x*(5./7. - ((1/r + r)*x)/2. + 2*x*x/7.) + 
				((-1 + r*x)*(-3*(fk + fp)*r + 7*(fk + fp*r2)*x - 
				4*(fk + fp)*r*x2)) / (14.*kmp2)  ) ;

			KI2uud1a = (fp*r*(-1 + x2)*(3*(fk + fp)*r - 7*(fk + fp*r2)*x + 
				4*(fk + fp)*r*x2)) / (14.*kmp2);

			KI2uud2a = (fk*fp*x*(10*r - 7*(1 + r2)*x + 4*r*x2))/7. - 
				((3*(fk + fp)*r - 7*(fk + fp*r2)*x + 
				4*(fk + fp)*r*x2)*(2*fk*r*(-1 + r*x) + fp*(-2*x + r*(-1 + 3*x2))))
				/(14.*kmp2) ; 

			KI3uuu2a =  fk * KI2uud1a ;

			KI3uuu3a =  (fk*fp*r*(3*(fk + fp) - 7*(fk/r + fp*r)*x + 4*(fk + fp)*x2) 
						* (r + 2*x - 3*r*x2))/(14.*kmp2) ;            
            
            P13dd_B +=   w*KP13dd*psl;
            P13du_B +=   w*KP13du*psl;
            P13uu_B +=   w*KP13uu*psl;
            
            sigma32PSL_B +=   w*Ksigma32PSL*psl;
            
            I1udd1a_B +=   w*KI1udd1a*psl;
            I2uud1a_B +=   w*KI2uud1a*psl;
            I2uud2a_B +=   w*KI2uud2a*psl;
            I3uuu2a_B +=   w*KI3uuu2a*psl;
            I3uuu3a_B +=   w*KI3uuu3a*psl;           
            
        }
        
        
        
        P13dd_p   += dkk[i]*(P13dd_A + P13dd_B) /  (2.0*ki);
        P13du_p   += dkk[i]*(P13du_A + P13du_B) /  (2.0*ki);
        P13uu_p   += dkk[i]*(P13uu_A + P13uu_B) /  (2.0*ki);
        sigma32PSL_p   += dkk[i]*(sigma32PSL_A + sigma32PSL_B) /  (2.0*ki);
        
        I1udd1a_p   += dkk[i]*(I1udd1a_A + I1udd1a_B) /  (2.0*ki);
        I2uud1a_p   += dkk[i]*(I2uud1a_A + I2uud1a_B) /  (2.0*ki);
        I2uud2a_p   += dkk[i]*(I2uud2a_A + I2uud2a_B) /  (2.0*ki);
        I3uuu2a_p   += dkk[i]*(I3uuu2a_A + I3uuu2a_B) /  (2.0*ki);
        I3uuu3a_p   += dkk[i]*(I3uuu3a_A + I3uuu3a_B) /  (2.0*ki);
        
              
        P13dd_A =   P13dd_B;      P13dd_B = 0.0;
        P13du_A =   P13du_B;      P13du_B = 0.0;
        P13uu_A =   P13uu_B;      P13uu_B = 0.0;
        sigma32PSL_A =   sigma32PSL_B;      sigma32PSL_B = 0.0;
        
        I1udd1a_A =   I1udd1a_B;      I1udd1a_B = 0.0;
        I2uud1a_A =   I2uud1a_B;      I2uud1a_B = 0.0;
        I2uud2a_A =   I2uud2a_B;      I2uud2a_B = 0.0;
        I3uuu2a_A =   I3uuu2a_B;      I3uuu2a_B = 0.0;
        I3uuu3a_A =   I3uuu3a_B;      I3uuu3a_B = 0.0;
      
    }
    //~ pkl_k = psInterpolation_nr(ki, kPS, pPS, nPSLT);
    pkl_k = Interpolation_nr(ki, kPKL, pPKL, nPKLT, pPKL2);
    P13dd_p *= (rpow(ki,3.0)/FOURPI2)*pkl_k;
    P13du_p *= (rpow(ki,3.0)/FOURPI2)*pkl_k;
    P13uu_p *= (rpow(ki,3.0)/FOURPI2)*pkl_k;
    sigma32PSL_p *= (rpow(ki,3.0)/FOURPI2)*pkl_k;
    I1udd1a_p *= (rpow(ki,3.0)/FOURPI2)*pkl_k;
    I2uud1a_p *= (rpow(ki,3.0)/FOURPI2)*pkl_k;
    I2uud2a_p *= (rpow(ki,3.0)/FOURPI2)*pkl_k;
    I3uuu2a_p *= (rpow(ki,3.0)/FOURPI2)*pkl_k;
    I3uuu3a_p *= (rpow(ki,3.0)/FOURPI2)*pkl_k;



    //~ etaQRs(QRstmp) = eta;
    //
	kFs(QRstmp)    = ki;
		
	P22dd(  QRstmp)      = P22dd_p;
	P22du(  QRstmp)      = P22du_p;
	P22uu(  QRstmp)      = P22uu_p;
	// A 
	I1udd1A(  QRstmp)      = I1udd1tA_p + 2.0*I1udd1a_p;
	I2uud1A(  QRstmp)      = I2uud1tA_p + 2.0*I2uud1a_p;
	I2uud2A(  QRstmp)      = I2uud2tA_p + 2.0*I2uud2a_p;
	I3uuu2A(  QRstmp)      = I3uuu2tA_p + 2.0*I3uuu2a_p;
	I3uuu3A(  QRstmp)      = I3uuu3tA_p + 2.0*I3uuu3a_p;
	//  B plus C   
	I2uudd1BpC(  QRstmp)   = I2uudd1BpC_p  - ki*ki*gd.sigma2v*pkl_k;
	I2uudd2BpC(  QRstmp)   = I2uudd2BpC_p;
	I3uuud2BpC(  QRstmp)   = I3uuud2BpC_p - 2.0*ki*ki*gd.sigma2v*fk*pkl_k;
	I3uuud3BpC(  QRstmp)   = I3uuud3BpC_p;
	I4uuuu2BpC(  QRstmp)   = I4uuuu2BpC_p;
	I4uuuu3BpC(  QRstmp)   = I4uuuu3BpC_p - ki*ki*gd.sigma2v*fk*fk*pkl_k;
	I4uuuu4BpC(  QRstmp)   = I4uuuu4BpC_p;
	//  Bias
	Pb1b2(    QRstmp) = Pb1b2_p;
	Pb1bs2(   QRstmp) = Pb1bs2_p;
	Pb22(     QRstmp) = Pb22_p;
	Pb2s2(    QRstmp) = Pb2s2_p;
	Ps22(     QRstmp) = Ps22_p;
	Pb2theta( QRstmp) = Pb2theta_p;
	Pbs2theta(QRstmp) = Pbs2theta_p;
	//
	P13dd(  QRstmp)      = P13dd_p;
	P13du(  QRstmp)      = P13du_p;
	P13uu(  QRstmp)      = P13uu_p;
	
	sigma32PSL(QRstmp)   = sigma32PSL_p;
    

    free_dvector(dkk,1,cmd.nquadSteps);
    free_dvector(kk,1,cmd.nquadSteps);
    
    return *QRstmp;
}
*/







//Modificacion-LCDMfk
// Q and R functions quadrature
// kk is the inner integration moment: 
// kk = k * r, so usual notation kk = p
//~ global_kFs ki_functions(real eta, real ki)
global_kFs ki_functions(real ki, double kPKL[], double pPKL[], int nPKLT, double pPKL2[])
{
    int i, j;
    real pkl_k;
     
    real PSLA, PSLB, psl;
	real fk, fp, fkmp, pklp, pklkmp;
    real rmin, rmax;
    real r, deltar, r2, y, y2;
    real mumin, mumax;
    real x, w, x2;
    real psl1;
    int Nx, nGL;
    real ypi, dk;
    real *xxGL, *wwGL, *xGL, *wGL;
    real kmin, kmax, pmin, pmax;
	real AngleEvQ, S2evQ, G2evQ, F2evQ, G2evR, F2evR; 
	real Gamma2evR, Gamma2fevR, C3Gamma3, C3Gamma3f, G3K, F3K, AngleEvR; 
    real A, ApOverf0, CFD3, CFD3p;
    
    
    if (_KERNELS_LCDMfk_==1) {
		A=KA_LCDM; ApOverf0 = KAp_LCDM/gd.f0;
		CFD3 = KR1_LCDM;  CFD3p = KR1p_LCDM;
	} else {		
		A=1; ApOverf0 = 0;
		CFD3 = 1;  CFD3p = 1;
	}
    
	real P22dd_p =0.0, P22dd_A = 0.0, P22dd_B = 0.0;
    real P22du_p =0.0, P22du_A = 0.0, P22du_B = 0.0;
    real P22uu_p =0.0, P22uu_A = 0.0, P22uu_B = 0.0;
    
    real I1udd1tA_p =0.0, I1udd1tA_A = 0.0, I1udd1tA_B = 0.0;
    real I2uud1tA_p =0.0, I2uud1tA_A = 0.0, I2uud1tA_B = 0.0;
    real I2uud2tA_p =0.0, I2uud2tA_A = 0.0, I2uud2tA_B = 0.0;
    real I3uuu2tA_p =0.0, I3uuu2tA_A = 0.0, I3uuu2tA_B = 0.0;
    real I3uuu3tA_p =0.0, I3uuu3tA_A = 0.0, I3uuu3tA_B = 0.0;
    
    real I2uudd1BpC_p =0.0, I2uudd1BpC_A =0.0, I2uudd1BpC_B =0.0;
    real I2uudd2BpC_p =0.0, I2uudd2BpC_A =0.0, I2uudd2BpC_B =0.0 ;
    real I3uuud2BpC_p =0.0, I3uuud2BpC_A =0.0, I3uuud2BpC_B =0.0 ;
    real I3uuud3BpC_p =0.0, I3uuud3BpC_A =0.0, I3uuud3BpC_B =0.0 ;
    real I4uuuu2BpC_p =0.0, I4uuuu2BpC_A =0.0, I4uuuu2BpC_B =0.0 ;
    real I4uuuu3BpC_p =0.0, I4uuuu3BpC_A =0.0, I4uuuu3BpC_B =0.0 ;
    real I4uuuu4BpC_p =0.0, I4uuuu4BpC_A =0.0, I4uuuu4BpC_B =0.0 ;
      
    real  Pb1b2_p=0.0,     Pb1b2_A =0.0,     Pb1b2_B =0.0;
    real  Pb1bs2_p=0.0,    Pb1bs2_A =0.0,    Pb1bs2_B =0.0; 
    real  Pb22_p=0.0,      Pb22_A =0.0,      Pb22_B =0.0; 
    real  Pb2s2_p=0.0,     Pb2s2_A =0.0,     Pb2s2_B =0.0; 
    real  Ps22_p=0.0,      Ps22_A =0.0,      Ps22_B =0.0; 
    real  Pb2theta_p=0.0,  Pb2theta_A =0.0,  Pb2theta_B =0.0; 
    real  Pbs2theta_p=0.0, Pbs2theta_A =0.0, Pbs2theta_B =0.0; 

    real KP22dd, KP22du, KP22uu;
    real KI1udd1tA, KI2uud1tA, KI2uud2tA, KI3uuu2tA, KI3uuu3tA;
    real KI2uudd1BpC, KI2uudd2BpC, KI3uuud2BpC, KI3uuud3BpC; 
    real KI4uuuu2BpC, KI4uuuu3BpC, KI4uuuu4BpC;
    real KPb1b2, KPb1bs2, KPb22, KPb2s2, KPs22, KPb2theta, KPbs2theta;
    
	real P13dd_p =0.0, P13dd_A = 0.0, P13dd_B = 0.0;
    real P13du_p =0.0, P13du_A = 0.0, P13du_B = 0.0;
    real P13uu_p =0.0, P13uu_A = 0.0, P13uu_B = 0.0;

    real sigma32PSL_p =0.0, sigma32PSL_A = 0.0, sigma32PSL_B = 0.0;
    
    real I1udd1a_p =0.0, I1udd1a_A = 0.0, I1udd1a_B = 0.0;
    real I2uud1a_p =0.0, I2uud1a_A = 0.0, I2uud1a_B = 0.0;
    real I2uud2a_p =0.0, I2uud2a_A = 0.0, I2uud2a_B = 0.0;
    real I3uuu2a_p =0.0, I3uuu2a_A = 0.0, I3uuu2a_B = 0.0;
    real I3uuu3a_p =0.0, I3uuu3a_A = 0.0, I3uuu3a_B = 0.0;
    
    real KP13dd, KP13du, KP13uu, Ksigma32PSL;
    real KI1udd1a, KI2uud1a, KI2uud2a, KI3uuu2a, KI3uuu3a;
    
    real *kk, *dkk;
    //
    pointPSTableptr p;
    //
    global_kFs_ptr QRstmp;
    
    QRstmp = (global_kFs_ptr) allocate(1 * sizeof(global_kFs));
    
    
 
    kmin = kPS[1];
    kmax = kPS[nPSLT];
    pmin = MAX(kmin,0.01*cmd.kmin);
    pmax = MIN(kmax,16.0*cmd.kmax);
   
    
    dk = (rlog10(pmax) - rlog10(pmin))/((real)(cmd.nquadSteps - 1));
    kk=dvector(1,cmd.nquadSteps);
    dkk=dvector(1,cmd.nquadSteps);
    kk[1] = rpow(10.0,rlog10(pmin));
    for (i=2; i<cmd.nquadSteps; i++) {
        ypi = rlog10(pmin) + dk*((real)(i - 1));
        kk[i] = rpow(10.0,ypi);
        dkk[i] = (kk[i]-kk[i-1]);
    }

// Q functions

    PSLA = 0.0;
    rmax = kmax/ki;
    rmin = rmin/ki;
    Nx=10;
    xxGL=dvector(1,Nx);
    wwGL=dvector(1,Nx);
    
	for (i=2; i<cmd.nquadSteps; i++) {
		r = kk[i]/ki;
        r2= r*r;
		PSLB = Interpolation_nr(kk[i], kPKL, pPKL, nPKLT, pPKL2);
		
		pklp=PSLB;
		fp = Interpolation_nr(kk[i], kPS, fkT, nPSLT, fkT2);
		fp /= gd.f0;
		mumin = MAX( -1.0, (1.0 + rsqr(r) - rsqr(rmax)) / (2.0*r)  );
		mumax = MIN(   1.0, (1.0  + rsqr(r) - rsqr(rmin)) / (2.0*r)  );
        
		if (r>=0.5)
			mumax = 0.5/r;
			gauleg(mumin,mumax,xxGL,wwGL,Nx);
        for (j=1; j<=Nx; j++) {
            x = xxGL[j];
            w = wwGL[j];
            x2= x*x;
            
            y2=1.0 + r2 - 2.0 * r * x;
            y = rsqrt(y2);
            psl = Interpolation_nr(ki * y, kPKL, pPKL, nPKLT, pPKL2); 
            pklkmp=psl;
			fkmp = Interpolation_nr(ki * y, kPS, fkT, nPSLT, fkT2);
			fkmp /= gd.f0;
			 
			AngleEvQ = (x - r)/y;
			S2evQ = AngleEvQ *AngleEvQ - 1./3.;
			F2evQ = 1./2. + 3./14. * A + (1./2. - 3./14. * A) * AngleEvQ*AngleEvQ + 
					AngleEvQ / 2. * (y/r + r/y);
			G2evQ = 3./14. * A * (fp + fkmp) + 3./14. * ApOverf0 
					+ (1./2. * (fp + fkmp) - 3./14. *  A *  (fp + fkmp) - 
					3./14. * ApOverf0)*AngleEvQ*AngleEvQ  
					+ AngleEvQ/2. *  (fkmp * y/r + fp * r/y);
			

			KP22dd = 2*r2*F2evQ*F2evQ;
			KP22du = 2*r2*F2evQ*G2evQ;
			KP22uu = 2*r2*G2evQ*G2evQ;
// A


			KI1udd1tA = 2.* (fp * r * x + fkmp * r2 * (1. - r * x)/y2 ) * F2evQ ;

			KI2uud1tA = - fp*fkmp * r2 *(1. - x2)/y2 * F2evQ;

			KI2uud2tA = 2.* (fp* r * x + fkmp *  r2*(1. - r * x)/y2) * G2evQ + 
						fp * fkmp * ( r2*(1. - 3.* x2) + 2.* r * x )/y2 * F2evQ;
						
			KI3uuu2tA = fp * fkmp *  r2 * (x2 - 1.)/y2 * G2evQ ;

			KI3uuu3tA = fp * fkmp * ( r2*(1. - 3.* x2) + 2.* r * x )/y2 * G2evQ;


     
// B+C     
			//~ KI2uudd1BpC =  fp*( fp * (1.-x2) + fkmp * r2 * (-1. + x2) / y2 ) / 2.;
			KI2uudd1BpC = 1/4. * (1.-x2)*(fp*fp + fkmp*fkmp*r2*r2/y2/y2) 
						+ fkmp*fp *r2 *(-1.+x2)/y2/2.; //Modificacion_Julio1
			
			KI2uudd2BpC = ( fp*fp*(-1. + 3.*x2) + 
				 2 * fkmp * fp*r * (r + 2*x - 3*r*x2) / y2 + 
				 fkmp * fkmp * r2*
				  ( 2 - 4*r*x + r2*(-1 + 3*x2) )/(y2*y2)    )/4. ;
  
			KI3uuud2BpC = -(   fkmp*fp*(fkmp*(-2 + 3*r*x)*r2 - 
				 fp*(-1 + 3*r*x)*(1 - 2*r*x + r2))*(-1 + x2 )    )/ (2.* y2*y2)  ; 
  
  			KI3uuud3BpC =   (  fkmp*fp*( -(fp*(1 - 2*r*x + r2)*(1 - 3*x2 + r*x*(-3 + 5*x2))) + 
				   fkmp*r*(2*x + r*(2 - 6*x2 + r*x*(-3 + 5*x2))) ) ) / (2.* y2*y2)  ;
				    
			KI4uuuu2BpC =  (3*rpow(fkmp,2)*rpow(fp,2)*r2*rpow(-1 + x2,2)) / (16. *y2*y2 )  ;         
				   
			KI4uuuu3BpC = -(rpow(fkmp,2)*rpow(fp,2)*(-1 + x2)*
				  (2 + 3*r*(-4*x + r*(-1 + 5*x2))))  / (8. *y2*y2 );
				   
			KI4uuuu4BpC =  (rpow(fkmp,2)*rpow(fp,2)*(-4 + 8*r*x*(3 - 5*x2) + 
				   12*x2 + r2*(3 - 30*x2 + 35*rpow(x,4))))  / (16. *y2*y2 );			
			
			       
  // Bias    
       
			KPb1b2  = r2 * F2evQ;
			KPb1bs2 = r2 * F2evQ*S2evQ;
			KPb22   = 1./2. * r2 * (1./2. * (1. - pklp/pklkmp) 
							+ 1./2. * (1. - pklkmp/pklp));
			KPb2s2  = 1./2. * r2 * (  1./2.* (S2evQ - 2./3. * pklp/pklkmp) 
							        + 1./2.* (S2evQ - 2./3. * pklkmp/pklp));
			KPs22   = 1./2. * r2 * (1./2. * (S2evQ * S2evQ - 4./9. * pklp/pklkmp) 
							+ 1./2.* (S2evQ * S2evQ - 4./9. * pklkmp/pklp));
			KPb2theta  = r2 * G2evQ;
			KPbs2theta = r2 * S2evQ * G2evQ;      
       
            //
    
            P22dd_B +=   w*KP22dd*psl;
            P22du_B +=   w*KP22du*psl;
            P22uu_B +=   w*KP22uu*psl;
            
            I1udd1tA_B +=   w*KI1udd1tA*psl;
            I2uud1tA_B +=   w*KI2uud1tA*psl;
            I2uud2tA_B +=   w*KI2uud2tA*psl;
            I3uuu2tA_B +=   w*KI3uuu2tA*psl;
            I3uuu3tA_B +=   w*KI3uuu3tA*psl;
          
            I2uudd1BpC_B +=   w*KI2uudd1BpC*psl;
            I2uudd2BpC_B +=   w*KI2uudd2BpC*psl;
            I3uuud2BpC_B +=   w*KI3uuud2BpC*psl;
            I3uuud3BpC_B +=   w*KI3uuud3BpC*psl;
            I4uuuu2BpC_B +=   w*KI4uuuu2BpC*psl;
            I4uuuu3BpC_B +=   w*KI4uuuu3BpC*psl;
            I4uuuu4BpC_B +=   w*KI4uuuu4BpC*psl;

            
            Pb1b2_B     +=   w*KPb1b2*psl;
            Pb1bs2_B    +=   w*KPb1bs2*psl;
            Pb22_B      +=   w*KPb22*psl;
            Pb2s2_B     +=   w*KPb2s2*psl;
            Ps22_B      +=   w*KPs22*psl;
            Pb2theta_B  +=   w*KPb2theta*psl;
            Pbs2theta_B +=   w*KPbs2theta*psl;

        }
        
        P22dd_p   += dkk[i]*(P22dd_A*PSLA + P22dd_B*PSLB)/2.0;
        P22du_p   += dkk[i]*(P22du_A*PSLA + P22du_B*PSLB)/2.0;
        P22uu_p   += dkk[i]*(P22uu_A*PSLA + P22uu_B*PSLB)/2.0;
        
        I1udd1tA_p   += dkk[i]*(I1udd1tA_A*PSLA + I1udd1tA_B*PSLB)/2.0;
        I2uud1tA_p   += dkk[i]*(I2uud1tA_A*PSLA + I2uud1tA_B*PSLB)/2.0;
        I2uud2tA_p   += dkk[i]*(I2uud2tA_A*PSLA + I2uud2tA_B*PSLB)/2.0;
        I3uuu2tA_p   += dkk[i]*(I3uuu2tA_A*PSLA + I3uuu2tA_B*PSLB)/2.0;
        I3uuu3tA_p   += dkk[i]*(I3uuu3tA_A*PSLA + I3uuu3tA_B*PSLB)/2.0;
        
        I2uudd1BpC_p   += dkk[i]*(I2uudd1BpC_A*PSLA + I2uudd1BpC_B*PSLB)/2.0;
        I2uudd2BpC_p   += dkk[i]*(I2uudd2BpC_A*PSLA + I2uudd2BpC_B*PSLB)/2.0;
        I3uuud2BpC_p   += dkk[i]*(I3uuud2BpC_A*PSLA + I3uuud2BpC_B*PSLB)/2.0;
        I3uuud3BpC_p   += dkk[i]*(I3uuud3BpC_A*PSLA + I3uuud3BpC_B*PSLB)/2.0;
        I4uuuu2BpC_p   += dkk[i]*(I4uuuu2BpC_A*PSLA + I4uuuu2BpC_B*PSLB)/2.0;
        I4uuuu3BpC_p   += dkk[i]*(I4uuuu3BpC_A*PSLA + I4uuuu3BpC_B*PSLB)/2.0;
        I4uuuu4BpC_p   += dkk[i]*(I4uuuu4BpC_A*PSLA + I4uuuu4BpC_B*PSLB)/2.0;

        Pb1b2_p       += dkk[i]*(Pb1b2_A    *PSLA + Pb1b2_B    *PSLB)/2.0;
        Pb1bs2_p      += dkk[i]*(Pb1bs2_A   *PSLA + Pb1bs2_B   *PSLB)/2.0;
        Pb22_p        += dkk[i]*(Pb22_A     *PSLA + Pb22_B     *PSLB)/2.0;
        Pb2s2_p       += dkk[i]*(Pb2s2_A    *PSLA + Pb2s2_B    *PSLB)/2.0;
        Ps22_p        += dkk[i]*(Ps22_A     *PSLA + Ps22_B     *PSLB)/2.0;
        Pb2theta_p    += dkk[i]*(Pb2theta_A *PSLA + Pb2theta_B *PSLB)/2.0;
        Pbs2theta_p   += dkk[i]*(Pbs2theta_A*PSLA + Pbs2theta_B*PSLB)/2.0;

     
        P22dd_A =   P22dd_B;      P22dd_B = 0.0;
        P22du_A =   P22du_B;      P22du_B = 0.0;
        P22uu_A =   P22uu_B;      P22uu_B = 0.0;
        
        I1udd1tA_A   =   I1udd1tA_B;        I1udd1tA_B   = 0.0;
        I2uud1tA_A   =   I2uud1tA_B;        I2uud1tA_B   = 0.0;
        I2uud2tA_A   =   I2uud2tA_B;        I2uud2tA_B   = 0.0;
        I3uuu2tA_A   =   I3uuu2tA_B;        I3uuu2tA_B   = 0.0;
        I3uuu3tA_A   =   I3uuu3tA_B;        I3uuu3tA_B   = 0.0;
         
        I2uudd1BpC_A =   I2uudd1BpC_B;      I2uudd1BpC_B = 0.0;
        I2uudd2BpC_A =   I2uudd2BpC_B;      I2uudd2BpC_B = 0.0;
        I3uuud2BpC_A =   I3uuud2BpC_B;      I3uuud2BpC_B = 0.0;
        I3uuud3BpC_A =   I3uuud3BpC_B;      I3uuud3BpC_B = 0.0;
        I4uuuu2BpC_A =   I4uuuu2BpC_B;      I4uuuu2BpC_B = 0.0;
        I4uuuu3BpC_A =   I4uuuu3BpC_B;      I4uuuu3BpC_B = 0.0;
        I4uuuu4BpC_A =   I4uuuu4BpC_B;      I4uuuu4BpC_B = 0.0;
        
        Pb1b2_A     =   Pb1b2_B;        Pb1b2_B     = 0.0;
        Pb1bs2_A    =   Pb1bs2_B;       Pb1bs2_B    = 0.0;
        Pb22_A      =   Pb22_B;         Pb22_B      = 0.0;
        Pb2s2_A     =   Pb2s2_B;        Pb2s2_B     = 0.0;
        Ps22_A      =   Ps22_B;         Ps22_B      = 0.0;
        Pb2theta_A  =   Pb2theta_B;     Pb2theta_B  = 0.0;
        Pbs2theta_A =   Pbs2theta_B;    Pbs2theta_B = 0.0;
       
 
        
        
        PSLA = PSLB;
    }
		
	P22dd_p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	P22du_p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	P22uu_p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
    
	I1udd1tA_p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I2uud1tA_p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I2uud2tA_p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I3uuu2tA_p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I3uuu3tA_p   *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
        
    I2uudd1BpC_p    *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I2uudd2BpC_p    *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I3uuud2BpC_p    *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I3uuud3BpC_p    *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I4uuuu2BpC_p    *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I4uuuu3BpC_p    *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	I4uuuu4BpC_p    *= 2.0*(rpow(ki,3)/FOURPI2)/ki;

	Pb1b2_p        *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	Pb1bs2_p       *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	Pb22_p         *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	Pb2s2_p        *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	Ps22_p         *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	Pb2theta_p     *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
	Pbs2theta_p    *= 2.0*(rpow(ki,3)/FOURPI2)/ki;
    
    free_dvector(wwGL,1,Nx);
    free_dvector(xxGL,1,Nx);

//  R functions



    nGL=10;
    xGL=dvector(1,nGL);
    wGL=dvector(1,nGL);
    gauleg(-1.0,1.0,xGL,wGL,nGL);
    
    //~ fk = 1.0; // This is f(k)/f0 we need to interpolate f(k) 
	fk = Interpolation_nr(ki, kPS, fkT, nPSLT, fkT2);
	fk /= gd.f0;    
    for (i=2; i<cmd.nquadSteps; i++) {
        r = kk[i]/ki;
        r2= r*r;
        //~ psl = psInterpolation_nr(kk[i], kPS, pPS, nPSLT);
        psl = Interpolation_nr(kk[i], kPKL, pPKL, nPKLT, pPKL2); ;
        fp = Interpolation_nr(kk[i], kPS, fkT, nPSLT, fkT2);
		fp /= gd.f0; 
        for (j=1; j<=nGL; j++) {
            x = xGL[j];
            w = wGL[j];
            x2 =x*x;     
            y2=1.0 + r2 - 2.0 * r * x;
            
            
			Gamma2evR  = A *(1. - x2);
			Gamma2fevR = A *(1. - x2)*(fk + fp)/2. + 1./2. * ApOverf0 *(1 - x2);

			C3Gamma3  = 2.*5./21. * CFD3  *(1 - x2)*(1 - x2)/y2;
			C3Gamma3f = 2.*5./21. * CFD3p *(1 - x2)*(1 - x2)/y2 *(fk + 2 * fp)/3.;

			G3K = C3Gamma3f/ 2. + (2 * Gamma2fevR * x)/(7. * r) - (fk  * x2)/(6 * r2) 
				+ fp * Gamma2evR*(1 - r * x)/(7 * y2) 
				- 1./7.*(fp * Gamma2evR + 2 *Gamma2fevR) * (1. - x2)/y2;
			F3K = C3Gamma3/6. - x2/(6 * r2) + (Gamma2evR * x *(1 - r * x))/(7. *r *y2);
			   
			AngleEvR = -x;
			F2evR = 1./2. + 3./14. * A + (1./2. - 3./14. * A) * AngleEvR*AngleEvR 
				+ AngleEvR/2. *(1./r + r);

			G2evR = 3./14.* A *(fp + fk) + 3./14.* ApOverf0 
				+ ((fp + fk)/2. - 3./14.* A *(fp + fk) - 3./14.* ApOverf0)*AngleEvR*AngleEvR 
				+ AngleEvR/2. * (fk/r + fp * r);             
            
			KP13dd = 6.* r2 * F3K;
			KP13du = 3.* r2 * G3K + 3.* r2 * F3K * fk;
			KP13uu = 6.* r2 * G3K * fk;

			Ksigma32PSL = ( 5.0* r2 * (7. - 2*r2 + 4*r*x + 6*(-2 + r2)*x2 - 
				12*r*x2*x + 9*x2*x2)) / (24.0 * y2 ) ;
				
			KI1udd1a = 2.*r2*(1 - r * x)/y2 *G2evR + 2*fp *r* x* F2evR ;

			KI2uud1a = -fp * r2 * (1 - x2)/y2*G2evR ;

			KI2uud2a = ( (r2 *(1 - 3.*x2) + 2.* r* x) /y2*fp + 
				  fk*2.*r2*(1 - r * x)/y2)* G2evR + 2*x*r*fp*fk*F2evR;


			KI3uuu2a = fk * KI2uud1a;

			KI3uuu3a = (r2 *(1 - 3.*x2) + 2.* r* x)/y2 * fp*fk * G2evR;              

    
            P13dd_B +=   w*KP13dd*psl;
            P13du_B +=   w*KP13du*psl;
            P13uu_B +=   w*KP13uu*psl;
            
            sigma32PSL_B +=   w*Ksigma32PSL*psl;
            
            I1udd1a_B +=   w*KI1udd1a*psl;
            I2uud1a_B +=   w*KI2uud1a*psl;
            I2uud2a_B +=   w*KI2uud2a*psl;
            I3uuu2a_B +=   w*KI3uuu2a*psl;
            I3uuu3a_B +=   w*KI3uuu3a*psl;           
            
        }
        
        
        
        P13dd_p   += dkk[i]*(P13dd_A + P13dd_B) /  (2.0*ki);
        P13du_p   += dkk[i]*(P13du_A + P13du_B) /  (2.0*ki);
        P13uu_p   += dkk[i]*(P13uu_A + P13uu_B) /  (2.0*ki);
        
        sigma32PSL_p   += dkk[i]*(sigma32PSL_A + sigma32PSL_B) /  (2.0*ki);
        
        I1udd1a_p   += dkk[i]*(I1udd1a_A + I1udd1a_B) /  (2.0*ki);
        I2uud1a_p   += dkk[i]*(I2uud1a_A + I2uud1a_B) /  (2.0*ki);
        I2uud2a_p   += dkk[i]*(I2uud2a_A + I2uud2a_B) /  (2.0*ki);
        I3uuu2a_p   += dkk[i]*(I3uuu2a_A + I3uuu2a_B) /  (2.0*ki);
        I3uuu3a_p   += dkk[i]*(I3uuu3a_A + I3uuu3a_B) /  (2.0*ki);
        
              
        P13dd_A =   P13dd_B;      P13dd_B = 0.0;
        P13du_A =   P13du_B;      P13du_B = 0.0;
        P13uu_A =   P13uu_B;      P13uu_B = 0.0;
        sigma32PSL_A =   sigma32PSL_B;      sigma32PSL_B = 0.0;
        
        I1udd1a_A =   I1udd1a_B;      I1udd1a_B = 0.0;
        I2uud1a_A =   I2uud1a_B;      I2uud1a_B = 0.0;
        I2uud2a_A =   I2uud2a_B;      I2uud2a_B = 0.0;
        I3uuu2a_A =   I3uuu2a_B;      I3uuu2a_B = 0.0;
        I3uuu3a_A =   I3uuu3a_B;      I3uuu3a_B = 0.0;
      
    }

    pkl_k = Interpolation_nr(ki, kPKL, pPKL, nPKLT, pPKL2);
    P13dd_p      *= (rpow(ki,3.0)/FOURPI2)*pkl_k;
    P13du_p      *= (rpow(ki,3.0)/FOURPI2)*pkl_k;
    P13uu_p      *= (rpow(ki,3.0)/FOURPI2)*pkl_k;
    sigma32PSL_p *= (rpow(ki,3.0)/FOURPI2)*pkl_k;
    I1udd1a_p    *= (rpow(ki,3.0)/FOURPI2)*pkl_k;
    I2uud1a_p    *= (rpow(ki,3.0)/FOURPI2)*pkl_k;
    I2uud2a_p    *= (rpow(ki,3.0)/FOURPI2)*pkl_k;
    I3uuu2a_p    *= (rpow(ki,3.0)/FOURPI2)*pkl_k;
    I3uuu3a_p    *= (rpow(ki,3.0)/FOURPI2)*pkl_k;



	kFs(QRstmp)    = ki;
		
	P22dd(  QRstmp)      = P22dd_p;
	P22du(  QRstmp)      = P22du_p;
	P22uu(  QRstmp)      = P22uu_p;
	// A TNS
	I1udd1A(  QRstmp)      = I1udd1tA_p + 2.0*I1udd1a_p;
	I2uud1A(  QRstmp)      = I2uud1tA_p + 2.0*I2uud1a_p;
	I2uud2A(  QRstmp)      = I2uud2tA_p + 2.0*I2uud2a_p;
	I3uuu2A(  QRstmp)      = I3uuu2tA_p + 2.0*I3uuu2a_p;
	I3uuu3A(  QRstmp)      = I3uuu3tA_p + 2.0*I3uuu3a_p;
	// D function: B + C - G  
	I2uudd1BpC(  QRstmp)   = I2uudd1BpC_p 
								- ki*ki*gd.sigma2v*pkl_k;
	//~ I2uudd1BpC(  QRstmp)   = I2uudd1BpC_p;
	I2uudd2BpC(  QRstmp)   = I2uudd2BpC_p;
	I3uuud2BpC(  QRstmp)   = I3uuud2BpC_p
								- 2.0*ki*ki*gd.sigma2v*fk*pkl_k;
	//~ I3uuud2BpC(  QRstmp)   = I3uuud2BpC_p;
	I3uuud3BpC(  QRstmp)   = I3uuud3BpC_p;
	I4uuuu2BpC(  QRstmp)   = I4uuuu2BpC_p;
	I4uuuu3BpC(  QRstmp)   = I4uuuu3BpC_p
								- ki*ki*gd.sigma2v*fk*fk*pkl_k;
	//~ I4uuuu3BpC(  QRstmp)   = I4uuuu3BpC_p;
	I4uuuu4BpC(  QRstmp)   = I4uuuu4BpC_p;
	

	
	//  Bias
	Pb1b2(    QRstmp) = Pb1b2_p;
	Pb1bs2(   QRstmp) = Pb1bs2_p;
	Pb22(     QRstmp) = Pb22_p;
	Pb2s2(    QRstmp) = Pb2s2_p;
	Ps22(     QRstmp) = Ps22_p;
	Pb2theta( QRstmp) = Pb2theta_p;
	Pbs2theta(QRstmp) = Pbs2theta_p;
	//
	P13dd(  QRstmp)      = P13dd_p;
	P13du(  QRstmp)      = P13du_p;
	P13uu(  QRstmp)      = P13uu_p;
	
	sigma32PSL(QRstmp)   = sigma32PSL_p;
    

    free_dvector(dkk,1,cmd.nquadSteps);
    free_dvector(kk,1,cmd.nquadSteps);
    
    return *QRstmp;
}









































// END Qs and Rs

//~ global_kFs qrs;
global_kFs kfunctions;
//~ global_kFs qrs_nw;

//~ global global_kFs ki_functions_driver(real eta, real ki)
//~ global global_kFs ki_functions_driver(real ki)
global global_kFs ki_functions_driver(real ki, double kPKL[], double pPKL[], int nPKLT, double pPKL2[])
{
    //~ quadrature(ki);
    kfunctions = ki_functions(ki,kPKL,pPKL,nPKLT,pPKL2);
    //~ return qrs;
    return kfunctions;
}


#define ROMO 1
#define NULLMETHOD 0
#define TRAPEZOID 2
#define TRAPEZOID3 5





void quadraturemethod_string_to_int(string method_str,int *method_int)
{
    *method_int=-1;
    if (strcmp(method_str,"romberg") == 0) {
        *method_int = ROMO;
        strcpy(gd.quadraturemethod_comment, "romberg open quadrature method");
    }
//
    if (strcmp(method_str,"trapezoid") == 0) {
        *method_int = TRAPEZOID;
        strcpy(gd.quadraturemethod_comment, "trapezoid quadrature method");
    }
//
    if (strcmp(method_str,"trapezoid3") == 0) {
        *method_int = TRAPEZOID3;
        strcpy(gd.quadraturemethod_comment, "trapezoid3 quadrature method");
    }
//
    if (strnull(method_str)) {
        *method_int = NULLMETHOD;
        strcpy(gd.quadraturemethod_comment,
               "given null quadrature method ... running deafult (trapezoid)");
        fprintf(stdout,"\n\tintegration: default integration method (trapezoid)...\n");
    }
//
    if (*method_int == -1) {
        *method_int = TRAPEZOID;
        strcpy(gd.quadraturemethod_comment,
               "Unknown quadrature method ... running deafult (trapezoid)");
        fprintf(stdout,"\n\tquadrature: Unknown method... %s ",cmd.quadratureMethod);
        fprintf(stdout,
                "\n\trunning default quadrature method (trapezoid)...\n");
    }
}

#undef ROMO
#undef TRAPEZOID
#undef TRAPEZOID3
#undef NULLMETHOD



#undef KK
#undef QROMBERG




local  real Interpolation_nr(real k, double kPS[], double pPS[], int nPS, double pPS2[])
{
    real psftmp;   
    splint(kPS,pPS,pPS2,nPS,k,&psftmp);  
    return (psftmp);
}

local real sigma2L_function_int(real y)
{
    real p;
    real PSL;

    p = rpow(10.0,y);
    PSL = psInterpolation_nr(p, kPS, pPS, nPSLT);

    return p*PSL;
}


local real sigma2v_function_int(real y)
{
    real p;
    real PSL,fk;

    p = rpow(10.0,y);
    PSL = psInterpolation_nr(p, kPS, pPS, nPSLT);
    fk = Interpolation_nr(p, kPS, fkT, nPSLT, fkT2);
	fk /= gd.f0;
    return p*PSL*fk*fk;
}

local real Sigma2_int(real y)
{
    real p, PSL_nw, kosc;
    
    kosc=1.0/104.;
   

    p = rpow(10.0,y);
    PSL_nw = Interpolation_nr(p, kPS, pPS_nw, nPSLT,pPS2_nw);
    return p * PSL_nw * (1- rj0Bessel(p/kosc) + 2.*rj2Bessel(p/kosc) );
}



local real deltaSigma2_int(real y)
{
    real p, PSL_nw, kosc;
    
    kosc=1.0/104.;
   

    p = rpow(10.0,y);
    PSL_nw = Interpolation_nr(p, kPS, pPS_nw, nPSLT,pPS2_nw);
    return 3.0 * p * PSL_nw * rj2Bessel(p/kosc) ;
}




local real sigma_constants(void)
{

    real sigma2v, sigma2L, Sigma2, deltaSigma2;
    real kmin, kmax;
    real ymin, ymax, ymaxSigma;
    real EPSQ = 0.000001;
    int KK = 5;
	real ks; 


    kmin = kPS[1];
    kmax = kPS[nPSLT];
    ymin = rlog10(kmin);
    ymax = rlog10(kmax);

    sigma2v= (1.0/SIXPI2)*rlog(10.0)*
				qromo(sigma2v_function_int,ymin,ymax,midpnt,EPSQ,KK);
    gd.sigma2v = sigma2v;


    sigma2L= (1.0/SIXPI2)*rlog(10.0)*
				qromo(sigma2L_function_int,ymin,ymax,midpnt,EPSQ,KK);	
    gd.sigma2L = sigma2L;	
    

	ks = 0.4;
    ymaxSigma = rlog10(ks);
    
    Sigma2 = (1.0/SIXPI2)*rlog(10.0)*
				qromo(Sigma2_int,ymin,ymaxSigma,midpnt,EPSQ,KK);	
    gd.Sigma2 = Sigma2;
    
    deltaSigma2 = (1.0/SIXPI2)*rlog(10.0)*
				qromo(deltaSigma2_int,ymin,ymaxSigma,midpnt,EPSQ,KK);	
    gd.deltaSigma2 = deltaSigma2;
    
};




















