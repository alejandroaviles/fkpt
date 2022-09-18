/*==============================================================================
 NAME: rsd.c				[code for redshift space correlation function - GSM]
 Alejandro Aviles (avilescervantes@gmail.com)
 ================================================================================ 
*/
#include "rsd.h"

 

global void compute_rsdmultipoles(void)
{
    stream outstr;
    //~ pointrfunctionsTableptr p;
    real aTime;
    int counter, i;
    real k, kmin, kmax;
    int Nk;
	real b1, b2, bs2, b3nl, alpha0, alpha2, alpha4, ctilde;
	real PshotP, alpha0shot, alpha2shot;
		

	b1= cmd.b1;
	b2= cmd.b2;
	bs2= cmd.bs2;
	b3nl= cmd.b3nl;
	alpha0= cmd.alpha0;
	alpha2= cmd.alpha2;
	alpha4= cmd.alpha4;
	ctilde= cmd.ctilde;
	PshotP= cmd.PshotP;
    alpha0shot = cmd.alpha0shot;
    alpha2shot = cmd.alpha2shot;
 
	
	kmin = cmd.kmin;
	kmax = cmd.kmax;
	Nk = cmd.Nk;


	
	//~ /* TESTS */
	//~ real pkEFTIR, pkloop, pkloop_nw, pkl, mu, kk,pkkaiser,fk, pkME;
	//~ int ii=40;	
	//~ kk       = kFArrays.kT[ii]; 
	//~ mu=0.5	;
	//~ pkl     = kFArrays.pklT[ii];
	//~ fk	    = kFArrays.fkT[ii];	
	//~ pkEFTIR = pk_EFT_IR(ii, mu, b1, b2, bs2, b3nl, 
			//~ alpha0, alpha2, alpha4, ctilde, 
			//~ PshotP, alpha0shot, alpha2shot);
	//~ pkloop = pk_loop(ii, mu, b1, b2, bs2, b3nl, 
			//~ alpha0, alpha2, alpha4, ctilde, 
			//~ PshotP, alpha0shot, alpha2shot);
	//~ pkloop_nw = pk_loop_nw(ii, mu, b1, b2, bs2, b3nl, 
			//~ alpha0, alpha2, alpha4, ctilde, 
			//~ PshotP, alpha0shot, alpha2shot);
	//~ pkkaiser = rpow(b1 + fk * mu*mu,2) * pkl;
	//~ pkME = pk_ME(ii, mu, b1, b2, bs2, b3nl);
	//~ fprintf(stdout, "\n\nk=%g, mu=%g, pkl=%g, pkKaiser=%g, pk_ME =%g, pkloop=%g, pkIR=%g\n"
	//~ ,kk,mu, pkl,pkkaiser,pkME, pkloop,pkEFTIR)	;	
	//~ /* END TESTS */



    global_rsdmultipoles pk_ells; 
    
	outstr = stropen(gd.fpfnamersd,"w!");


                
	fprintf(outstr,"# InputPklFile=%s, redshift z=%g, OmegaM=%g, h = %g, f0=%g, \
	b1=%g, b2=%g, bs2=%g, b3nl=%g, alpha0=%g, alpha2=%g, alpha4=%g, c~=%g, \
	PshotP=%g, alpha0shot=%g, alpha4shot=%g\n",
                cmd.fnamePS, cmd.xstop, cmd.om, cmd.h, gd.f0, b1, b2, 
				bs2, b3nl, alpha0, alpha2, alpha4, ctilde, PshotP, 
				alpha0shot, alpha2shot );       
                                  

    fprintf(outstr,"%1s%12s%16s%16s%16s",
            "#","1.s[Mpc/h]","2.monopole","3.quadrupole","4.hexadecapole\n");
    
        
    fprintf(stdout,"\nrsd: ");
    fprintf(stdout," computing... ");
	aTime = second();	
    
    for(i=0; i<cmd.Nk; i++) {
	
		pk_ells = pk_rsd_multipoles(i, b1, b2, bs2, b3nl, 
			alpha0, alpha2, alpha4, ctilde, PshotP, alpha0shot, alpha2shot);  
	
        fprintf(outstr,"%e %e %e %e\n",
                pk_ells.k,  // [Mpc/h]
                pk_ells.P0,  // [monopole]
                pk_ells.P2,  // [quadrupole]
                pk_ells.P4  // [hexadecapole]
                );        
	};	
	fclose(outstr);
    fprintf(stdout,"...time = %g seconds \n",second()-aTime); 
 
 
}   





global_rsdmultipoles pk_rsd_multipoles(int i, real b1, real b2, 
			real bs2, real b3nl, 
			real alpha0, real alpha2, real alpha4, real ctilde, 
			real PshotP, real alpha0shot, real alpha2shot)
{
	
    real k, pk_IR; 
    real P0=0.0, P2=0.0, P4=0.0;		
	
    global_rsdmultipoles_ptr rsdfunp;	
    rsdfunp = (global_rsdmultipoles_ptr) allocate(1 * sizeof(global_rsdmultipoles)); 	
	
	int NGL, integer;
    real *muGL,*wGL, mu, mu2, w, legP2, legP4;


	NGL=16;

	muGL=dvector(1,NGL);
	wGL=dvector(1,NGL);
	gauleg(-1.0,1.0,muGL,wGL,NGL);
		

	for(int muii=1; muii<=NGL; muii++)
	{	
			
		mu = muGL[muii];
		w = wGL[muii];
		mu2 = mu*mu;
		
		legP2 = 0.5 * (3.0 * mu2-1.0);
		legP4 = (35.0*mu2*mu2 - 30.0* mu2 +3.0 )/8.0;
			
		pk_IR = pk_EFT_IR(i, mu, b1, b2, bs2, b3nl, 
			alpha0, alpha2, alpha4, ctilde, 
			PshotP, alpha0shot, alpha2shot);
		
		P0 += w*pk_IR* 0.5;
		P2 += w*pk_IR* 5./2.* legP2;
		P4 += w*pk_IR* 9./2.* legP4; //Modificacion-LCDMfk	

	}; 
					
	k_rsdmultipoles( rsdfunp) = kFArrays.kT[i]; 
	P0_rsdmultipoles(rsdfunp) = P0;
	P2_rsdmultipoles(rsdfunp) = P2;
	P4_rsdmultipoles(rsdfunp) = P4;
     
    return *rsdfunp;	
}




local real pk_EFT_IR(int i, real mu, real b1, real b2, real bs2, real b3nl, 
			real alpha0, real alpha2, real alpha4, real ctilde, 
			real PshotP, real alpha0shot, real alpha2shot)
{
    real mu2, f0, k, fk;
    real Sigma2Total, ExpDamp;
    real pk_IR, pkloop_nw, pkloop, pkl, pkl_nw;
    real pkKaiser;


	k       = kFArrays.kT[i];
	pkl     = kFArrays.pklT[i];
	pkl_nw	= kFArrays_nw.pklT[i];
	fk	    = kFArrays.fkT[i]; 
	f0      = gd.f0;
	

	mu2=mu*mu;
	
	Sigma2Total = (1 + f0 * mu2* (2. + f0)) * gd.Sigma2 
		+ f0*f0 * mu2 * (mu2 - 1.0) * gd.deltaSigma2;

	ExpDamp = rexp(-k*k * Sigma2Total);
		
	pkloop = pk_loop(i, mu, b1, b2, bs2, b3nl, 
			alpha0, alpha2, alpha4, ctilde, 
			PshotP, alpha0shot, alpha2shot);
    
	pkloop_nw = pk_loop_nw(i, mu, b1, b2, bs2, b3nl, 
			alpha0, alpha2, alpha4, ctilde, 
			PshotP, alpha0shot, alpha2shot);
    
    
	pk_IR = rpow(b1 + fk * mu2,2) * 
			( pkl_nw + ExpDamp * (pkl - pkl_nw) * (1 + k*k*Sigma2Total) ) 
			+ ExpDamp * pkloop + (1. - ExpDamp ) * pkloop_nw;
			
			
	
	//~ pkKaiser = rpow(b1 + fk * mu2,2) * pkl;		
	
	//~ fprintf(stdout, "\n\nk=%g,  pkl=%g, pkKaiser=%g, pkloop=%g, pkIR=%g, inside\n",k,pkl,pkKaiser,pkloop, pk_IR);		
	return pk_IR;			
}


	
local real pk_loop(int i, real mu, real b1, real b2, real bs2, real b3nl, 
			real alpha0, real alpha2, real alpha4, real ctilde, 
			real PshotP, real alpha0shot, real alpha2shot)
{
	real k, pkl, P22dd, P22du, P22uu, P13dd, P13du, P13uu, fk;
	real I1udd1A, I2uud1A, I2uud2A, I3uuu2A, I3uuu3A;
	real I2uudd1D,I2uudd2D,I3uuud2D,I3uuud3D,I4uuuu2D,I4uuuu3D,I4uuuu4D;
	real Pb1b2, Pb1bs2, Pb22, Pb2s2, Pbs22, Pb2theta, Pbs2theta, sigma32pkl; 
	real ATNS, D, PddXloop, PduXloop, PuuXloop, GTNS; 
	real f0, mu2, sigma2v;
	real pPMEloop, Pctilde, PKaiser, Peft, Pshot, pkloop;
	
	f0 =gd.f0;
	mu2 = mu*mu;
	sigma2v = gd.sigma2v;
		
		k     = kFArrays.kT[i];
		pkl   = kFArrays.pklT[i];
		fk    = kFArrays.fkT[i];
		P22dd = kFArrays.P22ddT[i];
		P22du = kFArrays.P22duT[i];
		P22uu = kFArrays.P22uuT[i];
		P13dd = kFArrays.P13ddT[i];
		P13du = kFArrays.P13duT[i];
		P13uu = kFArrays.P13uuT[i];
		// A 
		I1udd1A  = kFArrays.I1udd1AT[i];
		I2uud1A  = kFArrays.I2uud1AT[i];
		I2uud2A  = kFArrays.I2uud2AT[i];
		I3uuu2A  = kFArrays.I3uuu2AT[i];
		I3uuu3A  = kFArrays.I3uuu3AT[i];
		//  D   
		I2uudd1D = kFArrays.I2uudd1BpCT[i];
		I2uudd2D = kFArrays.I2uudd2BpCT[i];
		I3uuud2D = kFArrays.I3uuud2BpCT[i];
		I3uuud3D = kFArrays.I3uuud3BpCT[i];
		I4uuuu2D = kFArrays.I4uuuu2BpCT[i];
		I4uuuu3D = kFArrays.I4uuuu3BpCT[i];
		I4uuuu4D = kFArrays.I4uuuu4BpCT[i];
		//  Bias
		Pb1b2      = kFArrays.Pb1b2T[i];
		Pb1bs2     = kFArrays.Pb1bs2T[i];
		Pb22       = kFArrays.Pb22T[i];
		Pb2s2      = kFArrays.Pb2s2T[i];
		Pbs22      = kFArrays.Ps22T[i];
		Pb2theta   = kFArrays.Pb2thetaT[i];
		Pbs2theta  = kFArrays.Pbs2thetaT[i];
		sigma32pkl = kFArrays.sigma32PSLT[i];
	
	PddXloop = b1*b1 * (P22dd+P13dd) + 2.*b1*b2*Pb1b2 + 2*b1*bs2 * Pb1bs2 
			+ b2*b2 * Pb22 + bs2 * bs2 * Pbs22 + 2.0*b2*bs2 * Pb2s2 
			+ 2.0*b1*b3nl * sigma32pkl;

	PduXloop = b1 * (P22du+P13du) + b2 * Pb2theta + bs2 * Pbs2theta 
			+ b3nl * fk/f0 * sigma32pkl;
	PuuXloop = P22uu+P13uu;
	
	ATNS = rpow(b1,3) * A_TNS(f0/b1, mu, I1udd1A, I2uud1A, I2uud2A, I3uuu2A, I3uuu3A);
	
	D = rpow(b1,4) * D_function(f0/b1, mu, I2uudd1D,I2uudd2D,I3uuud2D,I3uuud3D,
						I4uuuu2D,I4uuuu3D,I4uuuu4D);

    //~ GTNS = rpow(k*mu*f0,2) * sigma2v * (b1*b1  + 2*b1*mu2 * fk 
			//~ +  mu2*mu2*fk*fk)* pkl;                               
    GTNS = 0;	                               
                                   					
						
	pPMEloop = PddXloop + 2.0 * mu2 * f0 * PduXloop 
		+ mu2*mu2 * f0*f0 * PuuXloop + ATNS + D - GTNS;


	PKaiser = rpow(b1 + mu2 * fk,2 ) * pkl ;
	Pctilde = ctilde * rpow(mu*k*f0,4) *sigma2v*sigma2v * PKaiser;
	Peft =(alpha0 + alpha2*mu2 + alpha4*mu2*mu2 ) *k*k*pkl + Pctilde;
	
	Pshot= PshotP*(alpha0shot + alpha2shot*k*k*mu2);

	pkloop= pPMEloop+Peft+Pshot;
	
	//~ fprintf(stdout,"running %g", pkloop);
	
	return pkloop;					
							
};			



	
local real pk_loop_nw(int i, real mu, real b1, real b2, real bs2, real b3nl, 
			real alpha0, real alpha2, real alpha4, real ctilde, 
			real PshotP, real alpha0shot, real alpha2shot)
{
	real k, pkl_nw, P22dd, P22du, P22uu, P13dd, P13du, P13uu, fk;
	real I1udd1A, I2uud1A, I2uud2A, I3uuu2A, I3uuu3A;
	real I2uudd1D,I2uudd2D,I3uuud2D,I3uuud3D,I4uuuu2D,I4uuuu3D,I4uuuu4D;
	real Pb1b2, Pb1bs2, Pb22, Pb2s2, Pbs22, Pb2theta, Pbs2theta, sigma32pkl; 
	real ATNS, D, PddXloop, PduXloop, PuuXloop, GTNS; 
	real f0, mu2, sigma2v;
	real pPMEloop, Pctilde, PKaiser, Peft, Pshot, pkloop_nw;
	
	f0 =gd.f0;
	mu2 = mu*mu;
	sigma2v = gd.sigma2v;
		
		k        = kFArrays_nw.kT[i];
		pkl_nw   = kFArrays_nw.pklT[i];
		fk    = kFArrays_nw.fkT[i];
		P22dd = kFArrays_nw.P22ddT[i];
		P22du = kFArrays_nw.P22duT[i];
		P22uu = kFArrays_nw.P22uuT[i];
		P13dd = kFArrays_nw.P13ddT[i];
		P13du = kFArrays_nw.P13duT[i];
		P13uu = kFArrays_nw.P13uuT[i];
		// A 
		I1udd1A  = kFArrays_nw.I1udd1AT[i];
		I2uud1A  = kFArrays_nw.I2uud1AT[i];
		I2uud2A  = kFArrays_nw.I2uud2AT[i];
		I3uuu2A  = kFArrays_nw.I3uuu2AT[i];
		I3uuu3A  = kFArrays_nw.I3uuu3AT[i];
		//  D   
		I2uudd1D = kFArrays_nw.I2uudd1BpCT[i];
		I2uudd2D = kFArrays_nw.I2uudd2BpCT[i];
		I3uuud2D = kFArrays_nw.I3uuud2BpCT[i];
		I3uuud3D = kFArrays_nw.I3uuud3BpCT[i];
		I4uuuu2D = kFArrays_nw.I4uuuu2BpCT[i];
		I4uuuu3D = kFArrays_nw.I4uuuu3BpCT[i];
		I4uuuu4D = kFArrays_nw.I4uuuu4BpCT[i];
		//  Bias
		Pb1b2      = kFArrays_nw.Pb1b2T[i];
		Pb1bs2     = kFArrays_nw.Pb1bs2T[i];
		Pb22       = kFArrays_nw.Pb22T[i];
		Pb2s2      = kFArrays_nw.Pb2s2T[i];
		Pbs22      = kFArrays_nw.Ps22T[i];
		Pb2theta   = kFArrays_nw.Pb2thetaT[i];
		Pbs2theta  = kFArrays_nw.Pbs2thetaT[i];
		sigma32pkl = kFArrays_nw.sigma32PSLT[i];
	
	
	PddXloop = b1*b1 * (P22dd+P13dd) + 2.*b1*b2*Pb1b2 + 2*b1*bs2 * Pb1bs2 
			+ b2*b2 * Pb22 + bs2 * bs2 * Pbs22 + 2.0*b2*bs2 * Pb2s2 
			+ 2.0*b1*b3nl * sigma32pkl;

	PduXloop = b1 * (P22du+P13du) + b2 * Pb2theta + bs2 * Pbs2theta 
			+ b3nl * fk/f0 * sigma32pkl;
	PuuXloop = P22uu+P13uu;
	
	ATNS = rpow(b1,3) * A_TNS(f0/b1, mu, I1udd1A, I2uud1A, I2uud2A, I3uuu2A, I3uuu3A);
	
	D = rpow(b1,4) * D_function(f0/b1, mu, I2uudd1D,I2uudd2D,I3uuud2D,I3uuud3D,
						I4uuuu2D,I4uuuu3D,I4uuuu4D);
    //~ GTNS = rpow(k*mu*f0,2) * sigma2v * (b1*b1  + 2*b1*mu2 * fk 
			//~ +  mu2*mu2*fk*fk)* pkl_nw; 		 SAME RESULT	
	GTNS = 0;					
						
	pPMEloop = PddXloop + 2.0 * mu2 * f0 * PduXloop 
		+ mu2*mu2 * f0*f0 * PuuXloop + ATNS + D - GTNS;


	PKaiser = rpow(b1 + mu2 * fk,2 ) * pkl_nw ;
	Pctilde = ctilde * rpow(mu*k*f0,4) *sigma2v*sigma2v * PKaiser;
	Peft =(alpha0 + alpha2*mu2 + alpha4*mu2*mu2 ) *k*k*pkl_nw + Pctilde;
	
	Pshot= PshotP*(alpha0shot + alpha2shot*k*k*mu2);

	pkloop_nw= pPMEloop+Peft+Pshot;
	
	//~ fprintf(stdout,"running %g", pkloop);
	
	return pkloop_nw;					
							
};			







local real A_TNS(real f0, real mu, real I1udd1A, real I2uud1A, 
		real I2uud2A, real I3uuu2A, real I3uuu3A)
{
	real result;
	real mu2, mu4, mu6, f02, f03;
	mu2=mu*mu;
	mu4=mu2*mu2;
	mu6=mu4*mu2;
	f02=f0*f0;
	f03=f02*f0;
	
    result = f0  *  I1udd1A*mu2 
           + f02 * (I2uud1A*mu2 + I2uud2A*mu4)  
		   + f03 * (I3uuu2A*mu4 + I3uuu3A*mu6);
			
	return result;
}


local real D_function(real f0, real mu, real I2uudd1D,real I2uudd2D,
	real I3uuud2D,real I3uuud3D,real I4uuuu2D,real I4uuuu3D,real I4uuuu4D)
{
	real result;
	real mu2, mu4, mu6, mu8, f02,f03,f04;
	mu2=mu*mu;
	mu4=mu2*mu2;
	mu6=mu4*mu2;
	mu8=mu6*mu2;
	f02=f0*f0;
	f03=f02*f0;
	f04=f03*f0;
	
    result = f02*(I2uudd1D*mu2 + I2uudd2D*mu4) + 
		f03*(I3uuud2D*mu4 + I3uuud3D*mu6) + 
		f04*(I4uuuu2D*mu4 + I4uuuu3D*mu6 + I4uuuu4D*mu8);
			
	return result;
}







// NOt used, but is good to have it here
local real pk_ME(int i, real mu, real b1, real b2, real bs2, real b3nl)
{
	real k, pkl, P22dd, P22du, P22uu, P13dd, P13du, P13uu, fk;
	real I1udd1A, I2uud1A, I2uud2A, I3uuu2A, I3uuu3A;
	real I2uudd1D,I2uudd2D,I3uuud2D,I3uuud3D,I4uuuu2D,I4uuuu3D,I4uuuu4D;
	real Pb1b2, Pb1bs2, Pb22, Pb2s2, Pbs22, Pb2theta, Pbs2theta, sigma32pkl; 
	real ATNS, D, PddXloop, PduXloop, PuuXloop; 
	real f0, mu2, sigma2v;
	real pPMEloop;
	
	f0 =gd.f0;
	mu2 = mu*mu;
	sigma2v = gd.sigma2v;
		
		k     = kFArrays.kT[i];
		pkl   = kFArrays.pklT[i];
		fk    = kFArrays.fkT[i];
		P22dd = kFArrays.P22ddT[i];
		P22du = kFArrays.P22duT[i];
		P22uu = kFArrays.P22uuT[i];
		P13dd = kFArrays.P13ddT[i];
		P13du = kFArrays.P13duT[i];
		P13uu = kFArrays.P13uuT[i];
		// A 
		I1udd1A  = kFArrays.I1udd1AT[i];
		I2uud1A  = kFArrays.I2uud1AT[i];
		I2uud2A  = kFArrays.I2uud2AT[i];
		I3uuu2A  = kFArrays.I3uuu2AT[i];
		I3uuu3A  = kFArrays.I3uuu3AT[i];
		//  D   
		I2uudd1D = kFArrays.I2uudd1BpCT[i];
		I2uudd2D = kFArrays.I2uudd2BpCT[i];
		I3uuud2D = kFArrays.I3uuud2BpCT[i];
		I3uuud3D = kFArrays.I3uuud3BpCT[i];
		I4uuuu2D = kFArrays.I4uuuu2BpCT[i];
		I4uuuu3D = kFArrays.I4uuuu3BpCT[i];
		I4uuuu4D = kFArrays.I4uuuu4BpCT[i];
		//  Bias
		Pb1b2      = kFArrays.Pb1b2T[i];
		Pb1bs2     = kFArrays.Pb1bs2T[i];
		Pb22       = kFArrays.Pb22T[i];
		Pb2s2      = kFArrays.Pb2s2T[i];
		Pbs22       = kFArrays.Ps22T[i];
		Pb2theta   = kFArrays.Pb2thetaT[i];
		Pbs2theta  = kFArrays.Pbs2thetaT[i];
		sigma32pkl = kFArrays.sigma32PSLT[i];
	
	
	PddXloop = b1*b1 * (P22dd+P13dd) + 2.*b1*b2*Pb1b2 + 2*b1*bs2 * Pb1bs2 
			+ b2*b2 * Pb22 + bs2 * bs2 * Pbs22 + 2.0*b2*bs2 * Pb2s2 
			+ 2.0*b1*b3nl * sigma32pkl;

	PduXloop = b1 * (P22du+P13du) + b2 * Pb2theta + bs2 * Pbs2theta 
			+ b3nl * fk/f0 * sigma32pkl;
	PuuXloop = P22uu+P13uu;
	
	ATNS = rpow(b1,3) * A_TNS(f0/b1, mu, I1udd1A, I2uud1A, I2uud2A, I3uuu2A, I3uuu3A);
	
	D = rpow(b1,4) * D_function(f0/b1, mu, I2uudd1D,I2uudd2D,I3uuud2D,I3uuud3D,
						I4uuuu2D,I4uuuu3D,I4uuuu4D);
						
						
	pPMEloop = PddXloop + 2.0 * mu2 * f0 * PduXloop 
		+ mu2*mu2 * f0*f0 * PuuXloop + ATNS + D;
		
	fprintf(stdout, "\n\nk=%g, mu=%g, PddXloop=%g, PduXloop=%g, PuuXloop =%g, ATNS=%g, D=%g\n"
	,k,mu, PddXloop,PduXloop,PuuXloop, ATNS,D)	;			

	return pPMEloop;
};

