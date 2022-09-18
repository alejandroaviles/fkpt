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
		
	//~ b1= cmd.b1;
	//~ b2= cmd.b2;
	//~ bs= cmd.bs;
	b1= 1.7;
	b2= 0.2;
	bs2= 0.2;
	b3nl= 0.1;
	alpha0= -50;
	alpha2= -40;
	alpha4= -30;
	ctilde= -0.1;
		
	PshotP= 500;
    alpha0shot = 1.5;
    alpha2shot = 2.2;
    
	
	
	kmin = cmd.kmin;
	kmax = cmd.kmax;
	Nk = cmd.Nk;

	real pkEFTIR, pkloop, pkl, mu, kk,pkkaiser,fk, pkME;
	int ii=40;
	
	kk       = kFArrays.kT[ii]; 
	mu=0.5	;
	pkl     = kFArrays.pklT[ii];
	fk	    = kFArrays.fkT[ii];	
	
	
	
	pkEFTIR = pk_EFT_IR(ii, mu, b1, b2, bs2, b3nl, 
			alpha0, alpha2, alpha4, ctilde, 
			PshotP, alpha0shot, alpha2shot);
	pkloop = pk_loop(ii, mu, b1, b2, bs2, b3nl, 
			alpha0, alpha2, alpha4, ctilde, 
			PshotP, alpha0shot, alpha2shot);
	pkkaiser = rpow(b1 + fk * mu*mu,2) * pkl;
	pkME = pk_ME(ii, mu, b1, b2, bs2, b3nl);
	fprintf(stdout, "\n\nk=%g, mu=%g, pkl=%g, pkKaiser=%g, pk_ME =%g, pkloop=%g, pkIR=%g\n"
	,kk,mu, pkl,pkkaiser,pkME, pkloop,pkEFTIR)	;	


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


	NGL=6;

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
		P4 += w*pk_IR* 7./2.* legP4;	

	}; 
					
	k_rsdmultipoles( rsdfunp) = kFArrays.kT[i]; 
	P0_rsdmultipoles(rsdfunp) = P0;
	P2_rsdmultipoles(rsdfunp) = P2;
	P4_rsdmultipoles(rsdfunp) = P4;
     
    return *rsdfunp;	
}




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
	
local real pk_loop(int i, real mu, real b1, real b2, real bs2, real b3nl, 
			real alpha0, real alpha2, real alpha4, real ctilde, 
			real PshotP, real alpha0shot, real alpha2shot)
{
	real k, pkl, P22dd, P22du, P22uu, P13dd, P13du, P13uu, fk;
	real I1udd1A, I2uud1A, I2uud2A, I3uuu2A, I3uuu3A;
	real I2uudd1D,I2uudd2D,I3uuud2D,I3uuud3D,I4uuuu2D,I4uuuu3D,I4uuuu4D;
	real Pb1b2, Pb1bs2, Pb22, Pb2s2, Pbs22, Pb2theta, Pbs2theta, sigma32pkl; 
	real ATNS, D, PddXloop, PduXloop, PuuXloop; 
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


	PKaiser = rpow(b1 + mu2 * fk,2 ) * pkl ;
	Pctilde = ctilde * rpow(mu*k*f0,4) *sigma2v*sigma2v * PKaiser;
	Peft =(alpha0 + alpha2*mu2 + alpha4*mu2*mu2 ) *k*k*pkl + Pctilde;
	
	Pshot= PshotP*(alpha0shot + alpha2shot*k*k*mu2);

	pkloop= pPMEloop+Peft+Pshot;
	
	//~ fprintf(stdout,"running %g", pkloop);
	
	return pkloop;					
							
};			





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
	//~ Sigma2Total=0.0;
	ExpDamp = rexp(-k*k * Sigma2Total);
		
	pkloop = pk_loop(i, mu, b1, b2, bs2, b3nl, 
			alpha0, alpha2, alpha4, ctilde, 
			PshotP, alpha0shot, alpha2shot);
    
	pkloop_nw = pk_loop(i, mu, b1, b2, bs2, b3nl, 
			alpha0, alpha2, alpha4, ctilde, 
			PshotP, alpha0shot, alpha2shot);
    
    //~ pkloop=0.0;
    
	pk_IR = rpow(b1 + fk * mu2,2) * 
			( pkl_nw + ExpDamp * (pkl - pkl_nw) * (1 + k*k*Sigma2Total) ) 
			+ ExpDamp * pkloop + (1. - ExpDamp ) * pkloop_nw;
	
	pkKaiser = rpow(b1 + fk * mu2,2) * pkl;		
	
	//~ fprintf(stdout, "\n\nk=%g,  pkl=%g, pkKaiser=%g, pkloop=%g, pkIR=%g, inside\n",k,pkl,pkKaiser,pkloop, pk_IR);		
	return pk_IR;			
}








		//


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










/*
global_rsdmultipoles pk_rsd_multipoles(int i, real b1, real b2, 
			real bs2, real b3nl, 
			real alpha0, real alpha2, real alpha4, real ctilde, 
			real PshotP, real alpha0shot, real alpha2shot)
{
	
	real k, pkl, pkl_nw, pkloop, pkloop_nw, fk, f0, pk_IR;	
	real Sigma2Total, ExpDamp;
    real P0=0.0, P2=0.0, P4=0.0;		
	
    global_rsdmultipoles_ptr rsdfunp;	
    rsdfunp = (global_rsdmultipoles_ptr) allocate(1 * sizeof(global_rsdmultipoles)); 	
	
	int NGL, integer;
    real *muGL,*wGL, mu, mu2, w, legP2, legP4;

	k       = kFArrays.kT[i];
	pkl     = kFArrays.pklT[i];
	pkl_nw	= kFArrays_nw.pklT[i];
	fk	    = kFArrays.fkT[i]; 
	f0      = gd.f0;

	NGL=6;

	muGL=dvector(1,NGL);
	wGL=dvector(1,NGL);
	gauleg(-1.0,1.0,muGL,wGL,NGL);
		

	for(int muii=1; muii<=NGL; muii++)
	{	
			
		mu = muGL[muii];
		w = wGL[muii];
		mu2 = mu*mu;
		
		//~ if(i==1){
			//~ fprintf(stdout,"\nmuii=%d, mu=%g, w=%g", muii,mu,w);
		//~ };
		
		
		legP2 = 0.5 * (3.0 * mu2-1.0);
		legP4 = (35.0*mu2*mu2 - 30.0* mu2 +3.0 )/8.0;
		
		
		Sigma2Total = (1 + f0 * mu2* (2. + f0)) * gd.Sigma2 
			+ f0*f0 * mu2 * (mu2 - 1.0) * gd.deltaSigma2;

		ExpDamp = rexp(-k*k * Sigma2Total);
		
		pkloop = pk_loop(i, mu, b1, b2, bs2, b3nl, 
			alpha0, alpha2, alpha4, ctilde, 
			PshotP, alpha0shot, alpha2shot);
    
		pkloop_nw = pk_loop(i, mu, b1, b2, bs2, b3nl, 
			alpha0, alpha2, alpha4, ctilde, 
			PshotP, alpha0shot, alpha2shot);
			
    
		pk_IR = rpow(b1 + fk * mu2,2) * 
			( pkl_nw + ExpDamp * (pkl - pkl_nw) * (1 + k*k*Sigma2Total) ) 
			+ ExpDamp * pkloop + (1. - ExpDamp ) * pkloop_nw;
		

			
			P0 += w*pk_IR* 0.5;
			P2 += w*pk_IR* 5./2.* legP2;
			P4 += w*pk_IR* 7./2.* legP4;	

		}; 
			
		
	k_rsdmultipoles( rsdfunp) = k; 
	P0_rsdmultipoles(rsdfunp) = P0;
	P2_rsdmultipoles(rsdfunp) = P2;
	P4_rsdmultipoles(rsdfunp) = P4;
     
    return *rsdfunp;	
}
*/




 

/*
we integrate over the variable y = r_parallel - s mu, 
 mu is the angle betwen vec-s and the line of sight 
 s^2 = r_parallel^2 + s_perpendicular^2 

global_rsdmultipoles CF_rsd_multipoles(real s, real b1, real b2, real bs, real c1eft, real c2eft, real sigma2eft)
{
    int i, j;
    real *muGL, *wGL;
    int gsm_NGL;	
    real mu, w;
	
	int gsm_sizeyT;
	real gsm_width;  
	real y, ymin,deltay,delta,yprev;
	real legP2, legP4;
	real rpar,spar,sperp,r,sigma2,intvars,integrand;
	real sigma2parVal, sigma2perpVal, v12Val, xiVal, eftFoG;	
	real nearlyzero = 0.00000000001;


    real xis0p=0.0, xis0A=0.0, xis0B=0.0;	
    real xis2p=0.0, xis2A=0.0, xis2B=0.0;	
    real xis4p=0.0, xis4A=0.0, xis4B=0.0;	
	
    global_rsdmultipoles_ptr rsdfunp;	
    rsdfunp = (global_rsdmultipoles_ptr) allocate(1 * sizeof(global_rsdmultipoles));   

	gsm_width=cmd.gsm_width;
	gsm_sizeyT =cmd.gsm_sizeyT;
	gsm_NGL=cmd.gsm_NGL;
	
	//~ gsm_width=30.0;
	//~ gsm_sizeyT=120;
	//~ gsm_NGL=16;	
	
	delta = gsm_width/( (real)gsm_sizeyT - 1);
	ymin=0.0;


	muGL=dvector(1,gsm_NGL);
	wGL=dvector(1,gsm_NGL);
	gauleg(-1.0,1.0,muGL,wGL,gsm_NGL);

	yprev=0.0;
	for(int yii=1;yii<=gsm_sizeyT; yii++)
	{
		y = ymin + (real)(yii-1)*delta;
		deltay= y-yprev;
		
		//~ if (s ==10.){
			//~ fprintf(stdout,"\n%g, %g ",y,deltay);
		//~ }

		for(int muii=1; muii<=gsm_NGL; muii++)
		{	
			
			mu = muGL[muii];
			w = wGL[muii];
			legP2 = 0.5 * (3.0 * mu*mu-1.0);
			legP4 = (35.0*rpow(mu,4.0) - 30.0* rpow(mu,2.0) +3.0 )/8.0;
			
			
		//~ if (s ==10.){
			//~ if(y==ymin){
			//~ fprintf(stdout,"\n%g, %g, ",mu,w);
		    //~ }
		//~ }			
			
		
			
			rpar= y + s * mu;
			spar = mu * s;
			sperp = rsqrt(1.0 - mu*mu) * s;
			r = rsqrt(rpar*rpar + sperp*sperp);
			
			eftFoG = sigma2eft * (1.0 + xi_zaF(r))/(1.0 + xiVal);
			
			xiVal= xi12F(r, b1, b2, bs, c1eft, c2eft, sigma2eft);	
			v12Val= v12F(r, b1, b2, bs, c1eft, c2eft, sigma2eft) / (1+xiVal);		
			sigma2parVal =sigma2par12F(r, b1, b2, bs, c1eft, c2eft, sigma2eft) / (1+xiVal)
											- v12Val*v12Val + eftFoG;
			sigma2perpVal= sigma2perp12F(r, b1, b2, bs, c1eft, c2eft, sigma2eft)/ (1+xiVal) 
											+ eftFoG;
			
			sigma2 = rpow(rpar/r,2.0) * sigma2parVal + (1- rpow(rpar/r,2)) * sigma2perpVal;
			
			if (sigma2<nearlyzero) sigma2=nearlyzero;
			
			intvars = (1. + xiVal) / rsqrt(TWO_PI * sigma2 ) * rexp( - rpow(spar - rpar - v12Val*rpar/r,2) / (2.0*sigma2)  );
			integrand = intvars;

			
			xis0B += w*integrand*0.5;
			xis2B += w*integrand*5./2.* legP2;
			xis4B += w*integrand*7./2.* legP4;	
		}; 
		
		if (yii==1) xis0A=xis0B;
		
		xis0p += (xis0A+xis0B)/2. * deltay;
		xis2p += (xis2A+xis2B)/2. * deltay;
		xis4p += (xis4A+xis4B)/2. * deltay;
		
		xis0A=xis0B;xis0B=0;
		xis2A=xis2B;xis2B=0;
		xis4A=xis4B;xis4B=0;
		
		yprev = y;
	};
    
    xis0p *=   2.0;
    xis0p +=  -1.0;
    xis2p *=   2.0;
    xis4p *=   2.0;
    
     s_rsdmultipoles(rsdfunp) = s; 
     xi_0_rsdmultipoles(rsdfunp) = xis0p;
     xi_2_rsdmultipoles(rsdfunp) = xis2p;
     xi_4_rsdmultipoles(rsdfunp) = xis4p;
     
    return *rsdfunp;     
	
};

 */





// (density weighted) velocity moments
/*
local real xi12F(real r, real b1, real b2, real bs, real c1eft, real c2eft, real sigma2eft)
{
    real func;
    func = xi_00F(r) + b1* xi_10F(r) + b1*b1 * xi_20F(r)
                 + b2*xi_01F(r)+ b2*b2*xi_02F(r)+ b1*b2*xi_11F(r)
                 + bs*xi_001F(r)+ bs*bs*xi_002F(r)+ b1*bs*xi_101F(r)
                 + b2*bs*xi_011F(r) + c1eft * xi_eftF(r);
    return (func);	
}	



local real v12F(real r, real b1, real b2, real bs, real c1eft, real c2eft, real sigma2eft)
{
    real func;
    func = v_00F(r) + b1* v_10F(r) + b1*b1 * v_20F(r)
                 + b2*v_01F(r) + b1*b2*v_11F(r)
                 + bs*v_001F(r)+  b1*bs*v_101F(r)
                 + c2eft * v_eftF(r);
    return (func);		
}	


local real sigma2par12F(real r, real b1, real b2, real bs, real c1eft, real c2eft, real sigma2eft)
{
    real func;
    func = spar_00F(r) + b1* spar_10F(r) + b1*b1 * spar_20F(r)
                 + b2*spar_01F(r)  + bs*spar_001F(r);
    return (func);	
}	


local real sigma2perp12F(real r, real b1, real b2, real bs, real c1eft, real c2eft, real sigma2eft)
{
    real func;
    func = sperp_00F(r) + b1* sperp_10F(r) + b1*b1 * sperp_20F(r)
                 + b2*sperp_01F(r)  + bs*sperp_001F(r);
    return (func);		
}	





local  real Interpolation_nr(real k, double kPS[], double pPS[], int nPS, double pPS2[])
{
    real psftmp;   
    splint(kPS,pPS,pPS2,nPS,k,&psftmp);  
    return (psftmp);
}
*/
