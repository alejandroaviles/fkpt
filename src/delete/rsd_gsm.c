/*==============================================================================
 NAME: rsd.c				[code for redshift space correlation function - GSM]
 Alejandro Aviles (avilescervantes@gmail.com)
 ================================================================================ 
*/
#include "rsd.h"

#define USENEW 1

global void compute_gsm(void)
{
    stream outstr;
    pointrfunctionsTableptr p;
    real aTime;
    int counter;
    real r, nearlyzero, ds, si;
    real smin, smax;
    int Ns;
	real b1, b2, bs, c2eft, c1eft, s2eft;
		
	b1= cmd.b1;
	b2= cmd.b2;
	bs= cmd.bs;
	c1eft= cmd.c1eft;
	c2eft= cmd.c2eft;
	s2eft= cmd.s2eft;
	
	smin = golists.smin;
	smax = golists.smax;
	Ns = golists.Ns;

    global_rsdmultipoles rsd; 
    
    
    
    outstr = stropen(gd.fpfnamersd,"w!");


                
   fprintf(outstr,"# InputPklFile=%s, redshift z=%g, OmegaM=%g, h = %g, f=%g, b1=%g, b2=%g, bs=%g, c1eft=%g, c2eft=%g, s2FoG=%g\n",
                cmd.fnamePS, cmd.xstop, cmd.om, cmd.h, gd.f, cmd.b1, cmd.b2, cmd.bs, cmd.c1eft, cmd.c2eft,cmd.s2eft);       
               
    fprintf(outstr,"# Precision: q-functions: NqperLogDecade=%d, Nk_qFunctionsQuad=%d. gsm: gsm_width=%g, gsm_sizeyT=%d, gsm_NGL=%d\n",
                cmd.NqperLogDecade, cmd.Nk_qFunctionsQuad, cmd.gsm_width, cmd.gsm_sizeyT, cmd.gsm_NGL);                    

    fprintf(outstr,"%1s%12s%16s%16s%16s",
            "#","1.s[Mpc/h]","2.monopole","3.quadrupole","4.hexadecapole\n");
    
        
    fprintf(stdout,"\nrsd: ");
    fprintf(stdout," Ns=%d values from smin=%g to smax=%g ",
            Ns, smin, smax);
	aTime = second();	
    ds = (smax - smin)/((real)(Ns - 1));    
    for (int i=1; i<=Ns; i++) {
        si = smin + ds*((real)(i - 1));       
        rsd =CF_rsd_multipoles(si, b1, b2, bs, c1eft, c2eft, s2eft);  
		
		
        fprintf(outstr,"%e %e %e %e\n",
                rsd.s,  // [Mpc/h]
                rsd.xi_0,  // [monopole]
                rsd.xi_2,  // [quadrupole]
                rsd.xi_4  // [hexadecapole]
                );        
	};	
	fclose(outstr);
    fprintf(stdout,"...time = %g seconds \n",second()-aTime); 
 
}    

/*
we integrate over the variable y = r_parallel - s mu, 
 mu is the angle betwen vec-s and the line of sight 
 s^2 = r_parallel^2 + s_perpendicular^2 
 */
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







// (density weighted) velocity moments

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




////////////////////////////////////////////////////////////
///////   AUX FUNCTIONS
////////////////////////////////////////////////////////////



//2
local real xi_00F(real r)
{
    real func;
    //~ if (USENEW==0){
		//~ func = Interpolation_nr(r, rTab, xi_00T, nrfunctionsTable, xi_00T2);    
    //~ } else {
    //~ func = Interpolation_nr(r, rArrays.rTab, rArrays.xi_00T, golists.Nr, rArraysd.xi_00T); 
    func = Interpolation_nr(r, rArrays.rTab, rArrays.xi_00T, golists.Nr, rArraysd.xi_00T); 
	//~ };
    return (func);			
};
//3
local real xi_10F(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.xi_10T, golists.Nr, rArraysd.xi_10T);
    return (func);			
};
//4
local real xi_20F(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.xi_20T, golists.Nr, rArraysd.xi_20T);
    return (func);			
};
//5
local real xi_01F(real r)
{
    real func;
    //~ func = Interpolation_nr(r, rTab, xi_01T, golists.Nr, xi_01T2);
    func = Interpolation_nr(r, rArrays.rTab, rArrays.xi_01T, golists.Nr, rArrays.xi_01T);
    return (func);			
};
//6
local real xi_02F(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.xi_02T, golists.Nr, rArrays.xi_02T);
    return (func);			
};
//7
local real xi_11F(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.xi_11T, golists.Nr, rArrays.xi_11T);
    return (func);			
};
//8
local real xi_eftF(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.xi_eftT, golists.Nr, rArrays.xi_eftT);
    return (func);			
};
//9
local real xi_001F(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.xi_001T, golists.Nr, rArrays.xi_001T);
    return (func);			
};
//10
local real xi_002F(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.xi_002T, golists.Nr, rArrays.xi_002T);
    return (func);			
};
//11
local real xi_101F(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.xi_101T, golists.Nr, rArrays.xi_101T);
    return (func);			
};
//12
local real xi_011F(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.xi_011T, golists.Nr, rArrays.xi_011T);
    return (func);			
};

//v12:
// 17  
local real v_00F(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.v_00T, golists.Nr, rArraysd.v_00T);
    return (func);			
};
// 18  
local real v_10F(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.v_10T, golists.Nr, rArraysd.v_10T);
    return (func);			
};
// 19  
local real v_20F(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.v_20T, golists.Nr, rArraysd.v_20T);
    return (func);			
};
// 20  
local real v_01F(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.v_01T, golists.Nr, rArraysd.v_01T);
    return (func);			
};
// 21  
local real v_11F(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.v_11T, golists.Nr, rArraysd.v_11T);
    return (func);			
};
// 22  
local real v_001F(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.v_001T, golists.Nr, rArraysd.v_001T);
    return (func);			
};
// 23  
local real v_101F(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.v_101T, golists.Nr, rArraysd.v_101T);
    return (func);			
};
// 24  
local real v_eftF(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.v_eftT, golists.Nr, rArraysd.v_eftT);
    return (func);			
};

//sigma12 parallel:
// 25  
local real spar_00F(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.spar_00T, golists.Nr, rArraysd.spar_00T);
    return (func);			
};
// 26  
local real spar_10F(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.spar_10T, golists.Nr, rArraysd.spar_10T);
    return (func);			
};
// 27  
local real spar_20F(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.spar_20T, golists.Nr, rArraysd.spar_20T);
    return (func);			
};
// 28 
local real spar_01F(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.spar_01T, golists.Nr, rArraysd.spar_01T);
    return (func);			
};
// 29  
local real spar_001F(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.spar_001T, golists.Nr, rArraysd.spar_001T);
    return (func);			
};

//sigma12 perpendicular:
// 30  
local real sperp_00F(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.sperp_00T, golists.Nr, rArraysd.sperp_00T);
    return (func);			
};
// 31 
local real sperp_10F(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.sperp_10T, golists.Nr, rArraysd.sperp_10T);
    return (func);			
};
// 32 
local real sperp_20F(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.sperp_20T, golists.Nr, rArraysd.sperp_20T);
    return (func);			
};
// 33 
local real sperp_01F(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.sperp_01T, golists.Nr, rArraysd.sperp_01T);
    return (func);			
};
// 34 
local real sperp_001F(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.sperp_001T, golists.Nr, rArraysd.sperp_001T);
    return (func);			
};
// 35
local real xi_zaF(real r)
{
    real func;
    func = Interpolation_nr(r, rArrays.rTab, rArrays.xi_zaT, golists.Nr, rArraysd.xi_zaT);
    return (func);			
};



local  real Interpolation_nr(real k, double kPS[], double pPS[], int nPS, double pPS2[])
{
    real psftmp;   
    splint(kPS,pPS,pPS2,nPS,k,&psftmp);  
    return (psftmp);
}

