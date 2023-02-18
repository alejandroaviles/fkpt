/*==============================================================================
 NAME: clpt.c				[code for redshift space correlation function - GSM]
 Alejandro Aviles (avilescervantes@gmail.com), ...
 * 
================================================================================ 
 * 
 * We integrate over y and \alpha
 * 
 * with y=|vec-y| :
 * vec-y= vec-q - vec-r;
 * alpha = hat-y . hat-r; cosine angle between vec-y and vec-r
 * 
 * relations to more usual q and mu (mu = hat-q . hat-r :   cosine angle between vec-q and vec-r ):
 *  q       = sqrt(y^2 + r^2 + 2.*r*y*alpha )
 *  mu    = (r + alpha * y) / q
 *          and
 *  y       = sqrt( q^2 + r^2  - 2*q*r*mu )
 *  alpha = ( - r + mu*q )   / y
 * 
 ================================================================================ 
*/

#include "clpt.h"


#define CLPT_OUT_FMT    \
"%e %e %e %e %e %e %e %e \
 %e %e %e %e %e %e %e %e \
 %e %e %e %e %e %e %e %e \
 %e %e %e %e %e %e %e %e %e %e\n"


global void compute_clpt(void)
{
    stream outstr;
    //~ pointqfunctionsTableptr p;
    real aTime;
    real dr, ri;
    real placeholder;
    int i, Nr;
    real sigma2v,sigmav, f, f2;
    real rmin, rmax, nr;

    sigmav = rsqrt(sigma2L_function());    
    gd.particles_meanPath = 2.0*sigmav;
    
    rmin=golists.rmin;
    rmax= golists.rmax;
	Nr = golists.Nr;
	
    
    global_zacorrfunctions zacorrfun;
    global_clptcorrfunctions clptcorrfun;
    
   

    
    fprintf(stdout,"\nclpt: ");

    
    f = f_growth_LCDM(); 
    gd.f=f;
    f2=f*f;   
    
    fprintf(stdout," Nr=%d values from rmin=%g to rmax=%g ", Nr, rmin, rmax);

    outstr = stropen(gd.fpfnameclptfunctions,"w!");

   fprintf(outstr,"# InputPklFile=%s, redshift z=%g, OmegaM=%g, h=%g, Linear growth rate  f=%g\n",
                cmd.fnamePS, cmd.xstop, cmd.om, cmd.h, f);       
               
    fprintf(outstr,"# Precision: q-functions: NqperLogDecade=%d, Nk_qFunctionsQuad=%d\n",
                cmd.NqperLogDecade,cmd.Nk_qFunctionsQuad); 


    fprintf(outstr,"%3s%12s%16s%16s%16s%18s%18s%18s%28s%20s%20s%20s%20s%20s%16s%16s%16s \
    %19s%19s%19s%19s%19s%19s%24s%24s%24s%24s%24s%24s%24s%24s%24s%24s%24s%24s",
            "#","1.r[Mpc/h]","2.xiCLPT[1]","3.xi_10[b1]","4.xi_20[b1^2]","5.xi_01[b2] ",
            "6.xi_02[b2^2]","7.xi_11[b1*b2]","8.xi_eft[alpha1eft]",
            "9.xi_001[bs]","10.xi_002[bs^2] ","11.xi_101[b1*bs]","12.xi_011[b2*bs]",
            "13.xi_L(linearCF)  ","14.xi_ZA","15.xi_A", "16.xi_W","17.v12_00[1]","18.v12_10[b1]",
            "19.v12_20[b1^2]","20.v12_01[b2]","21.v12_11[b1*b2]","22.v12_001[bs]","23.v12_101[b1*bs]","24.v12_eft[alpha2eft]",
            "25.s12par_00[1]","26.s12par_10[b1]","27.s12par_20[b1^2]","28.s12par_01[b2]","29.s12par_001[bs]",
            "30.s12perp_00[1]","31.s12perp_10[b1]","32.s12perp_20[b1^2]","33.s12perp_01[b2]","34.s12perp_001[bs]\n");

//~ 8.xi_eft[alphaeft1] = 8.xi_nabla[2(1+b1)bnabla]


    placeholder = 0.0;
    aTime = second();
    


    
    dr = (rmax - rmin)/((real)(Nr - 1));
    
    for (i=1; i<=Nr; i++) {
        ri = rmin + dr*((real)(i - 1));
        zacorrfun = zacorrelation_functions(ri, sigmav);
        clptcorrfun = clptcorrelation_functions(ri, sigmav);
        
        //~ One_plus_xiZA=1+ zacorrfun.xi;            

        rArrays.rTab[i-1]=zacorrfun.r;
		rArrays.xi_00T[i-1]=zacorrfun.xi+clptcorrfun.xiA+clptcorrfun.xiW;
		rArrays.xi_10T[i-1]=clptcorrfun.xi10;
		rArrays.xi_20T[i-1]=clptcorrfun.xi20;
		rArrays.xi_01T[i-1]=clptcorrfun.xi01;
		rArrays.xi_02T[i-1]=clptcorrfun.xi02;
		rArrays.xi_11T[i-1]=clptcorrfun.xi11;
		rArrays.xi_eftT[i-1]=clptcorrfun.nabla2xi;
		rArrays.xi_001T[i-1]=clptcorrfun.xi001;
		rArrays.xi_002T[i-1]=clptcorrfun.xi002;
		rArrays.xi_101T[i-1]=clptcorrfun.xi101;
		rArrays.xi_011T[i-1]=clptcorrfun.xi011;         
		rArrays.xi_LT[i-1]=iBessel0F(ri);
		rArrays.xi_zaT[i-1]=zacorrfun.xi;
		rArrays.xi_AT[i-1]=clptcorrfun.xiA; 
		rArrays.xi_WT[i-1]=clptcorrfun.xiW; 
		rArrays.v_00T[i-1]= clptcorrfun.v12_00 * f;
		rArrays.v_10T[i-1]= clptcorrfun.v12_10 * f;
		rArrays.v_20T[i-1]= clptcorrfun.v12_20 * f;
		rArrays.v_01T[i-1]= clptcorrfun.v12_01 * f;
		rArrays.v_11T[i-1]= clptcorrfun.v12_11 * f;
		rArrays.v_001T[i-1]=  clptcorrfun.v12_001 * f;
		rArrays.v_101T[i-1]=  clptcorrfun.v12_101 * f;
		rArrays.v_eftT[i-1]=clptcorrfun.v12_eft * f;         
		rArrays.spar_00T[i-1]=clptcorrfun.s12par_00 * f2;
		rArrays.spar_10T[i-1]= clptcorrfun.s12par_10 * f2;
		rArrays.spar_20T[i-1]= clptcorrfun.s12par_20 * f2;
		rArrays.spar_01T[i-1]=clptcorrfun.s12par_01 * f2;
		rArrays.spar_001T[i-1]=clptcorrfun.s12par_001 * f2;           
		rArrays.sperp_00T[i-1]=clptcorrfun.s12perp_00 * f2;
		rArrays.sperp_10T[i-1]= clptcorrfun.s12perp_10 * f2;
		rArrays.sperp_20T[i-1]= clptcorrfun.s12perp_20 * f2;
		rArrays.sperp_01T[i-1]=clptcorrfun.s12perp_01 * f2;
		rArrays.sperp_001T[i-1]=clptcorrfun.s12perp_001 * f2;  
		
		
		
        fprintf(outstr,CLPT_OUT_FMT,
                zacorrfun.r,
                //xi
                zacorrfun.xi+clptcorrfun.xiA+clptcorrfun.xiW,   // \xi_matter CLPT [1]
                clptcorrfun.xi10,  // [b1]
                clptcorrfun.xi20,  // [b1^2]
                clptcorrfun.xi01,  // [b2]
                clptcorrfun.xi02, // [b2^2]
                clptcorrfun.xi11,  // [b1*b2]
                clptcorrfun.nabla2xi,  // [alpha1eft]~[2(1+b1) bnabla]
                clptcorrfun.xi001,    // [bs]
                clptcorrfun.xi002,   // [bs^2]
                clptcorrfun.xi101,     // [b1*bs]
                clptcorrfun.xi011,     // [b2 * bs]
                // 
                iBessel0F(ri),  // This is \xi_L
                zacorrfun.xi,   // Zeldovich approximation
                clptcorrfun.xiA,
                clptcorrfun.xiW,  
                //v12
                clptcorrfun.v12_00 * f,  // [1]
                clptcorrfun.v12_10 * f,   // [b1]
                clptcorrfun.v12_20 * f,   // [b1^2]
                clptcorrfun.v12_01 * f,  // [b2]
                clptcorrfun.v12_11 * f,    // [b1*b2]
                clptcorrfun.v12_001 * f,    // [bs]
                clptcorrfun.v12_101 * f,     // [b1*bs]
                clptcorrfun.v12_eft * f,  // [alpha2eft]
                //s12 parallel
                clptcorrfun.s12par_00 * f2,   // [1]
                clptcorrfun.s12par_10 * f2,  // [b1]
                clptcorrfun.s12par_20 * f2,   // [b1^2]
                clptcorrfun.s12par_01 * f2,   // [b1*b2]
                clptcorrfun.s12par_001 * f2,     // [bs]
                //s12 perpendicular
                clptcorrfun.s12perp_00 * f2,   // [1]
                clptcorrfun.s12perp_10 * f2,    // [b1]
                clptcorrfun.s12perp_20 * f2,   // [b1^2]
                clptcorrfun.s12perp_01 * f2,   // [b1*b2]
                clptcorrfun.s12perp_001 * f2     // [bs]
                );
    }
    fclose(outstr);
    fprintf(stdout,"...time = %g seconds",second()-aTime); 
    
    
}
#undef CLPT_OUT_FMT


// correlation functions

/*
 * We integrate over d^3 y  with vec y= vec q- vec r
 * y = |vec y| =  sqrt( q^2 + r^2 - 2 r q mu), 
 * with mu the angle between vec q and vec r
 * alpha is the angle between vec y and vec r, 
 * 
 *  q = sqrt(y^2 + r^2 + 2.*r*y*alpha )
 * mu = (r + alpha * y) / q 
 * 
 * We integrate over y and \alpha
*/
global_zacorrfunctions zacorrelation_functions(real r, real sigmav)
{
    int i, j;
//
    real *alphaGL, *wGL;
    int Nangle;
    real ymin, ymax, dy, y, q;
    int Ny;
    real xip, xiaA, xiaB;
    real mu, w, alpha;
    real fxf,hyf;
    real qminusmur, y2;
    real sigmapsi;
//
    global_zacorrfunctions_ptr zacorrfuncp;
    
    zacorrfuncp = (global_zacorrfunctions_ptr) allocate(1 * sizeof(global_zacorrfunctions));
    
    ymin = 0.001;
    ymax = 8*sigmav;
    Ny = 100;
    dy = (ymax - ymin)/((real)(Ny - 1));
    
    Nangle=12;
    alphaGL=dvector(1,Nangle);
    wGL=dvector(1,Nangle);
    gauleg(-1.0,1.0,alphaGL,wGL,Nangle);
    
    xip = 0.0; xiaA = 0.0; xiaB = 0.0;
    for (i=1; i<=Ny; i++) {
        y  = ymin + dy*((real)(i - 1));
        y2 = y*y;
        
        for (j=1; j<=Nangle; j++) {         
           alpha = alphaGL[j];
           w = wGL[j];
           
           q = rsqrt(y2 + r*r + 2.*r*y*alpha );
           mu = (r + alpha * y) / q ;
           qminusmur = q - mu * r;
           
           fxf = fXF(q); 
           hyf = hYF(q);
           
           xiaB += w * MZAF(y2, qminusmur, fxf, hyf)*y2 ;
        }
        if(i==1) xiaA = xiaB;
        xip += dy*(xiaA + xiaB)/2.0;
        xiaA = xiaB;
        xiaB = 0.0;
//
    }
    xip += -1.0;
//
    rzacorrfun(zacorrfuncp)    = r;
    xizacorrfun(zacorrfuncp)   = xip;
    
    free_dvector(wGL,1,Nangle);
    free_dvector(alphaGL,1,Nangle);
    
    return *zacorrfuncp;
}

/*
 * We integrate over d^3 y  with vec y= vec q- vec r
 * y = |vec y| =  sqrt( q^2 + r^2 - 2 r q mu), 
 * with mu the angle between vec q and vec r
 * alpha is the angle between vec y and vec r, 
 * 
 *  q = sqrt(y^2 + r^2 + 2.*r*y*alpha )
 * mu = (r + alpha * y) / q 
 * 
 * We integrate over y and \alpha
*/
global_clptcorrfunctions clptcorrelation_functions(real r, real sigmav)
{
    int i, j;
//
    real *alphaGL, *wGL;
    int NGL;
    real ymin, ymax, dy, y;
    int Ny;
    // xi
    real xiAp=0.0, xiAA=0.0, xiAB=0.0;
    real xiWp=0.0, xiWA=0.0, xiWB=0.0;
    real xi10p = 0.0, xi10A = 0.0, xi10B = 0.0;
    real xi20p = 0.0, xi20A = 0.0, xi20B = 0.0;
    real xi01p = 0.0, xi01A = 0.0, xi01B = 0.0;
    real xi02p = 0.0, xi02A = 0.0, xi02B = 0.0;
    real xi11p = 0.0, xi11A = 0.0, xi11B = 0.0;
    real nabla2xip = 0.0, nabla2xiA = 0.0, nabla2xiB = 0.0;    
    real xi001p=0.0, xi001A=0.0, xi001B=0.0;
    real xi002p=0.0, xi002A=0.0, xi002B=0.0;
    real xi101p=0.0, xi101A=0.0, xi101B=0.0;
    real xi011p=0.0, xi011A=0.0, xi011B=0.0;
    //v12
	real v12_00p=0.0, v12_00A=0.0, v12_00B=0.0;
	real v12_10p=0.0, v12_10A=0.0, v12_10B=0.0;
	real v12_20p=0.0, v12_20A=0.0, v12_20B=0.0;
	real v12_01p=0.0, v12_01A=0.0, v12_01B=0.0;
	real v12_11p=0.0, v12_11A=0.0, v12_11B=0.0;   
	real v12_001p=0.0, v12_001A=0.0, v12_001B=0.0;
	real v12_101p=0.0, v12_101A=0.0, v12_101B=0.0; 
	real v12_eftp=0.0, v12_eftA=0.0, v12_eftB=0.0;
	// sigma12 perpendicular
	real s12perp_00p=0.0,s12perp_00A=0.0,s12perp_00B=0.0;	
	real s12perp_10p=0.0,s12perp_10A=0.0,s12perp_10B=0.0;	
	real s12perp_20p=0.0,s12perp_20A=0.0,s12perp_20B=0.0;	
	real s12perp_01p=0.0,s12perp_01A=0.0,s12perp_01B=0.0;	
	real s12perp_001p=0.0,s12perp_001A=0.0,s12perp_001B=0.0;
	// sigma12 parallel
	real s12par_00p=0.0,s12par_00A=0.0,s12par_00B=0.0;	
	real s12par_10p=0.0,s12par_10A=0.0,s12par_10B=0.0;	
	real s12par_20p=0.0,s12par_20A=0.0,s12par_20B=0.0;	
	real s12par_01p=0.0,s12par_01A=0.0,s12par_01B=0.0;	
	real s12par_001p=0.0,s12par_001A=0.0,s12par_001B=0.0;
		
    
    real mu, alpha, w, mza, wmza;
    
    
    real dotXLf, dotXloopf, dotYLf, dotYloopf, ddotXLf, ddotXloopf, ddotYLf, ddotYloopf;
	real dotULf, dotUloopf, dotU20f, dotU11f, dotX10f, dotY10f, ddotX10f,ddotY10f; 
	real dotTf, dotV1f, dotV3f, ddotTf, ddotV1f, ddotV3f, dotV10f, ddotVbf;
    
   

    real y2, q, qminusmur, qigif, qiqjGijf, AijGij_xf, AijGij_yf;
    real fplush, fplushmu2, rigif, qirjGijf, rirjGijf, aux1, aux3, TrG;
    
    real XLf, YLf, Xloopf, Yloopf, Vf, Tf, ULf, Uloopf, X10f, Y10f, xiLf, U11f, U20f, nabla2xif;
    real dotVaf, dotVbf,ddotVaf,V10f,ij0f,ij2f,ij4f;
    real fxf, hyf;
    real V12f, Xupsilonf, Yupsilonf, chi12f, zetaf;
    real rcut=210.0; // r> rcut some functions are zero

	real auxddotAL, auxddotAloop, auxdotALdotAL, auxdotALdotUL, auxddotW, \
            auxULddotAL, auxddotA10, auxdotULdotUL, auxxiLddotAL;
    
    global_clptcorrfunctions_ptr clptcorrfunp;
    //~ global_v12_ptr v12funp;
    
    clptcorrfunp = (global_clptcorrfunctions_ptr) allocate(1 * sizeof(global_clptcorrfunctions));
    //~ v12funp =  (global_v12_ptr) allocate(1 * sizeof(global_v12));
 
 
    
    ymin = 0.001;
    ymax = 8*sigmav;
    Ny = 50;
	dy = (ymax - ymin)/((real)(Ny - 1));
        

    NGL=12;
    alphaGL=dvector(1,NGL);
    wGL=dvector(1,NGL);
    gauleg(-1.0,1.0,alphaGL,wGL,NGL);
    
    //~ suffix=_Ny200_NGL12
    //~ suffix=_Ny50_NGL12
    //~ suffix=_Ny800_NGL64    
    
    for (i=1; i<Ny; i++) {
        y  = ymin + dy*((real)(i - 1));
        y2 = y*y;
       

        
        for (j=1; j<=NGL; j++) {      
           alpha = alphaGL[j];
           w = wGL[j];
           
           q = rsqrt(y2 + r*r + 2.*r*y*alpha );
           mu = (r + alpha * y) / q ;
           qminusmur = q - mu * r;
                    
			XLf    = XLF(q); 
			YLf    = YLF(q);
			Xloopf    = XloopF(q); 
			Yloopf    = YloopF(q);
			Tf        = TF(q); 
			Vf        = preVF(q)-1.0/5.0 * Tf;
			ULf  = ULF(q);
			Uloopf    = UloopF(q);
			X10f      = X10F(q); 
			Y10f      = Y10F(q);
			xiLf = iBessel0F(q);
			U11f      = U11F(q);        
			U20f      = U20F(q);
			nabla2xif    = nabla2xiF(q);  // nabla2 bias
			
			dotVaf = dotVaF(q);
			dotVbf = dotVbF(q);
			ddotVaf = ddotVaF(q);
			V10f = V10F(q);
			ij0f=xiLf;
			ij2f=iBessel2F(q);
			ij4f=iBessel4F(q);
			
			// Derived for tidal
			
			zetaf = 4.* ( 7.*rpow(ij0f,2.) + 10.*rpow(ij2f,2.) + 18.*rpow(ij4f,2.) ) / 315. ;  
			Xupsilonf  = q*q * 4. *  rpow(7.*ij0f + 16.*ij2f + 9.*ij4f, 2.) / 99225.;    
			Yupsilonf  =  q*q * 4.* ( -28.* ij0f * ij2f - 18.0 * ( 7. * ij0f + ij2f ) * ij4f + 49. * rpow(ij0f,2.) 
			                      - 23. * rpow(ij2f,2.) + 54. * rpow(ij4f,2.) )  /  33075.;   
			V12f    = q * 4.* ij2f * ( 14.*ij0f + 5.*ij2f - 9.*ij4f ) / 315.;  // WRONG! CHECK!
			chi12f = 4.*rpow(ij2f,2.) / 3.;  
			
			
			
			
			// Derived for s12 and v12
			
			dotXLf    = XLf;
			dotYLf    = YLf;
			ddotXLf = XLf;
			ddotYLf = YLf;
			dotXloopf     = 2.0*Xloopf ;
			dotYloopf     = 2.0*Yloopf;
			ddotXloopf  = 4.0*Xloopf ;
			ddotYloopf  = 4.0*Yloopf;

			dotULf       = ULf;
			dotUloopf = 3.0*Uloopf;
			dotU20f    = 2.0*U20f;
			dotU11f    = 2.0*U11f;

			dotX10f     = 3.0/2.0 * X10f; 
			dotY10f     = 3.0/2.0 * Y10f;
			ddotX10f  = 2.0 * X10f;
			ddotY10f  = 2.0 * Y10f;

			dotTf       = 4.0/3.0 * Tf;
			dotV1f    = dotVaf - 1.0/5.0 * dotTf;
			dotV3f    = dotVaf + dotVbf - 1.0/5.0 * dotTf;
			ddotTf    = 5.0/3.0 * Tf;
			ddotV1f = ddotVaf - 1.0/5.0 * ddotTf;
			ddotVbf = dotVbf;
			ddotV3f = ddotVaf + ddotVbf - 1.0/5.0 * ddotTf;

			dotV10f = -2.0 * V10f;
			
			
			
			//	
            fxf= 1.0/XLf;
		    hyf = -YLf / ( XLf*XLf+ XLf*YLf );
            fplush= fxf + hyf;
            fplushmu2= fxf + hyf*mu*mu;
            qigif = qminusmur* fplush;
            qiqjGijf = fplush - qigif * qigif ; 
			rigif = fplush * mu * q  -  fplushmu2 * r ;
			qirjGijf =fplush *mu- qigif * rigif;
			rirjGijf =fplushmu2    -     rigif * rigif;
			aux1= q * fplush - r * mu * hyf;
			TrG = 3.0 * fxf + hyf - (r * fxf)*(r * fxf) - rpow( aux1 ,2.0 ) + 2.0 * r * fxf* aux1 * mu;        
            AijGij_xf  = -(rsqr(fxf)*y2 + hyf*(-1.0 + hyf*rsqr(qminusmur))
								+ fxf*(-3.0 + 2.0*hyf*rsqr(qminusmur)));
            AijGij_yf = - fplush*(-1.0 + fplush*rsqr(qminusmur) );
            
            
            
            wmza = w * MZAF(y2, qminusmur, fxf, hyf) * y2;            
            
            
            
            if (r<rcut){ 
            	xiAB  += wmza *  ( - 0.5  * ( AijGij_xf*Xloopf + AijGij_yf*Yloopf) );
            	xiWB += wmza *  ( - 1./6. *  fplush*(qminusmur) * (
   				          fplush*(-3.0 + fplush*rsqr(qminusmur))*Tf
										+ 3.0*(rsqr(fxf)*y2 + hyf*(-3.0 + hyf*rsqr(qminusmur))
										+ fxf*(-5.0 + 2.0*hyf*rsqr(qminusmur)))*Vf)    
                                );
            };
            xi10B += wmza  *  (-2.0 * (ULf +Uloopf) * qigif -   AijGij_xf*X10f - AijGij_yf*Y10f ); 	
            xi20B += wmza  *  ( xiLf - ULf*ULf * qiqjGijf - U11f * qigif    );
            if (r<rcut){            
            	xi01B     += wmza  *  ( -ULf*ULf * qiqjGijf - U20f * qigif    );  
            	xi02B     += wmza  *  ( 0.5 * xiLf * xiLf     );      
            	xi11B     += wmza  *  ( -2.0 * ULf * xiLf *qigif     );  
            	nabla2xiB    += wmza  *  ( nabla2xif    );   
            };
			
			// Tidal bias
				xi001B += wmza  * (- AijGij_xf*Xupsilonf - AijGij_yf*Yupsilonf  - 2.0 * V10f * qigif );  // WRONG! CHECK!
				xi002B += wmza  * zetaf;
				xi101B += wmza  * ( -2.0 * V12f * qigif);
				xi011B += wmza  * chi12f;
			
		//v12
			aux3          =  -2.0 * ULf * dotULf * qigif * mu;

			v12_00B +=  wmza * ( -(XLf +dotXloopf)   * rigif      - (YLf + dotYloopf) * qigif *  mu  
										-dotV1f * qirjGijf - 0.5 * dotV3f * mu * TrG  -  0.5* dotTf * qiqjGijf * mu);
			v12_10B +=  wmza * ( 2.0* (dotULf + dotUloopf ) * mu + 2.0 * (-dotX10f * rigif  -   dotY10f * qigif * mu) 
										-2.0*ULf  *  (dotXLf * qirjGijf + dotYLf * qiqjGijf* mu) );
			v12_20B +=  wmza * ( dotU11f * mu  +aux3 + xiLf * (-dotXLf* rigif - dotYLf * qigif * mu) );   
			v12_01B +=  wmza * ( dotU20f * mu  + aux3);    
			v12_11B +=  wmza * ( 2.0 * xiLf * dotULf * mu );    
		
		
			v12_001B +=  wmza * ( 2.0 *  dotV10f * mu  -Xupsilonf * rigif  -   Yupsilonf * qigif * mu);  //Check!    
			v12_101B +=  wmza * ( 2.0 *  V12f * mu );  //Check!   
			
		// *** sigma12 perpendicular ***

		//CDM
		auxddotAL         = ddotXLf + 0.5* ddotYLf*(1. - mu*mu);
		auxddotAloop   = ddotXloopf + 0.5* ddotYloopf*(1. - mu*mu);
		auxddotW          = -ddotV1f * qigif - ddotV3f*( qigif - mu * rigif) - 0.5*ddotTf*qigif*(1. - mu*mu);
		auxdotALdotAL = -0.5*dotXLf*dotXLf *(TrG - rirjGijf) -dotXLf*dotYLf*(qiqjGijf - qirjGijf*mu) 
										-0.5*dotYLf*dotYLf*qiqjGijf *(1. - mu*mu);

		//b1
		auxdotALdotUL = -2.0* dotULf *(dotXLf*(  qigif - mu * rigif) +  dotYLf *(1. - mu*mu) * qigif);
		auxULddotAL = auxddotAL*(-2.0*dotULf*qigif);
		auxddotA10 = 2.0*ddotX10f + ddotY10f * (1.0 - mu*mu);

		//b1^2
		auxdotULdotUL = dotULf *dotULf*(1. - mu*mu); //this is also for b2
		auxxiLddotAL = auxddotAL*xiLf;


		s12perp_00B += wmza * (auxddotAL + auxddotAloop + auxddotW + auxdotALdotAL);
		s12perp_10B += wmza * (auxdotALdotUL + auxULddotAL + auxddotA10);
		s12perp_20B += wmza * (auxdotULdotUL + auxxiLddotAL);     //[b1^2]
		s12perp_01B += wmza * (auxdotULdotUL);	
		s12perp_001B += wmza * (Xupsilonf + 0.5* Yupsilonf*(1. - mu*mu) );

			
		// *** sigma12 parallel ***


		// CDM
		auxddotAL = ddotXLf + ( ddotYLf)*mu*mu;
		auxddotAloop = ddotXloopf + ( ddotYloopf)*mu*mu;
		auxddotW = -ddotV1f * qigif -  2.0*ddotV3f*mu * rigif  -  ddotTf*qigif*mu*mu;
		auxdotALdotAL = -dotXLf*dotXLf * rirjGijf   - 2.0*dotXLf*dotYLf*qirjGijf*mu 
										- dotYLf*dotYLf*qiqjGijf* mu*mu;

		//b1
		auxdotALdotUL = -4.0* dotULf *(dotXLf * mu * rigif + dotYLf *  mu*mu * qigif);
		auxULddotAL = auxddotAL*(-2.0*dotULf*qigif);
		auxddotA10 = 2.0*(ddotX10f + ddotY10f * mu*mu);

		//b1^2
		auxdotULdotUL = 2.0*dotULf *dotULf*mu*mu;  //this is also for b2
		auxxiLddotAL = auxddotAL*xiLf;


		s12par_00B += wmza * (auxddotAL + auxddotAloop + auxddotW + auxdotALdotAL);
		s12par_10B += wmza * (auxdotALdotUL + auxULddotAL + auxddotA10);
		s12par_20B += wmza * (auxdotULdotUL + auxxiLddotAL);
		s12par_01B += wmza * (auxdotULdotUL);
		s12par_001B += wmza * (Xupsilonf +  Yupsilonf*mu*mu);
			
        }
        
        xiAp  += dy*(xiAA + xiAB)/2.0; xiAA = xiAB; xiAB = 0.0;
        xiWp += dy*(xiWA + xiWB)/2.0; xiWA = xiWB; xiWB = 0.0;
        xi10p += dy*(xi10A + xi10B)/2.0; xi10A = xi10B; xi10B = 0.0;
        xi20p += dy*(xi20A + xi20B)/2.0; xi20A = xi20B; xi20B = 0.0;
        xi01p += dy*(xi01A + xi01B)/2.0; xi01A = xi01B; xi01B = 0.0;
        xi02p += dy*(xi02A + xi02B)/2.0; xi02A = xi02B; xi02B = 0.0;
        xi11p += dy*(xi11A + xi11B)/2.0; xi11A = xi11B; xi11B = 0.0;
        nabla2xip += dy*(nabla2xiA + nabla2xiB)/2.0; nabla2xiA = nabla2xiB; nabla2xiB = 0.0;
        
        v12_00p += dy*(v12_00A + v12_00B)/2.0; v12_00A = v12_00B; v12_00B = 0.0;
        v12_10p += dy*(v12_10A + v12_10B)/2.0; v12_10A = v12_10B; v12_10B = 0.0;
        v12_20p += dy*(v12_20A + v12_20B)/2.0; v12_20A = v12_20B; v12_20B = 0.0;
        v12_01p += dy*(v12_01A + v12_01B)/2.0; v12_01A = v12_01B; v12_01B = 0.0;
        v12_11p += dy*(v12_11A + v12_11B)/2.0; v12_11A = v12_11B; v12_11B = 0.0;
        v12_eftp += dy*(v12_eftA + v12_eftB)/2.0; v12_eftA = v12_eftB; v12_eftB = 0.0;
        
         s12perp_00p += dy*(s12perp_00A + s12perp_00B)/2.0; s12perp_00A = s12perp_00B; s12perp_00B = 0.0;
         s12perp_10p += dy*(s12perp_10A + s12perp_10B)/2.0; s12perp_10A = s12perp_10B; s12perp_10B = 0.0;
         s12perp_20p += dy*(s12perp_20A + s12perp_20B)/2.0; s12perp_20A = s12perp_20B; s12perp_20B = 0.0;
         s12perp_01p += dy*(s12perp_01A + s12perp_01B)/2.0; s12perp_01A = s12perp_01B; s12perp_01B = 0.0;
         s12perp_001p += dy*(s12perp_001A + s12perp_001B)/2.0; s12perp_001A = s12perp_001B; s12perp_001B = 0.0;
         
         s12par_00p += dy*(s12par_00A + s12par_00B)/2.0; s12par_00A = s12par_00B; s12par_00B = 0.0;
         s12par_10p += dy*(s12par_10A + s12par_10B)/2.0; s12par_10A = s12par_10B; s12par_10B = 0.0;
         s12par_20p += dy*(s12par_20A + s12par_20B)/2.0; s12par_20A = s12par_20B; s12par_20B = 0.0;
         s12par_01p += dy*(s12par_01A + s12par_01B)/2.0; s12par_01A = s12par_01B; s12par_01B = 0.0;
         s12par_001p += dy*(s12par_001A + s12par_001B)/2.0; s12par_001A = s12par_001B; s12par_001B = 0.0;
       

			xi001p += dy*(xi001A + xi001B)/2.0; xi001A = xi001B; xi001B = 0.0;
			xi002p += dy*(xi002A + xi002B)/2.0; xi002A = xi002B; xi002B = 0.0;
			xi101p += dy*(xi101A + xi101B)/2.0; xi101A = xi101B; xi101B = 0.0;
			xi011p += dy*(xi011A + xi011B)/2.0; xi011A = xi011B; xi011B = 0.0;
			v12_001p += dy*(v12_001A + v12_001B)/2.0; v12_001A = v12_001B; v12_001B = 0.0;
            v12_101p += dy*(v12_101A + v12_101B)/2.0; v12_101A = v12_101B; v12_101B = 0.0;
			s12par_001p += dy*(s12par_001A + s12par_001B)/2.0; s12par_001A = s12par_001B; s12par_001B = 0.0;
			s12perp_001p += dy*(s12perp_001A + s12perp_001B)/2.0; s12perp_001A = s12perp_001B; s12perp_001B = 0.0;
             
        
        
//
    }
    
    rclptcorrfun(clptcorrfunp) = r;
    xiAclptcorrfun(clptcorrfunp) = xiAp;
    xiWclptcorrfun(clptcorrfunp) = xiWp;
    xi10clptcorrfun(clptcorrfunp) = xi10p;
    xi20clptcorrfun(clptcorrfunp) = xi20p;
    xi01clptcorrfun(clptcorrfunp) = xi01p;
    xi02clptcorrfun(clptcorrfunp) = xi02p;
    xi11clptcorrfun(clptcorrfunp) = xi11p;
    nabla2xiclptcorrfun(clptcorrfunp) = nabla2xip;
    xi001clptcorrfun(clptcorrfunp) = xi001p;
    xi002clptcorrfun(clptcorrfunp) = xi002p;
    xi101clptcorrfun(clptcorrfunp) = xi101p;
    xi011clptcorrfun(clptcorrfunp) = xi011p;
      
    v12_00_vfun(clptcorrfunp)    = v12_00p;
	v12_10_vfun(clptcorrfunp)    = v12_10p;
	v12_20_vfun(clptcorrfunp)    = v12_20p;
	v12_01_vfun(clptcorrfunp)    = v12_01p;
	v12_11_vfun(clptcorrfunp)    = v12_11p;
	v12_001_vfun(clptcorrfunp)  = v12_001p;
	v12_101_vfun(clptcorrfunp)  = v12_101p;
	v12_eft_vfun(clptcorrfunp)    = v12_eftp;
     
    s12perp_00_sfun(clptcorrfunp)    = s12perp_00p;
    s12perp_10_sfun(clptcorrfunp)    = s12perp_10p;
    s12perp_20_sfun(clptcorrfunp)    = s12perp_20p;
    s12perp_01_sfun(clptcorrfunp)    = s12perp_01p;
    s12par_00_sfun(clptcorrfunp)    = s12par_00p;
    s12par_10_sfun(clptcorrfunp)    = s12par_10p;
    s12par_20_sfun(clptcorrfunp)    = s12par_20p;
    s12par_01_sfun(clptcorrfunp)    = s12par_01p;
    
    free_dvector(wGL,1,NGL);
    free_dvector(alphaGL,1,NGL);
    
    return *clptcorrfunp;
}



//// q-functions

local real XLF(real q)
{
    real func;
    func = Interpolation_nr(q, qArrays.qTab, qArrays.XLT, golists.Nq, qArraysd.XLT);
    return (func);
}

local real YLF(real q)
{
    real func;
    func = Interpolation_nr(q, qArrays.qTab, qArrays.YLT, golists.Nq, qArraysd.YLT);
    return (func);
}

local real XloopF(real q)
{
    real func;
    func = Interpolation_nr(q, qArrays.qTab, qArrays.XloopT, golists.Nq, qArraysd.XloopT);
    return (func);
}

local real YloopF(real q)
{
    real func;
    func = Interpolation_nr(q, qArrays.qTab, qArrays.YloopT, golists.Nq, qArraysd.YloopT);
    return (func);
}

local real preVF(real q)
{
    real func;
    func = Interpolation_nr(q, qArrays.qTab, qArrays.preVT, golists.Nq, qArraysd.preVT);
    return (func);
}

local real TF(real q)
{
    real func;
    func = Interpolation_nr(q, qArrays.qTab, qArrays.TT, golists.Nq, qArraysd.TT);
    return (func);
}

local real X10F(real q)
{
    real func;
    func = Interpolation_nr(q, qArrays.qTab, qArrays.X10T, golists.Nq, qArraysd.X10T);
    return (func);
}

local real Y10F(real q)
{
    real func;
    func = Interpolation_nr(q, qArrays.qTab, qArrays.Y10T, golists.Nq, qArraysd.Y10T);
    return (func);
}

local real ULF(real q)
{
    real func;
    func = Interpolation_nr(q, qArrays.qTab, qArrays.ULT, golists.Nq, qArraysd.ULT);
    return (func);
}

local real UloopF(real q)
{
    real func;
    func = Interpolation_nr(q, qArrays.qTab, qArrays.UloopT, golists.Nq, qArraysd.UloopT);
    return (func);
}

local real U11F(real q)
{
    real func;
    func = Interpolation_nr(q, qArrays.qTab, qArrays.U11T, golists.Nq, qArraysd.U11T);
    return (func);
}

local real U20F(real q)
{
    real func;
    func = Interpolation_nr(q, qArrays.qTab, qArrays.U20T, golists.Nq, qArraysd.U20T);
    return (func);
}

local real dotVaF(real q)
{
    real func;
    func = Interpolation_nr(q, qArrays.qTab, qArrays.dotVaT, golists.Nq, qArraysd.dotVaT);
    return (func);
}


local real ddotVaF(real q)
{
    real func;
    func = Interpolation_nr(q, qArrays.qTab, qArrays.ddotVaT, golists.Nq, qArraysd.ddotVaT);
    return (func);
}

local real dotVbF(real q)
{
    real func;
    func = Interpolation_nr(q, qArrays.qTab, qArrays.dotVbT, golists.Nq, qArraysd.dotVbT);
    return (func);
}


local real V10F(real q)
{
    real func;
    func = Interpolation_nr(q, qArrays.qTab, qArrays.V10T, golists.Nq, qArraysd.V10T);
    return (func);
}

local real iBessel0F(real q)
{
    real func;
    func = Interpolation_nr(q, qArrays.qTab, qArrays.iBessel0T, golists.Nq, qArraysd.iBessel0T);
    return (func);
}

local real iBessel2F(real q)
{
    real func;
    func = Interpolation_nr(q, qArrays.qTab, qArrays.iBessel2T, golists.Nq, qArraysd.iBessel2T);
    return (func);
}


local real iBessel4F(real q)
{
    real func;
    func = Interpolation_nr(q, qArrays.qTab, qArrays.iBessel4T, golists.Nq, qArraysd.iBessel4T);
    return (func);
}

local real nabla2xiF(real q)
{
    real func;
    func = Interpolation_nr(q, qArrays.qTab, qArrays.nabla2xiT, golists.Nq, qArraysd.nabla2xiT);
    return (func);
}

// DERIVED FUNCTIONS:

local real UF(real q)
{
    real func;
    func = ULF(q) + UloopF(q);
    return (func);
}

local real fXF(real q)
{
    real func;
    func = 1.0/XLF(q);
    return (func);
}

local real hYF(real q)
{
    real func;
    func = -YLF(q)/(XLF(q)*XLF(q)+ XLF(q)*YLF(q));
    return (func);
}

local real MZAF(real y2, real qminusmur, real fx, real hy)
{
    real func;
    func = 0.39894228040143267794 * fx * rsqrt(fx+hy)  
    * rexp( -0.5*( hy* rsqr(qminusmur) + fx*y2) );
     // 0.39894228040143267794=INVSQRTDTWOPI = 1/sqrt(2*pi)
    return (func);
}





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


local real sigma2L_function(void)
{
    real result;
    real kmin, kmax;
    real ymin, ymax;
    real EPSQ = 0.000001;
    int KK = 5;

    kmin = kPS[1];
    kmax = kPS[nPSLT];
    ymin = rlog10(kmin);
    ymax = rlog10(kmax);

    result= (1.0/SIXPI2)*rlog(10.0)
    *qromo(sigma2L_function_int,ymin,ymax,midpnt,EPSQ,KK);
    
    
    //~ fprintf(stdout,"\n sigma2v= %g \n",result);

    return result;

}

