/*==============================================================================
 NAME: qfunctions.c				[code for redshift space correlation function - GSM]
 Alejandro Aviles (avilescervantes@gmail.com), ...
 ================================================================================ 
*/

#include "qfunctions.h"



global void compute_qfunctions(void)
{
    stream outstr;
    pointQsRsTableptr p;
    real aTime;
    real qval;
    real qi;
    real dq;
    real qmin, qmax, qsperLogDecade;
    int i, Nq;
     
    
    
    qmin = golists.qmin;     
    qmax = golists.qmax;    
    Nq = golists.Nq;   
       
    global_qfunctions qfun;
    global_corrfunctions corrfun;
    
    fprintf(stdout,"\nq-functions:");
    fprintf(stdout," Nq=%d values from qmin=%g to qmax=%g ",
            Nq, qmin, qmax);


    aTime = second();
    dq = (rlog10(qmax) - rlog10(qmin))/((real)(Nq - 1));      
    for (i=1; i<=Nq; i++) {
        qval = rlog10(qmin) + dq*((real)(i - 1));
        qi = rpow(10.0,qval);
        fflush(stdout);
        qfun = qfunctions(qi);
        //~ corrfun = correlation_functions(qi);

		qArrays.qTab[i-1] = qfun.q;
		qArrays.XLT[i-1] = qfun.XL;
		qArrays.YLT[i-1] = qfun.YL;
		qArrays.XloopT[i-1] = qfun.Xloop;
		qArrays.YloopT[i-1] = qfun.Yloop;
		qArrays.preVT[i-1] = qfun.preV;
		qArrays.TT[i-1] = qfun.T;
		qArrays.X10T[i-1] = qfun.X10;
		qArrays.Y10T[i-1] = qfun.Y10;
		qArrays.ULT[i-1] = qfun.UL;
		qArrays.UloopT[i-1] = qfun.Uloop;
		qArrays.U11T[i-1] = qfun.U11;
		qArrays.U20T[i-1] = qfun.U20;
		qArrays.dotVaT[i-1] = qfun.dotVa;
		qArrays.ddotVaT[i-1] = qfun.ddotVa;
		qArrays.dotVbT[i-1] = qfun.dotVb; 
		qArrays.V10T[i-1] = qfun.V10;
		qArrays.iBessel0T[i-1] = qfun.iBessel0;
		qArrays.iBessel2T[i-1] = qfun.iBessel2;
		qArrays.iBessel4T[i-1] = qfun.iBessel4;
		qArrays.nabla2xiT[i-1] = qfun.nabla2xi;                           
                             
    }
    
    fprintf(stdout,"...time = %g seconds ",second()-aTime);    
 
}






// q functions
global_qfunctions qfunctions(real qi)
{
    global_qfunctions_ptr qfunp;
    int i, Nk;
    real dk, kvali, kvalim1, ki, kim1;
    real kk, k2;
    real deltak;
    real kmaxh, kminh;
    
//
    real ULp, ULA, ULB;
    real Uloopp, UloopA, UloopB;
    real U11p, U11A, U11B;
    real U20p, U20A, U20B;
//
    real XLp, XLA, XLB;
    real Xloopp, XloopA, XloopB;
    real X10p, X10A, X10B;
    real YLp, YLA, YLB;
    real Yloopp, YloopA, YloopB;
    real Y10p, Y10A, Y10B;
//
    real Vp, preVp, preVA, preVB;
    real Tp, TA, TB;
    
    real dotVap, dotVaA, dotVaB;
    real ddotVap, ddotVaA, ddotVaB;
    real dotVbp, dotVbA, dotVbB;
    
    real V10p, V10A, V10B;    
    real xiell0p, xiell0A, xiell0B;
    real xiell2p, xiell2A, xiell2B;
    real xiell4p, xiell4A, xiell4B;
    
    real nabla2xip, nabla2xiA, nabla2xiB;
//    
    real j0, j1, j2, j3, j4, kq, damping;   
    real R1val, R2val, Q1val, Q2val, Q5val, Q8val, Qs2val, pklval; 

    
    qfunp = (global_qfunctions_ptr) allocate(1 * sizeof(global_qfunctions));

    //~ Nk = 1200; //original
    //~ Nk = 800;    //     
    Nk = cmd.Nk_qFunctionsQuad;


    //~ kmaxh=20.0;
    kmaxh = cmd.kmax;
    kminh = cmd.kmin; 
    dk = (rlog10(kmaxh) - rlog10(kminh))  /  ((real)(Nk - 1)) ;
    
    ULp = 0.;ULA = 0.; 
    Uloopp = 0.; UloopA = 0.;
    U11p = 0.;U11A = 0.; 
    U20p = 0.;U20A = 0;
//
    XLp = 0.;XLA = 0.;
    Xloopp = 0.;XloopA = 0.;
    X10p = 0.;X10A = 0;
    YLp = 0.;YLA = 0;
    Yloopp = 0.;YloopA = 0.;
    Y10p = 0.;Y10A = 0.;
//
    preVp = 0.;preVA = 0;
    Tp = 0.;TA = 0;
    
    dotVap=0, dotVaA=0;
    ddotVap=0, ddotVaA=0;
    dotVbp=0, dotVbA=0;
    
    V10p=0, V10A=0;
    xiell0p=0, xiell0A=0;
    xiell2p=0, xiell2A=0;
    xiell4p=0, xiell4A=0;
    
    nabla2xip=0, nabla2xiA=0;
//
    for (i=2; i<=Nk; i++) {
        kvali = rlog10(cmd.kmin) + dk*((real)(i - 1));
        kvalim1 = rlog10(cmd.kmin) + dk*((real)(i - 2));
        ki = rpow(10.0,kvali);
        kim1 = rpow(10.0,kvalim1);
        deltak = (ki - kim1);
        kk = rpow(10.0,kvali);
        
        k2=kk*kk;
        kq  = kk * qi;
        damping = rexp(-k2); // anti-aliasing
        
        j0 =rj0Bessel(kq);
        j1 =rj1Bessel(kq);
        j2 =rj2Bessel(kq);
        j3 =rj3Bessel(kq);
        j4 =rj4Bessel(kq);
        
        R1val  =R1F(kk);
        R2val  =R2F(kk);
        Q1val  =Q1F(kk);
        Q2val  =Q2F(kk);
        Q5val  =Q5F(kk);
        Q8val  =Q8F(kk);
        Qs2val  =Qs2F(kk);
        pklval =PSLF(kk);     

        ULB = -0.5* kk * pklval * j1; // * damping ;   //WARNING: REMOVE DAMPING!!!
        ULp = ULp + (ULA + ULB)*deltak/2.0;
        ULA = ULB;

        UloopB = -0.5* damping *  kk*( (5./21.)*R1val ) * j1;
        Uloopp = Uloopp + (UloopA + UloopB)*deltak/2.0;
        UloopA = UloopB;

        U11B = -0.5* damping * kk* 6./7. * ( R1val + R2val ) * j1;  
        U11p = U11p + (U11A + U11B)*deltak/2.0;
        U11A = U11B;

        U20B = -0.5* damping * kk*( (3./7.)*Q8val ) * j1;
        U20p = U20p + (U20A + U20B)*deltak/2.0;
        U20A = U20B;

        XLB = pklval*(1./3. - j1 / kq );
        XLp = XLp + (XLA + XLB)*deltak/2.0;
        XLA = XLB;

        XloopB = ( (9./98.)*Q1val + (10./21.)*R1val )*(1./3. - j1 / kq );
        Xloopp = Xloopp + (XloopA + XloopB)*deltak/2.0;
        XloopA = XloopB;

        X10B = (1./14.)*(
                    R1val * (2. + 3.* j0 ) - 2.*R2val
                    -3.*( 3.*R1val + 4.*R2val + 2.*Q5val ) * j1 / kq 
                    );
        X10p = X10p + (X10A + X10B)*deltak/2.0;
        X10A = X10B;

        YLB = pklval * j2;
        YLp = YLp + (YLA + YLB)*deltak/2.0;
        YLA = YLB;

        YloopB = ( (9./98.)*Q1val + (10./21.)*R1val ) * j2;
        Yloopp = Yloopp + (YloopA + YloopB)*deltak/2.0;
        YloopA = YloopB;

        Y10B = (-3./14.) * (3.* R1val + 4.* R2val + 2.* Q5val )
				*(j0 - 3. * j1 / kq );
        Y10p = Y10p + (Y10A + Y10B)*deltak/2.0;
        Y10A = Y10B;

        preVB = -(3./35.)*( Q1val - 3.*Q2val + 2.*R1val - 6.*R2val ) * j1 / kk;
        preVp = preVp + (preVA + preVB)*deltak/2.0;
        preVA = preVB;

        TB = -(9./14.)*(Q1val + 2.*Q2val + 2.*R1val + 4.*R2val )* j3 / kk;
        Tp = Tp + (TA + TB)*deltak/2.0;
        TA = TB;
        
        dotVaB = 3.0 * j1 * (-Q1val + 8.0 * Q2val - 7.0 * R1val + 16.0* R2val) / (70.  * kk) ;
        dotVap = dotVap + (dotVaA + dotVaB)*deltak/2.0;
        dotVaA = dotVaB;
        
        ddotVaB = 3.0 * j1 * (5.0 * Q2val - 5.0 * R1val + 10 * R2val ) / (35.0 * kk);
        ddotVap = ddotVap + (ddotVaA + ddotVaB)*deltak/2.0;
        ddotVaA = ddotVaB;
        
        dotVbB = - 3.0 * j1 * (Q1val - R1val ) / (14.0 * kk);
        dotVbp = dotVbp + (dotVbA + dotVbB)*deltak/2.0;
        dotVbA = dotVbB;
        
        V10B =  -1.0 / 7.0  * ki * Qs2val * j1 * damping;    //Mine
        //~ V10B = -2.0 / 7.0  * ki * Qs2val * j1 * damping;//;   //  (1609.02908)
        V10p = V10p + (V10A + V10B)*deltak/2.0;
        V10A = V10B;
        
        xiell0B = 0.5 *  k2* pklval *j0 * damping;  // This is also xiL
        xiell0p = xiell0p + (xiell0A + xiell0B)*deltak/2.0;
        xiell0A = xiell0B;
        
        xiell2B = 0.5 *  k2* pklval *j2 * damping; 
        xiell2p = xiell2p + (xiell2A + xiell2B)*deltak/2.0;
        xiell2A = xiell2B;
        
        xiell4B = 0.5 *  k2* pklval * j4 * damping;
        xiell4p = xiell4p + (xiell4A + xiell4B)*deltak/2.0;
        xiell4A = xiell4B;
       
        nabla2xiB = - k2 * xiell0B * damping * damping;
        nabla2xip = nabla2xip + (nabla2xiA + nabla2xiB)*deltak/2.0;
        nabla2xiA = nabla2xiB;
        


    }

    ULp /= PI2;    Uloopp /= PI2;
    U11p /= PI2; U20p /= PI2;
    XLp /= PI2;    Xloopp /= PI2;
    X10p /= 2.0*PI2; Y10p /= 2.0*PI2;
    YLp /= PI2; Yloopp /= PI2; 
    preVp = preVp/PI2; Tp = Tp/PI2;

    dotVap   /= PI2; ddotVap /= PI2;dotVbp   /= PI2;
    V10p /= PI2;
	xiell0p /= PI2; xiell2p /= PI2; xiell4p /= PI2;
	nabla2xip /= PI2;

    qqfun(qfunp) = qi;
    ULqfun(qfunp) = ULp;
    Uloopqfun(qfunp) = Uloopp;
    U11qfun(qfunp) = U11p;
    U20qfun(qfunp) = U20p;
    XLqfun(qfunp) = XLp;
    Xloopqfun(qfunp) = Xloopp;
    X10qfun(qfunp) = X10p;
    YLqfun(qfunp) = YLp;
    Yloopqfun(qfunp) = Yloopp;
    Y10qfun(qfunp) = Y10p;
    preVqfun(qfunp) = preVp;
    Tqfun(qfunp) = Tp;
    dotVaqfun(qfunp) =  dotVap;
    ddotVaqfun(qfunp) =  ddotVap;
    dotVbqfun(qfunp) =  dotVbp;
    V10qfun(qfunp) =  V10p;
    iBessel0qfun(qfunp) = xiell0p;
    iBessel2qfun(qfunp) = xiell2p;
    iBessel4qfun(qfunp) = xiell4p;
    nabla2xiqfun(qfunp) = nabla2xip;
                    
    return *qfunp;
}








local real PSLF(real k)
{
    real func;
    func = Interpolation_nr(k, kArrays.kT, kArrays.pklT, golists.Nk, kArraysd.pklT);
    return (func);
}

local real Q1F(real k)
{
    real func;
    func = Interpolation_nr(k, kArrays.kT, kArrays.Q1T, golists.Nk, kArraysd.Q1T);
    return (func);
}

local real Q2F(real k)
{
    real func;
    func = Interpolation_nr(k, kArrays.kT, kArrays.Q2T, golists.Nk, kArraysd.Q2T);
    return (func);
}

local real Q3F(real k)
{
    real func;
    func = Interpolation_nr(k, kArrays.kT, kArrays.Q3T, golists.Nk, kArraysd.Q3T);
    return (func);
}

local real Q5F(real k)
{
    real func;
    func = Interpolation_nr(k, kArrays.kT, kArrays.Q5T, golists.Nk, kArraysd.Q5T);
    return (func);
}

local real Q8F(real k)
{
    real func;
    func = Interpolation_nr(k, kArrays.kT, kArrays.Q8T, golists.Nk, kArraysd.Q8T);
    return (func);
}

local real Qs2F(real k)
{
    real func;
    func = Interpolation_nr(k, kArrays.kT, kArrays.Qs2T, golists.Nk, kArraysd.Qs2T);
    return (func);
}

local real R1F(real k)
{
    real func;
    func = Interpolation_nr(k, kArrays.kT, kArrays.R1T, golists.Nk, kArraysd.R1T);
    return (func);
}

local real R2F(real k)
{
    real func;
    func = Interpolation_nr(k, kArrays.kT, kArrays.R2T, golists.Nk, kArraysd.R2T);
    return (func);
}



local  real Interpolation_nr(real k, double kPS[], double pPS[], int nPS, double pPS2[])
{
    real psftmp;
    
//    if ( k < kPS[1] || k > kPS[nPS] )
//        fprintf(gd.outlog,"\n\nInterpolation_nr: warning! :: k is out of range... %g\n",k);
    
    splint(kPS,pPS,pPS2,nPS,k,&psftmp);
    
    return (psftmp);
}






/*
// correlation functions  
* It can be used for more precision in xiL
* Above we use a damping factor to avoid aliasing
global_corrfunctions correlation_functions(real qi)
{
    global_corrfunctions_ptr corrfunp;
    int i, Nk;
    real kmin, kmax, dk, kvali, kvalim1, ki, kim1;
    real kk;
    real deltak;
    
    real xip, xiA, xiB;
    real Lapxip, LapxiA, LapxiB;
    real nabla4xip, nabla4xiA, nabla4xiB;
//
    corrfunp = (global_corrfunctions_ptr) allocate(1 * sizeof(global_corrfunctions));
    
    kmin = 0.0001;
    kmax = 20.0;
//    Nk = 12000;
    Nk = 1200;
    if (Nk==1)
        dk = 0.;
    else
        dk = (rlog10(kmax) - rlog10(kmin))/((real)(Nk - 1));
    
    xip = 0.;
    xiA = rsqr(kmin)*PSLF(kmin)*rj0Bessel(kmin*qi);
    Lapxip = 0.;
    LapxiA = rpow(kmin,4.0)*PSLF(kmin)*rj0Bessel(kmin*qi)*rexp(-rsqr(kmin));
    nabla4xip = 0.;
    nabla4xiA = 0.; //aa-modificacion
//    nabla4xiA = rpow(kmin,6.0)*PSLF(kmin)*rj0Bessel(kmin*qi)*rexp(-rsqr(2.0*kmin));
//
    for (i=2; i<=Nk; i++) {
        kvali = rlog10(kmin) + dk*((real)(i - 1));
        kvalim1 = rlog10(kmin) + dk*((real)(i - 2));
        ki = rpow(10.0,kvali);
        kim1 = rpow(10.0,kvalim1);
        deltak = (ki - kim1);
        kk = rpow(10.0,kvali);
//
        xiB = rsqr(kk)*PSLF(kk)*rj0Bessel(kk*qi)*rexp(-rsqr(kk));
        xip = xip + (xiA + xiB)*deltak/2.0;
        xiA = xiB;
//
        LapxiB = rsqr(kk)*xiB*rexp(-rsqr(kk));
        Lapxip = Lapxip + (LapxiA + LapxiB)*deltak/2.0;
        LapxiA = LapxiB;
//
//        nabla4xiB = rpow(kk,4.0)*xiB*rexp(-rsqr(2*kk));
        nabla4xiB = 0; //aa-modificacion
        nabla4xip = nabla4xip + (nabla4xiA + nabla4xiB)*deltak/2.0;
        nabla4xiA = nabla4xiB;
    }
//
    xip /= TWOPI2;
    Lapxip /= -TWOPI2;
    nabla4xip /= TWOPI2;
//
    qcorrfun(corrfunp) = qi;
    xicorrfun(corrfunp) = xip;
    Lapxicorrfun(corrfunp) = Lapxip;
    nabla4xicorrfun(corrfunp) = nabla4xip;
    
    return *corrfunp;
}
*/



