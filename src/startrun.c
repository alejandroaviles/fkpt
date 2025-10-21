/*==============================================================================
 MODULE: startrun.c
 ==============================================================================*/

#include "globaldefs.h"
#include "protodefs.h"


local void ReadParameterFile(char *);
local void PrintParameterFile(char *);

local void startrun_parameterfile(void);
local void startrun_cmdline(void);
local void ReadParametersCmdline(void);
local void startrun_Common(void);
local void startrun_ParamStat(void);
local void CheckParameters(void);

local void InputPSTable(void);
local void PSLTable(void);
local void PSLTableNW(void);
local void GaussLegendrePoints(void);

void StartRun(string head0, string head1, string head2, string head3)
{
    real aTime;
    aTime = second();

    gd.headline0 = head0; gd.headline1 = head1;
    gd.headline2 = head2; gd.headline3 = head3;
    //~ printf("\n%s\n%s: %s\n\t %s\n",
		//~ gd.headline0, gd.headline1, gd.headline2, gd.headline3);

		//~ printf("\n \t\t fkpt - compute power spectrum in modified gravity models \n\n");
		//~ gd.headline0, gd.headline1, gd.headline2, gd.headline3);

    cmd.paramfile = GetParam("paramfile");
    if (!strnull(cmd.paramfile))
		startrun_parameterfile();
	else
		startrun_cmdline();

    StartOutput();

    fprintf(gd.outlog,"\nStartRun elapsed time: %g sec.\n\n",second()-aTime);
    fflush(gd.outlog);

    //~ stropen(gd.fpfnameTables,"w!");

}

local void startrun_parameterfile(void)
{
	ReadParameterFile(cmd.paramfile);
    	startrun_ParamStat();
	startrun_Common();
	PrintParameterFile(cmd.paramfile);
}

#define parameter_null	"parameters_null-mgpt"

local void startrun_cmdline(void)
{
	ReadParametersCmdline();
	startrun_Common();
	PrintParameterFile(parameter_null);
	//PrintParameterFile(cmd.paramfile);
}

local void ReadParametersCmdline(void)
{
stream outparams;


    cmd.chatty = GetiParam("chatty");

// bias and counterterms
    cmd.b1 = GetdParam("b1");
    cmd.b2 = GetdParam("b2");
    cmd.bs2 = GetdParam("bs2");
    cmd.b3nl = GetdParam("b3nl");
    cmd.alpha0 = GetdParam("alpha0");
    cmd.alpha2 = GetdParam("alpha2");
    cmd.alpha4 = GetdParam("alpha4");
    cmd.ctilde = GetdParam("ctilde");
    cmd.PshotP = GetdParam("pshotp");
    cmd.alpha0shot = GetdParam("alpha0shot");
    cmd.alpha2shot = GetdParam("alpha2shot");

    //~ cmd.c1eft = GetdParam("c1eft");
    //~ cmd.c2eft = GetdParam("c2eft");
    //~ cmd.s2eft = GetdParam("sigma2eft");
    cmd.c1eft = 0.0;
    cmd.c2eft = 0.0;
    cmd.s2eft = 0.0;
// output array params:
    //~ cmd.smin = GetdParam("smin");
    //~ cmd.smax = GetdParam("smax");
    //~ cmd.Ns = GetiParam("Ns");
    cmd.smin = 1;
    cmd.smax = 130;
    cmd.Ns = 100;

//  q functions:
    //~ cmd.NqperLogDecade = GetiParam("NqperLogDecade");
    //~ cmd.Nk_qFunctionsQuad = GetiParam("Nk_qFunctionsQuad");
    cmd.NqperLogDecade = 100;
    cmd.Nk_qFunctionsQuad = 1200;
//  For GSM integration:
    //~ cmd.gsm_width = GetdParam("gsm_width");  //Width of Gaussian in GSM integration
    //~ cmd.gsm_sizeyT = GetiParam("gsm_sizeyT"); //Number of points in 1-d GSM integration
    //~ cmd.gsm_NGL = GetiParam("gsm_NGL"); //Number of GL points in angular integration for rsd multipoles
    cmd.gsm_width = 100.;  //Width of Gaussian in GSM integration
    cmd.gsm_sizeyT = 160; //Number of points in 1-d GSM integration
    cmd.gsm_NGL = 16; //Number of GL points in angular integration for rsd multipoles
// Modified gravity model parameters:
    cmd.mgmodel = GetParam("model");   //";Modified gravity model to study (HS or DGP), default is LCDM",
    cmd.suffixModel = GetParam("suffixModel"); //";Suffix model to add to output filenames", ":suffix",
    cmd.model_paramfile = GetParam("modelParamfile");
    cmd.fR0 = GetdParam("fR0");
    cmd.screening =1.0;  ";set to =0 if you want no screenings", ":sc",
// DGP:
    cmd.eps_DGP = -1.0;
    cmd.rc_DGP = 1.0;
//
// Power spectrum table:
    cmd.fnamePS = GetParam("fnamePS");
    cmd.is_PS_input_LCDM = GetiParam("is_PS_input_LCDM");
    cmd.kmin = GetdParam("kmin");
    cmd.kmax = GetdParam("kmax");
    cmd.Nk = GetiParam("Nk");
//
// CLPT correlation functions table:
    cmd.rmin = 1.;
    cmd.rmax = 210.;
    cmd.Nr = 210;
// Background cosmology:
    cmd.om = GetdParam("Om");
    cmd.olstr ="1 - Om";  // This is for DGP only
    cmd.h = GetdParam("h");
//
// Differential equations evolution parameters:
    cmd.x =-4.0;    //";Initial conformal time value :: Log[1/(1 + zini)]",
    cmd.dxstr = "2/5";  //";Conformal time integration step",
    cmd.dxmin = 0.;  // ";Min conformal time integration step size",
    cmd.eps = 0.0001;  //";Differential equations solver tolerance error parameter",":epssolver",
    cmd.xstop = GetdParam("zout");
    cmd.maxnsteps = 10000;  //";Maximum number of integration steps", ":maxn",
    cmd.integration_method = "rkqs";  //";Integration method to use", ":solver",
//
// Quadrature parameters:
    cmd.quadratureMethod = "trapezoid3"; //  Quadrature method to use", ":quadm
    cmd.nquadSteps = GetiParam("nquadSteps"); // Number of k´s from the power spectrum table to integrate (trapezoid)",":nquad",
    cmd.ngausslegpoints = 16;  //NOT USED   ;Number of Gauss-Legendre of integration points", ":nglpts
    cmd.epsquad = 1.0e-6;    //NOT USED  Quadrature tolerance error parameter (open Romberg method: romberg)
// Post processing parameters:
    cmd.postprocessing = FALSE;
    cmd.options = "";



    if(cmd.chatty==1){
		printf("\n \t\t fkpt - compute power spectrum in modified gravity models \n\n");

		fprintf(stdout,"Reading pkl from file %s in the form column1: k[h/Mpc], 2: pkl[Mpc/h]^3 \n", cmd.fnamePS);
		fprintf(stdout,"redshift: z=%g\n",    cmd.xstop);
		fprintf(stdout,"OmegaM=%g\n",         cmd.om);
		fprintf(stdout,"h=%g\n",              cmd.h);
		fprintf(stdout,"b1=%g\n",             cmd.b1);
		fprintf(stdout,"b2=%g\n",             cmd.b2);
		fprintf(stdout,"bs2=%g\n",            cmd.bs2);
		fprintf(stdout,"b3nl=%g\n",           cmd.b3nl);
		fprintf(stdout,"alpha0=%g\n",         cmd.alpha0);
		fprintf(stdout,"alpha2=%g\n",         cmd.alpha2);
		fprintf(stdout,"alpha4=%g\n",         cmd.alpha4);
		fprintf(stdout,"ctilde=%g\n",         cmd.ctilde);
		fprintf(stdout,"Pshotp=%g\n",         cmd.PshotP);
		fprintf(stdout,"alpha0shot=%g\n",     cmd.alpha0shot);
		fprintf(stdout,"alpha2shot=%g\n",     cmd.alpha2shot);
	}

}

#undef parameter_null

local void startrun_Common(void)
{
	real dx1, dx2;
    char *ep;

    setFilesDirs_log();
    strcpy(gd.mode,"w");
	if(!(gd.outlog=fopen(gd.logfilePath, gd.mode)))
		error("\nstart_Common: error opening file '%s' \n",gd.logfilePath);

//
// Background cosmology:
    ep = strchr(cmd.olstr, '-');
    if (ep == NULL) {
        gd.ol = GetdParam("OL");
        fprintf(gd.outlog,"\nOL string input without '-' :: %s\n",cmd.olstr);
        fprintf(gd.outlog,"\nOLambda, Om and sum  : %g %g %g\n",gd.ol, cmd.om, gd.ol+cmd.om);
    } else {
        ep = strchr(cmd.olstr, '1');
        if (ep == NULL)
            error("\nstart_Common: OL not in the format '1 - Om' \n");
        else {
            ep = strchr(cmd.olstr, 'O');
            if (ep == NULL)
                error("\nstart_Common: OL not in the format '1 - Om' \n");
            else
                if (*(ep+1) == 'm') {
                    fprintf(gd.outlog,"\nFound Om\n");
                    gd.ol = 1. - cmd.om;
                    fprintf(gd.outlog,"\nOLambda and Om : %g %g\n",gd.ol, cmd.om);
                } else {
                    fprintf(gd.outlog,"\nNot found Om\n");
                    error("\nstart_Common: OL not in the format (1 - Om) \n");
                }
        }

    }


    gd.dx = (sscanf(cmd.dxstr, "%lf/%lf", &dx1, &dx2) == 2 ?
				dx1/dx2 : atof(cmd.dxstr));
    if ( dx2 == 0. )
        error("\n\nstartrun_Common: dx : dx2 must be finite\n");

    CheckParameters();
    GaussLegendrePoints();
    quadraturemethod_string_to_int(cmd.quadratureMethod, &gd.quadmethod_int);
	integration_method_string_to_int(cmd.integration_method, &gd.method_int);
    gd.xnow = cmd.x;
    gd.xout = gd.xnow;
    gd.xoutinfo = gd.xnow;

    gd.xstop = rlog(1.0/(1.0+cmd.xstop));

    set_model();
    setFilesDirs();

    if (!strnull(cmd.fnamePS)) {
        sprintf(gd.fnamePSPath,"Input/%s",cmd.fnamePS);
        fprintf(gd.outlog,"PS file Path and file name: %s\n",gd.fnamePSPath);
        fflush(gd.outlog);
        InputPSTable();
        PSLTable();
        PSLTableNW();
    }
    if (cmd.nquadSteps > nPSLT)
        error("CheckParameters: nquadSteps > nPSLT\n");
    /*
#define usrmodel_parameter_null    "usrmodel_parameters_null"

    if ( (strcmp(cmd.mgmodel,"USER") == 0) || (strcmp(cmd.mgmodel,"user") == 0))
        if (!strnull(cmd.model_paramfile)) {
            fprintf(stdout,"\n\nUser model :: using parameter file: %s\n",cmd.model_paramfile);
            ReadMGModelParameterFile(cmd.model_paramfile);
            PrintMGModelParameterFile(cmd.model_paramfile);
        } else
            PrintMGModelParameterFile(usrmodel_parameter_null);
#undef usrmodel_parameter_null
*/
}

local void startrun_ParamStat(void)
{
	real dx1, dx2;
// output array params:
    if (GetParamStat("smin") & ARGPARAM)
        cmd.smin = GetdParam("smin");
    if (GetParamStat("smax") & ARGPARAM)
        cmd.smax = GetdParam("smax");
    if (GetParamStat("Nk") & ARGPARAM)
        cmd.smin = GetiParam("Nk");
// bias and counterterms
    if (GetParamStat("b1") & ARGPARAM)
        cmd.b1 = GetdParam("b1");
    if (GetParamStat("b2") & ARGPARAM)
        cmd.b2 = GetdParam("b2");
    if (GetParamStat("bs") & ARGPARAM)
        cmd.bs2 = GetdParam("bs2");
//    cmd.bs = GetdParam("bs");
    if (GetParamStat("c1eft") & ARGPARAM)
        cmd.c1eft = GetdParam("c1eft");
    if (GetParamStat("c2eft") & ARGPARAM)
        cmd.c2eft = GetdParam("c2eft");
    if (GetParamStat("sigma2eft") & ARGPARAM)
        cmd.s2eft = GetdParam("sigma2eft");
//  q functions:
    if (GetParamStat("NqperLogDecade") & ARGPARAM)
        cmd.NqperLogDecade = GetiParam("NqperLogDecade");
    if (GetParamStat("Nk_qFunctionsQuad") & ARGPARAM)
        cmd.Nk_qFunctionsQuad = GetiParam("Nk_qFunctionsQuad");
//  For GSM integration:
    if (GetParamStat("gsm_width") & ARGPARAM)
        cmd.gsm_width = GetdParam("gsm_width");
    if (GetParamStat("gsm_sizeyT") & ARGPARAM)
        cmd.gsm_sizeyT = GetiParam("gsm_sizeyT");
    if (GetParamStat("gsm_NGL") & ARGPARAM)
        cmd.gsm_NGL = GetiParam("gsm_NGL");
// Power spectrum table:
    if (GetParamStat("fnamePS") & ARGPARAM)
        cmd.fnamePS = GetParam("fnamePS");
    if (GetParamStat("is_PS_input_LCDM") & ARGPARAM)
        cmd.is_PS_input_LCDM = GetiParam("is_PS_input_LCDM");
    if (GetParamStat("kmin") & ARGPARAM)
        cmd.kmin = GetdParam("kmin");
    if (GetParamStat("kmax") & ARGPARAM)
        cmd.kmax = GetdParam("kmax");
    if (GetParamStat("Nk") & ARGPARAM)
        cmd.Nk = GetiParam("Nk");
// CLPT correlation functions table:
    cmd.rmin = 1.;
    cmd.rmax = 210.;
    cmd.Nr = 210;
// Modified gravity model parameters:
    cmd.mgmodel = "LCDM";
    if (GetParamStat("suffixModel") & ARGPARAM)
    cmd.suffixModel = GetParam("suffixModel");
    cmd.fR0 = 1.0e-10;
    cmd.screening =1.0;
// DGP:
    cmd.eps_DGP = -1.0;
    cmd.rc_DGP = 1.0;
// Background cosmology:
    if (GetParamStat("Om") & ARGPARAM)
        cmd.om = GetdParam("Om");
    cmd.olstr ="1 - Om";  // This is for DGP only
    if (GetParamStat("h") & ARGPARAM)
        cmd.h = GetdParam("h");
// Differential equations evolution parameters:
    cmd.x =-4.0;
    cmd.dxstr = "2/5";
	gd.dx = (sscanf(cmd.dxstr, "%lf/%lf", &dx1, &dx2) == 2 ?
                    dx1/dx2 : atof(cmd.dxstr));
	if ( dx2 == 0. )
		error("\n\nstartrun_ParamStat: deta : deta2 must be finite\n");

	if (GetParamStat("zout") & ARGPARAM)
		cmd.xstop = GetdParam("zout");
    cmd.eps = 0.0001;
    cmd.maxnsteps = 10000;

    cmd.integration_method = "rkqs";
	//fprintf(gd.outlog,"\n\nrunning now %s integration method ...\n",
		//		cmd.integration_method);
// Quadrature parameters:
    cmd.quadratureMethod = "trapezoid3";
        //fprintf(gd.outlog,"\n\nrunning instead %s quadrature method ...\n",
             //   cmd.quadratureMethod);
    if (GetParamStat("nquadSteps") & ARGPARAM)
        cmd.nquadSteps = GetiParam("nquadSteps");
    cmd.ngausslegpoints = 16;  //NOT USED   ;Number of Gauss-Legendre of integration points", ":nglpts
    cmd.epsquad = 1.0e-6;

// Post processing parameters:
    cmd.postprocessing = FALSE;
    cmd.options = "";
}

local void CheckParameters(void)
{
// Power spectrum table:
    if (strnull(cmd.fnamePS))
        error("CheckParameters: You should give a power spectrum filename\n");
    if (cmd.kmin < 0.0)
        error("CheckParameters: absurd value for kmin\n");
    if (cmd.kmax < 0.0)
        error("CheckParameters: absurd value for kmax\n");
    if (cmd.kmin > cmd.kmax)
        error("CheckParameters: kmin can not be greater than kmax\n");
    if (cmd.Nk < 0)
        error("CheckParameters: absurd value for Nk\n");
//
// CLPT correlation functions table:
    if (cmd.rmin < 0.0)
        error("CheckParameters: absurd value for rmin\n");
    if (cmd.rmax < 0.0)
        error("CheckParameters: absurd value for rmax\n");
    if (cmd.rmin > cmd.rmax)
        error("CheckParameters: rmin can not be greater than rmax\n");
    if (cmd.Nr < 0)
        error("CheckParameters: absurd value for Nr\n");
// Background cosmology:
    if (cmd.om > 1.0 || cmd.om < 0.0)
        error("CheckParameters: absurd value for om\n");
    if ( gd.ol < 0. )
        error("\n\nstartrun_ParamStat: OL (=%g) : must be positive\n",gd.ol);
    if (cmd.h < 0.0)
        error("CheckParameters: absurd value for h\n");

// Differential equations evolution parameters:
    if (gd.dx == 0)
        error("CheckParameters: absurd value for deta\n");
    if(cmd.x == rlog(1.0/(1.0+cmd.xstop)) )
        error("\n\nstartrun_Common: etaini and etaout=exp(-zout)-1 must be different\n");

    if (cmd.eps > 1.0e-4 || cmd.eps <= 0)
        error("CheckParameters: inapropriate or absurd value for epsSolver\n");
    if (cmd.maxnsteps < 1)
        error("CheckParameters: absurd value for maxnsteps\n");

// Quadrature parameters:
    if (cmd.nquadSteps <= 1)
        error("CheckParameters: absurd value for nquadSteps\n");
    if (cmd.ngausslegpoints <= 1)
        error("CheckParameters: absurd value for ngausslegpoints\n");
    if (cmd.epsquad >= 1.0e-1 || cmd.epsquad <= 0)
        error("CheckParameters: absurd value for epsquad\n");
}

local void ReadParameterFile(char *fname)
{
#define DOUBLE 1
#define STRING 2
#define INT 3
#define BOOLEAN 4
#define MAXTAGS 300

  FILE *fd,*fdout;

  char buf[200],buf1[200],buf2[200],buf3[200];
  int  i,j,nt;
  int  id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int  errorFlag=0;

  nt=0;
    SPName(cmd.fnamePS,"fnamePS",100);
    IPName(cmd.is_PS_input_LCDM,"is_PS_input_LCDM");
    RPName(cmd.xstop,"zout");
    RPName(cmd.om,"Om");
    RPName(cmd.h,"h");
// bias and counterterms
    RPName(cmd.b1,"b1");
    RPName(cmd.b2,"b2");
    RPName(cmd.bs2,"bs2");
//    RPName(cmd.bs,"bs");
    RPName(cmd.c1eft,"c1eft");
    RPName(cmd.c2eft,"c2eft");
    RPName(cmd.s2eft,"sigma2eft");
// output array params:
    RPName(cmd.smin,"smin");
    RPName(cmd.smax,"smax");
    IPName(cmd.Ns,"Ns");
//  q functions:
    IPName(cmd.NqperLogDecade,"NqperLogDecade");
    IPName(cmd.Nk_qFunctionsQuad,"Nk_qFunctionsQuad");
//  For GSM integration:
    RPName(cmd.gsm_width,"gsm_width");  //Width of Gaussian in GSM integration
    IPName(cmd.gsm_sizeyT,"gsm_sizeyT"); //Number of points in 1-d GSM integration
    IPName(cmd.gsm_NGL,"gsm_NGL"); //Number of GL points in angular integration for rsd multipoles
// Power spectrum table:
    RPName(cmd.kmin,"kmin");
    RPName(cmd.kmax,"kmax");
    IPName(cmd.Nk,"Nk");
//
// CLPT correlation functions table:
    cmd.rmin = 1.;
    cmd.rmax = 210.;
    cmd.Nr = 210;
// Modified gravity model parameters:
    SPName(cmd.suffixModel,"suffixModel",100);
    SPName(cmd.model_paramfile,"modelParamfile",100);
    cmd.mgmodel = "LCDM";
    cmd.fR0 = 1.0e-10;
    cmd.screening =1.0;
// DGP:
    cmd.eps_DGP = -1.0;
    cmd.rc_DGP = 1.0;
    cmd.olstr ="1 - Om";
//
// Differential equations evolution parameters:

    cmd.dxstr = "2/5";
    cmd.eps = 0.0001;
    cmd.maxnsteps = 10000;
// Quadrature parameters:

    IPName(cmd.nquadSteps,"nquadSteps");
    cmd.integration_method = "rkqs";
    cmd.quadratureMethod = "trapezoid3";
    cmd.ngausslegpoints = 16;
    cmd.epsquad = 1.0e-6;
// Post processing parameters:
    cmd.postprocessing = FALSE;
    SPName(cmd.options,"options",100);
//

	if((fd=fopen(fname,"r"))) {
		while(!feof(fd)) {
			fgets(buf,200,fd);
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<1)
                continue;
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<2)
                *buf2='\0';
            if(buf1[0]=='%')
				continue;
            for(i=0,j=-1;i<nt;i++)
                if(strcmp(buf1,tag[i])==0) {
                    j=i;
					tag[i][0]=0;
					break;
				}
			if(j>=0) {
                switch(id[j]) {
					case DOUBLE:
						*((double*)addr[j])=atof(buf2);
						break;
					case STRING:
						strcpy(addr[j],buf2);
						break;
					case INT:
						*((int*)addr[j])=atoi(buf2);
						break;
					case BOOLEAN:
						if (strchr("tTyY1", *buf2) != NULL) {
							*((bool*)addr[j])=TRUE;
                        } else
                            if (strchr("fFnN0", *buf2) != NULL)  {
                                *((bool*)addr[j])=FALSE;
                            } else {
                                error("getbparam: %s=%s not bool\n",buf1,buf2);
                            }
						break;
                }
            } else {
                fprintf(stdout, "Error in file %s: Tag '%s' %s.\n",
					fname, buf1, "not allowed or multiple defined");
                errorFlag=1;
            }
        }
        fclose(fd);
    } else {
        fprintf(stdout,"Parameter file %s not found.\n", fname);
        errorFlag=1;
        exit(1);
    }

    for(i=0;i<nt;i++) {
        if(*tag[i]) {
            fprintf(stdout,
                "Error. I miss a value for tag '%s' in parameter file '%s'.\n",
                tag[i],fname);
            exit(0);
        }
    }
#undef DOUBLE
#undef STRING
#undef INT
#undef BOOLEAN
#undef MAXTAGS
}

#define FMTT	"%-35s%s\n"
#define FMTI	"%-35s%d\n"
#define FMTR	"%-35s%g\n"

local void PrintParameterFile(char *fname)
{
    FILE *fdout;
    char buf[200];

    sprintf(buf,"%s/%s%s%s",gd.tmpDir,fname,cmd.suffixModel,"-usedvalues");
    if(!(fdout=fopen(buf,"w"))) {
        fprintf(stdout,"error opening file '%s' \n",buf);
        exit(0);
    } else {
        fprintf(fdout,"%s\n",
		"%-------------------------------------------------------------------");
        fprintf(fdout,"%s %s\n","% Parameter input file for:",gd.headline0);
        fprintf(fdout,"%s\n","%");
        fprintf(fdout,"%s %s: %s\n%s\t    %s\n","%",gd.headline1,gd.headline2,"%",
            gd.headline3);
        fprintf(fdout,"%s\n%s\n",
		"%-------------------------------------------------------------------",
		"%");



        fprintf(fdout,FMTT,"fnamePS",cmd.fnamePS);
        fprintf(fdout,FMTR,"zout",cmd.xstop);
        fprintf(fdout,FMTR,"Om",cmd.om);
        fprintf(fdout,FMTR,"h",cmd.h);
        fprintf(fdout,FMTR,"b1",cmd.b1);
        fprintf(fdout,FMTR,"b2",cmd.b2);
//        fprintf(fdout,FMTR,"bs",cmd.bs);
        fprintf(fdout,FMTR,"bs2",cmd.bs2);
        fprintf(fdout,FMTR,"c1eft",cmd.c1eft);
        fprintf(fdout,FMTR,"c2eft",cmd.c2eft);
        fprintf(fdout,FMTR,"sigma2eft",cmd.s2eft);


	fprintf(fdout,FMTR,"smin",cmd.smin);
	fprintf(fdout,FMTR,"smax",cmd.smax);
	fprintf(fdout,FMTI,"Ns",cmd.Ns);

	fprintf(fdout,FMTR,"kmin",cmd.kmin);
	fprintf(fdout,FMTR,"kmax",cmd.kmax);
	fprintf(fdout,FMTI,"Nk",cmd.Nk);
	fprintf(fdout,FMTT,"quadratureMethod",cmd.quadratureMethod);

	fprintf(fdout,FMTI,"NqperLogDecade",cmd.NqperLogDecade);
	fprintf(fdout,FMTI,"Nk_qFunctionsQuad",cmd.Nk_qFunctionsQuad);

	fprintf(fdout,FMTR,"gsm_width",cmd.gsm_width);
	fprintf(fdout,FMTI,"gsm_sizeyT",cmd.gsm_sizeyT);
	fprintf(fdout,FMTI,"gsm_NGL",cmd.gsm_NGL);
	fprintf(fdout,FMTT,"suffixModel",cmd.suffixModel);
	fprintf(fdout,FMTT,"modelParamfile",cmd.model_paramfile);
// Power spectrum table:


//
// CLPT correlation functions table:
       // fprintf(fdout,FMTR,"rmin",cmd.rmin);
       // fprintf(fdout,FMTR,"rmax",cmd.rmax);
       // fprintf(fdout,FMTI,"Nr",cmd.Nr);
// Modified gravity model parameters:
      //  fprintf(fdout,FMTT,"mgModel",cmd.mgmodel);
//        fprintf(fdout,FMTI,"nHS",cmd.nHS);
      //  fprintf(fdout,FMTR,"fR0",cmd.fR0);
      //  fprintf(fdout,FMTR,"screening",cmd.screening);
// DGP:
      //  fprintf(fdout,FMTR,"epsDGP",cmd.eps_DGP);
      //  fprintf(fdout,FMTR,"rcDGP",cmd.rc_DGP);
      //  fprintf(fdout,FMTT,"OL",cmd.olstr);
//
// Differential equations evolution parameters:
     //   fprintf(fdout,FMTR,"etaini",cmd.x);
     //   fprintf(fdout,FMTR,"detamin",cmd.dxmin);
     //   fprintf(fdout,FMTR,"epsSolver",cmd.eps);
    //    fprintf(fdout,FMTT,"deta",cmd.dxstr);
     //   fprintf(fdout,FMTI,"maxnsteps",cmd.maxnsteps);
      //  fprintf(fdout,FMTT,"solverMethod",cmd.integration_method);
//
// Quadrature parameters:
    //    fprintf(fdout,FMTI,"nquadSteps",cmd.nquadSteps);
    //    fprintf(fdout,FMTI,"ngausslegpoints",cmd.ngausslegpoints);
    //    fprintf(fdout,FMTR,"epsquad",cmd.epsquad);
// Post processing parameters:
    //    fprintf(fdout,FMTT,"postprocessing",cmd.postprocessing ? "true" : "false");
    //    fprintf(fdout,FMTT,"options",cmd.options);
//
        fprintf(fdout,"\n\n");
    }
    fclose(fdout);
}

#undef FMTT
#undef FMTI
#undef FMTR


#define NPT 10
#define SPREAD 1.0
local void InputPSTable(void)
{
    stream outstr;
    pointPSTableptr p, plog, pn;
    pointPSTableptr PSLCDMtabtmp;
    int nPSTabletmp;
    int i;
    real dk, kval, PSval, kmin, kmax;
    int mwt;
    double al,bl,chi2,q,siga,sigb,*x,*y,*sig;
    double au, bu;
    real *kPStmp;
    real *pPStmp;
    real *pPS2tmp;
    char namebuf[256];
    real kminext, kmaxext, dktmp, kmn, kmx;
    int Nkext=800, NkL=50, NkU=50;
    real kminT=1.0e-5, kmaxT=400.0;

    fprintf(gd.outlog,"\n\nReading power spectrum from file %s...\n",gd.fnamePSPath);
    inout_InputData(gd.fnamePSPath, 1, 2, &nPSTabletmp);

    if (nPSTabletmp < 1)
        error("\n\nInputPSTable: nPSTable = %d is absurd\n\n", nPSTabletmp);

    PSLCDMtabtmp = (pointPSTableptr) allocate(nPSTabletmp * sizeof(pointPSTable));

    fprintf(gd.outlog,"nPSTable : %d\n", nPSTabletmp);

    i = 0;
    for (p=PSLCDMtabtmp; p<PSLCDMtabtmp+nPSTabletmp; p++) {
        kPos(p) = inout_xval[i];
        PS(p) = inout_yval[i];
        ++i;
    }

    fprintf(gd.outlog,"\n\nCreating log power spectrum...\n");

    PSLCDMLogtab = (pointPSTableptr) allocate(nPSTabletmp * sizeof(pointPSTable));
    plog = PSLCDMLogtab;
    nPSLogT=0;
    for (p=PSLCDMtabtmp; p<PSLCDMtabtmp+nPSTabletmp; p++) {
        kPos(plog) = rlog10(kPos(p));
        PS(plog) = rlog10(PS(p));
        plog++;
        nPSLogT++;
    }

    fprintf(gd.outlog,"\nTotal numbers in Log PS: %d %ld\n",nPSLogT,plog-PSLCDMLogtab);
    fprintf(gd.outlog,"Total numbers in Normal PS: %d %ld\n\n",nPSTabletmp,p-PSLCDMtabtmp);

    fprintf(gd.outlog,"\n\nLinear fit (a + b x) to log-log power spectrum at minset and maxset...\n");
    fprintf(gd.outlog,"\n\nLinear fit to log-log power spectrum at minset and maxset...\n");

    x=dvector(1,NPT);
    y=dvector(1,NPT);
    sig=dvector(1,NPT);

// Low-ks of the PS
    fprintf(gd.outlog,"\nAt low-ks of the spectrum...\n");

    plog = PSLCDMLogtab;
    for (i=1;i<=NPT;i++) {
        x[i]=kPos(plog);
        y[i]=PS(plog);
        sig[i]=SPREAD;
        plog++;
    }
    for (mwt=0;mwt<=1;mwt++) {
        fit(x,y,NPT,sig,mwt,&al,&bl,&siga,&sigb,&chi2,&q);
        if (mwt == 0)
            fprintf(gd.outlog,"\nIgnoring standard deviations\n");
        else
            fprintf(gd.outlog,"\nIncluding standard deviations\n");
        fprintf(gd.outlog,"%12s %9.6f %18s %9.6f \n",
               "a  =  ",al,"uncertainty:",siga);
        fprintf(gd.outlog,"%12s %9.6f %18s %9.6f \n",
               "b  =  ",bl,"uncertainty:",sigb);
        fprintf(gd.outlog,"%19s %14.6f \n","chi-squared: ",chi2);
        fprintf(gd.outlog,"%23s %10.6f \n","goodness-of-fit: ",q);
    }

// High-ks of the PS
    fprintf(gd.outlog,"\nAt high-ks of the spectrum...\n");

    plog = PSLCDMLogtab+nPSLogT-1;
    for (i=1;i<=NPT;i++) {
        x[i]=kPos(plog);
        y[i]=PS(plog);
        sig[i]=SPREAD;
        plog--;
    }
    for (mwt=0;mwt<=1;mwt++) {
        fit(x,y,NPT,sig,mwt,&au,&bu,&siga,&sigb,&chi2,&q);
        if (mwt == 0)
            fprintf(gd.outlog,"\nIgnoring standard deviations\n");
        else
            fprintf(gd.outlog,"\nIncluding standard deviations\n");
        fprintf(gd.outlog,"%12s %9.6f %18s %9.6f \n",
               "a  =  ",au,"uncertainty:",siga);
        fprintf(gd.outlog,"%12s %9.6f %18s %9.6f \n",
               "b  =  ",bu,"uncertainty:",sigb);
        fprintf(gd.outlog,"%19s %14.6f \n","chi-squared: ",chi2);
        fprintf(gd.outlog,"%23s %10.6f \n","goodness-of-fit: ",q);
    }

// Extending power spectrum
    kPStmp = dvector(1,nPSLogT);
    pPStmp = dvector(1,nPSLogT);
    pPS2 = dvector(1,nPSLogT);
    pn = PSLCDMLogtab;
    i=1;
    for (pn = PSLCDMLogtab; pn<PSLCDMLogtab+nPSLogT; pn++) {
        kPStmp[i] = kPos(pn);
        pPStmp[i] = PS(pn);
        i++;
    }
//
    spline(kPStmp,pPStmp,nPSLogT,1.0e30,1.0e30,pPS2);

    kmin = kPos(PSLCDMtabtmp);
    kmax = kPos(PSLCDMtabtmp+nPSTabletmp-1);
    fprintf(gd.outlog,"\nkmin, kmax of the given power spectrum (with %f values): %g %d",kmin, kmax, nPSTabletmp);
    dktmp = (rlog10(kmax) - rlog10(kmin))/((real)(nPSTabletmp - 1));
    kminext = rpow(10.0, rlog10(kmin)-((real)NkL)*dktmp);
    kmaxext = rpow(10.0, rlog10(kmax)+((real)NkU)*dktmp);
//
    fprintf(gd.outlog,"\n\nNkL, NkU: %d %d\n",NkL, NkU);
    NkL = rlog10(kmin/kminT)/dktmp;
    NkU = rlog10(kmaxT/kmax)/dktmp;
    fprintf(gd.outlog,"\n\nNkL, NkU targets: %d %d\n",NkL,NkU);
    kminext = rpow(10.0, rlog10(kmin)-((real)NkL)*dktmp);
    kmaxext = rpow(10.0, rlog10(kmax)+((real)NkU)*dktmp);
//

    fprintf(gd.outlog,"\nkmin, kmax of the extended power spectrum (first try): %g %g",
            kminext, kmaxext);

    kmn = MIN(kminext,cmd.kmin);
    kmx = MAX(kmaxext,cmd.kmax);
    fprintf(gd.outlog,"\nkmin, kmax of the extended power spectrum (second try): %g %g\n",
            kmn, kmx);

    nPSTable = Nkext;
    fprintf(gd.outlog,"\n\nCreating new PSTable with %d values\n",nPSTable);
    PSLCDMtab = (pointPSTableptr) allocate(nPSTable * sizeof(pointPSTable));
    dk = (rlog10(kmx) - rlog10(kmn))/((real)(nPSTable - 1));
    p = PSLCDMtab;

    for (i=1; i<=nPSTable; i++) {
        kval = rlog10(kmn) + dk*((real)(i - 1));
        if (rpow(10.0,kval) >= kmin && rpow(10.0,kval) <= kmax)
            PSval = psInterpolation_nr(kval, kPStmp, pPStmp, nPSLogT);
        else
            if (rpow(10.0,kval) < kmin)
                PSval = al + bl*kval;
            else
                if (rpow(10.0,kval) > kmax)
                    PSval = au + bu*kval;
                else
                    error("\n\nError: InputPSTable :: kmin, kmax, kval: %g %g %g", kmin, kmax, kval);
//
        kPos(p) = rpow(10.0,kval);
        PS(p) = rpow(10.0,PSval);
        p++;
    }

    sprintf(namebuf,"%s/%s%s_%s",gd.tmpDir,cmd.fnamePS,cmd.suffixModel,"ext.dat");
    outstr = stropen(namebuf,"w!");
    for (p=PSLCDMtab; p<PSLCDMtab+nPSTable; p++) {
        fprintf(outstr,"%g %g\n",
                kPos(p),PS(p));
    }
    fclose(outstr);
//

    free_dvector(pPS2,1,nPSLogT);
    free_dvector(pPStmp,1,nPSLogT);
    free_dvector(kPStmp,1,nPSLogT);

    free_dvector(sig,1,NPT);
    free_dvector(y,1,NPT);
    free_dvector(x,1,NPT);
}
#undef NPT
#undef SPREAD



local void PSLTable(void)
{
    char namebuf[256];
    stream outstr;
    real kmin, Dpkmin, Dpk, f0h;
    pointPSTableptr p, pn;
    int i;
    int rescalePS = 0;

//
    real xstoptmp, Dp0, Dpzout, fac;

    xstoptmp = gd.xstop;
    gd.xstop = 0.;

    Dp0 = DpFunction_LCDM(0.); // LCDM
    fprintf(gd.outlog,"\n\n Dp(0) = %g",Dp0);
    gd.xstop = xstoptmp;
    Dpzout = DpFunction(0.); // aa: This was for rescale the input pk redshift, but I removed it
    fprintf(gd.outlog,"\n Dp(%g) = %g\n",cmd.xstop,Dpzout);

    //~ fprintf(stdout,"\nzout  = %g", cmd.xstop);
    gd.Dplus=Dpzout/Dp0;
    //~ fprintf(stdout,"\nDplus= %g\n\n", gd.Dplus);//Dpzout/Dp0);


    if (rescalePS !=0){  // aa: This was for rescale the input pk redshift, but I removed it
		fac = rsqr(Dpzout/Dp0);
        fprintf(stdout,"Warning: rescaling the redshift with rsqr(Dpzout/Dp0) = %g", fac);
		for (p = PSLCDMtab; p<PSLCDMtab+nPSTable; p++) {
			PS(p) *= fac;
		}
	}

    sprintf(namebuf,"%s/%s%s_%s",gd.tmpDir,cmd.fnamePS,cmd.suffixModel,"ext2.dat");
    outstr = stropen(namebuf,"w!");
    for (p=PSLCDMtab; p<PSLCDMtab+nPSTable; p++) {
        fprintf(outstr,"%g %g\n",
                kPos(p),PS(p));
    }
    fclose(outstr);

    kmin = kPos(PSLCDMtab);
    Dpkmin = DpFunction(kmin);
    fprintf(gd.outlog,"\n\n Dpkmin = %g\n",Dpkmin);

    PSLT = (pointPSTableptr) allocate(nPSTable * sizeof(pointPSTable));
    nPSLT = 0;
    pn = PSLT;
    if (cmd.is_PS_input_LCDM==1) {
        printf("\nRescaling input power spectrum from LCDM to MG. Change to is_PS_input_LCDM=0 if dont want to\n");
        for (p = PSLCDMtab; p<PSLCDMtab+nPSTable; p++) {
            kPos(pn) = kPos(p);
            Dpk = DpFunction(kPos(p));
            kPos(pn) = kPos(p);
            Dpk = DpFunction(kPos(p));
            // printf("%g ", rsqr(Dpk/Dpkmin));
            PS(pn) = rsqr(Dpk/Dpkmin)*PS(p);
            pn++;
            nPSLT++;
        }
    }

    kPS = dvector(1,nPSLT);
    pPS = dvector(1,nPSLT);
    pPS2 = dvector(1,nPSLT);
    pn = PSLT;
    i=1;
    for (pn = PSLT; pn<PSLT+nPSLT; pn++) {
        kPS[i] = kPos(pn);
        pPS[i] = PS(pn);
        i++;
    }
//

    fkT=dvector(1,nPSLT);
    fkT2 = dvector(1,nPSLT);
    DplusT=dvector(1,nPSLT);
    DplusT2 = dvector(1,nPSLT);

    for (int i = 1; i<=nPSLT; i++) {
        fkT[i] = f_growth(kPS[i]);
        DplusT[i] = DpFunction(kPS[i]);
    }

    gd.f0=fkT[1];


    spline(kPS,pPS,nPSLT,1.0e30,1.0e30,pPS2);
    spline(kPS,fkT,nPSLT,1.0e30,1.0e30,fkT2);

    sprintf(namebuf,"%s/%s%s%s",gd.clptDir,"linear",cmd.suffixModel,".dat");
    outstr = stropen(namebuf,"w!");

    fprintf(outstr,"# InputPklFile=%s\n", cmd.fnamePS);
    fprintf(outstr,"# z=%g, OmegaM=%g, h=%g\n",
		cmd.xstop, cmd.om, cmd.h);
    fprintf(outstr,"%1s%11s%17s%20s%22s",
            "#"," 1.k[h/Mpc]","2.pkl","3.f(k)","4.D+(k)\n");

    for (i=1; i<=nPSLT; i++) {
        fprintf(outstr,"%-19.9e %-19.9e %-19.9e %-19.9e\n",
            kPS[i],pPS[i],fkT[i],DplusT[i]);
    }
    fclose(outstr);
}

real Interpolation_local(real k, double kPS[], double pPS[], int nPS, double pPS2[]);



local void PSLTableNW(void)
{
    char namebuf[256];
    stream outstr;
    double ksmin;
    double ksmax;
    double dk;
    int Nks;
    int cutmin;
    int cutmax;
    int i;
    double *kPSNW;
    double *pPSNW;
    double *pPSNW2;

    double *kPSNWtmp;
    double *pPSNWtmp;
    double *pPSNWtmp2;

    double *logkpkT;
    double *logkpkTeven;
    double *logkpkTodd;

    double *mlogkpkTevencutted;
    double *logkpkTevencutted;
    double *logkpkTevencutted2;

    double *mlogkpkToddcutted;
    double *logkpkToddcutted;
    double *logkpkToddcutted2;

    double *preT;

    double *logkpkTsave;
    double *logkpkTsave2;
    double *logkpkTOut;
    double *logkpkTIm;
    double *logkpkTImOut;
    double PSL;
    int FSTtype;
    int FSTsizestep;
// THESE VECTOR ARE NOT USED. REMOVE
//    int nPSTabletmp;
//    real *kPStmp;
//    real *pPStmp;
//    real *pPStmp2;
// ALSO THESE LINES. REMOVE
//    char fnamePSPath[100];
//    string fnamePS;
//
    int j;
    int Nksevencutted;
    int Nksoddcutted;

    double outkmin, outkmax, doutk;
    int outsize;
    double kval;
    double ki;
    double kPSNW8;
    double PNWval;
    double PNW;
    double PSLkval;
    double DeltaAppf;
    double fac1;
    real aTime;
    aTime = second();

    outkmin = 0.00001;
    outkmax = 10;
    outsize = nPSLT; //800;

    ksmin = 7.0e-5 / cmd.h;
    ksmax = 7.0 / cmd.h;
    Nks = 65536;            // 2^16
    cutmin = 120;
    cutmax = 240;

    if (Nks==1)
        dk = 0.;
    else
        dk = (ksmax - ksmin)/((double)(Nks - 1));

    fprintf(gd.outlog," -> PSLTableNW: dk value: %g\n",dk);

    logkpkT = dvector(1,Nks);

// THESE VECTOR ARE NOT USED. REMOVE
//    kPStmp = dvector(1,pmgpt->nPSLT);
//    pPStmp = dvector(1,pmgpt->nPSLT);
//    pPStmp2 = dvector(1,pmgpt->nPSLT);

    logkpkTsave = dvector(1,Nks);
    logkpkTsave2 = dvector(1,Nks);
    logkpkTeven = dvector(1,Nks/2);
    logkpkTodd = dvector(1,Nks/2);

    Nksevencutted = Nks/2 - (cutmax-cutmin)-4;
    Nksoddcutted = Nks/2 - (cutmax-cutmin)-1;

    fprintf(gd.outlog," -> PSTableNW : numbers discarded %d %d\n", Nksevencutted, Nksoddcutted);

    mlogkpkTevencutted = dvector(1,Nksevencutted);
    logkpkTevencutted = dvector(1,Nksevencutted);
    logkpkTevencutted2 = dvector(1,Nksevencutted);

    mlogkpkToddcutted = dvector(1,Nksoddcutted);
    logkpkToddcutted = dvector(1,Nksoddcutted);
    logkpkToddcutted2 = dvector(1,Nksoddcutted);

    preT = dvector(1,Nks);

    kPSNWtmp = dvector(1,Nks);
    pPSNWtmp = dvector(1,Nks);
    pPSNWtmp2 = dvector(1,Nks);

    kPSNW = dvector(1,outsize);
    pPSNW = dvector(1,outsize);
    pPSNW2 = dvector(1,outsize);


    for (i=1; i<=Nks; i++) {
        kval = ksmin + (double)(i-1)*dk;
        PSL = Interpolation_local(kval, kPS, pPS, nPSLT, pPS2);
        logkpkT[i] = rlog10(kval * PSL);
        logkpkTsave[i] = logkpkT[i];
    }

    dsinft(logkpkT, Nks);
        for (j=1; j<=Nks/2; j++) {
            logkpkTeven[j] = logkpkT[2*j];
            logkpkTodd[j] = logkpkT[2*j-1];
        }

    i=1;
    for (j=2; j<=cutmin-2; j++) {
        mlogkpkTevencutted[i] = (double)j;
        logkpkTevencutted[i] = logkpkTeven[j];
        ++i;
    }
    for (j=cutmax+2; j<=Nks/2; j++) {
        mlogkpkTevencutted[i] = (double)j;
        logkpkTevencutted[i] = logkpkTeven[j];
        ++i;
    }
    fprintf(gd.outlog," -> PSTableNW : even cutted %d\n", i-1);

// El primero es cero:
    i=1;
    for (j=1; j<=cutmin-1; j++) {
        mlogkpkToddcutted[j] = (double)j;
        logkpkToddcutted[j] = logkpkTodd[j];
        ++i;
    }
    for (j=cutmax+1; j<=Nks/2; j++) {
        mlogkpkToddcutted[i] = (double)j;
        logkpkToddcutted[i] = logkpkTodd[j];
        ++i;
    }
    fprintf(gd.outlog," -> PSTableNW : odd cutted %d\n", i-1);


    spline(mlogkpkTevencutted,logkpkTevencutted,Nksevencutted,1.0e30,1.0e30,logkpkTevencutted2);
    spline(mlogkpkToddcutted,logkpkToddcutted,Nksoddcutted,1.0e30,1.0e30,logkpkToddcutted2);

    double yy, dyy;
    unsigned long jj;

    for (i=1; i<=Nks/2; i++) {
        if(cutmin < i && i < cutmax) {
            locate(mlogkpkTevencutted, Nksevencutted, (double)(i+1), &jj);
            polint(&mlogkpkTevencutted[jj-2],&logkpkTevencutted[jj-2],5,(double)(i+1),&yy,&dyy);
            preT[2*i] = yy;
            locate(mlogkpkToddcutted, Nksoddcutted, (double)(i), &jj);
            polint(&mlogkpkToddcutted[jj-2],&logkpkToddcutted[jj-2],5,(double)(i),&yy,&dyy);
            preT[2*i-1] = yy;
        } else{
            preT[2*i] = logkpkT[2*i];
            preT[2*i-1] = logkpkT[2*i-1];
        }
    }

    dsinft(preT, Nks);

    fac1 = (2.0/(double) Nks);
    for (i=1; i<=Nks; i++) {
        preT[i] *= fac1;
    }

    for (i=1; i<=Nks; i++) {
        kPSNWtmp[i] = ksmin + (double)(i-1)*dk;
        pPSNWtmp[i] = rpow(10.0,preT[i])/kPSNWtmp[i];
    }
    spline(kPSNWtmp,pPSNWtmp,Nks,1.0e30,1.0e30,pPSNWtmp2);

    doutk = (rlog10(outkmax) - rlog10(outkmin))/((double)(outsize - 1));

    kPSNW8 = kPSNWtmp[8];

    //~ fprintf(stdout, "\n Nks=%d, nPSLT=%d, outsize=%d\n",Nks,nPSLT,outsize);

    pPS_nw = dvector(1,nPSLT);   //nuevo
    pPS2_nw = dvector(1,nPSLT); //nuevo


    for (i=1; i<=outsize; i++) {

        ki= kPS[i];

		if (ki <= 0.001) {
			PSL       = Interpolation_local(kPSNW8, kPS, pPS, nPSLT, pPS2);
			PNWval    = Interpolation_local(kPSNW8, kPSNWtmp, pPSNWtmp, Nks, pPSNWtmp2);
			PSLkval   = Interpolation_local(ki, kPS, pPS, nPSLT, pPS2);
            DeltaAppf = (kPS[i]*(PSL - PNWval) / PNWval) /kPSNW8;
            pPS_nw[i] = PSLkval/(DeltaAppf + 1);
		} else if (ki > 0.001 & ki< ksmax){
                pPS_nw[i] = Interpolation_local(ki, kPSNWtmp, pPSNWtmp, Nks, pPSNWtmp2);
        } else {
                pPS_nw[i] = Interpolation_local(ki, kPS, pPS, nPSLT, pPS2);
        }


    }





    free_dvector(pPSNW2,1,outsize);
    free_dvector(pPSNW ,1,outsize);
    free_dvector(kPSNW ,1,outsize);

    free_dvector(pPSNWtmp2,1,Nks);
    free_dvector(pPSNWtmp,1,Nks);
    free_dvector(kPSNWtmp,1,Nks);

    free_dvector(preT,1,Nks);

    free_dvector(logkpkToddcutted2,1,Nksoddcutted);
    free_dvector(logkpkToddcutted,1,Nksoddcutted);
    free_dvector(mlogkpkToddcutted,1,Nksoddcutted);

    free_dvector(logkpkTevencutted2,1,Nksevencutted);
    free_dvector(logkpkTevencutted,1,Nksevencutted);
    free_dvector(mlogkpkTevencutted,1,Nksevencutted);

    free_dvector(logkpkTodd,1,Nks/2);
    free_dvector(logkpkTeven,1,Nks/2);
    free_dvector(logkpkTsave2,1,Nks);
    free_dvector(logkpkTsave,1,Nks);

// THESE VECTOR ARE NOT USED. REMOVE
//    free_dvector(pPStmp2,1,pmgpt->nPSTable);
//    free_dvector(pPStmp,1,pmgpt->nPSTable);
//    free_dvector(kPStmp,1,pmgpt->nPSTable);

    free_dvector(logkpkT,1,Nks);

    fprintf(gd.outlog," -> PSLTableNW CPU time %g...\n",second()-aTime);
}




real Interpolation_local(real k, double kPS[], double pPS[], int nPS, double pPS2[])
{
    real psftmp;

//    if ( k < kPS[1] || k > kPS[nPS] )
//        fprintf(stdout,"\n\nInterpolation_nr_local: warning! :: k is out of range... %g %g %g\n",k, kPS[1], kPS[nPS]);

    splint(kPS,pPS,pPS2,nPS,k,&psftmp);

    return (psftmp);
}












local void GaussLegendrePoints(void)
{
    double dx,ss,xval;
    real x1=-1.0, x2=1.0;

    int i;
    double xx=0.0;

    pGL = (global_GL_ptr) allocate(sizeof(global_GL));

    nGL(pGL) = cmd.ngausslegpoints;
    xGL(pGL)=dvector(1,cmd.ngausslegpoints);
    wGL(pGL)=dvector(1,cmd.ngausslegpoints);
    x1GL(pGL) = -1.0;
    x2GL(pGL) = 1.0;

    fprintf(gd.outlog,"\nComputation of GLs...\n");
    gauleg(x1GL(pGL),x2GL(pGL),xGL(pGL),wGL(pGL),nGL(pGL));


    fprintf(gd.outlog,"\nEnd of Gauss-Legendre computing (StartRun)\n");

}





global  real psInterpolation_nr(real k, double kPS[], double pPS[], int nPS)
{
//    pointPSTableptr pf, pi;
    real psftmp;
//    real dps;

//    pi = PSLCDMtab;
//    pf = PSLCDMtab+nPSTable-1;

    splint(kPS,pPS,pPS2,nPS,k,&psftmp);

    return (psftmp);
}
