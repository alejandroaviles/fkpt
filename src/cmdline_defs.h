/*==============================================================================
 HEADER: cmdline_defs.h		[gsm]    
 Alejandro Aviles (other collaborators: Mario A. Rodriguez-Meza ...)
 * ==============================================================================
*/

#ifndef _cmdline_defs_h
#define _cmdline_defs_h

#define HEAD1	""
#define HEAD2	"GSM code."
#define HEAD3	"..."
   //~ fprintf(outstr,"InputPklFile=%s, redshift z=%g, OmegaM=%g, h=%g, b1=%g, \
                         //~ b2=%g, bs=%g, c1eft=%g, c2eft=%g, s2FoG=%g,\
                         //~ For Precision: q-funtcionts: NqperLogDecade=%g, Nk_qFunctionsQuad=%g, \
                         //~ gsm: gsm_width=%g,gsm_sizeyT=%d, gsm_NGL=%d\n",
                //~ cmd.fnamePS, cmd.xstop, cmd.om, cmd.h, cmd.b1, cmd.b2, cmd.bs, cmd.c1eft, cmd.c2eft,cmd.s2eft,
               //~ cmd.NqperLogDecade,cmd.Nk_qFunctionsQuad, cmd.gsm_width,cmd.gsm_sizeyT,cmd.gsm_NGL);
string defv[] = {  ";"HEAD1": " HEAD2 "\n\t " HEAD3,
    "paramfile=",                   ";Parameter input file. Overwritten by what follows",
//
// Power spectrum table:
    "fnamePS=pkl_z05.dat",    ";Input filename power spectrum table (k,P(k)). (At redshift zout)",
    "zout=0.5",                         ";Output redshift value",":z",
    //~ "Om=0.281",                     ";Omega matter value (z=0)",":OmegaM",
    "Om=0.317501",                     ";Omega matter value (z=0)",":OmegaM",
    //~ "h=0.697",                         ";Hubble parameter",
    "h=0.6711",                         ";Hubble parameter",
// bias parameters:
    "b1=1.7167",                     ";Linear local bias",
    "b2=-0.3546",                     ";Second order local bias",
    "bs2=-0.409543",                     ";tidal bias",
    "b3nl=0.0728076",                     ";b3nl bias",
    "alpha0=-6.6232",                     ";alpha0 EFT parameter",
    "alpha2=-16.3094",                     ";alpha2 EFT parameter",
    "alpha4=-20.6674",                     ";alpha4 EFT parameter",
    "ctilde=-0.0",                     ";NLO EFT parameter. It is degenerate with alpha2shot, hence set one of them to zero",
    "pshotp=4809.502114",                     "; Poissonian Shot Noise. Pshot = pshotp *(alpha0shot + alpha2shot mu^2 k^2). You can fix it to 1 and vary alpha0 shot noise",
    "alpha0shot=0.0805",                     ";alpha0 shot noise Pshot = pshotp *(alpha0shot + alpha2shot mu^2 k^2)",
    "alpha2shot=-4.3887",                     ";alpha2 shot noise Pshot = pshotp *(alpha0shot + alpha2shot mu^2 k^2)",
// output gsm multipoles
    //~ "smin=1",                     ";Output rsd 2pcf multipoles s minimum",
    //~ "smax=130",                     ";Output rsd 2pcf multipoles s maximum",
    //~ "Ns=100",                     ";number of s in rsd 2pcf multipoles output",
    //~ "c1eft=0.0",                      ";eft counterterm in xi(r), (almost) degenerate with nabla2 bias",
    //~ "c2eft=0.0",                      ";eft counterterm in pairwise velocity infall (not implemented yet)",
    //~ "sFoG=-2.0",                     ";eft counterterm in the pairwise velocity dispersion (similar to sigma_FoG)",":sigma2eft",
//  k-functions
    "kmin=1e-3",                    ";kmin in output",
    "kmax=0.5",                     ";kmax in output",
    "Nk=120",                       ";Total number of output",":nk",
    "nquadSteps=300",               ";Number of internal momenta p loop quadratures (trapezoid)",":nquad",
// model
    "model=LCDM",                  ";LCDM, HS, DGP",":m",    
    "fR0=1.0e-10",                 "; HS fR0 parameter",":fr0",
    "suffixModel=",                 ";Suffix model to add to output filenames", ":suffix",
    "modelParamfile=",              ";If mgmodel=USER, to use the model in models_user.h", ":mpf",
    NULL,
};

#endif // ! _cmdline_defs_h





/*
//~ // CLPT correlation functions table:
    "rmin=1",                      "; NOT USED rmin of the range for CLPT correlation functions NOT USED",
    "rmax=210",                     "; NOT USED rmax of the range for CLPT correlation functions NOT USED",
    "Nr=210",                       "; NOT USED Total number of r´s in the CLPT correlation function",":nr",
//~ // Modified gravity model parameters:
    "mgModel=LCDM",                 ";Modified gravity model to study (HS or DGP), default is LCDM", ":mgm",
    //~ "suffixModel=",                 ";Suffix model to add to output filenames", ":suffix",
    "fR0=1.0e-5",                   ";Hu-Sawicky f_R0",
    "screening=1.0",                ";set to =0 if you want no screenings", ":sc",
//~ // DGP:
    "epsDGP=-1.0",                  ";epsilon DGP parameter, =-1 normal branch, =1 self-accelerating branch",":epsdgp",
    "rcDGP=1.0",                    ";crossover scale parameter in DGP, (in units of 1/H0)",":rcdgp",
//~ //
    //~ "modelParamfile=",              ";If mgmodel=USER, to use the model in models_user.h", ":mpf",
//~ //
//~ // Background cosmology:
    "Om=0.281",                     ";Omega matter value (z=0)",":om",
    "OL= 1 - Om",                   ";Omega Lambda value (z=0). Only works for DGP!",":ol",
    "h=0.697",                      ";Hubble parameter",
//~ //
//~ // Differential equations evolution parameters:
    "etaini=-4.0",                  ";Initial conformal time value :: Log[1/(1 + zini)]",
    "deta=2/5",                     ";Conformal time integration step",
    "detamin=0.",                   ";Min conformal time integration step size",
    "epsSolver=1.0e-4",             ";Differential equations solver tolerance error parameter",":epssolver",
    "maxnsteps=10000",              ";Maximum number of integration steps", ":maxn",
    "solverMethod=rkqs",	        ";Integration method to use", ":solver",
//~ //
//~ // Quadrature parameters:
    "quadratureMethod=trapezoid3",   ";Quadrature method to use", ":quadm",
    "ngausslegpoints=16",           ";Number of Gauss-Legendre of integration points", ":nglpts",
    "epsquad=1.0e-6",               ";Quadrature tolerance error parameter (open Romberg method: romberg)",
//~ //
//~ // Post processing parameters:
    "postprocessing=false",			";Post processing options", ":pp",
    "options=",                     ";Various control options", ":opt",
*/
