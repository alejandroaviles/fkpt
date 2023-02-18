/*==============================================================================
 MODULE: gsm_io.c		[fkpt]
 ==============================================================================*/

#include "globaldefs.h"
#include "protodefs.h"


void StartOutput(void)
{


    fprintf(gd.outlog,"\n%8s%8s%8s", "maxnsteps", "etaini", "deta");
    fprintf(gd.outlog,"%8s\n","etaout");
    fprintf(gd.outlog,"%8d%8.2f%8.4f%8.4f\n",cmd.maxnsteps,cmd.x,gd.dx,gd.xstop);
    if (! strnull(cmd.options))
        fprintf(stdout,"\n\toptions: %s\n", cmd.options);
    
}

// I/O directories:
global void setFilesDirs_log(void)
{
    char buf[200];
    
    
    sprintf(gd.tmpDir,"tmp");    
    sprintf(buf,"if [ ! -d %s ]; then mkdir %s; fi",gd.tmpDir,gd.tmpDir);
    system(buf);    
    sprintf(gd.logfilePath,"%s/gsm%s.log",gd.tmpDir,cmd.suffixModel);
    return;
}

global void setFilesDirs(void)
{
    char buf[200];
    
    sprintf(gd.clptDir,"Output");
    sprintf(buf,"if [ ! -d %s ]; then mkdir %s; fi",gd.clptDir,gd.clptDir);
    fprintf(gd.outlog,"system: %s\n",buf);
    system(buf);

    //~ sprintf(gd.fpfnameTables,"Output/tables%s.dat",cmd.suffixModel); 
    //~ sprintf(gd.fpfnameParams,"Output/params%s.dat",cmd.suffixModel);   
    sprintf(gd.fpfnamekfun,"Output/kfunctions%s.dat",cmd.suffixModel);  
    sprintf(gd.fpfnamekfunnw,"Output/kfunctions_nw%s.dat",cmd.suffixModel);
    //~ sprintf(gd.fpfnamekfun2,"kfunctions%s_2.dat",cmd.suffixModel);
    sprintf(gd.fpfnameSPTPowerSpectrum,"SPTPowerSpectrum%s.dat",cmd.suffixModel);
    sprintf(gd.fpfnameqfunctions,"Output/qfunctions%s.dat",cmd.suffixModel);
    //~ sprintf(gd.fpfnameqfunctions2,"qfunctions%s_2.dat",cmd.suffixModel);
    sprintf(gd.fpfnameclptfunctions,"Output/clpt%s.dat",cmd.suffixModel);
    //~ sprintf(gd.fpfnameclptfunctions2,"clpt%s_2.dat",cmd.suffixModel);    
    sprintf(gd.fpfnamev12,"v12%s.dat",cmd.suffixModel);
    sprintf(gd.fpfnamesigma12_parallel,"sigma12_parallel%s.dat",cmd.suffixModel);
    sprintf(gd.fpfnamesigma12_perp,"sigma12_perp%s.dat",cmd.suffixModel);
    sprintf(gd.fpfnamersd,"rsd_multipoles%s.dat",cmd.suffixModel);
    
    return;
}

void EndRun(void)
{
    char   buf[200];
	FILE *fd;

	fclose(gd.outlog);
    
    if(cmd.chatty==1 || cmd.chatty==2) printf("\nTotal time : %g seconds\n\n", second() - gd.cpuinit);
    return;
}



