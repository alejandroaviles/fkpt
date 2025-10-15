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
    /* No on-disk logging: send logs to the bit bucket,
       and don’t create any directories. */
    gd.tmpDir[0] = '\0';
    snprintf(gd.logfilePath, sizeof(gd.logfilePath), "/dev/null");
}

/* Also make setFilesDirs() a no-op that disables outputs.
   We’ll guard the actual writes next (in startrun.c etc.). */
global void setFilesDirs(void)
{
    /* Tell downstream code “don’t write anywhere”. */
    gd.tmpDir[0]  = '\0';
    gd.clptDir[0] = '\0';

    /* If your build defines any precomputed output-path buffers,
       blank them so accidental writers won’t find paths. It’s safe to
       leave this block out if these fields don’t exist in your gd. */
#ifdef HAVE_GD_FILEPATHS
    gd.fpfnamekfun[0]            = '\0';
    gd.fpfnamekfunnw[0]          = '\0';
    gd.fpfnameqfunctions[0]      = '\0';
    gd.fpfnameclptfunctions[0]   = '\0';
    gd.fpfnameSPTPowerSpectrum[0]= '\0';
    gd.fpfnamersd[0]             = '\0';
#endif
}

void EndRun(void)
{
    char   buf[200];
	FILE *fd;

	fclose(gd.outlog);
    
    if(cmd.chatty==1 || cmd.chatty==2) printf("\nTotal time : %g seconds\n\n", second() - gd.cpuinit);
    return;
}



