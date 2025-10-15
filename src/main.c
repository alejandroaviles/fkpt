/*==============================================================================
 NAME: main.c				[code for redshift space correlation function - GSM]
 Alejandro Aviles (avilescervantes@gmail.com), ...
 *  (other people who collaborated: Mario A. Rodriguez-Meza ...)
 ================================================================================ 
 Use: ./fkpt -help
 References:  arXiv:...
 ==============================================================================*/

#ifndef FKPT_PYEXT   // building the CLI => own the globals
#define global
#endif

#include "globaldefs.h"
#include "cmdline_defs.h"
#include "protodefs.h"
#include "models.h"

int main(int argc, string argv[])
{
    gd.cpuinit = second();
    InitParam(argv, defv);
    StartRun(argv[0], HEAD1, HEAD2, HEAD3);
    MainLoop();
    EndRun();
    return 0;
}



void MainLoop(void)
{
    global_variables();
    compute_kfunctions();
    compute_rsdmultipoles();
    write();
    free_variables();	
}
