
// Alejandro Aviles (avilescervantes@gmail.com)


#include "globaldefs.h"
#include "protodefs.h"

#define FMTQFUNCTIONSDAT    \
"%e %e %e %e %e %e %e \
 %e %e %e %e %e %e %e \
 %e %e %e %e %e %e %e\n"

local void print_kfunctions(void);
local void print_kfunctions_nw(void);
local void print_linear(void);


global void write(void){
	


	
	    fprintf(stdout,"\nz = %g\n",cmd.xstop);
	    fprintf(stdout,"Dplus = %g\n",gd.Dplus);
	    fprintf(stdout,"f0 = %g\n",gd.f0);    
	    fprintf(stdout,"2*sigma_psi = %g  Mpc/h (Lagrangian particles mean displacement)\n",
							gd.particles_meanPath); 
	
	print_linear();
										
	print_kfunctions();										
	print_kfunctions_nw();
	
	



	
	
	
}
#undef FMTQFUNCTIONSDAT





#define FMTKOUT    \
"%e %e %e %e %e %e %e \
%e %e %e %e %e %e %e \
%e %e %e %e %e %e %e \
%e %e %e %e %e %e %e \
%e %e\n"

void print_kfunctions(void){

	stream outstr;
	outstr = stropen(gd.fpfnamekfun,"w!");

    fprintf(outstr,"# InputPklFile=%s, redshift z=%g\n",
                cmd.fnamePS, cmd.xstop);                     
    fprintf(outstr,"# Precision: k-functions: nquadStetps=%d\n",
                cmd.nquadSteps); 
	fprintf(outstr,"%1s%8s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s \
		%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s",
            "#","1.k","2.P22dd","3.P22du","4.Pdduu",
            "5.P13dd","6.P13du","7.P13uu",
            "8.I1udd1A", "9.I2uud1A", "10.I2uud2A", "11.I3uuu2A", "12.I3uuu3A",
            "13.I2uudd1D", "14.I2uudd2D", "15.I3uuud2D", "16.I3uuud3D",
            "17.I4uuuu2", "18.I4uuuu3D", "19.I4uuuu4D",
            "20.Pb1b2", "21.Pb1bs2", "22.Pb22", "23.Pb2s2",
            "24.Ps22", "25.Pb2theta", "26.Pbs2theta",
            "27.sigma32pkl","28.pklinear", "29.fk", "30.f0\n");


   for (int i=0; i<golists.Nk; i++) {
	    fprintf(outstr,FMTKOUT ,
			kFArrays.kT[i],
			kFArrays.P22ddT[i],
			kFArrays.P22duT[i],
			kFArrays.P22uuT[i],
			kFArrays.P13ddT[i],
			kFArrays.P13duT[i],
			kFArrays.P13uuT[i],
			kFArrays.I1udd1AT[i],
			kFArrays.I2uud1AT[i],
			kFArrays.I2uud2AT[i],
			kFArrays.I3uuu2AT[i],
			kFArrays.I3uuu3AT[i],  
			kFArrays.I2uudd1BpCT[i],
			kFArrays.I2uudd2BpCT[i],
			kFArrays.I3uuud2BpCT[i],
			kFArrays.I3uuud3BpCT[i],
			kFArrays.I4uuuu2BpCT[i],
			kFArrays.I4uuuu3BpCT[i],
			kFArrays.I4uuuu4BpCT[i],
			kFArrays.Pb1b2T[i],
			kFArrays.Pb1bs2T[i],
			kFArrays.Pb22T[i],
			kFArrays.Pb2s2T[i],
			kFArrays.Ps22T[i],
			kFArrays.Pb2thetaT[i],
			kFArrays.Pbs2thetaT[i],
			kFArrays.sigma32PSLT[i],
			kFArrays.pklT[i],
			kFArrays.fkT[i],
			gd.f0
		);	
	};
    fclose(outstr);	
	
};






void print_kfunctions_nw(void){

	stream outstr;

	outstr = stropen(gd.fpfnamekfunnw,"w!");

    fprintf(outstr,"# InputPklFile=%s, redshift z=%g\n",
                cmd.fnamePS, cmd.xstop);                     
    fprintf(outstr,"# Precision: k-functions: nquadStetps=%d\n",
                cmd.nquadSteps); 
	fprintf(outstr,"%1s%8s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s \
		%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s",
            "#","1.k","2.P22dd","3.P22du","4.Pdduu",
            "5.P13dd","6.P13du","7.P13uu",
            "8.I1udd1A", "9.I2uud1A", "10.I2uud2A", "11.I3uuu2A", "12.I3uuu3A",
            "13.I2uudd1D", "14.I2uudd2D", "15.I3uuud2D", "16.I3uuud3D",
            "17.I4uuuu2", "18.I4uuuu3D", "19.I4uuuu4D",
            "20.Pb1b2", "21.Pb1bs2", "22.Pb22", "23.Pb2s2",
            "24.Ps22", "25.Pb2theta", "26.Pbs2theta",
            "27.sigma32pkl","28.pklinear", "29.fk", "30.f0\n");


   for (int i=0; i<golists.Nk; i++) {
	    fprintf(outstr,FMTKOUT ,
			kFArrays_nw.kT[i],
			kFArrays_nw.P22ddT[i],
			kFArrays_nw.P22duT[i],
			kFArrays_nw.P22uuT[i],
			kFArrays_nw.P13ddT[i],
			kFArrays_nw.P13duT[i],
			kFArrays_nw.P13uuT[i],
			kFArrays_nw.I1udd1AT[i],
			kFArrays_nw.I2uud1AT[i],
			kFArrays_nw.I2uud2AT[i],
			kFArrays_nw.I3uuu2AT[i],
			kFArrays_nw.I3uuu3AT[i],  
			kFArrays_nw.I2uudd1BpCT[i],
			kFArrays_nw.I2uudd2BpCT[i],
			kFArrays_nw.I3uuud2BpCT[i],
			kFArrays_nw.I3uuud3BpCT[i],
			kFArrays_nw.I4uuuu2BpCT[i],
			kFArrays_nw.I4uuuu3BpCT[i],
			kFArrays_nw.I4uuuu4BpCT[i],
			kFArrays_nw.Pb1b2T[i],
			kFArrays_nw.Pb1bs2T[i],
			kFArrays_nw.Pb22T[i],
			kFArrays_nw.Pb2s2T[i],
			kFArrays_nw.Ps22T[i],
			kFArrays_nw.Pb2thetaT[i],
			kFArrays_nw.Pbs2thetaT[i],
			kFArrays_nw.sigma32PSLT[i],
			kFArrays_nw.pklT[i],
			kFArrays_nw.fkT[i],
			gd.f0
		);	
	};
    fclose(outstr);	
	
};






void free_variables(void){

	free(kFArrays.kT);
	free(kFArrays.P22ddT);
	free(kFArrays.P22duT);
	free(kFArrays.P22uuT);
		// A 
	free(kFArrays.I1udd1AT);
	free(kFArrays.I2uud1AT);
	free(kFArrays.I2uud2AT);
	free(kFArrays.I3uuu2AT);
	free(kFArrays.I3uuu3AT);
	//  B plus C   
	free(kFArrays.I2uudd1BpCT);
	free(kFArrays.I2uudd2BpCT);
	free(kFArrays.I3uuud2BpCT);
	free(kFArrays.I3uuud3BpCT);
	free(kFArrays.I4uuuu2BpCT);
	free(kFArrays.I4uuuu3BpCT);
	free(kFArrays.I4uuuu4BpCT);
	//  Bias
	free(kFArrays.Pb1b2T);
	free(kFArrays.Pb1bs2T);
	free(kFArrays.Pb22T);
	free(kFArrays.Pb2s2T);
	free(kFArrays.Ps22T);
	free(kFArrays.Pb2thetaT);
	free(kFArrays.Pbs2thetaT);
	//
	free(kFArrays.P13ddT);
	free(kFArrays.P13duT);
	free(kFArrays.P13uuT);
	free(kFArrays.sigma32PSLT);
	free(kFArrays.pklT);
	free(kFArrays.fkT);
	
	
	free(kFArraysd.kT);
	free(kFArraysd.P22ddT);
	free(kFArraysd.P22duT);
	free(kFArraysd.P22uuT);
		// A 
	free(kFArraysd.I1udd1AT);
	free(kFArraysd.I2uud1AT);
	free(kFArraysd.I2uud2AT);
	free(kFArraysd.I3uuu2AT);
	free(kFArraysd.I3uuu3AT);
	//  B plus C   
	free(kFArraysd.I2uudd1BpCT);
	free(kFArraysd.I2uudd2BpCT);
	free(kFArraysd.I3uuud2BpCT);
	free(kFArraysd.I3uuud3BpCT);
	free(kFArraysd.I4uuuu2BpCT);
	free(kFArraysd.I4uuuu3BpCT);
	free(kFArraysd.I4uuuu4BpCT);
	//  Bias
	free(kFArraysd.Pb1b2T);
	free(kFArraysd.Pb1bs2T);
	free(kFArraysd.Pb22T);
	free(kFArraysd.Pb2s2T);
	free(kFArraysd.Ps22T);
	free(kFArraysd.Pb2thetaT);
	free(kFArraysd.Pbs2thetaT);
	//
	free(kFArraysd.P13ddT);
	free(kFArraysd.P13duT);
	free(kFArraysd.P13uuT);
	free(kFArraysd.sigma32PSLT);
	free(kFArraysd.pklT);
	free(kFArraysd.fkT);	
	
	
	
	free(kFArrays_nw.kT);
	free(kFArrays_nw.P22ddT);
	free(kFArrays_nw.P22duT);
	free(kFArrays_nw.P22uuT);
		// A 
	free(kFArrays_nw.I1udd1AT);
	free(kFArrays_nw.I2uud1AT);
	free(kFArrays_nw.I2uud2AT);
	free(kFArrays_nw.I3uuu2AT);
	free(kFArrays_nw.I3uuu3AT);
	//  B plus C   
	free(kFArrays_nw.I2uudd1BpCT);
	free(kFArrays_nw.I2uudd2BpCT);
	free(kFArrays_nw.I3uuud2BpCT);
	free(kFArrays_nw.I3uuud3BpCT);
	free(kFArrays_nw.I4uuuu2BpCT);
	free(kFArrays_nw.I4uuuu3BpCT);
	free(kFArrays_nw.I4uuuu4BpCT);
	//  Bias
	free(kFArrays_nw.Pb1b2T);
	free(kFArrays_nw.Pb1bs2T);
	free(kFArrays_nw.Pb22T);
	free(kFArrays_nw.Pb2s2T);
	free(kFArrays_nw.Ps22T);
	free(kFArrays_nw.Pb2thetaT);
	free(kFArrays_nw.Pbs2thetaT);
	//
	free(kFArrays_nw.P13ddT);
	free(kFArrays_nw.P13duT);
	free(kFArrays_nw.P13uuT);
	free(kFArrays_nw.sigma32PSLT);
	free(kFArrays_nw.pklT);
	free(kFArrays_nw.fkT);
	
	
	free(kFArraysd_nw.kT);
	free(kFArraysd_nw.P22ddT);
	free(kFArraysd_nw.P22duT);
	free(kFArraysd_nw.P22uuT);
		// A 
	free(kFArraysd_nw.I1udd1AT);
	free(kFArraysd_nw.I2uud1AT);
	free(kFArraysd_nw.I2uud2AT);
	free(kFArraysd_nw.I3uuu2AT);
	free(kFArraysd_nw.I3uuu3AT);
	//  B plus C   
	free(kFArraysd_nw.I2uudd1BpCT);
	free(kFArraysd_nw.I2uudd2BpCT);
	free(kFArraysd_nw.I3uuud2BpCT);
	free(kFArraysd_nw.I3uuud3BpCT);
	free(kFArraysd_nw.I4uuuu2BpCT);
	free(kFArraysd_nw.I4uuuu3BpCT);
	free(kFArraysd_nw.I4uuuu4BpCT);
	//  Bias
	free(kFArraysd_nw.Pb1b2T);
	free(kFArraysd_nw.Pb1bs2T);
	free(kFArraysd_nw.Pb22T);
	free(kFArraysd_nw.Pb2s2T);
	free(kFArraysd_nw.Ps22T);
	free(kFArraysd_nw.Pb2thetaT);
	free(kFArraysd_nw.Pbs2thetaT);
	//
	free(kFArraysd_nw.P13ddT);
	free(kFArraysd_nw.P13duT);
	free(kFArraysd_nw.P13uuT);
	free(kFArraysd_nw.sigma32PSLT);
	free(kFArraysd_nw.pklT);
	free(kFArraysd_nw.fkT);
	
};


void print_linear(void){
	
	stream outstr;
    char namebuf[256];


    sprintf(namebuf,"%s/%s%s%s",gd.clptDir,"linear",cmd.suffixModel,".dat");
    outstr = stropen(namebuf,"w!");
 
    fprintf(outstr,"# InputPklFile=%s\n", cmd.fnamePS); 
    fprintf(outstr,"# z=%g, OmegaM=%g, h=%g\n",
		cmd.xstop, cmd.om, cmd.h);                    
    fprintf(outstr,"%1s%11s%17s%20s%22s%22s",
            "#"," 1.k[h/Mpc]","2.pkl","3.f(k)","4.D+(k)","5.pkl_nw\n");
   
    for (int i=1; i<=nPSLT; i++) {
        fprintf(outstr,"%-19.9e %-19.9e %-19.9e %-19.9e %-19.9e\n",
            kPS[i],
            pPS[i],
            fkT[i],
            DplusT[i],
            pPS_nw[i]);
    }
    
    fclose(outstr);	


	
};
