#!/usr/bin/env python
# coding: utf-8

# In[ ]:


### Package to parse the fkpt code from Alejandro Aviles
import numpy as np
import sys, platform, os, subprocess
#from classy import Class

path_fkpt='../' 

#Definitions
def generate_ps(proc = os.getpid(), h = 0.6711, ombh2 = 0.022, omch2 = 0.122, omnuh2 = 0.0006442, 
                As = 2e-9, ns = 0.965, z = 0.97, N_ur = 2.0328,
                khmin = 0.0001, khmax = 2.0, nbk = 1000):
    '''Generates the linear (cb) power spectrum using Class.
    
    Args:
        proc: number of process,
        h = H0/100, with H0 the Hubble constant,
        omXh2: Omega_X h², where X = b (baryons), c (CDM), nu (neutrinos),
        As: amplitude of primordial curvature fluctuations,
        ns: spectral index,
        z: redshift,
        khmin, khmax: minimal and maximal wave-number,
        nbk: number of points in [khmin, khmax].
        
    Rertuns:
        kh: vector of wave-number,
        z: redshift,
        pk: linear (cb) power spectrum,
        Om: Omega matter value.
        proc: number of process.
    '''
    
    params = {
             'output':'mPk',
             'omega_b':ombh2,
             'omega_cdm':omch2,
             'omega_ncdm':omnuh2, 
             'h':h,
             'A_s':As,
             'n_s':ns,
             'P_k_max_1/Mpc':khmax,
             'z_max_pk':10.,#Default value is 10 
             #'tau_reio':0.09, #0.066  
             #'N_eff':Neff,
             'N_ur':N_ur,      #massless neutrinos 
             'N_ncdm':1 #massive neutrinos species
             }
    
    cosmo = Class()
    cosmo.set(params)
    cosmo.compute()
    
    #Specify k
    k = np.logspace(np.log10(khmin*h), np.log10(khmax*h), num = nbk) #Mpc^-1
    
    #Computes the linear (cb) power spectrum
    Plin = np.array([cosmo.pk_cb(ki, z) for ki in k])
    
    #Tranforming to h/Mpc and (Mpc/h)^3
    k /= h
    Plin *= h**3

    np.savetxt(path_fkpt+'Python/Input/ps_'+str(proc)+'_z'+str(z)+'.txt', np.array( (k, Plin)).T)
    
    #Omega_m
    Om = (ombh2 + omch2 + omnuh2)/h**2
    
    return({'kh':k, 'z':z, 'pk':Plin, 'Om':Om, 'proc':proc})




def run_fkpt(proc = os.getpid(), pk_name='', Om = 0.317501, h = 0.6711, zout = 0.5, 
             b1 = 1.7167, b2 = -0.3546, kmin = 1.0e-3, kmax = 0.5, Nk = 120,
             bs2 = -0.409543, b3nl = 0.0728076, alpha0 = -6.6232, alpha2 = -16.3094, alpha4 = -20.6674, 
             ctilde =-0.0, pshotp = 4809.502114, alpha0shot = 0.0805, alpha2shot = -4.3887,
             nquadSteps = 300, model='LCDM', fR0=1.0e-10, suffixModel= '', remove=True):
    
    
    args_exe = [path_fkpt+'fkpt', 'Om=%f'%Om ,'h=%f'%h  , 'b1=%f'%b1, 'b2=%f'%b2, 'bs2=%f'%bs2, 
                'b3nl=%f'%b3nl ,'alpha0=%f'%alpha0 , 'alpha2=%f'%alpha2, 'alpha4=%f'%alpha4, 
                'ctilde=%f'%ctilde, 'pshotp=%f'%pshotp, 'alpha0shot=%f'%alpha0shot, 'alpha2shot=%f'%alpha2shot,
                'zout=%f'%zout, 'kmin=%f'%kmin, 'Nk=%f'%Nk,
                'kmax=%f'%kmax, 'suffixModel=_%d'%proc,
                'model='+str(model)+'', 'fR0=%1.15f'%fR0]
    #print(args_exe)
    
    if pk_name != '': args_exe.append('fnamePS=%s'%pk_name)
    res_proc = subprocess.call(args_exe)
    
    ### Part to read the data and then remove the file 
    name_read = path_fkpt+'Python/rsd_multipoles'+'_%d'%proc+'.dat'
    data = np.loadtxt(name_read, skiprows=3).T
    
    
    ###create folder for multipoles
    os.makedirs(path_fkpt+'Python/Output/multipoles', exist_ok=True)
    
    ###rename outputs files: to be identified by its pid process and redshift
    os.rename(path_fkpt+'Python/rsd_multipoles'+'_%d'%proc+'.dat', 
              path_fkpt+'Python/Output/multipoles/rsd_multipoles_'+str(proc)+'_z'+str(zout)+'.dat')
    
    os.rename(path_fkpt+'Python/Output/linear_%d'%proc+'.dat',
              path_fkpt+'Python/Output/linear_'+str(proc)+'_z'+str(zout)+'.dat')
    
    os.rename(path_fkpt+'Python/Output/kfunctions_%d'%proc+'.dat',
              path_fkpt+'Python/Output/kfunctions_'+str(proc)+'_z'+str(zout)+'.dat')
    
    os.rename(path_fkpt+'Python/Output/kfunctions_nw_%d'%proc+'.dat',
              path_fkpt+'Python/Output/kfunctions_nw_'+str(proc)+'_z'+str(zout)+'.dat')
    
    #if remove==True:
        ## Remove the rsd results
        #args_rm = ['rm', 'rsd_multipoles'+'_%d'%proc+'.dat']
        #subprocess.call(args_rm)
        ## Remove all the intermediate results in Output
        #args_rm = ['rm', 'Output/linear_%d'%proc+'.dat', 'Output/kfunctions_%d'%proc+'.dat','Output/kfunctions_nw_%d'%proc+'.dat' ]
        #subprocess.call(args_rm)
        ## Remove the linear power spectrum in Input
        #args_rm = ['rm', 'Input/ps_%d'%proc+'.txt' ]
        #subprocess.call(args_rm)
    
    
    return({'k':data[0], 'mono':data[1], 'quad':data[2], 'hexa':data[3]})