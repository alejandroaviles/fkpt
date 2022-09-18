# fkpt
Perturbation theory for Modified Theory using "fk"-Kernels


Alejandro Aviles (Conacyt/ININ, Mexico),
avilescervantes@gmail.com

#

fkpt is a code that computes the 1-loop redshift space power spectrum for tracers. 


Compile:

```
/MGPT/src$ make
```

Run: in the parent directory

```
/MGPT$ ./mgpt
```
This will compute the LCDM power spectra and correlation functions.

For help:

```
/MGPT$ ./mgpt -help
```


In help you can see how to change parameters, in the form [option]=[value], for example:

```
/MGPT$ ./mgpt om=0.3 h=0.7 mgm=HS fR0=1.0e-6 suffix=_F6z05 zout=0.5 fnamePS=pklin.dat
```

computes Hu-Sawicky f_R0 = -10^-6 with background cosmology h=0.7, Omega_m = 0.3, at z=0.5, for the input real space linear power spectrum pklin.dat. The output files will have a suffix _F6z05. The input pkl should be given in a two column (k,pkl) file in Mpc/h units at the desire output redshift. 

## References

The fkpt theory is based on 

https://arxiv.org/abs/2208.02791


