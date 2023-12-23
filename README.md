# fkpt
Perturbation theory for LCDM and Modified Gravity theories using "fk"-Kernels.

fkpt is a code that computes the 1-loop redshift space power spectrum for tracers. 


## Authors: 

- Alejandro Aviles (ICF-UNAM, Mexico), avilescervantes@gmail.com
- Mario A. Rodriguez-Meza (ININ, Mexico)


Special thanks to Sebastien Fromenteau and [Hernán E. Noriega](https://github.com/henoriega) for helping with the Python interface.
#

Hu-Sawicky f(R) is the only modified gravity model implemented so far. It is straightforward to do it for other models, hope I can find the time to do it soon. Contact me if you have questions about how to do it.

If you run the LCDM model, it will use the exact time dependence for the kernels. 



## Instructions for C

First use git clone:

```
git clone https://github.com/alejandroaviles/fkpt.git
```

**Compile:**

```
/fkpt$ make
```

**Run:** 

```
/fkpt$ ./fkpt
```
This will compute the LCDM redshift space multipoles of the power spectrum

For help:

```
/fkpt$ ./fkpt -help
```


In help you can see how to change parameters, in the form [option]=[value], for example:

```
/fkpt$ ./fkpt Om=0.3 h=0.7 model=HS fR0=1.0e-6 suffix=_F6z05 zout=0.5 fnamePS=pklin.dat
```

computes Hu-Sawicky f_R0 = -10^-6 with background cosmology h=0.7, Omega_m = 0.3, at z=0.5, for the input real space linear power spectrum pklin.dat. The output files will have a suffix _F6z05. 

The input linear power spectrum should be the LCDM one and given in a two column (k,pkl) file in Mpc/h units at the desire output redshift.

## Instructions for Python

**Dependences**

The code employs the standard libraries:
- **[NumPy](https://numpy.org/)**
- **[SciPy](https://scipy.org/)**


Once compiled with everything ready, check the [Jupyter Notebook](https://github.com/alejandroaviles/fkpt/blob/main/Python/run_fkpt.ipynb) which contains a helpful example. 

AP is implemented for Python only


## References

Please cite https://arxiv.org/abs/2312.10510 if you use this code for you research. 





