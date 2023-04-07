# fkpt
Perturbation theory for Modified Gravity Theory using "fk"-Kernels.

fkpt is a code that computes the 1-loop redshift space power spectrum for tracers. 


## Author: 

Alejandro Aviles (Conacyt/ININ, Mexico),
avilescervantes@gmail.com

**Other people who contributed to this code**:
- [Hernán E. Noriega](https://github.com/henoriega)
- Mario A. Rodriguez-Meza
- Sebastien Fromenteau
#



HS is the only modified gravity model implemented so far. It is straightforward to do it for other models, hope I can find the time to do it soon. Contact me if you have questions about how to do it.



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
- NumPy 
- SciPy

Once compiled with everything ready, check the [Jupyter Notebook](https://github.com/alejandroaviles/fkpt/blob/main/Python/run_fkpt.ipynb) which contains a helpful example. 


## References

The fkpt theory is based on: 

https://arxiv.org/abs/2012.05077

https://arxiv.org/abs/2208.02791

Please cite these papers, at least the first one, if you use this code for you research.


