# ---- imports (standard packages) ----
import os, tempfile
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy

import pyfkpt, inspect, importlib
print(pyfkpt.__file__)              # path to the loaded extension
import os; print(os.path.getmtime(pyfkpt.__file__))  # timestamp

# ---- pyfkpt (new fkpt wrapper) ----
import pyfkpt.rsd as pyfkpt
    
# ---- ISiTGR: linear matter P(k) at z_pk ----
import isitgr
from isitgr import model

# ---- cosmology & nuisance parameters ----
# cosmology (same as your example)
h      = 0.6711
ombh2  = 0.022
omch2  = 0.122
omnuh2 = 0.0 #0.0006442
As     = 2e-9
ns     = 0.965
z_pk   = 0.3
N_ur   = 1.0 + 2.0328
Om     = (ombh2 + omch2 + omnuh2) / (h**2)
print("Omega_m:", Om)

# k-range for CAMB
khmin, khmax, nbk = 1.0e-4, 2.0, 1000

# nuisance vector (b1, b2, bs2, b3nl, EFT, stoch, shot)
b1 = 1.70               
b2 = -0.45
bs2  = -4/7*(b1 - 1)
b3nl = 32/315*(b1 - 1)

alpha0, alpha2, alpha4, ctilde = 3.0, -29.0, 0.0, 0.0
alpha0shot, alpha2shot, pshotp = 0.08, -8.0, 5000.0

nuis = [b1,b2,bs2,b3nl,alpha0,alpha2,alpha4,ctilde,alpha0shot,alpha2shot,pshotp]

# --- ISiTGR setup: MG results ---
mu0=0.5
pars = isitgr.CAMBparams()
pars.set_cosmology(H0=h*100, ombh2=ombh2, omch2=omch2, mnu=93.14*omnuh2,
                   MG_parameterization="muSigma", mu0=mu0)
pars.InitPower.set_params(As=As, ns=ns)

# finite-difference step in ln a around z_pk
dln_a = 1e-3
a0 = 1.0/(1.0 + z_pk)
a1 = a0 * np.exp(-dln_a)
a2 = a0 * np.exp(+dln_a)
z1 = 1.0/a1 - 1.0
z2 = 1.0/a2 - 1.0

# include z=0 because your f_of_k uses P(z=0,k) for normalization
z_list = sorted([z1, z_pk, z2, 0.0])
pars.set_matter_power(redshifts=z_list, kmax=2.0)
pars.NonLinear = model.NonLinear_none

# --- run ISiTGR as you already do (with z_list including z1, z_pk, z2, 0.0) ---
results = isitgr.get_results(pars)

# 1) Use the cached interpolator to get the internal k-grid and P_lin(k,z_pk)
PK, _, ks = results.get_matter_power_interpolator(
    nonlinear=False, hubble_units=True, k_hunit=True, return_z_k=True, silent=True
)
pk_lin = PK.P(z_pk, ks)        # shape (nk_internal,)

# 2) f(k,z_pk) uses the same internal k-grid (ks)
kh_fk, fk_lin = results.f_of_k(z_pk, hubble_units=True, k_hunit=True)
# sanity: these should now match exactly
assert np.allclose(kh_fk, ks)

# 3) large-scale limit f0 (choose your small-k cut)
smallk = ks < 3e-3
f0 = float(np.mean(fk_lin[smallk])) if np.any(smallk) else float(fk_lin[0])

# Base FKPT params (shared)
base = dict(
    z=z_pk, Om=Om, h=h,
    b1=b1, b2=b2, bs2=bs2, b3nl=b3nl,
    alpha0=alpha0, alpha2=alpha2, alpha4=alpha4, ctilde=ctilde,
    PshotP=pshotp, alpha0shot=alpha0shot, alpha2shot=alpha2shot,
    kmin=float(max(1e-3, ks.min())),
    kmax=float(min(ks.max(), 0.5)),
    Nk=int(min(ks.size, 240)),
    nquadSteps=300, chatty=0,
    model="HDKI",    # must match models.c
    mg_variant="mu_OmDE",
    mu0=mu0
)

# evaluation grid
k_eval = np.linspace(base["kmin"], base["kmax"], base["Nk"])

print("nPSLT from pyfkpt:", ks.shape)
# ==========================
# (4) MG (unrescaled PS, external f(k), EdS kernels)
# ==========================
p_EdS = deepcopy(base)
p_EdS["rescale_PS"] = False
p_EdS["use_beyond_eds_kernels"] = False  # <---- EdS kernels
tables_MG_EdS = pyfkpt.compute_multipoles(k=ks, pk=pk_lin, fk=fk_lin.copy(), f0=f0, **p_EdS)
k_MG_EdS, P0_MG_EdS, P2_MG_EdS, P4_MG_EdS = pyfkpt.rsd_multipoles(
    k=k_eval, nuis=nuis, z=z_pk, Om=Om, ap=False, tables=tables_MG_EdS
)

# ==========================
# (5) MG (unrescaled PS, external f(k), beyond-EdS kernels)
# ==========================
p_bEdS = deepcopy(base)
p_bEdS["rescale_PS"] = False
p_bEdS["use_beyond_eds_kernels"] = True  # <---- beyond EdS kernels
tables_MG_beyond = pyfkpt.compute_multipoles(k=ks, pk=pk_lin, fk=fk_lin.copy(), f0=f0, **p_bEdS)
k_MG_beyond, P0_MG_beyond, P2_MG_beyond, P4_MG_beyond = pyfkpt.rsd_multipoles(
    k=k_eval, nuis=nuis, z=z_pk, Om=Om, ap=False, tables=tables_MG_beyond
)