# pyfkpt/pyfkpt/test.py  (minimal smoke test)

from __future__ import annotations
import os, numpy as np
from pathlib import Path
import pyfkpt  # compute_multipoles, store_tables, load_fkpt_outputs, rsd_multipoles

# --- quick debug: print what pyfkpt loads from the legacy tables ---
def dump_pyfkpt_tables(out_dir, proc, z, suffix, Om_for_f0):
    import os, numpy as np, pyfkpt

    # load (this also fills the module globals used during integration)
    TW, TN, f0 = pyfkpt.load_fkpt_outputs(
        path_out=out_dir, proc=proc, z=z, suffix=suffix, Om_for_f0=Om_for_f0
    )
    labels = [
        "k", "pklin", "Fk_over_f0",
        "Ploop_dd", "Ploop_dt", "Ploop_tt",
        "Pb1b2", "Pb1bs2", "Pb22", "Pb2bs2", "Pb2s2",
        "sigma23pkl", "Pb2t", "Pbs2t",
        "I1udd_1", "I2uud_1", "I2uud_2", "I3uuu_2", "I3uuu_3",
        "I2uudd_1D", "I2uudd_2D", "I3uuud_2D", "I3uuud_3D",
        "I4uuuu_2D", "I4uuuu_3D", "I4uuuu_4D",
        "sigma2w (const)"
    ]

   # def brief(tag, arr):
   #     a = np.asarray(arr)
   #     if np.ndim(a) == 0:
   #         print(f"{tag:>14s}: {a: .6e}")
   #     else:
   #         print(f"{tag:>14s}: shape={a.shape}   min={a.min(): .6e}   max={a.max(): .6e}   first5={a[:5]}")
#
   # print("\n=== pyfkpt table snapshot ===")
   # print(f"f0 = {f0:.6f}")
   # for tag, col in zip(labels, TW):
   #     brief("W."+tag, col)
   # for tag, col in zip(labels, TN):
   #     brief("NW."+tag, col)

    # --- raw-file cross checks (are the sums & columns consistent?) ---
    KF   = np.loadtxt(os.path.join(out_dir, f"kfunctions_{proc}_z{z}.dat"))
    KFNW = np.loadtxt(os.path.join(out_dir, f"kfunctions_nw_{proc}_z{z}.dat"))
    LIN  = np.loadtxt(os.path.join(out_dir, f"linear_{proc}_z{z}.dat"))

    # columns by construction in the legacy files:
    # 0:k, 1:P22dd, 2:P22du, 3:P22uu, 4:P13dd, 5:P13du, 6:P13uu,
    # 7..18: I*-blocks, 19..26: bias blocks, 27:pklin, 28:Fk, 29:f0

    def relerr(name, a, b):
        denom = np.maximum(1.0, np.max(np.abs(b)))
        print(f"max rel.err {name:>10s}: {np.max(np.abs(a-b))/denom: .3e}")

    # check that the loop sums match what rsd.py builds
    relerr("Ploop_dd", TW[3], KF[1] + KF[4])
    relerr("Ploop_dt", TW[4], KF[2] + KF[5])
    relerr("Ploop_tt", TW[5], KF[3] + KF[6])

    # check bias block mapping
    relerr("Pb1b2",   TW[6],  KF[19]);  relerr("Pb1bs2",  TW[7],  KF[20])
    relerr("Pb22",    TW[8],  KF[21]);  relerr("Pb2bs2",  TW[9],  KF[22])
    relerr("Pb2s2",   TW[10], KF[23]);  relerr("Pb2t",     TW[12], KF[24])
    relerr("Pbs2t",   TW[13], KF[25]);  relerr("sigma23",  TW[11], KF[26])

    # recompute sigma_w^2 from linear like use_fkpt.py / your legacy code
    k_lin, pklin, Fk = LIN[0], LIN[1], LIN[2]
    f0_from_file = KF[29,0]
    w = (Fk / f0_from_file)**2
    sigW  = (1.0/(6.0*np.pi**2))*np.trapz(pklin * w, x=k_lin)
    sigNW = (1.0/(6.0*np.pi**2))*np.trapz(np.loadtxt(os.path.join(out_dir, f"linear_{proc}_z{z}.dat"))[4] * w, x=k_lin)
    print(f"sigma2w recomputed:   W={sigW:.6e}   NW={sigNW:.6e}")
    print(f"sigma2w in tables:    W={TW[-1]:.6e}   NW={TN[-1]:.6e}")

def get_matter_power_spectrum_from_file(path, z):
    data1 = np.loadtxt(path+'ps_k.txt')
    data2 = np.loadtxt(path+'ps.txt')
    k = data1             # h/Mpc
    pk1d = data2              # (Mpc/h)^3
    # mimic CAMB's (k, z_array, pk[z_index, k_index]) signature
    return k, np.array([z]), pk1d[np.newaxis, :]

# -------------------- parameters (as requested) --------------------
#bias parameters
b1 = 1
b2 = 0
bs2 = 0 #0-4/7*(b1 - 1)
b3nl = 0 #32/315*(b1 - 1)

#EFT parameters
alpha0 = 0.0               #units: [Mpc/h]^2              
alpha2 = 0.0             #units: [Mpc/h]^2
alpha4 = 0.0               #units: [Mpc/h]^2
ctilde = 0.0               #units: [Mpc/h]^4

#Stochatics parameters
alpha0shot = 0.00
alpha2shot = 0.0          #units: [Mpc/h]^2
pshotp = 5000              #units: [Mpc/h]^3

# cosmology
h      = 0.6711
ombh2  = 0.022
omch2  = 0.122
omnuh2 = 0.0006442
As     = 2e-9
ns     = 0.965
z_pk   = 0.5
N_ur   = 1+2.0328

Om = (ombh2 + omch2 + omnuh2) / (h**2)
print("Omega_m =", Om)

# consistent file-tag everywhere
_SUFFIX = "pyfkpt"

# -------------------- main --------------------
def main():
    # repo paths
    root = Path(__file__).resolve().parent.parent
    (root / "Input").mkdir(parents=True, exist_ok=True)
    (root / "Output").mkdir(parents=True, exist_ok=True)

    # ---- CAMB: total matter P(k) at z=z_pk, matching the numbers above ----
    import camb
    from camb import model
    
    # one massive nu with Σmν = 93.14 * Ων h^2; N_ur massless d.o.f.
    mnu_eV = 93.14 * omnuh2
    pars = camb.set_params(
        H0=100*h, ombh2=ombh2, omch2=omch2,
        mnu=mnu_eV, nnu=N_ur, num_massive_neutrinos=1, As=As, ns=ns,
)
    pars.set_matter_power(redshifts=[z_pk], kmax=2.05)
    pars.NonLinear = model.NonLinear_none
    results = camb.get_results(pars)

    k_camb, _, pk_camb = results.get_matter_power_spectrum(minkh=1.0e-4, maxkh=2.0, npoints=1000)
    pk_camb = pk_camb[0]  # total matter P(k) [ (Mpc/h)^3 ]

    #k, _, pk = get_matter_power_spectrum_from_file('/n/home12/cgarciaquintero/MG/codes/pyfkpt/tests/', z_pk)
    #pk = pk[0]  # keep downstream code unchanged
    
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.interpolate import interp1d

    # --- your inputs (already computed) ---
    # k_camb, _, pk_camb
    # k, _, pk

    # 1) Overlay the two P(k)
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, height_ratios=[3, 1], figsize=(7, 6))

    ax1.loglog(k_camb, pk_camb, lw=2, label='CAMB')
    ax1.loglog(k,       pk,       lw=2, ls='--', label='CLASS file')
    ax1.set_ylabel(r'$P_{\rm m}(k)\;[(\mathrm{Mpc}/h)^3]$')
    ax1.set_xlim(1e-4, 2.0)
    ax1.grid(True, which='both', alpha=0.3)
    ax1.legend()

    # 2) Ratio CAMB / CLASS on the common k-range
    kmin = max(k_camb.min(), k.min())
    kmax = min(k_camb.max(), k.max())
    m = (k_camb >= kmin) & (k_camb <= kmax)

    interp_class = interp1d(k, pk, kind='linear', bounds_error=False, fill_value=np.nan)
    ratio = pk_camb[m] / interp_class(k_camb[m])

    ax2.semilogx(k_camb[m], ratio, lw=1.5)
    ax2.set_xlabel(r'$k\,[h\,\mathrm{Mpc}^{-1}]$')
    ax2.set_ylabel('CAMB / CLASS')
    ax2.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    fig.savefig('plot.png', dpi=150)

    # ensure the C side also sees an Input/pyfkpt_ps_<pid> file (FKPT’s convention)
    fname = f"pyfkpt_ps_{os.getpid()}"
    np.savetxt(root / "Input" / fname, np.column_stack([k, pk]))

    # ---- run pyfkpt fast path (C extension), then write & reload legacy tables ----
    params = dict(
        z=z_pk, Om=Om, h=h,
        b1=b1, b2=b2, bs2=bs2, b3nl=b3nl,
        alpha0=alpha0, alpha2=alpha2, alpha4=alpha4, ctilde=ctilde,
        PshotP=pshotp, alpha0shot=alpha0shot, alpha2shot=alpha2shot,
        kmin=1.0e-3, kmax=0.5, Nk=120,
        nquadSteps=300, chatty=0, model="LCDM", fR0=1e-15,
    )

    out = pyfkpt.compute_multipoles(k=k, pk=pk, **params)
    print(sorted(out.keys()))
    # write full 30-col tables that loader expects
    pyfkpt.store_tables_full(out, out=str(root / "Output"), proc=os.getpid(), z=z_pk, suffix=_SUFFIX)
    dump_pyfkpt_tables(out_dir=str(root / "Output"), proc=os.getpid(), z=z_pk, suffix=_SUFFIX, Om_for_f0=Om)

    # ---- assemble multipoles ----
    nuis = [b1, b2, bs2, b3nl, alpha0, alpha2, alpha4, ctilde, alpha0shot, alpha2shot, pshotp]
    k_eval = np.asarray(out["k"], float)
    k_eval = k_eval[(k_eval > 0) & (k_eval <= params["kmax"])]
    kgrid, P0, P2, P4 = pyfkpt.rsd_multipoles(
        k=k_eval, nuis=nuis, z=z_pk, Om=Om, ap=False,
        path_out=str(root / "Output"), proc=os.getpid(), suffix=_SUFFIX
    )

    # For LCDM (scale-independent growth), this should be ~1 everywhere:
    T, _ = pyfkpt.rsd._table_interp_bundle(k_eval)   # or however you reach it

    # ---- optional fkpt overlay from tests/ ----
    tests_dir = root / "tests"
    fk_k, kP0_fk, kP2_fk = None, None, None
    try:
        fk_k  = np.loadtxt(tests_dir / "fkpt_k.txt")
        fk_p0 = np.loadtxt(tests_dir / "fkpt_P0.txt")
        fk_p2 = np.loadtxt(tests_dir / "fkpt_P2.txt")
        kP0_fk, kP2_fk = fk_k * fk_p0, fk_k * fk_p2
    except Exception:
        pass

    # ---- plot ----
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.8, 5.0))
    ax.plot(kgrid, kgrid*P0, label=r"pyfkpt $\ell=0$", lw=2.2)
    ax.plot(kgrid, kgrid*P2, label=r"pyfkpt $\ell=2$", lw=2.2)

    if fk_k is not None and kP0_fk is not None:
        ax.plot(fk_k, kP0_fk, "--", lw=1.8, label=r"fkpt $\ell=0$")
        ax.plot(fk_k, kP2_fk, "--", lw=1.8, label=r"fkpt $\ell=2$")

    xmax = float(np.nanmax(kgrid if fk_k is None else np.concatenate([kgrid, fk_k])))
    ax.set_xlim(0.0, xmax)
    ax.set_xlabel(r"$k\,[h\,\mathrm{Mpc}^{-1}]$")
    ax.set_ylabel(r"$k\,P_\ell(k)\,[h^{-1}\,\mathrm{Mpc}]^2$")
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()
    out_png = root / "pyfkpt_smoke_multipoles.png"
    fig.savefig(out_png, dpi=150)
    print("Saved:", out_png)

if __name__ == "__main__":
    main()