#!/usr/bin/env python3
from __future__ import annotations
import os
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pyfkpt.rsd as pyfkpt  # legacy-style API lives here

# -------------------- optional: quick table dump --------------------
def get_matter_power_spectrum_from_file(path, z):
    data1 = np.loadtxt(path+'ps_k.txt')
    data2 = np.loadtxt(path+'ps.txt')
    k = data1             # h/Mpc
    pk1d = data2              # (Mpc/h)^3
    # mimic CAMB's (k, z_array, pk[z_index, k_index]) signature
    return k, np.array([z]), pk1d[np.newaxis, :]
    
def dump_pyfkpt_tables(out_dir: str, proc: int, z: float, suffix: str, Om_for_f0: float) -> None:
    TW, TN, f0 = pyfkpt.load_fkpt_outputs(
        path_out=out_dir, proc=proc, z=z, suffix=suffix, Om_for_f0=Om_for_f0
    )
    labels = [
        "k","pklin","Fk_over_f0",
        "Ploop_dd","Ploop_dt","Ploop_tt",
        "Pb1b2","Pb1bs2","Pb22","Pb2bs2","Pb2s2",
        "sigma23pkl","Pb2t","Pbs2t",
        "I1udd_1","I2uud_1","I2uud_2","I3uuu_2","I3uuu_3",
        "I2uudd_1D","I2uudd_2D","I3uuud_2D","I3uuud_3D",
        "I4uuuu_2D","I4uuuu_3D","I4uuuu_4D",
        "sigma2w (const)"
    ]
    def brief(tag, arr):
        a = np.asarray(arr)
        if a.ndim == 0:
            print(f"{tag:>14s}: {a: .6e}")
        else:
            print(f"{tag:>14s}: shape={a.shape} min={a.min(): .6e} max={a.max(): .6e} first5={a[:5]}")
    print("\n=== pyfkpt table snapshot ===")
    print(f"f0 = {f0:.6f}")
    for tag, col in zip(labels, TW):  brief("W."+tag, col)
    for tag, col in zip(labels, TN):  brief("NW."+tag, col)

# -------------------- user params --------------------
#bias parameters
b1 = 1.70               
b2 = -0.45
bs2 = -4/7*(b1 - 1)
b3nl = 32/315*(b1 - 1)

#EFT parameters
alpha0 = 3.0               #units: [Mpc/h]^2              
alpha2 = -29.0             #units: [Mpc/h]^2
alpha4 = 0.0               #units: [Mpc/h]^2
ctilde = 0.0               #units: [Mpc/h]^4

#Stochatics parameters
alpha0shot = 0.08
alpha2shot = -8.0          #units: [Mpc/h]^2
pshotp = 5000              #units: [Mpc/h]^3

# cosmology (matches your earlier test)
h      = 0.6711
ombh2  = 0.022
omch2  = 0.122
omnuh2 = 0.0006442
As     = 2e-9
ns     = 0.965
z_pk   = 0.5
N_ur   = 1 + 2.0328
Om     = (ombh2 + omch2 + omnuh2) / (h**2)

_SUFFIX = "pyfkpt"   # consistent file tag

def main():
    root = Path(__file__).resolve().parents[1] if "__file__" in globals() else Path.cwd()
    (root / "Input").mkdir(parents=True, exist_ok=True)
    (root / "Output").mkdir(parents=True, exist_ok=True)

    # ---------- 1) CAMB: total matter P(k) at z = z_pk ----------
    import camb
    from camb import model
    mnu_eV = 93.14 * omnuh2  # 1 massive nu with sum mnu = 93.14 * Ων h^2
    pars = camb.set_params(
        H0=100*h, ombh2=ombh2, omch2=omch2,
        mnu=mnu_eV, nnu=N_ur, num_massive_neutrinos=1, As=As, ns=ns,
    )
    pars.set_matter_power(redshifts=[z_pk], kmax=2.05)
    pars.NonLinear = model.NonLinear_none
    results = camb.get_results(pars)
    k_camb, _, pk_camb = results.get_matter_power_spectrum(
        minkh=1.0e-4, maxkh=2.0, npoints=1000
    )
    pk_camb = pk_camb[0]  # (Mpc/h)^3 at z=z_pk

    # ---------- 1b) Optional overlay from tests/ (CLASS file) ----------
    tests_dir = root / "tests"
    pk_path_k  = tests_dir / "ps_k.txt"
    pk_path_pk = tests_dir / "ps.txt"
    have_class = pk_path_k.exists() and pk_path_pk.exists()
    if have_class:
        k_cls  = np.loadtxt(pk_path_k)
        pk_cls = np.loadtxt(pk_path_pk)

    # quick overlay/ratio figure (CAMB vs CLASS-file if present)
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, height_ratios=[3,1], figsize=(7,6))
    ax1.loglog(k_camb, pk_camb, lw=2, label='CAMB')
    if have_class:
        ax1.loglog(k_cls, pk_cls, lw=2, ls='--', label='CLASS file')
    ax1.set_ylabel(r'$P_{\rm m}(k)\;[(\mathrm{Mpc}/h)^3]$')
    ax1.set_xlim(1e-4, 2.0); ax1.grid(True, which='both', alpha=0.3); ax1.legend()

    if have_class:
        kmin = max(k_camb.min(), k_cls.min())
        kmax = min(k_camb.max(), k_cls.max())
        m = (k_camb >= kmin) & (k_camb <= kmax)
        cls_on_camb = np.interp(k_camb[m], k_cls, pk_cls, left=np.nan, right=np.nan)
        ratio = pk_camb[m] / cls_on_camb
        ax2.semilogx(k_camb[m], ratio, lw=1.5)
        ax2.set_ylabel('CAMB / CLASS')
    else:
        ax2.text(0.02, 0.5, "No tests/ CLASS file\nfor ratio", transform=ax2.transAxes)
    ax2.set_xlabel(r'$k\,[h\,\mathrm{Mpc}^{-1}]$'); ax2.grid(True, which='both', alpha=0.3)
    fig.tight_layout(); fig.savefig(root / "plot.png", dpi=150)

    # Working (k, pk) for FKPT:
    k, pk = k_camb, pk_camb
    k, _, pk = get_matter_power_spectrum_from_file('/n/home12/cgarciaquintero/MG/codes/pyfkpt/tests/', z_pk)
    pk = pk[0]  # keep downstream code unchanged

    # ensure the C side also sees an Input/pyfkpt_ps_<pid> file (FKPT convention)
    fname = f"pyfkpt_ps_{os.getpid()}"
    np.savetxt(root / "Input" / fname, np.column_stack([k, pk]))

    # ---------- 2) compute FKPT tables in-memory (in-process) ----------
    params = dict(
        z=z_pk, Om=Om, h=h,
        b1=b1, b2=b2, bs2=bs2, b3nl=b3nl,
        alpha0=alpha0, alpha2=alpha2, alpha4=alpha4, ctilde=ctilde,
        PshotP=pshotp, alpha0shot=alpha0shot, alpha2shot=alpha2shot,
        kmin=float(max(1.0e-3, k.min())), kmax=min(float(k.max()), 0.5), Nk=min(k.size, 240),
        nquadSteps=300, chatty=0, model="LCDM", fR0=1e-15,
    )
    out = pyfkpt.compute_multipoles(k=k, pk=pk, **params)
    print("compute_multipoles keys:", sorted(out.keys()))

    # ---------- 3) write full tables, reload, (optional) dump ----------
    proc, suffix = os.getpid(), _SUFFIX
    out_dir = str(root / "Output")
    pyfkpt.store_tables_full(out, out=out_dir, proc=proc, z=z_pk, suffix=suffix)
    TW, TN, f0 = pyfkpt.load_fkpt_outputs(
        path_out=out_dir, proc=proc, z=z_pk, suffix=suffix, Om_for_f0=Om
    )
    print("loaded tables: len(W), len(NW), f0 =", len(TW), len(TN), f0)
    # dump_pyfkpt_tables(out_dir, proc, z_pk, suffix, Om_for_f0=Om)  # uncomment if you want the snapshot

    # ---------- 4) assemble multipoles (pyfkpt) ----------
    nuis = [b1,b2,bs2,b3nl,alpha0,alpha2,alpha4,ctilde,alpha0shot,alpha2shot,pshotp]
    k_eval = np.asarray(out["k"], float)
    k_eval = k_eval[(k_eval > 0) & (k_eval <= params["kmax"])]
    kgrid, P0, P2, P4 = pyfkpt.rsd_multipoles(
        k=k_eval, nuis=nuis, z=z_pk, Om=Om, ap=False,
        path_out=out_dir, proc=proc, suffix=suffix
    )
    print("kgrid size:", kgrid.size, "P0[:3]:", np.asarray(P0[:3]))

    # ---------- 5) OLD FKPT overlay & metrics (if tests/fkpt_*.txt exist) ----------
    fk_k, fk_P0, fk_P2, fk_P4 = None, None, None, None
    try:
        fk_k  = np.loadtxt(tests_dir / "fkpt_k.txt")
        fk_P0 = np.loadtxt(tests_dir / "fkpt_P0.txt")
        fk_P2 = np.loadtxt(tests_dir / "fkpt_P2.txt")
        fk_P4 = np.loadtxt(tests_dir / "fkpt_P4.txt")
        p4_path = tests_dir / "fkpt_P4.txt"
        if p4_path.exists():
            fk_P4 = np.loadtxt(p4_path)
    except Exception as e:
        print("No old-FKPT reference files found or failed to load:", e)

    def compare_and_print(new_k, new_P, old_k, old_P, ell):
        if old_k is None or old_P is None:
            return None
        # restrict to overlapping k-range
        lo = max(new_k.min(), old_k.min())
        hi = min(new_k.max(), old_k.max())
        m = (new_k >= lo) & (new_k <= hi)
        if not np.any(m):
            print(f"[ℓ={ell}] no overlap in k-range")
            return None
        old_on_new = np.interp(new_k[m], old_k, old_P, left=np.nan, right=np.nan)
        mask = np.isfinite(old_on_new) & np.isfinite(new_P[m])
        if not np.any(mask):
            print(f"[ℓ={ell}] no finite overlap")
            return None
        diff = new_P[m][mask] - old_on_new[mask]
        rel  = np.abs(diff) / np.maximum(1e-12, np.abs(old_on_new[mask]))
        print(f"[ℓ={ell}]  max |Δ| = {np.max(np.abs(diff)):.3e}   "
              f"max rel = {np.max(rel):.3e}   rms rel = {np.sqrt(np.mean(rel**2)):.3e}   "
              f"n={mask.sum()}")
        return (new_k[m][mask], old_on_new[mask])

    _ = compare_and_print(kgrid, P0, fk_k, fk_P0, 0)
    _ = compare_and_print(kgrid, P2, fk_k, fk_P2, 2)
    if fk_P4 is not None:
        _ = compare_and_print(kgrid, P4, fk_k, fk_P4, 4)

    # ---------- 6) plot k P_ell with old-FKPT (dashed) ----------
    fig, ax = plt.subplots(figsize=(7.8, 5.0))
    ax.plot(kgrid, kgrid*P0, label=r"pyfkpt $\ell=0$", lw=2.2)
    ax.plot(kgrid, kgrid*P2, label=r"pyfkpt $\ell=2$", lw=2.2)
    ax.plot(kgrid, kgrid*P4, label=r"pyfkpt $\ell=4$", lw=1.8)

    if fk_k is not None and fk_P0 is not None:
        ax.plot(fk_k, fk_k*fk_P0, "--", lw=1.7, label=r"fkpt (old) $\ell=0$")
    if fk_k is not None and fk_P2 is not None:
        ax.plot(fk_k, fk_k*fk_P2, "--", lw=1.7, label=r"fkpt (old) $\ell=2$")
    if fk_k is not None and fk_P4 is not None:
        ax.plot(fk_k, fk_k*fk_P4, "--", lw=1.7, label=r"fkpt (old) $\ell=4$")

    xmax = float(np.nanmax(kgrid if fk_k is None else np.concatenate([kgrid, fk_k])))
    ax.set_xlim(0.0, xmax)
    ax.set_xlabel(r"$k\,[h\,\mathrm{Mpc}^{-1}]$")
    ax.set_ylabel(r"$k\,P_\ell(k)\,[h^{-1}\,\mathrm{Mpc}]^2$")
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()
    out_png = root / "pyfkpt_smoke_multipoles.png"
    fig.savefig(out_png, dpi=150)
    print("Saved:", out_png)

    # (optional) LCDM sanity: Fk_over_f0 ~ 1
    try:
        T, _ = pyfkpt._table_interp_bundle(kgrid)  # private helper; OK for smoke
        print("mean Fk_over_f0 =", float(np.nanmean(T[1])))
    except Exception:
        pass

if __name__ == "__main__":
    main()