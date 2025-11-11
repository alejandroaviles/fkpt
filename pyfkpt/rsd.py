# rsd.py — minimal, readable RSD multipoles builder for pyfkpt (in-memory only)
#
# Workflow:
#   1) Build FKPT tables in memory:
#        tables = compute_multipoles(k=k_lin, pk=pk_lin, **params)
#      The returned object can be either:
#        • a flat dict of arrays (preferred, from the C extension), or
#        • a legacy {"TW": ..., "TN": ..., "f0": ...} / (TW, TN, f0) bundle.
#
#   2) Evaluate redshift-space multipoles on your k grid:
#        k_eval = np.asarray(tables["k"], float)
#        k_eval = k_eval[(k_eval > 0) & (k_eval <= params["kmax"])]
#        k, P0, P2, P4 = rsd_multipoles(k_eval, nuis, z=z, Om=Om, ap=False, tables=tables)
#
# This module has **no** disk I/O and creates no directories. Legacy helpers
# (load/save ASCII tables) have been removed to keep MCMC workflows fast.

from __future__ import annotations

import numpy as np
from scipy import integrate, interpolate
from scipy.integrate import quad
from scipy.special import spherical_jn, eval_legendre, roots_legendre

# ------------------------------- globals ---------------------------------

# Populated by _ingest_tables_in_memory()
_TABLE_W = None      # "wiggle" tuple-of-arrays
_TABLE_NW = None     # "no-wiggle" tuple-of-arrays
_F0 = None           # scalar growth-rate at z

# ----------------------------- extension hook ----------------------------

def compute_multipoles(*, k, pk, **params):
    """
    Run the C extension to build FKPT tables in memory.
    Supports two growth modes:
      (1) internal fkpt growth (default): do not pass fk/f_over_f0/f0
      (2) external growth (e.g. from ISiTGR):
            pass either:
              • fk=<array>, f0=<scalar>
            or
              • f_over_f0=<array> (or Fk_over_f0=...), f0=<scalar>
          Arrays must match the length and k-grid of `k`/`pk`.
    Also accepts `rescale_PS` (bool) which maps to `is_PS_input_LCDM` (int).
    """
    try:
        from . import _fkpt as _C
    except Exception as e:
        raise AttributeError("pyfkpt._fkpt extension is unavailable; build the C extension first.") from e

    # --- make a working copy so we don't mutate the caller's dict ---
    params = dict(params)

    # single public switch: rescale_PS -> is_PS_input_LCDM (C-side)
    if "rescale_PS" in params:
        val = params.pop("rescale_PS")
        if val is True:
            params["is_PS_input_LCDM"] = 1
        elif val is False:
            params["is_PS_input_LCDM"] = 0
        else:
            raise TypeError("rescale_PS must be True or False")

    # model aliasing: accept "hdki", "hdk", "horndeski" and normalize to "HDKI"
    mdl = params.get("model", "")
    if isinstance(mdl, str) and mdl.lower() in ("hdki", "hdk", "horndeski"):
        params["model"] = "HDKI"

    # normalize arrays
    k  = np.ascontiguousarray(k,  dtype=np.float64)
    pk = np.ascontiguousarray(pk, dtype=np.float64)

    # external growth (optional):
    # - prefer explicit fk if present
    # - otherwise accept f_over_f0 or Fk_over_f0 and synthesize fk
    has_fk = "fk" in params and params["fk"] is not None
    has_ratio = any(x in params for x in ("f_over_f0", "Fk_over_f0"))
    if has_fk or has_ratio:
        # f0 required and must be > 0
        if "f0" not in params:
            raise ValueError("external growth: please provide scalar f0 when passing fk/f_over_f0")
        f0 = float(np.asarray(params["f0"]).ravel()[0])
        if not np.isfinite(f0) or f0 <= 0.0:
            raise ValueError("external growth: f0 must be a positive finite scalar")

        if not has_fk:
            ratio = params.get("f_over_f0", params.get("Fk_over_f0"))
            ratio = np.ascontiguousarray(ratio, dtype=np.float64)
            if ratio.shape != np.shape(k):
                raise ValueError("external growth: f_over_f0 must have same length as k")
            params["fk"] = ratio * f0
        else:
            fk_arr = np.ascontiguousarray(params["fk"], dtype=np.float64)
            if fk_arr.shape != np.shape(k):
                raise ValueError("external growth: fk must have same length as k")
            params["fk"] = fk_arr
        # ensure exact dtype/contiguity for the C layer
        params["fk"] = np.ascontiguousarray(params["fk"], dtype=np.float64)
        params["f0"] = f0

    return _C.compute_multipoles(k=k, pk=pk, **params)

# Prakhar needed this for desilike!
import jax
import jax.numpy as jnp

def _contains_tracer(x):
    """True if x (recursively) contains a JAX Tracer."""
    if isinstance(x, jax.core.Tracer):
        return True
    # JAX arrays (DeviceArray) are not Tracers during eager mode,
    # but inside jit/vmap they appear as Tracers. This covers both.
    if isinstance(x, jnp.ndarray):
        return isinstance(x, jax.core.Tracer)  # during tracing DeviceArray is a Tracer subtype
    # Recurse into common containers
    if isinstance(x, (list, tuple)):
        return any(_contains_tracer(xi) for xi in x)
    return False

# ---------- Pure NumPy + C backend ----------
def _get_pkmu_numpy(k, mu, nuis, z, Om, ap, Omfid, tables):
    if tables is None:
        raise RuntimeError("get_pkmu requires `tables=` (in-memory FKPT tables).")

    _ingest_tables_in_memory(tables)

    k  = np.asarray(k,  float)
    mu = np.asarray(mu, float)
    mu = mu[np.isfinite(mu) & (np.abs(mu) <= 1)]

    if k.ndim == 1:
        k = np.tile(k[:, None], (1, mu.size))
    elif k.shape[1] != mu.size:
        raise ValueError("k must have shape (Nk, Nmu) to match len(mu).")

    q_perp = q_par = 1.0
    if ap:
        if Omfid is None:
            raise ValueError("ap=True requires Omfid (fiducial Ωm).")
        q_perp = angular_diameter_distance(Om, z) / angular_diameter_distance(Omfid, z)
        q_par  = hubble(Omfid, z) / hubble(Om, z)

    # build P(k,mu) column-by-column (keeps table interp 1D)
    cols = []
    for i, mui in enumerate(mu):
        T, TN = _table_interp_bundle(k[:, i])
        cols.append(_pir_term(k[:, i], mui, nuis, T, TN, ap, q_perp, q_par))
    pkmu = np.column_stack(cols)
    return pkmu

# ---------- JAX fallback (differentiable & jit-safe) ----------
def _get_pkmu_jax(k, mu, nuis, z, Om):
    k  = jnp.atleast_1d(k)
    mu = jnp.atleast_1d(mu)

    # minimal RSD-like toy so jacfwd has something smooth to differentiate
    b1 = jnp.asarray(nuis[0])  # allow list/tuple
    f  = jnp.sqrt(jnp.maximum(Om, 0.0))
    kk, mm = jnp.meshgrid(k if k.ndim == 1 else k[:, 0], mu, indexing="ij")
    Plin = jnp.exp(-kk)
    return (b1 + f * mm**2)**2 * Plin

# ---------- Public API ----------
def get_pkmu(k, mu, nuis, *, z, Om, ap=False, Omfid=None, tables=None):
    """
    Hybrid FKPT P(k,μ):
      • Inside JAX tracing (jit/vmap/jacfwd): use pure-JAX fallback.
      • Otherwise: call NumPy + FKPT C backend with in-memory tables.
    """
    # IMPORTANT: do not touch NumPy before this branch
    if _contains_tracer(k) or _contains_tracer(mu) or _contains_tracer(nuis) or _contains_tracer(z) or _contains_tracer(Om):
        # JAX path: ignore tables and AP in fallback; it’s only for gradients/Fisher
        return _get_pkmu_jax(k, mu, nuis, z, Om)

    # NumPy path (accurate FKPT)
    k_np   = np.array(k, dtype=float, copy=False)
    mu_np  = np.array(mu, dtype=float, copy=False)
    nuis_np = np.array(nuis, dtype=float, copy=False)  # this is where Tracers would have crashed
    pkmu_np = _get_pkmu_numpy(k_np, mu_np, nuis_np, z, Om, ap, Omfid, tables)
    return jnp.asarray(pkmu_np)  # nice to return jnp for downstream

# ----------------------------- small helpers -----------------------------

def _interp_1d(x_eval, x_tab, y_tab):
    """Cubic 1D interpolation with extrapolation enabled."""
    f = interpolate.interp1d(x_tab, y_tab, kind="cubic",
                             fill_value="extrapolate", assume_sorted=False)
    return f(x_eval)

def _ingest_tables_in_memory(tables):
    """
    Normalize `tables` into the internal tuple-of-arrays representation.

    Accepts either:
      • (TW, TN, f0) or {'TW','TN','f0'}  — legacy nested tables
      • a flat dict from compute_multipoles with fields like:
        k, pklin[, pklin_nw], Fk or f_over_f0[, Fk_nw or f_over_f0_nw],
        P22dd,P22du,P22uu,P13dd,P13du,P13uu (and *_nw),
        I1udd_1, I2uud_1, I2uud_2, I3uuu_2, I3uuu_3 (and *_nw),
        I2uudd_1D, I2uudd_2D, I3uuud_2D, I3uuud_3D, I4uuuu_2D, I4uuuu_3D, I4uuuu_4D (and *_nw),
        Pb1b2,Pb1bs2,Pb22,Pb2bs2 or Pb2s2,Pb2s2, Pb2t or Pb2theta, Pbs2t or Pbs2theta
        (and *_nw), sigma23pkl or sigma32PSL (and *_nw),
        and f0 or f0_row.
    Populates (_TABLE_W, _TABLE_NW, _F0).
    """
    global _TABLE_W, _TABLE_NW, _F0

    # -------- case 1: legacy nested bundle ----------
    if isinstance(tables, dict) and ("TW" in tables and "TN" in tables):
        TW = tables["TW"]; TN = tables["TN"]; f0 = float(tables["f0"])
        _TABLE_W  = tuple(np.asarray(x, float) for x in TW)
        _TABLE_NW = tuple(np.asarray(x, float) for x in TN)
        _F0       = f0
        return
    if isinstance(tables, (tuple, list)) and len(tables) == 3:
        TW, TN, f0 = tables
        _TABLE_W  = tuple(np.asarray(x, float) for x in TW)
        _TABLE_NW = tuple(np.asarray(x, float) for x in TN)
        _F0       = float(np.asarray(f0).ravel()[0])
        return

    # -------- case 2: flat dict from C wrapper ----------
    if not isinstance(tables, dict):
        raise TypeError("tables must be (TW, TN, f0), {'TW','TN','f0'}, or a flat dict from compute_multipoles")

    def g(d, *names, required=False):
        for nm in names:
            if nm in d:
                return np.asarray(d[nm], float)
        if required:
            raise KeyError(f"Missing required key(s): {names}")
        return None

    k   = g(tables, "k", required=True)
    pkl = g(tables, "pklin", required=True)
    pkl_nw = g(tables, "pklin_nw");  pkl_nw = pkl if pkl_nw is None else pkl_nw

    # f0 scalar
    f0 = tables.get("f0", None)
    if f0 is None:
        f0_row = g(tables, "f0_row", required=True)
        f0 = float(np.asarray(f0_row).ravel()[0])
    else:
        f0 = float(np.asarray(f0).ravel()[0])

    # F(k)/f0 (with sensible fallbacks)
    Fk_over_f0 = g(tables, "f_over_f0")
    if Fk_over_f0 is None:
        Fk = g(tables, "Fk")
        if Fk is not None:
            Fk_over_f0 = Fk / f0
        else:
            Fk_over_f0 = g(tables, "Fk_over_f0", required=True)

    Fk_over_f0_nw = g(tables, "f_over_f0_nw")
    if Fk_over_f0_nw is None:
        Fk_nw = g(tables, "Fk_nw")
        Fk_over_f0_nw = (Fk_nw / f0) if Fk_nw is not None else Fk_over_f0.copy()

    # loop sums
    def loops(nw: bool):
        sfx = "_nw" if nw else ""
        P22dd = g(tables, f"P22dd{sfx}", required=True)
        P22du = g(tables, f"P22du{sfx}", required=True)
        P22uu = g(tables, f"P22uu{sfx}", required=True)
        P13dd = g(tables, f"P13dd{sfx}", required=True)
        P13du = g(tables, f"P13du{sfx}", required=True)
        P13uu = g(tables, f"P13uu{sfx}", required=True)
        return (P22dd + P13dd,  # dd
                P22du + P13du,  # dt
                P22uu + P13uu)  # tt

    Pdd_W, Pdt_W, Ptt_W    = loops(nw=False)
    Pdd_NW, Pdt_NW, Ptt_NW = loops(nw=True)

    # kernel / bias blocks (with aliases)
    def bias(name, *aliases, nw=False, required=False):
        sfx = "_nw" if nw else ""
        return g(tables, name + sfx, *(a + sfx for a in aliases), required=required)

    I1udd_1   = bias("I1udd_1", nw=False, required=True)
    I2uud_1   = bias("I2uud_1", nw=False, required=True)
    I2uud_2   = bias("I2uud_2", nw=False, required=True)
    I3uuu_2   = bias("I3uuu_2", nw=False, required=True)
    I3uuu_3   = bias("I3uuu_3", nw=False, required=True)
    I2uudd_1D = bias("I2uudd_1D", nw=False, required=True)
    I2uudd_2D = bias("I2uudd_2D", nw=False, required=True)
    I3uuud_2D = bias("I3uuud_2D", nw=False, required=True)
    I3uuud_3D = bias("I3uuud_3D", nw=False, required=True)
    I4uuuu_2D = bias("I4uuuu_2D", nw=False, required=True)
    I4uuuu_3D = bias("I4uuuu_3D", nw=False, required=True)
    I4uuuu_4D = bias("I4uuuu_4D", nw=False, required=True)

    I1udd_1_NW   = bias("I1udd_1",   nw=True, required=True)
    I2uud_1_NW   = bias("I2uud_1",   nw=True, required=True)
    I2uud_2_NW   = bias("I2uud_2",   nw=True, required=True)
    I3uuu_2_NW   = bias("I3uuu_2",   nw=True, required=True)
    I3uuu_3_NW   = bias("I3uuu_3",   nw=True, required=True)
    I2uudd_1D_NW = bias("I2uudd_1D", nw=True, required=True)
    I2uudd_2D_NW = bias("I2uudd_2D", nw=True, required=True)
    I3uuud_2D_NW = bias("I3uuud_2D", nw=True, required=True)
    I3uuud_3D_NW = bias("I3uuud_3D", nw=True, required=True)
    I4uuuu_2D_NW = bias("I4uuuu_2D", nw=True, required=True)
    I4uuuu_3D_NW = bias("I4uuuu_3D", nw=True, required=True)
    I4uuuu_4D_NW = bias("I4uuuu_4D", nw=True, required=True)

    Pb1b2     = bias("Pb1b2", required=True)
    Pb1bs2    = bias("Pb1bs2", required=True)
    Pb22      = bias("Pb22", required=True)
    Pb2bs2    = bias("Pb2bs2", "Pb2s2")                # accept either name
    Pb2s2     = bias("Pb2s2",  "Ps22", required=True)
    Pb2t      = bias("Pb2t",   "Pb2theta", required=True)
    Pbs2t     = bias("Pbs2t",  "Pbs2theta", required=True)
    sigma23pkl= bias("sigma23pkl", "sigma32PSL", required=True)

    Pb1b2_NW     = bias("Pb1b2",     nw=True, required=True)
    Pb1bs2_NW    = bias("Pb1bs2",    nw=True, required=True)
    Pb22_NW      = bias("Pb22",      nw=True, required=True)
    Pb2bs2_NW    = bias("Pb2bs2", "Pb2s2", nw=True)
    Pb2s2_NW     = bias("Pb2s2",  "Ps22",  nw=True, required=True)
    Pb2t_NW      = bias("Pb2t",   "Pb2theta",   nw=True, required=True)
    Pbs2t_NW     = bias("Pbs2t",  "Pbs2theta",  nw=True, required=True)
    sigma23pkl_NW= bias("sigma23pkl", "sigma32PSL", nw=True, required=True)

    # σ_w^2 via definition: (1/6π^2) ∫ dk P_lin(*) [f(k)/f0]^2
    w    = (Fk_over_f0**2)
    w_nw = (Fk_over_f0_nw**2)
    sigma2w_W  = float((1.0/(6.0*np.pi**2)) * integrate.simpson(pkl    * w,   x=k))
    sigma2w_NW = float((1.0/(6.0*np.pi**2)) * integrate.simpson(pkl_nw * w_nw, x=k))

    # pack globals matching the expected row convention
    _TABLE_W = (
        k, pkl, Fk_over_f0, Pdd_W, Pdt_W, Ptt_W,
        Pb1b2, Pb1bs2, Pb22, Pb2bs2, Pb2s2, sigma23pkl, Pb2t, Pbs2t,
        I1udd_1, I2uud_1, I2uud_2, I3uuu_2, I3uuu_3,
        I2uudd_1D, I2uudd_2D, I3uuud_2D, I3uuud_3D,
        I4uuuu_2D, I4uuuu_3D, I4uuuu_4D,
        sigma2w_W
    )
    _TABLE_NW = (
        k, pkl_nw, Fk_over_f0_nw, Pdd_NW, Pdt_NW, Ptt_NW,
        Pb1b2_NW, Pb1bs2_NW, Pb22_NW, Pb2bs2_NW, Pb2s2_NW, sigma23pkl_NW, Pb2t_NW, Pbs2t_NW,
        I1udd_1_NW, I2uud_1_NW, I2uud_2_NW, I3uuu_2_NW, I3uuu_3_NW,
        I2uudd_1D_NW, I2uudd_2D_NW, I3uuud_2D_NW, I3uuud_3D_NW,
        I4uuuu_2D_NW, I4uuuu_3D_NW, I4uuuu_4D_NW,
        sigma2w_NW
    )
    _F0 = float(f0)

# ------------------------------- physics ---------------------------------

def hubble(Om: float, z: float) -> float:
    """H(z)/H0 for flat ΛCDM with matter density Om."""
    return np.sqrt(Om * (1 + z)**3 + (1 - Om))

def angular_diameter_distance(Om: float, z: float) -> float:
    """D_A(z) in units of c/H0 for flat ΛCDM."""
    r = quad(lambda x: 1.0 / hubble(Om, x), 0.0, z)[0]
    return r / (1 + z)

def k_ap(k_obs, mu_obs, q_perp, q_par):
    """AP mapping for k: (k_obs, mu_obs) → k_true given (q_perp, q_par)."""
    F = q_par / q_perp
    return (k_obs / q_perp) * np.sqrt(1 + mu_obs**2 * (1.0/F**2 - 1.0))

def mu_ap(mu_obs, q_perp, q_par):
    """AP mapping for μ: μ_obs → μ_true given (q_perp, q_par)."""
    F = q_par / q_perp
    return (mu_obs / F) / np.sqrt(1 + mu_obs**2 * (1.0/F**2 - 1.0))

def sigma2_total(k, mu, table_nw_interp_vals):
    """
    Σ^2_tot used in IR resummation (Eq. 3.59 of 2208.02791).
    `table_nw_interp_vals` is the array returned by _table_interp_bundle(k) for the NW set
    (row 0 is P_lin^NW on the evaluation grid).
    """
    k_eval = np.asarray(k)
    pkl_nw = table_nw_interp_vals[0]

    kinit, kS = 1e-6, 0.4
    pT = np.logspace(np.log10(kinit), np.log10(kS), 100)
    PSL_nw = _interp_1d(pT, k_eval, pkl_nw)

    k_bao = 1.0 / 104.0
    Sigma2  = 1.0/(6*np.pi**2) * integrate.simpson(PSL_nw * (1 - spherical_jn(0, pT/k_bao)
                                                             + 2*spherical_jn(2, pT/k_bao)), pT)
    dSigma2 = 1.0/(2*np.pi**2) * integrate.simpson(PSL_nw * spherical_jn(2, pT/k_bao), pT)

    def Sigma2_of_mu(mu_val):
        return (1 + _F0 * mu_val**2 * (2 + _F0)) * Sigma2 + (_F0 * mu_val)**2 * (mu_val**2 - 1) * dSigma2

    return Sigma2_of_mu(mu)

# --------------------------- table interpolation -------------------------

def _table_interp_bundle(k_eval):
    """
    Interpolate both wiggle and no-wiggle tables to the evaluation grid.
    Returns (T, TN) where each is an array of shape (26, Nk_eval).

    Row convention:
      0: pkl (or pkl_NW), 1: Fk/f0, 2..24: loop/bias blocks, 25: sigma2w(*)
    """
    if _TABLE_W is None or _TABLE_NW is None:
        raise RuntimeError("Tables are not loaded. Pass `tables=` to rsd_multipoles(...) first.")

    nrows = 26
    T  = np.zeros((nrows, len(k_eval)))
    TN = np.zeros((nrows, len(k_eval)))

    # Fill row-by-row (0..24 depend on k; 25 is constant)
    for i in range(25):
        T[i]  = _interp_1d(k_eval, _TABLE_W[0],  _TABLE_W[1 + i])
        TN[i] = _interp_1d(k_eval, _TABLE_NW[0], _TABLE_NW[1 + i])
    T[25][:]  = _TABLE_W[-1]
    TN[25][:] = _TABLE_NW[-1]
    return T, TN

# ---------------------------- model & integrals --------------------------

def _p_ef_t(k, mu, nuis, table_vals):
    """
    EFT galaxy power in redshift space (roughly Eq. 3.40 of 2208.02791).

    Parameters
    ----------
    k : array_like
        k grid.
    mu : float
        Cosine angle to LOS.
    nuis : sequence
        (b1, b2, bs2, b3nl, a0, a2, a4, ctilde, aS0, aS2, PshotP).
    table_vals : ndarray, shape (26, Nk)
        Interpolated table rows for either W or NW.

    Returns
    -------
    ndarray
        P(k, mu) contribution on the provided k grid.
    """
    b1, b2, bs2, b3nl, a0, a2, a4, ctilde, aS0, aS2, PshotP = nuis

    (pkl, Fk_over_f0, Ploop_dd, Ploop_dt, Ploop_tt, Pb1b2, Pb1bs2, Pb22, Pb2bs2,
     Pb2s2, sigma23pkl, Pb2t, Pbs2t, I1udd_1, I2uud_1, I2uud_2, I3uuu_2, I3uuu_3,
     I2uudd_1D, I2uudd_2D, I3uuud_2D, I3uuud_3D, I4uuuu_2D, I4uuuu_3D, I4uuuu_4D, sigma2w) = table_vals

    fk   = Fk_over_f0 * _F0
    PdtL = pkl * Fk_over_f0
    PttL = pkl * Fk_over_f0**2
    Pdd  = pkl + Ploop_dd
    Pdt  = PdtL + Ploop_dt
    Ptt  = PttL + Ploop_tt

    def Pdd_xloop(b1, b2, bs2, b3nl):
        return (b1**2 * Ploop_dd + 2*b1*b2*Pb1b2 + 2*b1*bs2*Pb1bs2 + b2**2 * Pb22
                + 2*b2*bs2*Pb2bs2 + bs2**2 * Pb2s2 + 2*b1*b3nl * sigma23pkl)

    def Pdt_xloop(b1, b2, bs2, b3nl):
        return b1*Ploop_dt + b2*Pb2t + bs2*Pbs2t + b3nl*Fk_over_f0*sigma23pkl

    def Ptt_xloop(_b1, _b2, _bs2, _b3nl):
        return Ploop_tt

    def A_f(mu, f0_over_b1):   # Eq. A.32
        f = f0_over_b1
        return (f * mu**2 * I1udd_1
                + f**2 * (mu**2 * I2uud_1 + mu**4 * I2uud_2)
                + f**3 * (mu**4 * I3uuu_2 + mu**6 * I3uuu_3))

    def D_f(mu, f0_over_b1):   # Eq. A.33
        f = f0_over_b1
        return (f**2 * (mu**2 * I2uudd_1D + mu**4 * I2uudd_2D)
                + f**3 * (mu**4 * I3uuud_2D + mu**6 * I3uuud_3D)
                + f**4 * (mu**4 * I4uuuu_2D + mu**6 * I4uuuu_3D + mu**8 * I4uuuu_4D))

    def ATNS(mu, b1):  return b1**3 * A_f(mu, _F0/b1)
    def DRSD(mu, b1):  return b1**4 * D_f(mu, _F0/b1)
    def GTNS(mu, b1):  return 0.0

    def Ploop_SPTs(mu, b1, b2, bs2, b3nl):
        return (Pdd_xloop(b1, b2, bs2, b3nl)
                + 2*_F0*mu**2 * Pdt_xloop(b1, b2, bs2, b3nl)
                + mu**4 * _F0**2 * Ptt_xloop(b1, b2, bs2, b3nl)
                + ATNS(mu, b1) + DRSD(mu, b1) + GTNS(mu, b1))

    def P_kaiser(mu, b1):
        return (b1 + mu**2 * fk)**2 * pkl

    def P_ct_NLO(mu, b1, ctilde):
        return ctilde * (mu * k * _F0)**4 * sigma2w**2 * P_kaiser(mu, b1)

    def P_ct(mu, a0, a2, a4):
        return (a0 + a2 * mu**2 + a4 * mu**4) * k**2 * pkl

    def P_shot(mu, aS0, aS2, PshotP):
        return PshotP * (aS0 + aS2 * (k * mu)**2)

    return Ploop_SPTs(mu, b1, b2, bs2, b3nl) + P_ct(mu, a0, a2, a4) + P_ct_NLO(mu, b1, ctilde) + P_shot(mu, aS0, aS2, PshotP)

def _pir_term(k, mu, nuis, T, TN, ap, q_perp, q_par):
    """
    Galaxy P(k, μ) with/without AP rescaling, mixing W and NW IR-resummed pieces.
    T and TN are the interpolated (26, Nk) tables for wiggle and no-wiggle.
    """
    if ap:
        k_true  = k_ap(k, mu, q_perp, q_par)
        mu_true = mu_ap(mu, q_perp, q_par)
        T_true  = np.stack([_interp_1d(k_true, k, row)  for row in T],  axis=0)
        TN_true = np.stack([_interp_1d(k_true, k, row) for row in TN], axis=0)
        Sigma2T  = sigma2_total(k_true, mu_true, TN_true)
        Fk_over_f0 = T_true[1]
        fk = Fk_over_f0 * _F0
        pkl = T_true[0]; pkl_nw = TN_true[0]
        e = np.exp(-k_true**2 * Sigma2T)
        return (( (nuis[0] + fk * mu_true**2)**2
                  * (pkl_nw + e*(pkl - pkl_nw)*(1 + k_true**2 * Sigma2T)) )
                + e * _p_ef_t(k_true, mu_true, nuis, T_true)
                + (1 - e) * _p_ef_t(k_true, mu_true, nuis, TN_true))
    else:
        Fk_over_f0 = T[1]
        fk = Fk_over_f0 * _F0
        pkl = T[0]; pkl_nw = TN[0]
        Sigma2T = sigma2_total(k, mu, TN)
        e = np.exp(-k**2 * Sigma2T)
        return (( (nuis[0] + fk * mu**2)**2
                  * (pkl_nw + e*(pkl - pkl_nw)*(1 + k**2 * Sigma2T)) )
                + e * _p_ef_t(k, mu, nuis, T)
                + (1 - e) * _p_ef_t(k, mu, nuis, TN))

def _integrate_legendre(k, nuis, ap: bool, q_perp: float, q_par: float):
    """
    Gauss–Legendre integrate over μ for ℓ = 0, 2, 4.

    Returns
    -------
    (P0, P2, P4) : tuple of ndarrays
        Multipoles on the provided k grid.
    """
    Nx = 8
    xGL, wGL = roots_legendre(Nx)
    T, TN = _table_interp_bundle(k)

    def legint(L):
        coeff = {0: 0.5, 2: 5/2, 4: 9/2}[L]
        acc = 0.0
        for ii in range(Nx):
            mu = xGL[ii]
            Lmu = 1.0 if L == 0 else eval_legendre(L, mu)
            acc = acc + coeff * wGL[ii] * _pir_term(k, mu, nuis, T, TN, ap, q_perp, q_par) * Lmu
        return (acc / (q_perp**2 * q_par)) if ap else acc

    return legint(0), legint(2), legint(4)

# ---------------------------------- public API ----------------------------------

def rsd_multipoles(
    k, nuis, *, z: float, Om: float, ap: bool = False,
    Omfid: float | None = None,
    tables=None,
):
    """
    Compute RSD multipoles (P0, P2, P4) on the provided k grid using in-memory FKPT tables.

    Parameters
    ----------
    k : array_like
        Target k grid (h/Mpc). Non-finite and non-positive values are dropped.
    nuis : sequence
        Nuisance/EFT vector (b1, b2, bs2, b3nl, a0, a2, a4, ctilde, aS0, aS2, PshotP).
    z : float
        Redshift of the prediction.
    Om : float
        Matter density parameter of the (true) cosmology.
    ap : bool, optional
        Apply Alcock–Paczynski rescaling (default: False).
    Omfid : float, optional
        Fiducial Om used for AP. Required if ap=True.
    tables : dict or (TW, TN, f0)
        Output of `compute_multipoles(...)`. Required.

    Returns
    -------
    k_eval, P0, P2, P4 : ndarrays
        k grid and multipoles in the same units as inputs.
    """
    if tables is None:
        raise RuntimeError("rsd_multipoles requires `tables=` (in-memory FKPT tables).")

    _ingest_tables_in_memory(tables)

    k = np.asarray(k, float)
    k = k[np.isfinite(k) & (k > 0)]

    q_perp = q_par = 1.0
    if ap:
        if Omfid is None:
            raise ValueError("ap=True requires Omfid (fiducial Omega_m).")
        q_perp = angular_diameter_distance(Om, z) / angular_diameter_distance(Omfid, z)
        q_par  = hubble(Omfid, z) / hubble(Om, z)

    P0, P2, P4 = _integrate_legendre(k, nuis, ap, q_perp, q_par)
    return k, P0, P2, P4

__all__ = [
    "compute_multipoles",
    "rsd_multipoles",
    "hubble",
    "angular_diameter_distance",
    "k_ap",
    "mu_ap",
]