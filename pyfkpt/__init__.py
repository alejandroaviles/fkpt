# pyfkpt/__init__.py
"""
pyfkpt: minimal in-memory interface.

Public API:
- compute_multipoles(k, pk, **params)  -> dict of FKPT tables (from C extension)
- rsd_multipoles(k, nuis, z=..., Om=..., ap=False, Omfid=None, tables=...)
- cosmology helpers: hubble, angular_diameter_distance, k_ap, mu_ap

Note: This package intentionally avoids creating or reading any on-disk tables.
"""

from .rsd import (
    compute_tables,
    rsd_multipoles,
    hubble,
    angular_diameter_distance,
    k_ap,
    mu_ap,
)

__all__ = [
    "compute_tables",
    "rsd_multipoles",
    "hubble",
    "angular_diameter_distance",
    "k_ap",
    "mu_ap",
]