from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext
from pathlib import Path
import os

ROOT    = Path(__file__).resolve().parent
PKG_DIR = ROOT / "pyfkpt"
SRC     = ROOT / "src"
CSRC    = ROOT / "csrc"
LIBS    = ROOT / "libs"

# --- exclude lists as Path objects ---
EXCLUDE_FILES = {
    LIBS / "inout.c",
    LIBS / "getparam.c",
    SRC  / "main.c",                # compile main.c only via the shim
}
EXCLUDE_DIRS = {
    SRC / "delete",
}

def rel(p: Path) -> str:
    """Return a POSIX-style path relative to ROOT (what setuptools likes)."""
    return p.relative_to(ROOT).as_posix()

def list_c(d: Path):
    if not d.exists():
        return []
    excluded_files = {f.resolve() for f in EXCLUDE_FILES}
    excluded_dirs  = [p.resolve() for p in EXCLUDE_DIRS]

    out = []
    for p in d.rglob("*.c"):
        pr = p.resolve()
        # skip files in excluded directories
        if any(pr.is_relative_to(ed) for ed in excluded_dirs):
            continue
        # skip explicitly excluded files
        if pr in excluded_files:
            continue
        out.append(rel(pr))
    return out

# ---- collect sources ----
sources = []

# Python extension glue (module)
mod_c = PKG_DIR / "_fkptmodule.c"
if mod_c.exists():
    sources.append(rel(mod_c.resolve()))

# the shim that provides fkpt_main
shim_c = CSRC / "pyext_link_shims.c"
if not shim_c.exists():
    raise RuntimeError("Missing csrc/pyext_link_shims.c (the fkpt_main shim).")
sources.append(rel(shim_c.resolve()))

# project C sources (minus excluded)
sources += list_c(SRC) + list_c(CSRC) + list_c(LIBS)

# de-duplicate while preserving order
seen = set()
sources = [s for s in sources if not (s in seen or seen.add(s))]

if not sources:
    raise RuntimeError("No C sources found. Check tree layout and setup.py.")

# include dirs (NumPy added at build time)
include_dirs = [rel(p.resolve()) for p in (PKG_DIR, SRC, CSRC, LIBS) if p.exists()]

extra_compile_args = ["-O3", "-DNDEBUG", "-std=c99"]
libraries = ["m"] if os.name != "nt" else []
extra_link_args = []  # optionally: ["-Wl,-z,defs"] to catch unresolved symbols at link time

class build_ext(_build_ext):
    def finalize_options(self):
        super().finalize_options()
        try:
            import numpy as _np
            self.include_dirs.append(_np.get_include())
        except Exception:
            pass

ext = Extension(
    name="pyfkpt._fkpt",
    sources=sources,
    include_dirs=include_dirs,
    extra_compile_args=extra_compile_args,
    define_macros=[
        ('FKPT_PYEXT','1'),
        ('_GNU_SOURCE','1'),
        ('STARTRUN_DEBUG','0'),   # enable very chatty C logs
        ('PYFKPT_DEBUG','0'),     # enable py wrapper logs
    ],
    libraries=libraries,
    extra_link_args=extra_link_args,
    language="c",
)

setup(
    name="pyfkpt",
    version="0.1.0",
    packages=["pyfkpt"],
    package_dir={"": "."},
    ext_modules=[ext],
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
)