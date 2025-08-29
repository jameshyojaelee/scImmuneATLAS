"""Single-cell Immune Atlas package.

Environment hardening:
- Disable user site-packages to avoid NumPy/Matplotlib version conflicts.
"""

# Hard-disable user site to prevent mixing with ~/.local packages
import os as _os
import sys as _sys
try:
    import site as _site  # type: ignore
    _usr = _site.getusersitepackages() if hasattr(_site, "getusersitepackages") else None
    if _usr and _usr in _sys.path:
        _sys.path.remove(_usr)
except Exception:
    pass
_os.environ.setdefault("PYTHONNOUSERSITE", "1")

__version__ = "0.1.0"
__author__ = "Single-cell Immune Atlas Contributors"

from . import annotate, doublets, export, integration, io, qc, utils, viz

__all__ = [
    "io",
    "qc", 
    "doublets",
    "integration",
    "annotate",
    "viz",
    "export",
    "utils",
]
