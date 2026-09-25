"""AlterSeeK-Path: general-k path generation for VASP/QE/ABINIT.

Public API:
    from alterseek import spin_symmetry, brillouin_zone, general_k_path
    from alterseek import run_workflow
    from alterseek import SpinSymmetryError

Everything else is internal and may change between versions.
"""

from .version import __version__

# Resolved on first use so that importing the package stays cheap.
_PUBLIC = {
    "spin_symmetry": ("alterseek.find_sf_operations", "run"),
    "brillouin_zone": ("alterseek.kpoints", "brillouin_zone"),
    "general_k_path": ("alterseek.kpoints", "general_k_path"),
    "run_workflow": ("alterseek.workflow", "run_workflow"),
    "SpinSymmetryError": ("alterseek.find_sf_operations", "SpinSymmetryError"),
}

__all__ = ["__version__", *_PUBLIC]


def __getattr__(name):
    if name in _PUBLIC:
        from importlib import import_module

        module_name, attribute = _PUBLIC[name]
        return getattr(import_module(module_name), attribute)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(__all__)
