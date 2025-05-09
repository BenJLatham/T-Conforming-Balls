"""
T-Conf: Symmetric Transfinite Mesh Generation

This package provides utilities to construct structured, T-conforming
meshes for concentric disk and spherical geometries using Gmsh.

Submodules:
- core_disk: High-level meshing for 2D annular disks.
- core_sphere: High-level meshing for 3D spherical domains.
- transfinite: Low-level transfinite mesh helpers.
- mesh: Sweep-based volume construction for spherical meshes.
- utils: Gmsh initialization and diagnostic tools.
"""

from .core_disk import Transfinite_Disk
from .core_sphere import Transfinite_Sphere
from .utils import ensure_gmsh_available

__all__ = [
    "Transfinite_Disk",
    "Transfinite_Sphere",
    "ensure_gmsh_available",
]
