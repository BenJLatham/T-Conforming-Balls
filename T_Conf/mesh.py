"""
T-Conf Sphere Mesh Helpers

This module defines functions to construct spherical mesh volumes by
revolving surface curves and assembling wedge-shaped regions, faithfully
replicating the logic in the original T_conf_sphere.py.
"""

from typing import Sequence, Optional
import gmsh
import numpy as np
from T_Conf.utils import ensure_gmsh_available, find_missing_numbers
from T_Conf.transfinite import set_transfinite


def revolve_surface(
    revolve_axis: Sequence[float],
    revolve_center: Sequence[float],
    angle_sweep: float,
    surface_loop: int,
    n_div_radial: int,
    n_div_curved: int,
    bumpcoef: float,
    surface_type: str,
    apply_transfinite: bool = True,
) -> (int, int):
    """
    Revolve a 2D surface loop to create positive and negative sweep volumes,
    and optionally apply transfinite meshing on the swept surfaces.

    Parameters
    ----------
    revolve_axis : Sequence[float]
        Direction vector of the revolution axis (x, y, z).
    revolve_center : Sequence[float]
        A point on the revolution axis (x, y, z).
    angle_sweep : float
        Sweep angle in radians for positive and negative revolutions.
    surface_loop : int
        ID of the 2D curve loop defining the surface to revolve.
    n_div_radial : int
        Number of divisions for radial (straight) edges.
    n_div_curved : int
        Number of divisions for curved (circular) edges.
    bumpcoef : float
        Coefficient to bias the transfinite curve progression.
    surface_type : {'inner', 'outer'}
        Flags arrangement logic ('Right' for inner, computed for outer).
    apply_transfinite : bool
        If True, apply surface transfinite meshing to new surfaces.

    Returns
    -------
    (vol_plus, vol_minus) : Tuple[int, int]
        The volume IDs from the positive and negative sweeps.
    """
    ensure_gmsh_available()

    # Positive sweep
    out_plus = gmsh.model.geo.revolve(
        [(2, surface_loop)],
        revolve_axis[0],
        revolve_axis[1],
        revolve_axis[2],
        revolve_center[0],
        revolve_center[1],
        revolve_center[2],
        +angle_sweep,
    )
    vol_plus = next(v for (dim, v) in out_plus if dim == 3)

    # Negative sweep
    out_minus = gmsh.model.geo.revolve(
        [(2, surface_loop)],
        revolve_axis[0],
        revolve_axis[1],
        revolve_axis[2],
        revolve_center[0],
        revolve_center[1],
        revolve_center[2],
        -angle_sweep,
    )
    vol_minus = next(v for (dim, v) in out_minus if dim == 3)

    # Apply transfinite meshing on surfaces
    if apply_transfinite:
        _apply_transfinite_sweep(out_plus, n_div_radial, n_div_curved, bumpcoef, surface_type)
        _apply_transfinite_sweep(out_minus, n_div_radial, n_div_curved, bumpcoef, surface_type)

    return vol_plus, vol_minus


def _apply_transfinite_sweep(
    sweep_results: Sequence[tuple],
    n_div_radial: int,
    n_div_curved: int,
    bumpcoef: float,
    surface_type: str,
) -> None:
    """
    Internal helper to process newly created surfaces from a revolve sweep
    and apply transfinite mesh settings.
    """
    # Identify new surface tags
    new_surfs = [tag for (dim, tag) in sweep_results if dim == 2]
    gmsh.model.geo.synchronize()

    for s_tag in new_surfs:
        # Determine arrangement
        if surface_type == "inner":
            use_arr = "Right"
        else:
            bbox = gmsh.model.getBoundingBox(2, s_tag)
            tol = 1e-6
            avg_z = (bbox[2] + bbox[5]) / 2.0
            use_arr = "Right" if (abs(avg_z) < tol and (bbox[0] >= 0 and bbox[3] >= 0)) else "Left"

        # Classify boundary edges
        bnd = gmsh.model.getBoundary([(2, s_tag)], combined=False, recursive=False)
        corner_points = []
        radial_curves = []
        curved_curves = []

        for dim_c, tag_c in bnd:
            if dim_c == 1:
                etype = gmsh.model.getType(1, abs(tag_c))
                if etype == "Circle":
                    curved_curves.append(tag_c)
                elif etype == "Line":
                    radial_curves.append(tag_c)
                # Collect endpoints
                endpts = gmsh.model.getBoundary([(1, tag_c)], combined=False, recursive=False)
                if len(endpts) == 2:
                    pt_a, pt_b = endpts[0][1], endpts[1][1]
                    for pt in (pt_a, pt_b):
                        if pt not in corner_points:
                            corner_points.append(pt)

        # Apply curve transfinite
        set_transfinite(
            radial_curves, curved_curves, [(s_tag, use_arr)], n_div_radial, n_div_curved, bumpcoef
        )

        # Apply surface transfinite arrangement and corner tags
        if len(corner_points) == 3:
            if 1 in corner_points:
                i = corner_points.index(1)
                corner_points = corner_points[i:] + corner_points[:i]
                gmsh.model.geo.mesh.setTransfiniteSurface(
                    s_tag, arrangement=use_arr, cornerTags=corner_points
                )
            else:
                _, missing = find_missing_numbers(corner_points, np.arange(2, 21))
                i = corner_points.index(missing[0])
                corner_points = corner_points[i:] + corner_points[:i]
                gmsh.model.geo.mesh.setTransfiniteSurface(s_tag, cornerTags=corner_points)
        else:
            gmsh.model.geo.mesh.setTransfiniteSurface(s_tag, arrangement=use_arr)


def create_wedge_volumes(
    inner_surf: int,
    outer_surf: int,
    inner_bulk_surf: Optional[int],
    outer_bulk_surf: Optional[int],
    revolve_axis: Sequence[float],
    revolve_center: Sequence[float],
    angle_sweep: float,
    n_div_radial: int,
    n_div_curved: int,
    bumpcoef: float,
) -> (Sequence[Optional[int]], Sequence[Optional[int]]):
    """
    Create inner and outer wedge volumes by revolving band and bulk surfaces.
    """
    ensure_gmsh_available()

    # Band volumes (transfinite)
    vol_i_plus, vol_i_minus = revolve_surface(
        revolve_axis,
        revolve_center,
        angle_sweep,
        inner_surf,
        n_div_radial,
        n_div_curved,
        bumpcoef,
        "inner",
        True,
    )
    gmsh.model.geo.synchronize()
    vol_o_plus, vol_o_minus = revolve_surface(
        revolve_axis,
        revolve_center,
        angle_sweep,
        outer_surf,
        n_div_radial,
        n_div_curved,
        bumpcoef,
        "outer",
        True,
    )
    gmsh.model.geo.synchronize()

    # Bulk volumes (unstructured)
    if inner_bulk_surf is not None:
        vol_ib_plus, vol_ib_minus = revolve_surface(
            revolve_axis,
            revolve_center,
            angle_sweep,
            inner_bulk_surf,
            n_div_radial,
            n_div_curved,
            bumpcoef,
            "inner",
            False,
        )
        gmsh.model.geo.synchronize()
    else:
        vol_ib_plus = vol_ib_minus = None

    if outer_bulk_surf is not None:
        vol_ob_plus, vol_ob_minus = revolve_surface(
            revolve_axis,
            revolve_center,
            angle_sweep,
            outer_bulk_surf,
            n_div_radial,
            n_div_curved,
            bumpcoef,
            "outer",
            False,
        )
    else:
        vol_ob_plus = vol_ob_minus = None

    gmsh.model.geo.synchronize()

    inner_vols = [vol_i_plus, vol_i_minus, vol_ib_plus, vol_ib_minus]
    outer_vols = [vol_o_plus, vol_o_minus, vol_ob_plus, vol_ob_minus]
    return inner_vols, outer_vols
