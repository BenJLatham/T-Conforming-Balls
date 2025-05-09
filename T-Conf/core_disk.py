"""
T-Conf Disk Core Module

High-level API for generating a symmetric, T-conforming annular disk mesh
using Gmsh. The mesh consists of a transfinite band of thickness 2*a
around radius rInner, and an outer annular bulk region.
"""

import math
import gmsh
from typing import List
from T_Conf.utils import ensure_gmsh_available
from T_Conf.transfinite import set_transfinite
import numpy as np


def Transfinite_Disk(
    a: float,
    rInner: float,
    rOuter: float,
    h_band: float = 0.1,
    h_outer: float = 0.5,
    n_sectors: int = 4,
    InnerMaterialName: str = "InnerMaterial",
    OuterMaterialName: str = "OuterMaterial",
    InterfaceName: str = "Interface",
    OuterBoundaryName: str = "Outer",
) -> None:
    """
    Build a 2D disk mesh with three regions.

    Raises
    ------
    ValueError
        If `a` is too large (rInner - a < 0 or rInner + a > rOuter).
    """
    ensure_gmsh_available()

    # Initialize Gmsh and add model if needed
    if not gmsh.isInitialized():
        gmsh.initialize()
    try:
        gmsh.model.getCurrent()
    except Exception:
        gmsh.model.add("TransfiniteDisk")

    # Parameter checks
    if rInner - a < 0:
        raise ValueError("rInner - a must be nonnegative.")
    if rInner + a > rOuter:
        raise ValueError("rInner + a must be <= rOuter.")

    # Transfinite parameters
    n_div_radial = max(4, int(np.ceil(2 * a / h_band)))
    n_div_angular = max(4, int(math.ceil((0.5 * math.pi) / h_band)))

    def meshSizeAt(x, y):
        r = math.sqrt(x * x + y * y)
        # if we’re in the band: fine
        if rInner - a <= r <= rInner + a:
            return h_band
        # otherwise coarser
        return h_outer

    inner_surfaces: List[int] = []
    outer_surfaces: List[int] = []
    p_center = gmsh.model.geo.addPoint(0, 0, 0, meshSizeAt(0, 0))

    # Helper to create a point
    def make_pt(r, theta):
        x, y = r * math.cos(theta), r * math.sin(theta)
        return gmsh.model.geo.addPoint(x, y, 0, meshSizeAt(x, y))

    n_sectors = 4
    for i in range(n_sectors):
        theta0 = 2 * math.pi * i / n_sectors
        theta1 = 2 * math.pi * (i + 1) / n_sectors
        # Points at key radii
        p = {
            "in0": make_pt(rInner - a, theta0),
            "in1": make_pt(rInner - a, theta1),
            "mid0": make_pt(rInner, theta0),
            "mid1": make_pt(rInner, theta1),
            "out0": make_pt(rInner + a, theta0),
            "out1": make_pt(rInner + a, theta1),
            "o20": make_pt(rOuter, theta0),
            "o21": make_pt(rOuter, theta1),
        }

        # inner bulk
        if rInner - a > 0:
            # inner band
            curves_ib = [
                gmsh.model.geo.addLine(p["in0"], p["mid0"]),
                gmsh.model.geo.addCircleArc(p["mid0"], p_center, p["mid1"]),
                gmsh.model.geo.addLine(p["mid1"], p["in1"]),
                gmsh.model.geo.addCircleArc(p["in1"], p_center, p["in0"]),
            ]
            arc_interface = curves_ib[1]
            loop_ib = gmsh.model.geo.addCurveLoop(curves_ib)
            surf_ib = gmsh.model.geo.addPlaneSurface([loop_ib])
            inner_surfaces.append(surf_ib)
            inner_arc = curves_ib[3]
            loop = gmsh.model.geo.addCurveLoop(
                [
                    gmsh.model.geo.addLine(p_center, p["in0"]),
                    -inner_arc,
                    gmsh.model.geo.addLine(p["in1"], p_center),
                ]
            )
            inner_surfaces.append(gmsh.model.geo.addPlaneSurface([loop]))
        elif rInner == a:
            curves_ib = [
                gmsh.model.geo.addLine(p_center, p["mid0"]),
                gmsh.model.geo.addCircleArc(p["mid0"], p_center, p["mid1"]),
                gmsh.model.geo.addLine(p["mid1"], p_center),
            ]
            arc_interface = curves_ib[1]
            loop_ib = gmsh.model.geo.addCurveLoop(curves_ib)
            surf_ib = gmsh.model.geo.addPlaneSurface([loop_ib])
            inner_surfaces.append(surf_ib)

        # outer band
        curves_ob = [
            gmsh.model.geo.addLine(p["mid0"], p["out0"]),
            gmsh.model.geo.addCircleArc(p["out0"], p_center, p["out1"]),
            gmsh.model.geo.addLine(p["out1"], p["mid1"]),
            -arc_interface,
        ]
        arc_outer = curves_ob[1]
        loop_ob = gmsh.model.geo.addCurveLoop(curves_ob)
        surf_ob = gmsh.model.geo.addPlaneSurface([loop_ob])
        outer_surfaces.append(surf_ob)
        if a + rInner < rOuter:
            # outer bulk
            l_b0 = gmsh.model.geo.addLine(p["out0"], p["o20"])
            arc_bound = gmsh.model.geo.addCircleArc(p["o20"], p_center, p["o21"])
            l_b1 = gmsh.model.geo.addLine(p["o21"], p["out1"])

            # Reuse the exact interface arc
            loop_bulk = gmsh.model.geo.addCurveLoop([l_b0, arc_bound, l_b1, -arc_outer])
            outer_surfaces.append(gmsh.model.geo.addPlaneSurface([loop_bulk]))

        set_transfinite(
            radial_curves=[curves_ib[0], curves_ib[2], curves_ob[0], curves_ob[2]],
            circular_curves=[curves_ib[1], curves_ib[3], curves_ob[1], curves_ob[3]],
            surfaces=[(surf_ib, "Right"), (surf_ob, "Left")],
            n_div_radial=n_div_radial,
            n_div_circular=n_div_angular,
        )

    # Tagging physical groups

    innerGroupTag = gmsh.model.addPhysicalGroup(2, inner_surfaces)
    gmsh.model.setPhysicalName(2, innerGroupTag, InnerMaterialName)

    outerGroupTag = gmsh.model.addPhysicalGroup(2, outer_surfaces)
    gmsh.model.setPhysicalName(2, outerGroupTag, OuterMaterialName)

    # Final sync
    gmsh.model.geo.synchronize()
