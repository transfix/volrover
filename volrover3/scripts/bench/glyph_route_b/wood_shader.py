"""Reusable wood-gradient glyph shader — the recovered per-vertex colour for route B.

Route C paints each cylinder vertex: C_WOOD_LIGHT on the axis (cap centres),
C_WOOD_DARK at the rim — a radial light-core/dark-ring gradient. Glyph instancing
drops the source's per-vertex colours, so we recompute the SAME radial shade in a
custom shader from the glyph-local radius length(vertexMC.xz) (unit cylinder radius
== 1, axis == Y), and set the material base colour BEFORE lighting so the wood is
lit exactly like route C's mesh.

WHY THE PRIOR SPIKE FAILED: actor shader-property replacements are applied BEFORE
the mapper's own //VTK::Color::Impl substitution, so replacing //VTK::Color::Impl
deletes the block that DECLARES ambientColor/diffuseColor -> "undefined variable".
We instead inject at //VTK::Normal::Impl (after the colour decls, before the light
math), where those variables are in scope. Verified against the dumped 9.5 shader.
"""
import math

import numpy as np
import vtk
from vtkmodules.util.numpy_support import numpy_to_vtk

C_WOOD_LIGHT = (0.6549, 0.4901, 0.2392)
C_WOOD_DARK = (0.3607, 0.2510, 0.2000)
BASE_TRI = 5  # pentagon cross-section — same as lsystem_forest's _CYL


def apply_wood_shader(actor, radius=1.0):
    """Install the radial wood-gradient replacements on `actor`'s ShaderProperty.

    The glyph helper reads the ACTOR's shader property (vtkGlyph3DMapper has no
    GetShaderProperty of its own), so this is where the replacements must go.
    """
    sp = actor.GetShaderProperty()

    # --- vertex: carry a per-vertex radial shade (1 on axis -> 0 at rim) ---
    sp.AddVertexShaderReplacement(
        "//VTK::Normal::Dec", True,
        "//VTK::Normal::Dec\nout float woodShade;", False)
    sp.AddVertexShaderReplacement(
        "//VTK::PositionVC::Impl", True,
        "//VTK::PositionVC::Impl\n"
        "  woodShade = 1.0 - clamp(length(vertexMC.xz)/%f, 0.0, 1.0);" % radius,
        False)

    # --- fragment: mix DARK<->LIGHT by shade, set base colour before lighting ---
    sp.AddFragmentShaderReplacement(
        "//VTK::Normal::Dec", True,
        "//VTK::Normal::Dec\nin float woodShade;", False)
    # //VTK::Normal::Impl survives to the final source as a no-op comment and sits
    # after the ambientColor/diffuseColor declarations, before the light math.
    sp.AddFragmentShaderReplacement(
        "//VTK::Normal::Impl", True,
        "//VTK::Normal::Impl\n"
        "  vec3 woodC = mix(vec3(%f,%f,%f), vec3(%f,%f,%f), woodShade);\n"
        "  ambientColor = ambientIntensity * woodC;\n"
        "  diffuseColor = diffuseIntensity * woodC;" % (
            C_WOOD_DARK + C_WOOD_LIGHT),
        False)


def _unit_cyl_arrays(base_tri=BASE_TRI):
    """The forest's _CYL as a UNIT cylinder (base y=0, top y=1, radius 1, ring in
    XZ), with explicit cap-CENTRE vertices at radius 0 and triangle-fan caps — the
    topology route C paints (light core on axis, dark rim). vtkCylinderSource has
    NO centre vertex (its cap is one n-gon over the rim), so it renders no light
    core; we build the real thing. Returns (points Nx3, tris flat)."""
    ring = np.column_stack((np.cos(np.arange(base_tri) * 2 * math.pi / base_tri),
                            np.sin(np.arange(base_tri) * 2 * math.pi / base_tri)))
    cyl_v = 2 * base_tri + 2
    pts = np.zeros((cyl_v, 3))
    pts[1:base_tri + 1, 0::2] = ring          # bottom ring (y=0)
    pts[base_tri + 2:, 0::2] = ring           # top ring
    pts[base_tri + 1:, 1] = 1.0               # top half y=1 (centre + top ring)
    tris = []
    for i in range(base_tri):
        b0, b1 = 1 + i, 1 + (i + 1) % base_tri
        t0, t1 = base_tri + 2 + i, base_tri + 2 + (i + 1) % base_tri
        tris += [0, b0, b1, base_tri + 1, t1, t0, b0, t1, b1, b0, t0, t1]
    return pts, np.array(tris)


def unit_wood_cylinder(base_tri=BASE_TRI):
    """vtkPolyData of the forest unit cylinder (see _unit_cyl_arrays)."""
    pts, tris = _unit_cyl_arrays(base_tri)
    pd = vtk.vtkPolyData()
    vpts = vtk.vtkPoints()
    vpts.SetData(numpy_to_vtk(np.ascontiguousarray(pts), deep=1))
    pd.SetPoints(vpts)
    cells = vtk.vtkCellArray()
    tri = tris.reshape(-1, 3)
    for a, b, c in tri:
        cells.InsertNextCell(3)
        cells.InsertCellPoint(int(a)); cells.InsertCellPoint(int(b)); cells.InsertCellPoint(int(c))
    pd.SetPolys(cells)
    return pd
