import os
"""Shadow interaction for route B (the doc's single most-likely silent break).

Two questions:
  1. Does the custom wood-gradient glyph shader COMPILE and render under the
     vtkShadowMapPass pipeline (baker depth pass + shadow pass)?  -> capture any
     'Could not set shader program' on stderr.
  2. Does a SWAYED instance's cast shadow MOVE?  Route B animates via the instance
     ARRAYS (data), not a vertex-shader palette, so the same glyph input drives the
     depth pass; the shadow should follow the sway for free. We prove it: put the
     forest above a ground plane, render at two very different sway phases, and
     confirm the ground shadow footprint changes.

Compared against route C (merged mesh) as the shadow-casting reference.
"""
import sys
import numpy as np
import vtk
from vtkmodules.util.numpy_support import numpy_to_vtk, vtk_to_numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import forest_geom as fg
import render_forest as rf

OUT = os.path.dirname(os.path.abspath(__file__))


def ground_actor(span):
    p = vtk.vtkPlaneSource(); p.SetOrigin(-span, -span, -0.2)
    p.SetPoint1(span, -span, -0.2); p.SetPoint2(-span, span, -0.2)
    p.SetResolution(1, 1); p.Update()
    m = vtk.vtkPolyDataMapper(); m.SetInputConnection(p.GetOutputPort())
    a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetColor(0.75, 0.72, 0.62)
    a.GetProperty().SetAmbient(0.3); a.GetProperty().SetDiffuse(0.8)
    return a


def render_at(build, t, tag, exaggerate=0.0):
    forest = fg.build_forest(9, seed=20260817, spacing=22.0)
    if exaggerate:
        for tr in forest:
            tr.sway = exaggerate  # big sway so the shadow motion is unmistakable
    actors, update = build(forest)
    update(t)
    span = 90.0
    actors = actors + [ground_actor(span)]
    ren, rw = rf.make_scene(actors, shadows=True, size=(700, 500))
    # a low sun-ish camera so long shadows fall across the ground toward us
    ren.ResetCamera()
    cam = ren.GetActiveCamera()
    b = [0] * 6; ren.ComputeVisiblePropBounds(b)
    cam.SetFocalPoint(0, 0, 15); cam.SetPosition(70, -120, 55); cam.SetViewUp(0, 0, 1)
    ren.ResetCameraClippingRange()
    rw.Render()
    img = rf.screenshot(rw, "%s/shadow_%s.png" % (OUT, tag))
    return img


def ground_shadow_mask(img):
    """Ground pixels (tan plane) that are darkened by shadow."""
    r, g, b = img[:, :, 0].astype(int), img[:, :, 1].astype(int), img[:, :, 2].astype(int)
    is_ground = (r > 90) & (r < 205) & (np.abs(r - g) < 40) & (b < r)
    lit = is_ground & (r > 165)
    shadowed = is_ground & (r <= 165)
    return shadowed, lit


print("=== route B: shadows, two sway phases ===")
b0 = render_at(rf.build_route_b, 0.0, "B_t0", exaggerate=0.25)
b1 = render_at(rf.build_route_b, 3.5, "B_t1", exaggerate=0.25)
s0, l0 = ground_shadow_mask(b0)
s1, l1 = ground_shadow_mask(b1)
moved = int(np.sum(s0 != s1))
print("route B ground-shadow px: t0=%d  t1=%d  changed=%d" % (s0.sum(), s1.sum(), moved))

print("=== route C reference: shadows, same two phases ===")
c0 = render_at(rf.build_route_c, 0.0, "C_t0", exaggerate=0.25)
c1 = render_at(rf.build_route_c, 3.5, "C_t1", exaggerate=0.25)
cs0, _ = ground_shadow_mask(c0)
cs1, _ = ground_shadow_mask(c1)
cmoved = int(np.sum(cs0 != cs1))
print("route C ground-shadow px: t0=%d  t1=%d  changed=%d" % (cs0.sum(), cs1.sum(), cmoved))

print()
print("route B casts shadows:", "YES" if s0.sum() > 200 else "NO")
print("route B shadow MOVES with sway:", "YES" if moved > 200 else "NO")
print("route C casts shadows:", "YES" if cs0.sum() > 200 else "NO")
