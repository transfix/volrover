import os
"""Disentangle the 20 shadow-pass shader errors: is it my WOOD glyph shader
failing in the baker, or the plain NEEDLE glyph (a VTK glyph+shadow issue)?

Renders four isolated cases under the shadow pipeline and reports whether each
compiled clean (no 'Could not set shader program' on stderr, captured by the
caller) and whether the beauty colours are right.
"""
import sys
import numpy as np
import vtk

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import forest_geom as fg
import render_forest as rf

OUT = os.path.dirname(os.path.abspath(__file__))
case = sys.argv[1]  # wood | needle | woodC | needleC

forest = fg.build_forest(9, seed=20260817, spacing=22.0)
bA, bU = rf.build_route_b(forest); bU(1.0)
cA, cU = rf.build_route_c(forest); cU(1.0)
wood_b, needle_b = bA[0], bA[1]
# route C: split wood actors (even idx) from needle actors (odd idx)
wood_c = [a for i, a in enumerate(cA) if i % 2 == 0]
needle_c = [a for i, a in enumerate(cA) if i % 2 == 1]

actors = {"wood": [wood_b], "needle": [needle_b],
          "woodC": wood_c, "needleC": needle_c}[case]

def gplane(span):
    p = vtk.vtkPlaneSource(); p.SetOrigin(-span, -span, -0.2)
    p.SetPoint1(span, -span, -0.2); p.SetPoint2(-span, span, -0.2); p.Update()
    m = vtk.vtkPolyDataMapper(); m.SetInputConnection(p.GetOutputPort())
    a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetColor(0.75, 0.72, 0.62)
    return a

ren, rw = rf.make_scene(actors + [gplane(90)], shadows=True, size=(500, 360))
ren.ResetCamera(); cam = ren.GetActiveCamera()
cam.SetFocalPoint(0, 0, 15); cam.SetPosition(70, -120, 55); cam.SetViewUp(0, 0, 1)
ren.ResetCameraClippingRange()
rw.Render()
img = rf.screenshot(rw, "%s/shiso_%s.png" % (OUT, case))
print("CASE %s rendered (see stderr for shader errors)" % case)
