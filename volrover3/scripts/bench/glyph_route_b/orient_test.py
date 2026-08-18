import os
"""Nail the vtkGlyph3DMapper QUATERNION orientation convention.

Source is a thin unit cylinder along +Y. We feed a quaternion that should rotate
+Y -> +X and render from +Z. If the cylinder lies HORIZONTAL (wide, short) the
convention is confirmed; if VERTICAL, the component order / matrix convention is
wrong. Tests (w,x,y,z) vs (x,y,z,w) array order and both matrix conventions.
"""
import sys
import numpy as np
import vtk
from vtkmodules.util.numpy_support import vtk_to_numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wood_shader import unit_wood_cylinder

print("QUATERNION mode available:", hasattr(vtk.vtkGlyph3DMapper(), "SetOrientationModeToQuaternion"))


def quat_wxyz_from_R(R):
    """Quaternion (w,x,y,z) from a 3x3 rotation matrix acting on COLUMN vectors."""
    t = np.trace(R)
    if t > 0:
        s = math.sqrt(t + 1.0) * 2
        w = 0.25 * s
        x = (R[2, 1] - R[1, 2]) / s
        y = (R[0, 2] - R[2, 0]) / s
        z = (R[1, 0] - R[0, 1]) / s
    elif R[0, 0] > R[1, 1] and R[0, 0] > R[2, 2]:
        s = math.sqrt(1.0 + R[0, 0] - R[1, 1] - R[2, 2]) * 2
        w = (R[2, 1] - R[1, 2]) / s; x = 0.25 * s
        y = (R[0, 1] + R[1, 0]) / s; z = (R[0, 2] + R[2, 0]) / s
    elif R[1, 1] > R[2, 2]:
        s = math.sqrt(1.0 + R[1, 1] - R[0, 0] - R[2, 2]) * 2
        w = (R[0, 2] - R[2, 0]) / s; x = (R[0, 1] + R[1, 0]) / s
        y = 0.25 * s; z = (R[1, 2] + R[2, 1]) / s
    else:
        s = math.sqrt(1.0 + R[2, 2] - R[0, 0] - R[1, 1]) * 2
        w = (R[1, 0] - R[0, 1]) / s; x = (R[0, 2] + R[2, 0]) / s
        y = (R[1, 2] + R[2, 1]) / s; z = 0.25 * s
    return np.array([w, x, y, z])


import math
# R rotating column vector +Y -> +X : rotation about Z by -90 deg
R = np.array([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
assert np.allclose(R @ np.array([0, 1.0, 0]), [1, 0, 0])
q_wxyz = quat_wxyz_from_R(R)
print("q(w,x,y,z) =", np.round(q_wxyz, 4))

cyl = unit_wood_cylinder(base_tri=12)


def render_extent(quat4, label):
    inp = vtk.vtkPolyData(); ip = vtk.vtkPoints(); ip.InsertNextPoint(0, 0, 0)
    inp.SetPoints(ip)
    sc = vtk.vtkFloatArray(); sc.SetNumberOfComponents(3); sc.SetName("scale")
    sc.InsertNextTuple3(0.3, 6.0, 0.3)  # thin + long along source +Y
    inp.GetPointData().AddArray(sc)
    q = vtk.vtkFloatArray(); q.SetNumberOfComponents(4); q.SetName("quat")
    q.InsertNextTuple4(*quat4)
    inp.GetPointData().AddArray(q)
    m = vtk.vtkGlyph3DMapper(); m.SetSourceData(cyl); m.SetInputData(inp)
    m.SetScaleArray("scale"); m.SetScaleModeToScaleByVectorComponents()
    m.SetOrientationArray("quat"); m.SetOrientationModeToQuaternion()
    m.OrientOn(); m.ScalarVisibilityOff()
    a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetColor(0.6, 0.4, 0.2)
    ren = vtk.vtkRenderer(); ren.AddActor(a); ren.SetBackground(0, 0, 0)
    rw = vtk.vtkRenderWindow(); rw.SetOffScreenRendering(1); rw.SetSize(200, 200)
    rw.AddRenderer(ren)
    cam = ren.GetActiveCamera(); cam.SetPosition(0, 0, 30); cam.SetFocalPoint(0, 0, 0)
    cam.SetViewUp(0, 1, 0); ren.ResetCamera(); rw.Render()
    w2i = vtk.vtkWindowToImageFilter(); w2i.SetInput(rw); w2i.ReadFrontBufferOff(); w2i.Update()
    arr = vtk_to_numpy(w2i.GetOutput().GetPointData().GetScalars()).reshape(200, 200, 3)
    mask = arr[:, :, 0] > 40
    ys, xs = np.where(mask)
    if len(xs) == 0:
        print("%-28s (nothing rendered)" % label); return
    w_ext = xs.max() - xs.min(); h_ext = ys.max() - ys.min()
    orient = "HORIZONTAL(+X) OK" if w_ext > 1.6 * h_ext else (
        "VERTICAL(+Y) wrong" if h_ext > 1.6 * w_ext else "square/ambiguous")
    print("%-28s width=%3d height=%3d -> %s" % (label, w_ext, h_ext, orient))


render_extent(tuple(q_wxyz), "quat as (w,x,y,z)")
render_extent((q_wxyz[1], q_wxyz[2], q_wxyz[3], q_wxyz[0]), "quat as (x,y,z,w)")
# sanity: identity quaternion should stay VERTICAL (+Y)
render_extent((1, 0, 0, 0), "identity (w,x,y,z) [expect VERT]")
