import os
"""Verify the wood-gradient glyph shader compiles AND produces the radial albedo
gradient. Viewed END-ON (cap facing camera) on a fat cylinder so the full radial
disc (light core -> dark rim) is visible. Ambient-only so red spread == albedo.

  1. route-C per-vertex : the look we must match
  2. glyph no shader     : flat (confirmed colour loss)
  3. glyph + wood shader : gradient recovered
Metric: fraction of wood pixels that are 'light core' (red > 130). Flat ~ 0.
"""
import sys
import numpy as np
import vtk
from vtkmodules.util.numpy_support import vtk_to_numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wood_shader import apply_wood_shader, unit_wood_cylinder, C_WOOD_LIGHT, C_WOOD_DARK

OUT = os.path.dirname(os.path.abspath(__file__))


def input_points():
    inp = vtk.vtkPolyData(); ip = vtk.vtkPoints(); ip.InsertNextPoint(0, 0, 0)
    inp.SetPoints(ip)
    scale = vtk.vtkFloatArray(); scale.SetNumberOfComponents(3); scale.SetName("scale")
    scale.InsertNextTuple3(3.0, 1.0, 3.0)  # fat, short disc
    inp.GetPointData().AddArray(scale)
    return inp


def flat(actor):
    p = actor.GetProperty()
    p.SetAmbient(1.0); p.SetDiffuse(0.0); p.SetSpecular(0.0)
    p.SetColor(0.36, 0.25, 0.20)


def shoot(actor, label, png, size=(240, 240)):
    ren = vtk.vtkRenderer(); ren.AddActor(actor); ren.SetBackground(0.1, 0.1, 0.1)
    rw = vtk.vtkRenderWindow(); rw.SetOffScreenRendering(1); rw.SetSize(*size)
    rw.AddRenderer(ren)
    cam = ren.GetActiveCamera()
    cam.SetPosition(0, 12, 0); cam.SetFocalPoint(0, 0.5, 0); cam.SetViewUp(0, 0, 1)
    ren.ResetCamera(); rw.Render()
    w2i = vtk.vtkWindowToImageFilter(); w2i.SetInput(rw); w2i.ReadFrontBufferOff()
    w2i.Update()
    wr = vtk.vtkPNGWriter(); wr.SetFileName("%s/%s" % (OUT, png))
    wr.SetInputConnection(w2i.GetOutputPort()); wr.Write()
    arr = vtk_to_numpy(w2i.GetOutput().GetPointData().GetScalars()).reshape(-1, 3).astype(int)
    wood = arr[(arr[:, 0] > 40) & (arr[:, 0] >= arr[:, 2])]
    light = int(np.sum(wood[:, 0] > 130))
    frac = light / max(1, len(wood))
    print("%-24s wood=%5d  light-core px=%4d (%.1f%%)  maxR=%3d" %
          (label, len(wood), light, 100 * frac, int(wood[:, 0].max()) if len(wood) else 0))
    return frac


cyl = unit_wood_cylinder(base_tri=40)
pts = vtk_to_numpy(cyl.GetPoints().GetData())
rad = np.sqrt(pts[:, 0] ** 2 + pts[:, 2] ** 2)
shade = 1.0 - np.clip(rad / 1.0, 0, 1)
cols = (np.array(C_WOOD_DARK)[None] * (1 - shade)[:, None]
        + np.array(C_WOOD_LIGHT)[None] * shade[:, None])

# 1. route-C reference
tf = vtk.vtkTransform(); tf.Scale(3.0, 1.0, 3.0)
tpf = vtk.vtkTransformPolyDataFilter(); tpf.SetInputData(cyl); tpf.SetTransform(tf); tpf.Update()
pc = vtk.vtkPolyData(); pc.DeepCopy(tpf.GetOutput())
ca = vtk.vtkUnsignedCharArray(); ca.SetNumberOfComponents(3); ca.SetName("Colors")
for c in cols:
    ca.InsertNextTypedTuple(tuple((c * 255).astype(int)))
pc.GetPointData().SetScalars(ca)
mm = vtk.vtkPolyDataMapper(); mm.SetInputData(pc)
mm.SetScalarModeToUsePointData(); mm.ScalarVisibilityOn()
ra = vtk.vtkActor(); ra.SetMapper(mm); flat(ra)
ref = shoot(ra, "1. route-C per-vertex", "grad_1_routeC.png")

# 2. glyph, no shader
m2 = vtk.vtkGlyph3DMapper(); m2.SetSourceData(cyl); m2.SetInputData(input_points())
m2.SetScaleArray("scale"); m2.SetScaleModeToScaleByVectorComponents(); m2.ScalarVisibilityOff()
a2 = vtk.vtkActor(); a2.SetMapper(m2); flat(a2)
noshader = shoot(a2, "2. glyph no shader", "grad_2_flat.png")

# 3. glyph + wood shader
m3 = vtk.vtkGlyph3DMapper(); m3.SetSourceData(cyl); m3.SetInputData(input_points())
m3.SetScaleArray("scale"); m3.SetScaleModeToScaleByVectorComponents(); m3.ScalarVisibilityOff()
a3 = vtk.vtkActor(); a3.SetMapper(m3); flat(a3)
apply_wood_shader(a3, radius=1.0)
b = shoot(a3, "3. glyph + wood shader", "grad_3_shader.png")

print()
ok = (noshader < 0.02) and (b > 0.5 * ref) and (b > 0.05)
print("VERDICT:", "GRADIENT RECOVERED (flat glyph had %.1f%%, shader %.1f%% vs route-C %.1f%%)"
      % (100 * noshader, 100 * b, 100 * ref) if ok else "FAILED")
