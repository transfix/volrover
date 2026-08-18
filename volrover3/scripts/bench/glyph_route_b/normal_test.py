import os
"""Is route B's wood darker because vtkGlyph3DMapper's per-instance normal matrix
mishandles SMOOTH (beveled) normals under non-uniform per-axis scale?

Test: a strongly non-uniform cylinder (scale 0.4 x 6 x 0.4), side-lit. Compare
  ref  = real scaled mesh, normals recomputed AFTER scaling (geometrically correct)
  glyphS = glyph, source normals SMOOTH (SplittingOff) — beveled ring normals
  glyphH = glyph, source normals HARD  (SplittingOn) — all normals axial/radial,
           i.e. eigenvectors of the scale, immune to inverse-transpose error
If glyphH matches ref but glyphS is darker, the cause is confirmed and the fix is
hard normals on the wood source.
"""
import sys
import numpy as np
import vtk
from vtkmodules.util.numpy_support import vtk_to_numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wood_shader import unit_wood_cylinder

SC = (0.4, 6.0, 0.4)


def normals(pd, split):
    nf = vtk.vtkPolyDataNormals(); nf.SetInputData(pd)
    (nf.SplittingOn if split else nf.SplittingOff)()
    nf.ConsistencyOn(); nf.ComputePointNormalsOn(); nf.Update()
    return nf.GetOutput()


def brightness(actor, label):
    ren = vtk.vtkRenderer(); ren.AddActor(actor); ren.SetBackground(0, 0, 0)
    ren.RemoveAllLights(); ren.SetLightFollowCamera(False)
    L = vtk.vtkLight(); L.SetPositional(False)
    L.SetPosition(1, 0, 0.3); L.SetFocalPoint(0, 0, 0); L.SetColor(1, 1, 1); L.SetIntensity(1.0)
    ren.AddLight(L)
    rw = vtk.vtkRenderWindow(); rw.SetOffScreenRendering(1); rw.SetSize(200, 200); rw.AddRenderer(ren)
    cam = ren.GetActiveCamera(); cam.SetPosition(0, 0, 30); cam.SetFocalPoint(0, 2.5, 0)
    cam.SetViewUp(0, 1, 0); ren.ResetCamera(); rw.Render()
    w2i = vtk.vtkWindowToImageFilter(); w2i.SetInput(rw); w2i.ReadFrontBufferOff(); w2i.Update()
    arr = vtk_to_numpy(w2i.GetOutput().GetPointData().GetScalars()).reshape(-1, 3).astype(int)
    fg = arr[arr.max(1) > 20]
    print("%-10s lit wood mean R=%.1f  (px %d)" % (label, fg[:, 0].mean(), len(fg)))
    return fg[:, 0].mean()


cyl = unit_wood_cylinder(base_tri=5)

# ref: real scaled mesh
tf = vtk.vtkTransform(); tf.Scale(*SC)
tpf = vtk.vtkTransformPolyDataFilter(); tpf.SetInputData(cyl); tpf.SetTransform(tf); tpf.Update()
refpd = normals(tpf.GetOutput(), split=False)
mm = vtk.vtkPolyDataMapper(); mm.SetInputData(refpd); mm.ScalarVisibilityOff()
ra = vtk.vtkActor(); ra.SetMapper(mm); ra.GetProperty().SetColor(0.65, 0.49, 0.24)
ra.GetProperty().SetAmbient(0.0); ra.GetProperty().SetDiffuse(1.0)
ref = brightness(ra, "ref-mesh")


def glyph_actor(src):
    inp = vtk.vtkPolyData(); ip = vtk.vtkPoints(); ip.InsertNextPoint(0, 0, 0); inp.SetPoints(ip)
    s = vtk.vtkFloatArray(); s.SetNumberOfComponents(3); s.SetName("scale"); s.InsertNextTuple3(*SC)
    inp.GetPointData().AddArray(s)
    m = vtk.vtkGlyph3DMapper(); m.SetSourceData(src); m.SetInputData(inp)
    m.SetScaleArray("scale"); m.SetScaleModeToScaleByVectorComponents(); m.ScalarVisibilityOff()
    a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetColor(0.65, 0.49, 0.24)
    a.GetProperty().SetAmbient(0.0); a.GetProperty().SetDiffuse(1.0)
    return a


gS = brightness(glyph_actor(normals(cyl, split=False)), "glyph-smooth")
gH = brightness(glyph_actor(normals(cyl, split=True)), "glyph-hard")

print()
print("ratio glyph-smooth/ref = %.3f   glyph-hard/ref = %.3f" % (gS / ref, gH / ref))
print("CONCLUSION:",
      "beveled-normal inverse-transpose error CONFIRMED; hard normals fix it"
      if abs(gH - ref) < 0.05 * ref and gS < 0.92 * ref else
      "not the (sole) cause — investigate further")
