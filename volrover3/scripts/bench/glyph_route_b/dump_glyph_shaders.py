"""Dump the ACTUAL generated glyph fragment+vertex shaders.

Technique: VTK prints the full shader source (with line numbers) to stderr
whenever vtkShaderProgram fails to compile. We inject ONE deliberately-invalid
line at the end of the fragment Impl anchor to force the dump; everything ELSE
in the printed source is the real, already-substituted generated shader — the
real variable names and the real //VTK anchors we need to match.
"""
import sys
import vtk

cyl = vtk.vtkCylinderSource()
cyl.SetResolution(12); cyl.SetRadius(1.0); cyl.SetHeight(1.0); cyl.Update()

inp = vtk.vtkPolyData(); ip = vtk.vtkPoints()
for x in (-3, 0, 3):
    ip.InsertNextPoint(x, 0, 0)
inp.SetPoints(ip)

m = vtk.vtkGlyph3DMapper()
m.SetSourceConnection(cyl.GetOutputPort())
m.SetInputData(inp)
m.ScalarVisibilityOff()

a = vtk.vtkActor(); a.SetMapper(m)
a.GetProperty().SetColor(0.36, 0.25, 0.20)

sp = a.GetShaderProperty()
# Force a compile failure at the very end of the fragment shader so VTK dumps
# the WHOLE generated fragment shader (with all real substitutions) to stderr.
sp.AddFragmentShaderReplacement(
    "//VTK::Color::Impl", True,
    "//VTK::Color::Impl\n  FORCE_A_DUMP_OF_THE_REAL_SHADER = 1.0;", False)

ren = vtk.vtkRenderer(); ren.AddActor(a); ren.SetBackground(0.1, 0.1, 0.1)
ren.AddLight(vtk.vtkLight())
rw = vtk.vtkRenderWindow(); rw.SetOffScreenRendering(1); rw.SetSize(240, 120)
rw.AddRenderer(ren)
cam = ren.GetActiveCamera(); cam.SetPosition(0, 0, 20); cam.SetFocalPoint(0, 0, 0)
ren.ResetCamera()
rw.Render()
print("=== render done ===", file=sys.stderr)
