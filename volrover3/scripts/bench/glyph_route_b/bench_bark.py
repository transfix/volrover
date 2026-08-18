import os
"""Performance of bark on the route-B wood glyph: fragment BUMP vs vertex
DISPLACEMENT vs both, at low vs high source tessellation.

Wood-only (needles excluded) so the number IS the wood-shader cost. Bump and
displacement change only the RENDER (the per-frame instance update is identical),
so we measure pure GPU render throughput: K back-to-back rw.Render() with one GPU
finish at the end — this divides out the ~1/60 s readback stall that otherwise
floors every synced frame.
"""
import sys
import time

import numpy as np
import vtk
from vtkmodules.util.numpy_support import numpy_to_vtk

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import forest_geom as fg
import render_forest as rf
from bark_shader import apply_bark_shader, cylinder_polydata

N_TREES = int(sys.argv[1]) if len(sys.argv) > 1 else 70
forest = fg.build_forest(N_TREES, seed=20260817)
nseg = sum(t.n_seg for t in forest)
wood, _ = fg.route_b_instances(forest, 1.0)
wpos, wquat, wscale = wood


def build_wood(sides, rings, bump, displace):
    src = cylinder_polydata(sides, rings)
    src_v = src.GetNumberOfPoints()
    inp = vtk.vtkPolyData(); inp.SetPoints(rf._points(wpos))
    sc = numpy_to_vtk(np.ascontiguousarray(wscale, np.float32), deep=1); sc.SetName("scale")
    q = numpy_to_vtk(np.ascontiguousarray(wquat, np.float32), deep=1); q.SetName("quat")
    inp.GetPointData().AddArray(sc); inp.GetPointData().AddArray(q)
    m = vtk.vtkGlyph3DMapper(); m.SetSourceData(src); m.SetInputData(inp)
    m.SetScaleArray("scale"); m.SetScaleModeToScaleByVectorComponents()
    m.SetOrientationArray("quat"); m.SetOrientationModeToQuaternion(); m.OrientOn()
    m.ScalarVisibilityOff()
    a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetColor(*fg.C_WOOD_DARK)
    apply_bark_shader(a, bump=bump, displace=displace)
    return a, src_v


def measure(label, sides, rings, bump, displace):
    a, src_v = build_wood(sides, rings, bump, displace)
    ren, rw = rf.make_scene([a], shadows=False, size=(900, 600))
    rf.frame_island_camera(ren, forest)
    zb = vtk.vtkFloatArray()
    for _ in range(6):
        rw.Render()
    rw.GetZbufferData(450, 300, 450, 300, zb)
    K = 120
    t0 = time.perf_counter()
    for _ in range(K):
        rw.Render()
    rw.GetZbufferData(450, 300, 450, 300, zb)  # single finish
    ms = (time.perf_counter() - t0) * 1e3 / K
    total_v = src_v * nseg
    print("%-14s src=%4d v  scene=%7.2fM v  render=%6.3f ms/frame  (%6.1f fps GPU-only)"
          % (label, src_v, total_v / 1e6, ms, 1000.0 / ms))
    return ms, src_v, a, ren, rw


print("wood-only, %d trees, %d segments\n" % (N_TREES, nseg))
base, _, _, _, _ = measure("base-lowres", 5, 1, False, False)
measure("base-hires", 24, 12, False, False)
measure("bump-lowres", 5, 1, True, False)
measure("bump-hires", 24, 12, True, False)
measure("disp-hires", 24, 12, False, True)
mbark, _, abark, renbark, rwbark = measure("bark-full", 24, 12, True, True)

# save a bark screenshot for the visual check
OUT = os.path.dirname(os.path.abspath(__file__))
rwbark.Render()
rf.screenshot(rwbark, "%s/bark_full.png" % OUT)
# and a close crop by moving the camera in
cam = renbark.GetActiveCamera()
b = [0] * 6; renbark.ComputeVisiblePropBounds(b)
cam.SetFocalPoint((b[0]+b[1])/2, (b[2]+b[3])/2, (b[4]+b[5])/2)
cam.Zoom(3.0); renbark.ResetCameraClippingRange(); rwbark.Render()
rf.screenshot(rwbark, "%s/bark_full_zoom.png" % OUT)
print("\nbark-full costs %.2f ms/frame over base (%.3f -> %.3f); on the ~68 ms route-B"
      " frame that is +%.1f%%." % (mbark - base, base, mbark, 100 * (mbark - base) / 68.0))
