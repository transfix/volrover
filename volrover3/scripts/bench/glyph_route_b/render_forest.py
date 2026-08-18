import os
"""Render the reproduced forest two ways in a shared VTK scene:

  route C : one merged wood polydata actor (per-vertex colour) + one needle line
            actor per tree; per frame overwrite points + Modified (full-VBO reupload).
  route B : ONE vtkGlyph3DMapper over a unit cylinder (all wood segments) + ONE
            over a unit needle star (all leaves); per frame rewrite the instance
            (pos/quat/scale) arrays. Wood colour recovered by the glyph shader.

Same renderer / lights / camera for both; only actor construction differs.
"""
import math
import sys

import numpy as np
import vtk
from vtkmodules.util.numpy_support import (numpy_to_vtk, numpy_to_vtkIdTypeArray,
                                            vtk_to_numpy)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import forest_geom as fg
from wood_shader import apply_wood_shader, unit_wood_cylinder, C_WOOD_DARK

SUN_AZ = -52.0


def _np3(a):
    return numpy_to_vtk(np.ascontiguousarray(a, dtype=np.float32), deep=1)


def _points(a):
    p = vtk.vtkPoints(); p.SetData(_np3(a)); return p


# ── needle star source (lines), unit radius/height ───────────────────────────
def unit_needle_star():
    stride = fg.NEEDLES + 1
    pts = np.zeros((stride, 3))
    pts[1:, 0::2] = fg._NRING
    pts[1:, 1] = 1.0
    pd = vtk.vtkPolyData(); pd.SetPoints(_points(pts))
    lines = vtk.vtkCellArray()
    for k in range(fg.NEEDLES):
        lines.InsertNextCell(2); lines.InsertCellPoint(0); lines.InsertCellPoint(1 + k)
    pd.SetLines(lines)
    return pd


def _with_normals(pd):
    nf = vtk.vtkPolyDataNormals(); nf.SetInputData(pd)
    nf.SplittingOff(); nf.ConsistencyOn(); nf.ComputePointNormalsOn(); nf.Update()
    return nf.GetOutput()


def _tris_cells(tris):
    """vtkCellArray of triangles from an (M,3) index array (legacy format)."""
    conn = np.column_stack((np.full(len(tris), 3), tris)).ravel().astype(np.int64)
    ca = vtk.vtkCellArray()
    ca.SetCells(len(tris), numpy_to_vtkIdTypeArray(conn, deep=1))
    return ca


def _lines_cells(lines):
    conn = np.column_stack((np.full(len(lines), 2), lines)).ravel().astype(np.int64)
    ca = vtk.vtkCellArray()
    ca.SetCells(len(lines), numpy_to_vtkIdTypeArray(conn, deep=1))
    return ca


# ── route C ──────────────────────────────────────────────────────────────────
def build_route_c(forest):
    actors = []
    wood_pds = []
    needle_pds = []
    for tree in forest:
        wpd = vtk.vtkPolyData(); wpd.SetPoints(_points(tree.wood_buf))
        wpd.SetPolys(_tris_cells(tree.wood_tris.reshape(-1, 3)))
        col = numpy_to_vtk(np.ascontiguousarray((tree.wood_colors * 255).astype(np.uint8)), deep=1)
        col.SetName("Colors")
        wpd.GetPointData().SetScalars(col)
        wpd2 = _with_normals(wpd)
        wm = vtk.vtkPolyDataMapper(); wm.SetInputData(wpd2)
        wm.SetScalarModeToUsePointData(); wm.ScalarVisibilityOn()
        wa = vtk.vtkActor(); wa.SetMapper(wm)
        actors.append(wa); wood_pds.append((wpd, wpd2, wm))
        if tree.needle_buf is not None:
            npd = vtk.vtkPolyData(); npd.SetPoints(_points(tree.needle_buf))
            npd.SetLines(_lines_cells(tree.needle_lines.reshape(-1, 2)))
            nm = vtk.vtkPolyDataMapper(); nm.SetInputData(npd); nm.ScalarVisibilityOff()
            na = vtk.vtkActor(); na.SetMapper(nm); na.GetProperty().SetColor(*fg.C_NEEDLE)
            actors.append(na); needle_pds.append((npd, nm))
        else:
            needle_pds.append(None)

    def update(t):
        wi = ni = 0
        for tree in forest:
            wb, nb = fg.repose_route_c(tree, t)
            wpd, wpd2, wm = wood_pds[wi]; wi += 1
            # overwrite the (normals-filtered) polydata points in place + Modified
            wpd2.GetPoints().SetData(_np3(wb)); wpd2.GetPoints().Modified(); wpd2.Modified()
            if needle_pds[ni] is not None:
                npd, nm = needle_pds[ni]
                npd.GetPoints().SetData(_np3(nb)); npd.GetPoints().Modified(); npd.Modified()
            ni += 1
    return actors, update


# ── route B ──────────────────────────────────────────────────────────────────
def build_route_b(forest, wood_res=fg.BASE_TRI):
    wood_src = _with_normals(unit_wood_cylinder(base_tri=wood_res))
    needle_src = unit_needle_star()

    wood_in = vtk.vtkPolyData()
    needle_in = vtk.vtkPolyData()

    def make_glyph(inp, src, is_line):
        m = vtk.vtkGlyph3DMapper(); m.SetSourceData(src); m.SetInputData(inp)
        m.SetScaleArray("scale"); m.SetScaleModeToScaleByVectorComponents()
        m.SetOrientationArray("quat"); m.SetOrientationModeToQuaternion(); m.OrientOn()
        m.ScalarVisibilityOff()
        return m

    wm = make_glyph(wood_in, wood_src, False)
    nm = make_glyph(needle_in, needle_src, True)
    wa = vtk.vtkActor(); wa.SetMapper(wm)
    wa.GetProperty().SetColor(*C_WOOD_DARK)
    apply_wood_shader(wa, radius=1.0)
    na = vtk.vtkActor(); na.SetMapper(nm); na.GetProperty().SetColor(*fg.C_NEEDLE)
    actors = [wa, na]

    def _fill(inp, pos, quat, scale):
        inp.SetPoints(_points(pos))
        sc = numpy_to_vtk(np.ascontiguousarray(scale, np.float32), deep=1); sc.SetName("scale")
        q = numpy_to_vtk(np.ascontiguousarray(quat, np.float32), deep=1); q.SetName("quat")
        inp.GetPointData().AddArray(sc); inp.GetPointData().AddArray(q)
        inp.Modified()

    def update(t):
        wood, needle = fg.route_b_instances(forest, t)
        _fill(wood_in, *wood)
        _fill(needle_in, *needle)
    update(0.0)
    return actors, update


# ── shared scene ─────────────────────────────────────────────────────────────
def sun_dir(az_deg, el_deg):
    az, el = math.radians(az_deg), math.radians(el_deg)
    return (math.cos(el) * math.sin(az), -math.cos(el) * math.cos(az), math.sin(el))


def make_scene(actors, shadows=False, size=(900, 600)):
    ren = vtk.vtkRenderer(); ren.SetBackground(0.62, 0.76, 0.92)
    for a in actors:
        ren.AddActor(a)
    ren.RemoveAllLights()
    ren.SetLightFollowCamera(False)
    for az, el, inten, col in ((SUN_AZ, 34.0, 1.0, (1.0, 0.94, 0.82)),
                               (128.0, 52.0, 0.55, (0.66, 0.85, 0.55))):
        L = vtk.vtkLight(); L.SetLightTypeToSceneLight(); L.SetPositional(False)
        d = sun_dir(az, el)
        L.SetPosition(d[0], d[1], d[2]); L.SetFocalPoint(0, 0, 0)
        L.SetColor(*col); L.SetIntensity(inten); ren.AddLight(L)
    rw = vtk.vtkRenderWindow(); rw.SetOffScreenRendering(1); rw.SetSize(*size)
    rw.AddRenderer(ren)
    if shadows:
        seq = vtk.vtkSequencePass(); passes = vtk.vtkRenderPassCollection()
        shadow = vtk.vtkShadowMapPass()
        passes.AddItem(shadow.GetShadowMapBakerPass()); passes.AddItem(shadow)
        seq.SetPasses(passes)
        cam = vtk.vtkCameraPass(); cam.SetDelegatePass(seq)
        ren.SetPass(cam)
    return ren, rw


def screenshot(rw, path):
    w2i = vtk.vtkWindowToImageFilter(); w2i.SetInput(rw); w2i.ReadFrontBufferOff(); w2i.Update()
    wr = vtk.vtkPNGWriter(); wr.SetFileName(path); wr.SetInputConnection(w2i.GetOutputPort()); wr.Write()
    return vtk_to_numpy(w2i.GetOutput().GetPointData().GetScalars()).reshape(rw.GetSize()[1], rw.GetSize()[0], 3)


def frame_island_camera(ren, forest):
    ren.ResetCamera()
    cam = ren.GetActiveCamera()
    b = [0, 0, 0, 0, 0, 0]; ren.ComputeVisiblePropBounds(b)
    cx, cy = (b[0] + b[1]) / 2, (b[2] + b[3]) / 2
    span = max(b[1] - b[0], b[3] - b[2])
    cam.SetFocalPoint(cx, cy, (b[4] + b[5]) / 2)
    cam.SetPosition(cx + 0.2 * span, cy - 1.1 * span, (b[4] + b[5]) / 2 + 0.7 * span)
    cam.SetViewUp(0, 0, 1)
    ren.ResetCameraClippingRange()
