"""Bark on the route-B wood glyph: procedural FRAGMENT bump mapping and/or VERTEX
displacement, layered on top of the recovered wood-colour gradient.

Both read the same cheap procedural bark height barkH() (vertical furrows + axial
breakup, a handful of sines — no texture, so nothing extra to upload).

  bump      : perturb the view normal per-FRAGMENT via the surface-gradient method
              (Mikkelsen 2010) — dFdx/dFdy of barkH and view position, no tangents
              needed. Fragment-shader cost, ~independent of tessellation.
  displace  : move each SOURCE vertex along its normal by barkH before the glyph
              transform — needs a finely-tessellated source (a pentagon has nothing
              to displace), so its cost is (source verts x instances) of GPU vertex
              work + the displacement ALU. Still ONE draw call, and the hi-res
              source is uploaded ONCE (route C would have to bake it into every
              tree and re-upload it every frame).

Colour injection stays at //VTK::Normal::Impl (see wood_shader.py for why).
"""
import math

import numpy as np
import vtk
from vtkmodules.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray

C_WOOD_LIGHT = (0.6549, 0.4901, 0.2392)
C_WOOD_DARK = (0.3607, 0.2510, 0.2000)

# procedural bark height, shared by vertex (displace) and fragment (bump).
BARK_GLSL = """
float barkH(vec3 p){
  float theta = atan(p.z, p.x);
  float y = p.y;
  float f = 0.0;
  f += 0.50*sin(theta*9.0  + 1.7*sin(y*3.0));
  f += 0.30*sin(theta*19.0 + 3.0);
  f += 0.20*sin(theta*37.0 + y*6.0);
  f += 0.15*sin(y*40.0 + theta*4.0);
  return f;
}
"""


def unit_wood_cylinder_res(sides, rings):
    """A tessellated unit cylinder (base y=0, top y=1, radius 1, ring in XZ) with
    `sides` around and `rings` segments along the length, plus cap-centre vertices
    at radius 0. `sides=5, rings=1` reproduces the forest _CYL topology."""
    ang = np.arange(sides) * 2 * math.pi / sides
    ring_xz = np.column_stack((np.cos(ang), np.sin(ang)))
    ys = np.linspace(0.0, 1.0, rings + 1)
    pts = [np.array([0.0, 0.0, 0.0])]          # 0: bottom centre
    for yi in ys:
        for c, s in ring_xz:
            pts.append(np.array([c, yi, s]))
    pts.append(np.array([0.0, 1.0, 0.0]))      # last: top centre
    pts = np.array(pts)
    top_c = len(pts) - 1

    def ring_idx(r, s):
        return 1 + r * sides + (s % sides)

    tris = []
    for s in range(sides):                      # bottom cap fan
        tris += [0, ring_idx(0, s + 1), ring_idx(0, s)]
    for s in range(sides):                       # top cap fan
        tris += [top_c, ring_idx(rings, s), ring_idx(rings, s + 1)]
    for r in range(rings):                        # walls
        for s in range(sides):
            a, b = ring_idx(r, s), ring_idx(r, s + 1)
            c, d = ring_idx(r + 1, s), ring_idx(r + 1, s + 1)
            tris += [a, c, d, a, d, b]
    return pts, np.array(tris)


def cylinder_polydata(sides, rings):
    pts, tris = unit_wood_cylinder_res(sides, rings)
    pd = vtk.vtkPolyData()
    vp = vtk.vtkPoints(); vp.SetData(numpy_to_vtk(np.ascontiguousarray(pts), deep=1))
    pd.SetPoints(vp)
    conn = np.column_stack((np.full(len(tris) // 3, 3), tris.reshape(-1, 3))).ravel().astype(np.int64)
    ca = vtk.vtkCellArray(); ca.SetCells(len(tris) // 3, numpy_to_vtkIdTypeArray(conn, deep=1))
    pd.SetPolys(ca)
    nf = vtk.vtkPolyDataNormals(); nf.SetInputData(pd)
    nf.SplittingOff(); nf.ConsistencyOn(); nf.ComputePointNormalsOn(); nf.Update()
    return nf.GetOutput()


def apply_bark_shader(actor, bump=False, displace=False,
                      bump_scale=0.35, disp_scale=0.06, radius=1.0):
    sp = actor.GetShaderProperty()
    # ---- vertex ----
    sp.AddVertexShaderReplacement(
        "//VTK::Normal::Dec", True,
        "//VTK::Normal::Dec\nout float woodShade;\nout vec3 barkCoordVS;\n" + BARK_GLSL, False)
    if displace:
        # displace the SOURCE vertex along its normal BEFORE the glyph transform.
        sp.AddVertexShaderReplacement(
            "vec4 vertex = GCMCMatrix * vertexMC;", True,
            "vec3 dispMC = vertexMC.xyz + normalMC * (%f) * barkH(vertexMC.xyz);\n"
            "  vec4 vertex = GCMCMatrix * vec4(dispMC, 1.0);" % disp_scale, False)
    sp.AddVertexShaderReplacement(
        "//VTK::PositionVC::Impl", True,
        "//VTK::PositionVC::Impl\n"
        "  woodShade = 1.0 - clamp(length(vertexMC.xz)/%f, 0.0, 1.0);\n"
        "  barkCoordVS = vertexMC.xyz;" % radius, False)
    # ---- fragment ----
    sp.AddFragmentShaderReplacement(
        "//VTK::Normal::Dec", True,
        "//VTK::Normal::Dec\nin float woodShade;\nin vec3 barkCoordVS;\n" + BARK_GLSL, False)
    impl = ("//VTK::Normal::Impl\n"
            "  vec3 woodC = mix(vec3(%f,%f,%f), vec3(%f,%f,%f), woodShade);\n"
            "  ambientColor = ambientIntensity * woodC;\n"
            "  diffuseColor = diffuseIntensity * woodC;\n" % (C_WOOD_DARK + C_WOOD_LIGHT))
    if bump:
        impl += (
            "  {\n"
            "    float h = barkH(barkCoordVS);\n"
            "    vec3 sS = dFdx(vertexVC.xyz);\n"
            "    vec3 sT = dFdy(vertexVC.xyz);\n"
            "    vec3 vn = normalVCVSOutput;\n"
            "    vec3 R1 = cross(sT, vn);\n"
            "    vec3 R2 = cross(vn, sS);\n"
            "    float det = dot(sS, R1);\n"
            "    vec3 sg = sign(det) * (dFdx(h)*R1 + dFdy(h)*R2);\n"
            "    normalVCVSOutput = normalize(abs(det)*vn - (%f)*sg);\n"
            "  }\n" % bump_scale)
    sp.AddFragmentShaderReplacement("//VTK::Normal::Impl", True, impl, False)
