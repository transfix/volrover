import os
"""Scaling bench B vs C from 70 to 300 trees. Vsync MUST be off
(__GL_SYNC_TO_VBLANK=0) or renders floor at 1/60 s. Reports, per tree count:
  update : median per-frame CPU pose (shared by both routes)
  render : pure GPU render (back-to-back renders, one finish) — C is O(actors), B O(1)
  fps    : sustained update+render (real frame)
"""
import sys
import time

import numpy as np
import vtk

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import forest_geom as fg
import render_forest as rf

COUNTS = [int(x) for x in sys.argv[1:]] or [70, 140, 200, 250, 300]
F = 20


def measure(build, forest):
    actors, update = build(forest)
    ren, rw = rf.make_scene(actors, shadows=False, size=(900, 600))
    rf.frame_island_camera(ren, forest)
    zb = vtk.vtkFloatArray()
    fin = lambda: rw.GetZbufferData(450, 300, 450, 300, zb)
    for i in range(6):
        update(0.05 * i); rw.Render()
    fin()
    ups = []
    for i in range(F):
        a = time.perf_counter(); update(1.0 + 0.05 * i); ups.append((time.perf_counter() - a) * 1e3)
    t0 = time.perf_counter()
    for _ in range(F):
        rw.Render()
    fin(); rend = (time.perf_counter() - t0) * 1e3 / F
    t0 = time.perf_counter()
    for i in range(F):
        update(2.0 + 0.05 * i); rw.Render()
    fin(); sus = (time.perf_counter() - t0) * 1e3 / F
    return float(np.median(ups)), rend, sus, len(actors)


print("trees  actors_C |    route C: upd  rend   fps  |    route B: upd  rend   fps  | rend C/B | fps B/C")
for n in COUNTS:
    forest = fg.build_forest(n, seed=20260817)
    cu, cr, cs, ca = measure(rf.build_route_c, forest)
    bu, br, bs, ba = measure(rf.build_route_b, forest)
    print("%5d  %6d   |  %6.1f %6.1f %5.1f  |  %6.1f %6.1f %5.1f  |  %6.1fx | %5.2fx"
          % (n, ca, cu, cr, 1000 / cs, bu, br, 1000 / bs, cr / max(br, 1e-6), (1000 / bs) / (1000 / cs)))
