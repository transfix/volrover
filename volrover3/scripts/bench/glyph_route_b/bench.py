import os
"""fps/cost comparison: route B (2 glyph draw calls) vs route C (2 actors/tree).

Trees only (no terrain/volumes), 900x600 offscreen, 70 trees like the real forest.
So absolute fps is HIGHER than the doc's whole-scene numbers — the comparison that
matters is B-vs-C on THIS identical harness, which isolates the draw-approach.

Per frame both routes re-pose EVERY tree (no stagger) so the animation path is
stressed equally, and we time two phases separately:
  update : CPU re-pose + push to GPU-side structures (route C: overwrite points +
           Modified per tree; route B: recompute instance arrays + set them)
  render : rw.Render() — the GL draw cost (140 draw calls vs 2)
Reports median over many frames; run repeatedly (the doc warns of ~2x variance).
"""
import sys
import time

import numpy as np
import vtk

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import forest_geom as fg
import render_forest as rf

N_TREES = int(sys.argv[1]) if len(sys.argv) > 1 else 70
FRAMES = 60
WARMUP = 8


def gl_info(rw):
    rw.Render()
    return rw.ReportCapabilities().split("OpenGL renderer string:")[-1].splitlines()[0].strip()


def median(xs):
    return float(np.median(xs))


def bench(build, label, forest):
    actors, update = build(forest)
    ren, rw = rf.make_scene(actors, shadows=False, size=(900, 600))
    rf.frame_island_camera(ren, forest)
    n_actors = len(actors)
    zb = vtk.vtkFloatArray()

    def render_sync():
        rw.Render()
        rw.GetZbufferData(450, 300, 450, 300, zb)  # 1-px readback -> blocks on GPU

    for i in range(WARMUP):
        update(0.05 * i); render_sync()
    up, rd = [], []
    for i in range(FRAMES):
        t = 1.0 + 0.05 * i
        a = time.perf_counter(); update(t); b = time.perf_counter()
        render_sync(); c = time.perf_counter()
        up.append((b - a) * 1e3); rd.append((c - b) * 1e3)
    mu, mr = median(up), median(rd)
    # sustained throughput: K frames back-to-back, ONE gpu finish at the end
    # (pipelines frames like a real render loop, no per-frame readback stall)
    K = 60
    t0 = time.perf_counter()
    for i in range(K):
        update(2.0 + 0.05 * i); rw.Render()
    rw.GetZbufferData(450, 300, 450, 300, zb)  # single finish
    sustained_ms = (time.perf_counter() - t0) * 1e3 / K
    fps = 1000.0 / sustained_ms
    print("%-9s actors=%3d  update=%6.2f  render+sync=%6.2f  sustained=%6.2f ms => %5.1f fps"
          % (label, n_actors, mu, mr, sustained_ms, fps))
    return mu, mr, fps, sustained_ms, n_actors


forest = fg.build_forest(N_TREES, seed=20260817)
nseg = sum(t.n_seg for t in forest)
nleaf = sum(t.n_leaf for t in forest)
_, rw0 = rf.make_scene([], size=(64, 64))
print("GL renderer:", gl_info(rw0))
print("forest: %d trees, %d wood segments, %d needle stars\n" % (N_TREES, nseg, nleaf))

cU = bench(rf.build_route_c, "route C", forest)
bU = bench(rf.build_route_b, "route B", forest)
print()
print("render+sync speedup (C/B): %.1fx   sustained fps B/C: %.2fx (%.1f vs %.1f fps)"
      % (cU[1] / bU[1], cU[3] / bU[3], bU[2], cU[2]))
