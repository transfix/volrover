import os
"""Where does the ~200-300 ms/frame Python pose actually go? Break route B's
per-frame update into three phases so we can estimate a pure-C++/cvcGL rewrite:

  cascade  : _world_transforms — the per-module Python for-loop of 4x4 matmuls
             (the wind cascade). Python-call-bound; C++ does this natively.
  instances: _tree_instances — batched numpy matmul + quaternion (already
             vectorised, near-C, but per-tree dispatch + temporaries).
  marshal  : concatenate + numpy_to_vtk into the glyph input arrays (Python->VTK
             copy; a C++ node writes straight into the VTK/cvcGL buffer).
"""
import sys
import time

import numpy as np
import vtk
from vtkmodules.util.numpy_support import numpy_to_vtk

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import forest_geom as fg

N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
FR = 25
forest = fg.build_forest(N, seed=20260817)
nseg = sum(t.n_seg for t in forest)
nleaf = sum(t.n_leaf for t in forest)


def timed_update(t, acc):
    wp, wq, ws, np_, nq, ns = [], [], [], [], [], []
    for tree in forest:
        a = time.perf_counter()
        world = fg._world_transforms(tree, t)
        b = time.perf_counter()
        inst = fg._tree_instances(tree, world)
        c = time.perf_counter()
        p, q, s = inst["seg"]; wp.append(p); wq.append(q); ws.append(s)
        p, q, s = inst["leaf"]; np_.append(p); nq.append(q); ns.append(s)
        acc[0] += b - a
        acc[1] += c - b
    a = time.perf_counter()
    wood = (np.concatenate(wp), np.concatenate(wq), np.concatenate(ws))
    needle = (np.concatenate(np_), np.concatenate(nq), np.concatenate(ns))
    # marshal into VTK arrays (what the glyph mapper consumes)
    for pos, quat, scale in (wood, needle):
        vp = vtk.vtkPoints(); vp.SetData(numpy_to_vtk(np.ascontiguousarray(pos, np.float32), deep=1))
        numpy_to_vtk(np.ascontiguousarray(quat, np.float32), deep=1)
        numpy_to_vtk(np.ascontiguousarray(scale, np.float32), deep=1)
    acc[2] += time.perf_counter() - a


acc = [0.0, 0.0, 0.0]
for i in range(FR):
    timed_update(1.0 + 0.05 * i, acc)
cascade, instances, marshal = [x * 1e3 / FR for x in acc]
total = cascade + instances + marshal
nmod = sum(len(t.mods) for t in forest)
print("route B pose, %d trees (%d modules, %d wood seg, %d needle):" % (N, nmod, nseg, nleaf))
print("  cascade  (Python module loop) : %6.1f ms  (%4.1f%%)" % (cascade, 100 * cascade / total))
print("  instances(numpy matmul+quat)  : %6.1f ms  (%4.1f%%)" % (instances, 100 * instances / total))
print("  marshal  (concat + numpy_to_vtk): %6.1f ms  (%4.1f%%)" % (marshal, 100 * marshal / total))
print("  TOTAL                          : %6.1f ms" % total)
print("  per-module cascade cost        : %.1f us" % (cascade * 1000 / nmod))
