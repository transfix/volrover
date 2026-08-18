import os
"""Numeric correctness: reconstruct route-B instanced geometry from (pos,quat,
scale) and confirm it matches route-C's merged per-vertex buffer, segment for
segment, at a swayed time t. If this passes, the quaternion/scale/compose math is
right and any later visual difference is shading, not geometry."""
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import forest_geom as fg
from wood_shader import _unit_cyl_arrays

U, _ = _unit_cyl_arrays(fg.BASE_TRI)          # unit _CYL points (base y0, radius 1)
NRING = fg._NRING
STRIDE = fg.NEEDLES + 1
UN = np.zeros((STRIDE, 3))                     # unit needle star (radius 1, height 1)
UN[1:, 0::2] = NRING
UN[1:, 1] = 1.0


def R_from_quat(q):
    w, x, y, z = q
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
        [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
        [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])


forest = fg.build_forest(3, seed=1)
t = 1.234  # a swayed moment
worst_w = worst_n = 0.0
for ti, tree in enumerate(forest):
    world = fg._world_transforms(tree, t)
    # route C merged buffers
    wb, nb = fg.repose_route_c(tree, t)
    # route B instances
    inst = fg._tree_instances(tree, world)
    # wood: reconstruct each segment and compare to the matching _CYL block
    wp, wq, wscale = inst["seg"]
    for s in range(tree.n_seg):
        recon = wp[s] + (R_from_quat(wq[s]) @ (U * wscale[s]).T).T
        ref = wb[s * fg._CYL_V:(s + 1) * fg._CYL_V]
        worst_w = max(worst_w, float(np.abs(recon - ref).max()))
    # needles
    lp, lq, lscale = inst["leaf"]
    for s in range(tree.n_leaf):
        recon = lp[s] + (R_from_quat(lq[s]) @ (UN * lscale[s]).T).T
        ref = nb[s * STRIDE:(s + 1) * STRIDE]
        worst_n = max(worst_n, float(np.abs(recon - ref).max()))
    print("tree %d: %d segs, %d leaves, %d modules" %
          (ti, tree.n_seg, tree.n_leaf, len(tree.mods)))

print()
print("max |route B - route C| wood  =", worst_w)
print("max |route B - route C| needle=", worst_n)
print("VERDICT:", "MATCH" if worst_w < 1e-6 and worst_n < 1e-6 else "MISMATCH")
