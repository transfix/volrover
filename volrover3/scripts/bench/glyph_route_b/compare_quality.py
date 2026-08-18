import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import forest_geom as fg
import render_forest as rf

OUT = os.path.dirname(os.path.abspath(__file__))
N = 12
T = 1.0

forest = fg.build_forest(N, seed=20260817)
nseg = sum(t.n_seg for t in forest)
nleaf = sum(t.n_leaf for t in forest)
print("forest: %d trees, %d wood segments, %d needle stars" % (N, nseg, nleaf))

# route C
cA, cU = rf.build_route_c(forest)
cU(T)
renC, rwC = rf.make_scene(cA, shadows=False)
rf.frame_island_camera(renC, forest)
rwC.Render()
imgC = rf.screenshot(rwC, "%s/cmp_routeC.png" % OUT)
print("route C: %d actors" % len(cA))

# route B (same camera params)
bA, bU = rf.build_route_b(forest)
bU(T)
renB, rwB = rf.make_scene(bA, shadows=False)
# copy route C's camera exactly for a fair image compare
renB.SetActiveCamera(renC.GetActiveCamera())
renB.ResetCameraClippingRange()
rwB.Render()
imgB = rf.screenshot(rwB, "%s/cmp_routeB.png" % OUT)
print("route B: %d actors (2 glyphs)" % len(bA))

# crude image similarity on the foreground (non-sky) pixels
sky = np.array([158, 194, 235])
fgmaskC = np.any(np.abs(imgC.astype(int) - sky) > 30, axis=2)
fgmaskB = np.any(np.abs(imgB.astype(int) - sky) > 30, axis=2)
inter = np.sum(fgmaskC & fgmaskB); union = np.sum(fgmaskC | fgmaskB)
print("foreground coverage: C=%d px  B=%d px  IoU=%.3f" %
      (fgmaskC.sum(), fgmaskB.sum(), inter / max(1, union)))
