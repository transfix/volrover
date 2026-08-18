"""Self-contained reproduction of lsystem_forest.py's tree grammar, emitting
either route C (one merged wood+needle actor per tree, per-vertex colour) or
route B (two global vtkGlyph3DMapper actors: unit-cylinder instanced over every
wood segment, unit-needle-star instanced over every leaf).

The tree math (grammar, expansion, module cascade, sway) is copied verbatim from
scripts/examples/lsystem_forest.py so the geometry is representative on every axis
that matters for the render comparison: module hierarchy, segment/needle counts,
per-axis-scaled segments, wind that re-poses the trunk. Terrain/volumes/sun are
omitted — this harness is about the two ways to draw the trees.
"""
import math
import random

import numpy as np

# ── constants (verbatim from lsystem_forest.py) ──────────────────────────────
TREE_RULES = (
    "FF[RL1][RR2][RRR3]F[RL3][RR1][RRR2]RFLR0",
    "FL[T[RF]2]R[TRFL]RTFL4",
    "FL[TRF3]RFLRTFL2",
    "FL[TFL2RFL]R[T[RFLF3]]RTFL2",
    "FL[TRFL4]RFLRTFL4",
)
YROTATE, TILT = 10.0, 120.0
MICRO_TILT = 1.0e-4
T_SCALE, T_RADSCALE = 0.9, 0.6
T_LENGTH, T_RADIUS = 5.0, 0.7
BASE_TRI, NEEDLES = 5, 9
LEAF_LEN, LEAF_RAD = 4.0, 1.0
MATURITY = (1, 2, 2, 3, 3, 3, 4)
TREE_SIZE = (0.32, 0.75)
SWAY_LEVELS = 2
C_WOOD_LIGHT = (0.6549, 0.4901, 0.2392)
C_WOOD_DARK = (0.3607, 0.2510, 0.2000)
C_NEEDLE = (0.1373, 0.5568, 0.1373)
_TREE_AXIS = (0.0, 1.0, 0.0)


def mat_rotate(angle, x, y, z):
    c, s = math.cos(angle), math.sin(angle)
    k = 1.0 - c
    return np.array([
        [c + k * x * x, k * x * y - s * z, k * x * z + s * y, 0.0],
        [k * x * y + s * z, c + k * y * y, k * y * z - s * x, 0.0],
        [k * x * z - s * y, k * y * z + s * x, c + k * z * z, 0.0],
        [0.0, 0.0, 0.0, 1.0]])


def mat_translate(x, y, z):
    m = np.identity(4)
    m[:3, 3] = (x, y, z)
    return m


def xform(m, pts):
    return pts @ m[:3, :3].T + m[:3, 3]


TREE_UP = mat_rotate(math.pi / 2.0, 1.0, 0.0, 0.0)

# unit cylinder (_CYL) local layout, RADIUS/HEIGHT applied per segment below
_RING = np.column_stack((np.cos(np.arange(BASE_TRI) * 2 * math.pi / BASE_TRI),
                         np.sin(np.arange(BASE_TRI) * 2 * math.pi / BASE_TRI)))
_NRING = np.column_stack((np.cos(np.arange(NEEDLES) * 2 * math.pi / NEEDLES),
                          np.sin(np.arange(NEEDLES) * 2 * math.pi / NEEDLES)))
_CYL_V = 2 * BASE_TRI + 2
_CYL_TRIS = []
for _i in range(BASE_TRI):
    _b0, _b1 = 1 + _i, 1 + (_i + 1) % BASE_TRI
    _t0, _t1 = BASE_TRI + 2 + _i, BASE_TRI + 2 + (_i + 1) % BASE_TRI
    _CYL_TRIS += [0, _b0, _b1, BASE_TRI + 1, _t1, _t0, _b0, _t1, _b1, _b0, _t0, _t1]
_CYL_TRIS = np.array(_CYL_TRIS)
_CYL_COLORS = np.array([C_WOOD_LIGHT] + [C_WOOD_DARK] * BASE_TRI
                       + [C_WOOD_LIGHT] + [C_WOOD_DARK] * BASE_TRI)
_TURN_MICRO = mat_rotate(TILT * MICRO_TILT, 0.0, 0.0, 1.0)
_TURN_TILT = mat_rotate(TILT, 0.0, 0.0, 1.0)
_TURN_ROLL = mat_rotate(YROTATE, 0.0, 1.0, 0.0)


class TreeModule(object):
    __slots__ = ("name", "parent", "level", "hang", "segs", "leaves")

    def __init__(self, name, parent, level, hang):
        self.name, self.parent, self.level, self.hang = name, parent, level, hang
        self.segs, self.leaves = [], []


def expand_tree(rule, depth, scale, radscale, name, parent, level, out):
    mod = TreeModule(name, parent, level, np.identity(4))
    out.append(mod)
    cur = np.identity(4)
    stack = []
    seg_len, seg_rad = T_LENGTH * scale, T_RADIUS * radscale
    step = mat_translate(0.0, seg_len, 0.0)
    for ch in rule:
        if ch == "F":
            cur = cur @ _TURN_MICRO
            mod.segs.append((cur, seg_len, seg_rad))
            cur = cur @ step
        elif ch == "[":
            stack.append(cur)
        elif ch == "]":
            cur = stack.pop()
        elif ch == "L":
            mod.leaves.append((cur, scale))
        elif ch == "R":
            cur = cur @ _TURN_ROLL
        elif ch == "T":
            cur = cur @ _TURN_TILT
        elif ch.isdigit() and depth > 1:
            child = expand_tree(TREE_RULES[int(ch)], depth - 1, scale * T_SCALE,
                                radscale * T_RADSCALE, "%s_%d" % (name, len(out)),
                                name, level + 1, out)
            child.hang = cur.copy()
    return mod


# ── the tree/module record used by both routes ───────────────────────────────
class Mod(object):
    __slots__ = ("parent", "hang", "swayer",
                 # route C local geometry
                 "local_wood", "w_off", "w_n", "local_needle", "n_off", "n_n",
                 # route B per-instance (segment transforms + leaf transforms, local)
                 "seg_M", "seg_scale", "leaf_M", "leaf_scale")


class Tree(object):
    __slots__ = ("mods", "px", "py", "pz", "phase", "sway",
                 # route C
                 "wood_verts0", "wood_tris", "wood_colors", "needle_verts0", "needle_lines",
                 "wood_buf", "needle_buf",
                 # route B: per-segment/leaf local matrices flattened for the tree
                 "seg_mod", "seg_local", "seg_scale_arr", "n_seg",
                 "leaf_mod", "leaf_local", "leaf_scale_arr", "n_leaf")


def _module_local_arrays(segs, leaves):
    """Route-C local wood+needle geometry for one module (verbatim math)."""
    wpts = np.empty((len(segs) * _CYL_V, 3))
    local = np.zeros((_CYL_V, 3))
    for s, (m, height, radius) in enumerate(segs):
        local[1:BASE_TRI + 1, 0::2] = _RING * radius
        local[BASE_TRI + 2:, 0::2] = _RING * radius
        local[BASE_TRI + 1:, 1] = height
        wpts[s * _CYL_V:(s + 1) * _CYL_V] = xform(m, local)
    wtris = (_CYL_TRIS + (np.arange(len(segs)) * _CYL_V)[:, None]).ravel()
    wcol = np.tile(_CYL_COLORS, (len(segs), 1))
    stride = NEEDLES + 1
    npts = np.empty((len(leaves) * stride, 3))
    nloc = np.zeros((stride, 3))
    for c, (m, sc) in enumerate(leaves):
        nloc[1:, 0::2] = _NRING * (LEAF_RAD * sc)
        nloc[1:, 1] = LEAF_LEN * sc
        npts[c * stride:(c + 1) * stride] = xform(m, nloc)
    root = (np.arange(len(leaves)) * stride)[:, None]
    nlines = np.empty((len(leaves), NEEDLES, 2), dtype=np.int64)
    nlines[:, :, 0] = root
    nlines[:, :, 1] = root + 1 + np.arange(NEEDLES)
    return wpts, wtris, wcol, npts, nlines.reshape(-1, 2)


def build_tree(px, py, pz, mods):
    """Prepare one tree for BOTH routes: route-C merged buffers + bind-pose verts,
    and route-B per-segment/leaf local transforms (composed with module world at
    pose time)."""
    name_idx = {m.name: i for i, m in enumerate(mods)}
    tree = Tree()
    tree.px, tree.py, tree.pz = px, py, pz
    tmods = []
    wpts_all, wtris_all, wcol_all, npts_all, nlines_all = [], [], [], [], []
    world = [None] * len(mods)
    w_cur = n_cur = 0
    for i, m in enumerate(mods):
        wpts, wtris, wcol, npts, nlines = _module_local_arrays(m.segs, m.leaves)
        hang = (mat_translate(px, py, pz) @ TREE_UP) if m.parent is None else m.hang
        parent = -1 if m.parent is None else name_idx[m.parent]
        W = hang if parent < 0 else world[parent] @ hang
        world[i] = W
        R, tt = W[:3, :3].T, W[:3, 3]
        wpts_all.append(wpts @ R + tt)
        wtris_all.append(wtris + w_cur)
        wcol_all.append(wcol)
        n_n = npts.shape[0]
        if n_n:
            npts_all.append(npts @ R + tt)
            nlines_all.append(nlines + n_cur)
        d = Mod()
        d.parent, d.hang, d.swayer = parent, hang, (m.level <= SWAY_LEVELS)
        d.local_wood, d.w_off, d.w_n = wpts, w_cur, wpts.shape[0]
        d.local_needle, d.n_off, d.n_n = (npts if n_n else None), n_cur, n_n
        # route B: local segment transforms (module-local) + scale, leaf transforms + scale
        d.seg_M = [m.segs[s][0] for s in range(len(m.segs))]
        d.seg_scale = [(m.segs[s][2], m.segs[s][1], m.segs[s][2]) for s in range(len(m.segs))]
        d.leaf_M = [lv[0] for lv in m.leaves]
        d.leaf_scale = [(LEAF_RAD * lv[1], LEAF_LEN * lv[1], LEAF_RAD * lv[1]) for lv in m.leaves]
        tmods.append(d)
        w_cur += wpts.shape[0]
        n_cur += n_n
    tree.mods = tmods
    # route B: flatten per-segment / per-leaf local matrices + scales for the tree
    seg_mod, seg_local, seg_scale = [], [], []
    leaf_mod, leaf_local, leaf_scale = [], [], []
    for i, d in enumerate(tmods):
        for M, sc in zip(d.seg_M, d.seg_scale):
            seg_mod.append(i); seg_local.append(M); seg_scale.append(sc)
        for M, sc in zip(d.leaf_M, d.leaf_scale):
            leaf_mod.append(i); leaf_local.append(M); leaf_scale.append(sc)
    tree.seg_mod = np.array(seg_mod, dtype=np.int64)
    tree.seg_local = np.array(seg_local) if seg_local else np.zeros((0, 4, 4))
    tree.seg_scale_arr = np.array(seg_scale) if seg_scale else np.zeros((0, 3))
    tree.n_seg = len(seg_mod)
    tree.leaf_mod = np.array(leaf_mod, dtype=np.int64)
    tree.leaf_local = np.array(leaf_local) if leaf_local else np.zeros((0, 4, 4))
    tree.leaf_scale_arr = np.array(leaf_scale) if leaf_scale else np.zeros((0, 3))
    tree.n_leaf = len(leaf_mod)
    tree.wood_verts0 = np.concatenate(wpts_all)
    tree.wood_tris = np.concatenate(wtris_all)
    tree.wood_colors = np.concatenate(wcol_all)
    tree.wood_buf = tree.wood_verts0.copy()
    if npts_all:
        tree.needle_verts0 = np.concatenate(npts_all)
        tree.needle_lines = np.concatenate(nlines_all)
        tree.needle_buf = tree.needle_verts0.copy()
    else:
        tree.needle_verts0 = np.zeros((0, 3))
        tree.needle_lines = np.zeros((0, 2), dtype=np.int64)
        tree.needle_buf = None
    return tree


def build_forest(n_trees, seed=20260817, spacing=26.0):
    """A grid of trees with the forest's per-tree maturity/size/phase/sway RNG."""
    random.seed(seed)
    cols = int(math.ceil(math.sqrt(n_trees)))
    forest = []
    for n in range(n_trees):
        gx = (n % cols - cols / 2.0) * spacing
        gy = (n // cols - cols / 2.0) * spacing
        size = random.uniform(*TREE_SIZE)
        mods = []
        expand_tree(TREE_RULES[0], random.choice(MATURITY), size, size,
                    "ftree%d" % n, None, 1, mods)
        ph = random.uniform(0.0, 2 * math.pi)
        sw = 0.020 + 0.016 * random.random()
        tree = build_tree(gx, gy, 0.0, mods)
        tree.phase, tree.sway = ph, sw
        forest.append(tree)
    return forest


# ── per-frame posing ─────────────────────────────────────────────────────────
def _world_transforms(tree, t):
    """The wind cascade: each module's world transform, sway baked into swayers."""
    a = tree.sway * math.sin(1.3 * t + tree.phase)
    sway = mat_rotate(a, *_TREE_AXIS)
    world = [None] * len(tree.mods)
    for i, m in enumerate(tree.mods):
        local = (m.hang @ sway) if m.swayer else m.hang
        world[i] = local if m.parent < 0 else world[m.parent] @ local
    return world


def repose_route_c(tree, t):
    """Route C: blit posed local verts into the merged buffers (== repose_tree)."""
    world = _world_transforms(tree, t)
    wb, nb = tree.wood_buf, tree.needle_buf
    for m, W in zip(tree.mods, world):
        R, tt = W[:3, :3].T, W[:3, 3]
        wb[m.w_off:m.w_off + m.w_n] = m.local_wood @ R + tt
        if m.n_n:
            nb[m.n_off:m.n_off + m.n_n] = m.local_needle @ R + tt
    return wb, nb


def quats_wxyz_from_Rs(R):
    """Batched quaternion (w,x,y,z) from (N,3,3) rotation matrices (column-vector
    convention). Vectorised branchless form of the standard trace method."""
    n = R.shape[0]
    q = np.empty((n, 4))
    m00, m11, m22 = R[:, 0, 0], R[:, 1, 1], R[:, 2, 2]
    tr = m00 + m11 + m22
    # case 0: tr > 0
    c0 = tr > 0
    # case 1: m00 largest; case 2: m11 largest; case 3: m22 largest
    c1 = (~c0) & (m00 >= m11) & (m00 >= m22)
    c2 = (~c0) & (~c1) & (m11 >= m22)
    c3 = (~c0) & (~c1) & (~c2)
    s = np.empty(n)
    s[c0] = np.sqrt(tr[c0] + 1.0) * 2
    q[c0, 0] = 0.25 * s[c0]
    q[c0, 1] = (R[c0, 2, 1] - R[c0, 1, 2]) / s[c0]
    q[c0, 2] = (R[c0, 0, 2] - R[c0, 2, 0]) / s[c0]
    q[c0, 3] = (R[c0, 1, 0] - R[c0, 0, 1]) / s[c0]
    s[c1] = np.sqrt(1.0 + m00[c1] - m11[c1] - m22[c1]) * 2
    q[c1, 0] = (R[c1, 2, 1] - R[c1, 1, 2]) / s[c1]
    q[c1, 1] = 0.25 * s[c1]
    q[c1, 2] = (R[c1, 0, 1] + R[c1, 1, 0]) / s[c1]
    q[c1, 3] = (R[c1, 0, 2] + R[c1, 2, 0]) / s[c1]
    s[c2] = np.sqrt(1.0 + m11[c2] - m00[c2] - m22[c2]) * 2
    q[c2, 0] = (R[c2, 0, 2] - R[c2, 2, 0]) / s[c2]
    q[c2, 1] = (R[c2, 0, 1] + R[c2, 1, 0]) / s[c2]
    q[c2, 2] = 0.25 * s[c2]
    q[c2, 3] = (R[c2, 1, 2] + R[c2, 2, 1]) / s[c2]
    s[c3] = np.sqrt(1.0 + m22[c3] - m00[c3] - m11[c3]) * 2
    q[c3, 0] = (R[c3, 1, 0] - R[c3, 0, 1]) / s[c3]
    q[c3, 1] = (R[c3, 0, 2] + R[c3, 2, 0]) / s[c3]
    q[c3, 2] = (R[c3, 1, 2] + R[c3, 2, 1]) / s[c3]
    q[c3, 3] = 0.25 * s[c3]
    return q


def _tree_instances(tree, world):
    """Per-segment and per-leaf (pos, quat, scale) for one posed tree.

    M = W_module @ local  (numpy matmul == the forest's xform composition), so
    pos = M[:3,3], rotation R = M[:3,:3] on column vectors, scale carried
    separately — exactly the T·R·S the glyph mapper reconstructs per instance."""
    Wm = np.array(world)  # (K,4,4)
    out = {}
    for kind, mod_idx, local, scale, n in (
            ("seg", tree.seg_mod, tree.seg_local, tree.seg_scale_arr, tree.n_seg),
            ("leaf", tree.leaf_mod, tree.leaf_local, tree.leaf_scale_arr, tree.n_leaf)):
        if n == 0:
            out[kind] = (np.zeros((0, 3)), np.zeros((0, 4)), np.zeros((0, 3)))
            continue
        M = Wm[mod_idx] @ local              # (n,4,4)
        pos = M[:, :3, 3]
        quat = quats_wxyz_from_Rs(M[:, :3, :3])
        out[kind] = (pos, quat, scale)
    return out


def route_b_instances(forest, t):
    """Global instance arrays for route B: all wood segments -> one glyph, all
    needle stars -> one glyph. Returns (wood_pos, wood_quat, wood_scale),
    (needle_pos, needle_quat, needle_scale)."""
    wp, wq, ws, np_, nq, ns = [], [], [], [], [], []
    for tree in forest:
        world = _world_transforms(tree, t)
        inst = _tree_instances(tree, world)
        p, q, s = inst["seg"]; wp.append(p); wq.append(q); ws.append(s)
        p, q, s = inst["leaf"]; np_.append(p); nq.append(q); ns.append(s)
    wood = (np.concatenate(wp), np.concatenate(wq), np.concatenate(ws))
    needle = (np.concatenate(np_), np.concatenate(nq), np.concatenate(ns))
    return wood, needle

