# lsystem_forest.py — an L-system island: terrain, forest, sea and sky.
#
# The companion to lsystem_tree.py. That example puts ONE tree in the scene graph
# as a live hierarchy, to show what a scene graph buys you over immediate mode.
# This one is about composition: a second L-system lays out the *world*, and the
# tree grammar is instanced across it at different maturities.
#
# WHAT THE GRAMMARS DO
#   * The terrain L-system walks a turtle over the ground plane. `F` paints a
#     patch of whatever material is current — `K` rock, `D` dirt, `G` grass — and
#     `S` drops a tree seed. Patches do not just colour the ground, they SHAPE it:
#     rock pushes the height up into ridges, dirt scoops it down into hollows. So
#     where the grammar wandered is where the island's relief, its materials and
#     its trees all come from — one walk, three outputs.
#   * Every seed that landed above the waterline sprouts the tree grammar from
#     lsystem_tree.py, at a maturity drawn per tree. A `1` is a seedling, a `4` is
#     a full canopy, so the forest reads as a population rather than a stamp.
#
# THE TREES ARE FULL HIERARCHIES
#   Every L-system module gets its own GeometryNode, exactly as lsystem_tree.py
#   does, so wind ACCUMULATES down each tree and the tips move more than the
#   trunk — the original demo's motion, not a rigid whole-tree lean.
#
#   This was not affordable until recently: posing a node was ~80% state-tree
#   write, and the write came back through handleStateChanged and ran the whole
#   transform cascade a second time, so 70 hierarchical trees cost ~300 ms/frame.
#   With cvc::gl::state_publisher batching those writes off the render path
#   (libcvc #193) the same scene poses in ~24 ms. REQUIRES that fix; against an
#   older cvcGL this example will crawl.
#
# THE TWO VOLUMES
#   * SEA — a cvc::volume over the island's footprint, filling the space between
#     the seabed and a travelling wave surface, so water sits exactly in the
#     hollows the dirt patches carved AND is only as thick as the water actually
#     is. That thickness is what sells it: the transfer function is translucent
#     enough that a knee-deep column barely tints the sand while open water
#     stacks up to solid navy, so the shallows show their bottom and the deeps do
#     not. Both the field and the transfer function animate — the field carries
#     the swell, the transfer function breathes the surface.
#   * SKY — a second volume slab overhead holding a band-limited noise field,
#     scrolled by the wind. Its transfer function pins alpha to exactly zero
#     below a threshold, so empty sky is empty rather than a faint grey box, and
#     the clouds read as discrete puffs.
#
# The script never touches the camera. It sets the world bounds once so volrover3
# frames the island; orbit, pan and zoom stay yours.
#
# LIVE CONTROLS (Python Console dock -> State tab, or from any script):
#     pycvc.state_set(app, "volrover3.lsystem_forest.speed",  "0.4")  # clock scale
#     pycvc.state_set(app, "volrover3.lsystem_forest.waves",  "1")    # 0/1
#     pycvc.state_set(app, "volrover3.lsystem_forest.clouds", "1")    # 0/1
#     pycvc.state_set(app, "volrover3.lsystem_forest.wind",   "1")    # 0/1
#
# HOW TO RUN (inside a running volrover3): Python Console dock -> "Jobs" tab ->
# "Load Script..." -> pick this file. Select it -> Stop to end.
#
# THE CONTRACT: define a module-level step(dt); import-time code runs once at
# submit, step(dt) runs every tick.

import math
import random

import numpy as np

import pycvc
import pycvc_gl
import vrhost

SEED = 20260817  # fixed so the island is reproducible run to run
STATE = "volrover3.lsystem_forest"

# ── the island ───────────────────────────────────────────────────────────────
HALF = 120.0  # terrain spans [-HALF, HALF] in x and y
TERRAIN_N = 96  # heightfield resolution (N x N vertices)
SEA_LEVEL = 0.0
PEAK = 34.0  # height of the central dome before patches
SHELF = -9.0  # offset that drowns the outer ring, making a coastline

# ── the terrain grammar ──────────────────────────────────────────────────────
#   F  step forward, painting a patch of the current material
#   + -  turn        [ ]  push/pop the turtle
#   K D G  set material to rock / dirt / grass
#   S  drop a tree seed here
#   A B C  non-terminals
# Four arms at ~108 deg to each other, so the walk covers the island radially
# instead of arcing across one side of it.
TERRAIN_AXIOM = "G[A][++++A][++++++++A][++++++++++++A]"
TERRAIN_RULES = {
    # A arcs as it recurses (the trailing +), so the main path curls around the
    # island instead of marching straight off the edge of it.
    "A": "FF[+GBS]F[-DC]F+SA",
    "B": "F[+FS]KFFG-B",
    "C": "FD[-FFS]F[+KFS]+C",
}
TERRAIN_DEPTH = 6
TURN = 27.0  # degrees per + / -
STEP0 = 7.0  # first step length, shrinking with depth
STEP_DECAY = 0.86
PATCH_R0 = 18.0  # first patch radius, shrinking likewise
PATCH_DECAY = 0.88

ROCK, DIRT, GRASS = 0, 1, 2
MAT_COLOR = np.array([
    [0.46, 0.45, 0.43],  # rock  — grey
    [0.42, 0.31, 0.20],  # dirt  — brown
    [0.27, 0.44, 0.19],  # grass — green
])
MAT_LIFT = (15.0, -11.0, 2.0)  # rock ridges up, dirt hollows down, grass ~flat
SAND = np.array([0.68, 0.62, 0.44])  # shoreline blend just above the waterline

# ── the tree grammar (verbatim from lsystem_tree.py / the 2004 original) ─────
TREE_RULES = (
    "FF[RL1][RR2][RRR3]F[RL3][RR1][RRR2]RFLR0",
    "FL[T[RF]2]R[TRFL]RTFL4",
    "FL[TRF3]RFLRTFL2",
    "FL[TFL2RFL]R[T[RFLF3]]RTFL2",
    "FL[TRFL4]RFLRTFL4",
)
# Turn amounts as the original *used* them: matrix.c fed these to cos()/sin()
# directly, so they are radians despite being written as degrees. See
# lsystem_tree.py — the tree's whole shape depends on it.
YROTATE, TILT = 10.0, 120.0
MICRO_TILT = 1.0e-4
T_SCALE, T_RADSCALE = 0.9, 0.6
T_LENGTH, T_RADIUS = 5.0, 0.7
BASE_TRI, NEEDLES = 5, 9
LEAF_LEN, LEAF_RAD = 4.0, 1.0
MATURITY = (1, 2, 2, 3, 3, 3, 4)  # drawn per tree — a population, not a stamp
TREE_SIZE = (0.32, 0.75)  # extra per-tree scale range
MAX_TREES = 70  # the grammar yields ~200 dry seeds; plant a sample of them
SWAY_LEVELS = 2  # pose modules this deep; sway accumulates below them for free
C_WOOD_LIGHT = (0.6549, 0.4901, 0.2392)
C_WOOD_DARK = (0.3607, 0.2510, 0.2000)
C_NEEDLE = (0.1373, 0.5568, 0.1373)

# ── the sea volume ───────────────────────────────────────────────────────────
SEA_N = 64  # x/y resolution
SEA_NZ = 20  # z layers
SEA_FLOOR = SEA_LEVEL - 20.0
SEA_TOP = SEA_LEVEL + 5.0
WAVE_AMP = 1.6
WAVE_LEN = 46.0
WAVE_SPEED = 7.0

# ── the cloud grammar ────────────────────────────────────────────────────────
# The clouds were sums of sine octaves, which is why they came out as round
# blobs: every octave is separable and axis-aligned, so its level sets are
# ellipses and no amount of stacking them makes an edge that turns a corner.
# A grammar does, because a turtle can branch. Same symbols as the terrain
# walk, deposited into a 2-D density map instead of onto the ground:
#   F  drift forward, depositing a puff        + -  turn
#   [ ]  push/pop        <  shrink the puff radius       A B  non-terminals
CLOUD_AXIOM = "[A][+++++A][-----A][++++++++++A][----------A][+++++++++++++++A]"
CLOUD_RULES = {
    "A": "FF[+<B]F[-<B]<FA",
    "B": "F[+<F]F<[-<F]B",
}
CLOUD_DEPTH = 6
CLOUD_TURN = 32.0     # degrees per + / -
CLOUD_STEP0 = 6.5     # first drift, in map cells
CLOUD_STEP_DECAY = 0.9
CLOUD_PUFF0 = 4.0     # first puff radius, in map cells
CLOUD_PUFF_DECAY = 0.88
CLOUD_MAPS = 2        # independent skies, crossfaded so cloud EVOLVES
CLOUD_MORPH_S = 26.0  # seconds for a full crossfade cycle

# ── the sky volume ───────────────────────────────────────────────────────────
SKY_N = 64
SKY_NZ = 14
SKY_BASE, SKY_TOP = 82.0, 104.0
CLOUD_DRIFT = 5.0  # world units per second
SKY_HALF = 150.0  # the cloud slab overhangs the island a little

SIM_DT = 1.0 / 60.0  # simulation quantum (see lsystem_tree.py on world_clock)
PRIME_DT = 0.25  # a tick slower than this at startup is scene setup, not lag
STALL_QUANTA = 60  # only shout about a stall once it costs ~half a second


# ── 4x4 transforms (row-major, as GraphicsNode.setTransform takes) ───────────
def mat_rotate(angle, x, y, z):
    c, s = math.cos(angle), math.sin(angle)
    k = 1.0 - c
    return np.array([
        [c + k * x * x, k * x * y - s * z, k * x * z + s * y, 0.0],
        [k * x * y + s * z, c + k * y * y, k * y * z - s * x, 0.0],
        [k * x * z - s * y, k * y * z + s * x, c + k * z * z, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ])


def mat_translate(x, y, z):
    m = np.identity(4)
    m[:3, 3] = (x, y, z)
    return m


def xform(m, pts):
    return pts @ m[:3, :3].T + m[:3, 3]


# The tree turtle stands on +Y; the world is Z-up.
TREE_UP = mat_rotate(math.pi / 2.0, 1.0, 0.0, 0.0)


# ── walk the terrain grammar ─────────────────────────────────────────────────
def walk_terrain():
    """Run the terrain turtle; return (patches, seeds).

    patches: list of (x, y, radius, material) — these both colour AND lift/carve
    seeds:   list of (x, y) tree sites
    """
    # Expanding to a fixed depth would blow up (the rules are all recursive), so
    # the turtle carries its own depth and simply stops recursing — the same
    # bounded-expansion trick drawPlant() used in the original.
    patches, seeds = [], []
    stack = []
    x, y, head, mat, d = 0.0, 0.0, 90.0, GRASS, 0
    step, radius = STEP0, PATCH_R0
    todo = list(TERRAIN_AXIOM)
    guard = 0
    while todo and guard < 20000:
        guard += 1
        c = todo.pop(0)
        if c == "F":
            nx = x + step * math.cos(math.radians(head))
            ny = y + step * math.sin(math.radians(head))
            patches.append((0.5 * (x + nx), 0.5 * (y + ny), radius, mat))
            x, y = nx, ny
        elif c == "+":
            head += TURN
        elif c == "-":
            head -= TURN
        elif c == "[":
            stack.append((x, y, head, mat, step, radius, d))
        elif c == "]":
            if stack:
                x, y, head, mat, step, radius, d = stack.pop()
        elif c == "K":
            mat = ROCK
        elif c == "D":
            mat = DIRT
        elif c == "G":
            mat = GRASS
        elif c == "S":
            seeds.append((x, y))
        elif c in TERRAIN_RULES and d < TERRAIN_DEPTH:
            d += 1
            step *= STEP_DECAY
            radius *= PATCH_DECAY
            todo = list(TERRAIN_RULES[c]) + todo
    return patches, seeds


# ── build the heightfield + per-vertex material colours ──────────────────────
def build_terrain(patches):
    ax = np.linspace(-HALF, HALF, TERRAIN_N)
    gx, gy = np.meshgrid(ax, ax, indexing="xy")  # gx[j,i] = x, gy[j,i] = y

    # Base island: a dome that falls below sea level toward the edges, so there
    # is always a coastline for the sea volume to fill.
    r2 = gx * gx + gy * gy
    h = PEAK * np.exp(-r2 / (0.34 * HALF * HALF)) + SHELF
    h += 4.5 * np.sin(gx * 0.045) * np.cos(gy * 0.041)
    h += 2.2 * np.sin(gx * 0.11 + 1.3) * np.sin(gy * 0.097)  # finer relief

    # Each patch lifts or carves within its radius, with a smooth falloff, and
    # stains the ground its material colour by the same weight.
    wsum = np.full(gx.shape, 0.35)  # a little default weight keeps bare ground grass
    csum = MAT_COLOR[GRASS] * 0.35
    csum = np.repeat(csum[None, None, :], gx.shape[0], 0).repeat(gx.shape[1], 1)
    lift = np.zeros(gx.shape)
    for px, py, pr, mat in patches:
        d2 = (gx - px) ** 2 + (gy - py) ** 2
        w = np.clip(1.0 - d2 / (pr * pr), 0.0, 1.0) ** 2
        lift += MAT_LIFT[mat] * w
        wsum += w
        csum += MAT_COLOR[mat] * w[:, :, None]

    # Both the relief and the colour are weighted AVERAGES over the patches
    # covering each point, not sums. Summing the lift makes overlapping patches
    # stack into needle-thin spires wherever the walk crossed itself; averaging
    # keeps the relief bounded by MAT_LIFT and reads as terrain.
    h += lift / wsum
    color = csum / wsum[:, :, None]
    # Blend to sand right at the waterline so the shore reads as a beach.
    shore = np.clip(1.0 - np.abs(h - SEA_LEVEL) / 4.5, 0.0, 1.0)[:, :, None]
    color = color * (1.0 - shore) + SAND * shore
    return gx, gy, h, color


def terrain_mesh(app, gx, gy, h, color):
    n = TERRAIN_N
    pts = np.stack((gx.ravel(), gy.ravel(), h.ravel()), axis=1)
    i = np.arange(n - 1)
    j = np.arange(n - 1)
    jj, ii = np.meshgrid(j, i, indexing="ij")
    v = (jj * n + ii).ravel()
    tris = np.column_stack((v, v + 1, v + n, v + 1, v + n + 1, v + n)).ravel()
    g = pycvc.geometry(app)
    g.add_vertices(pts.ravel().tolist())
    g.add_triangles(tris.tolist())
    g.set_colors(color.reshape(-1, 3).ravel().tolist())
    return g


def height_at(h, x, y):
    """Sample the heightfield at a world (x, y) by nearest grid vertex."""
    i = int(round((x + HALF) / (2 * HALF) * (TERRAIN_N - 1)))
    j = int(round((y + HALF) / (2 * HALF) * (TERRAIN_N - 1)))
    i = min(max(i, 0), TERRAIN_N - 1)
    j = min(max(j, 0), TERRAIN_N - 1)
    return float(h[j, i])


# ── bake one tree (the whole grammar flattened into two meshes) ──────────────
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
    """One rule expansion: its own geometry, and where it hangs off its parent."""

    __slots__ = ("name", "parent", "level", "hang", "segs", "leaves", "node", "needles",
                 "phase", "sway")

    def __init__(self, name, parent, level, hang):
        self.name, self.parent, self.level, self.hang = name, parent, level, hang
        self.segs, self.leaves, self.node, self.needles = [], [], None, None
        self.phase = self.sway = 0.0


def expand_tree(rule, depth, scale, radscale, name, parent, level, out):
    """Walk one rule; the module keeps only ITS OWN segments, in its local frame.

    Children record the turtle pose where they attach, which becomes their node
    transform — so the graph composes what a flattened bake would have baked in,
    and moving a module moves everything below it.
    """
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


def module_meshes(app, segs, leaves):

    pts = np.empty((len(segs) * _CYL_V, 3))
    local = np.zeros((_CYL_V, 3))
    for s, (m, height, radius) in enumerate(segs):
        local[1:BASE_TRI + 1, 0::2] = _RING * radius
        local[BASE_TRI + 2:, 0::2] = _RING * radius
        local[BASE_TRI + 1:, 1] = height
        pts[s * _CYL_V:(s + 1) * _CYL_V] = xform(m, local)
    tris = (_CYL_TRIS + (np.arange(len(segs)) * _CYL_V)[:, None]).ravel()
    wood = pycvc.geometry(app)
    wood.add_vertices(pts.ravel().tolist())
    wood.add_triangles(tris.tolist())
    wood.set_colors(np.tile(_CYL_COLORS, (len(segs), 1)).ravel().tolist())

    stride = NEEDLES + 1
    npts = np.empty((len(leaves) * stride, 3))
    nloc = np.zeros((stride, 3))
    for c, (m, sc) in enumerate(leaves):
        nloc[1:, 0::2] = _NRING * (LEAF_RAD * sc)
        nloc[1:, 1] = LEAF_LEN * sc
        npts[c * stride:(c + 1) * stride] = xform(m, nloc)
    root = (np.arange(len(leaves)) * stride)[:, None]
    lines = np.empty((len(leaves), NEEDLES, 2), dtype=np.int64)
    lines[:, :, 0] = root
    lines[:, :, 1] = root + 1 + np.arange(NEEDLES)
    needles = pycvc.geometry(app)
    needles.add_vertices(npts.ravel().tolist())
    needles.add_lines(lines.ravel().tolist())
    return wood, needles


# ── build the scene ──────────────────────────────────────────────────────────
_app = vrhost.app()
_sg = vrhost.scene()
random.seed(SEED)

_patches, _seeds = walk_terrain()
print("lsystem_forest: terrain grammar -> %d patches, %d seeds." % (len(_patches), len(_seeds)))

_gx, _gy, _H, _C = build_terrain(_patches)
_sg.addGraphics("forest_terrain", terrain_mesh(_app, _gx, _gy, _H, _C))
_terrain_node = _sg.geometry_node("forest_terrain")
_terrain_node.setUseSingleColor(False)

# Plant a tree at every seed that came down on dry land, each at its own maturity.
_dry = [(x, y) for (x, y) in _seeds if height_at(_H, x, y) >= SEA_LEVEL + 1.0]
if len(_dry) > MAX_TREES:
    # Node count is the frame budget, so thin the forest rather than plant every
    # seed the grammar dropped. Sampled, not truncated, to keep it spread out.
    _dry = random.sample(_dry, MAX_TREES)
print("lsystem_forest: %d seeds, %d on dry land, planting %d." %
      (len(_seeds), sum(1 for x, y in _seeds if height_at(_H, x, y) >= SEA_LEVEL + 1.0), len(_dry)))

_trees = []      # every module of every tree, flat, for the per-frame pose
_tree_roots = []  # the root module of each tree
for _n, (_sx, _sy) in enumerate(_dry):
    _hz = height_at(_H, _sx, _sy)
    _size = random.uniform(*TREE_SIZE)
    _mods = []
    expand_tree(TREE_RULES[0], random.choice(MATURITY), _size, _size,
                "ftree%d" % _n, None, 1, _mods)
    # One phase per TREE, not per module: the modules of a tree must lean
    # together or it reads as a bush in a blender rather than a tree in wind.
    _ph = random.uniform(0.0, 2 * math.pi)
    _sw = 0.010 + 0.008 * random.random()
    for _m in _mods:
        if not _m.segs:
            continue
        _wood, _needles = module_meshes(_app, _m.segs, _m.leaves)
        if _m.parent is None:
            _sg.addGraphics(_m.name, _wood)
            _m.node = _sg.geometry_node(_m.name)
            # The root carries the tree out to its seed and stands it upright;
            # every module below is expressed in its parent's frame.
            _m.hang = mat_translate(_sx, _sy, _hz) @ TREE_UP
            _tree_roots.append(_m)
        else:
            _m.node = _sg.add_child_geometry(_m.parent, _m.name, _wood)
        _m.node.setUseSingleColor(False)
        _m.phase, _m.sway = _ph, _sw
        _m.node.setTransform(_m.hang.ravel().tolist())
        if _m.leaves:
            _m.needles = _sg.add_child_geometry(_m.name, _m.name + "_n", _needles)
            _m.needles.setRenderMode(pycvc_gl.GeometryRenderMode_LINES)
            _m.needles.setUseSingleColor(True)
            _m.needles.setColor(*C_NEEDLE)
        _trees.append(_m)

# Only the top SWAY_LEVELS are re-posed each frame. Everything below inherits
# the motion through the graph, which is the whole point of the hierarchy — and
# it keeps the per-frame cost proportional to the trunks, not to the twigs.
_swayers = [m for m in _trees if m.level <= SWAY_LEVELS]
print("lsystem_forest: %d trees planted as %d modules (%d posed/frame), maturities %s." %
      (len(_tree_roots), len(_trees), len(_swayers), sorted(set(MATURITY))))


# ── the sea: a volume whose field is depth under a travelling wave ───────────
_sea_ax = np.linspace(-HALF, HALF, SEA_N)
_sea_z = np.linspace(SEA_FLOOR, SEA_TOP, SEA_NZ)
_sgx, _sgy = np.meshgrid(_sea_ax, _sea_ax, indexing="xy")
# Terrain height sampled on the sea grid, so the water can be masked off wherever
# the ground stands above it — water shows up only in the hollows.
_idx = np.clip(((_sea_ax + HALF) / (2 * HALF) * (TERRAIN_N - 1)).round().astype(int),
               0, TERRAIN_N - 1)
_sea_terrain = _H[np.ix_(_idx, _idx)]

_sea_vol = pycvc.volume(_app)
# Phase of the swell at each column, precomputed — only the time term changes.
_wave_phase = (2 * math.pi / WAVE_LEN) * (_sgx + 0.6 * _sgy)
_zcol = _sea_z[:, None, None].astype(np.float32)


def sea_field(t):
    """Depth below the wave surface, bounded by the seabed — the sea's field.

    The seabed clip is what makes shallow water actually look shallow. Masking
    per COLUMN (is there sea floor here at all?) is not enough: without the
    per-voxel bound, every wet column is filled from SEA_FLOOR up, so a sandbar
    under a foot of water carries the same 20-unit slab of water as open ocean
    and no transfer function can make the bottom show through it. Water only
    exists between the bed and the surface.
    """
    surf = SEA_LEVEL + WAVE_AMP * (np.sin(_wave_phase - WAVE_SPEED * t * 0.1)
                                   + 0.45 * np.sin(1.7 * _wave_phase + WAVE_SPEED * t * 0.13))
    below_surface = surf[None, :, :] - _zcol
    above_bed = _zcol - _sea_terrain[None, :, :]
    wet = (below_surface > 0.0) & (above_bed > 0.0)
    depth = np.clip(below_surface / 6.0, 0.0, 1.0)
    return (depth * wet).astype(np.float32)


def sea_transfer(t):
    """Transfer function for the sea: translucent in the shallows, opaque offshore.

    The scalar is depth below the surface, so 0 is the waterline and 1 is 6 m
    down. Colour therefore runs LIGHT at 0 and dark at 1 — the reverse looks
    like an X-ray of the sea.

    The alphas are in thousandths for the same reason the sky's are: VTK applies
    the opacity function once per ScalarOpacityUnitDistance, which it derives
    from the voxel spacing (~0.06 world units here), so a ray crossing 6 m of
    water applies it ~100 times. At the 0.3-per-sample this used to carry, a
    puddle was as opaque as the deep ocean. At ~0.01 a knee-deep column
    accumulates to roughly 0.15 and the sand reads straight through it, while
    open water still stacks to effectively solid.

    `k` drifts the mid alpha slightly so the surface keeps breathing between
    field updates.
    """
    k = 0.0100 + 0.0015 * math.sin(t * 0.9)
    color = [0.00, 0.42, 0.78, 0.74,   # waterline: pale turquoise over the sand
             0.25, 0.14, 0.55, 0.66,
             0.60, 0.04, 0.26, 0.46,
             1.00, 0.01, 0.09, 0.22]   # deep water: near-black navy
    opacity = [0.00, 0.0, 0.12, k * 0.45, 0.55, k, 1.00, k * 2.0]
    return color, opacity


# Load the real field BEFORE the node ever sees the volume. A VolumeNode built
# over an all-zero volume computes its scalar range as [0, 0] and installs a
# default transfer function against that degenerate range — after which you can
# update the voxels all you like and it still renders as one solid block.
_sea_vol.set_float_grid(sea_field(0.0).ravel().tolist(), SEA_N, SEA_N, SEA_NZ,
                        -HALF, -HALF, SEA_FLOOR, HALF, HALF, SEA_TOP)
_sea_grid = _sea_vol.grid()  # zero-copy (nz, ny, nx) view for the per-frame writes
_sg.addGraphics("forest_sea", _sea_vol)
_sea_node = _sg.volume_node("forest_sea")
_sea_node.setTransferFunction(*sea_transfer(0.0))

# ── the sky: a noise slab, scrolled by the wind ──────────────────────────────

# Band-limited noise: a few octaves of sine in x/y, tapered top and bottom so the
# slab has soft faces rather than a hard cut.
def walk_clouds(rng):
    """Run the cloud turtle; return an (SKY_N, SKY_N) density map.

    Deposits a soft radial puff per F. Because the turtle BRANCHES, the union of
    those puffs has concave, ragged outline — the thing summed sine octaves
    cannot produce. The map wraps in x so it can scroll forever without a seam.
    """
    field = np.zeros((SKY_N, SKY_N), dtype=np.float32)
    gy, gx = np.mgrid[0:SKY_N, 0:SKY_N].astype(np.float32)

    x, y = rng.uniform(0.25, 0.75) * SKY_N, rng.uniform(0.3, 0.7) * SKY_N
    head = rng.uniform(0.0, 360.0)
    step, puff, depth = CLOUD_STEP0, CLOUD_PUFF0, 0
    stack = []
    todo = list(CLOUD_AXIOM)
    guard = 0
    while todo and guard < 6000:
        guard += 1
        c = todo.pop(0)
        if c == "F":
            x += step * math.cos(math.radians(head))
            y += step * math.sin(math.radians(head))
            # Wrap in x (the scroll axis) and clamp in y.
            x %= SKY_N
            y = min(max(y, 0.0), float(SKY_N - 1))
            # Nearest image in x, so a puff near the seam deposits on both sides.
            dx = np.abs(gx - x)
            dx = np.minimum(dx, SKY_N - dx)
            d2 = dx * dx + (gy - y) ** 2
            field += np.exp(-d2 / (2.0 * puff * puff)).astype(np.float32)
        elif c == "+":
            head += CLOUD_TURN
        elif c == "-":
            head -= CLOUD_TURN
        elif c == "<":
            puff *= CLOUD_PUFF_DECAY
        elif c == "[":
            stack.append((x, y, head, step, puff, depth))
        elif c == "]":
            if stack:
                x, y, head, step, puff, depth = stack.pop()
        elif c in CLOUD_RULES and depth < CLOUD_DEPTH:
            depth += 1
            step *= CLOUD_STEP_DECAY
            todo = list(CLOUD_RULES[c]) + todo
    # Normalise on a high PERCENTILE, not the max: a single spot where several
    # branches overlap would otherwise set the scale and push the rest of the
    # sky under the threshold, leaving two lonely puffs.
    m = float(np.percentile(field, 99.0))
    return np.clip(field / m, 0.0, 1.0).astype(np.float32) if m > 0 else field


_rng = np.random.default_rng(SEED)
# Several independent skies. One would only ever translate; crossfading between
# them is what makes the cloud EVOLVE — puffs grow and dissolve in place rather
# than sliding past like a painted backdrop.
_cloud_maps = [walk_clouds(_rng) for _ in range(CLOUD_MAPS)]
_taper = np.sin(np.linspace(0, math.pi, SKY_NZ)).astype(np.float32)[:, None, None]
# The z faces are tapered above, but the x/y faces were a hard cut: cloud density
# ran right up to the slab boundary, so the volume ended in a straight vertical
# wall hanging in the sky. Fade to zero at every face instead.
_w = np.hanning(SKY_N).astype(np.float32) ** 0.5
_edge_fade = np.minimum.outer(_w, _w)[None, :, :]

# A high threshold is what makes this read as CLOUD rather than as overcast: only
# the top third of the noise becomes anything at all, so the slab is mostly holes
# and you can see the sky (and the island) through the gaps.
CLOUD_FLOOR = 0.22


def _sky_raw(shift, morph):
    """Density at a CONTINUOUS scroll offset and crossfade position.

    Two things make this fluid rather than steppy. The scroll is sub-cell: the
    slab used to jump a whole column at a time (np.roll takes an integer), which
    is exactly the snapping — here adjacent integer shifts are blended by the
    fractional part, so the drift is smooth at any speed. And the map itself is a
    crossfade between independent grammar-grown skies, so the cloud changes SHAPE
    as it travels instead of being a rigid pattern sliding past.
    """
    i = int(math.floor(morph)) % CLOUD_MAPS
    a = _cloud_maps[i]
    b = _cloud_maps[(i + 1) % CLOUD_MAPS]
    # Smoothstep the blend so the crossfade has no visible kick at either end.
    u = morph - math.floor(morph)
    u = u * u * (3.0 - 2.0 * u)
    base2d = (1.0 - u) * a + u * b

    c = int(math.floor(shift))
    f = shift - c
    base = (1.0 - f) * np.roll(base2d, c, axis=1) + f * np.roll(base2d, c + 1, axis=1)

    lumps = np.clip((base[None, :, :] - CLOUD_FLOOR) / (1.0 - CLOUD_FLOOR), 0.0, 1.0)
    # Square it: a linear ramp out of the threshold spreads thin haze over
    # everything above the floor, which reads as overcast. Squaring keeps the
    # dense cores and pushes the fringes towards the transparent band, so the
    # clouds have edges and the sky between them is genuinely empty.
    return lumps * lumps * _taper * _edge_fade


# Sampled across BOTH axes of variation — scroll offset and crossfade position —
# because the pattern rolls under a fixed edge window and blends between maps, so
# the peak moves on both. Normalising on one sample let the field drift past 1.0
# later and saturate silently against the top of the ramp.
_SKY_NORM = max(float(_sky_raw(c, m).max())
                for c in range(0, SKY_N, 8) for m in (0.0, 0.5, 1.0)) or 1.0


def sky_field(shift, morph):
    return (_sky_raw(shift, morph) / _SKY_NORM).astype(np.float32)


# Same ordering rule as the sea: real voxels first, then the node, then the TF.
_sky_vol = pycvc.volume(_app)
_sky_vol.set_float_grid(sky_field(0.0, 0.0).ravel().tolist(), SKY_N, SKY_N, SKY_NZ,
                        -SKY_HALF, -SKY_HALF, SKY_BASE, SKY_HALF, SKY_HALF, SKY_TOP)
_sky_grid = _sky_vol.grid()
_sg.addGraphics("forest_sky", _sky_vol)
_sky_node = _sg.volume_node("forest_sky")
# Cloud is not a lit surface. VolumeNode defaults to SetShade(1), and volume
# shading uses the scalar GRADIENT as its normal — on a soft noise field those
# gradients are weak and noisy, so it buys nothing and costs the colour 70% of
# its brightness to the 0.3 ambient term. Absorption/emission only, full
# brightness. (This is NOT why the clouds looked dark; see the opacity note
# below. It is just the right mode for cloud.)
_sky_node.setShading(False)
_sky_node.setAmbient(1.0)
_sky_node.setDiffuse(0.0)
_sky_node.setSpecular(0.0)
# Opacity has to be MINUTE here, and the reason is easy to get wrong: VTK applies
# the opacity function once per ScalarOpacityUnitDistance, which it derives from
# the voxel spacing — about 0.12 world units for this slab. A ray crossing 22
# units of sky therefore compounds the alpha ~180 times, so an innocent-looking
# 0.085 saturates to a solid grey lid. These values are chosen so a ray through
# the densest cloud lands near 0.6 total, and empty sky stays empty.
# Empty sky must be EXACTLY invisible, so the transparent band is wide and flat:
# alpha is pinned to 0 from 0 up to CLOUD_EMPTY, not ramped down towards it. A
# ramp that merely approaches zero still accumulates over the ~180 samples a ray
# takes through the slab, which is what turned clear sky into grey haze and made
# the volume read as a box.
CLOUD_EMPTY = 0.22


def sky_transfer():
    """Transfer function for the cloud slab.

    Two things have to hold at once. Empty sky must be EXACTLY invisible, so the
    transparent band is pinned flat at 0 rather than ramped down towards it — a
    ramp that merely approaches zero still accumulates over the ~180 samples a
    ray takes through the slab, which is what turned clear sky into grey haze
    and made the volume read as a box. And cloud must be OPAQUE where it does
    exist: composited over the sky, a puff that only reaches alpha 0.3 is 30%
    of the background, which reads as dirty smoke. The cores have to reach
    alpha ~1 to composite as white, which is only safe because the zeros are
    pinned and the field is sparse.
    """
    color = [0.0, 0.72, 0.76, 0.82, 0.45, 0.92, 0.94, 0.97, 1.0, 1.00, 1.00, 1.00]
    opacity = [0.0, 0.0, CLOUD_EMPTY, 0.0, 0.60, 0.075, 1.0, 0.200]
    return color, opacity


_sky_node.setTransferFunction(*sky_transfer())

# Bounds over everything, so the grid resizes and the camera's orbit centre lands
# on the island. This is the only thing the script does to the view.
vrhost.set_world_bounds(*_sg.compute_graphics_bounds())

for _k, _v in (("speed", 1), ("waves", 1), ("clouds", 1), ("wind", 1)):
    pycvc.state_set(_app, STATE + "." + _k, str(_v))

_clock = pycvc.world_clock(SIM_DT)
_t = 0.0
_primed = False
_stalled = False
_TREE_AXIS = (0.0, 1.0, 0.0)
_bucket = 0  # which slice of the forest gets re-posed this frame
TREE_STAGGER = 3

print("lsystem_forest: running — %d nodes. Set %s.speed / .waves / .clouds / .wind."
      % (_sg.num_graphics(), STATE))


def _state_float(key, default):
    try:
        return float(pycvc.state_get(_app, STATE + "." + key))
    except (ValueError, TypeError):
        return default


def step(dt):
    global _t, _primed, _stalled, _bucket

    if not _primed:
        # Scene setup (meshes, 2 volumes, node creation) blocks for seconds and
        # the scheduler measures dt from the tick before it — not simulation
        # time. See lsystem_tree.py.
        if dt > PRIME_DT:
            dt = SIM_DT
        else:
            _primed = True

    _clock.set_scale(_state_float("speed", 1.0))
    r = _clock.advance(dt)
    if r.dropped_steps > STALL_QUANTA and not _stalled:
        # A hitch of a frame or two is normal (a heavy first render, a GC
        # pause) and reporting it would just cry wolf; a sustained stall is
        # worth surfacing rather than silently skipping world time.
        _stalled = True
        print("lsystem_forest: stalled — world_clock dropped %d quanta." % r.dropped_steps)
    if not r.steps:
        return
    _t = _clock.t() + r.alpha * SIM_DT  # world seconds, smoothed into the frame

    if _state_float("waves", 1.0):
        _sea_grid[:] = sea_field(_t)  # in-place write into the voxel buffer
        _sea_node.setVolume(_sea_vol)
        _sea_node.setTransferFunction(*sea_transfer(_t))

    if _state_float("clouds", 1.0):
        # Continuous now, so this runs every frame rather than only when an
        # integer column ticked over. That gating is precisely what made the
        # drift snap; paying for it every frame is the cost of fluid motion.
        shift = _t * CLOUD_DRIFT * SKY_N / (2 * SKY_HALF)
        morph = _t / CLOUD_MORPH_S * CLOUD_MAPS
        _sky_grid[:] = sky_field(shift, morph)
        _sky_node.setVolume(_sky_vol)
        # setVolume RESETS the transfer function to VolumeNode's default
        # grayscale ramp, so it must be re-applied after every upload.
        _sky_node.setTransferFunction(*sky_transfer())

    if _state_float("wind", 1.0):
        # Pose only the upper modules; the rest of each tree follows through the
        # graph, so the sway ACCUMULATES and the tips travel further than the
        # trunk — the original demo's motion rather than a rigid lean. Spread
        # over TREE_STAGGER frames: the sway is a ~5 s cycle, so a two-frame lag
        # on part of the forest is invisible.
        _bucket = (_bucket + 1) % TREE_STAGGER
        for i in range(_bucket, len(_swayers), TREE_STAGGER):
            m = _swayers[i]
            a = m.sway * math.sin(1.3 * _t + m.phase)
            m.node.setTransform((m.hang @ mat_rotate(a, *_TREE_AXIS)).ravel().tolist())
