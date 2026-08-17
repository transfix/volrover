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
# WHY THE TREES ARE FLAT HERE
#   lsystem_tree.py gives every L-system module its own GeometryNode so wind can
#   accumulate down the hierarchy. That costs ~422 nodes for ONE tree, and node
#   count is what the frame budget actually goes on (setTransform cascades to
#   every descendant). A forest cannot afford it, so each tree here is BAKED into
#   one wood mesh + one needle mesh and sways as a whole from its root node. Same
#   grammar, opposite end of the granularity trade — that contrast is the point.
#
# THE TWO VOLUMES
#   * SEA — a cvc::volume over the island's footprint. Its scalar field is depth
#     below a travelling wave surface, and it is zeroed anywhere the terrain
#     stands above sea level, so water appears exactly in the hollows the dirt
#     patches carved. Both the field AND the transfer function animate: the field
#     carries the swell, the transfer function breathes the surface highlight.
#   * SKY — a second volume slab overhead holding a band-limited noise field,
#     scrolled by the wind. A soft ramp transfer function turns it into cloud.
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


def bake_tree(rule, depth, scale, radscale, m, segs, leaves):
    """Recursively flatten the tree grammar into segment/leaf placements."""
    stack = []
    seg_len, seg_rad = T_LENGTH * scale, T_RADIUS * radscale
    step = mat_translate(0.0, seg_len, 0.0)
    for ch in rule:
        if ch == "F":
            m = m @ _TURN_MICRO
            segs.append((m, seg_len, seg_rad))
            m = m @ step
        elif ch == "[":
            stack.append(m)
        elif ch == "]":
            m = stack.pop()
        elif ch == "L":
            leaves.append((m, scale))
        elif ch == "R":
            m = m @ _TURN_ROLL
        elif ch == "T":
            m = m @ _TURN_TILT
        elif ch.isdigit() and depth > 1:
            bake_tree(TREE_RULES[int(ch)], depth - 1, scale * T_SCALE,
                      radscale * T_RADSCALE, m, segs, leaves)
    return m


def tree_meshes(app, maturity, size):
    segs, leaves = [], []
    bake_tree(TREE_RULES[0], maturity, size, size, TREE_UP, segs, leaves)

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

_trees = []
for _n, (_sx, _sy) in enumerate(_dry):
    _hz = height_at(_H, _sx, _sy)
    _mat = random.choice(MATURITY)
    _size = random.uniform(*TREE_SIZE)
    _wood, _needles = tree_meshes(_app, _mat, _size)
    _wn, _nn = "forest_tree%d" % _n, "forest_tree%d_needles" % _n
    _sg.addGraphics(_wn, _wood)
    _node = _sg.geometry_node(_wn)
    _node.setUseSingleColor(False)
    _leaf = _sg.add_child_geometry(_wn, _nn, _needles)
    _leaf.setRenderMode(pycvc_gl.GeometryRenderMode_LINES)
    _leaf.setUseSingleColor(True)
    _leaf.setColor(*C_NEEDLE)
    _trees.append({"node": _node, "pos": (_sx, _sy, _hz),
                   "phase": random.uniform(0.0, 2 * math.pi),
                   "sway": 0.012 + 0.010 * random.random()})
print("lsystem_forest: %d trees planted, maturities %s." %
      (len(_trees), sorted(set(MATURITY))))


# ── the sea: a volume whose field is depth under a travelling wave ───────────
_sea_ax = np.linspace(-HALF, HALF, SEA_N)
_sea_z = np.linspace(SEA_FLOOR, SEA_TOP, SEA_NZ)
_sgx, _sgy = np.meshgrid(_sea_ax, _sea_ax, indexing="xy")
# Terrain height sampled on the sea grid, so the water can be masked off wherever
# the ground stands above it — water shows up only in the hollows.
_idx = np.clip(((_sea_ax + HALF) / (2 * HALF) * (TERRAIN_N - 1)).round().astype(int),
               0, TERRAIN_N - 1)
_sea_terrain = _H[np.ix_(_idx, _idx)]
_wet = (_sea_terrain < SEA_LEVEL).astype(np.float32)  # 1 where there is sea floor

_sea_vol = pycvc.volume(_app)
# Phase of the swell at each column, precomputed — only the time term changes.
_wave_phase = (2 * math.pi / WAVE_LEN) * (_sgx + 0.6 * _sgy)
_zcol = _sea_z[:, None, None].astype(np.float32)


def sea_field(t):
    """Depth below the wave surface, 0 on land — the sea's scalar field."""
    surf = SEA_LEVEL + WAVE_AMP * (np.sin(_wave_phase - WAVE_SPEED * t * 0.1)
                                   + 0.45 * np.sin(1.7 * _wave_phase + WAVE_SPEED * t * 0.13))
    depth = np.clip((surf[None, :, :] - _zcol) / 6.0, 0.0, 1.0)
    return (depth * _wet[None, :, :]).astype(np.float32)


def sea_transfer(t):
    """Transfer function for the sea. The knee where opacity climbs is what reads
    as the surface, so drifting it slightly makes the water look alive even
    between field updates."""
    k = 0.16 + 0.05 * math.sin(t * 0.9)
    color = [0.00, 0.02, 0.10, 0.22,
             0.35, 0.05, 0.28, 0.48,
             0.70, 0.10, 0.45, 0.62,
             1.00, 0.55, 0.80, 0.85]  # foam-ish at the very top
    opacity = [0.00, 0.00, max(k - 0.10, 0.001), 0.00, k, 0.32, 0.75, 0.72, 1.00, 0.90]
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
_kx = np.linspace(-3.0, 3.0, SKY_N)
_cx, _cy = np.meshgrid(_kx, _kx, indexing="xy")
_rng = np.random.default_rng(SEED)
_cloud2d = np.zeros((SKY_N, SKY_N), dtype=np.float32)
for _oct in (1.0, 2.3, 4.7, 9.1):  # enough octaves that the edges are ragged
    _px, _py = _rng.uniform(0, 2 * math.pi, 2)
    _cloud2d += (1.0 / _oct) * np.sin(_oct * _cx + _px) * np.cos(_oct * _cy * 1.3 + _py)
_cloud2d = (_cloud2d - _cloud2d.min()) / (np.ptp(_cloud2d) or 1.0)
_taper = np.sin(np.linspace(0, math.pi, SKY_NZ)).astype(np.float32)[:, None, None]

# A high threshold is what makes this read as CLOUD rather than as overcast: only
# the top third of the noise becomes anything at all, so the slab is mostly holes
# and you can see the sky (and the island) through the gaps.
CLOUD_FLOOR = 0.62


def sky_field(shift_cols):
    base = np.roll(_cloud2d, int(shift_cols) % SKY_N, axis=1)
    lumps = np.clip((base[None, :, :] - CLOUD_FLOOR) / (1.0 - CLOUD_FLOOR), 0.0, 1.0)
    return (lumps * _taper).astype(np.float32)


# Same ordering rule as the sea: real voxels first, then the node, then the TF.
_sky_vol = pycvc.volume(_app)
_sky_vol.set_float_grid(sky_field(0).ravel().tolist(), SKY_N, SKY_N, SKY_NZ,
                        -SKY_HALF, -SKY_HALF, SKY_BASE, SKY_HALF, SKY_HALF, SKY_TOP)
_sky_grid = _sky_vol.grid()
_sg.addGraphics("forest_sky", _sky_vol)
_sky_node = _sg.volume_node("forest_sky")
# Opacity has to be MINUTE here, and the reason is easy to get wrong: VTK applies
# the opacity function once per ScalarOpacityUnitDistance, which it derives from
# the voxel spacing — about 0.12 world units for this slab. A ray crossing 22
# units of sky therefore compounds the alpha ~180 times, so an innocent-looking
# 0.085 saturates to a solid grey lid. These values are chosen so a ray through
# the densest cloud lands near 0.6 total, and empty sky stays empty.
_sky_node.setTransferFunction(
    [0.0, 0.55, 0.60, 0.70, 0.5, 0.85, 0.88, 0.93, 1.0, 1.00, 1.00, 1.00],
    [0.0, 0.0, 0.30, 0.0, 0.65, 0.0035, 1.0, 0.0110])

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
_sky_col = -1  # last cloud column offset uploaded (see step)
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
    global _t, _primed, _stalled, _sky_col, _bucket

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
        # The slab scrolls under a whole column per second, so re-uploading it
        # every frame would push identical voxels 60 times a second. Only when
        # the integer offset actually changes is there anything new to send.
        col = int(_t * CLOUD_DRIFT * SKY_N / (2 * SKY_HALF))
        if col != _sky_col:
            _sky_col = col
            _sky_grid[:] = sky_field(col)
            _sky_node.setVolume(_sky_vol)

    if _state_float("wind", 1.0):
        # Each tree leans on its own phase. One transform per tree — the forest's
        # counterpart to lsystem_tree.py's per-module sway. Spread over
        # TREE_STAGGER frames: the sway is a ~5 s cycle, so a two-frame lag on
        # part of the forest is invisible and it keeps setTransform (which
        # cascades into each tree's needle child) off the critical path.
        _bucket = (_bucket + 1) % TREE_STAGGER
        for i in range(_bucket, len(_trees), TREE_STAGGER):
            tr = _trees[i]
            a = tr["sway"] * math.sin(1.3 * _t + tr["phase"])
            m = mat_translate(*tr["pos"]) @ mat_rotate(a, *_TREE_AXIS)
            tr["node"].setTransform(m.ravel().tolist())
