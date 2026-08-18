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
#   ^ v  raise / lower the billow height deposited by the puffs that follow
#   C    a turret: a compact vertical stack, which is what gives cumulus its
#        cauliflower top rather than a smooth dome
CLOUD_AXIOM = "[A][+++++A][-----A][++++++++++A][----------A][+++++++++++++++A]"
CLOUD_RULES = {
    # A is the anvil-ward drift; it throws off B fringes and C turrets, and the
    # ^/v keep the profile from being uniform along the run.
    "A": "FF[+<B]^F[-<C]<F[+<C]vFA",
    "B": "F[+<F]F<[-<F]vB",
    "C": "^<F[+<F][-<F]^<FC",
}
CLOUD_DEPTH = 6
CLOUD_TURN = 32.0     # degrees per + / -
CLOUD_STEP0 = 6.5     # first drift, in map cells
CLOUD_STEP_DECAY = 0.9
CLOUD_PUFF0 = 7.0     # first puff radius, in map cells
CLOUD_PUFF_DECAY = 0.88
CLOUD_MAPS = 2        # independent skies, crossfaded so cloud EVOLVES
CLOUD_MORPH_S = 26.0  # seconds for a full crossfade cycle

# ── the sky volume ───────────────────────────────────────────────────────────
SKY_N = 48
SKY_NZ = 22
SKY_BASE, SKY_TOP = 74.0, 122.0  # a deep slab, so a puff can be round rather than a disc
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
    _sw = 0.020 + 0.016 * random.random()  # ~1-2 deg per level; accumulates down the tree
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
# Water is a lit surface, unlike cloud: volume shading uses the scalar gradient
# as a normal, and the sea's gradient is strongest exactly at the wave surface,
# which is where a highlight belongs. A tight, bright specular turns a
# directional light into a sun glint that travels across the swell.
_sea_node.setShading(True)
_sea_node.setAmbient(0.18)
_sea_node.setDiffuse(0.72)
_sea_node.setSpecular(0.85)
_sea_node.setSpecularPower(70.0)

# ── the sky: a noise slab, scrolled by the wind ──────────────────────────────

# Band-limited noise: a few octaves of sine in x/y, tapered top and bottom so the
# slab has soft faces rather than a hard cut.
# Soft fades at the slab faces so a puff drifting against one is not sliced
# off square. Shape itself comes from the 3-D deposits, not from these.
_w = np.hanning(SKY_N).astype(np.float32) ** 0.5
_edge_fade = np.minimum.outer(_w, _w)[None, :, :]
_zfade = np.sin(np.linspace(0.12, math.pi - 0.12, SKY_NZ)).astype(np.float32)[:, None, None]


def walk_clouds(rng):
    """Run the cloud turtle; return a full 3-D (nz, ny, nx) density field.

    The turtle moves in 3-D and deposits a 3-D Gaussian per F, so a puff is a
    BALL, not a column. This is what rounds the undersides: the previous version
    built a 2-D map and extruded it under a per-column ceiling, which domes the
    top but leaves every cloud sitting on the slab floor with a flat bottom, and
    no amount of shading hides that. Branching in 3-D is also what lets a turret
    genuinely sit above and behind its parent rather than merely being taller.

    Wraps in x (the scroll axis) so the sky can drift forever without a seam.
    """
    field = np.zeros((SKY_NZ, SKY_N, SKY_N), dtype=np.float32)
    gz, gy, gx = np.mgrid[0:SKY_NZ, 0:SKY_N, 0:SKY_N].astype(np.float32)
    # z is squashed relative to x/y: the slab is much thinner than it is wide, so
    # a puff that is round in world units spans far fewer cells vertically.
    zscale = (SKY_N / float(SKY_NZ)) * ((SKY_TOP - SKY_BASE) / (2.0 * SKY_HALF))

    x, y = rng.uniform(0.25, 0.75) * SKY_N, rng.uniform(0.3, 0.7) * SKY_N
    z = SKY_NZ * 0.42
    head = rng.uniform(0.0, 360.0)
    step, puff, depth = CLOUD_STEP0, CLOUD_PUFF0, 0
    climb = 0.0
    stack = []
    todo = list(CLOUD_AXIOM)
    guard = 0
    while todo and guard < 4000:
        guard += 1
        c = todo.pop(0)
        if c == "F":
            x = (x + step * math.cos(math.radians(head))) % SKY_N
            y = min(max(y + step * math.sin(math.radians(head)), 0.0), float(SKY_N - 1))
            z = min(max(z + climb, 1.0), float(SKY_NZ - 2))
            dx = np.abs(gx - x)
            dx = np.minimum(dx, SKY_N - dx)  # nearest image across the seam
            dz = (gz - z) * zscale
            d2 = dx * dx + (gy - y) ** 2 + dz * dz
            field += np.exp(-d2 / (2.0 * puff * puff)).astype(np.float32)
        elif c == "+":
            head += CLOUD_TURN
        elif c == "-":
            head -= CLOUD_TURN
        elif c == "<":
            puff *= CLOUD_PUFF_DECAY
        elif c == "^":
            climb += 0.55          # a turret climbs as it goes
        elif c == "v":
            climb -= 0.45          # a fringe sags away underneath
        elif c == "[":
            stack.append((x, y, z, head, step, puff, depth, climb))
        elif c == "]":
            if stack:
                x, y, z, head, step, puff, depth, climb = stack.pop()
        elif c in CLOUD_RULES and depth < CLOUD_DEPTH:
            depth += 1
            step *= CLOUD_STEP_DECAY
            todo = list(CLOUD_RULES[c]) + todo

    # Normalise on a high percentile, not the max: one spot where several
    # branches overlap would otherwise set the scale and push the rest of the
    # sky under the threshold.
    m = float(np.percentile(field, 99.9))
    if m > 0:
        field = np.clip(field / m, 0.0, 1.0).astype(np.float32)
    # Fade at the slab faces so nothing is cut off square.
    return (field * _edge_fade * _zfade).astype(np.float32)


_rng = np.random.default_rng(SEED)
# Several independent skies. One would only ever translate; crossfading between
# them is what makes the cloud EVOLVE — puffs grow and dissolve in place rather
# than sliding past like a painted backdrop.
_cloud_maps = [walk_clouds(_rng) for _ in range(CLOUD_MAPS)]
# Soft fade at the top and bottom faces of the slab, so a puff that drifts
# against them is not sliced off square. Shape is now carried by the 3-D
# deposits themselves, not by this.
# The z faces are tapered above, but the x/y faces were a hard cut: cloud density
# ran right up to the slab boundary, so the volume ended in a straight vertical
# wall hanging in the sky. Fade to zero at every face instead.

# A high threshold is what makes this read as CLOUD rather than as overcast: only
# the top third of the noise becomes anything at all, so the slab is mostly holes
# and you can see the sky (and the island) through the gaps.
CLOUD_FLOOR = 0.10


def _sky_raw(shift, morph):
    """Density at a CONTINUOUS scroll offset and crossfade position.

    Sub-cell scroll (adjacent integer shifts blended by the fractional part, so
    the drift is smooth at any speed rather than jumping a whole column), and a
    smoothstepped crossfade between independently grown 3-D skies so the cloud
    changes SHAPE as it travels instead of sliding past rigidly.
    """
    i = int(math.floor(morph)) % CLOUD_MAPS
    j = (i + 1) % CLOUD_MAPS
    u = morph - math.floor(morph)
    u = u * u * (3.0 - 2.0 * u)
    vol = (1.0 - u) * _cloud_maps[i] + u * _cloud_maps[j]

    k = int(math.floor(shift))
    f = shift - k
    vol = (1.0 - f) * np.roll(vol, k, axis=2) + f * np.roll(vol, k + 1, axis=2)

    lumps = np.clip((vol - CLOUD_FLOOR) / (1.0 - CLOUD_FLOOR), 0.0, 1.0)
    # Squaring keeps dense cores and pushes fringes into the transparent band,
    # so clouds have edges and the sky between them is genuinely empty.
    return lumps * lumps


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
# Cloud shading was off while the field was a flat extruded slab — there were no
# gradients worth lighting, only noise. Now that the grammar gives each column a
# billow height the field HAS shape, so a directional light picks out the tops
# and leaves the undersides dim, which is most of what makes cloud read as
# volume rather than as fog. Ambient stays high so the shadowed side is grey-blue
# rather than black.
_sky_node.setShading(True)
_sky_node.setAmbient(0.55)
_sky_node.setDiffuse(0.85)
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

# -- the afternoon sun -------------------------------------------------------
# The SCENE says what time of day it is, instead of inheriting whatever light
# the host window happened to configure. That matters here because VTK's default
# is a HEADLIGHT: it rides the camera, so it lights every surface head-on and
# flattens exactly what this scene is made of -- the billow of a cloud, the roll
# of a swell, the depth of a canopy. Fly around with a headlight and nothing
# changes; fly around with a fixed sun and the island turns in the light.
#
# Elevation 34 degrees reads as mid-afternoon: high enough to catch the cloud
# tops, low enough that the water throws a specular track back toward the camera
# and the billows keep a shaded underside. The sun is warmed a little, and a dim
# COOL fill from the opposite side keeps shadowed faces blue rather than black --
# that fill stands in for sky light, it is not pretending to be a second sun.
SUN_AZ = -52.0
_sun = _sg.addDirectionalLight(SUN_AZ, 34.0, 1.0, 0.94, 0.82, 0.95)
_fill = _sg.addDirectionalLight(128.0, 52.0, 0.55, 0.66, 0.85, 0.55)

# Shadows are OFF by default here, and the reason is measured (docs/RENDER_PERFORMANCE.md
# fix #2): vtkShadowMapBakerPass re-renders every actor from every light, so at this
# scene's ~3776 actors x 2 lights it more than doubles the render (217 -> ~590 ms) and
# turns a live orbit into a slideshow. Shadows are not affordable until the actor count
# comes down (the skinning work, fix #3), so they are opt-in rather than on by default.
#
# Flip them on live to see them (accepting the frame cost):
#     pycvc.state_set(app, "volrover3.lsystem_forest.shadows", "1")
#
# setShadowsEnabled installs passes on the RENDER TARGET, so it can only succeed once the
# scene is attached to one. Under volrover3 that has happened and it returns True; in a
# headless harness that builds the scene before any renderer exists it returns False, and
# the scene is simply lit without shadows rather than failing to load.
_shadows_on = False
_sg.setShadowsEnabled(False)
print("lsystem_forest: sun at az %.0f el 34 (%d lights); shadows off by default "
      "(set %s.shadows=1 to enable -- costs ~2x render at this actor count)"
      % (SUN_AZ, _sg.numLights(), STATE))

# Bounds over everything, so the grid resizes and the camera's orbit centre lands
# on the island. This is the only thing the script does to the view.
vrhost.set_world_bounds(*_sg.compute_graphics_bounds())

# -- the sun itself, and an honest half of a lens flare -----------------------
# Added AFTER set_world_bounds deliberately. The disc sits ~430 units out, and
# folding that into compute_graphics_bounds would push the camera reset back far
# enough to turn the island into a speck. The sun is scenery, not world.
#
# Flat-lit rather than shaded: ambient 1, no diffuse, no specular. A shaded ball
# out there would be lit from BEHIND -- by its own light -- and render as a dark
# disc, which is the one thing a sun must not be.
SUN_DIST = 430.0
SUN_R = 13.0


def sun_dir(az_deg, el_deg):
    """Unit vector toward the sun, matching addDirectionalLight's az/el convention."""
    az, el = math.radians(az_deg), math.radians(el_deg)
    return np.array([math.cos(el) * math.sin(az), -math.cos(el) * math.cos(az), math.sin(el)])


def disc_mesh(app, centre, normal, radius, seg=48):
    """A triangle-fan disc facing `normal`, built in the plane perpendicular to it."""
    n = np.asarray(normal, dtype=float)
    n = n / np.linalg.norm(n)
    # Any vector not parallel to n works to start the basis; swap near the poles.
    up = np.array([0.0, 0.0, 1.0]) if abs(n[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    u = np.cross(n, up)
    u /= np.linalg.norm(u)
    v = np.cross(n, u)
    th = np.linspace(0.0, 2.0 * math.pi, seg, endpoint=False)
    c = np.asarray(centre, dtype=float)
    rim = c[None, :] + radius * (np.cos(th)[:, None] * u[None, :] + np.sin(th)[:, None] * v[None, :])
    pts = np.vstack((c[None, :], rim))
    i = np.arange(seg)
    tris = np.column_stack((np.zeros(seg, dtype=int), 1 + i, 1 + (i + 1) % seg)).ravel()
    g = pycvc.geometry(app)
    g.add_vertices(pts.ravel().tolist())
    g.add_triangles(tris.tolist())
    return g


def _sun_geoms(el_deg):
    d = sun_dir(SUN_AZ, el_deg)
    c = d * SUN_DIST
    face = -d  # face the origin, which is where the camera orbits
    # The halo sits fractionally FARTHER out so the disc always wins the depth test.
    return (disc_mesh(_app, c, face, SUN_R),
            disc_mesh(_app, c * 1.02, face, SUN_R * 3.2))


_disc_g, _halo_g = _sun_geoms(34.0)
_sun_disc = _sg.addGraphics("forest_sun", _disc_g)
_sun_disc.setColor(1.0, 0.97, 0.88)
_sun_disc.setAmbient(1.0)
_sun_disc.setDiffuse(0.0)
_sun_disc.setSpecular(0.0)

# A wider, fainter disc behind it: the corona/bloom you actually see around a
# bright source. This is the half of "lens flare" that can be done honestly in
# world space. The ghosts -- the chain of coloured discs along the axis from the
# sun through screen centre -- are a SCREEN-space effect: their positions depend
# on where the sun projects in the image, so they need a post-processing pass
# over the rendered frame, which cvcGL does not expose yet. Faking them with
# world-space quads would put them at fixed 3-D points that only line up from
# one camera angle, and this scene is meant to be flown around.
_sun_halo = _sg.addGraphics("forest_sun_halo", _halo_g)
_sun_halo.setColor(1.0, 0.90, 0.72)
_sun_halo.setAmbient(1.0)
_sun_halo.setDiffuse(0.0)
_sun_halo.setSpecular(0.0)
_sun_halo.setOpacity(0.22)


def place_sun(el_deg):
    """Aim the light and move the disc together, so they cannot disagree."""
    _sg.setLightDirection(_sun, SUN_AZ, el_deg)
    dg, hg = _sun_geoms(el_deg)
    _sun_disc.setGeometry(dg)
    _sun_halo.setGeometry(hg)

for _k, _v in (("speed", 1), ("waves", 1), ("clouds", 1), ("wind", 1),
               ("sun", 34), ("shadows", 0)):
    pycvc.state_set(_app, STATE + "." + _k, str(_v))

_clock = pycvc.world_clock(SIM_DT)
_t = 0.0
_primed = False
_stalled = False
_TREE_AXIS = (0.0, 1.0, 0.0)
_sun_el = 34.0  # last elevation pushed, so we only re-aim when it moves
_bucket = 0  # which slice of the forest gets re-posed this frame
TREE_STAGGER = 3

# ── decouple the volume uploads from the frame rate (fix #1) ─────────────────
# The field NUMPY is cheap (~1-2 ms for both volumes). The cost of a per-frame
# volume update is setVolume(): it re-imports the whole volume into VTK every call
# -- deep-copies the cvc::volume, re-runs the vtkImageData setup + a memcpy of every
# voxel, and RESETS the transfer function (which is why the TF has to be re-applied
# below) -- and a live renderer then re-uploads the field to the GPU.
# None of that has to happen at frame rate: water and cloud motion reads fine at
# ~30 Hz. So each volume is refreshed once every VOL_STRIDE sim-steps, and the two
# are offset so at most ONE volume re-uploads on any given frame. That cuts the
# per-frame volume upload work ~4x (two uploads/frame -> one every other frame) with
# no visible change to the motion. Set VOL_STRIDE=1 to restore per-frame updates.
# The real fix is a lightweight in-place VolumeNode refresh in cvcGL (see the doc);
# until that lands, striding the uploads is the reversible, C++-free mitigation.
VOL_STRIDE = 2
_vol_phase = 0

print("lsystem_forest: running — %d nodes. Set %s.speed / .waves / .clouds / .wind / .sun / .shadows."
      % (_sg.num_graphics(), STATE))


def _state_float(key, default):
    try:
        return float(pycvc.state_get(_app, STATE + "." + key))
    except (ValueError, TypeError):
        return default


def step(dt):
    global _t, _primed, _stalled, _bucket, _sun_el, _shadows_on, _vol_phase

    if not _primed:
        # Scene setup (meshes, 2 volumes, node creation) blocks for seconds and
        # the scheduler measures dt from the tick before it — not simulation
        # time. See lsystem_tree.py.
        if dt > PRIME_DT:
            dt = SIM_DT
        else:
            _primed = True

    el = _state_float("sun", 34.0)
    if abs(el - _sun_el) > 1e-3:
        # Sweep this down toward the horizon for evening light: the specular
        # track on the water stretches and the cloud undersides catch the warm.
        _sun_el = el
        place_sun(el)

    want_shadows = _state_float("shadows", 0.0) >= 0.5
    if want_shadows != _shadows_on:
        # Off by default (fix #2): shadow mapping re-renders every actor per light,
        # ~2x the frame at this actor count. Toggle it only when it changes, like
        # the sun elevation above -- installing/removing the passes is not free.
        _sg.setShadowsEnabled(want_shadows)
        _shadows_on = want_shadows

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

    # Refresh at most one volume per frame (fix #1): sea on one phase, sky on the
    # next, each every VOL_STRIDE sim-steps. The field is still sampled at the true
    # _t whenever it IS refreshed, so the motion keeps its phase -- it just updates
    # at ~30 Hz instead of frame rate, which water and cloud do not need to beat.
    _vol_phase += 1
    _sea_due = (_vol_phase % VOL_STRIDE) == 0
    _sky_due = (_vol_phase % VOL_STRIDE) == (1 % VOL_STRIDE)

    if _sea_due and _state_float("waves", 1.0):
        _sea_grid[:] = sea_field(_t)  # in-place write into the voxel buffer
        _sea_node.setVolume(_sea_vol)
        _sea_node.setTransferFunction(*sea_transfer(_t))

    if _sky_due and _state_float("clouds", 1.0):
        # Continuous drift (sub-cell scroll + crossfade), so the field is smooth at
        # any refresh rate -- unlike the old integer-column gate, which is what made
        # the drift snap. Striding the UPLOAD (not the phase) keeps it fluid.
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
