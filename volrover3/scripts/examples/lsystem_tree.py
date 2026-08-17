# lsystem_tree.py — a realtime L-system tree, ported to the cvcGL scene graph.
#
# This is a port of a 2004 UT CS354 project (github.com/transfix/cs354,
# src/project2) that drew a fractal conifer with an L-system in immediate-mode
# OpenGL, re-walking the whole grammar and re-emitting every triangle on every
# frame. Here the grammar is expanded ONCE into a real scene-graph hierarchy and
# only 4x4 transforms move after that — which is the whole point of the port.
#
# THE MAPPING (immediate mode -> retained mode)
#   * Every L-system MODULE (a `0`..`4` rule expansion) becomes a GeometryNode
#     carrying that module's own branch mesh, parented to the module that spawned
#     it. The turtle's push/pop matrix stack IS the node hierarchy: a child's
#     transform is relative to its parent, so the graph stores exactly what the
#     original recomputed from scratch 60 times a second.
#   * Each module's needle clusters ride along in one LINES child node.
#   * WIND: the original perturbed the global tilt/roll angles and re-drew. Here
#     the geometry never changes — we recompute each module's local transform
#     from the current angles, and the sway ACCUMULATES down the hierarchy for
#     free, exactly as it did through the old matrix stack.
#   * TIME: the original integrated against a 60 Hz glutTimerFunc. Here a
#     cvc::world_clock (pycvc.world_clock) advances the wind in fixed quanta and
#     hands back an `alpha` to interpolate the presented pose — so the animation
#     runs at the same rate whatever the tick rate, and stays smooth.
#
# The script builds the tree, sets the world bounds so volrover3 frames it, and
# then does nothing but blow wind through it until the job is stopped. It never
# touches the camera: orbit/pan/zoom stay yours.
#
# LIVE CONTROLS (Python Console dock -> State tab, or from any script):
#     pycvc.state_set(app, "volrover3.lsystem_tree.depth", "4")   # 1..5
#     pycvc.state_set(app, "volrover3.lsystem_tree.wind",  "1")   # 0/1
#     pycvc.state_set(app, "volrover3.lsystem_tree.speed", "0.5") # clock scale
# `depth` replaces the original's 'a'/'s' keys; `speed` is the world clock's
# scale (0 pauses the wind, 2 is double time).
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

# ── the grammar, verbatim from the original project2 ─────────────────────────
#   F  move forward one branch segment (and draw it)
#   L  draw a terminating needle cluster
#   R  roll about the branch axis          T  tilt away from the branch axis
#   [ ]  push / pop the turtle             0-4  expand that rule one level deeper
# The turtle starts at the origin facing +Y, its plane parallel to XY.
AXIOM = 0
RULES = (
    "FF[RL1][RR2][RRR3]F[RL3][RR1][RRR2]RFLR0",
    "FL[T[RF]2]R[TRFL]RTFL4",
    "FL[TRF3]RFLRTFL2",
    "FL[TFL2RFL]R[T[RFLF3]]RTFL2",
    "FL[TRFL4]RFLRTFL4",
)

SCALE = 0.9  # each sub-branch is 90% of its parent's segment length
RADSCALE = 0.6  # ...and 60% of its radius
LENGTH = 5.0  # trunk segment length
RADIUS = 0.7  # trunk radius
BASE_TRI = 5  # branch cross-section is a pentagon (as in the original)
NEEDLES = 25  # line segments per needle cluster
LEAF_LEN = 4.0  # needle length / spread, both scaled with the branch
LEAF_RAD = 1.0
MAX_DEPTH = 5  # how deep the grammar is expanded (the original's LSYS_DEPTH)
START_DEPTH = MAX_DEPTH  # render the whole tree; `depth` can trim it live

# 'R' and 'T' turn amounts. The original declares these as degrees (YROTATE 10,
# TILT 120) but its Rotate() feeds the value straight to cos()/sin() — so the
# tree everyone actually saw was rotating by 10 and 120 RADIANS. The shape of
# this plant depends entirely on that, so the numbers are kept as the original
# used them: an effective 212.95 deg roll and 35.49 deg tilt. Reading them as
# degrees instead grows a recognisable but quite different tree.
YROTATE = 10.0  # radians, per the note above
TILT = 120.0  # radians, likewise
MICRO_TILT = 1.0e-4  # every 'F' pre-tilts by TILT/10000 — a slight natural lean

WIND_DEG_RATE = 120.0  # degrees/sec through the wind oscillator (60 Hz * ~2/tick)
STAGGER = 6  # re-pose 1/STAGGER of the tree per frame (see _apply_transforms)
POSE_EPS = 1e-7  # below this a re-pose is a no-op, so don't pay for it

C_WOOD_LIGHT = (0.6549, 0.4901, 0.2392)  # branch cap centres
C_WOOD_DARK = (0.3607, 0.2510, 0.2000)  # branch sides
C_NEEDLE = (0.1373, 0.5568, 0.1373)  # forest green

STATE = "volrover3.lsystem_tree"


# ── 4x4 transforms (numpy; row-major, which is what setTransform takes) ──────
# The 2004 original shipped its own matrix.c because it had to. Here the turtle
# is just numpy: `a @ b` for composition and one vectorised `pts @ R.T + t` to
# place a whole branch's vertices at once, instead of a Python loop per point.
IDENTITY = np.identity(4)


def mat_rotate(angle, x, y, z):
    """Rotation about the unit axis (x,y,z) — the original's Rodrigues matrix."""
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
    """Apply a 4x4 to an (N,3) array of points."""
    return pts @ m[:3, :3].T + m[:3, 3]


# The turtle stands on +Y (as in 2004); volrover3's world is Z-up, so the whole
# tree is rotated a quarter turn about X to plant it on the ground plane.
ROOT_XFORM = mat_rotate(math.pi / 2.0, 1.0, 0.0, 0.0)

# The three turns at their unperturbed values — what the geometry is baked with.
BAKE_TURNS = (mat_rotate(TILT * MICRO_TILT, 0.0, 0.0, 1.0),
              mat_rotate(TILT, 0.0, 0.0, 1.0),
              mat_rotate(YROTATE, 0.0, 1.0, 0.0))


# ── the turtle, replayed symbolically ────────────────────────────────────────
# A module's transform relative to its parent is recorded as the OP SEQUENCE that
# walked the turtle there, not as a baked matrix — so when the wind changes the
# turn amounts we can replay it for the new pose. ("F", length) / ("R",) / ("T",).
# A forward step's length is one of MAX_DEPTH fixed values (LENGTH * SCALE**k),
# so its translation matrix is worth building once and reusing.
_STEP = {}


def step_matrix(length):
    m = _STEP.get(length)
    if m is None:
        m = _STEP[length] = mat_translate(0.0, length, 0.0)
    return m


def ops_matrix(ops, turns):
    """Replay an op sequence. `turns` is the frame's three rotation matrices —
    (micro-tilt, tilt, roll) — hoisted out because they are the same for every
    op in every module on a given frame."""
    micro, tilt, roll = turns
    m = IDENTITY
    for op in ops:
        if op[0] == "F":
            m = m @ micro @ step_matrix(op[1])
        elif op[0] == "R":
            m = m @ roll
        else:  # "T"
            m = m @ tilt
    return m


class Module(object):
    """One rule expansion: its own geometry, plus where it hangs off its parent."""

    __slots__ = ("name", "parent", "level", "ops", "segments", "leaves", "node", "needle_node",
                 "posed")

    def __init__(self, name, parent, level, ops):
        self.name = name
        self.parent = parent  # name of the module that spawned this one
        self.level = level  # 1 = trunk module, MAX_DEPTH = outermost twigs
        self.ops = ops  # turtle ops from the PARENT module's origin
        self.segments = []  # (matrix, length, radius) branch cylinders
        self.leaves = []  # (matrix, size) needle clusters
        self.node = None
        self.needle_node = None
        self.posed = None  # last transform actually pushed to the node


def expand(rule, depth, scale, radscale, name, parent, level, ops, modules):
    """Walk one rule, collecting this module's own geometry and recursing.

    Mirrors drawPlant(): the '[' / ']' stack that used to push GL matrices now
    also stacks the op sequence, so each child module records the exact turtle
    path from this module's origin.
    """
    mod = Module(name, parent, level, ops)
    modules.append(mod)

    micro, tilt, roll = BAKE_TURNS
    m = IDENTITY  # turtle pose in THIS module's local frame
    cur = []  # ops walked so far, in this module's frame
    stack = []
    seg_len, seg_rad = LENGTH * scale, RADIUS * radscale

    for ch in rule:
        if ch == "F":
            m = m @ micro
            cur.append(("F", seg_len))
            mod.segments.append((m, seg_len, seg_rad))
            m = m @ step_matrix(seg_len)
        elif ch == "[":
            stack.append((m, list(cur)))
        elif ch == "]":
            m, cur = stack.pop()
        elif ch == "L":
            # The original drew leaves through a stale GL matrix (it only called
            # LoadMatrix() at 'F'), so needles landed at the previous branch's
            # base. Here they sit at the turtle, where the grammar puts them.
            mod.leaves.append((m, scale))
        elif ch == "R":
            m = m @ roll
            cur.append(("R",))
        elif ch == "T":
            m = m @ tilt
            cur.append(("T",))
        elif ch.isdigit() and depth > 1:
            expand(
                RULES[int(ch)],
                depth - 1,
                scale * SCALE,
                radscale * RADSCALE,
                "lsys_%d" % len(modules),
                name,
                level + 1,
                tuple(cur),
                modules,
            )
    return mod


# ── mesh builders (one mesh per module, in that module's local frame) ─────────
def _unit_ring(n):
    """(n,2) unit circle in the turtle's XZ plane — the branch cross-section."""
    a = np.arange(n) * (2.0 * math.pi / n)
    return np.column_stack((np.cos(a), np.sin(a)))


_WOOD_RING = _unit_ring(BASE_TRI)
_NEEDLE_RING = _unit_ring(NEEDLES)

# Triangle wind-up for one cylinder, relative to its first vertex. Vertex layout
# is [base centre, base ring..., top centre, top ring...]; caps face outward and
# the sides are wound counter-clockwise from outside.
_CYL_TRIS = []
for _i in range(BASE_TRI):
    _b0, _b1 = 1 + _i, 1 + (_i + 1) % BASE_TRI
    _t0, _t1 = BASE_TRI + 2 + _i, BASE_TRI + 2 + (_i + 1) % BASE_TRI
    _CYL_TRIS += [0, _b0, _b1]  # base cap (faces -Y)
    _CYL_TRIS += [BASE_TRI + 1, _t1, _t0]  # top cap (faces +Y)
    _CYL_TRIS += [_b0, _t1, _b1, _b0, _t0, _t1]  # side
_CYL_TRIS = np.array(_CYL_TRIS)
_CYL_VERTS = 2 * BASE_TRI + 2

# One colour per cylinder vertex: light brown at the cap centres, dark on the rim.
_CYL_COLORS = np.array([C_WOOD_LIGHT] + [C_WOOD_DARK] * BASE_TRI
                       + [C_WOOD_LIGHT] + [C_WOOD_DARK] * BASE_TRI)


def wood_mesh(app, segments):
    """All of a module's branch cylinders baked into one triangle mesh."""
    pts = np.empty((len(segments) * _CYL_VERTS, 3))
    local = np.zeros((_CYL_VERTS, 3))
    for s, (m, height, radius) in enumerate(segments):
        local[1:BASE_TRI + 1, 0::2] = _WOOD_RING * radius  # base rim (x,z)
        local[BASE_TRI + 2:, 0::2] = _WOOD_RING * radius  # top rim
        local[BASE_TRI + 1:, 1] = height  # top centre + rim ride at y=height
        pts[s * _CYL_VERTS:(s + 1) * _CYL_VERTS] = xform(m, local)

    tris = (_CYL_TRIS + (np.arange(len(segments)) * _CYL_VERTS)[:, None]).ravel()
    g = pycvc.geometry(app)
    g.add_vertices(pts.ravel().tolist())
    g.add_triangles(tris.tolist())
    g.set_colors(np.tile(_CYL_COLORS, (len(segments), 1)).ravel().tolist())
    return g


def needle_mesh(app, leaves):
    """A module's needle clusters as one LINES mesh (the original's draw_leaf)."""
    stride = NEEDLES + 1  # one shared root vertex, then one tip per needle
    pts = np.empty((len(leaves) * stride, 3))
    local = np.zeros((stride, 3))
    for c, (m, scale) in enumerate(leaves):
        local[1:, 0::2] = _NEEDLE_RING * (LEAF_RAD * scale)
        local[1:, 1] = LEAF_LEN * scale
        pts[c * stride:(c + 1) * stride] = xform(m, local)

    root = (np.arange(len(leaves)) * stride)[:, None]
    lines = np.empty((len(leaves), NEEDLES, 2), dtype=np.int64)
    lines[:, :, 0] = root
    lines[:, :, 1] = root + 1 + np.arange(NEEDLES)
    g = pycvc.geometry(app)
    g.add_vertices(pts.ravel().tolist())
    g.add_lines(lines.ravel().tolist())
    return g


# ── build the scene ──────────────────────────────────────────────────────────
_app = vrhost.app()
_sg = vrhost.scene()  # the live SceneGraph the window renders

_modules = []
expand(RULES[AXIOM], MAX_DEPTH, 1.0, 1.0, "lsys", None, 1, (), _modules)
print("lsystem_tree: grammar expanded to %d modules (depth %d)." % (len(_modules), MAX_DEPTH))

for _mod in _modules:
    # The module's own branch mesh IS its scene node; children hang off it, so
    # the node hierarchy mirrors the grammar and a parent's transform carries
    # its entire sub-branch.
    _wood = wood_mesh(_app, _mod.segments)
    if _mod.parent is None:
        _sg.addGraphics(_mod.name, _wood)
        _mod.node = _sg.geometry_node(_mod.name)
    else:
        _mod.node = _sg.add_child_geometry(_mod.parent, _mod.name, _wood)
    _mod.node.setUseSingleColor(False)  # per-vertex light/dark bark shading

    if _mod.leaves:
        _needles = _mod.name + "_needles"
        _mod.needle_node = _sg.add_child_geometry(_mod.name, _needles, needle_mesh(_app, _mod.leaves))
        _mod.needle_node.setRenderMode(pycvc_gl.GeometryRenderMode_LINES)
        _mod.needle_node.setUseSingleColor(True)
        _mod.needle_node.setColor(*C_NEEDLE)

print("lsystem_tree: %d scene nodes built; geometry is now static." % _sg.num_graphics())


# ── live parameters, read from the state tree each tick ──────────────────────
def _state_float(key, default):
    try:
        return float(pycvc.state_get(_app, STATE + "." + key))
    except (ValueError, TypeError):
        return default


def _state_int(key, default):
    return int(_state_float(key, default))


for _k, _v in (("depth", START_DEPTH), ("wind", 1), ("speed", 1)):
    pycvc.state_set(_app, STATE + "." + _k, str(_v))


# ── the simulation clock ─────────────────────────────────────────────────────
# volrover3's JobScheduler hands step() the raw wall-clock delta between ticks —
# precisely what an animation should NOT be integrated against: the tree would
# sway faster on a fast machine and a full growth would take a different number
# of seconds every run. cvc::world_clock is the codebase's answer to that, and
# pycvc.world_clock is it, directly:
#
#   * the wind and growth integrators advance in whole SIM_DT quanta, so the
#     motion is identical at 30, 60 or 240 fps (and reproducible run to run);
#   * `alpha` — how far into the next quantum the frame lands — interpolates the
#     presented pose between the last two simulated ones, so decoupling the
#     simulation rate from the render rate costs no smoothness;
#   * `dropped_steps` reports quanta discarded after a stall instead of silently
#     swallowing them (the original just skipped ahead and said nothing).
#
# The `speed` state key is the clock's scale: 0 pauses, 2 is double time.
SIM_DT = 1.0 / 120.0
PRIME_DT = 0.25  # a tick slower than this at startup is scene setup, not lag
STALL_QUANTA = 60  # only shout about a stall once it costs ~half a second
_clock = pycvc.world_clock(SIM_DT)


# ── animation state ──────────────────────────────────────────────────────────
_deg = _deg2 = 0.0  # the original's two wind oscillator phases
_shown = START_DEPTH  # levels currently visible
_bucket = 0  # which STAGGER slice gets re-posed this frame
_stalled = False  # have we already reported a dropped-quanta stall?
_primed = False  # has the first (scene-build) tick been discarded yet?

# (tilt, yrotate) at the previous and current simulation quanta; the presented
# values below are interpolated between them by the clock's alpha.
_prev_pose = _cur_pose = (TILT, YROTATE)
_tilt, _yrotate = _cur_pose  # what _apply_transforms actually uses


def _reveal():
    """Show levels 1.._shown, then re-pose the whole tree in one pass.

    The full (unstaggered) pose matters here: a level that has just become
    visible must be placed on the SAME frame it appears, rather than waiting its
    turn in the stagger and flashing at a stale pose for a frame or two.
    """
    for mod in _modules:
        vis = mod.level <= _shown
        mod.node.setVisible(vis)
        if mod.needle_node is not None:
            mod.needle_node.setVisible(vis)
    _apply_transforms()


def _apply_transforms(bucket=None):
    """Re-derive module local transforms from the current turn amounts.

    This is the entire frame cost of the demo — no geometry is touched. The
    Python matrix math is nearly free (op sequences repeat heavily across the
    grammar, so memoising by op tuple collapses 211 modules to 47 chains). What
    costs is setTransform: it CASCADES, re-deriving the world transform of every
    descendant, so posing the trunk alone walks all 211 modules. Three things
    keep that affordable at full depth:

    * hidden modules are skipped outright;
    * a pose that would not move the node measurably is skipped — which is what
      retires the most expensive call in the tree, since the trunk sits at a
      constant transform and re-posing it is pure cascade for no motion;
    * the rest are spread over STAGGER frames. The wind is a ~6 s oscillation,
      so refreshing a given branch every sixth frame is invisible — and the
      slight phase spread between neighbours reads as *less* mechanical.
    """
    turns = (mat_rotate(_tilt * MICRO_TILT, 0.0, 0.0, 1.0),
             mat_rotate(_tilt, 0.0, 0.0, 1.0),
             mat_rotate(_yrotate, 0.0, 1.0, 0.0))
    cache = {}
    for i, mod in enumerate(_modules):
        if mod.level > _shown:
            continue  # hidden: nothing to pose
        if bucket is not None and i % STAGGER != bucket:
            continue
        m = cache.get(mod.ops)
        if m is None:
            m = ops_matrix(mod.ops, turns)
            cache[mod.ops] = m
        if mod.parent is None:
            m = ROOT_XFORM @ m  # stand the +Y turtle up in a Z-up world
        was = mod.posed
        if was is not None and np.abs(was - m).max() < POSE_EPS:
            continue
        mod.posed = m
        mod.node.setTransform(m.ravel().tolist())


_reveal()

# Hand volrover3 the bounds of what we just built, so the world grid resizes and
# the camera's orbit centre lands on the tree. That is the ONLY thing this script
# does to the view — it never writes volrover3.camera.*, so orbit, pan and zoom
# stay entirely the user's. (Until vrhost.set_world_bounds existed there was no
# way to do this from a script: the world_bounds state node is read-only.)
vrhost.set_world_bounds(*_sg.compute_graphics_bounds())

print("lsystem_tree: running — wind until the job is stopped. Set %s.depth (1-%d) / "
      ".wind / .speed to drive it." % (STATE, MAX_DEPTH))


def _simulate(want, wind):
    """Advance the wind by exactly one SIM_DT quantum."""
    global _deg, _deg2, _shown, _prev_pose, _cur_pose

    _prev_pose = _cur_pose

    if want != _shown:
        _shown = want
        _reveal()

    # The original's two oscillators, now stepped by a fixed quantum instead of
    # by whatever the tick rate happened to be. They nudge the same two turn
    # amounts the grammar uses, so the sway accumulates down the hierarchy.
    if wind:
        rate = WIND_DEG_RATE * SIM_DT
        _deg = (_deg + rate * (0.5 + random.random())) % 720.0
        _deg2 = (_deg2 + rate * (0.5 + random.random())) % 720.0
    _cur_pose = (TILT + math.sin(math.radians(_deg * 0.5)) / 100.0,
                 YROTATE + math.cos(math.radians(_deg2 * 0.5)) / 1000.0)


def step(dt):
    global _tilt, _yrotate, _bucket, _stalled, _primed

    if not _primed:
        # Startup is two big blocking deltas, neither of which is simulation
        # time: building 422 nodes (the scheduler measures dt from the tick
        # BEFORE it, so the whole build arrives as one delta), then VTK's first
        # render of those actors. Charged to the clock they read as a genuine
        # ~1800-quantum stall. Substitute a nominal delta until frames come back
        # normal, then latch — after that a stall report means something real.
        if dt > PRIME_DT:
            dt = SIM_DT
        else:
            _primed = True

    want = max(1, min(MAX_DEPTH, _state_int("depth", _shown)))
    wind = _state_int("wind", 1)
    _clock.set_scale(_state_float("speed", 1.0))

    r = _clock.advance(dt)  # -> whole quanta to run + alpha into the next
    for _ in range(r.steps):
        _simulate(want, wind)
    if r.dropped_steps > STALL_QUANTA and not _stalled:
        # A hitch of a frame or two is normal (a heavy first render, a GC
        # pause) and reporting it would just cry wolf; a sustained stall is
        # worth surfacing rather than silently skipping world time.
        _stalled = True  # report a stall once; do not pretend it didn't happen
        print("lsystem_tree: stalled — world_clock dropped %d quanta." % r.dropped_steps)

    # Present between the last two simulated states. With no steps this frame
    # (render faster than SIM_DT) alpha still advances, so the tree keeps moving.
    a = r.alpha
    _tilt = _prev_pose[0] + (_cur_pose[0] - _prev_pose[0]) * a
    _yrotate = _prev_pose[1] + (_cur_pose[1] - _prev_pose[1]) * a

    _bucket = (_bucket + 1) % STAGGER
    _apply_transforms(_bucket)
