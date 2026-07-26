# chase_cam_demo.py — a self-contained volrover3 example: build a live 3D scene
# (rolling terrain + a marker) and animate it, filmed by a third-person chase
# camera, all from the embedded Python interpreter. No external packages — only
# pycvc / pycvc_gl (the SDK bindings) and vrhost (the host shim).
#
# HOW TO RUN (inside a running volrover3): Python Console dock -> "Jobs" tab ->
# "Load Script..." -> pick this file. It appears in the jobs table; the marker
# walks a loop over the terrain and the camera follows. Select it -> Stop to end.
#
# THE CONTRACT: define a module-level step(dt); import-time code runs once at
# submit, step(dt) runs every tick. This is the generic "live render window
# changing over time" pattern — swap the marker's path for your own simulation.

import math

import pycvc
import pycvc_gl
import vrhost

# ── adopt the host's live app + scene (no singleton) ─────────────────────────
_app = vrhost.app()
_sg = vrhost.scene()  # the live SceneGraph the window renders

_CAM = "volrover3.camera"


# ── tiny scene builders (single-color materials render faithfully) ───────────
def _add_mesh(name, verts, tris, color):
    g = pycvc.geometry(_app)
    g.add_vertices(list(verts))
    g.add_triangles(list(tris))
    _sg.addGraphics(name, g)
    node = _sg.geometry_node(name)
    node.setUseSingleColor(True)
    node.setColor(*color)


def _terrain(n=48, half=100.0, amp=12.0):
    """An analytic rolling-hills heightfield as a triangle mesh over [-half,half]^2."""

    def h(x, y):
        return (
            amp * math.exp(-((x * x + y * y) / 2600.0))
            + 0.5 * amp * math.exp(-(((x - 45) ** 2 + (y + 30) ** 2) / 900.0))
            + 2.0 * math.sin(x * 0.06) * math.cos(y * 0.05)
        )

    step = 2 * half / (n - 1)
    verts = []
    for i in range(n):
        y = -half + i * step
        for j in range(n):
            x = -half + j * step
            verts += [x, y, h(x, y)]
    tris = []
    for i in range(n - 1):
        for j in range(n - 1):
            v = i * n + j
            tris += [v, v + 1, v + n, v + 1, v + n + 1, v + n]
    return verts, tris, h


_verts, _tris, _height = _terrain()
_add_mesh("terrain", _verts, _tris, (0.34, 0.42, 0.28))

_S = 5.0  # marker size; base at local z=0, apex at +S (sits on the ground)
_add_mesh(
    "agent",
    [0, 0, _S, -_S, -_S, 0, _S, -_S, 0, 0, _S, 0],
    [0, 1, 2, 0, 2, 3, 0, 3, 1, 1, 3, 2],
    (0.95, 0.15, 0.15),
)


# ── the agent's demo path: a wandering loop, draped onto the terrain ─────────
def _agent_pos(t, radius=70.0, wander=8.0, loop_s=22.0, lift=0.4):
    th = 2 * math.pi * ((t / loop_s) % 1.0)
    r = radius + wander * math.sin(3 * th)
    x, y = r * math.cos(th), r * math.sin(th)
    return (x, y, _height(x, y) + lift)


# ── a position-driven chase camera (heading from smoothed velocity + damping) ─
class _Chase:
    def __init__(self, back=55.0, height=40.0, up_look=3.0, vel_tau=0.4, cam_tau=0.55):
        self.back, self.height, self.up_look = back, height, up_look
        self.vel_tau, self.cam_tau = vel_tau, cam_tau
        self.p = self.pp = self.v = self.head = self.eye = self.tgt = None

    @staticmethod
    def _ema(prev, new, dt, tau):
        if prev is None:
            return list(new)
        a = 1.0 - math.exp(-dt / tau)
        return [prev[i] + (new[i] - prev[i]) * a for i in range(len(new))]

    def update(self, pos, dt):
        dt = max(dt, 1e-4)
        self.p = self._ema(self.p, pos, dt, 0.15)
        if self.pp is not None:
            self.v = self._ema(self.v, [(self.p[i] - self.pp[i]) / dt for i in range(3)], dt, self.vel_tau)
        self.pp = list(self.p)
        if self.v is not None:
            s = math.hypot(self.v[0], self.v[1])
            if s >= 0.05:
                self.head = (self.v[0] / s, self.v[1] / s)
        hx, hy = self.head if self.head else (1.0, 0.0)
        p = self.p
        te = [p[0] - hx * self.back, p[1] - hy * self.back, p[2] + self.height]
        tt = [p[0], p[1], p[2] + self.up_look]
        self.eye = self._ema(self.eye, te, dt, self.cam_tau)
        self.tgt = self._ema(self.tgt, tt, dt, self.cam_tau)
        return self.eye, self.tgt


def _drive_cam(eye, tgt):
    vx, vy, vz = (tgt[i] - eye[i] for i in range(3))
    m = math.sqrt(vx * vx + vy * vy + vz * vz) or 1.0

    def s(k, val):
        pycvc.state_set(_app, _CAM + "." + k, "%.6f" % float(val))

    s("position.x", eye[0]); s("position.y", eye[1]); s("position.z", eye[2])
    s("view_direction.x", vx / m); s("view_direction.y", vy / m); s("view_direction.z", vz / m)
    s("up_vector.x", 0.0); s("up_vector.y", 0.0); s("up_vector.z", 1.0)
    s("fov", 55.0)


_chase = _Chase()
_t = 0.0
for _i in range(30):  # prime the camera to a settled pose
    _e, _g = _chase.update(_agent_pos(_i / 30.0), 1 / 30.0)
_drive_cam(_e, _g)
_t = 1.0
_sg.getGraphics("agent").setPosition(*_agent_pos(_t))
print("chase_cam_demo: marker walks a loop over the terrain; camera follows. Stop from the Jobs tab.")


def step(dt):
    global _t
    _t += dt
    pos = _agent_pos(_t)
    _sg.getGraphics("agent").setPosition(*pos)  # move in place — no rebuild
    _drive_cam(*_chase.update(pos, dt))
