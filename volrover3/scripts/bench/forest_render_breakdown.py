"""Where does a lsystem_forest frame actually go?

Not a test. Run it to re-measure on any machine rather than trusting the numbers in
docs/RENDER_PERFORMANCE.md, which go stale the moment the scene changes.

It answers three questions the development-time profiling got wrong:

  1. How much of the frame is SIMULATION vs RENDER? The "22 fps" quoted while the example
     was being built was step() cost only -- the render was never timed, which is how a
     ~200 ms render went unnoticed.
  2. Which nodes cost the most? Measured by toggling visibility and re-rendering, so no
     scene rebuild is involved.
  3. Is the frame draw-call bound or fill-rate bound? Rendering the same scene at a
     quarter of the pixels settles it: ~1x means CPU/draw-call bound, ~4x means fill rate.

Run standalone (no volrover3 app needed -- it stubs vrhost and renders offscreen):

    python volrover3/scripts/bench/forest_render_breakdown.py

CAVEAT: setVisible() appears to cascade to children. Needle nodes are leaves, so the
needle figure is measured directly and is trustworthy; the wood figure is derived by
subtraction, because hiding wood nodes also hides their needle children.
"""

import os
import sys
import time
import types

import pycvc
import pycvc_gl

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.normpath(os.path.join(HERE, "..", "examples", "lsystem_forest.py"))
W, H = 900, 600
REPS = 8

# Frame the island from outside so terrain, trees, sea and cloud deck are all in shot --
# measuring a view that happens to cull most of the scene would prove nothing.
CAM = (250.0, -250.0, 96.0, 0.0, 0.0, 60.0, 0.0, 0.0, 1.0, 50.0, 1.0, 1e5)


def load_scene():
    """Exec the example against a stubbed vrhost, so no Qt app is required."""
    app = pycvc.make_app()
    sg = pycvc_gl.SceneGraph(app, "bench")
    vr = types.ModuleType("vrhost")
    vr.app = lambda: app
    vr.scene = lambda: sg
    vr.set_world_bounds = lambda *a: None  # the bench drives its own camera
    sys.modules["vrhost"] = vr

    mod = types.ModuleType("forest")
    mod.__file__ = SCRIPT
    quiet, real = open(os.devnull, "w"), sys.stdout
    sys.stdout = quiet  # the example is chatty on load
    try:
        exec(compile(open(SCRIPT).read(), SCRIPT, "exec"), mod.__dict__)
        for _ in range(3):
            mod.step(1.0 / 30.0)  # settle: first ticks prime the volumes
    finally:
        sys.stdout = real
    return app, sg, mod, quiet


def main():
    app, sg, mod, quiet = load_scene()

    trees = mod.__dict__.get("_trees", [])
    needles = [m.needles for m in trees if getattr(m, "needles", None) is not None]
    wood = [m.node for m in trees if getattr(m, "node", None) is not None]
    vols = [n for n in (mod.__dict__.get("_sea_node"), mod.__dict__.get("_sky_node")) if n]

    print("scene: %d nodes  (%d needle actors, %d wood actors, %d volumes)"
          % (sg.num_graphics(), len(needles), len(wood), len(vols)))

    r = pycvc_gl.SceneRenderer(sg, W, H, True)
    sg.setShadowsEnabled(False)
    sg.processEvents()
    r.setCamera(*CAM)
    r.render()

    def show(nodes, visible):
        for n in nodes:
            n.setVisible(visible)
        sg.processEvents()

    def render_ms(label, n=REPS):
        r.render()  # warm: first render after a visibility change rebuilds GL state
        t0 = time.time()
        for _ in range(n):
            r.render()
        ms = (time.time() - t0) / n * 1000.0
        print("  %-30s %8.1f ms" % (label, ms))
        return ms

    print("\nRENDER ONLY (no step(), shadows off, %dx%d):" % (W, H))
    base = render_ms("everything")
    show(needles, False)
    no_needles = render_ms("needles hidden")
    show(needles, True)
    show(vols, False)
    no_vols = render_ms("volumes hidden")
    show(vols, True)
    show(needles, False)
    show(wood, False)
    bare = render_ms("all trees hidden")
    show(needles, True)
    show(wood, True)

    needle_ms = base - no_needles
    wood_ms = no_needles - bare
    print("\nATTRIBUTED:")
    for label, ms, count in (("needles", needle_ms, len(needles)),
                             ("wood (by subtraction)", wood_ms, len(wood)),
                             ("volumes", base - no_vols, len(vols)),
                             ("terrain/grid/axis/sun", bare, 0)):
        per = ("  %6.1f us/actor" % (ms * 1000.0 / count)) if count else ""
        print("  %-24s %7.1f ms  (%2.0f%%)%s" % (label, ms, 100.0 * ms / base, per))

    # Simulation vs render: the distinction the original profiling missed.
    sys.stdout = quiet
    t0 = time.time()
    for _ in range(4):
        mod.step(1.0 / 30.0)
        sg.processEvents()
    sim = (time.time() - t0) / 4 * 1000.0
    sys.stdout = sys.__stdout__
    print("\nSIMULATION vs RENDER:")
    print("  step() + processEvents        %8.1f ms" % sim)
    print("  render                        %8.1f ms" % base)
    print("  => a live frame is about      %8.1f ms  (%.1f fps)" % (sim + base, 1000.0 / (sim + base)))

    # Shadows: affordable or not, at this actor count?
    #
    # Measured BEFORE the resolution test on purpose. A SceneGraph holds ONE renderer, so
    # constructing a second SceneRenderer re-targets the scene and silently detaches the
    # first -- rendering through the stale handle then costs 0.0 ms and looks like a
    # miraculous speedup. Anything using `r` must happen before `r2` exists.
    print("")
    sg.setShadowsEnabled(True)
    sg.processEvents()
    with_shadows = render_ms("SHADOWS ON")
    print("  shadow mapping costs %+.1f ms (%.1fx) -- it re-renders every actor per light"
          % (with_shadows - base, with_shadows / base))
    sg.setShadowsEnabled(False)
    sg.processEvents()

    # Draw-call bound or fill-rate bound? This invalidates `r` -- keep it last.
    r2 = pycvc_gl.SceneRenderer(sg, W // 2, H // 2, True)
    r2.setCamera(*CAM)
    r2.render()
    t0 = time.time()
    for _ in range(REPS):
        r2.render()
    small = (time.time() - t0) / REPS * 1000.0
    print("\nRESOLUTION TEST: %dx%d %.1f ms -> %dx%d %.1f ms  (%.2fx for 1/4 the pixels)"
          % (W, H, base, W // 2, H // 2, small, base / small))
    print("  ~1x => CPU/draw-call bound;  ~4x => fill-rate bound")
    print("\nNOTE: run this on an IDLE machine. These timings vary by 2x under concurrent")
    print("      load, and the volume attribution in particular is order-sensitive.")


if __name__ == "__main__":
    main()
