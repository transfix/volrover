# Render performance: the lsystem_forest scene

Measured 2026-08-18 on the sandipaws dev box (Windows, MSVC release build), against
`scripts/examples/lsystem_forest.py` at ~3781 scene-graph nodes.

This document exists because a number quoted in a commit message is not evidence. Every
figure below is reproducible with `scripts/bench/forest_render_breakdown.py`; re-run it
rather than trusting this file after the scene changes.

## Read this before quoting any number below

Two independent runs of the same measurement on the same binary disagree by roughly 2x.
Run 2 was taken while a multi-agent analysis job was loading the machine, which is exactly
the mistake this document warns against. **Re-run the bench on an idle box before treating
any absolute figure as settled.**

| | run 1 (quiet) | run 2 (loaded) |
|---|---:|---:|
| render, shadows off | 217 ms | 403 ms |
| needles | 144 ms (66%) | 249 ms (62%) |
| wood (by subtraction) | 60 ms (28%) | 138 ms (34%) |
| sea + sky volumes | 4 ms (2%) | 140 ms (35%) |
| terrain/grid/axis/sun | 13 ms (6%) | 17 ms (4%) |
| `step()` simulation | ~380 ms | 1634 ms |
| resolution 1/4 pixels | 1.01x | 1.35x |

What is **robust across both runs**, and safe to act on:

* The **needles are the single largest render component** (62-66%).
* **The trees together dominate** the render (~94-96%).
* **`step()` costs as much as, or more than, the entire render.** Whatever the absolute
  numbers, the simulation is not a rounding error next to drawing -- it is half the frame
  or worse.
* There are **3776 actors**, and the implied per-actor cost (57-107 us) is one to two
  orders of magnitude above what a draw call should cost.

What is **not settled**, and needs a clean re-run:

* **The volumes' share.** 2% in one run, 35% in the other. Volume raycasting cost is
  order-sensitive -- how much of the volume is occluded by already-drawn geometry changes
  how early rays terminate -- and the two runs toggled visibility in different orders. Do
  not conclude the volumes are cheap on the strength of run 1.
* **Whether the frame is purely draw-call bound.** 1.01x vs 1.35x for a quarter of the
  pixels. Run 1 says entirely CPU-bound; run 2 says mostly, with some fill-rate component
  (consistent with the volumes mattering more than run 1 suggested).

## Where the render time goes

Pure render cost, no simulation, shadows off, 900x600, measured by toggling node
visibility and re-rendering. Figures from run 1; see the caveat above.

| component | actors | ms | share |
|---|---:|---:|---:|
| needles | 1888 | 144.1 | 66% |
| wood (branch segments) | 1888 | 60.4 | 28% |
| sea + sky volumes | 2 | 4.2 | 2% |
| terrain, grid, axis, sun | ~5 | 12.6 | 6% |
| **total** | **3781** | **217.1** | |

### The needles are the biggest cost, but not because they are lines

The obvious hypothesis is that the needles are slow because they are a lot of `GL_LINES`.
The primitive counts say otherwise:

| | per actor | |
|---|---|---|
| needle cluster | 30 verts | 27 line segments |
| wood segment | 48 verts | 80 triangles |

Each needle cluster has **fewer** primitives than each wood segment, and still costs
**2.4x more per actor** (76 us vs 32 us). Primitive throughput is not the problem.

The problem is that there are **3776 actors**, hence 3776 draw calls per frame, at ~57 us
each. A draw call should cost single-digit microseconds. We are paying fixed per-actor
overhead -- state validation, shader program binding, matrix upload, per-prop culling --
3776 times a frame.

Two consequences worth stating plainly:

* **A shader on the needles alone will not fix this.** Expanding 27 segments into sprites
  or camera-facing quads optimises the part that is already cheap. The actor count, which
  is what actually costs, is unchanged.
* **The per-module hierarchy is the direct cause.** It was introduced so wind accumulates
  down each tree, and it does that well -- the tips move more than the trunk, which is the
  behaviour the original demo had. It also creates one actor per module. That tradeoff was
  never measured at the time.

## What is NOT the problem

Recording these because they were plausible, investigated, and ruled out:

* ~~**Volume rendering.**~~ **Withdrawn.** Run 1 measured the two volumes at 4.2 ms (2%)
  and this section originally declared them exonerated. Run 2 measured 140 ms (35%). The
  claim was made on a single sample and does not hold; volume cost is unresolved.
* **Fill rate / resolution.** Weakly held: 1.01x in run 1, 1.35x in run 2. The frame is
  clearly dominated by CPU-side per-actor work, but a fill-rate component cannot be ruled
  out until the volume question is settled.
* **The Python math in `step()`'s tree posing.** Posing was profiled separately and
  optimised earlier (state_publisher batching, cached paths, echo suppression); node posing
  is now ~2.8 us. The ~380 ms in `step()` is the volume field rebuild, not the trees.

## Two traps when measuring this scene

**A SceneGraph holds ONE renderer.** Constructing a second `SceneRenderer` re-targets the
scene and silently detaches the first. Rendering through the stale handle afterwards costs
0.0 ms and reads as a spectacular optimisation. The first version of the bench measured
shadows after building a second renderer and duly reported that shadow mapping was free.

**`setVisible()` appears to cascade to children.** Needle nodes are leaves, so the needle
figure (144.1 ms, measured directly by hiding only them) is solid. The wood figure is
derived by subtraction -- hiding the wood nodes also hid their needle children, so it
cannot be read directly. The three numbers are self-consistent:

```
needles hidden      -> 73.0 ms   =>  needles      = 217.1 - 73.0 = 144.1
all trees hidden    -> 12.6 ms   =>  wood         =  73.0 - 12.6 =  60.4
                                     trees total  = 144.1 + 60.4 = 204.5 = 217.1 - 12.6
```

## Fixes, ranked by win / effort

### 1. Stop rebuilding the volume fields every frame (no C++ needed)

The largest single win and the cheapest to make. `step()` recomputes both 3-D scalar
fields in numpy on every tick, ~380 ms, which is more than the entire render. The water and
cloud motion does not need to update at frame rate; decoupling it (update every N ticks, or
interpolate between two precomputed fields) roughly halves the frame time on its own.
Fully reversible.

### 2. Turn shadows off by default in this example

Shadow mapping more than doubles the render (217 -> ~590 ms) because
`vtkShadowMapBakerPass` re-renders every actor from every light, and this scene has 3776
actors and 2 lights. Given the actor count, shadows are not affordable here until the
actor count comes down. One line in the example.

### 3. Reduce actor count (the real fix)

3776 -> 70 actors, one per tree, is the change that addresses the actual bottleneck.
Merging each tree's modules into a single polydata is straightforward; the difficulty is
that per-module wind must survive it. The approach that preserves it is skinning: give each
vertex a module index, upload the tree's ~21 module matrices as a uniform array or texture,
and apply the right matrix in a vertex shader. One draw call per tree, wind intact.

This needs new cvcGL surface and is the largest piece of work. Design writeup pending --
see the section below.

## Shader / skinning design

_To be filled in: analysis of the VTK-side options (vtkOpenGLPolyDataMapper shader
replacements, glyph/instanced mappers, geometry shaders) and the minimum cvcGL API needed._

## How the offline films hid all of this

Worth recording so the next person does not repeat it. The film harness
(`make_clip2.py`, scratchpad) differed from the live path in ways that each concealed a
symptom:

* it stubbed `vrhost.set_world_bounds` to a no-op and drove the camera itself, so the
  live camera-reset path was never exercised;
* it predated the lighting work, so it never enabled shadows;
* it rendered 900x600 offscreen;
* and the "22 fps" figure quoted during development was **`step()` cost only** -- the
  render was never timed at all. That is how a 217 ms render went unnoticed for so long.
