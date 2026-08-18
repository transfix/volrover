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
  is now ~2.8 us.
* ~~**The numpy volume field rebuild.**~~ **Corrected.** This section originally said
  "the ~380 ms in `step()` is the volume field rebuild". Timed directly (STEP BREAKDOWN in
  the bench), the numpy that builds both fields is ~1-2 ms — it is not where `step()` goes.
  The cost is the per-frame `setVolume()` UPLOAD, not the field math; see fix #1.

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

### 1. Decouple the volume UPLOADS from the frame rate (done — no C++)

The largest cheap win, but not for the reason first recorded here. The claim was that
`step()` "recomputes both 3-D scalar fields in numpy, ~380 ms". Timed directly (STEP
BREAKDOWN in the bench), that numpy is ~1-2 ms for both volumes together. The field
rebuild is not the cost.

What costs is `_sea_node.setVolume()` / `_sky_node.setVolume()`, called every frame.
cvcGL's `VolumeNode::setVolume` (`src/cvcGL/VolumeNode.cpp`) does a full RE-IMPORT on every
call: it deep-copies the whole `cvc::volume` (`make_shared<cvc::volume>(vol)`), re-runs
`updateImageData` (SetDimensions/Spacing/Origin + `AllocateScalars` + a `memcpy` of every
voxel), and calls `setDefaultTransferFunction()` which RESETS the transfer function — which
is exactly why the example has to re-apply the TF after every upload — plus state-tree writes
and per-call logging. A live renderer then re-uploads the field to the GPU. All of that runs
even though the example already wrote the new voxels in place through the zero-copy `grid()`
view; the re-import throws that write away.

**Implemented mitigation** (reversible, C++-free): `VOL_STRIDE` in the example refreshes each
volume once every 2 sim-steps and offsets sea against sky, so at most one volume re-uploads
on any frame — two uploads/frame down to one every other frame, ~4x fewer, with no visible
change (water and cloud read fine at ~30 Hz). `VOL_STRIDE=1` restores per-frame updates.

**The real fix is in cvcGL** and retires the re-import entirely: a lightweight in-place refresh
on `VolumeNode` that memcpys the changed voxels into the already-allocated vtkImageData and
calls `Modified()`, skipping the deep copy, the reallocation, the TF reset, and the state
writes. Sketch: `VolumeNode::refreshData()` + a `pycvc_gl` binding, a few lines against the
existing `updateImageData` guarded on dims/type being unchanged. With it, per-frame volume
updates are cheap and striding becomes unnecessary — this is what makes the water and cloud
genuinely free rather than merely rationed.

### 2. Shadows: off, then back ON once the actors came down (done)

Shadow mapping more than doubles the render: `vtkShadowMapBakerPass` re-renders every actor
from every light, so its cost scales with the ACTOR count — at 3776 actors × 2 lights it
turned a live orbit into a slideshow, so shadows were made off-by-default. Fix #3 then cut the
actor count to ~140, which is exactly the quantity the shadow baker scales on, so **shadows
are now ON by default again** and affordable. A live toggle (`…lsystem_forest.shadows` = 0/1)
remains for the last few fps. The order matters: shadows became affordable *because* the actor
count came down, not on their own.

### 3. Reduce actor count — DONE (route C), ~13× on the measured box

3776 → ~140 actors, roughly one per tree, was the real bottleneck, and it is now implemented
as **route C** (the merged-mesh CPU re-pose from the design below). Each tree is built as ONE
merged wood actor (+ ONE merged needle actor); the per-module hierarchy still exists as data,
so the wind cascade is unchanged — it is applied to the merged vertices each frame via a new
`GeometryNode::updateVertices` fast path in cvcGL (points-only, no cell/normal rebuild) instead
of one `setTransform` per module. Per-vertex wood colour is kept natively (no shader).

Measured on this box (offscreen, 900×600), against the 3781-node baseline:

| | baseline (3776 actors) | route C (145 actors) |
|---|---:|---:|
| render only (shadows off) | 225 ms | **16.5 ms** |
| live frame, shadows **off** | ~1.6 fps | **~18 fps** |
| live frame, shadows **on**, bake every frame | (hopeless) | ~11 fps |
| live frame, shadows **on**, bake 1/3 (default) | (hopeless) | **~17.5 fps** |

The trees drop from ~194 ms of render (3776 draw calls) to ~0.4 ms; the terrain heightfield is
now the largest single render cost.

**Shadows are on by default and they render** (I was briefly wrong that they didn't — the shadow
baker prints benign `ReleaseGraphicsResources` teardown warnings offscreen, but it casts
correctly; verified by a shadows-on/off frame diff — 49.5% of pixels change — and a screenshot
with a visible cast tree-shadow on the ground). Two follow-on fixes made them nearly free:

* **Re-bake the shadow map on a stride, not every frame.** Because the wind moves geometry every
  frame, `vtkShadowMapBakerPass` re-bakes every frame — that per-frame bake, not the shadow
  sampling, is the whole cost (~1.69× the frame). But shadows from a slow ~5 s sway barely change
  frame to frame, so cvcGL's `setShadowUpdateInterval(N)` (a `StridedShadowBaker`) bakes every Nth
  frame and reuses the map between. At `SHADOW_STRIDE=3` the live frame goes from ~11 fps back to
  **~17.5 fps** — 93% of the shadow cost gone, right next to the shadows-off rate — and the shadows
  still track the sway (just every third frame). It is the shadow analogue of fix #1's `VOL_STRIDE`.
* **numpy-direct `updateVertices`.** The per-frame re-pose passes the merged vertex buffer straight
  through the buffer protocol instead of `.ravel().tolist()`, dropping the per-vertex Python-float
  allocations (~11× cheaper per call) — this is what lifted the shadows-off rate to ~18 fps.

At 3776 actors the per-frame re-bake — every actor, every light, every frame — was out of the
question, which is the whole point: **shadows are affordable because the actor count came down**,
and the stride then makes them nearly free.

Routes B (glyph instancing) and A (shader skinning) remain documented alternatives below. B has
since been **prototyped end-to-end** (2026-08-18) — the wood-colour loss is recovered with a
custom shader, shadows are safe, and it is measured against route C. The verdict is *do not build
it yet*: it optimises tree render, which route C already made a non-bottleneck. See **"Route B:
executed spike results"** below for the numbers; the source-read reasoning that follows is left
intact as the pre-spike record. Route C reaches the same actor-count number, keeps the colour for
free, and shipped.

## Design detail: the three routes (route C shipped; B/A for reference)

Fix #3 wants ~3776 actors -> one `vtkActor` per tree, without losing per-module wind. Read
this the way the rest of the doc asks to be read: everything below marked "from source" is a
**source read of `libcvc@origin/master`, not an executed spike**. The doc's own rule — a quoted
number is not evidence — applies to a quoted source path too. Treat the VTK-9.5 behaviours as
*load-bearing hypotheses to be proven by the spike at the end*, not as settled facts.

### Two facts that bound every option (from source)

* **The cascade is CPU-side already.** `GraphicsNode::updateTransform()` composes
  `vtkMatrix4x4::Multiply4x4(parent->m_worldMatrix, m_transform, m_worldMatrix)` top-down and
  hands each module's actor a `vtkTransform` via `SetUserTransform` (`applyWorldTransformToProps`).
  "Wind accumulates down the tree" is that recursion re-multiplying a swayed ancestor's whole
  subtree. There is no free VTK assembly hierarchy to lose — any flattened form must reproduce
  that matrix product itself (a transpose/aliasing bug in the reimplementation mis-poses every
  child, so match `Multiply4x4`'s semantics and row/col convention exactly).
* **Wind re-poses the trunk, so the whole tree moves.** `SWAY_LEVELS=2` re-poses trunk-side
  modules whose subtree is the entire tree, so in any merged form effectively all of a swaying
  tree's vertices move each frame. The "only pose 2 levels, children inherit for free" property
  does not survive flattening — but the work it saves is trivial (Mflops), so this costs nothing.

One more: `updatePolyData` is mode-switched (`SetVerts`/`SetLines`/`SetPolys`, one primitive
kind per node), which is *why* the scene is 1888 wood nodes + 1888 needle nodes. A single
`vtkPolyData` can legally hold polys **and** lines, but `ensureNormals()` is the catch (below).

### The three routes

**A. Shader-replacement skinning (`vtkOpenGLPolyDataMapper`).** Merge each tree into one
polydata; tag every vertex with a module index (point-data array); upload the tree's ~21 module
matrices as a `mat4[]` uniform palette; inject a vertex-shader prologue that applies
`palette[moduleIndex]` before the model-view matrix. Hits ~70 actors, wind intact.
*Source suggests* the key fast-path works — value-only uniform updates don't recompile the
shader (`vtkOpenGLUniforms::SetUniformValue` bumps only `Parent->Modified()`, the rebuild gate
stays closed, `UpdateShaders` still pushes custom uniforms each frame) — but this and the anchor
strings it edits (`//VTK::PositionVC::Dec`, `MCVCMatrix * vertexMC`, the normal line) are exactly
what a point release can shift, and `vtkShaderProgram::Substitute` **silently no-ops on a miss**,
so a wrong anchor degrades to a rest-pose render with no error. Complexity: highest — cvcGL has
never written GLSL or touched `vtkShaderProperty`, and the SWIG-director node it needs runs into
the pycvc-family SWIG-runtime landmine (must build with pycvc's 4.2.0 runtime, not cvcpkg's
4.4.1, or cross-module type sharing breaks). Two correctness traps the first draft glossed:
**normals under per-axis scale need the inverse-transpose of the palette's 3×3, not the palette**
(multiplying a normal by a non-uniform-scale matrix bends lighting); and `vtkPolyDataNormals` on
the merged rest mesh **smooths normals across trunk↔branch seams**, which are then wrong at the
joints under sway unless normals are split there (which defeats "run `ensureNormals()` unchanged").

**B. Glyph instancing (`vtkGlyph3DMapper`) — stock VTK does the skinning.** *From source*
(`vtkOpenGLGlyph3DMapper::RebuildStructures` / `vtkOpenGLGlyph3DHelper`), the mapper builds a
per-instance `T·R·S` model matrix + normal matrix from position/orientation(quaternion)/scale
arrays and issues one `glDrawElementsInstanced` per source. Every wood segment and needle is
rigid + axis-aligned scale, so it fits. Instanced *segment-level* → **one cylinder source over
all 1888 segment bases + one line source over all needles = 2 actors, 2 draw calls for the whole
forest**, well past the target, wind kept by rewriting the affected instances' arrays. Zero GLSL.
The verified loss: the glyph VS overwrites `//VTK::Color::Impl` with a per-instance color, so the
wood's **per-vertex** light-top/dark-ring gradient (`_CYL_COLORS`) collapses to one color per
segment — a product call that must be signed off before B is real (recovering it is a texture
path, new work). `vtkPointGaussianMapper`/`vtkTensorGlyph` don't fit and are rejected.

**C. CPU re-pose of one merged polydata per tree (the honest baseline).** Merge each tree in
numpy (`np.concatenate` + index offsets — the example already builds the per-module arrays);
each frame transform the swaying tree's vertices *and normals* (again: **inverse-transpose for
the normals**) in numpy and re-upload. Wind kept. Actor count: 70 for wood, plus needles as a
second per-tree LINES node (→140) or one scene-wide needle actor (→71); a *true* 70 needs one
node carrying both polys and lines, which trips `ensureNormals()` (next paragraph). Cost is
provably below the actor-count wall — but note the draft's "0.6 MB/frame, points-only" is
optimistic: swapping `GetPoints()` + `Modified()` bumps the polydata MTime, and
`vtkOpenGLPolyDataMapper::BuildBufferObjects` then **re-packs and re-uploads the whole
interleaved VBO** (points+normals+colors+tcoords), not just the moved points — there is no
partial-buffer update. Still cheap at ~70 trees, but measure the real bytes/frame. Complexity:
lowest — one contained points-only `GeometryNode` fast path plus numpy the example already writes.

**Mixed tris+lines caveat (the "true 70" that is really 140).** `ensureNormals()`'s own comment
says `vtkPolyDataNormals` "strips the lines out entirely," and it only runs the filter when
`GetNumberOfPolys() != 0`. A merged wood+needle polydata has polys > 0, so the filter runs and —
per that comment — drops the needle lines. So one actor holding both needs either normals
precomputed per-cell-type with the filter skipped, or two nodes per tree (~140 actors). **140 is
the safe, already-feasible number and is still a ~27× cut; do not quote 70 as if it were free.**

| route | actors | wind | draw calls | new surface | headline risk |
|---|---|---|---|---|---|
| A. shader skinning | ~70 | yes | ~70 | large (GLSL + `vtkShaderProperty` + director node) | silent substitution; shadow baker; inverse-transpose normals |
| B. glyph instancing | 2–70 | yes | 2 | moderate (`GlyphNode` + arrays) | ~~per-vertex wood color lost~~ **recovered by shader (spiked)**; wood *lighting* ~15% off (glyph normals under non-uniform scale) |
| C. CPU-repose merge | 70–140 | yes | 70–140 | small (points-only update) | full-VBO re-upload; mixed tris+lines |

### Recommendation: stage it, do not lead with the shader

1. **Ship C first (plan of record).** The measured bottleneck is *actor count*; C hits it
   (3776 → ~140) with the smallest new surface, no GLSL, no coupling to VTK shader internals,
   and a per-frame cost below concern. It is the change most in proportion to the evidence.
2. **Escalate to B** only if draw-call count or upload ever becomes the real bottleneck (it is
   not today). B collapses the forest to ~2 draw calls and gets the GPU-side skinning fix #3
   imagined *for free from stock VTK*, at the cost of a `GlyphNode` and — now that it is spiked —
   a ~15% wood-lighting delta rather than an outright colour loss (the gradient is recovered by a
   shader; see "Route B: executed spike results"). The spike confirms this escalation is not needed
   today: B's render win is real but redundant while the per-frame pose dominates.
3. **Keep A as a documented fallback, not the first build.** It is plausibly feasible, but it is
   the most code and the most fragile and buys nothing over B on the GPU.

This refines fix #3 rather than contradicting it: skinning does preserve wind — it is just not
the only, nor the cheapest, way to the actor-count number, and it is the most expensive.

### Minimum cvcGL / pycvc API

**For C (the first move) — genuinely small:**
* `GeometryNode::updateVertices(const cvc::geometry&)` — a **topology-preserving fast path**:
  overwrite `m_polyData->GetPoints()` (and pre-transformed normals) and `Modified()`, skipping
  cell rebuild, `ensureNormals()`, and colour/TCoord regeneration. This is the whole
  performance-critical addition; do NOT route the re-pose through `setGeometry()` (it re-runs
  `vtkPolyDataNormals` over the tree ~23×/frame).
* `pycvc.geometry.set_vertices(numpy)` + a `node.update_vertices(geom)` binding, marshaled to the
  owner thread like `setColor` — but note this pushes a full array *every frame*, not occasionally,
  so measure the marshaling cost and guard the buffer swap against an in-flight upload.
* Optional (for true-70): teach `updatePolyData`/`ensureNormals` to keep both `SetPolys` and
  `SetLines` and compute normals without the line-stripping filter.

**If/when A is built (fallback):** a `SkinnedGeometryNode : GeometryNode` with
`setBoneCount`/`setBoneParents`/`setBoneLocal`(+`setBonesLocal`/`resetBones`), a `setGeometry`
override that reads a per-vertex `ModuleIndex` array (reuse `cvc::geometry`'s unbound
`functions()` scalar channel), calls `MapDataArrayToVertexAttribute` +
`SetVBOShiftScaleMethod(DISABLE_SHIFT_SCALE)` on a `vtkOpenGLPolyDataMapper::SafeDownCast` of the
mapper (**hard-fail on a null downcast, don't silently skip skinning**), installs the shader
replacements, and uploads the palette via `GetVertexCustomUniforms()->SetUniformMatrix4x4v(...)`.
Promote `m_actor`/`m_mapper`/`m_polyData` to `protected`. Plus the `pycvc_gl.i` director/factory
wiring — subject to the SWIG-runtime constraint above.

**For B:** a separate `GlyphNode(GraphicsNode)` wrapping `vtkGlyph3DMapper` + N source polydatas
with numpy `set_instances(positions, quats, scales[, source_index])` / `update_instances(...)`.

All three **help** #193 (the per-frame path stops calling `setTransform`, so no `getState("matrix")`
write/echo) and #189 (fewer traversed props). Only A's shadow-pass behaviour is a genuine open
risk against #194's shadow passes.

### Effort, and the spike that must come first

*(The B spike below has since been executed — see "Route B: executed spike results". The checklist
here was written for it and for A; B's items 3–4 came back: normals under non-uniform scale ARE off
under glyph instancing, shadows move correctly with data-driven sway. A remains un-spiked.)*

Rough order of magnitude, explicitly *pre-spike* and optimistic: C ~1–2 days (for the 140-actor
version with a full-VBO re-upload deemed acceptable), B ~3–5 days (if the wood-color loss is
pre-accepted and the quaternion convention lands first try) — **now ~1–2 days post-spike, both
unknowns resolved**, A ~1–2+ weeks with the residual risk concentrated in the shadow pass and
substitution fragility. **Do not trust the A estimate until the shadow spike lands.**

Because C needs none of A's unknowns, it can land in parallel. The spike exists to prove the
things source-reading cannot, and it must use geometry representative on every risky axis at once
— **normals + lighting + shadows + a per-axis-scaled module + a line primitive**, not a bare quad:

1. Dump the *actual* generated vertex shader for real geometry with normals, lighting and shadows
   all on, and grep for each anchor string, asserting each occurs the expected number of times
   (VTK emits `MCDCMatrix * vertexMC` on the common path; `MCVCMatrix` only when view-space
   position is needed — the anchor may be absent or appear more than once).
2. Prove a `mat4[N]` custom uniform round-trips to the GPU and indexes correctly in GLSL (not
   merely that the call compiles), and instrument the shader compile count across many palette
   updates to confirm **zero** recompiles after frame 1. Keep a `vtkTextureObject`+`texelFetch`
   palette fallback live until this holds. Query real `GL_MAX_VERTEX_UNIFORM_COMPONENTS` headroom
   *with shadows+lighting on* (the 1024 floor is shared, not free).
3. Sway a per-axis-scaled module and verify lighting is correct (inverse-transpose normals) and
   seams don't crack; run a merged polys+lines tree through the real `ensureNormals()` and confirm
   whether lines are dropped; assert the mapped `ModuleIndex` length == point count and reads back
   right per vertex.
4. With shadows on, sway a module and confirm the **cast shadow moves** — i.e. the baker's
   depth-only program variant received the vertex replacements *and* the current palette. This is
   the single most likely silent break (canopy sways, shadow baked from rest pose). If it fails,
   that is the signal to commit to C+B rather than fight the baker.

## Route B: executed spike results (2026-08-18)

Everything above about route B is a source read. This section is the **executed spike** —
`scripts/bench/glyph_route_b/` — run against real VTK 9.5 on an NVIDIA GTX 1650. It is a pure-VTK
harness (no cvcGL): the forest tree grammar and per-module wind cascade are reproduced in numpy and
**verified against the exact `lsystem_forest.py` math to 7e-15** (`geom_check.py`), so every
route-B-vs-C difference below is *rendering*, not geometry. Re-run it rather than trusting these
numbers; the doc's own rule still applies.

**Verdict: B is real and de-risked, but not worth a cvcGL `GlyphNode` today.** It reproduces route
C's look almost exactly and slashes render cost, but the cost it slashes — tree *render* — is one
route C already brought under budget, while the cost that actually dominates now (the per-frame CPU
pose) is untouched by B. Build it only if the forest grows several-fold or a scene becomes
draw-call bound.

### The wood-colour loss is recoverable (the blocker is gone)

The glyph drops the source's per-vertex colours, as predicted. It is recovered with a custom shader
on the **actor** (`vtkGlyph3DMapper` has no `GetShaderProperty`; the glyph helper reads the actor's),
and the fix hinges on an anchor the first attempt got wrong:

* Injecting the wood colour at `//VTK::Color::Impl` **fails to compile** ("undefined variable
  ambientColor"). Actor shader-property replacements are applied *before* the mapper's own colour
  substitution, so replacing `//VTK::Color::Impl` deletes the very block that *declares*
  `ambientColor`/`diffuseColor`. Inject instead at **`//VTK::Normal::Impl`** — after the colour
  declarations, before the light math — where those variables are in scope. (Dump the real generated
  shader with `dump_glyph_shaders.py`; do not guess anchors.)
* The radial light-core/dark-rim gradient is recomputed in-shader from the glyph-local radius
  `length(vertexMC.xz)` (unit-cylinder radius 1), so it is independent of source tessellation.
* The glyph **source must be the forest's `_CYL`** (explicit cap-*centre* vertices at radius 0 + a
  triangle-fan cap). `vtkCylinderSource`'s cap is one n-gon over the rim with no centre vertex, so it
  renders *no* light core at all — an easy false negative when spiking.

`test_wood_gradient.py`: flat glyph shows 0 % light-core pixels; the shader recovers 23.5 %, matching
route C's 22.2 %. Under flat lighting the wood albedo matches route C to <1/255.

Orientation convention (nailed in `orient_test.py`, previously an open risk): `vtkGlyph3DMapper`
quaternion array order is **(w, x, y, z)**; per instance `pos = (W·m)[:3,3]`, `R = (W·m)[:3,:3]`
(column-vector), `scale = (seg_rad, seg_len, seg_rad)`.

### Shadows are safe — the doc's "single most likely silent break" does not apply to B

The feared failure (canopy sways, shadow baked from rest pose) is a *route A* risk, because A moves
vertices with a shader palette the depth pass might not receive. **B animates via the instance
arrays (data), so the same glyph input drives the baker's depth pass for free** — a swayed instance's
cast shadow moves (`shadow_test.py`: ~13.7k ground-shadow px change between two sway phases, on par
with route C). And the custom wood shader **compiles and renders correctly under `vtkShadowMapPass`**
(`shadow_isolate.py`: 0 compile failures for the wood glyph). The needle *lines* do throw shadow-pass
shader errors — but so do route C's needle-line polydata (more, in fact) — so that is a VTK
lines-under-shadow quirk, not a route-B problem.

### Quality: needles perfect, wood colour perfect, wood *lighting* ~15 % off

Same forest, same camera/lights, route B vs route C (`compare_quality.py`): silhouettes are
pixel-identical (foreground IoU 1.000), overall image diff 3.9/255. Needles are effectively perfect
(0.9/255). The one real gap is **wood lighting**: route B's lit trunks come out ~15 % darker in
aggregate. The albedo/gradient is identical (matches under flat light); the difference is normals.
`vtkGlyph3DMapper`'s per-instance normal matrix mishandles the beveled cap/ring normals under
**non-uniform per-axis scale** (rad ≠ len) — the same inverse-transpose subtlety flagged for route A,
here inside the glyph helper with **no clean knob** (hard/split normals make it worse, not better;
`normal_test.py`). Acceptable for a demo, but it is a genuine fidelity cost route C does not pay.

### Performance: B's render is O(1) in actors; but the per-frame *pose* dominates and B doesn't touch it

`bench.py`, 70 trees (140 actors) → collapsed to 2 glyph draw calls, 900×600, GPU-synced. Trees only
(no terrain/volumes), so absolute fps is higher than the whole-scene figures earlier — the point is
B-vs-C on one identical harness.

| | route C (140 actors) | route B (2 glyphs) |
|---|---:|---:|
| per-frame CPU pose (`update`) | ~68 ms | ~66 ms |
| render (GPU-inclusive) | ~48 ms | ~18 ms |
| sustained | ~9.3 fps | ~14.7 fps |

Two things the numbers say plainly:

* **B's render advantage is real and grows with scene size.** Route B's render is ~flat in actor
  count (2 draw calls): 16→19 ms as trees go 70→210. Route C's scales with actors: 49→143 ms. So the
  render speedup grows **3.0× → 7.6×** over that range. This is exactly fix #3's draw-call argument,
  now measured.
* **But render is no longer the bottleneck — the CPU pose is.** The per-frame wind cascade is ~66 ms
  for *both* routes (it is the same numpy either way), so it dominates the frame and route B does not
  reduce it. Net sustained gain is only ~1.6×, not the ~11× the raw draw-call collapse suggests. Once
  actor count is down, the next real target is the pose, not the draw path — and that is a shared
  cost neither B nor the actor-count fix addresses.

The honest reading: route C already took tree render from ~194 ms to ~0.4 ms of the *scene*. B would
take that ~0.4 ms to ~0.2 ms — a real multiple, an irrelevant absolute — while adding the wood-lighting
regression above and a new `GlyphNode` + instance-array binding surface. **Escalate to B only when
draw calls actually bite: a forest of several hundred trees, or a genuinely draw-call-bound scene.**
The spike removes B's unknowns (colour recovered, shadows safe, quaternion nailed), so if that day
comes it is ~1–2 days, not the pre-spike 3–5.

### Bonus finding: procedural bark is a route-B-native capability

Because B already renders the wood through a custom shader, bark (fragment **bump** and/or vertex
**displacement**) is a shader edit, and it is cheap — `bench_bark.py`, wood-only, 70 trees, vsync
disabled (offscreen renders otherwise floor at 1/60 s):

| variant | source verts/seg | wood render |
|---|---:|---:|
| base (pentagon) | 12 | 0.56 ms |
| + fragment bump | 12 | 0.45 ms |
| tessellated (24×12), no bark | 314 | 3.56 ms |
| + vertex displacement | 314 | 3.14 ms |
| bark: displacement + bump | 314 | 3.47 ms |

* **Fragment bump mapping is essentially free** (per-fragment ALU; ~0 ms here) — bark *texture* for
  nothing. Done with the surface-gradient method (`dFdx`/`dFdy` of a procedural height and view
  position), so no precomputed tangents.
* **Vertex displacement is also free — but it needs tessellation** (a pentagon has nothing to
  displace). Going 12→314 verts/segment costs +3 ms of GPU render — the whole cost is the extra
  vertices, not the displacement. On the pose-dominated ~68 ms frame that is **+4 %**.
* Crucially the tessellated source is uploaded **once** and the per-frame instance arrays are
  **unchanged**, so bark adds **zero** per-frame CPU cost. The equivalent on route C would bake a
  ~25× denser mesh into every tree's merged buffer and re-upload all of it every frame (route C's
  cost *is* that per-frame vertex upload) — so bark is a natural fit for B and a poor one for C. If
  bark ever becomes a requirement, that flips the B-vs-C calculus toward B.

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
