# Route B (glyph instancing) spike — lsystem_forest wood

Executed prototype of **route B** from [`docs/RENDER_PERFORMANCE.md`](../../../docs/RENDER_PERFORMANCE.md):
render the forest with `vtkGlyph3DMapper` (one unit cylinder instanced over every
wood segment, one needle star over every leaf) instead of route C's one merged
actor per tree. Answers the questions the doc left as source-read hypotheses:
recover the per-vertex wood colour a glyph drops, check the shadow interaction,
and compare fps/quality against route C.

Pure VTK 9.5 — **no cvcGL/pycvc needed**. The forest tree grammar/cascade is
reproduced in numpy (`forest_geom.py`) and verified against the real
`scripts/examples/lsystem_forest.py` math to 7e-15 (`geom_check.py`), so any
route-B-vs-C difference is rendering, not geometry.

## Running

Needs a Python with VTK 9.5 (e.g. a cvcpkg `deps-live` prefix):

```sh
P=../../../../deps-live            # a prefix with bin/python3.11 + VTK 9.5
export LD_LIBRARY_PATH=$P/lib PYTHONPATH=$P/lib/python3.11/site-packages
$P/bin/python3.11 bench.py 70
```

Offscreen rendering is **vsync-locked to 60 Hz** on this box — sub-16.6 ms render
costs read as a flat 16.667 ms until you disable it:

```sh
export __GL_SYNC_TO_VBLANK=0 vblank_mode=0    # required for bench_bark.py to be meaningful
```

## Files

**Modules**
- `forest_geom.py` — forest tree grammar + module cascade in numpy; builds both
  route-C merged buffers and route-B per-instance (pos, quat, scale) arrays.
- `render_forest.py` — VTK actors for both routes in one shared scene (lights,
  camera, optional `vtkShadowMapPass`).
- `wood_shader.py` — the recovered wood-gradient glyph shader (inject at
  `//VTK::Normal::Impl`, not `//VTK::Color::Impl`) + the forest `_CYL` source.
- `bark_shader.py` — procedural bark: fragment bump (surface-gradient method) and
  vertex displacement, on a tessellated cylinder source.

**Verification**
- `geom_check.py` — route B (pos/quat/scale) reconstructs route C to 7e-15.
- `orient_test.py` — nails the glyph quaternion convention: order **(w,x,y,z)**.
- `test_wood_gradient.py` — the wood gradient recovers (light-core % matches C).
- `normal_test.py` — glyph per-instance normals under non-uniform scale deviate
  from a correct mesh (the ~15 % wood-lighting delta; no clean fix).
- `dump_glyph_shaders.py` — dump the real generated glyph fragment shader.

**Comparison / benchmarks**
- `compare_quality.py` — route B vs C side-by-side (needles pixel-perfect, IoU 1.0).
- `bench.py [n_trees]` — fps/cost B vs C; render is O(1) in actors for B, O(actors)
  for C. The per-frame CPU pose dominates both and B doesn't reduce it.
- `bench_bark.py [n_trees]` — cost of bark: fragment bump ≈ free; vertex
  displacement free but needs tessellation (+~3 ms GPU, +~4 % of the frame).
- `shadow_test.py` / `shadow_isolate.py` — the wood shader survives the shadow
  baker (0 errors) and swayed instances cast moving shadows; needle-line shadow
  errors are a VTK lines+shadow quirk common to both routes.

## Bottom line

Route B works and is de-risked (colour recovered, shadows safe, quaternion nailed),
but it optimises tree render — which route C already made a non-bottleneck. Build a
cvcGL `GlyphNode` only if the forest scales past ~200–300 trees or a scene becomes
draw-call bound. See `docs/RENDER_PERFORMANCE.md` for the full write-up.
