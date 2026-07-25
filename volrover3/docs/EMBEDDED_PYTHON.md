# volrover3 Embedded Python Environment — Design

**Status:** design · initial branch `feat/embedded-python-interpreter`
**Goal:** bake a **first-class** CPython interpreter into volrover3 (not a plugin) so scripts can drive
the *running* application — its `cvc::app` state tree, its scene graph, and (see §7) its **Qt UI** — through
a robust, volrover3-specific binding surface. Conceptually modeled on the **verlihub `plugins/python`**
embedding, adapted to Python 3 + SWIG + the no-singleton architecture.

---

## 0. Requirements (verbatim intent)

1. **Model the embedding on verlihub's python3 plugin** — one interpreter, disciplined GIL bracketing on every
   crossing, name-based hook discovery, per-script exception isolation, a load/reload/unload control surface.
2. **Bake it in** as a core subsystem owned by the app — *not* a loadable plugin.
3. **Single interpreter by default**, multi offered but off — because pycvc / pycvc_gl / VTK / numpy call into C
   libraries that are **not sub-interpreter-safe** (single-phase-init C globals alias across sub-interpreters).
4. **SWIG bindings against volrover3-specific host objects.**
5. **A host-control module** that controls the running application's **state** and **delivers the live
   `cvc::app`** to pycvc scripts (so they operate on the *running* app, not a fresh one).
6. **No singletons** — own the `cvc::app` explicitly and thread it; remove `volrover3::app()` and
   `AppState::instance()`. (Authorized: edit volrover3 as needed.)
7. **Qt UI from Python** (§7) — expose the live `QMainWindow` so scripts can manipulate, extend, or **outright
   replace** the UI with PySide6 widgets defined in Python.

## 1. Roadmap alignment

This subsystem is the **native runtime form** of the embedded-interpreter goal already on the cvcpkg roadmap —
it is not new scope, it is the concrete first step:

- **CVCPKG-ROADMAP Phase 19 (Application Packaging & Desktop Delivery)** — the "Embedded-Python single binaries"
  note and `cvcpkg bake`. The interpreter designed here is later *baked* (single-interpreter, inittab-registered
  host module, embedded entry script — the exact stitching in
  [`libcvc-deps/docs/roadmap/static-single-binary-python.md`](../../../libcvc-deps/docs/roadmap/static-single-binary-python.md) §2).
- **static-single-binary-python.md Case B** is literally *"WASM VolRover with an embedded interpreter (libcvc +
  vtk-python)"* — this design's single-interpreter + host-module-inittab shape is what that bake consumes.
- **Phase 20** lists **verlihub** as a featured recipe — the app we model the embedding on.
- The **Python/Qt binding + hermeticity** architecture doc (engagement-docs `modernization/`) already scopes the
  `pycvc` / `pycvc_gl` / `pycvc.qt` surface and the **PySide6 + Shiboken6** hermetic recipes as the "linchpin
  heavy lift" — §7 here builds directly on that.

## 2. What we copy from verlihub (and what we drop)

Grounded in `/home/joe/src/verlihub/plugins/python`:

**Copy (the load-bearing patterns):**
- **Exactly-once-per-process interpreter lifecycle.** verlihub brings CPython up once (`w_Begin`:
  `PyEval_InitThreads` + `Py_Initialize`, then parks the GIL via `PyThreadState_Swap(NULL)` +
  `PyEval_ReleaseLock`, `wrapper.cpp:1275-1285`) and finalizes once (`w_End`, `wrapper.cpp:1302-1304`). **Never
  per-script.**
- **Disciplined GIL bracketing on every C++↔Python crossing** (`w_CallHook` Acquire/Release,
  `wrapper.cpp:1633/1957`; Python→host releases + reacquires the same state, `Call()` `wrapper.cpp:554-561`).
- **Name-based hook discovery** — scan a script's namespace for well-known functions and cache a presence
  bitmap (`wrapper.cpp:1475-1482`); scripts subscribe just by defining `OnFoo`, no register API.
- **Per-script exception isolation** — catch, `PyErr_Print`, clear, return a safe default, keep going
  (`wrapper.cpp:1903-1955`) so one bad script never takes down the host.
- **A load/reload/unload + filesystem-as-registry control surface** (`cconsole` `!pyload/!pylist/!pyunload/
  !pyreload/!pylog`) → ported to a Qt dock (§4.4).

**Drop (chat/plugin-specific or Py2-era):**
- The **MULTI sub-interpreter-per-script** model (`Py_NewInterpreter`, `wrapper.cpp:1363`) — unsafe for
  C-extension packages (§6). We default to **single**.
- The **flat scalar/string marshalling ABI** (`w_Targs`, `wrapper.h:182-185`) — **SWIG proxy objects** replace
  it entirely (the whole point of wrapping real volrover3 objects).
- The **`cpiPython::me` process-global singleton + `vh.myid` reflection identity** (`cpipython.cpp:60,68`;
  `wrapper.cpp:340-362`) — replaced by an **injected wrapped host handle** (§4/§5), and consistent with the
  no-singleton rule.
- The **separate `dlopen(RTLD_GLOBAL)` wrapper `.so`** (`cpipython.cpp:102-124`) — needed only because the
  verlihub plugin is itself a `.so` that must not link libpython. volrover3 is an **executable**, so it links
  libpython directly and uses **`-Wl,--export-dynamic`** to make libpython symbols visible to imported
  C-extensions (numpy/vtk) — the executable-world equivalent of `RTLD_GLOBAL`.
- Everything is **Python 2 C-API** (`Py_InitModule`, `PyString/PyInt`, `PyEval_ReleaseLock`) → ported to the
  Py3 C-API + SWIG-generated init.

## 3. No-singleton app ownership (a volrover3 refactor)

Today the process app is a **Meyers singleton**: `cvc::app& volrover3::app()` returns a function-local `static`
(`volrover3_app.cpp:4-7`), used at **155 sites** across 12 files; `AppState::instance()` (`AppState.h:15`) is a
second singleton, **45 sites**. Both violate the no-singleton rule and force the app-delivery workaround
(aliasing a static with a no-op deleter).

**Target model — own it, thread it:**
- `main()` (or a small `volrover3::Application` bootstrap object) constructs **one** `std::shared_ptr<cvc::app>`
  and hands it to `MainWindow`.
- `MainWindow` owns the `shared_ptr<cvc::app>` and threads it into: the **cvcGL** `SceneGraph` (via #136's
  **injected-app ctor** `SceneGraph(cvc::app&, prefix)`), `AppState` (now a plain object owned by `MainWindow`,
  constructed with the app), the SceneNodes (already thread the app as `SceneNode::_ctx`, `SceneNode.h:16,21` —
  good), and the new `EmbeddedInterpreter`.

**Extract ALL 3D/VTK into cvcGL — volrover3 becomes a thin Qt + Python shell (companion to Phase 0).**
The directive: *all 3D-graphics and VTK code used by volrover3 lives in cvcGL and is wrapped by pycvc_gl;
remove the redundant copies from volrover3.* Concretely, volrover3 today carries a **private fork** of the
entire scene graph — 9 classes (`SceneGraph`, `GraphicsNode`, `GeometryNode`, `VolumeNode`, `SceneNode`,
`GridNode`, `AxisNode`, `BBoxNode`, `NullGraphicNode`, under `inc/volrover3/`) — and links only `cvc::cvc`.

**cvcGL already has all 9** (`inc/cvc/gl/*.h`), and is the **superset**: byte-diffs are 0–2 lines for most
(only an include path differs), and cvcGL's `SceneGraph` (+21 lines) already carries the #136 injected-app ctor
+ extras; `SceneNode` similar. So this is overwhelmingly **delete-the-fork-and-repoint**, not a reconciliation:

1. **Delete** volrover3's 9 forked node headers/impls.
2. **Repoint** includes `<volrover3/XNode.h>` → `<cvc/gl/XNode.h>` and link **`cvc::cvcGL`** (the `cvcgl`
   package) alongside `cvc::cvc`.
3. **Reconcile** the few real deltas (adopt cvcGL's injected-app `SceneGraph` ctor; port any volrover3-only
   node tweaks *up into cvcGL* first if any survive the diff — the `SceneGraph`/`SceneNode` +lines are the only
   candidates).
4. **`VTKRenderWidget`** stays in volrover3 (it is a Qt widget) but becomes a **thin host**: it hands its
   `vtkRenderer` to cvcGL's `SceneGraph::setRenderer` and drives `processEvents()`/render — no scene logic of
   its own. Same for the camera/transfer-function *rendering* logic (the 3D parts move to cvcGL; the Qt dialogs
   stay).

**Result:** the C++ app, `pycvc_gl` scripts, and the state tree share **one** scene-graph engine over **one**
app. A script's `pycvc_gl.Scene(app)` and the running UI's scene are literally the same objects — so §7's
"Python manipulates/replaces the UI" is matched by "Python manipulates/replaces the **3D scene**" through the
same cvcGL that C++ uses. volrover3 C++ shrinks to: `QApplication` + the Qt widget/dialog shell + `AppState` +
the embedded interpreter + the app bootstrap. Everything 3D is cvcGL (reusable, wasm/bake-able, Python-wrapped).
- `volrover3::app()` and `AppState::instance()` are **deleted**; call sites take the app/appstate from their
  owner (widgets from `MainWindow`, nodes from the SceneGraph, dialogs constructed with a reference).

This is staged (§10 Phase 0) because it touches ~200 call sites; it is a mechanical "thread the handle that is
already reachable" refactor, and it is a prerequisite for a clean (non-aliased) app delivery in §5.

## 4. Components (new C++, all compiled into `volrover3_lib` so they are unit-testable)

`volrover3_lib` is the existing `STATIC` library of all sources minus `main.cpp` (`CMakeLists.txt:415`); the
tests link it. The interpreter subsystem compiles **into** it.

1. **`volrover3::EmbeddedInterpreter`** (`inc/volrover3/EmbeddedInterpreter.h` / `.cpp`)
   Owns the CPython lifecycle (`Py_InitializeEx(0)` once on the UI thread, `Py_FinalizeEx` once in the dtor), the
   GIL policy (`m_mainState = PyEval_SaveThread()` to park after init), script load/exec/unload, and the injected
   host handle. One instance, owned by `MainWindow` as `std::unique_ptr<EmbeddedInterpreter> m_interp`,
   constructed **right after** `m_sceneGraph` (`MainWindow.cpp:71`) and destroyed in `~MainWindow` **before**
   any exit-time static teardown. Holds `PyThreadState* m_mainState`, `std::shared_ptr<PyHost> m_host`, and a
   `std::map<std::string, PyObject*>` of per-script module namespaces (single-interpreter isolation).

2. **`volrover3::PyHost`** (`inc/volrover3/PyHost.h` / `.cpp`) — the host-control facade SWIG wraps and scripts
   drive. Constructed with the injected `std::shared_ptr<cvc::app>` + `std::shared_ptr<SceneGraph>` + (later) the
   `QMainWindow*`. **Not a singleton.** Surface:
   - `std::shared_ptr<cvc::app> app()` — **the delivery method**; returns the owned app handle verbatim (§5).
   - `std::string state_prefix() const` → `"volrover3"`.
   - Scene control: `scene()`, `request_render()` (posts via `SceneGraph::postEvent`, drained by the 16 ms
     `QTimer`, `VTKRenderWidget.cpp:25`), `list_nodes()`.
   - Typed host-state helpers over the same subtree `AppState` writes (camera fov/position, world bounds with
     the `readOnly(true)` toggle gotcha, show-fps) — ergonomic sugar; a script can equivalently call
     `pycvc.state_set(host.app(), "volrover3.camera.fov", …)`.
   - (§7) `quintptr main_window_ptr()` — the live `QMainWindow*` as an int, for the PySide6/Shiboken bridge.

3. **`vrhost` SWIG module** (`bindings/vrhost.i`) — `%module vrhost`, **`%import "pycvc.i"`**. Wraps `PyHost` and
   (read-mostly) `SceneGraph`. Co-resident in the one interpreter with `pycvc`/`pycvc_gl`.

4. **`volrover3::PyConsoleDock`** (Qt dock widget) — the control surface (verlihub's `cconsole` ported 1:1 to
   Qt): load/reload/unload/list scripts, a REPL line, log-level combo, output pane (redirected
   `sys.stdout`/`sys.stderr` + surfaced `PyErr`). Mirrors how `StateDashboardWidget` is created
   (`MainWindow.cpp:767`).

## 5. App delivery — the crux

pycvc (#136) crosses `cvc::app` **only** as `std::shared_ptr<cvc::app>` (`%shared_ptr(cvc::app)`,
`pycvc.i:99`), with **no module-global current app** — every op takes it explicitly (`state_set(app,…)`,
`volume(app)`, `sdf(app,…)`, `observer.watch(app)`). With the no-singleton refactor (§3) volrover3 **owns** a
`std::shared_ptr<cvc::app>`, so `PyHost::app()` simply returns it — real ownership, real control block, no alias
hack, no dangling risk (the interpreter is torn down in `~MainWindow` before the app shared_ptr drops).

**Delivery mechanism (no SWIG in volrover3, no `pycvc.i`).** volrover3 does **not** build a SWIG module and
does **not** carry `pycvc.i`. Instead `EmbeddedInterpreter::bind_host()` wraps the live app in a
`PyCapsule` named `"cvc.app"` (holding a heap `std::shared_ptr<cvc::app>` copy — shared ownership; the capsule
destructor frees the copy) and injects it into a **pure-Python `vrhost` shim** as `vrhost._app_capsule`.
`vrhost.app()` feeds that capsule to **`pycvc.app_from_capsule(cap)`**, which wraps it into **pycvc's OWN** app
proxy (the same `%shared_ptr(cvc::app)` typemap as `make_app`). The handle is therefore type-compatible with
pycvc **by construction** — no cross-module SWIG type sharing, and no SWIG-runtime-version coupling (the fragility
that shipping `pycvc.i` + `%import` used to impose). The script:

```python
import vrhost, pycvc
app = vrhost.app()                                   # THE live volrover3 app, as a pycvc app handle
pycvc.state_set(app, "volrover3.camera.fov", "42")   # writes the RUNNING app's state tree
print(pycvc.state_get(app, "volrover3.camera.fov"))  # -> 42; the FOV widget updates on the UI thread
```

Because `app` **is** the process app, `cvc::state::instance(*app)` inside pycvc *is* the tree the
`StateTreeWidget`/`StateDashboardWidget` and `SceneNode`s already watch — the round-trip is observable
end-to-end. A script that instead called `pycvc.make_app()` would get a **fresh, disconnected** tree the UI never
sees — exactly the bug #136 + this host module exist to prevent.

The app is injected at boot as `vrhost._app_capsule` (the `PyCapsule("cvc.app")` above), and the QMainWindow
address as `vrhost._main_window_ptr` (pushed by `EmbeddedInterpreter::set_main_window_ptr` once the window
exists) — the **injected-handle** pattern, no singleton, no reflection, no SWIG.

## 6. Single vs multi interpreter policy

**Default = SINGLE**; MULTI exposed as a `Mode` enum but **off by default and forbidden for any script importing
pycvc/pycvc_gl/vtk/numpy.** Rationale: `Py_NewInterpreter()` sub-interpreters share single-phase-init C
extensions' C-global state; a second sub-interpreter importing such a module **aliases and corrupts** those
globals. pycvc/pycvc_gl wrap libcvc/VTK (exactly this character) and are not sub-interpreter-safe; VTK also
self-initializes in `main.cpp` (`VTK_MODULE_INIT` + `vtk_module_autoinit`, `CMakeLists.txt:200`) so a second
interpreter re-initializing VTK-python risks double init. Per-interpreter GIL is only 3.12+ (PEP 684) and only
for modules that opt in via multi-phase init — which these don't.

**Single mechanics:** one `Py_InitializeEx(0)` for the process; per-script isolation at the Python level (a
fresh module/globals dict per script in the shared `sys.modules`), **not** `Py_NewInterpreter`. All C extensions
imported **once** and shared. **GIL:** `PyImport_AppendInittab` for `vrhost` (+ any baked pycvc) **before**
`Py_InitializeEx`; `PyEval_SaveThread()` to park; every C++→Python entry brackets `PyGILState_Ensure/Release`;
blocking cvc ops called from Python release the GIL (SWIG `-threads` / `Py_BEGIN_ALLOW_THREADS`) so
`app().startThread` workers aren't starved.

## 7. The Qt bridge — manipulating & replacing the UI from Python

**Yes — confirmed possible and tractable.** This is precisely the architecture of Qt-for-Python's official
**"Scriptable Application"** example: a C++ Qt app that owns the `QApplication` + `MainWindow`, embeds CPython,
and hands the live C++ `MainWindow*` to Python via shiboken (`pythonutils.cpp` binds it into the interpreter's
globals as `mainWindow`, after which Python calls dispatch into the real C++ object). Embedded-Python-in-a-
C++-QApp is a first-class, upstream-supported pattern — not the usual Python-drives-Qt.

**Mechanism.** Qt objects are *not* wrapped with SWIG — they use **PySide6 + Shiboken6** (Qt's official Python
bindings). The bridge between the two binding worlds is a **raw pointer exchanged as an integer** — the same
pattern as the already-shipped pycvc↔VTK `vtkPythonUtil` bridge (SWIG and Shiboken type tables are mutually
opaque, so the *only* safe handoff is an `int` address):

- `PyHost` exposes `quintptr main_window_ptr()` (the live `QMainWindow*` as an int).
- A Python helper `vrhost.main_window()` internally calls `shiboken6.wrapInstance(ptr, QtWidgets.QMainWindow)`
  → a **live PySide6 wrapper around the existing C++ QMainWindow**, in the one C++ `QApplication` event loop.
- Reverse direction: `shiboken6.getCppPointer(pyWidget)` hands a Python-created widget's `QWidget*` back to C++.

Once a script holds the `QMainWindow`, it can do **anything Qt exposes** — on the GUI thread:
- `setCentralWidget(w)` to swap the central widget; `addDockWidget/removeDockWidget`; rebuild menus/toolbars.
- Build whole `QWidget`/`QLayout` trees in PySide6 (incl. Python `QWidget` subclasses with signals/slots) that
  live natively in the C++ app.
- **Outright replace the UI:** tear down the C++ docks/central and install a fully Python-authored widget tree —
  leaving the C++ side as `QApplication` + interpreter + the backend objects (`vrhost`/`pycvc`). The one thing
  that stays C++ is the bootstrap (`QApplication`, interpreter init) and the VTK render widget's GL context
  (which can still be *re-parented* by Python).

**Hard constraints (confirmed against upstream):**
- **One Qt, one ABI.** PySide6 links whatever Qt it is built against; it **must** be the **same cvcpkg Qt6**
  volrover3 links, or QObject vtables mismatch → UB/crash. The Scriptable-Application docs state this outright
  ("use the same Qt version … to ensure binary compatibility"). This is the hermeticity linchpin.
- **One `QApplication`, one loop.** The embedded Python side must adopt the existing app
  (`QtWidgets.QApplication.instance()`) and **must not** create its own or call `app.exec()` — the C++ host owns
  the loop. All wrapped objects + Python-created widgets are driven by that one loop automatically.
- **GUI-thread only + GIL.** All Qt manipulation on the main (interpreter-calling) thread; background Python
  marshals via `QMetaObject::invokeMethod(QueuedConnection)` / queued signals; every Qt→Python slot fired from
  the loop must hold the GIL (composes with §6).
- **Ownership.** `shiboken6.wrapInstance(addr, T)` adopts an existing C++ object **without** taking ownership
  (Python GC won't delete the C++ `QMainWindow`) — but if C++ deletes it, the wrapper goes invalid (guard with
  `shiboken6.isValid`). Python-*created* widgets start Python-owned; **parenting them into the C++ tree**
  (`setCentralWidget`, `layout.addWidget`) transfers ownership to the Qt parent — so always parent injected
  widgets immediately (or hold a Python ref until you do). `setCentralWidget(new)` deletes the old central
  widget (Qt semantics), invalidating any wrapper of it.
- **Reverse handoff:** `shiboken6.getCppPointer(pyWidget)[0]` (returns a tuple — take element 0).

**Gating dependency — DONE.** Hermetic **`shiboken6` + `pyside6` 6.8.2 cvcpkg recipes** (libcvc-deps #373),
built via direct CMake against the cvcpkg `qt6` 6.8.2 (non-standalone → single `libQt6Core` in-process) and a
hermetic **`llvm18`** recipe for libclang (no system `libclang-18-dev`). One `std::uintptr_t`→64-bit SWIG fix
in `vrhost.i` was load-bearing (SWIG defaulted it to 32-bit, truncating the `QWidget*`). Optionally a
`volrover3.qt` Shiboken module wrapping volrover3's *own* widgets (`VTKRenderWidget`, dialogs) for typed
handles — additive; the generic `QMainWindow` handle already unlocks full control via Qt's own API.

### 7.1 Canonical example — add a menu item that pops a message box

Full runnable script: [`scripts/examples/menu_messagebox.py`](../scripts/examples/menu_messagebox.py). From the
Python Console dock's REPL: `exec(open("scripts/examples/menu_messagebox.py").read())`. The essence:

```python
import vrhost
from PySide6 import QtWidgets

window = vrhost.main_window()                      # the live C++ QMainWindow, adopted by shiboken
action = window.menuBar().addMenu("&Demo").addAction("Say &Hello")
action.triggered.connect(lambda: QtWidgets.QMessageBox.information(
    window, "volrover3", "Hello from embedded Python!"))
vrhost._demo_action = action                       # keep the slot's Python objects alive
```

`window.menuBar().addMenu(...)` mutates the **real** menu bar; the `QMessageBox` is a genuine modal parented to
the app window. Verified headless by `tests/VrHostQtBridgeTest.cpp`.

## 8. SWIG plan

Two binding families in **one** interpreter, sharing SWIG's runtime type table:
- **pycvc / pycvc_gl** — reused as-is from #136 (already a cvcpkg package). `%module(directors="1") pycvc` with
  `%shared_ptr(cvc::app)`; `pycvc_gl` already `%import`s `pycvc.i`.
- **`bindings/vrhost.i`** (new) — `%module vrhost`, **`%import "pycvc.i"`** (mandatory: otherwise `vrhost` mints
  a *distinct* `shared_ptr<cvc::app>` type and `pycvc.state_set(app)` rejects the handle). `%{ #include
  "PyHost.h" %}` + `%include "PyHost.h"`; a thin `%inline PyHost* _current_host();` plus a `.py`-level
  `host = _current_host()`. Add `directors="1"` + a director base if Python→C++ event hooks (`OnSceneChanged`,
  `OnTimer`) are needed. Unlike verlihub (hand-written `PyMethodDef`, `Py_InitModule`, no SWIG), we use
  SWIG-generated wrappers that yield **real proxy objects** — the reason to wrap real volrover3 objects.

## 9. Build plan

`volrover3/CMakeLists.txt` (no SWIG/Python today):
- `find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Embed)` + `find_package(SWIG REQUIRED)` +
  `include(UseSWIG)`.
- Link `Python3::Python` into `volrover3_lib` (so the subsystem is unit-testable); the executable inherits it.
- `target_link_options(volrover3 PRIVATE -Wl,--export-dynamic)` so libpython symbols reach imported
  C-extensions (numpy/vtk) — the executable equivalent of `RTLD_GLOBAL`.
- Generate `vrhost` and compile `vrhostPYTHON_wrap.cxx` **into** `volrover3_lib`, registered via
  `PyImport_AppendInittab("vrhost", PyInit_vrhost)` (no `.so` to locate on disk). SWIG needs pycvc's include dir
  on `-I`.
- Ship pycvc/pycvc_gl as cvcpkg runtime deps (already published); import from bundled site-packages. Keep the
  **same libcvc** (`cvc::cvc`) so `cvc::app` ABI matches. Ensure VTK-python matches the app's VTK.
- Ship a `scripts/` dir the console enumerates.
- **cvcpkg recipe:** volrover3 gains build deps `swig` + `python3-dev`; runtime deps `python3` + `pycvc` +
  `pycvc_gl` (+ transitive numpy/vtk-python); later `pyside6` + `shiboken6` (§7). Validate with `cvcpkg validate`.

## 10. Phased plan

- **Phase 0a — cvcGL extraction.** Delete volrover3's 9 forked node classes; repoint includes to `<cvc/gl/…>`
  and link `cvc::cvcGL`; reconcile the small `SceneGraph`/`SceneNode` deltas (port any volrover3-only tweaks up
  into cvcGL first); make `VTKRenderWidget` a thin host over cvcGL's `SceneGraph::setRenderer`. Verify the app
  renders identically. This shrinks volrover3 to the Qt+state shell and makes the 3D engine one shared,
  Python-wrapped (`pycvc_gl`) library.
- **Phase 0b — no-singleton app ownership.** Own one `shared_ptr<cvc::app>` in `main()`/`MainWindow`; thread it;
  delete `volrover3::app()` (155 sites) + de-singleton `AppState` (45 sites); construct the (now cvcGL)
  `SceneGraph` with the injected-app ctor.
- **Phase 0c — build scaffolding.** Add Python3/SWIG to CMake, `-Wl,--export-dynamic`; prove volrover3 still
  builds/links with embedded libpython.
- **Phase 1 — smallest end-to-end slice.** `EmbeddedInterpreter` boots once + parks the GIL; `PyHost::app()`
  returns the owned handle; `vrhost.i` compiled in + `vrhost.host` injected; pycvc importable; `run_string()`.
  **Acceptance gtest (in `volrover3_lib`):** run `import vrhost, pycvc; pycvc.state_set(vrhost.host.app(),
  'volrover3.camera.fov','42')` and assert C++-side that `cvc::state::instance(app)('volrover3')('camera.fov')
  == '42'` **and** the camera-changed handler fired — proving the script drove the *running* app.
- **Phase 2 — GIL + threading correctness.** Worker→Python and Python→UI through the existing channels; a Python
  `pycvc.state_observer` `watch(app)` whose `on_changed(path)` fires on cvc worker writes, marshaled to the UI
  via `SceneGraph::postEvent` / `QMetaObject::invokeMethod`. Add SWIG `-threads`.
- **Phase 3 — script lifecycle + isolation.** load/exec/unload `.py` from `scripts/` into per-script namespaces
  in the one interpreter; name-based hook discovery + exception isolation (ported from verlihub).
- **Phase 4 — PyConsoleDock Qt UI** (REPL + load/reload/unload/list + output pane).
- **Phase 5 — host-state control surface** (typed PyHost helpers; director hooks for host events; API docs).
- **Phase 6 — Qt bridge (PySide6/Shiboken).** `main_window_ptr()` + `vrhost.main_window()` wrapInstance helper;
  demo scripts that add a dock and replace the central widget from Python. **Gated on** the hermetic
  `pyside6`/`shiboken6` cvcpkg recipes.
- **Phase 7 — packaging + CI** (cvcpkg deps; Linux/macOS/Windows; ship `scripts/` + site-packages).

## 11. Open decisions (for review) & risks

**Decisions:**
- **AppState de-singletoning scope** — full removal now (Phase 0), or keep `AppState` as an owned object but
  retain a thin `instance()` shim temporarily to bound the diff? (Rule says remove; flagging the size.)
- **pycvc distribution into the interpreter** — bundle prebuilt `.so` on a private `PYTHONPATH`, or bake the
  wrapper into `volrover3_lib` via inittab (like `vrhost`)?
- **State-exposure boundary** — whole tree, just the `volrover3` subtree, or also `volrover3.graphics.root.*`
  (SceneNode-owned, different lifetime)?
- **Script capability boundary** — scripts currently get *full* host access (verlihub-style). Do we want an
  opt-in sandbox/capability layer, given §7 lets a script replace the entire UI?
- **When to pull PySide6/Shiboken** — the §7 Qt bridge is high-value but the recipes are the heavy lift; land
  Phases 0–5 (pycvc/app/state/console) first, or parallelize the recipe work?

**Risks:** static-destruction-order UB (mitigated by `~MainWindow` teardown before exit-time statics); sub-
interpreter corruption if MULTI ever enabled for C-ext scripts (mitigated by the hard policy); SWIG cross-module
type mismatch if `vrhost.i` fails to `%import pycvc.i` or is built against a different libcvc/SWIG (enforce in
CI); GIL deadlocks from a missing bracket (verlihub's load-bearing discipline); libpython symbol visibility
without `--export-dynamic`; VTK double-init between the app and pycvc_gl; the `readOnly(true)` world-bounds write
no-op gotcha; PySide6 ABI mismatch vs cvcpkg Qt6; cvcpkg closure completeness for the new dep set (bit libcvc
before — imagemagick/libxml2).
```

*Grounding: verlihub `plugins/python` (cpythoninterpreter/cpipython/wrapper/cconsole); volrover3
`volrover3_app.cpp`, `MainWindow.cpp`, `AppState.*`, `SceneNode.*`, `CMakeLists.txt`; pycvc `bindings/pycvc`
(#136); CVCPKG-ROADMAP Phase 19/20 + `static-single-binary-python.md`; the Python/Qt binding + hermeticity
architecture doc.*

---

## 12. Scheduler, interpreter modes, settings & console (design — grounded in verlihub's `dispatcher.py`)

Reference material vendored under `docs/reference/verlihub-scripts/` (verlihub `plugins/python/scripts/`:
`dispatcher.py` + `README.md`). verlihub's "scheduler" is a **host tick → `CallAll` fan-out → per-script
`OnTimer`** loop (`casyncsocketserver.cpp` 1 Hz clock → `cserverdc.cpp:1954` `CallAll` → `cpipython.cpp:1087`
`OnTimer` → `:465-519` in-order iterate `mPython`), and `dispatcher.py` is the Python-level **hook dispatcher**
that rides it — a script *registry* (`register_script`/`unregister_script`/`enable_script`/`disable_script`/
`list_scripts`) that fans a hook out to all registered jobs, solving single-interpreter hook collisions. That
registry maps 1:1 onto the jobs tab.

### 12.1 JobScheduler — the tick, ported to Qt
A thin UI-thread `volrover3::JobScheduler` driven by a `QTimer(tickMs)` (default `tick_ms: 100`), the direct
analog of verlihub's 1 Hz socket loop. `onTick()` = cooperative fan-out over an ordered `std::vector<Job>`:
skip `!online` / no-`step` jobs, call each `step(dt)` in order under one `PyGILState_Ensure/Release`, per-job
`try/catch` so one bad job neither kills the loop nor stalls the next tick (verlihub invariants:
`wrapper.cpp:1475-1483` cache the hook bitmap once, never re-`hasattr`; `cpythoninterpreter.cpp` guard). NOT
libcvc's `state_exec` engine — that steps a stackless DSL evaluator over `value_t` AST and can't host a CPython
program; its `register_fn` runs a Python callable only as a synchronous DSL leaf. We **borrow state_exec's
vocabulary** (`pid→job_id`, the `process_status` enum ready/running/paused/waiting/terminated/killed) so the
jobs tab renders uniformly, and optionally wrap `pycvc::Exec` as a `JobKind::Dsl` entry (Phase 7) that inherits
real pid/kill.

`JobInfo { int id; std::string name; JobStatus status; JobKind kind; InterpreterMode mode; double elapsed;
uint64_t steps; std::string lastError; }`. **Stop** = drop the job + clear its module namespace (single) /
`Py_EndInterpreter` (multi). **Interrupt** = `PyThreadState_SetAsyncExc(tid, KeyboardInterrupt)` — fires only at
bytecode boundaries, so a tight C-extension loop (VTK/pycvc) is **not** mid-call preemptible; the UI must say so
honestly. (Hard-kill of a runaway C loop would need a sacrificial worker thread — deferred.)

### 12.2 Single vs multi interpreter — a restart-applied settings flag (feasible; verlihub ships both)
verlihub proves both modes are production-ready (10/10 sub-interpreter, 11/11 single) — it selects at *compile*
time; volrover3 makes it a **`~/.volrover/settings.yaml` flag applied on restart** (never hot-switched:
`Py_Initialize`/`Py_Finalize` + the `vrhost` inittab registration + GIL park are once-per-process). The
once-per-process boot is identical in both modes; only the **job launcher** branches on the mode.
- **SINGLE** (default, the real product): all jobs share the one interpreter + `sys.modules`; cooperatively
  ticked; **full package compat** — pycvc/pycvc_gl/VTK/numpy work. The *only* mode with the live-app headline
  feature.
- **MULTI** (weaker sandbox): each job is a `Py_NewInterpreter` sub-interpreter, still ticked from the same
  loop (isolation ⟂ scheduling — verlihub already does exactly this). Such a job **cannot import
  pycvc/pycvc_gl/vrhost/vtk/vtkmodules/numpy** (single-phase-init C-globals corrupt across sub-interpreters —
  §6). So MULTI = pure-Python compute sandboxes with **no host/app/scene handle**.
- **Enforce the gate in code, not docs** (an errant `import pycvc` in a sub-interpreter is memory corruption,
  not a catchable error): a `PySys_AddAuditHook` + `sys.meta_path` denylist raising `ImportError` for
  `{pycvc, pycvc_gl, vrhost, vtk, vtkmodules, numpy}` in every MULTI job; plus a UI warning "app scripting
  (pycvc/vrhost) is disabled in multi-interpreter mode." `EmbeddedInterpreter` gains a `Config { InterpreterMode
  mode; std::string python_home; bool gate_multi_imports; }`, fed from `Settings` at the `MainWindow` boot site.

### 12.3 Settings — state-tree backed, persisted to `~/.volrover/settings.yaml` + `data.db`
**Settings are backed by the `cvc::state` tree, not a bare file.** Each volrover3 instance owns a dedicated
**section of the global state root** — the same `"volrover3"` prefix the SceneGraph + AppState already run under
(`SceneGraph(*m_app,"volrover3")`, `AppState(*m_app,"volrover3")`) — so **one instance section holds it all**:
`volrover3.settings.*` (interpreter mode, python home, tick, …), `volrover3.camera.*` (AppState),
`volrover3.graphics.*` (the scene graph), etc. Multiple instances get distinct sections
(`volrover3`, `volrover3.1`, … — the prefix is the instance id), so their settings + scene state never collide
in `cvc::state::instance(app)`. `volrover3::Settings` is a **typed facade over
`cvc::state::instance(app)("volrover3")("settings")`** (reads/writes state keys there via `state.value(...)`),
NOT a private field bag — so a Python script (`pycvc.state_set(app,"volrover3.settings.scheduler.tick_ms",…)`)
and the C++ settings UI see the same live values, and the reactive `childChanged` signals already drive updates.

**Persistence** = load/save the instance's state section to disk. `~/.volrover/settings.yaml` (yaml-cpp) mirrors
the `volrover3.settings.*` subtree (human-editable metadata: mode, python home, tick, scripts dir, history);
`~/.volrover/data.db` (Qt6Sql/WAL) holds bulkier/arbitrary persisted state (console history, job audit, KV, and
optionally scene snapshots). On startup `Settings::load()` reads yaml → `state_set` into
`volrover3.settings.*`; `save()` emits that subtree → yaml. A plain object (not a singleton),
**constructed before `EmbeddedInterpreter`** so `pythonHome()` feeds `PyConfig.home` and `mode()` selects the
launcher. (The prior "field bag" framing is superseded: the fields live in `cvc::state`.)
- **`settings.yaml`** (human-editable metadata): `interpreter: {mode: single, python_home: "",
  gate_multi_imports: true, scripts_dir: ~/.volrover/scripts}` · `console: {history_size: 500}` ·
  `scheduler: {tick_ms: 100}`. **YAML lib:** deps-vr3 has only libyaml (C); no yaml-cpp. → add a small
  **yaml-cpp cvcpkg recipe** (house style, `cvcpkg validate`) OR fall back to wrapping libyaml C for the fixed
  shallow schema. (Decision pending — see below.)
- **`data.db`** (arbitrary persisted state): **sqlite via Qt6Sql** (`QSQLITE`, already linked) with
  `PRAGMA journal_mode=WAL`. Tables `kv(key,value,updated_at)`, `console_history(id,ts,source)`,
  `job_runs(id,name,mode,started,ended,status,error)`; a raw-sqlite3-C second WAL connection is reserved for any
  future core-side `cvc::state` persistence to the same file.

### 12.4 PyConsoleDock — REPL + Jobs tabs
A `QDockWidget` cloned structurally from `StateDashboardWidget` (which already has an exec-console tab + a
process table with Pause/Resume/Kill + a refresh `QTimer`). Installed in `MainWindow::createDockWidgets()`.
- **REPL tab**: read-only mono output + input line with data.db-backed Up/Down history. Enter runs on the UI
  thread via a **new** `EmbeddedInterpreter::run_string_capture(src, out, err)` (swap `sys.stdout/stderr` to
  `io.StringIO` under the GIL, `PyRun_String(..., Py_single_input, ...)` for REPL echo, restore in a
  finally/RAII) — `run_string` is left unchanged.
- **Jobs tab**: a `QTableWidget` (Id/Name/Status/Mode/Steps/Elapsed) polled 1 Hz from `JobScheduler::listJobs()`;
  Interrupt → `interrupt(id)`, Stop → `kill(id)`. Python jobs + (optional) DSL jobs unify in one table; the live
  `state_exec` DSL scheduler keeps its own StateDashboard tab.

### 12.5 Phased plan
1. **Settings infra** (Settings class, `~/.volrover`, yaml + Qt6Sql data.db) — unblocked. · 2. **Mode plumbing**
(`Config` into `EmbeddedInterpreter`, restart-only `setMode`). · 3. **`run_string_capture`** (stdout/stderr
capture). · 4. **JobScheduler** (single mode: QTimer tick + registry). · 5. **PyConsoleDock** (REPL + Jobs). ·
6. **Multi mode + import gate** (`Py_NewInterpreter`/`Py_EndInterpreter` + audit-hook denylist + UI warning). ·
7. *(optional)* **DSL job kind** (wrap `pycvc::Exec`). Runs alongside the `vrhost` registration (Phase 1 core,
runtime-gated on imagemagick for `import pycvc`).

### 12.6 Decisions (locked)
- **YAML: add a `yaml-cpp` cvcpkg recipe** (libcvc-deps/recipes/, `cvcpkg validate`, built into deps-vr3) — the
  `Settings` class round-trips settings.yaml via `YAML::LoadFile`/`YAML::Emitter`.
- **Interrupt model: cooperative tick + a hard-kill worker path.** Normal jobs run cooperatively on the UI-thread
  tick (interrupt = `PyThreadState_SetAsyncExc`). A job may also be launched on a **sacrificial worker thread**
  so a runaway C-extension loop can be **force-terminated** (hard stop). Worker-thread jobs own their GIL via
  `PyGILState_Ensure` per step and marshal any scene/UI effect back through `SceneGraph::postEvent` /
  `QMetaObject::invokeMethod(QueuedConnection)` (never touch VTK/widgets off the UI thread). `JobKind`/launch
  carries a `runOnWorker` flag; the Jobs tab's Stop offers cooperative-stop for tick jobs and hard-kill (thread
  cancel + `Py_EndInterpreter` for its sub-interpreter, or interpreter-scoped teardown) for worker jobs.
- Tick **100 ms**; jobs submittable from **both** the REPL and `scripts_dir` (alphabetical load like verlihub);
  Python + DSL jobs **unified** in one table; `cvc::state`→data.db **out of scope** for now.
- **Risk to honor (hard-kill):** force-terminating a thread mid-C-call can leak/corrupt that job's interpreter
  state — hard-kill tears down the whole *job* (its sub-interpreter) and quarantines it, never the process; the
  UI marks a hard-killed job "terminated (unclean)". Cross-thread GIL hand-off + scene marshaling is the load-
  bearing discipline (§6 rules apply per worker).

**IMPLEMENTATION FINDING (the hard limit of "hard-kill"):** CPython gives no safe way to force-terminate a
thread that is **hung inside a C-extension call while holding the GIL**. `PyThreadState_SetAsyncExc` only fires
at Python **bytecode boundaries** — it breaks a pure-Python `while True: pass` but *cannot* interrupt a tight
C loop (a long VTK/pycvc call). And you must **never** `pthread_cancel`/detach-then-free a GIL-holding thread:
in single-interpreter mode the abandoned thread keeps the one GIL → the whole interpreter deadlocks; freeing the
job's `PyObject`s under it is a use-after-free. So the worker path is engineered as: (a) each worker steps on its
own thread and **releases the GIL between steps**, so a *yielding* slow job never freezes the app and can be
stopped cleanly (stop flag + join) or interrupted (`SetAsyncExc`, effective for Python loops); (b) a job that is
**genuinely hung in a C call** cannot be force-killed — `kill()` raises `SetAsyncExc`, then on a join timeout
**detaches the thread and moves the (still-referenced) Job to a zombie list** (never freed, never joined at
dtor) marking it `Killed (unclean)`, and stops managing it. True force-termination of such a job is only the
process boundary. This is honest: the worker path buys off-UI-thread execution + clean/interruptible stop for
the common case, and graceful quarantine (not a crash) for the pathological one. (Full force-kill would need
per-job **multi-mode** sub-interpreters *and* a cooperative job — still not a hung C-loop.)

## 13. As-built: the `vrhost` module (Phase 1 capstone) — deltas from §5/§8/§9

The design's app-delivery crux (§5) is implemented and verified end-to-end
(`tests/VrHostBindingTest.cpp`: a script calls `pycvc.state_set(vrhost.host.app(), …)`
and C++ reads the value back out of the **same** `cvc::state` tree). Two deviations
from the written plan, plus one hard operational constraint:

- **Bundled `_vrhost.so`, not baked-into-`volrover3_lib` + inittab.** §9 proposed
  compiling `vrhostPYTHON_wrap.cxx` into `volrover3_lib` and
  `PyImport_AppendInittab("vrhost", …)`. As built, `vrhost` is a normal
  `swig_add_library` module (`_vrhost.so` + `vrhost.py`) shipped in
  `share/volrover3/pymod`; `EmbeddedInterpreter` prepends that dir to `sys.path`
  (`Config::module_path`, else `$VOLROVER3_PYMODULE_PATH`) and `import`s it at boot.
  This reuses pycvc's proven UseSWIG path and keeps `PyHost.cpp` self-contained
  (it compiles into `_vrhost.so`, linking only `cvc::cvc`/`cvc::cvcGL`).
- **Host handed over by `PyCapsule`, not a direct `SWIG_NewPointerObj` in the host TU.**
  Because `_vrhost.so` is a separate object from `volrover3_lib`, the interpreter
  passes the live `PyHost*` across as a name-matched capsule
  (`PyCapsule_New(host, "volrover3.PyHost", …)` → `vrhost._bind_host_capsule`),
  then refreshes the module-level `vrhost.host`. No SWIG type is needed from the
  host TU. `EmbeddedInterpreter::bind_host()` does this once, GIL-held, right after
  init (and in the "already-initialized" reuse path via `PyGILState_Ensure`).
- **CRITICAL — SWIG runtime-version must match pycvc's.** Cross-module type sharing
  (the whole point of `%import "pycvc.i"`) works only if both modules register into
  the **same** SWIG runtime type table, whose capsule name is versioned:
  SWIG 4.2.x → `swig_runtime_data4`, SWIG ≥4.3 → `swig_runtime_data5`. The published
  **pycvc was built with SWIG 4.2.0**, but the cvcpkg `swig` recipe is **4.4.1**
  (`…data5`). Building `vrhost` with 4.4.1 makes `pycvc.state_set(vrhost.host.app())`
  fail with `TypeError: argument 1 of type 'std::shared_ptr< cvc::app > const &'`
  (+ a `"no destructor found"` leak) even though the mangled type names are identical
  — the tables simply never link. **vrhost must be generated with the same SWIG
  runtime version pycvc was.** Locally this is the system `/usr/bin/swig` 4.2.0
  (`-DSWIG_EXECUTABLE=/usr/bin/swig -DSWIG_DIR=$(swig -swiglib)`). The durable fix is
  to align the ecosystem: rebuild/republish pycvc (and pycvc_gl) with the cvcpkg
  4.4.1 swig, **or** pin the cvcpkg swig recipe to 4.2.x — a coordination decision,
  since the SWIG runtime bump is not ABI-compatible across the pycvc family.

Two `%import`-specific SWIG mechanics also had to be handled in `vrhost.i` (all
documented inline): `%import` does **not** re-emit the base module's `%{ … %}`
header block, its library `%include`s, or its `%shared_ptr` **typemaps** — only the
type *registrations*. So `vrhost.i` re-`#include`s `<cvc/core/exception.h>`, declares
its own self-contained `%exception` (via the always-emitted `SWIG_exception_fail`),
re-declares `%include <std_shared_ptr.i>` + `%shared_ptr(cvc::app)`, and forces
`import pycvc` **before** the `_vrhost` C-extension via `%pythonbegin` (SWIG links
type tables in import order — pycvc must register first).
