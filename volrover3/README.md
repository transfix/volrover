# VolRover3 - Volume Rover Version 3

A prototype visualization application built on libcvc for rendering volumetric data, surface meshes, and volumetric meshes.

## Features

- **Volume Rendering**: 3D texture-based volume rendering with GPU acceleration via VTK
- **Surface Mesh Visualization**: Triangle mesh rendering with normals and colors
- **Volumetric Mesh Support**: Tetrahedral and hexahedral mesh visualization
- **Transfer Functions**: Interactive color and opacity mapping for volume data
- **Quake-Style Camera**: First-person camera controls for intuitive navigation
- **Scene Elements**: Toggleable grid and coordinate axis for reference
- **File I/O Integration**: Support for CVC geometry and volume formats
- **Embedded Python**: script the *live* scene from a built-in CPython interpreter
  (Python Console / Jobs tab / `--run-job`) via `pycvc` / `pycvc_gl` / `vrhost`
  (see [Embedded Python & Scripting](#embedded-python--scripting))

## Building

VolRover3 requires:
- An installed [libcvc](https://github.com/transfix/libcvc) SDK
  (consumed externally via `find_package(cvc CONFIG)`; links `cvc::cvc`)
- Qt6 (Core, Widgets, OpenGL, OpenGLWidgets) — Qt5 fallback supported
- VTK (Visualization Toolkit) 9.0+

The libcvc SDK's own dependency closure (Boost, CGAL, HDF5, FFTW3,
ImageMagick, libiimod) must also be resolvable at configure time —
point `CMAKE_PREFIX_PATH` at a cvcpkg prefix or an extracted
libcvc-deps bundle alongside the SDK itself.

### Build Steps

From the volrover repository root (VolRover3 is disabled by default so
the legacy VolumeRover 2.0 build is unaffected):

```bash
cmake -B build -DBUILD_VOLROVER3=ON \
  -DCMAKE_PREFIX_PATH="/path/to/libcvc-sdk;/path/to/deps-prefix"
cmake --build build --target volrover3
```

Unit tests (Google Test) build by default with the app; disable with
`-DVOLROVER3_BUILD_TESTS=OFF` or run them with:

```bash
ctest --test-dir build/volrover3
```

## Usage

### Launch

```bash
./bin/volrover3
```

### Controls

**Camera Movement (Quake-Style)**:
- `W` - Move forward
- `S` - Move backward
- `A` - Strafe left
- `D` - Strafe right
- `E` or `Space` - Move up
- `Q` or `Ctrl` - Move down
- `Mouse drag` (left button) - Look around
- `Mouse wheel` - Zoom in/out

### Menu Options

**File Menu**:
- `Open Geometry...` - Load surface meshes (.off, .raw, .obj, etc.)
- `Open Volume...` - Load volume data (.rawiv, .mrc, .ccp4)

**View Menu**:
- `Show Grid` - Toggle ground grid display
- `Show Axis` - Toggle coordinate axis display

## Embedded Python & Scripting

VolRover3 embeds a first-class CPython interpreter that drives the **live** scene
through the `pycvc` / `pycvc_gl` bindings. The `vrhost` module hands a script the
running app and scene, so a script mutates what you see in the window — this is
how the GRL-SNAM lab and the Austin drive demo run inside the viewer.

### Install & activate (recommended)

Install into a cvcpkg prefix and use **that prefix's own activation script** — no
manual environment variables, because cvcpkg's `activate` already exports `PATH`,
`LD_LIBRARY_PATH`, `CMAKE_PREFIX_PATH`, and `PKG_CONFIG_PATH` for the prefix:

```bash
cmake --install build --prefix /path/to/cvcpkg-prefix --component volrover3
source /path/to/cvcpkg-prefix/bin/activate    # cvcpkg's script; sets up the env
volrover3                                      # just works
```

The installed binary self-locates its CPython home (`<prefix>/lib/pythonX.Y`) and
the `vrhost` module (`<prefix>/share/volrover3/pymod`) from its own location, and
its `RUNPATH` (`$ORIGIN/../lib`) finds the bundle's shared libs — so no
`VOLROVER3_*` variables are needed. Example scripts install to
`<prefix>/share/volrover3/examples/`.

### Running a script

Three ways, all against the same embedded interpreter:

- **Interactively** — open the **Python Console** dock, or its **Jobs** tab →
  *Load Script…* → *Run as Job*. A script that defines `step(dt)` is ticked
  cooperatively by the scheduler (it appears in the Jobs tab; pause/stop there).
- **At startup, as a job** — `volrover3 --run-job my_demo.py` submits the file as
  a scheduler job once the window is up.
- **At startup, once** — `volrover3 --exec-script setup.py` runs the file once (no
  scheduler). Add `--screenshot out.png --exit-after` for a headless capture.

### The live-scene API

```python
import vrhost, pycvc, pycvc_gl
app   = vrhost.app()        # the running cvc::app
scene = vrhost.scene()      # the live pycvc_gl.SceneGraph

g = pycvc.geometry(app)
g.add_vertices([0, 0, 0, 10, 0, 0, 0, 10, 0]); g.add_triangle(0, 1, 2)
scene.addGraphics("tri", g)                     # appears in the window immediately
scene.getGraphics("tri").setPosition(1, 2, 3)   # move it in place (no rebuild)
scene.getGraphics("tri").setTransform(m16)      # or a full 4x4 (row-major list)

def step(dt):                                   # optional: animate, ticked per frame
    scene.getGraphics("tri").setPosition(...)
    scene.processEvents()
```

Drive the camera and viewer settings through the state tree:

```python
# direct camera control (also: view_direction.{x,y,z}, up_vector.{x,y,z}, fov)
pycvc.state_set(app, "volrover3.camera.position.x", "10.0")
# tunable render/refresh rate, 1–240 fps (see State Management)
pycvc.state_set(app, "volrover3.viewer.max_fps", "60")
```

Hide reference gizmos for a clean demo:
`scene.setAxisVisible(False)` / `scene.setGridVisible(False)`.

## Supported File Formats

**Geometry**:
- `.off` - Object File Format
- `.raw`, `.rawn`, `.rawc`, `.rawnc` - CVC raw formats
- `.obj` - Wavefront OBJ (experimental via SDF)

**Volume**:
- `.rawiv` - RAWIV format
- `.mrc` - MRC/CCP4 format
- Other formats supported by libcvc

## Architecture

### Components

- **MainWindow**: Qt6 main application window with menus and docking
- **VTKRenderWidget**: VTK/OpenGL rendering widget with event handling
- **SceneGraph**: Scene management and traversal
- **SceneNode**: Base class for renderable objects
  - **GeometryNode**: Surface mesh rendering
  - **VolumeNode**: Volume rendering with transfer functions
  - **GridNode**: Reference grid
  - **AxisNode**: Coordinate axis
- **CameraController**: Quake-style first-person camera
- **TransferFunctionWidget**: Color and opacity mapping UI

### Rendering Pipeline

1. Load geometry/volume via libcvc file I/O
2. Convert to VTK data structures (vtkPolyData, vtkImageData)
3. Create appropriate mappers (vtkPolyDataMapper, vtkSmartVolumeMapper)
4. Add actors/volumes to VTK renderer
5. Scene graph manages visibility and updates
6. Camera controller handles user input
7. Transfer function widget controls volume appearance

### State Management

VolRover3 uses a reactive state management system built on `cvc::state`:

- **AppState**: Singleton managing application-wide state
  - Camera position, view direction, FOV
  - Geometry and volume data
  - World bounds and visibility flags
  - Transfer function parameters
  
- **State Tree**: All state stored in hierarchical tree at `volrover3.*`
  - Direct access: `cvc::state::instance()("volrover3")("camera_position_x")`
  - Bidirectional synchronization with AppState methods
  
- **Change Notifications**: Register callbacks for reactive updates
  - All callback methods return `boost::signals2::connection`
  - Disconnect when no longer needed for proper lifecycle management
  - See `docs/APPSTATE_CALLBACKS.md` for detailed API documentation

## API Documentation

- [AppState Callback System](docs/APPSTATE_CALLBACKS.md) - Reactive state change notifications
- [Multi-Object Graphics System](docs/GRAPHICS_SYSTEM.md) - SceneGraph / GraphicsNode architecture
- [GraphicsNode Data-Driven Updates](docs/GRAPHICS_DATA_DRIVEN_UPDATES.md) - State-tree-driven geometry reloads
- [libcvc Testing Guide](https://github.com/transfix/libcvc/blob/master/docs/TESTING.md) - Unit and integration test documentation (libcvc repo)

## Future Enhancements

- [ ] Isosurface extraction and rendering
- [ ] Multiple geometry/volume layers
- [ ] Advanced transfer function editor with histogram
- [ ] Screenshot and animation export
- [ ] Property inspector for loaded data
- [ ] Clipping planes
- [ ] Lighting controls
- [ ] Material editor
- [ ] Measurements and annotations

## License

Copyright © 2025 CVC (Computational Visualization Center)

See main project LICENSE for details.
