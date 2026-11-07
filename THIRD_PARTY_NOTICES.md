# Licensing and third-party notices

VolumeRover is licensed under the **GNU Lesser General Public License, version 2.1** (see [LICENSE](LICENSE)),
the same licence as libcvc. The source files carry a "GNU Lesser General Public License version 2.1" header,
with these copyright holders:

| Holder | Where | Files |
|---|---|---|
| The University of Texas at Austin | most of `src/` and `inc/` | 622 |
| Jesse Sweet (2010), Bharadwaj Subramanian (2009), written at CVC | `inc/SweetMesh/`, `src/SweetMesh/`, `inc/MMHLS/`, `src/MMHLS/` | 31 |
| The Technical University of Denmark (author Ojaswa Sharma) | `inc/MSLevelSet/`, `src/MSLevelSet/`, `inc/HigherOrderSegmentation/`, `src/HigherOrderSegmentation/`, `inc/MultiphaseSegmentation/`, `src/MultiphaseSegmentation/` | 19 |
| CyberPC Angel, LLC (2026) | `volrover3/` C++ and Python sources | 89 |

(Before this change the repository had no root licence file. 42 UT-held files under `inc/VolumeRover2/`,
`inc/ColorTable2/`, `src/ColorTable2/`, `inc/cvcalgo/ui/`, `src/cvcalgo/ui/`, `src/VolMagick/`,
`inc/RawIVEditor/`, `src/RawIVEditor/`, `inc/SimpleExample/`, `src/SimpleExample/` and `inc/MMHLS/cvcraw_geometry.h`,
and the 31 SweetMesh and MMHLS files, had GPL-2.0-or-later headers. 102 more UT files had an LGPL-2.1 grant
followed by GPL boilerplate. Those were mistakes. The DTU files and the `volrover3/` sources stated no licence.
All of them now carry the LGPL-2.1 header.)

## Files under the GNU General Public License

The files below carry GPL headers and are **not** covered by VolumeRover's LGPL grant; each remains under the
licence its header states. The table records whether the default build compiles them, so that anyone
distributing binaries knows which obligations apply.

| Component | Files | Copyright holder (per header) | Licence | In the default build? | Upstream |
|---|---|---|---|---|---|
| VRender (`inc/QGLViewer/VRender/`, `src/QGLViewer/VRender/`) | 28 | Cyril Soler (2005); each file also carries the libQGLViewer header (Gilles Debunne, LGPL-3.0-or-later) | GPL-2.0-or-later | **Yes**: compiled into the vendored QGLViewer library, which VolumeRover2 links, unless `-DNO_VECTORIAL_RENDER=ON` | VRender, http://artis.imag.fr/Members/Cyril.Soler/VRender, shipped inside libQGLViewer 3.0.0 |

**Consequence for binaries.** Because the default build links VRender into VolumeRover2, a VolumeRover2 binary
built with the default configuration is a combined work that must be distributed under the GPL's terms (and see
GPC below). Building with `-DNO_VECTORIAL_RENDER=ON` leaves VRender and GPC out. The VolumeRover sources
themselves remain LGPL-2.1. VolumeRover3 (`volrover3/`) links only libcvc, cvcGL, Python, yaml-cpp, Qt and VTK,
and none of the components above.

## Other third-party code

| Component | Path | Copyright holder | Licence | In the default build? |
|---|---|---|---|---|
| libQGLViewer 3.0.0 | `inc/QGLViewer/`, `src/QGLViewer/` | Gilles Debunne (2002–2025) | LGPL-3.0-or-later (https://github.com/GillesDebunne/libQGLViewer) | **Yes**: linked by VolumeRover2 |
| Generic Polygon Clipper (GPC) 2.32 | `inc/QGLViewer/VRender/gpc.h`, `src/QGLViewer/VRender/gpc.cpp` | Alan Murta, Advanced Interfaces Group, University of Manchester (1997–2004) | Free for non-commercial use only; commercial use needs the author's consent. Not an open-source licence, and not GPL-compatible. These two files also carry the VRender and libQGLViewer headers. | **Yes**: part of VRender |
| IMOD `libiimod` | `inc/VolMagick/libiimod/`, `src/VolMagick/libiimod/` | BL3DEMC and the Regents of the University of Colorado (1995–2009) | The headers refer to IMOD's `dist/COPYRIGHT`, which releases IMOD's C/C++ libraries under the LGPL (its `LGPL.txt` is version 3) and its programs under GPL-2.0 (https://bio3d.colorado.edu/imod/) | No: only with `-DVOLMAGICK_USING_IMOD_MRC=ON` |
| XmlRpc++ | `inc/XmlRPC/`, `src/XmlRPC/` | Chris Morley (2002–2003); CVC notice added at UT | Distributed here under the CVC LGPL notice; upstream XmlRpc++ (http://xmlrpcpp.sourceforge.net) is LGPL, version not checked | **Yes** |
| TinyXML | `src/ContourTiler/xml/tinyxml*.cpp`, `inc/ContourTiler/xml/tinyxml.h`, `inc/ContourTiler/xml/tinystr.h` | Lee Thomason | zlib | **Yes**: part of ContourTiler (`BUILD_TILING_LIB`, ON unless `DISABLE_CGAL`) |

## Other notices

- The `Segmentation` and `XmlRPC` sources carry the older CVC notice, which grants the code under "a Lesser General
  Public License (LGPL)" with no version number and asks that published work cite the listed references.
- `CMake/cuda/FindCUDA/` (NVIDIA, MIT) and `CMake/FindPETSc.cmake` (BSD) are build scripts under their own
  licences and are not part of any binary.
