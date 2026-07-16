#!/usr/bin/env bash
# cvcpkg/recipes/volrover3/build.sh - build the VolumeRover3 application.
#
# volrover3 consumes libcvc as an EXTERNAL SDK (find_package(cvc CONFIG)
# -> cvc::cvc), so CVC_DEPS_PREFIX must contain an installed libcvc SDK
# plus Qt6 and VTK (e.g. a cvcpkg prefix with @cvc/libcvc, qt6, vtk).
#
# The top-level VolumeRover project is configured with BUILD_VOLROVER3=ON
# (strictly additive — the VolumeRover 2.0 targets are untouched), only
# the volrover3 target is built, and only the `volrover3` install
# component is staged, so the packaged prefix contains exactly the
# volrover3 app tree.
set -euo pipefail

: "${CVC_SOURCE_DIR:?CVC_SOURCE_DIR must be set}"
: "${CVC_BUILD_DIR:?CVC_BUILD_DIR must be set}"
: "${CVC_INSTALL_DIR:?CVC_INSTALL_DIR must be set}"

CVC_BUILD_TYPE="${CVC_BUILD_TYPE:-Release}"
CVC_JOBS="${CVC_JOBS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)}"

case "$(echo "$CVC_BUILD_TYPE" | tr '[:upper:]' '[:lower:]')" in
  debug) CMAKE_BUILD_TYPE=Debug ;;
  *) CMAKE_BUILD_TYPE=Release ;;
esac

CMAKE_ARGS=(
  -G Ninja
  -S "$CVC_SOURCE_DIR"
  -B "$CVC_BUILD_DIR"
  -DCMAKE_INSTALL_PREFIX="$CVC_INSTALL_DIR"
  -DCMAKE_BUILD_TYPE="$CMAKE_BUILD_TYPE"
  -DVOLROVER_RELEASE=ON
  -DUSE_QT6=ON
  -DDISABLE_CGAL=ON
  -DBUILD_VOLROVER3=ON
  -DVOLROVER3_BUILD_TESTS=OFF
)

if [[ -n "${CVC_DEPS_PREFIX:-}" ]]; then
  CMAKE_ARGS+=("-DCMAKE_PREFIX_PATH=$CVC_DEPS_PREFIX")
fi

cmake "${CMAKE_ARGS[@]}"
cmake --build "$CVC_BUILD_DIR" -j "$CVC_JOBS" --target volrover3

# Stage only the volrover3 component. The install-time --prefix override
# also sidesteps the top-level macOS quirk that appends
# `VolumeRover-<ver>.app/Contents` to the configured prefix, so
# VolumeRover3.app lands directly at the prefix root.
cmake --install "$CVC_BUILD_DIR" --component volrover3 --prefix "$CVC_INSTALL_DIR"
