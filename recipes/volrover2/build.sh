#!/usr/bin/env bash
# recipes/volrover2/build.sh - build VolumeRover 2.0 from the
# in-repo source tree. Mirrors the flags proven by the repo's CI
# (.github/workflows/ci.yml / release.yml).
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
  # Clean version string (no git hash) for published packages
  -DVOLROVER_RELEASE=ON
  -DUSE_QT6=ON
  -DDISABLE_CGAL=ON
  # volrover3 is a separate recipe (recipes/volrover3); the
  # VolumeRover 2.0 package never enables it.
  -DBUILD_VOLROVER3=OFF
)

if [[ -n "${CVC_DEPS_PREFIX:-}" ]]; then
  CMAKE_ARGS+=("-DCMAKE_PREFIX_PATH=$CVC_DEPS_PREFIX")
  # Magick++-config (queried at configure time for ImageMagick 7
  # defines) resolves its flags through pkg-config.
  export PKG_CONFIG_PATH="$CVC_DEPS_PREFIX/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
fi

cmake "${CMAKE_ARGS[@]}"
cmake --build "$CVC_BUILD_DIR" -j "$CVC_JOBS"
cmake --install "$CVC_BUILD_DIR"
