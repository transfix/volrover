# cvcpkg/recipes/volrover2/build.ps1 - build VolumeRover 2.0 from the
# in-repo source tree. Mirrors the flags proven by the repo's CI
# (.github/workflows/ci.yml / release.yml).
$ErrorActionPreference = 'Stop'

if (-not $env:CVC_SOURCE_DIR) { throw 'CVC_SOURCE_DIR must be set' }
if (-not $env:CVC_BUILD_DIR) { throw 'CVC_BUILD_DIR must be set' }
if (-not $env:CVC_INSTALL_DIR) { throw 'CVC_INSTALL_DIR must be set' }
if (-not $env:CVC_BUILD_TYPE) { $env:CVC_BUILD_TYPE = 'Release' }
if (-not $env:CVC_JOBS) { $env:CVC_JOBS = [Environment]::ProcessorCount }

$cmakeBuildType = if ($env:CVC_BUILD_TYPE.ToLower() -eq 'debug') { 'Debug' } else { 'Release' }

$args = @(
  '-G', 'Ninja',
  '-S', $env:CVC_SOURCE_DIR,
  '-B', $env:CVC_BUILD_DIR,
  "-DCMAKE_INSTALL_PREFIX=$env:CVC_INSTALL_DIR",
  "-DCMAKE_BUILD_TYPE=$cmakeBuildType",
  # Clean version string (no git hash) for published packages
  '-DVOLROVER_RELEASE=ON',
  '-DUSE_QT6=ON',
  '-DDISABLE_CGAL=ON',
  # volrover3 is a separate recipe (cvcpkg/recipes/volrover3); the
  # VolumeRover 2.0 package never enables it.
  '-DBUILD_VOLROVER3=OFF'
)

if ($env:CVC_DEPS_PREFIX) {
  $args += "-DCMAKE_PREFIX_PATH=$env:CVC_DEPS_PREFIX"
}

& cmake @args
if ($LASTEXITCODE -ne 0) { throw 'cmake configure failed' }

& cmake --build $env:CVC_BUILD_DIR -j $env:CVC_JOBS
if ($LASTEXITCODE -ne 0) { throw 'cmake build failed' }

& cmake --install $env:CVC_BUILD_DIR
if ($LASTEXITCODE -ne 0) { throw 'cmake install failed' }
