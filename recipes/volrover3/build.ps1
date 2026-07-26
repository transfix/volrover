# recipes/volrover3/build.ps1 - build the VolumeRover3 application.
#
# volrover3 consumes libcvc as an EXTERNAL SDK (find_package(cvc CONFIG)
# -> cvc::cvc), so CVC_DEPS_PREFIX must contain an installed libcvc SDK
# plus Qt6 and VTK. Only the volrover3 target is built and only the
# `volrover3` install component is staged; the VolumeRover 2.0 targets
# are untouched.
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
  '-DVOLROVER_RELEASE=ON',
  '-DUSE_QT6=ON',
  '-DDISABLE_CGAL=ON',
  '-DBUILD_VOLROVER3=ON',
  '-DVOLROVER3_BUILD_TESTS=OFF'
)

if ($env:CVC_DEPS_PREFIX) {
  $args += "-DCMAKE_PREFIX_PATH=$env:CVC_DEPS_PREFIX"
}

& cmake @args
if ($LASTEXITCODE -ne 0) { throw 'cmake configure failed' }

& cmake --build $env:CVC_BUILD_DIR -j $env:CVC_JOBS --target volrover3
if ($LASTEXITCODE -ne 0) { throw 'cmake build failed' }

# Stage only the volrover3 component into the install prefix.
& cmake --install $env:CVC_BUILD_DIR --component volrover3 --prefix $env:CVC_INSTALL_DIR
if ($LASTEXITCODE -ne 0) { throw 'cmake install failed' }
