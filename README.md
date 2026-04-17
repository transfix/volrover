# VolumeRover 2.1

[![CI/CD](https://github.com/transfix/volrover/actions/workflows/ci.yml/badge.svg)](https://github.com/transfix/volrover/actions/workflows/ci.yml)

VolumeRover (a.k.a. VolRover) is interactive multi-purpose 3D image processing and visualization software developed at the [Computational Visualization Center (CVC)](https://cvc-lab.github.io/), Oden Institute for Computational Engineering and Sciences, The University of Texas at Austin.

**Project page:** https://cvc-lab.github.io/software/volumerover/

## Features

- Visualization of three-dimensional volumetric data of arbitrary size on commodity hardware
- Image contrast enhancement
- Filtering and noise reduction (bilateral, anisotropic diffusion, GDTV)
- Image segmentation
- Isocontouring (LBIE, FastContouring)
- Symmetry detection for virus maps
- Boundary-free image skeletonization
- Mesh generation and quality improvement
- Multi-tile server support for large-scale rendering
- Contour tiling

VolumeRover provides a unified interface to a number of CVC software packages including segmentation, contrast enhancement, secondary structure detection, and reconstruction.

## Building

### Requirements

- CMake 3.16+
- C++17 compiler
- Qt 6
- Boost (thread, filesystem, regex, date_time, program-options, and header-only libraries)
- log4cplus 2.x
- GLEW
- FFTW3
- GSL
- OpenGL
- ImageMagick (Magick++ 7.x)

### Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
cmake --install build
```

### Platform Notes

- **Linux:** Install dependencies via your package manager (e.g. `apt`, `dnf`).
- **macOS:** Install dependencies via [Homebrew](https://brew.sh/).
- **Windows:** Install dependencies via [vcpkg](https://vcpkg.io/). Qt 6 can be installed separately. Pass `-DCMAKE_TOOLCHAIN_FILE=[vcpkg-root]/scripts/buildsystems/vcpkg.cmake` to CMake.

## CI/CD

The project has CI/CD via GitHub Actions with matrix builds across Ubuntu, macOS, and Windows in both Debug and Release configurations. Packages are produced as DEB (Linux), DMG (macOS), and NSIS installer + ZIP (Windows).

## References

- Q. Zhang, R. Bettadapura, C. Bajaj. "Macromolecular Structure Modeling from 3DEM using VOLROVER 2.0." *Biopolymers*, 2012.
- C. Bajaj, S. Goswami, Q. Zhang. "Detection of Secondary and Supersecondary Structures of Proteins for Cryo-Electron Microscopy." *Journal of Structural Biology*, 177(2), 2012.
- C. Bajaj, Z. Yu, M. Auer. "Volumetric Feature Extraction and Visualization of Tomographic Molecular Imaging." *Journal of Structural Biology*, 144(1-2), 2003.

## License

Originally developed at CVC, The University of Texas at Austin. See individual source files for license details.

## Contact

- Project page: https://cvc-lab.github.io/software/volumerover/
- CVC: https://cvc-lab.github.io/
