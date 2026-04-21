/*
  Compatibility header: forwards to libcvc volmagick aggregate header
  and provides VolMagick namespace aliases for existing code.
*/

#ifndef __VOLMAGICK_H__
#define __VOLMAGICK_H__

#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif

#include <cvc/volmagick.h>
#include <cvc/app.h>
#include <CVC/Namespace.h>
#include <VolMagick/Exceptions.h>

namespace VolMagick
{
  // Type aliases
  using cvc::int64;
  using cvc::uint64;
  typedef cvc::data_type VoxelType;

  static const unsigned int *VoxelTypeSizes = cvc::data_type_sizes;
  static const char **VoxelTypeStrings = cvc::data_type_strings;

  // Enum value aliases
  using cvc::Undefined;
  using cvc::UChar;
  using cvc::UShort;
  using cvc::UInt;
  using cvc::Float;
  using cvc::Double;
  using cvc::UInt64;
  using cvc::Char;

  // Type aliases
  typedef cvc::dimension      Dimension;
  typedef cvc::bounding_box   BoundingBox;
  typedef cvc::index_bounding_box IndexBoundingBox;
  typedef cvc::volume_file_io VolumeFile_IO;

  // Voxels subclass adding dimension() compat method
  class Voxels : public cvc::voxels
  {
  public:
    using cvc::voxels::voxels; // inherit constructors
    using cvc::voxels::voxel_dimensions;

    Voxels() : cvc::voxels(cvc::app::instance()) {}

    Dimension& dimension() { return voxel_dimensions(); }
    const Dimension& dimension() const { return voxel_dimensions(); }
    void dimension(const Dimension& d) { voxel_dimensions(d); }
  };

  // Volume subclass adding dimension() compat method  
  class Volume : public cvc::volume
  {
  public:
    using cvc::volume::volume; // inherit constructors
    using cvc::volume::voxel_dimensions;

    Volume() : cvc::volume(cvc::app::instance()) {}
    Volume(const cvc::volume& v) : cvc::volume(v) {}
    Volume(const Volume& v) : cvc::volume(v) {}
    Volume(const Dimension& d, cvc::data_type vt = cvc::UChar)
      : cvc::volume(cvc::app::instance(), d, vt) {}
    Volume(const Dimension& d, cvc::data_type vt, const BoundingBox& bb)
      : cvc::volume(cvc::app::instance(), d, vt, bb) {}

    Dimension& dimension() { return voxel_dimensions(); }
    const Dimension& dimension() const { return voxel_dimensions(); }
    void dimension(const Dimension& d) { voxel_dimensions(d); }
  };

  // VolumeFileInfo subclass adding dimension() compat method
  class VolumeFileInfo : public cvc::volume_file_info
  {
  public:
    using cvc::volume_file_info::volume_file_info;
    using cvc::volume_file_info::voxel_dimensions;

    Dimension& dimension() { return voxel_dimensions(); }
    const Dimension& dimension() const { return voxel_dimensions(); }
    void dimension(const Dimension& d) { voxel_dimensions(d); }
  };

  // Free function aliases
  using cvc::readVolumeFile;
  using cvc::writeVolumeFile;
  using cvc::createVolumeFile;
  using cvc::calcGradient;

  // Overload for std::vector<VolMagick::Volume> (subclass of cvc::volume)
  inline void writeVolumeFile(const std::vector<Volume>& vols,
                              const std::string& filename)
  {
    std::vector<cvc::volume> base_vols(vols.begin(), vols.end());
    cvc::writeVolumeFile(base_vols, filename);
  }

  // ctx-aware overload for std::vector<VolMagick::Volume>
  inline void writeVolumeFile(cvc::app& ctx,
                              const std::vector<Volume>& vols,
                              const std::string& filename)
  {
    std::vector<cvc::volume> base_vols(vols.begin(), vols.end());
    cvc::writeVolumeFile(ctx, base_vols, filename);
  }
}

#endif
