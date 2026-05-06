/*
 * Tests for VolMagick compatibility layer.
 * Verifies that VolMagick:: types, subclasses, exception classes,
 * and free functions all work correctly through the compat headers.
 */

#include <gtest/gtest.h>

#include <VolMagick/VolMagick.h>
#include <VolMagick/Types.h>
#include <VolMagick/Dimension.h>
#include <VolMagick/BoundingBox.h>
#include <VolMagick/Voxels.h>
#include <VolMagick/Volume.h>
#include <VolMagick/Exceptions.h>
#include <VolMagick/endians.h>

#include <cvc/types.h>
#include <cvc/volume.h>
#include <cvc/voxels.h>

#include <type_traits>

// ---------------------------------------------------------------------------
// Type aliases
// ---------------------------------------------------------------------------
TEST(VolMagick_Types, VoxelTypeAliasMatchesCvc)
{
    static_assert(std::is_same<VolMagick::VoxelType, cvc::data_type>::value,
                  "VolMagick::VoxelType must alias cvc::data_type");
}

TEST(VolMagick_Types, EnumValuesAccessible)
{
    EXPECT_EQ(VolMagick::UChar,     cvc::UChar);
    EXPECT_EQ(VolMagick::UShort,    cvc::UShort);
    EXPECT_EQ(VolMagick::UInt,      cvc::UInt);
    EXPECT_EQ(VolMagick::Float,     cvc::Float);
    EXPECT_EQ(VolMagick::Double,    cvc::Double);
    EXPECT_EQ(VolMagick::Char,      cvc::Char);
    EXPECT_EQ(VolMagick::Undefined, cvc::Undefined);
}

TEST(VolMagick_Types, VoxelTypeSizesCorrect)
{
    ASSERT_NE(VolMagick::VoxelTypeSizes, nullptr);
    EXPECT_EQ(VolMagick::VoxelTypeSizes[VolMagick::UChar],  1u);
    EXPECT_EQ(VolMagick::VoxelTypeSizes[VolMagick::UShort], 2u);
    EXPECT_EQ(VolMagick::VoxelTypeSizes[VolMagick::Float],  4u);
    EXPECT_EQ(VolMagick::VoxelTypeSizes[VolMagick::Double], 8u);
}

TEST(VolMagick_Types, DimensionAlias)
{
    static_assert(std::is_same<VolMagick::Dimension, cvc::dimension>::value, "");
}

TEST(VolMagick_Types, BoundingBoxAlias)
{
    static_assert(std::is_same<VolMagick::BoundingBox, cvc::bounding_box>::value, "");
}

// ---------------------------------------------------------------------------
// VolMagick::Voxels subclass (dimension() compat)
// ---------------------------------------------------------------------------
TEST(VolMagick_Voxels, InheritsFromCvcVoxels)
{
    static_assert(std::is_base_of<cvc::voxels, VolMagick::Voxels>::value,
                  "VolMagick::Voxels must inherit from cvc::voxels");
}

TEST(VolMagick_Voxels, DefaultConstruction)
{
    VolMagick::Voxels v;
    // Default-constructed voxels have dimensions (4,4,4)
    EXPECT_EQ(v.voxel_dimensions()[0], 4u);
    EXPECT_EQ(v.voxel_dimensions()[1], 4u);
    EXPECT_EQ(v.voxel_dimensions()[2], 4u);
}

TEST(VolMagick_Voxels, DimensionCompatMethod)
{
    VolMagick::Voxels v;
    VolMagick::Dimension dim(8, 16, 32);
    v.dimension(dim);

    EXPECT_EQ(v.dimension()[0], 8u);
    EXPECT_EQ(v.dimension()[1], 16u);
    EXPECT_EQ(v.dimension()[2], 32u);

    // dimension() and voxel_dimensions() must return the same data
    EXPECT_EQ(v.dimension()[0], v.voxel_dimensions()[0]);
    EXPECT_EQ(v.dimension()[1], v.voxel_dimensions()[1]);
    EXPECT_EQ(v.dimension()[2], v.voxel_dimensions()[2]);
}

// ---------------------------------------------------------------------------
// VolMagick::Volume subclass
// ---------------------------------------------------------------------------
TEST(VolMagick_Volume, InheritsFromCvcVolume)
{
    static_assert(std::is_base_of<cvc::volume, VolMagick::Volume>::value,
                  "VolMagick::Volume must inherit from cvc::volume");
}

TEST(VolMagick_Volume, DefaultConstruction)
{
    VolMagick::Volume vol;
    // Default volume has dimensions (4,4,4)
    EXPECT_EQ(vol.voxel_dimensions()[0], 4u);
    EXPECT_EQ(vol.voxel_dimensions()[1], 4u);
    EXPECT_EQ(vol.voxel_dimensions()[2], 4u);
}

TEST(VolMagick_Volume, DimensionCompatMethod)
{
    VolMagick::Volume vol;
    VolMagick::Dimension dim(4, 8, 16);
    vol.dimension(dim);

    EXPECT_EQ(vol.dimension()[0], 4u);
    EXPECT_EQ(vol.dimension()[1], 8u);
    EXPECT_EQ(vol.dimension()[2], 16u);
    EXPECT_EQ(vol.dimension()[0], vol.voxel_dimensions()[0]);
}

TEST(VolMagick_Volume, ConstructFromCvcVolume)
{
    static cvc::app test_app;
    cvc::volume cv(test_app);
    cv.voxelType(cvc::Float);
    cv.voxel_dimensions(cvc::dimension(4, 4, 4));

    VolMagick::Volume v(cv);
    EXPECT_EQ(v.dimension()[0], 4u);
    EXPECT_EQ(v.voxelType(), cvc::Float);
}

TEST(VolMagick_Volume, CopyConstruction)
{
    VolMagick::Volume v1;
    v1.voxelType(cvc::UChar);
    v1.dimension(VolMagick::Dimension(2, 2, 2));

    VolMagick::Volume v2(v1);
    EXPECT_EQ(v2.dimension()[0], 2u);
    EXPECT_EQ(v2.voxelType(), cvc::UChar);
}

// ---------------------------------------------------------------------------
// VolMagick::VolumeFileInfo subclass
// ---------------------------------------------------------------------------
TEST(VolMagick_VolumeFileInfo, InheritsFromCvc)
{
    static_assert(
        std::is_base_of<cvc::volume_file_info, VolMagick::VolumeFileInfo>::value,
        "VolMagick::VolumeFileInfo must inherit from cvc::volume_file_info");
}

// ---------------------------------------------------------------------------
// VolMagick exceptions
// ---------------------------------------------------------------------------
TEST(VolMagick_Exceptions, ReadErrorMessage)
{
    VolMagick::ReadError err("file not found");
    std::string what = err.what_str();
    EXPECT_TRUE(what.find("ReadError") != std::string::npos)
        << "what_str() = " << what;
    EXPECT_TRUE(what.find("file not found") != std::string::npos)
        << "what_str() = " << what;
}

TEST(VolMagick_Exceptions, WriteErrorMessage)
{
    VolMagick::WriteError err("disk full");
    EXPECT_TRUE(err.what_str().find("disk full") != std::string::npos);
}

TEST(VolMagick_Exceptions, DefaultMessage)
{
    VolMagick::ReadError err;
    EXPECT_TRUE(err.what_str().find("ReadError") != std::string::npos);
}

TEST(VolMagick_Exceptions, CatchAsBaseException)
{
    try {
        throw VolMagick::IndexOutOfBounds("idx 99");
    } catch (const VolMagick::Exception& e) {
        EXPECT_TRUE(e.what_str().find("idx 99") != std::string::npos);
    }
}

TEST(VolMagick_Exceptions, AllExceptionTypesInstantiate)
{
    // Verify all 10 exception types compile and have what_str()
    EXPECT_FALSE(VolMagick::ReadError("x").what_str().empty());
    EXPECT_FALSE(VolMagick::WriteError("x").what_str().empty());
    EXPECT_FALSE(VolMagick::MemoryAllocationError("x").what_str().empty());
    EXPECT_FALSE(VolMagick::SubVolumeOutOfBounds("x").what_str().empty());
    EXPECT_FALSE(VolMagick::UnsupportedVolumeFileType("x").what_str().empty());
    EXPECT_FALSE(VolMagick::IndexOutOfBounds("x").what_str().empty());
    EXPECT_FALSE(VolMagick::NullDimension("x").what_str().empty());
    EXPECT_FALSE(VolMagick::VolumePropertiesMismatch("x").what_str().empty());
    EXPECT_FALSE(VolMagick::VolumeCacheDirectoryFileError("x").what_str().empty());
    EXPECT_FALSE(VolMagick::InvalidBoundingBox("x").what_str().empty());
}

// ---------------------------------------------------------------------------
// VolMagick endians
// ---------------------------------------------------------------------------
TEST(VolMagick_Endians, BigEndianReturnsInt)
{
    int result = VolMagick::big_endian();
    // On x86/x86_64 this should be 0 (little endian)
    EXPECT_EQ(result, 0);
}
