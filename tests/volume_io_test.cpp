/*
 * Tests for in-memory volume operations through the VolMagick compat layer.
 * Verifies construction, voxel access, dimension/bbox manipulation.
 *
 * NOTE: Direct file I/O tests (writeVolumeFile/readVolumeFile) are not included
 * here because libcvc's volume_file_io handler registration uses static
 * initializers that deadlock when loaded via shared library in test binaries.
 * Volume I/O is tested end-to-end through the volutils CLI tests instead.
 */

#include <gtest/gtest.h>

#include <VolMagick/VolMagick.h>
#include <CVC/BoundingBox.h>
#include <CVC/Dimension.h>

#include <cmath>

// ---------------------------------------------------------------------------
// Volume construction and property tests
// ---------------------------------------------------------------------------

TEST(VolumeOpsTest, DefaultVolumeProperties)
{
    VolMagick::Volume vol;

    // Default: 4x4x4 UChar with bbox [-0.5, 0.5]^3
    EXPECT_EQ(vol.dimension()[0], 4u);
    EXPECT_EQ(vol.dimension()[1], 4u);
    EXPECT_EQ(vol.dimension()[2], 4u);
    EXPECT_EQ(vol.voxelType(), VolMagick::UChar);
    EXPECT_DOUBLE_EQ(vol.boundingBox().minx, -0.5);
    EXPECT_DOUBLE_EQ(vol.boundingBox().maxx,  0.5);
}

TEST(VolumeOpsTest, SetVoxelTypeAndDimension)
{
    VolMagick::Volume vol;
    vol.voxelType(VolMagick::Float);
    vol.dimension(VolMagick::Dimension(8, 16, 32));

    EXPECT_EQ(vol.voxelType(), VolMagick::Float);
    EXPECT_EQ(vol.dimension()[0], 8u);
    EXPECT_EQ(vol.dimension()[1], 16u);
    EXPECT_EQ(vol.dimension()[2], 32u);
}

TEST(VolumeOpsTest, VoxelReadWrite)
{
    VolMagick::Volume vol;
    vol.voxelType(VolMagick::Float);
    vol.dimension(VolMagick::Dimension(4, 4, 4));

    // Write ascending values
    for (unsigned k = 0; k < 4; k++)
        for (unsigned j = 0; j < 4; j++)
            for (unsigned i = 0; i < 4; i++)
                vol(i, j, k, double(i + j * 4 + k * 16));

    // Read them back
    for (unsigned k = 0; k < 4; k++)
        for (unsigned j = 0; j < 4; j++)
            for (unsigned i = 0; i < 4; i++)
            {
                double expected = double(i + j * 4 + k * 16);
                EXPECT_FLOAT_EQ(vol(i, j, k), expected)
                    << "at (" << i << "," << j << "," << k << ")";
            }
}

TEST(VolumeOpsTest, BoundingBoxSetGet)
{
    VolMagick::Volume vol;
    VolMagick::BoundingBox bb(-10.5, -20.5, -30.5, 40.5, 50.5, 60.5);
    vol.boundingBox(bb);

    auto got = vol.boundingBox();
    EXPECT_DOUBLE_EQ(got.minx, -10.5);
    EXPECT_DOUBLE_EQ(got.miny, -20.5);
    EXPECT_DOUBLE_EQ(got.minz, -30.5);
    EXPECT_DOUBLE_EQ(got.maxx,  40.5);
    EXPECT_DOUBLE_EQ(got.maxy,  50.5);
    EXPECT_DOUBLE_EQ(got.maxz,  60.5);
}

TEST(VolumeOpsTest, MinMaxTracking)
{
    VolMagick::Volume vol;
    vol.voxelType(VolMagick::Float);
    vol.dimension(VolMagick::Dimension(4, 4, 4));

    for (unsigned k = 0; k < 4; k++)
        for (unsigned j = 0; j < 4; j++)
            for (unsigned i = 0; i < 4; i++)
                vol(i, j, k, double(i + j + k));

    // min should be 0 (at 0,0,0), max should be 9 (at 3,3,3)
    EXPECT_NEAR(vol.min(), 0.0, 1e-6);
    EXPECT_NEAR(vol.max(), 9.0, 1e-6);
}

TEST(VolumeOpsTest, VoxelTypeConversionPreservesData)
{
    VolMagick::Volume vol;
    vol.voxelType(VolMagick::Float);
    vol.dimension(VolMagick::Dimension(4, 4, 4));

    for (unsigned k = 0; k < 4; k++)
        for (unsigned j = 0; j < 4; j++)
            for (unsigned i = 0; i < 4; i++)
                vol(i, j, k, 100.0);

    // Convert to Double
    vol.voxelType(VolMagick::Double);
    EXPECT_EQ(vol.voxelType(), VolMagick::Double);
    EXPECT_NEAR(vol(0, 0, 0), 100.0, 1e-3);
}

TEST(VolumeOpsTest, ResizeClearsPreviousData)
{
    VolMagick::Volume vol;
    vol.voxelType(VolMagick::Float);
    vol.dimension(VolMagick::Dimension(4, 4, 4));

    vol(0, 0, 0, 42.0);
    EXPECT_FLOAT_EQ(vol(0, 0, 0), 42.0);

    // Resize should allocate new storage
    vol.dimension(VolMagick::Dimension(8, 8, 8));
    EXPECT_EQ(vol.dimension()[0], 8u);
}
