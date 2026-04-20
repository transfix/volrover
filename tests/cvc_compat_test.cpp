/*
 * Tests for CVC compatibility headers.
 * Verifies that CamelCase type aliases, namespace aliasing, and
 * exception types all work correctly through the compat layer.
 */

#include <gtest/gtest.h>

#include <CVC/Namespace.h>
#include <CVC/Types.h>
#include <CVC/BoundingBox.h>
#include <CVC/Dimension.h>
#include <CVC/Exception.h>

#include <cvc/types.h>
#include <cvc/bounding_box.h>
#include <cvc/dimension.h>
#include <cvc/exception.h>

#include <type_traits>

// ---------------------------------------------------------------------------
// Namespace alias: CVC == cvc
// ---------------------------------------------------------------------------
TEST(CVC_Namespace, CVC_AliasesCvc)
{
    // CVC::BoundingBox and cvc::bounding_box should be the exact same type
    static_assert(std::is_same<CVC::BoundingBox, cvc::bounding_box>::value,
                  "CVC::BoundingBox must alias cvc::bounding_box");
}

// ---------------------------------------------------------------------------
// CVC/Types.h CamelCase aliases
// ---------------------------------------------------------------------------
TEST(CVC_Types, DataTypeAliasMatchesCvcDataType)
{
    static_assert(std::is_same<CVC::DataType, cvc::data_type>::value,
                  "CVC::DataType must alias cvc::data_type");
}

TEST(CVC_Types, EnumValuesAccessibleViaCVC)
{
    EXPECT_EQ(CVC::UChar,   cvc::UChar);
    EXPECT_EQ(CVC::UShort,  cvc::UShort);
    EXPECT_EQ(CVC::Float,   cvc::Float);
    EXPECT_EQ(CVC::Double,  cvc::Double);
    EXPECT_EQ(CVC::Undefined, cvc::Undefined);
}

TEST(CVC_Types, DataTypeSizesNotNull)
{
    ASSERT_NE(CVC::DataTypeSizes, nullptr);
    // UChar should be 1 byte
    EXPECT_EQ(CVC::DataTypeSizes[CVC::UChar], 1u);
    // Float should be 4 bytes
    EXPECT_EQ(CVC::DataTypeSizes[CVC::Float], 4u);
    // Double should be 8 bytes
    EXPECT_EQ(CVC::DataTypeSizes[CVC::Double], 8u);
}

TEST(CVC_Types, DataTypeStringsNotNull)
{
    ASSERT_NE(CVC::DataTypeStrings, nullptr);
    EXPECT_STREQ(CVC::DataTypeStrings[CVC::UChar], cvc::data_type_strings[cvc::UChar]);
}

TEST(CVC_Types, AppFrameworkTypesExist)
{
    // Just verify these typedefs compile and refer to the correct types
    static_assert(std::is_same<CVC::DataMap, cvc::data_map>::value, "");
    static_assert(std::is_same<CVC::PropertyMap, cvc::property_map>::value, "");
    static_assert(std::is_same<CVC::ThreadPtr, cvc::thread_ptr>::value, "");
    static_assert(std::is_same<CVC::MutexPtr, cvc::mutex_ptr>::value, "");
    static_assert(std::is_same<CVC::Signal, cvc::signal>::value, "");
}

// ---------------------------------------------------------------------------
// CVC/BoundingBox.h
// ---------------------------------------------------------------------------
TEST(CVC_BoundingBox, TypeAliasMatchesCvc)
{
    static_assert(std::is_same<CVC::BoundingBox, cvc::bounding_box>::value, "");
    static_assert(std::is_same<CVC::IndexBoundingBox, cvc::index_bounding_box>::value, "");
}

TEST(CVC_BoundingBox, DefaultConstruction)
{
    CVC::BoundingBox bb;
    EXPECT_DOUBLE_EQ(bb.minx, 0.0);
    EXPECT_DOUBLE_EQ(bb.miny, 0.0);
    EXPECT_DOUBLE_EQ(bb.minz, 0.0);
    EXPECT_DOUBLE_EQ(bb.maxx, 0.0);
    EXPECT_DOUBLE_EQ(bb.maxy, 0.0);
    EXPECT_DOUBLE_EQ(bb.maxz, 0.0);
}

TEST(CVC_BoundingBox, ParameterizedConstruction)
{
    CVC::BoundingBox bb(-1.0, -2.0, -3.0, 4.0, 5.0, 6.0);
    EXPECT_DOUBLE_EQ(bb.minx, -1.0);
    EXPECT_DOUBLE_EQ(bb.miny, -2.0);
    EXPECT_DOUBLE_EQ(bb.minz, -3.0);
    EXPECT_DOUBLE_EQ(bb.maxx, 4.0);
    EXPECT_DOUBLE_EQ(bb.maxy, 5.0);
    EXPECT_DOUBLE_EQ(bb.maxz, 6.0);
}

TEST(CVC_BoundingBox, GenericBoundingBoxIsTemplateAlias)
{
    CVC::GenericBoundingBox<float> bb;
    static_assert(
        std::is_same<CVC::GenericBoundingBox<float>,
                     cvc::generic_bounding_box<float>>::value, "");
    EXPECT_FLOAT_EQ(bb.minx, 0.0f);
}

// ---------------------------------------------------------------------------
// CVC/Dimension.h
// ---------------------------------------------------------------------------
TEST(CVC_Dimension, TypeAliasMatchesCvc)
{
    static_assert(std::is_same<CVC::Dimension, cvc::dimension>::value, "");
}

TEST(CVC_Dimension, DefaultConstruction)
{
    CVC::Dimension dim;
    EXPECT_EQ(dim[0], 0u);
    EXPECT_EQ(dim[1], 0u);
    EXPECT_EQ(dim[2], 0u);
}

TEST(CVC_Dimension, ValueConstruction)
{
    CVC::Dimension dim(64, 128, 256);
    EXPECT_EQ(dim[0], 64u);
    EXPECT_EQ(dim[1], 128u);
    EXPECT_EQ(dim[2], 256u);
}

// ---------------------------------------------------------------------------
// CVC/Exception.h
// ---------------------------------------------------------------------------
TEST(CVC_Exception, BaseExceptionTypeAlias)
{
    static_assert(std::is_same<CVC::Exception, cvc::exception>::value, "");
}

TEST(CVC_Exception, ReadErrorTypeAlias)
{
    static_assert(std::is_same<CVC::ReadError, cvc::read_error>::value, "");
}

TEST(CVC_Exception, WriteErrorTypeAlias)
{
    static_assert(std::is_same<CVC::WriteError, cvc::write_error>::value, "");
}

TEST(CVC_Exception, VolroverSpecificExceptions)
{
    // These are defined locally via CVC_DEF_EXCEPTION, not aliases
    try {
        throw CVC::SubVolumeOutOfBounds("test msg");
    } catch (const CVC::Exception& e) {
        // Should be caught as CVC::Exception (base class)
        std::string what = e.what_str();
        EXPECT_TRUE(what.find("test msg") != std::string::npos)
            << "what_str() = " << what;
    }

    try {
        throw CVC::NetworkError("net error");
    } catch (const CVC::Exception& e) {
        std::string what = e.what_str();
        EXPECT_TRUE(what.find("net error") != std::string::npos)
            << "what_str() = " << what;
    }
}
