/*
 * Integration tests for the volutils CLI.
 * Runs volutils as a subprocess and verifies output/exit codes.
 *
 * Requires VOLUTILS_BINARY to be defined at compile time,
 * pointing to the volutils executable path.
 */

#include <gtest/gtest.h>

#include <cstdlib>
#include <cstdio>
#include <string>
#include <array>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace fs = std::filesystem;

#ifndef VOLUTILS_BINARY
#error "VOLUTILS_BINARY must be defined"
#endif

// Run a shell command and capture stdout+stderr. Returns exit code.
static int run(const std::string& cmd, std::string& output)
{
    std::string full_cmd = cmd + " 2>&1";
    FILE* pipe = popen(full_cmd.c_str(), "r");
    if (!pipe) return -1;

    std::array<char, 4096> buf;
    output.clear();
    while (fgets(buf.data(), buf.size(), pipe))
        output += buf.data();

    int status = pclose(pipe);
    return WEXITSTATUS(status);
}

class VolutilsCLITest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        volutils_ = VOLUTILS_BINARY;
        tmpdir_ = fs::temp_directory_path() / "volutils_cli_test";
        fs::create_directories(tmpdir_);
    }

    void TearDown() override
    {
        fs::remove_all(tmpdir_);
    }

    std::string tmpFile(const std::string& name)
    {
        return (tmpdir_ / name).string();
    }

    // Helper: create a small test rawiv via "fill" command
    std::string createTestVolume(const std::string& name,
                                 unsigned int dim = 8,
                                 double value = 1.0)
    {
        std::string path = tmpFile(name);
        std::ostringstream cmd;
        cmd << volutils_ << " fill -o " << path
            << " --xdim " << dim
            << " --ydim " << dim
            << " --zdim " << dim
            << " --value " << value;
        std::string out;
        int rc = run(cmd.str(), out);
        EXPECT_EQ(rc, 0) << "fill failed: " << out;
        return path;
    }

    std::string volutils_;
    fs::path tmpdir_;
};

// --- Help/usage ---

TEST_F(VolutilsCLITest, NoArgsShowsHelp)
{
    std::string out;
    int rc = run(volutils_, out);
    EXPECT_NE(rc, 0); // should fail with no subcommand
}

TEST_F(VolutilsCLITest, HelpFlag)
{
    std::string out;
    int rc = run(volutils_ + " --help", out);
    EXPECT_EQ(rc, 0);
    EXPECT_TRUE(out.find("volutils") != std::string::npos)
        << "Help output should mention volutils";
}

// --- fill command ---

TEST_F(VolutilsCLITest, FillCreatesFile)
{
    std::string path = tmpFile("filled.rawiv");
    std::string out;
    int rc = run(volutils_ + " fill -o " + path +
                 " --xdim 4 --ydim 4 --zdim 4 --value 42.0", out);
    EXPECT_EQ(rc, 0) << out;
    EXPECT_TRUE(fs::exists(path));
    EXPECT_GT(fs::file_size(path), 0u);
}

// --- info command ---

TEST_F(VolutilsCLITest, InfoShowsMetadata)
{
    std::string vol = createTestVolume("info_test.rawiv", 8, 0.0);
    std::string out;
    int rc = run(volutils_ + " info " + vol, out);
    EXPECT_EQ(rc, 0) << out;
    // Should show dimension info
    EXPECT_TRUE(out.find("8") != std::string::npos)
        << "Info should show dimension 8. Output: " << out;
}

// --- stats command ---

TEST_F(VolutilsCLITest, StatsShowsStatistics)
{
    std::string vol = createTestVolume("stats_test.rawiv", 4, 5.0);
    std::string out;
    int rc = run(volutils_ + " stats " + vol, out);
    EXPECT_EQ(rc, 0) << out;
    // Should show min/max that include 5.0
    EXPECT_TRUE(out.find("5") != std::string::npos)
        << "Stats should show value 5. Output: " << out;
}

// --- convert command ---

TEST_F(VolutilsCLITest, ConvertChangesType)
{
    std::string vol = createTestVolume("conv_src.rawiv", 4, 1.0);
    std::string dst = tmpFile("conv_dst.rawiv");
    std::string out;
    int rc = run(volutils_ + " convert -i " + vol + " -o " + dst +
                 " --type UChar", out);
    EXPECT_EQ(rc, 0) << out;
    EXPECT_TRUE(fs::exists(dst));
}

// --- compare command ---

TEST_F(VolutilsCLITest, CompareSameVolumes)
{
    std::string vol = createTestVolume("cmp.rawiv", 4, 3.14);
    std::string out;
    int rc = run(volutils_ + " compare -i " + vol + " " + vol, out);
    EXPECT_EQ(rc, 0) << out;
    // Identical volumes should have 0 max difference
    EXPECT_TRUE(out.find("0") != std::string::npos)
        << "Compare of identical volumes should show 0 diff. Output: " << out;
}

// --- add command ---

TEST_F(VolutilsCLITest, AddTwoVolumes)
{
    std::string v1 = createTestVolume("add1.rawiv", 4, 2.0);
    std::string v2 = createTestVolume("add2.rawiv", 4, 3.0);
    std::string dst = tmpFile("add_out.rawiv");
    std::string out;
    int rc = run(volutils_ + " add -i " + v1 + " " + v2 + " -o " + dst, out);
    EXPECT_EQ(rc, 0) << out;
    EXPECT_TRUE(fs::exists(dst));
}

// --- scale command ---

TEST_F(VolutilsCLITest, ScaleVolume)
{
    std::string vol = createTestVolume("scale_src.rawiv", 4, 1.0);
    std::string dst = tmpFile("scale_out.rawiv");
    std::string out;
    int rc = run(volutils_ + " scale -i " + vol + " -o " + dst +
                 " --factor 2.5", out);
    EXPECT_EQ(rc, 0) << out;
    EXPECT_TRUE(fs::exists(dst));
}

// --- negate command ---

TEST_F(VolutilsCLITest, NegateVolume)
{
    std::string vol = createTestVolume("neg_src.rawiv", 4, 1.0);
    std::string dst = tmpFile("neg_out.rawiv");
    std::string out;
    int rc = run(volutils_ + " negate -i " + vol + " -o " + dst, out);
    EXPECT_EQ(rc, 0) << out;
    EXPECT_TRUE(fs::exists(dst));
}

// --- bbox-shift command ---

TEST_F(VolutilsCLITest, BboxShift)
{
    std::string vol = createTestVolume("bbox_src.rawiv", 4, 0.0);
    std::string dst = tmpFile("bbox_shifted.rawiv");
    std::string out;
    int rc = run(volutils_ + " bbox-shift -i " + vol + " -o " + dst +
                 " --dx 1.0 --dy 2.0 --dz 3.0", out);
    EXPECT_EQ(rc, 0) << out;
    EXPECT_TRUE(fs::exists(dst));
}

// --- edge command ---

TEST_F(VolutilsCLITest, EdgeDetection)
{
    std::string vol = createTestVolume("edge_src.rawiv", 8, 1.0);
    std::string dst = tmpFile("edge_out.rawiv");
    std::string out;
    int rc = run(volutils_ + " edge -i " + vol + " -o " + dst, out);
    EXPECT_EQ(rc, 0) << out;
    EXPECT_TRUE(fs::exists(dst));
}

// --- bunny command ---

TEST_F(VolutilsCLITest, BunnyGeneratesGeometry)
{
    std::string dst = tmpFile("bunny.rawc");
    std::string out;
    int rc = run(volutils_ + " bunny -o " + dst, out);
    EXPECT_EQ(rc, 0) << out;
    EXPECT_TRUE(fs::exists(dst));
}

TEST_F(VolutilsCLITest, BunnyGeneratesVolume)
{
    std::string dst = tmpFile("bunny_vol.rawiv");
    std::string out;
    int rc = run(volutils_ + " bunny --volume -d 16 -o " + dst, out);
    EXPECT_EQ(rc, 0) << out;
    EXPECT_TRUE(fs::exists(dst));
    EXPECT_GT(fs::file_size(dst), 0u);
}

// --- error handling ---

TEST_F(VolutilsCLITest, MissingInputFileErrors)
{
    std::string out;
    int rc = run(volutils_ + " info /nonexistent/file.rawiv", out);
    EXPECT_NE(rc, 0);
}

TEST_F(VolutilsCLITest, UnknownSubcommandErrors)
{
    std::string out;
    int rc = run(volutils_ + " bogus_command", out);
    EXPECT_NE(rc, 0);
}

// --- downsample command ---

TEST_F(VolutilsCLITest, Downsample)
{
    std::string vol = createTestVolume("ds_src.rawiv", 8, 1.0);
    std::string dst = tmpFile("ds_out.rawiv");
    std::string out;
    int rc = run(volutils_ + " downsample -i " + vol + " -o " + dst, out);
    EXPECT_EQ(rc, 0) << out;
    EXPECT_TRUE(fs::exists(dst));
}

// --- clamp-min command ---

TEST_F(VolutilsCLITest, ClampMin)
{
    std::string vol = createTestVolume("clamp_src.rawiv", 4, -5.0);
    std::string dst = tmpFile("clamped.rawiv");
    std::string out;
    int rc = run(volutils_ + " clamp-min -i " + vol + " -o " + dst +
                 " --value 0.0", out);
    EXPECT_EQ(rc, 0) << out;
    EXPECT_TRUE(fs::exists(dst));
}
