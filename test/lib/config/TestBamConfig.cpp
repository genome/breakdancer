#include "config/BamConfig.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/assign/list_of.hpp>
#include <gtest/gtest.h>

#include <sstream>
#include <string>
#include <set>

using boost::assign::map_list_of;
using namespace std;

namespace {
    const string _cfg_str =
        "readgroup:rg1	platform:illumina	map:x.bam	readlen:90.00	lib:lib1	num:10001	lower:277.03	upper:525.50	mean:467.59	std:31.91	SWnormality:minus infinity	exe:samtools view\n"
        "readgroup:rg2	platform:illumina	map:x.bam	readlen:90.00	lib:lib1	num:10001	lower:277.03	upper:525.50	mean:467.59	std:31.91	SWnormality:minus infinity	exe:samtools view\n"
        "readgroup:rg3	mapqual:10	platform:illumina	map:y.bam	readlen:90.00	lib:lib2	num:10001	lower:311.36	upper:532.53	mean:475.76	std:28.67	SWnormality:minus infinity	exe:samtools view\n"
        "readgroup:rg4	mapqual:10	platform:illumina	map:y.bam	readlen:90.00	lib:lib2	num:10001	lower:311.36	upper:532.53	mean:475.76	std:28.67	SWnormality:minus infinity	exe:samtools view\n"
    ;
}

class TestConfig : public ::testing::Test {
public:
    void SetUp() {
        _cfg_stream << _cfg_str;
        _cut_sd = 3;
    }

protected:
    std::stringstream _cfg_stream;
    int _cut_sd;
};

TEST_F(TestConfig, legacyParse) {
    BamConfig cfg(_cfg_stream, _cut_sd);

    // test "fmaps", mapping input files -> first library they contain?
    ASSERT_EQ(2u, cfg.bam_files().size());
    EXPECT_EQ("x.bam", cfg.bam_files()[0]);
    EXPECT_EQ("y.bam", cfg.bam_files()[1]);

    EXPECT_EQ("lib1", cfg.readgroup_library("rg1"));
    EXPECT_EQ("lib1", cfg.readgroup_library("rg2"));
    EXPECT_EQ("lib2", cfg.readgroup_library("rg3"));
    EXPECT_EQ("lib2", cfg.readgroup_library("rg4"));

    // test libmaps, mapping library names -> input files?
    ASSERT_EQ(2u, cfg.num_libs());
    EXPECT_EQ(0u, cfg.library_config("lib1").bam_file_index);
    EXPECT_EQ(0u, cfg.library_config("lib1").index);
    EXPECT_EQ(1u, cfg.library_config("lib2").bam_file_index);
    EXPECT_EQ(1u, cfg.library_config("lib2").index);
    EXPECT_THROW(cfg.library_config("lib3"), out_of_range);

    LibraryConfig const& lc1 = cfg.library_config(0);
    EXPECT_EQ(0u, lc1.index);
    EXPECT_EQ("lib1", lc1.name);
    EXPECT_EQ("x.bam", lc1.bam_file);
    EXPECT_FLOAT_EQ(90.00, lc1.readlens);
    EXPECT_FLOAT_EQ(277.03, lc1.lowercutoff);
    EXPECT_FLOAT_EQ(525.50, lc1.uppercutoff);
    EXPECT_FLOAT_EQ(467.59, lc1.mean_insertsize);
    EXPECT_FLOAT_EQ(31.91, lc1.std_insertsize);
    EXPECT_EQ(-1, lc1.min_mapping_quality);

    LibraryConfig const& lc2 = cfg.library_config(1);
    EXPECT_EQ(1u, lc2.index);
    EXPECT_EQ("lib2", lc2.name);
    EXPECT_EQ("y.bam", lc2.bam_file);
    EXPECT_FLOAT_EQ(90.00, lc2.readlens);
    EXPECT_FLOAT_EQ(311.36, lc2.lowercutoff);
    EXPECT_FLOAT_EQ(532.53, lc2.uppercutoff);
    EXPECT_FLOAT_EQ(475.76, lc2.mean_insertsize);
    EXPECT_FLOAT_EQ(28.67, lc2.std_insertsize);
    EXPECT_EQ(10, lc2.min_mapping_quality);

    EXPECT_THROW(cfg.library_config(2), out_of_range);
}

TEST_F(TestConfig, noPlatformOrExe) {
    stringstream cfgss(
        "readgroup:rg1"
        "	map:x.bam"
        "	readlen:90.00"
        "	lib:lib1"
        "	num:10001"
        "	lower:277.03"
        "	upper:525.50"
        "	mean:467.59"
        "	std:31.91\n"
        );

    ASSERT_NO_THROW(BamConfig(cfgss, _cut_sd));
}

TEST_F(TestConfig, missingBam) {
    stringstream cfgss(
        "readgroup:rg1"
        "	readlen:90.00"
        "	lib:lib1"
        "	num:10001"
        "	lower:277.03"
        "	upper:525.50"
        "	mean:467.59"
        "	std:31.91\n"
        );

    ASSERT_THROW(BamConfig(cfgss, _cut_sd), runtime_error);
}
