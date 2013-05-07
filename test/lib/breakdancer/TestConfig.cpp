#include "breakdancer/Options.hpp"
#include "breakdancer/BDConfig.hpp"
#include "breakdancer/BamConfig.hpp"

#include <sstream>
#include <string>
#include <set>

#include <gtest/gtest.h>

using namespace std;

namespace {
    // FIXME: figure out how to generate these or something...
    string _cfg_str =
        "readgroup:rg1	platform:illumina	map:x.bam	readlen:90.00	lib:lib1	num:10001	lower:277.03	upper:525.50	mean:467.59	std:31.91	SWnormality:minus infinity	exe:samtools view\n"
        "readgroup:rg2	platform:illumina	map:x.bam	readlen:90.00	lib:lib1	num:10001	lower:277.03	upper:525.50	mean:467.59	std:31.91	SWnormality:minus infinity	exe:samtools view\n"
        "readgroup:rg3	platform:illumina	map:x.bam	readlen:90.00	lib:lib1	num:10001	lower:277.03	upper:525.50	mean:467.59	std:31.91	SWnormality:minus infinity	exe:samtools view\n"
        "readgroup:rg4	platform:illumina	map:x.bam	readlen:90.00	lib:lib1	num:10001	lower:277.03	upper:525.50	mean:467.59	std:31.91	SWnormality:minus infinity	exe:samtools view\n"
        "readgroup:rg5	platform:illumina	map:x.bam	readlen:90.00	lib:lib1	num:10001	lower:277.03	upper:525.50	mean:467.59	std:31.91	SWnormality:minus infinity	exe:samtools view\n"
        "readgroup:rg6	platform:illumina	map:x.bam	readlen:90.00	lib:lib1	num:10001	lower:277.03	upper:525.50	mean:467.59	std:31.91	SWnormality:minus infinity	exe:samtools view\n"
        "readgroup:rg7	platform:illumina	map:x.bam	readlen:90.00	lib:lib1	num:10001	lower:277.03	upper:525.50	mean:467.59	std:31.91	SWnormality:minus infinity	exe:samtools view\n"
        "readgroup:rg8	platform:illumina	map:y.bam	readlen:90.00	lib:lib2	num:10001	lower:311.36	upper:532.53	mean:475.76	std:28.67	SWnormality:minus infinity	exe:samtools view\n"
        "readgroup:rg9	platform:illumina	map:y.bam	readlen:90.00	lib:lib2	num:10001	lower:311.36	upper:532.53	mean:475.76	std:28.67	SWnormality:minus infinity	exe:samtools view\n"
        "readgroup:rg10	platform:illumina	map:y.bam	readlen:90.00	lib:lib2	num:10001	lower:311.36	upper:532.53	mean:475.76	std:28.67	SWnormality:minus infinity	exe:samtools view\n"
        "readgroup:rg11	platform:illumina	map:y.bam	readlen:90.00	lib:lib2	num:10001	lower:311.36	upper:532.53	mean:475.76	std:28.67	SWnormality:minus infinity	exe:samtools view\n"
        "readgroup:rg12	platform:illumina	map:y.bam	readlen:90.00	lib:lib2	num:10001	lower:311.36	upper:532.53	mean:475.76	std:28.67	SWnormality:minus infinity	exe:samtools view\n"
        "readgroup:rg13	platform:illumina	map:y.bam	readlen:90.00	lib:lib2	num:10001	lower:311.36	upper:532.53	mean:475.76	std:28.67	SWnormality:minus infinity	exe:samtools view\n"
        "readgroup:rg14	platform:illumina	map:y.bam	readlen:90.00	lib:lib2	num:10001	lower:311.36	upper:532.53	mean:475.76	std:28.67	SWnormality:minus infinity	exe:samtools view\n"
    ;
}

class TestConfig : public ::testing::Test {
public:
    void SetUp() {
        _cfg_stream << _cfg_str;
    }

protected:
    std::stringstream _cfg_stream;
    Options _opts;
};

class TestableBDConfig : public BDConfig {
public:
    TestableBDConfig(istream& in)
        : BDConfig(in)
    {
    }

    using BDConfig::_entries;
};

TEST_F(TestConfig, legacyParse) {
    BamConfig cfg(_cfg_stream, _opts);

    // test "fmaps", mapping input files -> first library they contain?
    ASSERT_EQ(2, cfg.fmaps.size());
    ASSERT_EQ(1, cfg.fmaps.count("x.bam"));
    ASSERT_EQ(1, cfg.fmaps.count("y.bam"));
    ASSERT_EQ("lib1", cfg.fmaps.find("x.bam")->second);
    ASSERT_EQ("lib2", cfg.fmaps.find("y.bam")->second);

    // test "exes", the executables perl would use to view the input file?
    ASSERT_EQ(2, cfg.exes.size());
    ASSERT_EQ(1, cfg.exes.count("x.bam"));
    ASSERT_EQ(1, cfg.exes.count("y.bam"));
    ASSERT_EQ("samtools view", cfg.exes.find("x.bam")->second);
    ASSERT_EQ("samtools view", cfg.exes.find("y.bam")->second);

    // test libmaps, mapping library names -> input files?
    ASSERT_EQ(2, cfg.library_info.size());
    ASSERT_EQ(1, cfg.library_info.count("lib1"));
    ASSERT_EQ(1, cfg.library_info.count("lib2"));
    ASSERT_EQ("x.bam", cfg.library_info.find("lib1")->second.bam_file);
    ASSERT_EQ("y.bam", cfg.library_info.find("lib2")->second.bam_file);

    // ...
}

TEST_F(TestConfig, bdConfigRawParse) {
    TestableBDConfig cfg(_cfg_stream);
    ASSERT_EQ(14, cfg._entries.size());

    set<string> expected_readgroups;
    for (size_t i = 0; i < 14; ++i) {
        stringstream rgname;
        rgname << "rg" << i+1;
        expected_readgroups.insert(rgname.str());
        ASSERT_EQ(rgname.str(), cfg._entries[i].readgroup);
    }

    ASSERT_TRUE(expected_readgroups == cfg.readgroups())
        << "Parsed read groups are not as expected";
}
