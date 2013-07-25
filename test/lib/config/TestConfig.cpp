#include "config/BamConfig.hpp"

#include "common/Options.hpp"

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
    }

protected:
    std::stringstream _cfg_stream;
    Options _opts;
};

TEST(ConfigEntry, translate_token) {
    map<string, ConfigField> expected_translations = map_list_of
        ("map", BAM_FILE)

        ("lib", LIBRARY_NAME)
        ("libname", LIBRARY_NAME)
        ("library_name", LIBRARY_NAME)

        ("groUp", READ_GROUP)
        ("ReadgroUp", READ_GROUP)
        ("Read_groUp", READ_GROUP)

        ("mean", INSERT_SIZE_MEAN)
        ("mean_insert", INSERT_SIZE_MEAN)
        ("mean_insert_size", INSERT_SIZE_MEAN)

        ("std", INSERT_SIZE_STDDEV)
        ("stddev", INSERT_SIZE_STDDEV)
        ("insert_stddev", INSERT_SIZE_STDDEV)
        ("insert_size_stddev", INSERT_SIZE_STDDEV)
        ("stddev_insert", INSERT_SIZE_STDDEV)
        ("stddev_insert_size", INSERT_SIZE_STDDEV)

        ("readlen", READ_LENGTH)
        ("rEaDlEnGtH", READ_LENGTH)
        ("average_readlen", READ_LENGTH)
        ("average_readlength", READ_LENGTH)

        ("upp", INSERT_SIZE_UPPER_CUTOFF)
        ("upper", INSERT_SIZE_UPPER_CUTOFF)
        ("uppEr_cutOff", INSERT_SIZE_UPPER_CUTOFF)
        ("inseRt_size_uPper_cutoff", INSERT_SIZE_UPPER_CUTOFF)

        ("low", INSERT_SIZE_LOWER_CUTOFF)
        ("lower", INSERT_SIZE_LOWER_CUTOFF)
        ("lower_cuToff", INSERT_SIZE_LOWER_CUTOFF)
        ("insert_size_lower_cutoff", INSERT_SIZE_LOWER_CUTOFF)

        ("mapqual", MIN_MAP_QUAL)
        ("mapPing_quAlity", MIN_MAP_QUAL)

        ("samp", SAMPLE_NAME)
        ("sample", SAMPLE_NAME)
        ("samplename", SAMPLE_NAME)
        ("sample_name", SAMPLE_NAME)
        ;

    EXPECT_EQ(UNKNOWN, ConfigEntry::translate_token("ZIOJFksfjlaiaowinfd"));

    typedef map<string, ConfigField>::const_iterator TIter;
    for (TIter i = expected_translations.begin(); i != expected_translations.end(); ++i) {
        string src = i->first;
        EXPECT_EQ(i->second, ConfigEntry::translate_token(src))
            << "Translating " << src;
        boost::to_upper(src);
        EXPECT_EQ(i->second, ConfigEntry::translate_token(src))
            << "Translating " << src;
        boost::to_lower(src);
        EXPECT_EQ(i->second, ConfigEntry::translate_token(src))
            << "Translating " << src;
    }
}

TEST_F(TestConfig, legacyParse) {
    BamConfig cfg(_cfg_stream, _opts);

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
    EXPECT_EQ(0u, cfg.library_config_by_name("lib1").index);
    EXPECT_EQ(1u, cfg.library_config_by_name("lib2").index);
    EXPECT_THROW(cfg.library_config_by_name("lib3"), out_of_range);

    LibraryConfig const& lc1 = cfg.library_config_by_index(0);
    EXPECT_EQ(0u, lc1.index);
    EXPECT_EQ("lib1", lc1.name);
    EXPECT_EQ("x.bam", lc1.bam_file);
    EXPECT_FLOAT_EQ(90.00, lc1.readlens);
    EXPECT_FLOAT_EQ(277.03, lc1.lowercutoff);
    EXPECT_FLOAT_EQ(525.50, lc1.uppercutoff);
    EXPECT_FLOAT_EQ(467.59, lc1.mean_insertsize);
    EXPECT_FLOAT_EQ(31.91, lc1.std_insertsize);
    EXPECT_EQ(-1, lc1.min_mapping_quality);

    LibraryConfig const& lc2 = cfg.library_config_by_index(1);
    EXPECT_EQ(1u, lc2.index);
    EXPECT_EQ("lib2", lc2.name);
    EXPECT_EQ("y.bam", lc2.bam_file);
    EXPECT_FLOAT_EQ(90.00, lc2.readlens);
    EXPECT_FLOAT_EQ(311.36, lc2.lowercutoff);
    EXPECT_FLOAT_EQ(532.53, lc2.uppercutoff);
    EXPECT_FLOAT_EQ(475.76, lc2.mean_insertsize);
    EXPECT_FLOAT_EQ(28.67, lc2.std_insertsize);
    EXPECT_EQ(10, lc2.min_mapping_quality);

    EXPECT_THROW(cfg.library_config_by_index(2), out_of_range);
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

    ASSERT_NO_THROW(BamConfig(cfgss, _opts));
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

    ASSERT_THROW(BamConfig(cfgss, _opts), runtime_error);
}
