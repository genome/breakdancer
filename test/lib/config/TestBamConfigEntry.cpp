#include "config/BamConfigEntry.hpp"

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

TEST(BamConfigEntry, translate_token) {
    typedef BamConfigEntry BCE;
    typedef BamConfigEntry::Field ConfigField;

    map<string, ConfigField> expected_translations = map_list_of
        ("map", BCE::BAM_FILE)

        ("lib", BCE::LIBRARY_NAME)
        ("libname", BCE::LIBRARY_NAME)
        ("library_name", BCE::LIBRARY_NAME)

        ("groUp", BCE::READ_GROUP)
        ("ReadgroUp", BCE::READ_GROUP)
        ("Read_groUp", BCE::READ_GROUP)

        ("mean", BCE::INSERT_SIZE_MEAN)
        ("mean_insert", BCE::INSERT_SIZE_MEAN)
        ("mean_insert_size", BCE::INSERT_SIZE_MEAN)

        ("std", BCE::INSERT_SIZE_STDDEV)
        ("stddev", BCE::INSERT_SIZE_STDDEV)
        ("insert_stddev", BCE::INSERT_SIZE_STDDEV)
        ("insert_size_stddev", BCE::INSERT_SIZE_STDDEV)
        ("stddev_insert", BCE::INSERT_SIZE_STDDEV)
        ("stddev_insert_size", BCE::INSERT_SIZE_STDDEV)

        ("readlen", BCE::READ_LENGTH)
        ("rEaDlEnGtH", BCE::READ_LENGTH)
        ("average_readlen", BCE::READ_LENGTH)
        ("average_readlength", BCE::READ_LENGTH)

        ("upp", BCE::INSERT_SIZE_UPPER_CUTOFF)
        ("upper", BCE::INSERT_SIZE_UPPER_CUTOFF)
        ("uppEr_cutOff", BCE::INSERT_SIZE_UPPER_CUTOFF)
        ("inseRt_size_uPper_cutoff", BCE::INSERT_SIZE_UPPER_CUTOFF)

        ("low", BCE::INSERT_SIZE_LOWER_CUTOFF)
        ("lower", BCE::INSERT_SIZE_LOWER_CUTOFF)
        ("lower_cuToff", BCE::INSERT_SIZE_LOWER_CUTOFF)
        ("insert_size_lower_cutoff", BCE::INSERT_SIZE_LOWER_CUTOFF)

        ("mapqual", BCE::MIN_MAP_QUAL)
        ("mapPing_quAlity", BCE::MIN_MAP_QUAL)

        ("samp", BCE::SAMPLE_NAME)
        ("sample", BCE::SAMPLE_NAME)
        ("samplename", BCE::SAMPLE_NAME)
        ("sample_name", BCE::SAMPLE_NAME)
        ;

    EXPECT_EQ(BCE::UNKNOWN, BamConfigEntry::translate_token("ZIOJFksfjlaiaowinfd"));

    typedef map<string, ConfigField>::const_iterator TIter;
    for (TIter i = expected_translations.begin(); i != expected_translations.end(); ++i) {
        string src = i->first;
        EXPECT_EQ(i->second, BamConfigEntry::translate_token(src))
            << "Translating " << src;
        boost::to_upper(src);
        EXPECT_EQ(i->second, BamConfigEntry::translate_token(src))
            << "Translating " << src;
        boost::to_lower(src);
        EXPECT_EQ(i->second, BamConfigEntry::translate_token(src))
            << "Translating " << src;
    }
}
