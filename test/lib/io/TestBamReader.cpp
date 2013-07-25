#include "io/BamReader.hpp"
#include "io/AlignmentFilter.hpp"

#include "TestData.hpp"

#include <boost/filesystem.hpp>
#include <unistd.h>
#include <cstdlib>
#include <gtest/gtest.h>

using namespace std;
namespace bfs = boost::filesystem;
namespace bdaf = breakdancer::alnfilter;

class TestBamReader : public ::testing::TestWithParam<BamInfo> {
protected:
    void SetUp() {
        _bam_path = TEST_DATA_DIRECTORY;
        _bam_path /= GetParam().path;
    }

    bfs::path _bam_path;
};

TEST_P(TestBamReader, read_count) {
    BamReader<bdaf::True> reader(_bam_path.string());
    EXPECT_EQ(reader.path(), _bam_path.string());

    bam1_t* b = bam_init1();
    size_t n_reads = 0;
    while (reader.next(b) > 0) {
        ++n_reads;
    }
    bam_destroy1(b);

    ASSERT_EQ(GetParam().n_reads, n_reads);
}

INSTANTIATE_TEST_CASE_P(RC, TestBamReader, ::testing::ValuesIn(TEST_BAMS));
