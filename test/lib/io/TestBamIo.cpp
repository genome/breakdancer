#include "common/Options.hpp"
#include "io/AlignmentFilter.hpp"
#include "io/BamMerger.hpp"
#include "io/BamReader.hpp"

#include "TestData.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>

#include <gtest/gtest.h>

namespace bdaf = breakdancer::alnfilter;
namespace bfs = boost::filesystem;
using boost::shared_ptr;
using namespace std;

class TestBamReader : public ::testing::TestWithParam<BamInfo> {
};

TEST_P(TestBamReader, read_count) {
    BamReader<bdaf::True> reader(GetParam().path);
    EXPECT_EQ(reader.path(), GetParam().path);

    bam1_t* b = bam_init1();
    size_t n_reads = 0;
    while (reader.next(b) > 0) {
        ++n_reads;
    }
    bam_destroy1(b);

    ASSERT_EQ(GetParam().n_reads, n_reads);
}

INSTANTIATE_TEST_CASE_P(RC, TestBamReader, ::testing::ValuesIn(TEST_BAMS));


TEST(TestBamMerger, read_count) {
    vector<string> paths;
    size_t expected = 0;

    vector< shared_ptr<IBamReader> > spReaders;
    vector<IBamReader*> readers;

    for (size_t i = 0; i < TEST_BAMS.size(); ++i) {
        paths.push_back(TEST_BAMS[i].path);
        expected += TEST_BAMS[i].n_reads;
        shared_ptr<IBamReader> p(new BamReader<bdaf::True>(TEST_BAMS[i].path));
        spReaders.push_back(p);
        readers.push_back(p.get());
    }

    BamMerger reader(readers);

    int last_tid = -1;
    long last_pos = 0;

    bam1_t* b = bam_init1();
    size_t n_reads = 0;
    while (reader.next(b) > 0) {
        ASSERT_LE(last_tid, b->core.tid);
        ASSERT_LE(last_pos, b->core.pos);
        last_tid = b->core.tid;
        last_pos = b->core.pos;
        ++n_reads;
    }
    bam_destroy1(b);

    ASSERT_EQ(expected, n_reads);
}
