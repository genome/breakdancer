#include "io/AlignmentFilter.hpp"
#include "io/BamReader.hpp"
#include "io/BamMerger.hpp"

#include "TestData.hpp"

#include <boost/shared_ptr.hpp>

#include <gtest/gtest.h>

namespace bdaf = breakdancer::alnfilter;
using boost::shared_ptr;
using namespace std;

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


