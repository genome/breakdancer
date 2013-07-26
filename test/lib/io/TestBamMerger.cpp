#include "io/AlignmentFilter.hpp"
#include "io/BamReader.hpp"
#include "io/BamMerger.hpp"
#include "io/RawBamEntry.hpp"

#include "TestData.hpp"

#include <boost/shared_ptr.hpp>

#include <gtest/gtest.h>

namespace bdaf = breakdancer::alnfilter;
using namespace std;

TEST(TestBamMerger, read_count) {
    vector<string> paths;
    size_t expected = 0;

    vector< boost::shared_ptr<BamReaderBase> > spReaders;
    vector<BamReaderBase*> readers;

    for (size_t i = 0; i < TEST_BAMS.size(); ++i) {
        paths.push_back(TEST_BAMS[i].path);
        expected += TEST_BAMS[i].n_reads;
        boost::shared_ptr<BamReaderBase> p(new BamReader<bdaf::True>(TEST_BAMS[i].path));
        spReaders.push_back(p);
        readers.push_back(p.get());
    }

    BamMerger reader(readers);

    int last_tid = -1;
    long last_pos = 0;

    RawBamEntry b;
    size_t n_reads = 0;
    while (reader.next(b) > 0) {
        ASSERT_LE(last_tid, b->core.tid);
        ASSERT_LE(last_pos, b->core.pos);
        last_tid = b->core.tid;
        last_pos = b->core.pos;
        ++n_reads;
    }

    ASSERT_EQ(expected, n_reads);
}


