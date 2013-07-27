#include "io/BamIo.hpp"

#include "io/BamReaderBase.hpp"
#include "io/RawBamEntry.hpp"

#include "TestData.hpp"

#include <boost/shared_ptr.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <iterator>
#include <vector>
#include <set>
#include <string>
#include <utility>

namespace bdaf = breakdancer::alnfilter;
using namespace std;

class TestBamIo : public ::testing::Test {
protected:
    void SetUp() {
        for (size_t i = 0; i < TEST_BAMS.size(); ++i) {
            _bam_paths.push_back(TEST_BAMS[i].path);
        }
    }

    std::vector<std::string> _bam_paths;
};

TEST_F(TestBamIo, openBam) {
    int tid = -1;
    pair<int, int> extents(-1, -1);
    set<int> positions;

    for (size_t i = 0; i < _bam_paths.size(); ++i) {
        boost::shared_ptr<BamReaderBase> reader(openBam(_bam_paths[i]));
        EXPECT_EQ(_bam_paths[i], reader->path());

        // Get the first record
        RawBamEntry entry;
        ASSERT_GT(reader->next(entry), 0);
        ASSERT_NE(-1, entry->core.tid);
        tid = entry->core.tid;

        while (reader->next(entry) > 0 && entry->core.tid == tid) {
            // keep track of all unique positions we have seen
            positions.insert(entry->core.pos);
        }

        // this won't work if there's only one starting position
        ASSERT_GT(positions.size(), 1u);

        // pick a position near the middle
        set<int>::const_iterator iter = positions.begin();
        advance(iter, positions.size() / 2);
        int halfway = *iter;

        // now let's try to read just that region
        stringstream region;
        // + 1 to correct for 0-based reader api vs 1-based region api
        // in samtools.
        region << reader->sequence_name(tid) << ":" << (halfway + 1);
        reader.reset(openBam(_bam_paths[i], region.str()));
        EXPECT_EQ(_bam_paths[i], reader->path());

        ASSERT_GT(reader->next(entry), 0);
        EXPECT_EQ(tid, entry->core.tid);
        ASSERT_GT(entry->core.pos, extents.first);
        ASSERT_LE(entry->core.pos, halfway);
    }
}

TEST_F(TestBamIo, openBams) {
    vector<boost::shared_ptr<BamReaderBase> > readers = openBams(_bam_paths);
    EXPECT_EQ(_bam_paths.size(), readers.size());
    for (size_t i = 0; i < TEST_BAMS.size(); ++i) {
        EXPECT_EQ(_bam_paths[i], readers[i]->path());
    }
}
