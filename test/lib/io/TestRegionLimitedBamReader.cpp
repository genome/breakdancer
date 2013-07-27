#include "io/RegionLimitedBamReader.hpp"

#include "io/AlignmentFilter.hpp"
#include "io/BamReader.hpp"
#include "io/RawBamEntry.hpp"

#include "TestData.hpp"

#include <boost/lexical_cast.hpp>

#include <gtest/gtest.h>

#include <string>
#include <sstream>
#include <stdexcept>

namespace bdaf = breakdancer::alnfilter;
using boost::lexical_cast;
using namespace std;

namespace {
    struct RegionCount {
        string seq;
        int begin;
        int end;

        size_t region_read_count;
        vector<string> read_names;
    };

    string make_region(string const& seq, int begin, int end) {
        stringstream ss;
        ss << seq << ":" << begin << "-" << end;
        return ss.str();
    }

    RegionCount getTestRegion(string const& bam_path, size_t max_reads) {
        BamReader<bdaf::True> in(bam_path);
        RawBamEntry b;

        // Get sequence id (tid) of the first read in the bam
        if (in.next(b) <= 0)
            throw std::runtime_error("Failed to read from " + bam_path);

        // So far we have processed one read
        size_t n_reads = 1;

        RegionCount rv;

        int tid = b->core.tid;

        // set sequence name, begin and initial end position in return struct.
        rv.seq = in.sequence_name(tid);
        rv.begin = b->core.pos;
        rv.end = rv.begin;

        // don't forget to collect the name of the first read!
        rv.read_names.push_back(bam1_qname(b));

        // XXX: We're assuming that we have at least > max_reads reads on the
        // first tid in the bam.
        while (n_reads < max_reads && in.next(b) > 0 && b->core.tid == tid) {
            ++n_reads;
            rv.end = b->core.pos;
            rv.read_names.push_back(bam1_qname(b));
        }

        // seek forward until we find a read that starts beyond our ending
        // position. again, we're assuming that such a read exists.
        while (in.next(b) > 0 && b->core.tid == tid && b->core.pos == rv.end) {
            ++n_reads;
            rv.end = b->core.pos;
            rv.read_names.push_back(bam1_qname(b));
        }

        // Now, the idea is that we should have n_reads == max_reads collected
        rv.region_read_count = n_reads;

        return rv;
    }
}

class TestRegionLimitedBamReader : public ::testing::TestWithParam<BamInfo> {
};

TEST_P(TestRegionLimitedBamReader, read_count) {
    string const& path = GetParam().path;
    RegionCount rc = getTestRegion(path, 104);

    // Samtools region parsing uses one based coordinates, but the C api
    // (which produced rc.begin, rc.end) uses one based coordinates.
    // We correct for this by adding one to begin and end.
    string region = make_region(rc.seq, rc.begin + 1, rc.end + 1);

    vector<string> observed_read_names;

    RegionLimitedBamReader<bdaf::True> reader(path, region.c_str());
    RawBamEntry b;

    size_t observed_read_count = 0;
    int first_pos = -1;
    int last_pos = -1;
    while (reader.next(b) > 0) {
        if (first_pos == -1)
            first_pos = b->core.pos;

        last_pos = b->core.pos;

        observed_read_names.push_back(bam1_qname(b));
        ++observed_read_count;
    }

    EXPECT_EQ(rc.region_read_count, observed_read_count)
        << "bam/region " << path << " " << region;

    EXPECT_EQ(rc.begin, first_pos);
    EXPECT_EQ(rc.end, last_pos);
    EXPECT_EQ(rc.read_names.size(), observed_read_names.size());
    EXPECT_EQ(rc.read_names, observed_read_names);
}

INSTANTIATE_TEST_CASE_P(RC, TestRegionLimitedBamReader,
    ::testing::ValuesIn(TEST_BAMS));
