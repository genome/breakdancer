#include "io/RegionLimitedBamReader.hpp"

#include "io/AlignmentFilter.hpp"
#include "io/BamReader.hpp"
#include "io/RawBamEntry.hpp"

#include "TestData.hpp"

#include <boost/lexical_cast.hpp>

#include <gtest/gtest.h>

#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

using boost::lexical_cast;
using namespace std;

namespace {
    struct RegionCount {
        string seq;
        int begin;
        int end;

        size_t read_count;
        vector<string> read_names;
    };

    string make_region(string const& seq, int begin, int end) {
        stringstream ss;
        ss << seq << ":" << begin << "-" << end;
        return ss.str();
    }

    RegionCount getTestRegion(string const& bam_path, size_t max_reads) {
        BamReader<AlignmentFilter::True> in(bam_path);
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
        rv.read_count = n_reads;

        return rv;
    }
}

class TestRegionLimitedBamReader : public ::testing::TestWithParam<BamInfo> {
};

TEST_P(TestRegionLimitedBamReader, read_count) {
    string const& path = GetParam().path;
    RegionCount rc = getTestRegion(path, 104); // take at most 104 reads

    // Samtools region parsing uses one based coordinates, but the C api
    // (which produced rc.begin, rc.end) uses one based coordinates.
    // We correct for this by adding one to begin and end.
    string region = make_region(rc.seq, rc.begin + 1, rc.end + 1);

    vector<string> observed_read_names;

    RegionLimitedBamReader<AlignmentFilter::True> reader(path, region.c_str());
    EXPECT_EQ(path, reader.path());
    EXPECT_EQ(path + " (region: " + region + ")", reader.description());

    size_t observed_read_count = 0;
    int first_pos = -1;
    int last_pos = -1;
    RawBamEntry b;
    while (reader.next(b) > 0) {
        if (first_pos == -1)
            first_pos = b->core.pos;

        last_pos = b->core.pos;

        observed_read_names.push_back(bam1_qname(b));
        ++observed_read_count;
    }

    EXPECT_EQ(rc.read_count, observed_read_count)
        << "bam/region " << path << " " << region;

    EXPECT_EQ(rc.begin, first_pos);
    EXPECT_EQ(rc.end, last_pos);
    EXPECT_EQ(rc.read_names.size(), observed_read_names.size());
    EXPECT_EQ(rc.read_names, observed_read_names);
}


TEST_P(TestRegionLimitedBamReader, set_region) {
    string const& path = GetParam().path;
    RegionCount rc = getTestRegion(path, 10); // take at most 10 reads

    string region = make_region(rc.seq, rc.begin + 1, rc.end + 1);

    vector<string> observed_read_names;

    RegionLimitedBamReader<AlignmentFilter::True> reader(path, region.c_str());
    EXPECT_EQ(path, reader.path());
    EXPECT_EQ(path + " (region: " + region + ")", reader.description());

    size_t observed_read_count = 0;
    int first_pos = -1;
    int last_pos = -1;
    RawBamEntry b;
    std::string first_name;
    while (reader.next(b) > 0) {
        if (first_pos == -1) {
            first_pos = b->core.pos;
            first_name = bam1_qname(b);
        }

        last_pos = b->core.pos;

        observed_read_names.push_back(bam1_qname(b));
        ++observed_read_count;
    }

    EXPECT_EQ(rc.read_count, observed_read_count)
        << "bam/region " << path << " " << region;

    reader.set_region(region.c_str());
    int rv = reader.next(b);
    EXPECT_GT(rv, 0);
    EXPECT_EQ(first_name, bam1_qname(b));
}

TEST_P(TestRegionLimitedBamReader, multi_region) {
    string const& path = GetParam().path;
    RegionCount rc1 = getTestRegion(path, 10); // take at most 10 reads
    RegionCount rc2 = getTestRegion(path, 200); // take at most 10 reads
    std::string region1 = make_region(rc1.seq, rc1.begin + 1, rc1.end + 1);
    std::string region2 = make_region(rc2.seq, rc2.begin + 1, rc2.end + 1);

    std::size_t expected_count = rc1.read_count + rc2.read_count;

    std::vector<std::string> regions{region1, region2};
    MultiRegionLimitedBamReader<AlignmentFilter::True> reader(path, regions);
    RawBamEntry b;
    std::size_t observed_count = 0;

    std::multiset<std::string> expected(rc1.read_names.begin(), rc1.read_names.end());
    expected.insert(rc2.read_names.begin(), rc2.read_names.end());

    std::multiset<std::string> observed;
    while (reader.next(b) > 0) {
        char const* name = bam1_qname(b);
        observed.insert(name);
        ++observed_count;
    }

    EXPECT_EQ(expected_count, observed_count);
    EXPECT_EQ(expected, observed);
}

INSTANTIATE_TEST_CASE_P(RC, TestRegionLimitedBamReader,
    ::testing::ValuesIn(TEST_BAMS));
