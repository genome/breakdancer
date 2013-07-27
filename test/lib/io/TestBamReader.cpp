#include "io/BamReader.hpp"

#include "io/AlignmentFilter.hpp"
#include "io/RawBamEntry.hpp"

#include "TestData.hpp"

#include <gtest/gtest.h>

namespace bdaf = breakdancer::alnfilter;

class TestBamReader : public ::testing::TestWithParam<BamInfo> {
};

INSTANTIATE_TEST_CASE_P(RC, TestBamReader, ::testing::ValuesIn(TEST_BAMS));

TEST_P(TestBamReader, read_count) {
    std::string const& path = GetParam().path;
    size_t expected_count = GetParam().n_reads;

    BamReader<bdaf::True> reader(path);
    EXPECT_EQ(path, reader.path());

    RawBamEntry b;
    size_t n_reads = 0;
    while (reader.next(b) > 0) {
        ++n_reads;
    }

    ASSERT_EQ(expected_count, n_reads);
}
