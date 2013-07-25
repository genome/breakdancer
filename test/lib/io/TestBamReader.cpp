#include "common/Options.hpp"
#include "io/AlignmentFilter.hpp"
#include "io/BamReader.hpp"

#include "TestData.hpp"

#include <gtest/gtest.h>

namespace bdaf = breakdancer::alnfilter;

namespace {
    Options const default_options;
}

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


