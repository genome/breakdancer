#include "io/LibraryFlagDistribution.hpp"

#include <gtest/gtest.h>

#include <sstream>

TEST(TestLibraryFlagDistribution, merge) {
    LibraryFlagDistribution lfd1;
    LibraryFlagDistribution lfd2;

    lfd1.read_count = 100;
    lfd2.read_count = 2000;

    // the read counts don't necessarily have to be the sum of
    // the counts within each lib.
    lfd1.read_counts_by_flag[0] = 10;
    lfd1.read_counts_by_flag[1] = 20;

    lfd2.read_counts_by_flag[1] = 30;
    lfd2.read_counts_by_flag[2] = 40;

    lfd1.merge(lfd2);

    EXPECT_EQ(2100u, lfd1.read_count);
    EXPECT_EQ(10u, lfd1.read_counts_by_flag[0]);
    EXPECT_EQ(50u, lfd1.read_counts_by_flag[1]);
    EXPECT_EQ(40u, lfd1.read_counts_by_flag[2]);
}
