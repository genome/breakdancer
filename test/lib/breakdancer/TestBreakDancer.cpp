#include "breakdancer/BreakDancer.hpp"

#include "breakdancer/BamConfig.hpp"
#include "breakdancer/ReadCountsByLib.hpp"

#include <map>
#include <memory>

#include <gtest/gtest.h>

using namespace std;
using namespace breakdancer;

TEST(BreakDancer, region_lib_read_counts) {
    BreakDancer bd;
    typedef ReadCountsByLib MapType;
    MapType counts;
    counts["x"] = 11;
    counts["y"] = 22;

    bd.add_per_lib_read_counts_to_region(0, counts);
    MapType const* c0 = bd.region_read_counts_by_library(0);
    ASSERT_TRUE(c0 != 0);
    ASSERT_TRUE(counts == *c0);

    counts["z"] = 33;
    bd.add_per_lib_read_counts_to_region(0, counts);
    ASSERT_EQ(3, c0->size());

    MapType expected;
    expected["x"] = 22;
    expected["y"] = 44;
    expected["z"] = 33;

    ASSERT_TRUE(expected == *c0);

    ASSERT_EQ(22, bd.region_lib_read_count(0, "x"));
    ASSERT_EQ(44, bd.region_lib_read_count(0, "y"));
    ASSERT_EQ(33, bd.region_lib_read_count(0, "z"));
}
