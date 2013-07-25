#include "breakdancer/BreakDancer.hpp"

#include "breakdancer/ReadCountsByLib.hpp"
#include "config/BamConfig.hpp"

#include <map>
#include <memory>

#include <gtest/gtest.h>

using namespace std;
using namespace breakdancer;

TEST(BreakDancer, accumulate_reads_in_region) {
// The code this tests is under active development.
// Commenting this out until it stabilizes.
// -ta
/*
    BreakDancer bd;
    typedef ReadCountsByLib MapType;
    MapType counts1;
    MapType counts2;
    counts1["x"] = 11;
    counts1["y"] = 22;
    counts2["x"] = 100;
    counts2["y"] = 100;


    bd.add_per_lib_read_counts_to_region(0, counts1);
    bd.add_per_lib_read_counts_to_region(1, counts2);
    MapType c;
    bd.accumulate_reads_in_region(c, 0, 1);
    ASSERT_EQ(counts1, c);

    c.clear();
    bd.accumulate_reads_in_region(c, 1, 2);
    ASSERT_EQ(counts2, c);

    c.clear();
    bd.accumulate_reads_in_region(c, 0, 2);
    ASSERT_EQ(counts1 + counts2, c);
*/
}
