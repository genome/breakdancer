#include "breakdancer/Region.hpp"
#include <boost/assign/list_of.hpp>
#include <cmath>
#include <gtest/gtest.h>

using namespace std;
using namespace region;
using namespace boost::assign;

TEST(TestRegion, construction) {
    int chrom_tid = 1;
    int start = 500000;
    int end = 500100;
    int num_normal_reads = 20;
    Read r1 = list_of("r1")("some_val");
    Read r2 = list_of("r2")("other_val");
    vector<Read> fake_reads = list_of(r1)(r2);

    Region test(chrom_tid, start, end, num_normal_reads, fake_reads);
    ASSERT_EQ(1, test.begins);
    ASSERT_EQ(500000, test.beginc);
    ASSERT_EQ(500100, test.lastc);
    ASSERT_EQ(20, num_normal_reads);
    ASSERT_EQ("r1",test.reads[0][0]);
    ASSERT_EQ("r2",test.reads[1][0]);

}
