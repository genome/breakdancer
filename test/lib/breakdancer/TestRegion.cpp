#include "breakdancer/Region.hpp"
#include <boost/assign/list_of.hpp>
#include <cmath>
#include <gtest/gtest.h>

using namespace std;
using namespace region;
using namespace boost::assign;

class TestRegion : public ::testing::Test {
    protected:
        void SetUp() {
            int chrom_tid = 1;
            int start = 500000;
            int end = 500100;
            int num_normal_reads = 20;
            Read r1 = list_of("r1")("some_val");
            Read r2 = list_of("r2")("other_val");
            vector<Read> fake_reads = list_of(r1)(r2);
            test_region = new Region(chrom_tid, start, end, num_normal_reads, fake_reads);
        }
        void TearDown() {
            delete test_region;
        }

        Region *test_region;
};

TEST_F(TestRegion, construction) {

    ASSERT_EQ(1, test_region->begins);
    ASSERT_EQ(500000, test_region->beginc);
    ASSERT_EQ(500100, test_region->lastc);
    ASSERT_EQ(20, test_region->num_normal_reads);
    ASSERT_EQ("r1", test_region->reads[0][0]);
    ASSERT_EQ("r2", test_region->reads[1][0]);

}

TEST_F(TestRegion, set_read_counts) {
    map<string, uint32_t> num_ROI_reads;
    num_ROI_reads["lib1"] = 2;
    num_ROI_reads["lib2"] = 1;

    map<string, uint32_t> num_FR_reads;
    num_FR_reads["lib1"] = 12;
    num_FR_reads["lib2"] = 5;
    num_FR_reads["lib3"] = 1;

    test_region->set_read_counts(num_ROI_reads, num_FR_reads);
    ASSERT_EQ(2,test_region->read_count_ROI_map["lib1"]);
    ASSERT_EQ(10,test_region->read_count_FR_map["lib1"]);
    ASSERT_EQ(1,test_region->read_count_FR_map["lib3"]);
#ifndef NDEBUG
    ASSERT_DEATH(test_region->set_read_counts(num_FR_reads, num_ROI_reads),".*Assertion.*failed.*");
#endif
}

