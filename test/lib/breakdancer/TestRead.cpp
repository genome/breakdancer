#include "breakdancer/Read.hpp"
#include <map>
#include <gtest/gtest.h>

using namespace std;
using namespace breakdancer;

namespace {
    uint8_t data[23] ={0x6a, 0x75, 0x6e, 0x6b, 0x0, 0x20, 0x0, 0x0, 0x0, 0x28, 0x27, 0x21, 0x52, 0x47, 0x5a, 0x72, 0x67, 0x33, 0x0, 0x41, 0x4d, 0x43, 0x25}; 
}

class TestRead : public ::testing::Test {
    protected:
        void SetUp() {
            core.tid = 22;
            core.pos = 29184911;
            core.bin = 6462;
            core.qual = 60;
            core.l_qname = 5;
            core.flag = 163;
            core.n_cigar = 1;
            core.l_qseq = 2;
            core.mtid = 22;
            core.mpos = 29185299;
            core.isize = 478;
            bam_record.core = core;
            bam_record.l_aux = 11;
            bam_record.data_len = 23;
            bam_record.m_data = 32;
            bam_record.data = &data[0];
            readgroup_platform["rg3"] = "helicos";
            readgroup_library["rg3"] = "some_lib";
            test_read = new Read(&bam_record, "sam", readgroup_platform, readgroup_library);
        }
        void TearDown() {
            delete test_read;
        }
        bam1_core_t core;
        bam1_t bam_record;
        Read* test_read;
        map<string, string> readgroup_platform;
        map<string, string> readgroup_library;
};

TEST_F(TestRead, readgroup) {
    ASSERT_EQ(test_read->readgroup, "rg3");
}

TEST_F(TestRead, queryname) {
    ASSERT_EQ(test_read->queryname(), "junk");
}

TEST_F(TestRead, query_sequence) {
    ASSERT_EQ(test_read->query_sequence(), "CT");
}

TEST_F(TestRead, quality_string) {
    ASSERT_EQ(test_read->quality_string(), "HB");
}

TEST_F(TestRead, ori) {
    ASSERT_EQ(test_read->ori(), "+");
}

TEST_F(TestRead, indexing) {
    ASSERT_EQ((*test_read)[0], "junk");
    ASSERT_EQ((*test_read)[1], "22");
    ASSERT_EQ((*test_read)[2], "29184911");
    ASSERT_EQ((*test_read)[3], "+");
    ASSERT_EQ((*test_read)[4], "478");
    ASSERT_EQ((*test_read)[5], "18");
    ASSERT_EQ((*test_read)[6], "37");
    ASSERT_EQ((*test_read)[7], "2");
    ASSERT_EQ((*test_read)[8], "some_lib");
    ASSERT_EQ((*test_read)[9], "CT");
    ASSERT_EQ((*test_read)[10], "HB");
}
