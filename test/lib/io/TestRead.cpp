#include "io/Read.hpp"

#include <map>
#include <memory>

#include <gtest/gtest.h>

using namespace std;
using namespace breakdancer;

namespace {
    uint8_t data[23] ={
        0x6a, 0x75, 0x6e, 0x6b, 0x0, 0x20, 0x0, 0x0, 0x0, 0x28, 0x27, 0x21,
        0x52, 0x47, 0x5a, 0x72, 0x67, 0x33, 0x0, 0x41, 0x4d, 0x43, 0x25
    };
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
            core.isize = -478;
            bam_record.core = core;
            bam_record.l_aux = 11;
            bam_record.data_len = 23;
            bam_record.m_data = 32;
            bam_record.data = &data[0];
            test_read.reset(new Read(&bam_record));
        }

        bam1_core_t core;
        bam1_t bam_record;
        auto_ptr<Read> test_read;
};

TEST_F(TestRead, readgroup) {
    ASSERT_EQ(test_read->readgroup(), "rg3");
}

TEST_F(TestRead, query_name) {
    ASSERT_EQ(test_read->query_name(), "junk");
}

TEST_F(TestRead, query_sequence) {
    ASSERT_EQ(test_read->query_sequence(), "CT");
}

TEST_F(TestRead, quality_string) {
    ASSERT_EQ(test_read->quality_string(), "HB");
}

TEST_F(TestRead, ori) {
    ASSERT_EQ(test_read->ori(), FWD);
}

TEST_F(TestRead, tid) {
    ASSERT_EQ(test_read->tid(), 22);
}

TEST_F(TestRead, pos) {
    ASSERT_EQ(test_read->pos(), 29184911);
}

TEST_F(TestRead, query_length) {
    ASSERT_EQ(test_read->query_length(), 2);
}

TEST_F(TestRead, isize) {
    ASSERT_EQ(test_read->isize(), -478);
}

TEST_F(TestRead, abs_isize) {
    ASSERT_EQ(test_read->abs_isize(), 478);
}

TEST_F(TestRead, set_bdflag) {
    test_read->set_bdflag(breakdancer::ARP_CTX);
    ASSERT_EQ(test_read->bdflag(),breakdancer::ARP_CTX);
    test_read->set_bdflag(breakdancer::NORMAL_FR);
}

TEST_F(TestRead, copy_constructor) {
    Read test_copy(*test_read);
    ASSERT_EQ(test_copy.ori(), FWD);
    ASSERT_EQ(test_copy.quality_string(), "HB");
    ASSERT_EQ(test_copy.query_sequence(), "CT");
    ASSERT_EQ(test_copy.query_name(), "junk");
    ASSERT_EQ(test_copy.readgroup(), "rg3");
}

TEST_F(TestRead, assignment) {
    Read test_copy = *test_read;
    ASSERT_EQ(test_copy.ori(), FWD);
    ASSERT_EQ(test_copy.quality_string(), "HB");
    ASSERT_EQ(test_copy.query_sequence(), "CT");
    ASSERT_EQ(test_copy.query_name(), "junk");
    ASSERT_EQ(test_copy.readgroup(), "rg3");
}
