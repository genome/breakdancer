#include "common/Options.hpp"
#include "io/AlignmentFilter.hpp"
#include "io/BamReaderBase.hpp"
#include "io/BamIo.hpp"

#include "TestData.hpp"

#include <boost/shared_ptr.hpp>

#include <gtest/gtest.h>

#include <vector>
#include <string>

namespace bdaf = breakdancer::alnfilter;
using namespace std;

class TestBamIo : public ::testing::Test {
protected:
    void SetUp() {
        for (size_t i = 0; i < TEST_BAMS.size(); ++i) {
            _bam_paths.push_back(TEST_BAMS[i].path);
        }
    }

    Options const default_options;
    std::vector<std::string> _bam_paths;
};

TEST_F(TestBamIo, openBam) {
    for (size_t i = 0; i < _bam_paths.size(); ++i) {
        boost::shared_ptr<BamReaderBase> reader(openBam(_bam_paths[i], default_options));
        EXPECT_EQ(_bam_paths[i], reader->path());
    }
}

TEST_F(TestBamIo, openBams) {
    vector<boost::shared_ptr<BamReaderBase> > readers = openBams(_bam_paths, default_options);
    EXPECT_EQ(_bam_paths.size(), readers.size());
    for (size_t i = 0; i < TEST_BAMS.size(); ++i) {
        EXPECT_EQ(_bam_paths[i], readers[i]->path());
    }
}
