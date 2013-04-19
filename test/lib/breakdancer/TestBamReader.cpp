#include "breakdancer/BamReader.hpp"

#include <boost/filesystem.hpp>
#include <unistd.h>
#include <cstdlib>
#include <gtest/gtest.h>


using namespace std;
namespace bfs = boost::filesystem;

class TestRegionLimitedBamReader : public ::testing::Test {
    protected:
        void SetUp() {
            _bam_path = getenv("BD_TEST_BAM_FILE");
            _temp_dir = "/tmp/BD_TestBamReader.XXXXXX";
            mkdtemp(&_temp_dir[0]);
        }


        void TearDown() {
            bfs::remove_all(_temp_dir);
        }

        char const* _bam_path;
        string _temp_dir;
};

TEST_F(TestRegionLimitedBamReader, readChromosome) {
    if (!_bam_path) {
        return;
    }

    //RegionLimitedBamReader reader(_bam_path, "21");
    BamReader reader(_bam_path);

    bfs::path output(_temp_dir);
    output /= "result.sam";
    samfile_t* sam_out = samopen(output.string().c_str(), "w", reader.header());
    ASSERT_TRUE(sam_out);

    bam1_t* b = bam_init1();
    size_t n = 0;
    while (reader.next(b) > 0) {
        ++n;
        samwrite(sam_out, b);
    }

    cerr << "Read " << n << " records\n";

    samclose(sam_out);
    bam_destroy1(b);
}
