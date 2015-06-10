#include "io/BamWriter.hpp"
#include "io/Alignment.hpp"
#include "io/BamIo.hpp"
#include "io/RawBamEntry.hpp"

#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>

#include <gtest/gtest.h>

#include <fstream>
#include <string>

namespace bfs = boost::filesystem;

namespace {
    std::string tempPath(std::string const& prefix, std::string const& suffix) {
        bfs::path tmpdir = bfs::temp_directory_path();
        std::string tmpl(prefix);
        tmpl += "%%%%-%%%%-%%%%-%%%%";
        tmpl += suffix;
        bfs::path rv = tmpdir / bfs::unique_path(tmpl);
        return rv.native();
    }

    // Our test sam data comprises three pairs aligned like this:
    // P1
    // GTTTT
    //      gccccttttt
    // P2
    // TGTTTT
    //      cgcccttttt
    // P3
    // TTGTTTTTTT
    //      ccgccttttt
    std::string samData =
        "@HD\tVN:1.0\tGO:none\tSO:coordinate\n"
        "@SQ\tSN:21\tLN:46944323\tUR:internet\tAS:spec\tM5:NA\tSP:unknown\n"

        // adjacent, but no overlap due to soft clipping
        "P1\t99\t21\t10\t60\t5M5S\t=\t15\t15\tGTTTTTTTTT\tHHHHHHHHHH\n"
        "P1\t147\t21\t15\t60\t10M\t=\t10\t-15\tGCCCCTTTTT\tHHHHHHHHHH\n"

        // 1bp overlap (involving soft clipping)
        "P2\t99\t21\t10\t60\t6M4S\t=\t15\t15\tTGTTTTTTTT\tHHHHHHHHHH\n"
        "P2\t147\t21\t15\t60\t10M\t=\t10\t-15\tCGCCCTTTTT\tHHHHHHHHHH\n"

        // 5bp overlap (no clipping)
        "P3\t99\t21\t10\t60\t10M\t=\t15\t15\tTTGTTTTTTT\tHHHHHHHHHH\n"
        "P3\t147\t21\t15\t60\t10M\t=\t10\t-15\tCCGCCTTTTT\tHHHHHHHHHH\n"
        ;
}

class TestBam : public ::testing::Test {
public:
    void SetUp() {
        samPath_ = tempPath("breakdancer-unit-test", ".sam");
        std::ofstream out(samPath_.c_str());
        out << samData;
        out.close();

        reader.reset(openBam(samPath_));

        boost::shared_ptr<RawBamEntry> entry(new RawBamEntry);
        while (reader->next(*entry) > 0) {
            reads.emplace_back(new Alignment(*entry));
            rawEntries.push_back(entry);
            entry.reset(new RawBamEntry);
        }

        EXPECT_EQ(6u, reads.size());
    }

    void TearDown() {
        bfs::remove(samPath_);
    }

protected:
    std::string samPath_;
    std::vector<Alignment::Ptr> reads;
    std::vector<boost::shared_ptr<RawBamEntry> > rawEntries;
    boost::shared_ptr<BamReaderBase> reader;
};

TEST_F(TestBam, leftmost) {

    // We only test for overlap on the first (leftmost) read
    EXPECT_TRUE (reads[0]->leftmost());
    EXPECT_FALSE(reads[1]->leftmost());
    EXPECT_TRUE (reads[2]->leftmost());
    EXPECT_FALSE(reads[3]->leftmost());
    EXPECT_TRUE (reads[4]->leftmost());
    EXPECT_FALSE(reads[5]->leftmost());
}


TEST_F(TestBam, bamWriter) {
    std::string tmpOut = tempPath("breakdancer-unit-test", ".sam");
    BamWriter writer(tmpOut, reader->header(), true);
    for (size_t i = 0; i < rawEntries.size(); ++i) {
        writer.write(*rawEntries[i]);
    }
    writer.close();

    std::ifstream in(tmpOut.c_str());
    in.seekg(0, std::ios::end);
    size_t size = in.tellg();
    in.seekg(0, std::ios::beg);
    std::vector<char> buf(size);
    ASSERT_TRUE((bool)in.read(buf.data(), size));
    buf.push_back(0); // make sure buffer is null terminated
    EXPECT_STREQ(samData.c_str(), buf.data());

    bfs::remove(tmpOut);
}
