#pragma once

#include "BasicRegion.hpp"
#include "BedWriter.hpp" // FIXME: try to move this to io lib
#include "ReadCountsByLib.hpp"
#include "ReadRegionData.hpp"
#include "common/Timer.hpp"
#include "io/FastqWriter.hpp"

#include <boost/scoped_ptr.hpp>
#include <boost/unordered_map.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <iterator>
#include <map>
#include <numeric>
#include <set>
#include <stdint.h>
#include <string>
#include <vector>

class Options;
class BamConfig; //FIXME we probably won't need to store this once done stubbing in LibraryInfo
class LibraryInfo;
class BamReaderBase;

class BreakDancer {
public:
    typedef BasicRegion::ReadVector ReadVector;
    typedef std::vector<BasicRegion*> RegionData;
    typedef std::vector<ReadCountsByLib> RoiReadCounts;
    typedef ReadRegionData::ReadsToRegionsMap ReadsToRegionsMap;

    BreakDancer(
        Options const& opts,
        BamConfig const& cfg,
        LibraryInfo const& lib_info,
        ReadRegionData& read_regions,
        BamReaderBase& merged_reader,
        int max_read_window_size
        );

    void push_read(Read& aln);
    void build_connection();


    void process_sv(std::vector<int> const& snodes);
    void set_max_read_window_size(int val) {
        _max_read_window_size = val;
    }

    void process_breakpoint();
    void process_final_region();

    void dump_fastq(ReadFlag const& flag, std::vector<Read> const& support_reads);

    void run();

    void set_read_density(std::string const& libName, float density);

private:
    uint32_t _region_lib_counts(size_t region_idx, std::string const& lib, RoiReadCounts const& x) const {
        if (region_idx >= x.size())
            return 0;
        RoiReadCounts::value_type::const_iterator found = x[region_idx].find(lib);
        if (found != x[region_idx].end())
            return found->second;
        return 0;
    }

private: // data
    Options const& _opts;
    BamConfig const& _cfg;
    LibraryInfo const& _lib_info;
    ReadRegionData& _rdata;
    BamReaderBase& _merged_reader;
    int _max_read_window_size;

    bool _collecting_normal_reads;
    int _nnormal_reads;
    int _ntotal_nucleotides;
    int _max_readlen;
    int _buffer_size;

    int _region_start_tid;
    int _region_start_pos;
    int _region_end_tid; // global (chr, should be int in samtools)
    int _region_end_pos; // global

    ReadVector reads_in_current_region;
    boost::scoped_ptr<FastqWriter> _fastq_writer;
    boost::scoped_ptr<std::ofstream> _bed_stream;
    boost::scoped_ptr<BedWriter> _bed_writer;

    std::map<std::string, float> _read_density;
};
