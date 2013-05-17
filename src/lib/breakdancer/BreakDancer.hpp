#pragma once

#include "ReadCountsByLib.hpp"
#include "BasicRegion.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <iterator>
#include <map>
#include <numeric>
#include <stdint.h>
#include <string>
#include <vector>

class Options;
class BamConfig;
class IBamReader;

class BreakDancer {
public:
    typedef breakdancer::Read ReadType;
    typedef BasicRegion::ReadVector ReadVector;
    typedef std::vector<BasicRegion*> RegionData;
    typedef std::vector<ReadCountsByLib> RoiReadCounts;

    BreakDancer(
        Options const& opts,
        BamConfig const& cfg,
        IBamReader& merged_reader,
        int max_read_window_size
        );

    ~BreakDancer();

    void push_read(ReadType& aln, bam_header_t const* bam_header);
    void build_connection(bam_header_t const* bam_header);

    ReadCountsByLib nread_ROI; // global
    ReadCountsByLib nread_FR;    // global
    std::map<std::string, float> read_density;
    ReadVector reads_in_current_region;

    void add_per_lib_read_counts_to_last_region(ReadCountsByLib const& counts) {
        assert(num_regions() > 0);

        size_t region_idx = num_regions()-1;
        if (region_idx >= _read_count_ROI_map.size())
            _read_count_ROI_map.resize(2*(region_idx+1));

        _read_count_ROI_map[region_idx] +=  counts;
    }

    void add_current_read_counts_to_last_region() {
        assert(num_regions() > 0);

        size_t region_idx = num_regions() - 1;
        if (region_idx >= _read_count_ROI_map.size())
            _read_count_ROI_map.resize(2*(region_idx+1));

        _read_count_ROI_map[region_idx] = nread_ROI;

        if (region_idx >= _read_count_FR_map.size())
            _read_count_FR_map.resize(2*(region_idx+1));

        _read_count_FR_map[region_idx] = nread_FR - nread_ROI;
    }

    void accumulate_reads_in_region(ReadCountsByLib& acc, size_t begin, size_t end) {
        for(size_t i = begin; i < std::min(end, _read_count_ROI_map.size()); i++){
            acc += _read_count_ROI_map[i];

            // flanking region doesn't contain the first node
            if(i > begin && i < _read_count_FR_map.size())
                acc += _read_count_FR_map[i];
        }
    }

    uint32_t region_lib_read_count(size_t region_idx, std::string const& lib) const {
        return _region_lib_counts(region_idx, lib, _read_count_ROI_map);
    }

    void swap_reads_in_region(size_t region_idx, ReadVector& reads) {
        _regions[region_idx]->swap_reads(reads);
    }

    ReadVector const& reads_in_region(size_t region_idx) const {
        return _regions[region_idx]->reads();
    }

    bool region_exists(size_t region_idx) const {
        return (region_idx < _regions.size()) && _regions[region_idx];
    }

    void clear_region(size_t region_idx) {
        delete _regions[region_idx];
        _regions[region_idx] = 0;
    }

    size_t num_regions() const {
        return _regions.size();
    }

    size_t add_region(BasicRegion* r) {
        _regions.push_back(r);
        return _regions.size() - 1;
    }

    BasicRegion const& get_region_data(size_t region_idx) const {
        return *_regions[region_idx];
    }

    void set_max_read_window_size(int val) {
        _max_read_window_size = val;
    }

    void process_breakpoint(bam_header_t const* bam_header);
    void process_final_region(bam_header_t const* bam_header);

    void run();

private:
    uint32_t _region_lib_counts(size_t region_idx, std::string const& lib, RoiReadCounts const& x) const {
        if (region_idx >= x.size())
            return 0;
        RoiReadCounts::value_type::const_iterator found = x[region_idx].find(lib);
        if (found != x[region_idx].end())
            return found->second;
        return 0;
    }

private:
    Options const& _opts;
    BamConfig const& _cfg;
    IBamReader& _merged_reader;
    int _max_read_window_size;

    RoiReadCounts _read_count_ROI_map;
    RoiReadCounts _read_count_FR_map;
    RegionData _regions;

    bool _normal_switch;
    int _nnormal_reads;
    int _ntotal_nucleotides;
    int _max_readlen;
    int _buffer_size;

public:
    std::map<std::string, std::vector<int> > _read_regions;
    int begins; // global (chr)
    int beginc; // global
    int lasts; // global (chr, should be int in samtools)
    int lastc; // global

};

// choose the predominant type of read in a region
inline
breakdancer::pair_orientation_flag choose_sv_flag(const int num_readpairs, const std::map<breakdancer::pair_orientation_flag, int> reads_per_type) {
    using namespace std;
    breakdancer::pair_orientation_flag flag = breakdancer::NA;
    float max_type_pct = 0.0;
    for(map<breakdancer::pair_orientation_flag, int>::const_iterator type_iter = reads_per_type.begin(); type_iter != reads_per_type.end(); ++type_iter){
        breakdancer::pair_orientation_flag current_flag = (*type_iter).first;
        float type_pct = static_cast<float>((*type_iter).second) / static_cast<float>(num_readpairs);
        if(max_type_pct < type_pct) {
            max_type_pct = type_pct;
            flag = current_flag;
        }
    }
    return flag;
}


