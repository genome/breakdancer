#pragma once

#include "BasicRegion.hpp"
#include "Read.hpp"
#include "ReadCountsByLib.hpp"

#include <boost/unordered_map.hpp>

#include <cassert>
#include <map>
#include <string>
#include <vector>

class ReadRegionData {
public:
    typedef breakdancer::Read ReadType;
    typedef BasicRegion::ReadVector ReadVector;
    typedef std::vector<BasicRegion*> RegionData;
    typedef std::vector<ReadCountsByLib> RoiReadCounts;
    typedef boost::unordered_map<std::string, std::vector<int> > ReadsToRegionsMap;
    typedef BasicRegion::const_read_iterator const_read_iterator;

    ~ReadRegionData();

    void accumulate_reads_between_regions(ReadCountsByLib& acc, size_t begin, size_t end) const;
    uint32_t region_lib_read_count(size_t region_idx, std::string const& lib) const;
    void swap_reads_in_region(size_t region_idx, ReadVector& reads);

    size_t num_reads_in_region(size_t region_idx) const;
    bool region_exists(size_t region_idx) const;

    size_t add_region(int start_tid, int start_pos, int end_pos, int normal_reads,
            ReadVector& reads);
    int sum_of_region_sizes(std::vector<int> const& region_ids) const;

    void clear_region(size_t region_idx);
    size_t num_regions() const;
    BasicRegion const& region(size_t region_idx) const;

    void incr_normal_read_count(ReadCountsByLib::LibId const& key);
    void clear_region_accumulator();
    void clear_flanking_region_accumulator();
    void collapse_accumulated_data_into_last_region(ReadVector const& reads);
    ReadsToRegionsMap const& read_regions() const;
    void erase_read(std::string const& read_name);
    bool read_exists(ReadType const& read) const;

    const_read_iterator region_read_begin(size_t region_idx) const;
    const_read_iterator region_read_end(size_t region_idx) const;

private:
    void _add_current_read_counts_to_region(size_t region_idx);
    void _add_per_lib_read_counts_to_last_region(ReadCountsByLib const& counts);
    ReadVector const& _reads_in_region(size_t region_idx) const;


private:
    RoiReadCounts _read_count_ROI_map;
    RoiReadCounts _read_count_FR_map;
    RegionData _regions;

    ReadCountsByLib nread_ROI;
    ReadCountsByLib nread_FR;

    ReadsToRegionsMap _read_regions;
};

inline
void ReadRegionData::swap_reads_in_region(size_t region_idx, ReadVector& reads) {
    assert(_regions[region_idx] != 0);
    _regions[region_idx]->swap_reads(reads);
}

inline
size_t ReadRegionData::num_reads_in_region(size_t region_idx) const {
    return _reads_in_region(region_idx).size();
}

inline
ReadRegionData::ReadVector const& ReadRegionData::_reads_in_region(size_t region_idx) const {
    return _regions[region_idx]->reads();
}

inline
bool ReadRegionData::region_exists(size_t region_idx) const {
    return (region_idx < _regions.size()) && _regions[region_idx];
}

inline
size_t ReadRegionData::num_regions() const {
    return _regions.size();
}

inline
BasicRegion const& ReadRegionData::region(size_t region_idx) const {
    assert(_regions[region_idx] != 0);
    return *_regions[region_idx];
}

inline
void ReadRegionData::incr_normal_read_count(ReadCountsByLib::LibId const& key) {
    ++nread_ROI[key];
    ++nread_FR[key];
}

inline
void ReadRegionData::clear_region_accumulator() {
    nread_ROI.clear();
}

inline
void ReadRegionData::clear_flanking_region_accumulator() {
    nread_FR.clear();
}

inline
ReadRegionData::ReadsToRegionsMap const&
ReadRegionData::read_regions() const {
    return _read_regions;
}

inline
void ReadRegionData::erase_read(std::string const& read_name) {
    _read_regions.erase(read_name);
}

inline
bool ReadRegionData::read_exists(ReadType const& read) const {
    return _read_regions.count(read.query_name()) > 0;
}
