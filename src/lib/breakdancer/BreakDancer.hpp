#pragma once

#include "BasicRegion.hpp"

#include <cstddef>
#include <map>
#include <stdint.h>
#include <string>
#include <vector>

class BreakDancer {
public:
    typedef BasicRegion::ReadVector ReadVector;
    typedef std::vector<BasicRegion*> RegionData;
    typedef std::vector<std::map<std::string, uint32_t> > RoiReadCounts;

    BreakDancer()
        : begins(-1)
        , beginc(-1)
        , lasts(-1)
        , lastc(-1)
    {
    }

    ~BreakDancer() {
        for (size_t i = 0; i < _regions.size(); ++i) {
            delete _regions[i];
        }
    }

    int begins; // global (chr)
    int beginc; // global
    int lasts; // global (chr, should be int in samtools)
    int lastc; // global
    std::map<std::string, uint32_t> nread_ROI; // global
    std::map<std::string, uint32_t> nread_FR;    // global
    std::map<std::string, float> read_density;
    ReadVector reads_in_current_region;

    void increment_region_lib_read_count(size_t region_idx, std::string const& lib, uint32_t nreads) {
        if (region_idx >= _read_count_ROI_map.size())
            _read_count_ROI_map.resize(2*(region_idx+1));
        _read_count_ROI_map[region_idx][lib] += nreads;
    }

    void set_region_lib_FR_count(size_t region_idx, std::string const& lib, uint32_t nreads) {
        if (region_idx >= _read_count_FR_map.size())
            _read_count_FR_map.resize(2*(region_idx+1));
        _read_count_FR_map[region_idx][lib] = nreads;
    }

    RoiReadCounts::value_type const* region_read_counts_by_library(size_t region_idx) {
        if (region_idx >= _read_count_ROI_map.size())
            return 0;
        return &_read_count_ROI_map[region_idx];
    }

    RoiReadCounts::value_type const* region_FR_counts_by_library(size_t region_idx) {
        if (region_idx >= _read_count_FR_map.size())
            return 0;
        return &_read_count_FR_map[region_idx];
    }

    uint32_t region_lib_read_count(size_t region_idx, std::string const& lib) const {
        return _region_lib_counts(region_idx, lib, _read_count_ROI_map);
    }

    uint32_t region_lib_FR_count(size_t region_idx, std::string const& lib) const {
        return _region_lib_counts(region_idx, lib, _read_count_FR_map);
    }


    // Ultimately, we'll want to be pushing onto this vector and
    // returning new region indices as appropriate.
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

    RoiReadCounts const& read_count_ROI_map() const {
        return _read_count_ROI_map;
    }

    RoiReadCounts const& read_count_FR_map() const {
        return _read_count_FR_map;
    }

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
    RoiReadCounts _read_count_ROI_map;
    RoiReadCounts _read_count_FR_map;
    RegionData _regions;
};


