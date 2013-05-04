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
    typedef std::map<int, std::map<std::string, uint32_t> > RoiReadCounts;

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
    RoiReadCounts read_count_ROI_map; // global
    RoiReadCounts read_count_FR_map; // global

    void increment_region_lib_read_count(int region_idx, std::string const& lib, uint32_t nreads) {
        read_count_ROI_map[region_idx][lib] += nreads;
    }

    void set_region_lib_FR_count(int region_idx, std::string const& lib, uint32_t nreads) {
        read_count_FR_map[region_idx][lib] = nreads;
    }

    RoiReadCounts::mapped_type const* region_read_counts_by_library(int region_idx) {
        RoiReadCounts::const_iterator rfound = read_count_ROI_map.find(region_idx);
        if (rfound == read_count_ROI_map.end())
            return 0;
        return &rfound->second;
    }

    RoiReadCounts::mapped_type const* region_FR_counts_by_library(int region_idx) {
        RoiReadCounts::const_iterator rfound = read_count_FR_map.find(region_idx);
        if (rfound == read_count_FR_map.end())
            return 0;
        return &rfound->second;
    }

    uint32_t region_lib_read_count(int region_idx, std::string const& lib) const {
        RoiReadCounts::const_iterator rfound = read_count_ROI_map.find(region_idx);
        if (rfound != read_count_ROI_map.end()) {
            std::map<std::string, uint32_t>::const_iterator lfound = rfound->second.find(lib);
            if (lfound != rfound->second.end())
                return lfound->second;
        }

        return 0;
    }

    uint32_t region_lib_FR_count(int region_idx, std::string const& lib) const {
        RoiReadCounts::const_iterator rfound = read_count_FR_map.find(region_idx);
        if (rfound != read_count_FR_map.end()) {
            std::map<std::string, uint32_t>::const_iterator lfound = rfound->second.find(lib);
            if (lfound != rfound->second.end())
                return lfound->second;
        }

        return 0;
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

private:
    RegionData _regions;
};


