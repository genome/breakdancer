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

    BreakDancer()
        : begins(-1)
        , beginc(-1)
        , lasts(-1)
        , lastc(-1)
    {
    }

    ~BreakDancer() {
        for (size_t i = 0; i < region_data.size(); ++i) {
            delete region_data[i];
        }
    }

    int begins; // global (chr)
    int beginc; // global
    int lasts; // global (chr, should be int in samtools)
    int lastc; // global
    std::map<std::string, uint32_t> nread_ROI; // global
    std::map<int, std::map<std::string, uint32_t> > read_count_ROI_map; // global
    std::map<std::string, uint32_t> nread_FR;    // global
    std::map<int, std::map<std::string, uint32_t> > read_count_FR_map; // global


    // Ultimately, we'll want to be pushing onto this vector and
    // returning new region indices as appropriate.
    void swap_reads_in_region(size_t region_idx, ReadVector& reads) {
        region_data[region_idx]->swap_reads(reads);
    }

    ReadVector const& reads_in_region(size_t region_idx) const {
        return region_data[region_idx]->reads();
    }

    bool region_exists(size_t region_idx) const {
        return (region_idx < region_data.size()) && region_data[region_idx];
    }

    void clear_region(size_t region_idx) {
        delete region_data[region_idx];
        region_data[region_idx] = 0;
    }

    void add_region(size_t region_idx, BasicRegion* r) {
        // FIXME: ultimately, we'll want to push_back and return the region_idx
        if (region_idx >= region_data.size())
            region_data.resize((region_idx+1)*2);
        else if (region_data[region_idx])
            delete region_data[region_idx];

        region_data[region_idx] = r;
    }

    BasicRegion const& get_region_data(size_t region_idx) const {
        using boost::format;
        if (!region_exists(region_idx)) {
            throw std::runtime_error(str(format("Unknown region %1%")
                % region_idx));
        }
        return *region_data[region_idx];
    }

private:
    RegionData region_data;
};


