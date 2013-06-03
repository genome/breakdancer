#pragma once

#include "BasicRegion.hpp"
#include "Read.hpp"
#include "ReadCountsByLib.hpp"

#include <map>
#include <vector>
#include <string>

class ReadRegionData {
public:
    typedef breakdancer::Read ReadType;
    typedef BasicRegion::ReadVector ReadVector;
    typedef std::vector<BasicRegion*> RegionData;
    typedef std::vector<ReadCountsByLib> RoiReadCounts;
    typedef std::map<std::string, std::vector<int> > ReadsToRegionsMap;

    ~ReadRegionData();

    void accumulate_reads_between_regions(ReadCountsByLib& acc, size_t begin, size_t end) const {
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
        BasicRegion::ReadVector const& reads = reads_in_region(region_idx);

        for(ReadVector::const_iterator i = reads.begin(); i != reads.end(); ++i)
            erase_read(i->query_name());

        delete _regions[region_idx];
        _regions[region_idx] = 0;
    }

    size_t num_regions() const {
        return _regions.size();
    }

    size_t add_region(int start_tid, int start_pos, int end_pos, int normal_reads,
            ReadVector& reads)
    {
        _regions.push_back(new BasicRegion(start_tid, start_pos, end_pos, normal_reads));
        size_t region_idx = _regions.size() - 1;
        _add_current_read_counts_to_region(region_idx);

        // This adds the region id to an array of region ids
        for(ReadVector::const_iterator iter = reads.begin(); iter != reads.end(); ++iter) {
            _read_regions[iter->query_name()].push_back(region_idx);
        }

        // we're essentially destroying reads_in_current_region here by swapping it with whatever
        //reads this region had (probably none) this is ok because it is just about to be cleared anyway.
        swap_reads_in_region(region_idx, reads);

        return region_idx;
    }

    BasicRegion const& region(size_t region_idx) const {
        return *_regions[region_idx];
    }

    void incr_normal_read_count(ReadCountsByLib::LibId const& key) {
        ++nread_ROI[key];
        ++nread_FR[key];
    }

    void clear_region_accumulator() {
        nread_ROI.clear();
    }

    void clear_flanking_region_accumulator() {
        nread_FR.clear();
    }

    void collapse_accumulated_data_into_last_region(ReadVector const& reads) {
        if(num_regions() > 0) {
            _add_per_lib_read_counts_to_last_region(nread_FR);
        }

        // remove any reads that are linking the last region with this new, merged in region
        for(ReadVector::const_iterator it_reg_seq = reads.begin(); it_reg_seq != reads.end(); ++it_reg_seq) {
            _read_regions.erase(it_reg_seq->query_name());
        }
    }

    ReadsToRegionsMap const& read_regions() const {
        return _read_regions;
    }

    void erase_read(std::string const& read_name) {
        _read_regions.erase(read_name);
    }

    bool read_exists(std::string const& read_name) {
        return _read_regions.count(read_name) > 0;
    }

private:
    void _add_current_read_counts_to_region(size_t region_idx) {
        if (region_idx >= _read_count_ROI_map.size())
            _read_count_ROI_map.resize(2*(region_idx+1));

        _read_count_ROI_map[region_idx] = nread_ROI;

        if (region_idx >= _read_count_FR_map.size())
            _read_count_FR_map.resize(2*(region_idx+1));

        _read_count_FR_map[region_idx] = nread_FR - nread_ROI;
    }

    uint32_t _region_lib_counts(size_t region_idx, std::string const& lib, RoiReadCounts const& x) const {
        if (region_idx >= x.size())
            return 0;
        RoiReadCounts::value_type::const_iterator found = x[region_idx].find(lib);
        if (found != x[region_idx].end())
            return found->second;
        return 0;
    }

private:
    void _add_per_lib_read_counts_to_last_region(ReadCountsByLib const& counts) {
        assert(num_regions() > 0);

        size_t region_idx = num_regions()-1;
        if (region_idx >= _read_count_ROI_map.size())
            _read_count_ROI_map.resize(2*(region_idx+1));

        _read_count_ROI_map[region_idx] +=  counts;
    }


private:
    RoiReadCounts _read_count_ROI_map;
    RoiReadCounts _read_count_FR_map;
    RegionData _regions;

    ReadCountsByLib nread_ROI;
    ReadCountsByLib nread_FR;

    ReadsToRegionsMap _read_regions;
};
