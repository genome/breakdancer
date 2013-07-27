#pragma once

#include "BasicRegion.hpp"
#include "ReadCountsByLib.hpp"
#include "common/Graph.hpp"
#include "common/Options.hpp"
#include "io/Read.hpp"

#include <boost/function.hpp>
#include <boost/range/algorithm/remove_copy_if.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <cassert>
#include <iterator>
#include <map>
#include <ostream>
#include <string>
#include <vector>

class ReadRegionData {
public:
    typedef breakdancer::Read ReadType;
    typedef BasicRegion::ReadVector ReadVector;
    typedef BasicRegion::const_read_iterator const_read_iterator;
    typedef BasicRegion::iterator_range read_iter_range;
    typedef std::vector<BasicRegion*> RegionData;
    typedef std::vector<ReadCountsByLib> RoiReadCounts;
    typedef boost::unordered_map<std::string, std::vector<int> > ReadsToRegionsMap;
    typedef UndirectedWeightedGraph<int, int> Graph; // tmpl params=vertex type, weight type.

public:
    ReadRegionData(Options const& opts)
        : _opts(opts)
    {
    }

    ~ReadRegionData();

    void summary(std::ostream& s) const;

//    Graph region_graph() const;

    void accumulate_reads_between_regions(ReadCountsByLib& acc, size_t begin, size_t end) const;
    uint32_t region_lib_read_count(size_t region_idx, std::string const& lib) const;
    void swap_reads_in_region(size_t region_idx, ReadVector& reads);

    void remove_reads_in_region_if(size_t region_idx, boost::function<bool(ReadType const&)> pred) {
        BasicRegion::ReadVector filtered;
        boost::remove_copy_if(region_reads_range(region_idx),
            std::back_inserter(filtered), pred);
        swap_reads_in_region(region_idx, filtered);
    }

    size_t num_reads_in_region(size_t region_idx) const;
    bool region_exists(size_t region_idx) const;

    size_t add_region(int start_tid, int start_pos, int end_pos, int normal_reads,
            ReadVector& reads);
    bool is_region_final(size_t region_idx) const;

    int sum_of_region_sizes(std::vector<int> const& region_ids) const;

    void clear_region(size_t region_idx);
    size_t num_regions() const;
    size_t last_region_idx() const;
    BasicRegion const& region(size_t region_idx) const;

    void incr_normal_read_count(ReadCountsByLib::LibId const& key);
    void clear_region_accumulator();
    void clear_flanking_region_accumulator();
    void collapse_accumulated_data_into_last_region(ReadVector const& reads);
    ReadsToRegionsMap const& read_regions() const;
    void erase_read(std::string const& read_name);
    bool read_exists(ReadType const& read) const;

    read_iter_range region_reads_range(size_t region_idx) const;

    Graph& persistent_graph() {
        return _persistent_graph;
    }

    void incr_region_access_counter(size_t region_idx) {
        if (!region_exists(region_idx))
            return;
        ++_regions[region_idx]->times_accessed;
    }

private:
    void _add_current_read_counts_to_region(size_t region_idx);
    void _add_per_lib_read_counts_to_last_region(ReadCountsByLib const& counts);
    ReadVector const& _reads_in_region(size_t region_idx) const;

    size_t DEBUG_unpaired_reads(size_t region_idx) const {
        if (!region_exists(region_idx))
            return 0u;

        size_t rv = 0;
        ReadVector const& v = _reads_in_region(region_idx);
        for (ReadVector::const_iterator i = v.begin(); i != v.end(); ++i) {
            ReadsToRegionsMap::const_iterator found = _read_regions.find(i->query_name());
            if (found == _read_regions.end() || found->second.size() != 2)
                ++rv;
        }
        return rv;
    }

private:
    Options const& _opts;
    RoiReadCounts _read_count_ROI_map;
    RoiReadCounts _read_count_FR_map;
    RegionData _regions;

    ReadCountsByLib nread_ROI;
    ReadCountsByLib nread_FR;

    ReadsToRegionsMap _read_regions;

    Graph _persistent_graph;
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
size_t ReadRegionData::last_region_idx() const {
    return num_regions() - 1;
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
