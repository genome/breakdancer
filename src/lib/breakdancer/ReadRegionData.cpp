#include "ReadRegionData.hpp"
#include "Timer.hpp"

#include <boost/format.hpp>
#include <boost/bind.hpp>
#include <iostream>

using std::cerr;

ReadRegionData::~ReadRegionData() {
    if (getenv("BD_DUMP_REGION_SUMMARY"))
        summary(cerr);

    for (size_t i = 0; i < _regions.size(); ++i) {
        delete _regions[i];
    }
}

void ReadRegionData::summary(std::ostream& out) const {
    out << "Region summary:\n";
    size_t n_total(0);
    size_t n_active(0);
    size_t n_deleted(0);
    for (size_t i = 0; i < _regions.size(); ++i) {
        bool active = _regions[i] != 0;
        ++n_total;
        if (active)
            ++n_active;
        else
            ++n_deleted;

        out << i << "\t" << (active ? "ACTIVE" : "DELETED")
            << "\t" << (active ? _regions[i]->size() : 0)
            << "\t" << __DEBUG_unpaired_reads(i)
            << "\n";
    }
    out << "Total regions: " << n_total << "\n";
    out << "Active regions: " << n_active << "\n";
    out << "Deleted regions: " << n_deleted << "\n";
}

ReadRegionData::Graph ReadRegionData::region_graph() const {
    Graph graph;
    typedef ReadsToRegionsMap::const_iterator IterType;

    using namespace boost::chrono;
    using std::cerr;
    using boost::format;

    for(IterType i = read_regions().begin(); i != read_regions().end(); i++) {
        // test
        std::vector<int> const& p = i->second;
        assert(p.size() < 3);
        if(p.size() != 2) // skip singleton read (non read pairs)
            continue;

        int const& r1 = p[0];
        int const& r2 = p[1];

        if (!region_exists(r1) || !region_exists(r2))
            continue;

        graph.increment_edge_weight(r1, r2);

    }

    return graph;
}

void ReadRegionData::accumulate_reads_between_regions(ReadCountsByLib& acc, size_t begin, size_t end) const {
    for(size_t i = begin; i < std::min(end, _read_count_ROI_map.size()); i++){
        acc += _read_count_ROI_map[i];

        // flanking region doesn't contain the first node
        if(i > begin && i < _read_count_FR_map.size())
            acc += _read_count_FR_map[i];
    }
}

uint32_t ReadRegionData::region_lib_read_count(size_t region_idx, std::string const& lib) const {
    if (region_idx >= _read_count_ROI_map.size())
        return 0;
    RoiReadCounts::value_type::const_iterator found = _read_count_ROI_map[region_idx].find(lib);
    if (found != _read_count_ROI_map[region_idx].end())
        return found->second;
    return 0;
}

size_t ReadRegionData::add_region(int start_tid, int start_pos, int end_pos, int normal_reads,
        ReadVector& reads)
{
    size_t region_idx = _regions.size();
    _regions.push_back(new BasicRegion(region_idx, start_tid, start_pos, end_pos, normal_reads));
    _add_current_read_counts_to_region(region_idx);

    // This adds the region id to an array of region ids
    for(ReadVector::const_iterator iter = reads.begin(); iter != reads.end(); ++iter) {
        std::vector<int>& regions = _read_regions[iter->query_name()];
        regions.push_back(region_idx);
    }

    // we're essentially destroying reads_in_current_region here by swapping it with whatever
    //reads this region had (probably none) this is ok because it is just about to be cleared anyway.
    swap_reads_in_region(region_idx, reads);

    return region_idx;
}

bool ReadRegionData::is_region_final(size_t region_idx) const {
    if (!region_exists(region_idx))
        return false;

    ReadVector const& v = _reads_in_region(region_idx);
    for (ReadVector::const_iterator i = v.begin(); i != v.end(); ++i) {
        if (_opts.chr != "0" && i->bdflag() == breakdancer::ARP_CTX)
            continue;

        ReadsToRegionsMap::const_iterator found = _read_regions.find(i->query_name());
        if (found == _read_regions.end() || found->second.size() != 2)
            return false;
    }

    return true;
}

int ReadRegionData::sum_of_region_sizes(std::vector<int> const& region_ids) const {
    typedef std::vector<int>::const_iterator IterType;
    int size(0);
    for (IterType i = region_ids.begin(); i != region_ids.end(); ++i)
        size += region(*i).size();
    return size;
}

void ReadRegionData::clear_region(size_t region_idx) {
    if (!region_exists(region_idx))
        return;

    BasicRegion::ReadVector const& reads = _reads_in_region(region_idx);
    for(ReadVector::const_iterator i = reads.begin(); i != reads.end(); ++i) {
        ReadsToRegionsMap::iterator found = _read_regions.find(i->query_name());
        if (found != _read_regions.end()) {
            std::vector<int> new_regions;
            std::vector<int>& old_regions(found->second);
            for (size_t j = 0; j != old_regions.size(); ++j) {
                if (old_regions[j] != int(region_idx))
                    new_regions.push_back(old_regions[j]);
            }
            if (!new_regions.empty())
                found->second.swap(new_regions);
            else
                _read_regions.erase(found);
        }
    }

    delete _regions[region_idx];
    _regions[region_idx] = 0;
}

void ReadRegionData::collapse_accumulated_data_into_last_region(ReadVector const& reads) {
    if(num_regions() > 0) {
        _add_per_lib_read_counts_to_last_region(nread_FR);
    }

    // remove any reads that are linking the last region with this new, merged in region
    for(ReadVector::const_iterator it_reg_seq = reads.begin(); it_reg_seq != reads.end(); ++it_reg_seq) {
        _read_regions.erase(it_reg_seq->query_name());
    }
}

ReadRegionData::read_iter_range ReadRegionData::region_reads_range(size_t region_idx) const {
    return _regions[region_idx]->reads_range(boost::bind(&ReadRegionData::read_exists, this, _1));
}

void ReadRegionData::_add_current_read_counts_to_region(size_t region_idx) {
    if (region_idx >= _read_count_ROI_map.size())
        _read_count_ROI_map.resize(2*(region_idx+1));

    _read_count_ROI_map[region_idx] = nread_ROI;

    if (region_idx >= _read_count_FR_map.size())
        _read_count_FR_map.resize(2*(region_idx+1));

    _read_count_FR_map[region_idx] = nread_FR - nread_ROI;
}

void ReadRegionData::_add_per_lib_read_counts_to_last_region(ReadCountsByLib const& counts) {
    assert(num_regions() > 0);

    size_t region_idx = last_region_idx();
    if (region_idx >= _read_count_ROI_map.size())
        _read_count_ROI_map.resize(2*(region_idx+1));

    _read_count_ROI_map[region_idx] +=  counts;
}
