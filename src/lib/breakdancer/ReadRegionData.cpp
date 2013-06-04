#include "ReadRegionData.hpp"

#include <boost/bind.hpp>

ReadRegionData::~ReadRegionData() {
    for (size_t i = 0; i < _regions.size(); ++i) {
        delete _regions[i];
    }
}

ReadRegionData::Graph ReadRegionData::region_graph() const {
    Graph graph;
    typedef ReadsToRegionsMap::const_iterator IterType;
    for(IterType i = read_regions().begin(); i != read_regions().end(); i++){
        // test
        std::vector<int> const& p = i->second;
        assert(p.size() < 3);
        if(p.size() != 2) // skip singleton read (non read pairs)
            continue;

        int const& r1 = p[0];
        int const& r2 = p[1];

        //track the number of links between two nodes
        //
        // This doesn't make a lot of sense to me. When r1 == r2 and r1 is not
        // in the map, both are set to one. If r1 is in the map, then we increment
        // twice. We should either double count or not. Doing a mixture of both is
        // silly. -ta
        if(graph.find(r1) != graph.end() && graph[r1].find(r2) != graph[r1].end()){
            ++graph[r1][r2];
            ++graph[r2][r1];
        }
        else{
            graph[r1][r2] = 1;
            graph[r2][r1] = 1;
        }
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
        _read_regions[iter->query_name()].push_back(region_idx);
    }

    // we're essentially destroying reads_in_current_region here by swapping it with whatever
    //reads this region had (probably none) this is ok because it is just about to be cleared anyway.
    swap_reads_in_region(region_idx, reads);

    return region_idx;
}

int ReadRegionData::sum_of_region_sizes(std::vector<int> const& region_ids) const {
    typedef std::vector<int>::const_iterator IterType;
    int size(0);
    for (IterType i = region_ids.begin(); i != region_ids.end(); ++i)
        size += region(*i).size();
    return size;
}

void ReadRegionData::clear_region(size_t region_idx) {
    BasicRegion::ReadVector const& reads = _reads_in_region(region_idx);

    for(ReadVector::const_iterator i = reads.begin(); i != reads.end(); ++i)
        erase_read(i->query_name());

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

    size_t region_idx = num_regions()-1;
    if (region_idx >= _read_count_ROI_map.size())
        _read_count_ROI_map.resize(2*(region_idx+1));

    _read_count_ROI_map[region_idx] +=  counts;
}
