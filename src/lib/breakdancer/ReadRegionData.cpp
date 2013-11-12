#include "ReadRegionData.hpp"
#include "common/Timer.hpp"

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
    out << "Number of tracked reads: " << _read_regions.size() << "\n";
    out << "Active region summary:\n";
    size_t n_total(0);
    size_t n_active(0);
    size_t n_deleted(0);
    out <<
        "region_id"
        "\tregion_tid"
        "\tregion_start"
        "\tregion_end"
        "\tnreads"
        "\tn_unpaired_reads"
        "\ttimes_processed"
        "\ttimes_collapsed\n"
        ;

    for (size_t i = 0; i < _regions.size(); ++i) {
        ++n_total;
        if (_regions[i]) {
            ++n_active;

            BasicRegion const& r = *_regions[i];
            out << i
                << "\t" << r.chr
                << "\t" << r.start
                << "\t" << r.end
                << "\t" << r.reads().size()
                << "\t" << DEBUG_unpaired_reads(i)
                << "\t" << r.times_accessed
                << "\t" << r.times_collapsed
                << "\n"
                ;

            size_t num_reads = std::min(size_t(10ull), r.reads().size());
            if (num_reads) {
                out << "\tfirst " << num_reads << " read names:\n";
                for (size_t j = 0; j < num_reads; ++j) {
                    out << "\t\t" << r.reads()[j].query_name() << "\n";
                }
            }
        }
        else {
            ++n_deleted;
        }
    }
    out << "Total regions: " << n_total << "\n";
    out << "Active regions: " << n_active << "\n";
    out << "Deleted regions: " << n_deleted << "\n";
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

    int non_ctx_reads(0);

    // This adds the region id to an array of region ids
    for(ReadVector::const_iterator iter = reads.begin(); iter != reads.end(); ++iter) {
        if (iter->bdflag() != breakdancer::ARP_CTX)
            ++non_ctx_reads;

        if (iter->ori() == FWD)
            ++_regions.back()->fwd_read_count;
        else
            ++_regions.back()->rev_read_count;

        std::vector<int>& regions = _read_regions[iter->query_name()];
        regions.push_back(region_idx);
        if (regions.size() == 2) {
            _persistent_graph.increment_edge_weight(regions[0], regions[1]);
        }
    }

    // we're essentially destroying reads_in_current_region here by swapping it with whatever
    //reads this region had (probably none) this is ok because it is just about to be cleared anyway.
    int valid_reads = _opts.chr.empty() ? reads.size() : non_ctx_reads;
    if (valid_reads >= _opts.min_read_pair) {
        swap_reads_in_region(region_idx, reads);
    }

    return region_idx;
}

bool ReadRegionData::is_region_final(size_t region_idx) const {
    if (!region_exists(region_idx) || region_idx == last_region_idx())
        return false;

    ReadVector const& v = _reads_in_region(region_idx);
    for (ReadVector::const_iterator i = v.begin(); i != v.end(); ++i) {
        if (!_opts.chr.empty() && i->bdflag() == breakdancer::ARP_CTX)
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
        ++_regions[last_region_idx()]->times_collapsed;
    }

    // remove any reads that are linking the last region with this new, merged in region
    for(ReadVector::const_iterator iter = reads.begin(); iter != reads.end(); ++iter) {
        std::string const& read_name = iter->query_name();
/*
        ReadsToRegionsMap::const_iterator found = _read_regions.find(read_name);
        if (found != _read_regions.end()) {
            cerr << "Perhaps we shouldn't drop this read (collapsing to " << last_region_idx() << "):\n";
            cerr << "\t" << read_name << " in " << found->second.size() << " regions:";
            for (size_t i = 0; i < found->second.size(); ++i) {
                cerr << " " << found->second[i];
            }
            cerr << "\n";
        }
*/
        _read_regions.erase(read_name);
    }
}

ReadRegionData::read_iter_range ReadRegionData::region_reads_range(size_t region_idx) const {
    return _regions[region_idx]->reads_range(
        boost::bind(&ReadRegionData::read_exists, this, _1)
        );
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
