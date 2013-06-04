#include "SvBuilder.hpp"
#include "LibraryInfo.hpp"

#include <algorithm>
#include <stdexcept>

using namespace std;
namespace bd = breakdancer;

namespace {
    const static bd::PerFlagArray<int>::type ZEROS = {{0}};
}

SvBuilder::SvBuilder()
    : current_region(-1)
    , num_regions(0)
    , num_pairs(0)
    , flag_counts(ZEROS)
    , type_orient_counts(2)
{
    for (int i = 0; i < 2; ++i) {
        type_orient_counts[i][FWD] = 0;
        type_orient_counts[i][REV] = 0;
    }
}

// choose the predominant type of read in a region
bd::pair_orientation_flag SvBuilder::choose_sv_flag() {
    bd::pair_orientation_flag flag = bd::NA;
    int const* max_ptr = max_element(flag_counts.begin(), flag_counts.end());
    if (max_ptr != flag_counts.end() && *max_ptr > 0) {
        flag = bd::pair_orientation_flag(max_ptr - flag_counts.begin());
    }
    return flag;
}

void SvBuilder::observe_read(Read const& read, int region) {
    if (region != current_region) {
        if (++num_regions > 2) {
            throw runtime_error("Attempted to build sv with more than 2 regions");
        }
        current_region = region;
    }
    int region_idx = num_regions - 1;

    ++type_orient_counts[region_idx][read.ori()];

    typedef ObservedReads::iterator IterType;
    pair<IterType, bool> inserted = observed_reads.insert(make_pair(read.query_name(), read));
    if(inserted.second) {
        // This is the first time we have seen a read with this name.
        observed_reads[read.query_name()] = read;
    }
    else {
        // We just found an existing read's mate. Good for him/her.
        bd::pair_orientation_flag bdflag = read.bdflag();
        string const& libname = read.lib_info().name;
        ++flag_counts[bdflag];
        ++type_library_readcount[bdflag][libname];
        type_library_meanspan[bdflag][libname] += read.abs_isize();

        ++num_pairs;
        reads_to_free.push_back(read.query_name());
        support_reads.push_back(read);
        support_reads.push_back(inserted.first->second);
        observed_reads.erase(inserted.first);
    }
}
