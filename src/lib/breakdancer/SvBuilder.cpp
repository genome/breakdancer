#include "SvBuilder.hpp"
#include "ReadCountsByLib.hpp"
#include "config/LibraryInfo.hpp"

#include <boost/bind.hpp>
#include <boost/range/algorithm/for_each.hpp>

#include <algorithm>
#include <cassert>
#include <stdexcept>

using namespace std;
namespace bd = breakdancer;

namespace {
    const static bd::PerFlagArray<int>::type ZEROS = {{0}};
}

SvBuilder::SvBuilder(Options const& opts, int n, BasicRegion const* regions[2],
        ReadsRange read_ranges[2], int max_readlen)
    : current_region(-1)
    , num_regions(n)
    , num_pairs(0)
    , diffspan(0)
    , flag_counts(ZEROS)
    , chr(_init_zero())
    , pos(_init_zero())
    , fwd_read_count(_init_zero())
    , rev_read_count(_init_zero())
    , allele_frequency(0.0f)
    , _opts(opts)
{
    assert(n == 1 || n == 2);

    for (int i = 0; i < n; ++i) {
        boost::range::for_each(read_ranges[i],
            boost::bind(&SvBuilder::_observe_read, this, _1, i));

        fwd_read_count[i] = regions[i]->fwd_read_count;
        rev_read_count[i] = regions[i]->rev_read_count;
    }

    flag = choose_sv_flag();

    chr[0] = regions[0]->chr;
    pos[0] = regions[0]->start;
    pos[1] = regions[0]->end;

    // the sv may be contained in a single region, in which case
    // regions[1] will be null.
    if (n == 2) {
        if(flag == bd::ARP_RF) {
            pos[1] = regions[1]->end + max_readlen - 5;
        }
        else if(flag == bd::ARP_FF) {
            pos[0] = pos[1];
            pos[1] = regions[1]->end + max_readlen - 5;
        }
        else if(flag == bd::ARP_RR) {
            pos[1] = regions[1]->start;
        }
        else {
            pos[0] = pos[1];
            pos[1] = regions[1]->start;
        }
        chr[1] = regions[1]->chr;
    }
    else {
        fwd_read_count[1] = fwd_read_count[0];
        rev_read_count[1] = rev_read_count[0];
        chr[1] = regions[0]->chr;
        pos[1] = regions[0]->end;
    }
}

void SvBuilder::compute_copy_number(ReadCountsByLib const& counts,
    std::map<std::string, float> const& read_density)
{
    typedef ReadCountsByLib::const_iterator IterType;
    float copy_number_sum = 0.0f;
    for(IterType iter = counts.begin(); iter != counts.end(); ++iter) {
        string const& lib = iter->first;
        copy_number[lib] = iter->second/(read_density.at(lib) * float(pos[1] - pos[0]))*2.0f;
        copy_number_sum += copy_number[lib];
    }
    copy_number_sum /= 2.0f * counts.size();
    allele_frequency = 1 - copy_number_sum;
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

void SvBuilder::_observe_read(Read const& read, int region_idx) {
    typedef ObservedReads::iterator IterType;
    pair<IterType, bool> inserted = observed_reads.insert(make_pair(read.query_name(), read));
    if(inserted.second) {
        // This is the first time we have seen a read with this name.
        observed_reads[read.query_name()] = read;
    }
    else {
        // We just found an existing read's mate. Good for him/her.
        bd::pair_orientation_flag bdflag = read.bdflag();
        size_t const& index = read.lib_index();
        ++flag_counts[bdflag];
        ++type_library_readcount[bdflag][index];
        type_library_meanspan[bdflag][index] += read.abs_isize();

        ++num_pairs;
        reads_to_free.push_back(read.query_name());
        support_reads.push_back(read);
        support_reads.push_back(inserted.first->second);
        observed_reads.erase(inserted.first);
    }
}
