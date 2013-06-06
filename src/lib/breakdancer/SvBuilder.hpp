#pragma once

#include "BasicRegion.hpp"
#include "ReadCountsByLib.hpp"
#include "Read.hpp"

#include <boost/array.hpp>

#include <map>
#include <string>
#include <vector>

struct SvBuilder {
    typedef breakdancer::Read Read;
    typedef std::map<std::string, Read> ObservedReads;
    typedef BasicRegion::iterator_range ReadsRange;

    SvBuilder(int n, BasicRegion const* regions[2], ReadsRange read_ranges[2], int max_readlen);

    void compute_copy_number(ReadCountsByLib const& counts,
        std::map<std::string, float> const& read_density);

// data
    int current_region;
    size_t num_regions;
    int num_pairs;

    breakdancer::PerFlagArray<int>::type flag_counts;
    std::vector<std::string> reads_to_free;

    // number of readpairs per each type/flag (first key) then library (second key)
    breakdancer::PerFlagArray<std::map<std::size_t, int> >::type type_library_readcount;
    // average ISIZE from BAM records
    breakdancer::PerFlagArray<std::map<std::size_t, int> >::type type_library_meanspan;

    ObservedReads observed_reads; //unpaired reads

    std::vector<Read> support_reads; //reads supporting the SV

    breakdancer::pair_orientation_flag flag;
    boost::array<int, 2> chr;
    boost::array<int, 2> pos;
    boost::array<int, 2> fwd_read_count;
    boost::array<int, 2> rev_read_count;

    std::map<std::string, float> copy_number;
    float allele_frequency;

private:
    breakdancer::pair_orientation_flag choose_sv_flag();
    void _observe_read(Read const& read, int region_idx);

    static boost::array<int, 2> _init_zero() {
        static boost::array<int, 2> zeros = {{0, 0}};
        return zeros;
    }
};
