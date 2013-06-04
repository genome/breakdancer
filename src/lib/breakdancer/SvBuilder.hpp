#pragma once

#include "BasicRegion.hpp"
#include "Read.hpp"

#include <boost/array.hpp>

#include <map>
#include <string>
#include <vector>

struct SvBuilder {
    typedef breakdancer::Read Read;
    typedef std::map<std::string, Read> ObservedReads;

    SvBuilder();

    void observe_read(Read const& read, BasicRegion const& region);
    breakdancer::pair_orientation_flag choose_sv_flag();

// data
    int current_region;
    size_t num_regions;
    int num_pairs;

    breakdancer::PerFlagArray<int>::type flag_counts;
    std::vector<std::string> reads_to_free;

    // number of readpairs per each type/flag (first key) then library (second key)
    breakdancer::PerFlagArray<std::map<std::string, int> >::type type_library_readcount;
    // average ISIZE from BAM records
    breakdancer::PerFlagArray<std::map<std::string, int> >::type type_library_meanspan;

    // vector of readcounts for each type/flag
    std::vector<boost::array<int, 2> > type_orient_counts;

    ObservedReads observed_reads; //unpaired reads

    std::vector<Read> support_reads; //reads supporting the SV
};
