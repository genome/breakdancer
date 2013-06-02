#pragma once

#include "ReadFlags.hpp"

#include <boost/array.hpp>

struct SvEntry {
    typedef breakdancer::pair_orientation_flag FlagType;

    SvEntry(int num_regions, FlagType flag)
        : num_regions(num_regions)
        , flag(flag)
        , chr(_init_minus_one())
        , pos(_init_zero())
        , fwd_read_count(_init_zero())
        , rev_read_count(_init_zero())
    {
    }

    int num_regions;
    FlagType flag;
    boost::array<int, 2> chr;
    boost::array<int, 2> pos;
    boost::array<int, 2> fwd_read_count;
    boost::array<int, 2> rev_read_count;

private:
    static boost::array<int, 2> _init_zero() {
        static boost::array<int, 2> zeros = {{0, 0}};
        return zeros;
    }

    static boost::array<int, 2> _init_minus_one() {
        static boost::array<int, 2> zeros = {{-1, -1}};
        return zeros;
    }
};
