#pragma once

#include "ReadFlags.hpp"

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/vector.hpp>

#include <map>
#include <stdint.h>
#include <string>
#include <vector>

struct LibraryFlagDistribution {
    size_t read_count;
    std::vector<uint32_t> read_counts_by_flag;

    bool operator==(LibraryFlagDistribution const& rhs) const {
        return read_count == rhs.read_count
            && read_counts_by_flag == rhs.read_counts_by_flag;
    }

    LibraryFlagDistribution()
        : read_count(0)
        , read_counts_by_flag(breakdancer::NUM_ORIENTATION_FLAGS, 0u)
    {
    }

    template<typename Archive>
    void serialize(Archive& arch, const unsigned int version) {
        arch
            & BOOST_SERIALIZATION_NVP(read_count)
            & BOOST_SERIALIZATION_NVP(read_counts_by_flag)
            ;
    }
};
