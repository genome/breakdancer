#pragma once

#include "common/ReadFlags.hpp"

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/vector.hpp>

#include <map>
#include <stdint.h>
#include <string>
#include <vector>

class LibraryFlagDistribution {
public:
    friend class boost::serialization::access;

public:
    size_t read_count;
    std::vector<uint32_t> read_counts_by_flag;

    LibraryFlagDistribution();
    bool operator==(LibraryFlagDistribution const& rhs) const;

private:
    template<typename Archive>
    void serialize(Archive& arch, const unsigned int version);
};

template<typename Archive>
void LibraryFlagDistribution::serialize(Archive& arch, const unsigned int version) {
    namespace bs = boost::serialization;
    arch
        & bs::make_nvp("readCount", read_count)
        & bs::make_nvp("readCountsByFlag", read_counts_by_flag)
        ;
}
