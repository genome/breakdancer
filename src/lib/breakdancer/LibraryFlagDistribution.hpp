#pragma once

#include <stdint.h>
#include <map>
#include <string>
#include <vector>
#include "ReadFlags.hpp"

struct LibraryFlagDistribution {
    LibraryFlagDistribution()
        : read_count(0)
        , read_counts_by_flag(breakdancer::NUM_ORIENTATION_FLAGS, 0u)
    {
    }

    size_t read_count;
    std::vector<uint32_t> read_counts_by_flag;
};

