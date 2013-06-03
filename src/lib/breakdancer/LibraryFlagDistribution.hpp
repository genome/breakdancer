#pragma once

#include <stdint.h>
#include <map>
#include <string>
#include <vector>

struct LibraryFlagDistribution {
    LibraryFlagDistribution()
        : index(0)
        , read_count(0)
        , read_counts_by_flag(breakdancer::NUM_ORIENTATION_FLAGS, 0u)
    {
    }

    size_t index;
    std::string name;
    std::string bam_file;
    size_t read_count;
    std::vector<uint32_t> read_counts_by_flag;
};

