#include "LibraryFlagDistribution.hpp"

LibraryFlagDistribution::LibraryFlagDistribution()
    : read_count(0)
    , read_counts_by_flag(breakdancer::NUM_ORIENTATION_FLAGS, 0u)
{
}

bool LibraryFlagDistribution::operator==(LibraryFlagDistribution const& rhs) const {
    return read_count == rhs.read_count
        && read_counts_by_flag == rhs.read_counts_by_flag;
}


