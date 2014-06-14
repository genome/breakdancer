#include "LibraryFlagDistribution.hpp"

LibraryFlagDistribution::LibraryFlagDistribution()
    : read_count(0)
    , read_counts_by_flag(ReadFlag::NUM_ORIENTATION_FLAGS, 0u)
{
}

bool LibraryFlagDistribution::operator==(LibraryFlagDistribution const& rhs) const {
    return read_count == rhs.read_count
        && read_counts_by_flag == rhs.read_counts_by_flag;
}


void LibraryFlagDistribution::merge(LibraryFlagDistribution const& other) {
    assert(read_counts_by_flag.size() == other.read_counts_by_flag.size());
    for (size_t i = 0; i < read_counts_by_flag.size(); ++i) {
        read_counts_by_flag[i] += other.read_counts_by_flag[i];
    }

    read_count += other.read_count;
}
