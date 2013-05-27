#pragma once

#include <stdint.h>
#include <map>
#include <string>
#include <vector>

struct LibraryInfo {
    LibraryInfo()
        : index(0)
        , mean_insertsize(0)
        , std_insertsize(0)
        , uppercutoff(0)
        , lowercutoff(0)
        , readlens(0)
        , min_mapping_quality(-1)
        , read_count(0)
        , read_counts_by_flag(breakdancer::NUM_ORIENTATION_FLAGS, 0u)
    {
    }

    size_t index;
    std::string name;
    std::string bam_file;
    float mean_insertsize;
    float std_insertsize;
    float uppercutoff;
    float lowercutoff;
    float readlens;
    int min_mapping_quality;
    size_t read_count;

    // FIXME: ultimately we'll want to renumber the bd flag types and
    // use a simple array for this.
    std::vector<uint32_t> read_counts_by_flag;

    struct Wrapper {
        Wrapper(LibraryInfo const& lib_info)
            : lib_info(lib_info)
        {
        }

        bool operator<(Wrapper const& rhs) const {
            return lib_info.name < rhs.lib_info.name;
        }

        LibraryInfo const& lib_info;
    };

    Wrapper wrap() const {
        return Wrapper(*this);
    }
};
