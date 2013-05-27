#pragma once

#include <stdint.h>
#include <string>
#include <map>

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
    std::map<breakdancer::pair_orientation_flag, uint32_t> read_counts_by_flag;

    uint32_t get_read_counts_by_flag(breakdancer::pair_orientation_flag flag) const {
        std::map<breakdancer::pair_orientation_flag, uint32_t>::const_iterator fiter = read_counts_by_flag.find(flag);
        if (fiter == read_counts_by_flag.end())
            return 0u;

        return fiter->second;
    }

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
