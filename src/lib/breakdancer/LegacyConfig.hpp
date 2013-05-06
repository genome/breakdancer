#pragma once

#include "ReadFlags.hpp"

#include <istream>
#include <map>
#include <stdint.h> // cstdint requires c++11
#include <boost/container/flat_map.hpp>

class Options;

template<typename K, typename V>
struct ConfigMap {
    typedef boost::container::flat_map<K, V> type;
};

struct LibraryInfo {
    LibraryInfo()
        : mean_insertsize(0)
        , std_insertsize(0)
        , uppercutoff(0)
        , lowercutoff(0)
        , readlens(0)
        , min_mapping_quality(-1)
    {
    }

    std::string bam_file;
    float mean_insertsize;
    float std_insertsize;
    float uppercutoff;
    float lowercutoff;
    float readlens;
    int min_mapping_quality;

    // FIXME: ultimately we'll want to renumber the bd flag types and
    // use a simple array for this.
    std::map<breakdancer::pair_orientation_flag, uint32_t> read_counts_by_flag;

    uint32_t get_read_counts_by_flag(breakdancer::pair_orientation_flag flag) const {
        std::map<breakdancer::pair_orientation_flag, uint32_t>::const_iterator fiter = read_counts_by_flag.find(flag);
        if (fiter == read_counts_by_flag.end())
            return 0u;

        return fiter->second;
    }
};

class LegacyConfig {
public:
    LegacyConfig();
    LegacyConfig(std::istream& in, Options const& opts);

    ConfigMap<std::string, std::string>::type exes;
    ConfigMap<std::string, std::string>::type fmaps;
    ConfigMap<std::string, LibraryInfo>::type library_info;
    ConfigMap<std::string, std::string>::type readgroup_library;
    ConfigMap<std::string, std::string>::type readgroup_platform;
    ConfigMap<std::string, std::string>::type ReadsOut;

    int max_read_window_size;
    int max_readlen;

    std::string const& platform_for_readgroup(std::string const& readgroup) const {
        using std::string;
        static string illumina("illumina");
        ConfigMap<string, string>::type::const_iterator it = readgroup_platform.find(readgroup);
        if(it != readgroup_platform.end()) {
            return it->second;
        }
        else {
            return illumina;
        }
    }
};
