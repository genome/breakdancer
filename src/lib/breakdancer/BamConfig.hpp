#pragma once

#include "ReadFlags.hpp"

#include <boost/container/flat_map.hpp>

#include <istream>
#include <map>
#include <stdint.h> // cstdint requires c++11
#include <vector>

class Options;
class BamSummary;
class IBamReader;

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
        , read_count(0)
    {
    }

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
};

class BamConfig {
public:
    BamConfig();
    BamConfig(std::istream& in, Options const& opts);

    ConfigMap<std::string, std::string>::type exes;
    ConfigMap<std::string, std::string>::type fmaps;
    ConfigMap<std::string, LibraryInfo>::type library_info;
    ConfigMap<std::string, std::string>::type readgroup_library;
    ConfigMap<std::string, std::string>::type readgroup_platform;
    ConfigMap<std::string, std::string>::type ReadsOut;

    int max_read_window_size;
    int max_readlen;

    std::vector<std::string> const& bam_files() const {
        return _bam_files;
    }

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

    uint32_t covered_reference_length() const {
        return _covered_ref_len;
    }

    uint32_t read_count_in_bam(std::string const& key) const {
        return _read_count_per_bam.at(key);
    }

private:
    void _analyze_bam(IBamReader& reader, Options const& opts);

private:
    // FIXME: we don't need separate vectors for filenames and file summaries
    // we can do with one if the summary can report the name.
    std::vector<std::string> _bam_files;
    std::vector<BamSummary*> _bam_summaries;
    ConfigMap<std::string, uint32_t>::type _read_count_per_bam;
    uint32_t _covered_ref_len;
};
