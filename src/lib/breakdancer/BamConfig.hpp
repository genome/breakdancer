#pragma once

#include "ReadFlags.hpp"
#include "LibraryInfo.hpp"
#include "utility.hpp"

#include <istream>
#include <map>
#include <stdint.h> // cstdint requires c++11
#include <vector>

class Options;
class BamSummary;
class IBamReader;

class BamConfig {
public:
    BamConfig();
    BamConfig(std::istream& in, Options const& opts);

    ConfigMap<std::string, std::string>::type exes;
    ConfigMap<std::string, std::string>::type fmaps;
    ConfigMap<std::string, std::string>::type readgroup_platform;
    ConfigMap<std::string, std::string>::type ReadsOut;

    int max_read_window_size;
    int max_readlen;

    size_t num_libs() const {
        return _library_info.size();
    }

    size_t num_bams() const {
        return _bam_files.size();
    }

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

    LibraryInfo const& library_info_by_index(size_t idx) const {
        return _library_info[idx];
    }

    LibraryInfo const& library_info_by_name(std::string const& lib) const {
        return _library_info[_lib_names_to_indices.at(lib)];
    }

    uint32_t covered_reference_length() const {
        return _covered_ref_len;
    }

    uint32_t read_count_in_bam(std::string const& key) const {
        return _read_count_per_bam.at(key);
    }

    std::string const& readgroup_library(std::string const& rg) const {
        ConfigMap<std::string, std::string>::type::const_iterator lib = _readgroup_library.find(rg);
        if(lib != _readgroup_library.end())
            return lib->second;
        else
            return fmaps.begin()->second;
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
    ConfigMap<std::string, size_t>::type _lib_names_to_indices;
    std::vector<LibraryInfo> _library_info;
    ConfigMap<std::string, std::string>::type _readgroup_library;
};
