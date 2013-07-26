#pragma once

#include "LibraryConfig.hpp"
#include "common/ReadFlags.hpp"
#include "common/ConfigMap.hpp"

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

#include <cassert>
#include <istream>
#include <map>
#include <sstream>
#include <stdint.h> // cstdint requires c++11
#include <vector>

class Options;
class BamReaderBase;


class BamConfig {
public:
    static const int DEFAULT_MAX_READ_WINDOW_SIZE;

public:
    friend class boost::serialization::access;

public:
    BamConfig();
    BamConfig(std::istream& in, Options const& opts);

    ConfigMap<std::string, std::string>::type ReadsOut;

    int max_read_window_size() const;

    size_t num_libs() const;
    size_t num_bams() const;
    std::vector<std::string> const& bam_files() const;

    LibraryConfig const& library_config_by_index(size_t idx) const;
    LibraryConfig const& library_config_by_name(std::string const& lib) const;
    std::string const& readgroup_library(std::string const& rg) const;

private:
    template<typename Archive>
    void serialize(Archive& arch, const unsigned int version);

private:
    ConfigMap<std::string, std::string>::type _bam_library;

    // FIXME: we don't need separate vectors for filenames and file summaries
    // we can do with one if the summary can report the name.
    std::vector<std::string> _bam_files;
    ConfigMap<std::string, size_t>::type _lib_names_to_indices;
    std::vector<LibraryConfig> _library_config;
    ConfigMap<std::string, std::string>::type _readgroup_library;

    int _max_read_window_size;
};

inline
LibraryConfig const& BamConfig::library_config_by_index(size_t idx) const {
    if (idx >= _library_config.size())
        throw std::out_of_range("library index out of range");
    return _library_config[idx];
}

inline
LibraryConfig const& BamConfig::library_config_by_name(std::string const& lib) const {
    return _library_config[_lib_names_to_indices.at(lib)];
}

inline
std::string const& BamConfig::readgroup_library(std::string const& rg) const {
    ConfigMap<std::string, std::string>::type::const_iterator lib = _readgroup_library.find(rg);
    if(lib != _readgroup_library.end()) {
        return lib->second;
    }
    else {
        assert(!_bam_library.empty());
        return _bam_library.begin()->second;
    }
}

template<typename Archive>
void BamConfig::serialize(Archive& arch, const unsigned int version) {
    namespace bs = boost::serialization;

    arch
        & bs::make_nvp("bamLibrary", _bam_library)
        & bs::make_nvp("bamFiles", _bam_files)
        & bs::make_nvp("libsToIndices", _lib_names_to_indices)
        & bs::make_nvp("libraryConfig", _library_config)
        & bs::make_nvp("readgroupToLibrary", _readgroup_library)
        ;
}
