#pragma once

#include "LibraryConfig.hpp"
#include "common/ReadFlags.hpp"
#include "common/ConfigMap.hpp"

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

#include <istream>
#include <map>
#include <sstream>
#include <stdint.h> // cstdint requires c++11
#include <vector>

class Options;
class IBamReader;

enum ConfigField {
    BAM_FILE,
    LIBRARY_NAME,
    READ_GROUP,
    INSERT_SIZE_MEAN,
    INSERT_SIZE_STDDEV,
    READ_LENGTH,
    INSERT_SIZE_UPPER_CUTOFF,
    INSERT_SIZE_LOWER_CUTOFF,
    MIN_MAP_QUAL,
    SAMPLE_NAME,
    UNKNOWN
};

class ConfigEntry {
public:
    static ConfigField translate_token(std::string const& tok);
    static std::string const& token_string(ConfigField tok);

    explicit ConfigEntry(std::string const& line);

    template<typename T>
    bool set_value(ConfigField const& f, T& value) const {
        std::map<ConfigField, std::string>::const_iterator found = _directives.find(f);
        if (found != _directives.end()) {
            value = boost::lexical_cast<T>(found->second);
            return true;
        }

        return false;
    }

    template<typename T>
    void set_required_value(ConfigField const& f, T& value, size_t line_num) {
        using boost::format;

        if (!set_value(f, value))
            // FIXME: print which field, or the regex
            throw std::runtime_error(str(format(
                "Required field '%1%' not found in config at line %2%!"
                ) % token_string(f) % line_num));
    }

private: // data
    std::map<ConfigField, std::string> _directives;
};


class BamConfig {
public:
    BamConfig();
    BamConfig(std::istream& in, Options const& opts);

    ConfigMap<std::string, std::string>::type ReadsOut;

    int max_read_window_size;
    int max_readlen;

    size_t num_libs() const {
        return _library_config.size();
    }

    size_t num_bams() const {
        return _bam_files.size();
    }

    std::vector<std::string> const& bam_files() const {
        return _bam_files;
    }

    LibraryConfig const& library_config_by_index(size_t idx) const {
        if (idx >= _library_config.size())
            throw std::out_of_range("library index out of range");
        return _library_config[idx];
    }

    LibraryConfig const& library_config_by_name(std::string const& lib) const {
        return _library_config[_lib_names_to_indices.at(lib)];
    }

    std::string const& readgroup_library(std::string const& rg) const {
        ConfigMap<std::string, std::string>::type::const_iterator lib = _readgroup_library.find(rg);
        if(lib != _readgroup_library.end())
            return lib->second;
        else
            return _bam_library.begin()->second;
    }

    template<typename Archive>
    void serialize(Archive& arch, const unsigned int version) {
        arch
            & BOOST_SERIALIZATION_NVP(_bam_library)
            & BOOST_SERIALIZATION_NVP(_bam_files)
            & BOOST_SERIALIZATION_NVP(_lib_names_to_indices)
            & BOOST_SERIALIZATION_NVP(_library_config)
            & BOOST_SERIALIZATION_NVP(_readgroup_library)
            ;
    }

private:
    ConfigMap<std::string, std::string>::type _bam_library;

    // FIXME: we don't need separate vectors for filenames and file summaries
    // we can do with one if the summary can report the name.
    std::vector<std::string> _bam_files;
    ConfigMap<std::string, size_t>::type _lib_names_to_indices;
    std::vector<LibraryConfig> _library_config;
    ConfigMap<std::string, std::string>::type _readgroup_library;
};
