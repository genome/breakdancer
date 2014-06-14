#pragma once

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include <map>
#include <string>

class BamConfigEntry {
public:
    enum Field {
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


public:
    static Field translate_token(std::string const& tok);
    static std::string const& token_string(Field tok);

    explicit BamConfigEntry(std::string const& line);

    template<typename T>
    bool set_value(Field const& f, T& value) const;

    // throws std::runtime error on failure
    template<typename T>
    void set_required_value(Field const& f, T& value, size_t line_num);

private: // data
    std::map<Field, std::string> _directives;
};

template<typename T>
inline
bool BamConfigEntry::set_value(Field const& f, T& value) const {
    std::map<Field, std::string>::const_iterator found = _directives.find(f);
    if (found != _directives.end()) {
        value = boost::lexical_cast<T>(found->second);
        return true;
    }

    return false;
}

template<typename T>
inline
void BamConfigEntry::set_required_value(Field const& f, T& value, size_t line_num) {
    using boost::format;

    if (!set_value(f, value))
        throw std::runtime_error(str(format(
            "Required field '%1%' not found in config at line %2%!"
            ) % token_string(f) % line_num));
}
