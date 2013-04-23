#pragma once

#include <boost/format.hpp>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <sstream>

class ConfigEntry {
public:
    explicit ConfigEntry(std::string const& line)
        : _line(line)
    {
        parse();
    }

    std::string readgroup;
    std::string platform;
    std::string bam_file;
    std::string library_name;
    double readlen;
    double mean;
    double stddev;
    double lower;
    double upper;
    std::string view_executable;

private: // data
    std::string _line;

private: // functions
    void parse();

    template<typename MapType>
    std::string const& find_required(MapType const& m, char const* k) {
        return find_required(m, std::vector<std::string>(1, k));
    }

    template<typename MapType, typename KeySequenceType>
    std::string const& find_required(MapType const& m, KeySequenceType const& k) {
        std::stringstream tried;

        typedef typename KeySequenceType::const_iterator KIter;
        for (KIter iter = k.begin(); iter != k.end(); ++iter) {
            typename MapType::const_iterator x = m.find(*iter);
            if (x != m.end())
                return x->second;
            if (iter != k.begin())
                tried << ", " << *iter;
        }

        using boost::format;
        throw std::runtime_error(str(format(
            "While parsing line '%1': failed to find a required value. Needed one of %2%"
            ) % _line % tried.str()));
    }
};

class BDConfig {
public:
    explicit BDConfig(std::istream& cfg_stream);

private:
    std::vector<ConfigEntry> _entries;
};
