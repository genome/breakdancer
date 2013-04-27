#pragma once

#include <istream>
#include <map>
#include <stdint.h> // cstdint requires c++11
#include <boost/container/flat_map.hpp>

class Options;

template<typename K, typename V>
struct ConfigMap {
    typedef boost::container::flat_map<K, V> type;
};

class LegacyConfig {
public:
    LegacyConfig(std::istream& in, Options const& opts);

    ConfigMap<std::string, std::string>::type exes;
    ConfigMap<std::string, std::string>::type fmaps;
    ConfigMap<std::string, std::string>::type libmaps;
    ConfigMap<std::string, float>::type mean_insertsize;
    ConfigMap<std::string, float>::type std_insertsize;
    ConfigMap<std::string, float>::type uppercutoff;
    ConfigMap<std::string, float>::type lowercutoff;
    ConfigMap<std::string, float>::type readlens;
    ConfigMap<std::string, int>::type mapQual;
    std::map<uint32_t, std::map<std::string,int> > x_readcounts;
    ConfigMap<std::string, std::string>::type readgroup_library;
    ConfigMap<std::string, std::string>::type readgroup_platform;
    std::map<std::string, std::string> ReadsOut;

    int d;
    int max_readlen;
};
