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
    std::map<std::string, float> std_insertsize;
    std::map<std::string, float> uppercutoff;
    std::map<std::string, float> lowercutoff;
    std::map<std::string, float> readlens;
    std::map<std::string, int> mapQual;
    std::map<uint32_t, std::map<std::string,int> > x_readcounts;
    std::map<std::string,std::string> readgroup_library;
    std::map<std::string, std::string> readgroup_platform;
    std::map<std::string, std::string> ReadsOut;

    int d;
    int max_readlen;
};
