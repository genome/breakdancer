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
    std::map<std::string, float> mean_insertsize;//global
    std::map<std::string, float> std_insertsize;//global
    std::map<std::string, float> uppercutoff;//global
    std::map<std::string, float> lowercutoff;//global
    std::map<std::string, float> readlens;//global
    std::map<std::string, int> mapQual;// global
    std::map<uint32_t, std::map<std::string,int> > x_readcounts;
    std::map<std::string,std::string> readgroup_library;
    std::map<std::string, std::string> readgroup_platform;
    std::map<std::string, std::string> ReadsOut;

    int d;
    int max_readlen;
};
